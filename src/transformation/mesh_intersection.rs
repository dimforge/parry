use crate::math::{Isometry, Point, Real, Vector};
use crate::query::{self, visitors::BoundingVolumeIntersectionsSimultaneousVisitor, PointQuery};
use crate::shape::trimesh::{TopoFace, TopoHalfEdge, TopoVertex, TriMeshTopology};
use crate::shape::{FeatureId, Segment, TriMesh, Triangle};
use core::fmt;
use std::collections::{HashMap, HashSet, VecDeque};

/// Error indicating that a query is not supported between certain shapes
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum MeshIntersectionError {
    MissingTopology,
    TriTriError,
}

impl fmt::Display for MeshIntersectionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingTopology => {
                f.pad("at least on of the meshes is missing its topology information")
            }
            Self::TriTriError => f.pad("internal failure while intersecting two triangles"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for MeshIntersectionError {}

const EPS: Real = 1.0e-5;

/// Computes the intersection of two meshes.
///
/// The meshes must be oriented, have their half-edge topology computed, and must not be self-intersecting.
/// The result mesh vertex coordinates are given in the local-space of `mesh1`.
pub fn intersect_meshes(
    pos1: &Isometry<Real>,
    mesh1: &TriMesh,
    pos2: &Isometry<Real>,
    mesh2: &TriMesh,
) -> Result<Option<TriMesh>, MeshIntersectionError> {
    if mesh1.topology().faces.is_empty() || mesh2.topology().faces.is_empty() {
        return Err(MeshIntersectionError::MissingTopology);
    }

    // NOTE: remove this, used for debugging only.
    mesh1.assert_half_edge_topology_is_valid();
    mesh2.assert_half_edge_topology_is_valid();

    let pos12 = pos1.inv_mul(pos2);
    // 1: collect all the potential triangle-triangle intersections.
    let mut intersections = vec![];
    let mut visitor =
        BoundingVolumeIntersectionsSimultaneousVisitor::with_relative_pos(pos12, |tri1, tri2| {
            intersections.push((*tri1, *tri2));
            true
        });
    mesh1.qbvh().traverse_bvtt(mesh2.qbvh(), &mut visitor);

    let mut deleted_faces1: HashSet<u32> = HashSet::default();
    let mut deleted_faces2: HashSet<u32> = HashSet::default();
    let mut new_indices1 = vec![];
    let mut new_indices2 = vec![];

    for (fid1, fid2) in &intersections {
        let tri1 = mesh1.triangle(*fid1);
        let tri2 = mesh2.triangle(*fid2).transformed(&pos12);

        if triangle_triangle_intersection(&tri1, &tri2).is_some() {
            let _ = deleted_faces1.insert(*fid1);
            let _ = deleted_faces2.insert(*fid2);
        }
    }

    extract_connected_components(&pos12, &mesh1, &mesh2, &deleted_faces1, &mut new_indices1);
    extract_connected_components(
        &pos12.inverse(),
        &mesh2,
        &mesh1,
        &deleted_faces2,
        &mut new_indices2,
    );

    println!(
        "indices1: {}, indices2: {}",
        new_indices1.len(),
        new_indices2.len()
    );

    let mut modif_mesh1 = ModifiedTriMesh::new(mesh1, 1);
    let mut modif_mesh2 = ModifiedTriMesh::new(mesh2, 2);

    // new_indices1.extend(
    //     deleted_faces1
    //         .iter()
    //         .map(|fid| mesh1.indices()[*fid as usize]),
    // );
    // new_indices2.extend(
    //     deleted_faces2
    //         .iter()
    //         .map(|fid| mesh2.indices()[*fid as usize]),
    // );

    cut_and_triangulate_intersections(
        &pos12,
        &deleted_faces1,
        &deleted_faces2,
        &mut modif_mesh1,
        &mut modif_mesh2,
        &mut new_indices1,
        &mut new_indices2,
        &mut intersections,
    );

    // At this point, we know what triangles we want from the first mesh,
    // and the ones we want from the second mesh. Now we need to build the
    // vertex buffer and adjust the indices accordingly.
    let mut new_vertices = vec![];

    // TODO: use Vec instead?
    let mut index_map1 = HashMap::new();
    let mut index_map2 = HashMap::new();
    for idx1 in &mut new_indices1 {
        for k in 0..3 {
            let new_id = *index_map1.entry(idx1[k]).or_insert_with(|| {
                let vtx = modif_mesh1.vertex(idx1[k], &Isometry::identity());
                new_vertices.push(vtx);
                new_vertices.len() - 1
            });
            idx1[k] = new_id as u32;
        }
    }

    for idx2 in &mut new_indices2 {
        for k in 0..3 {
            let new_id = *index_map2.entry(idx2[k]).or_insert_with(|| {
                let vtx = modif_mesh2.vertex(idx2[k], &pos12);
                new_vertices.push(vtx);
                new_vertices.len() - 1
            });
            idx2[k] = new_id as u32;
        }
    }

    new_indices1.append(&mut new_indices2);

    dbg!("Intersection done");

    if !new_indices1.is_empty() {
        Ok(Some(TriMesh::new(new_vertices, new_indices1)))
    } else {
        Ok(None)
    }
}

fn extract_connected_components(
    pos12: &Isometry<Real>,
    mesh1: &TriMesh,
    mesh2: &TriMesh,
    deleted_faces1: &HashSet<u32>,
    new_indices1: &mut Vec<[u32; 3]>,
) {
    let topo1 = mesh1.topology();
    let mut visited: HashSet<u32> = HashSet::default();
    let mut to_visit = vec![];
    dbg!(deleted_faces1.len());

    for face in deleted_faces1 {
        let eid = topo1.faces[*face as usize].half_edge;
        let edge_a = &topo1.half_edges[eid as usize];
        let edge_b = &topo1.half_edges[edge_a.next as usize];
        let edge_c = &topo1.half_edges[edge_b.next as usize];
        let edges = [edge_a, edge_b, edge_c];

        for edge in edges {
            let twin = &topo1.half_edges[edge.twin as usize];
            if !deleted_faces1.contains(&twin.face) {
                let tri1 = mesh1.triangle(twin.face as u32);

                if mesh2.contains_local_point(&pos12.inverse_transform_point(&tri1.a)) {
                    to_visit.push(twin.face);
                }
            }
        }
    }

    // Propagate.
    while let Some(face) = to_visit.pop() {
        if !visited.insert(face) {
            continue; // Already visited.
        }

        new_indices1.push(mesh1.indices()[face as usize]);

        let eid = topo1.faces[face as usize].half_edge;
        let edge_a = &topo1.half_edges[eid as usize];
        let edge_b = &topo1.half_edges[edge_a.next as usize];
        let edge_c = &topo1.half_edges[edge_b.next as usize];
        let edges = [edge_a, edge_b, edge_c];

        for edge in edges {
            let twin = &topo1.half_edges[edge.twin as usize];
            if !deleted_faces1.contains(&twin.face) {
                to_visit.push(twin.face);
            }
        }
    }
}

fn cut_and_triangulate_intersections(
    pos12: &Isometry<Real>,
    deleted_faces1: &HashSet<u32>,
    deleted_faces2: &HashSet<u32>,
    modif_mesh1: &mut ModifiedTriMesh,
    modif_mesh2: &mut ModifiedTriMesh,
    new_indices1: &mut Vec<[u32; 3]>,
    new_indices2: &mut Vec<[u32; 3]>,
    intersections: &mut Vec<(u32, u32)>,
) {
    let mut tris1 = vec![];
    let mut tris2 = vec![];
    let mut intersection_curve = vec![];

    while let Some((i1, i2)) = intersections.pop() {
        tris1.clear();
        tris2.clear();
        modif_mesh1.resolve_triangles(i1, &Isometry::identity(), &mut tris1);
        modif_mesh2.resolve_triangles(i2, &pos12, &mut tris2);

        if tris1.len() > 1 || tris2.len() > 1 {
            // Push the new, resolved pairs.
            for (i1, _, _) in &tris1 {
                for (i2, _, _) in &tris2 {
                    intersections.push((*i1, *i2));
                }
            }
        } else {
            // If we resolved to single triangles, we can compute the intersection.
            let (i1, tri1, idx1) = tris1[0];
            let (i2, tri2, idx2) = tris2[0];

            if let Some(mut inter) = triangle_triangle_intersection(&tri1, &tri2) {
                // First find the half-edge ids for each Feature::Edge.
                // By using these absolute half-edge ids, the splitting process is much simpler.
                if let FeatureId::Edge(eid) = &mut inter.a.f1 {
                    *eid = modif_mesh1.face_half_edges_ids(i1)[*eid as usize];
                }
                if let FeatureId::Edge(eid) = &mut inter.b.f1 {
                    *eid = modif_mesh1.face_half_edges_ids(i1)[*eid as usize];
                }
                if let FeatureId::Edge(eid) = &mut inter.a.f2 {
                    *eid = modif_mesh2.face_half_edges_ids(i2)[*eid as usize];
                }
                if let FeatureId::Edge(eid) = &mut inter.b.f2 {
                    *eid = modif_mesh2.face_half_edges_ids(i2)[*eid as usize];
                }

                let mut inter1 = IntersectionHalfEdge::new();
                let mut inter2 = IntersectionHalfEdge::new();

                // println!(
                //     "Intersection features 1: {:?}, {:?}",
                //     inter.a.f1, inter.b.f1
                // );
                // println!(
                //     "Intersection features 2: {:?}, {:?}",
                //     inter.a.f2, inter.b.f2
                // );
                match (inter.a.f1, inter.b.f1) {
                    (FeatureId::Vertex(vid_a), FeatureId::Vertex(vid_b)) => {
                        // We intersected at two vertices, so the intersection edge is located
                        // on one of the already existing sides of the triangle.
                        let face_eids = modif_mesh1.face_half_edges_ids(i1);
                        inter1.eid = match (vid_a, vid_b) {
                            (0, 1) | (1, 0) => face_eids[0],
                            (1, 2) | (2, 1) => face_eids[1],
                            (2, 0) | (0, 2) => face_eids[2],
                            _ => u32::MAX,
                        };
                    }
                    (FeatureId::Face(_), FeatureId::Face(_)) => {
                        let normal = tri1.scaled_normal();
                        let (new_tris, new_edges) = modif_mesh1.split_face(
                            i1,
                            inter.a.pt,
                            &mut inter1,
                            &Isometry::identity(),
                        );
                        // Now that we split the triangle, we need to figure out on which
                        // triangle piece the second point ends up.
                        let ia_ta = tri1.a - inter.a.pt;
                        let ia_tb = tri1.b - inter.a.pt;
                        let ia_tc = tri1.c - inter.a.pt;
                        let dir_b = inter.b.pt - inter.a.pt;
                        let sgn1 = normal.cross(&ia_ta).dot(&dir_b);
                        let sgn2 = normal.cross(&ia_tb).dot(&dir_b);
                        let sgn3 = normal.cross(&ia_tc).dot(&dir_b);

                        if sgn1.abs() < EPS && dir_b.dot(&ia_ta) >= 0.0 {
                            let _ = modif_mesh1.split_edge(
                                new_edges[0],
                                inter.b.pt,
                                &mut inter1,
                                &Isometry::identity(),
                            );
                        } else if sgn2.abs() < EPS && dir_b.dot(&ia_tb) >= 0.0 {
                            let _ = modif_mesh1.split_edge(
                                new_edges[1],
                                inter.b.pt,
                                &mut inter1,
                                &Isometry::identity(),
                            );
                        } else if sgn3.abs() < EPS && dir_b.dot(&ia_tc) >= 0.0 {
                            let _ = modif_mesh1.split_edge(
                                new_edges[2],
                                inter.b.pt,
                                &mut inter1,
                                &Isometry::identity(),
                            );
                        } else if sgn1 > 0.0 && sgn2 < 0.0 {
                            let _ = modif_mesh1.split_face(
                                new_tris[0],
                                inter.b.pt,
                                &mut inter1,
                                &Isometry::identity(),
                            );
                        } else if sgn2 > 0.0 && sgn3 < 0.0 {
                            let _ = modif_mesh1.split_face(
                                new_tris[1],
                                inter.b.pt,
                                &mut inter1,
                                &Isometry::identity(),
                            );
                        } else {
                            let _ = modif_mesh1.split_face(
                                new_tris[2],
                                inter.b.pt,
                                &mut inter1,
                                &Isometry::identity(),
                            );
                        }
                    }
                    (FeatureId::Edge(eid_a), FeatureId::Edge(eid_b)) if eid_a == eid_b => {
                        // We hit the same edge twice. Therefore, after the first edge split, we need
                        // to determine on which split edge the second intersection point lies.
                        let he = modif_mesh1.half_edge(eid_a);
                        let first_vertex = modif_mesh1.vertex(he.vertex, &Isometry::identity());

                        let new_edges = modif_mesh1.split_edge(
                            eid_a,
                            inter.a.pt,
                            &mut inter1,
                            &Isometry::identity(),
                        );

                        let to_split =
                            if (first_vertex - inter.a.pt).dot(&(inter.b.pt - inter.a.pt)) > 0.0 {
                                new_edges[0]
                            } else {
                                new_edges[1]
                            };

                        let _ = modif_mesh1.split_edge(
                            to_split,
                            inter.b.pt,
                            &mut inter1,
                            &Isometry::identity(),
                        );
                    }
                    _ => {
                        if let FeatureId::Vertex(vid) = inter.a.f1 {
                            inter1.first_vertex = idx1[vid as usize];
                        }
                        if let FeatureId::Vertex(vid) = inter.b.f1 {
                            assert_eq!(inter1.first_vertex, u32::MAX);
                            inter1.first_vertex = idx1[vid as usize];
                        }

                        // Always split the face first, because splitting the edge first would
                        // invalidate `i1`.
                        if matches!(inter.a.f1, FeatureId::Face(_)) {
                            let _ = modif_mesh1.split_face(
                                i1,
                                inter.a.pt,
                                &mut inter1,
                                &Isometry::identity(),
                            );
                        }
                        if matches!(inter.b.f1, FeatureId::Face(_)) {
                            let _ = modif_mesh1.split_face(
                                i1,
                                inter.b.pt,
                                &mut inter1,
                                &Isometry::identity(),
                            );
                        }

                        if let FeatureId::Edge(eid) = inter.a.f1 {
                            let _ = modif_mesh1.split_edge(
                                eid,
                                inter.a.pt,
                                &mut inter1,
                                &Isometry::identity(),
                            );
                        }

                        if let FeatureId::Edge(eid) = inter.b.f1 {
                            let _ = modif_mesh1.split_edge(
                                eid,
                                inter.b.pt,
                                &mut inter1,
                                &Isometry::identity(),
                            );
                        }
                    }
                }

                match (inter.a.f2, inter.b.f2) {
                    (FeatureId::Vertex(vid_a), FeatureId::Vertex(vid_b)) => {
                        // We intersected at two vertices, so the intersection edge is located
                        // on one of the already existing sides of the triangle.
                        let face_eids = modif_mesh2.face_half_edges_ids(i2);
                        inter2.eid = match (vid_a, vid_b) {
                            (0, 1) | (1, 0) => face_eids[0],
                            (1, 2) | (2, 1) => face_eids[1],
                            (2, 0) | (0, 2) => face_eids[2],
                            _ => u32::MAX,
                        };
                    }
                    (FeatureId::Face(_), FeatureId::Face(_)) => {
                        let normal = tri2.scaled_normal();
                        let (new_tris, new_edges) =
                            modif_mesh2.split_face(i2, inter.a.pt, &mut inter2, &pos12);
                        // Now that we split the triangle, we need to figure out on which
                        // triangle piece the second point ends up.
                        let ia_ta = tri2.a - inter.a.pt;
                        let ia_tb = tri2.b - inter.a.pt;
                        let ia_tc = tri2.c - inter.a.pt;
                        let dir_b = inter.b.pt - inter.a.pt;
                        let sgn1 = normal.cross(&ia_ta).dot(&dir_b);
                        let sgn2 = normal.cross(&ia_tb).dot(&dir_b);
                        let sgn3 = normal.cross(&ia_tc).dot(&dir_b);

                        if sgn1.abs() < EPS && dir_b.dot(&ia_ta) >= 0.0 {
                            let _ = modif_mesh2.split_edge(
                                new_edges[0],
                                inter.b.pt,
                                &mut inter2,
                                &pos12,
                            );
                        } else if sgn2.abs() < EPS && dir_b.dot(&ia_tb) >= 0.0 {
                            let _ = modif_mesh2.split_edge(
                                new_edges[1],
                                inter.b.pt,
                                &mut inter2,
                                &pos12,
                            );
                        } else if sgn3.abs() < EPS && dir_b.dot(&ia_tc) >= 0.0 {
                            let _ = modif_mesh2.split_edge(
                                new_edges[2],
                                inter.b.pt,
                                &mut inter2,
                                &pos12,
                            );
                        } else if sgn1 > 0.0 && sgn2 < 0.0 {
                            let _ = modif_mesh2.split_face(
                                new_tris[0],
                                inter.b.pt,
                                &mut inter2,
                                &pos12,
                            );
                        } else if sgn2 > 0.0 && sgn3 < 0.0 {
                            let _ = modif_mesh2.split_face(
                                new_tris[1],
                                inter.b.pt,
                                &mut inter2,
                                &pos12,
                            );
                        } else {
                            let _ = modif_mesh2.split_face(
                                new_tris[2],
                                inter.b.pt,
                                &mut inter2,
                                &pos12,
                            );
                        }
                    }
                    (FeatureId::Edge(eid_a), FeatureId::Edge(eid_b)) if eid_a == eid_b => {
                        // We hit the same edge twice. Therefore, after the first edge split, we need
                        // to determine on which split edge the second intersection point lies.
                        let he = modif_mesh2.half_edge(eid_a);
                        let first_vertex = modif_mesh2.vertex(he.vertex, &pos12);

                        let new_edges =
                            modif_mesh2.split_edge(eid_a, inter.a.pt, &mut inter2, &pos12);

                        let to_split =
                            if (first_vertex - inter.a.pt).dot(&(inter.b.pt - inter.a.pt)) > 0.0 {
                                new_edges[0]
                            } else {
                                new_edges[1]
                            };

                        let _ = modif_mesh2.split_edge(to_split, inter.b.pt, &mut inter2, &pos12);
                    }
                    _ => {
                        if let FeatureId::Vertex(vid) = inter.a.f2 {
                            inter2.first_vertex = idx2[vid as usize];
                        }
                        if let FeatureId::Vertex(vid) = inter.b.f2 {
                            assert_eq!(inter2.first_vertex, u32::MAX);
                            inter2.first_vertex = idx2[vid as usize];
                        }

                        // Always split the face first, because splitting the edge first would
                        // invalidate `i2`.
                        if matches!(inter.a.f2, FeatureId::Face(_)) {
                            let _ = modif_mesh2.split_face(i2, inter.a.pt, &mut inter2, &pos12);
                        }
                        if matches!(inter.b.f2, FeatureId::Face(_)) {
                            let _ = modif_mesh2.split_face(i2, inter.b.pt, &mut inter2, &pos12);
                        }

                        if let FeatureId::Edge(eid) = inter.a.f2 {
                            let _ = modif_mesh2.split_edge(eid, inter.a.pt, &mut inter2, &pos12);
                        }

                        if let FeatureId::Edge(eid) = inter.b.f2 {
                            let _ = modif_mesh2.split_edge(eid, inter.b.pt, &mut inter2, &pos12);
                        }
                    }
                }

                if inter1.eid != u32::MAX && inter2.eid != u32::MAX {
                    intersection_curve.push((inter1.eid, inter2.eid));
                }
            }
        }
    }

    extract_result(
        &pos12,
        &modif_mesh1,
        &modif_mesh2,
        &deleted_faces1,
        &deleted_faces2,
        &intersection_curve,
        new_indices1,
        new_indices2,
    );
}

fn extract_result(
    pos12: &Isometry<Real>,
    mesh1: &ModifiedTriMesh,
    mesh2: &ModifiedTriMesh,
    deleted_faces1: &HashSet<u32>,
    deleted_faces2: &HashSet<u32>,
    intersection_curve: &[(u32, u32)],
    new_indices1: &mut Vec<[u32; 3]>,
    new_indices2: &mut Vec<[u32; 3]>,
) {
    for (fid, added) in mesh1.added_faces.iter().enumerate() {
        if mesh1.face_splits.contains_key(&(fid as u32)) {
            continue;
        }

        if added[0] != u32::MAX {
            let tri = Triangle::new(
                mesh1.vertex(added[0], &Isometry::identity()),
                mesh1.vertex(added[1], &Isometry::identity()),
                mesh1.vertex(added[2], &Isometry::identity()),
            )
            .transformed(&pos12.inverse());

            let base_len = mesh1.trimesh.vertices().len() as u32;

            if (added[0] >= base_len && added[1] >= base_len && added[1] >= base_len)
                || (added[0] < base_len && added[2] < base_len && added[2] < base_len)
            {
                if mesh2.trimesh.contains_local_point(&tri.center()) {
                    new_indices1.push(*added);
                }
            } else if (added[0] >= base_len || mesh2.trimesh.contains_local_point(&tri.a))
                && (added[1] >= base_len || mesh2.trimesh.contains_local_point(&tri.b))
                && (added[2] >= base_len || mesh2.trimesh.contains_local_point(&tri.c))
            {
                new_indices1.push(*added);
            }
        }
    }

    for (fid, added) in mesh2.added_faces.iter().enumerate() {
        if mesh2.face_splits.contains_key(&(fid as u32)) {
            continue;
        }

        if added[0] != u32::MAX {
            let tri = Triangle::new(
                mesh2.vertex(added[0], pos12),
                mesh2.vertex(added[1], pos12),
                mesh2.vertex(added[2], pos12),
            );

            let base_len = mesh2.trimesh.vertices().len() as u32;

            if (added[0] >= base_len && added[1] >= base_len && added[1] >= base_len)
                || (added[0] < base_len && added[2] < base_len && added[2] < base_len)
            {
                if mesh1.trimesh.contains_local_point(&tri.center()) {
                    new_indices2.push(*added);
                }
            } else if (added[0] >= base_len || mesh1.trimesh.contains_local_point(&tri.a))
                && (added[1] >= base_len || mesh1.trimesh.contains_local_point(&tri.b))
                && (added[2] >= base_len || mesh1.trimesh.contains_local_point(&tri.c))
            {
                new_indices2.push(*added);
            }
        }
    }
    return;
    // Identify all the propagation starting faces.
    for (eid1, eid2) in intersection_curve {
        let edge1 = mesh1.half_edge(*eid1);
        let twin_edge1 = mesh1.half_edge(edge1.twin);

        let edge2 = mesh2.half_edge(*eid2);
        let twin_edge2 = mesh2.half_edge(edge2.twin);

        let tri1 = mesh1.original_face_vertices(edge1.face, &Isometry::identity());
        let tri2 = mesh2.original_face_vertices(edge2.face, pos12);
        let normal1 = tri1.scaled_normal();
        let normal2 = tri2.scaled_normal();

        let edge_dir1 = mesh1.vertex(twin_edge1.vertex, &Isometry::identity())
            - mesh1.vertex(edge1.vertex, &Isometry::identity());
        let edge_dir2 = mesh2.vertex(twin_edge2.vertex, pos12) - mesh2.vertex(edge2.vertex, pos12);

        let ref_dir1 = normal1.cross(&edge_dir1);
        let ref_dir2 = normal2.cross(&edge_dir2);

        let dot1 = ref_dir1.dot(&normal2);
        let dot2 = ref_dir2.dot(&normal1);

        {
            let to_add1 = if dot1 < 0.0 { *eid1 } else { edge1.twin };
            let edge = mesh1.half_edge(to_add1);
            if mesh1.face_splits.contains_key(&edge.face) {
                continue;
            }
            new_indices1.push(mesh1.face_indices(edge.face));
        }

        {
            let to_add2 = if dot2 < 0.0 { *eid2 } else { edge2.twin };
            let edge = mesh2.half_edge(to_add2);
            if mesh2.face_splits.contains_key(&edge.face) {
                continue;
            }
            new_indices2.push(mesh2.face_indices(edge.face));
        }
    }
}

#[derive(Copy, Clone, Debug)]
struct IntersectionHalfEdge {
    first_vertex: u32,
    second_vertex: u32,
    eid: u32,
}

impl IntersectionHalfEdge {
    fn new() -> Self {
        Self {
            first_vertex: u32::MAX,
            second_vertex: u32::MAX,
            eid: u32::MAX,
        }
    }
}

struct ModifiedTriMesh<'a> {
    tag: u8, // NOTE: only for debugging.
    trimesh: &'a TriMesh,
    // NOTE: the added vertices are in the local-space of the first trimesh.
    added_vertices: Vec<Point<Real>>,
    added_faces: Vec<[u32; 3]>,
    // The original triangle this face is split from.
    original_faces: Vec<u32>,
    added_topo: TriMeshTopology,
    half_edge_substitutes: HashMap<u32, u32>,
    // TODO: we need a hash-map only for the first time a face
    //       is split. If it’s split more than once, the subsequent splits
    //       could be stored in a Vec with direct indexed access.
    face_splits: HashMap<u32, u32>,
    edge_splits: HashSet<u32>, // For debugging
}

impl<'a> ModifiedTriMesh<'a> {
    fn new(trimesh: &'a TriMesh, tag: u8) -> Self {
        Self {
            tag,
            trimesh,
            added_vertices: vec![],
            added_faces: vec![],
            original_faces: vec![],
            added_topo: TriMeshTopology::default(),
            half_edge_substitutes: HashMap::default(),
            face_splits: HashMap::default(),
            edge_splits: HashSet::default(),
        }
    }

    // NOTE: this is mostly used for debugging, in order to visualize the whay faces have been cut.
    fn into_trimesh(mut self, pos: &Isometry<Real>) -> TriMesh {
        let new_indices: Vec<_> = self
            .trimesh
            .indices()
            .into_iter()
            .copied()
            .chain(self.added_faces.into_iter())
            .enumerate()
            .filter(|(i, idx)| !self.face_splits.contains_key(&(*i as u32)) && idx[0] != u32::MAX)
            .map(|(_, idx)| idx)
            .collect();

        let mut new_vertices: Vec<_> = self.trimesh.vertices().iter().map(|pt| pos * pt).collect();
        new_vertices.append(&mut self.added_vertices);

        TriMesh::new(new_vertices, new_indices)
    }

    fn resolve_triangles(
        &self,
        id: u32,
        pos: &Isometry<Real>,
        out: &mut Vec<(u32, Triangle, [u32; 3])>,
    ) {
        if let Some(face_splits) = self.face_splits.get(&id) {
            self.resolve_triangles(*face_splits, pos, out);
            self.resolve_triangles(*face_splits + 1, pos, out);
            self.resolve_triangles(*face_splits + 2, pos, out);
        } else if let Some((tri, idx)) = self.face_vertices(id, pos) {
            out.push((id, tri, idx));
        }
    }

    fn vertex(&self, vid: u32, pos: &Isometry<Real>) -> Point<Real> {
        let vertices = self.trimesh.vertices();
        if (vid as usize) < vertices.len() {
            pos * vertices[vid as usize]
        } else {
            self.added_vertices[vid as usize - vertices.len()]
        }
    }

    fn original_face_id(&self, fid: u32) -> u32 {
        let base_len = self.trimesh.indices().len();
        if fid >= base_len as u32 {
            self.original_faces[fid as usize - base_len]
        } else {
            fid
        }
    }

    fn original_face_vertices(&self, fid: u32, pos: &Isometry<Real>) -> Triangle {
        self.trimesh
            .triangle(self.original_face_id(fid))
            .transformed(pos)
    }

    fn face_vertices(&self, fid: u32, pos: &Isometry<Real>) -> Option<(Triangle, [u32; 3])> {
        let idx = self.face_indices(fid);
        // If this is not an invalid triangle (added by edge-splitting), push it.
        if idx[0] != u32::MAX {
            Some((
                Triangle::new(
                    self.vertex(idx[0], pos),
                    self.vertex(idx[1], pos),
                    self.vertex(idx[2], pos),
                ),
                idx,
            ))
        } else {
            None
        }
    }

    fn half_edge_base_id(&self) -> u32 {
        (self.trimesh.topology().half_edges.len() + self.added_topo.half_edges.len()) as u32
    }

    fn vertex_base_id(&self) -> u32 {
        (self.trimesh.vertices().len() + self.added_vertices.len()) as u32
    }

    fn face_base_id(&self) -> u32 {
        (self.trimesh.topology().faces.len() + self.added_topo.faces.len()) as u32
    }

    fn face_half_edges_ids(&self, fid: u32) -> [u32; 3] {
        assert!(!self.face_splits.contains_key(&fid));
        let prev_face_len = self.trimesh.topology().faces.len();
        let first_half_edge = if let Some(face) = self.trimesh.topology().faces.get(fid as usize) {
            face.half_edge
        } else {
            self.added_topo.faces[fid as usize - prev_face_len].half_edge
        };

        let mut result = [first_half_edge; 3];
        for k in 1..3 {
            let half_edge = self.half_edge(result[k - 1]);
            result[k] = half_edge.next;
        }
        result
    }

    fn face_indices_and_edge_index(&self, half_edge: &TopoHalfEdge) -> ([u32; 3], u8) {
        let indices = self.face_indices(half_edge.face);
        for i in 0..3 {
            if indices[i] == half_edge.vertex {
                return (indices, i as u8);
            }
        }

        unreachable!();
    }

    fn face_indices(&self, fid: u32) -> [u32; 3] {
        assert!(!self.face_splits.contains_key(&fid));
        if let Some(idx) = self.trimesh.indices().get(fid as usize) {
            *idx
        } else {
            let prev_face_len = self.trimesh.topology().faces.len();
            self.added_faces[fid as usize - prev_face_len]
        }
    }

    fn half_edge(&self, id: u32) -> &TopoHalfEdge {
        let id = self.half_edge_substitutes.get(&id).copied().unwrap_or(id);
        if let Some(hedge) = self.trimesh.topology().half_edges.get(id as usize) {
            hedge
        } else {
            let prev_half_edges = self.trimesh.topology().half_edges.len();
            &self.added_topo.half_edges[id as usize - prev_half_edges]
        }
    }

    fn substitute_half_edge(&mut self, id_to_substitute: u32) -> (&mut TopoHalfEdge, u32) {
        let mut id = self
            .half_edge_substitutes
            .get(&id_to_substitute)
            .copied()
            .unwrap_or(id_to_substitute);
        let base_len = self.trimesh.topology().half_edges.len();

        if let Some(hedge) = self.trimesh.topology().half_edges.get(id as usize) {
            id = self.half_edge_base_id();
            self.added_topo.half_edges.push(*hedge);
            let _ = self.half_edge_substitutes.insert(id_to_substitute, id);
        }

        (&mut self.added_topo.half_edges[id as usize - base_len], id)
    }

    fn assert_face_isnt_degenerate(&self, idx: [u32; 3], pos: &Isometry<Real>) {
        let tri = Triangle::new(
            self.vertex(idx[0], pos),
            self.vertex(idx[1], pos),
            self.vertex(idx[2], pos),
        );
        if tri.area() == 0.0 {
            // dbg!(tri);
        }
        // assert_ne!(tri.area(), 0.0);
    }

    fn push_original_face(&mut self, fid: u32) {
        let base_len = self.trimesh.indices().len();
        if (fid as usize) < base_len {
            self.original_faces.push(fid);
        } else {
            self.original_faces
                .push(self.original_faces[fid as usize - base_len]);
        }
    }

    fn split_edge(
        &mut self,
        half_edge_id: u32,
        new_pt: Point<Real>,
        intersection: &mut IntersectionHalfEdge,
        pos: &Isometry<Real>,
    ) -> [u32; 2] {
        // println!("Splitting half-egde: {:?}", half_edge_id);
        let vertex_base_id = self.vertex_base_id();
        let face_base_id = self.face_base_id();
        let half_edge_base_id = self.half_edge_base_id();
        let new_pt_id = self.vertex_base_id();

        let half_edge = *self.half_edge(half_edge_id);
        let fid = half_edge.face;
        let (face_indices, edge) = self.face_indices_and_edge_index(&half_edge);

        let twin_half_edge = *self.half_edge(half_edge.twin);
        let twin_fid = twin_half_edge.face;
        let (twin_face_indices, twin_edge) = self.face_indices_and_edge_index(&twin_half_edge);

        assert!(self.edge_splits.insert(half_edge_id));
        assert!(self.edge_splits.insert(half_edge.twin));

        // Grab the indices of the half-edges that will survive the split (one half-edge per triangle after-split).
        let next_half_edge = half_edge.next;
        let next_next_half_edge = self.half_edge(next_half_edge).next;

        let twin_next_half_edge = twin_half_edge.next;
        let twin_next_next_half_edge = self.half_edge(twin_next_half_edge).next;

        // Add 1 vertex.
        self.added_vertices.push(new_pt);
        self.added_topo.vertices.push(TopoVertex {
            half_edge: half_edge_base_id,
        });

        // Add 4 faces.
        assert_eq!(
            face_indices[(edge as usize + 0) % 3],
            twin_face_indices[(twin_edge as usize + 1) % 3]
        );
        assert_eq!(
            face_indices[(edge as usize + 1) % 3],
            twin_face_indices[(twin_edge as usize + 0) % 3]
        );

        let new_face_indices = [
            [
                face_indices[(edge as usize + 2) % 3],
                face_indices[(edge as usize + 0) % 3],
                new_pt_id,
            ],
            [
                face_indices[(edge as usize + 1) % 3],
                face_indices[(edge as usize + 2) % 3],
                new_pt_id,
            ],
            [u32::MAX; 3], // We always split triangles into 3 triangles. So, for the third tri, use invalid indices.
        ];
        let new_twin_face_indices = [
            [
                twin_face_indices[(twin_edge as usize + 2) % 3],
                twin_face_indices[(twin_edge as usize + 0) % 3],
                new_pt_id,
            ],
            [
                twin_face_indices[(twin_edge as usize + 1) % 3],
                twin_face_indices[(twin_edge as usize + 2) % 3],
                new_pt_id,
            ],
            [u32::MAX; 3],
        ];

        // println!("[{}] Splitting edge face {fid}", self.tag);
        // println!("[{}] Splitting edge face {twin_fid}", self.tag);
        assert!(
            self.face_splits.insert(fid, face_base_id).is_none(),
            "Found duplicate face split."
        );
        assert!(
            self.face_splits
                .insert(twin_fid, face_base_id + 3)
                .is_none(),
            "Found duplicate face split."
        );

        // self.assert_face_isnt_degenerate(new_face_indices[0], pos);
        // self.assert_face_isnt_degenerate(new_face_indices[1], pos);
        // self.assert_face_isnt_degenerate(new_twin_face_indices[0], pos);
        // self.assert_face_isnt_degenerate(new_twin_face_indices[1], pos);

        for k in 0..3 {
            self.added_faces.push(new_face_indices[k]);
            self.push_original_face(fid);
        }
        for k in 0..3 {
            self.added_faces.push(new_twin_face_indices[k]);
            self.push_original_face(twin_fid);
        }
        self.added_topo.faces.push(TopoFace {
            half_edge: next_next_half_edge,
        });
        self.added_topo.faces.push(TopoFace {
            half_edge: next_half_edge,
        });
        self.added_topo.faces.push(TopoFace {
            half_edge: u32::MAX,
        });
        self.added_topo.faces.push(TopoFace {
            half_edge: twin_next_next_half_edge,
        });
        self.added_topo.faces.push(TopoFace {
            half_edge: twin_next_half_edge,
        });
        self.added_topo.faces.push(TopoFace {
            half_edge: u32::MAX,
        });

        // Add the 8 new half-edges (we add them face-after-face, in CCW order).
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: half_edge_base_id + 1,
            twin: half_edge_base_id + 7,
            vertex: new_face_indices[0][1],
            face: face_base_id + 0,
        });
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: next_next_half_edge,
            twin: half_edge_base_id + 2,
            vertex: new_pt_id,
            face: face_base_id + 0,
        });
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: half_edge_base_id + 3,
            twin: half_edge_base_id + 1,
            vertex: new_face_indices[1][1],
            face: face_base_id + 1,
        });
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: next_half_edge,
            twin: half_edge_base_id + 4,
            vertex: new_pt_id,
            face: face_base_id + 1,
        });
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: half_edge_base_id + 5,
            twin: half_edge_base_id + 3,
            vertex: new_twin_face_indices[0][1],
            face: face_base_id + 3,
        });
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: twin_next_next_half_edge,
            twin: half_edge_base_id + 6,
            vertex: new_pt_id,
            face: face_base_id + 3,
        });
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: half_edge_base_id + 7,
            twin: half_edge_base_id + 5,
            vertex: new_twin_face_indices[1][1],
            face: face_base_id + 4,
        });
        self.added_topo.half_edges.push(TopoHalfEdge {
            next: twin_next_half_edge,
            twin: half_edge_base_id + 0,
            vertex: new_pt_id,
            face: face_base_id + 4,
        });

        // Modify the surviving adjascent half-edges:
        // - update `next`.
        // - update `face`.
        let (next_next, _) = self.substitute_half_edge(next_next_half_edge);
        next_next.next = half_edge_base_id + 0;
        next_next.face = face_base_id + 0;

        let (next, _) = self.substitute_half_edge(next_half_edge);
        next.next = half_edge_base_id + 2;
        next.face = face_base_id + 1;

        let (twin_next_next, _) = self.substitute_half_edge(twin_next_next_half_edge);
        twin_next_next.next = half_edge_base_id + 4;
        twin_next_next.face = face_base_id + 3;

        let (twin_next, _) = self.substitute_half_edge(twin_next_half_edge);
        twin_next.next = half_edge_base_id + 6;
        twin_next.face = face_base_id + 4;

        if intersection.first_vertex == u32::MAX {
            intersection.first_vertex = vertex_base_id;
        } else {
            intersection.second_vertex = vertex_base_id;
            let base = self.added_topo.half_edges.len() - 8;
            for k in 0..8 {
                // if self.added_topo.half_edges[base + k].vertex == intersection.first_vertex {
                if self.half_edge(half_edge_base_id + k).vertex == intersection.first_vertex {
                    intersection.eid = half_edge_base_id + k as u32;
                }
            }
        }

        [half_edge_base_id, half_edge_base_id + 3]
    }

    fn split_face(
        &mut self,
        fid: u32,
        new_pt: Point<Real>,
        intersection: &mut IntersectionHalfEdge,
        pos: &Isometry<Real>,
    ) -> ([u32; 3], [u32; 3]) {
        // println!("[{}] Splitting face {fid}", self.tag);
        let vertex_base_id = self.vertex_base_id();
        let face_base_id = self.face_base_id();
        let half_edge_base_id = self.half_edge_base_id();
        let idx = self.face_indices(fid);
        let new_pt_id = self.vertex_base_id();
        let half_edge_ids = self.face_half_edges_ids(fid);

        // Add 1 vertex.
        self.added_vertices.push(new_pt);
        self.added_topo.vertices.push(TopoVertex {
            half_edge: half_edge_base_id,
        });

        // Add 3 faces.
        assert!(
            self.face_splits.insert(fid, face_base_id).is_none(),
            "Found duplicate face split."
        );

        let tri1 = [idx[0], idx[1], new_pt_id];
        let tri2 = [idx[1], idx[2], new_pt_id];
        let tri3 = [idx[2], idx[0], new_pt_id];
        self.added_faces.push(tri1);
        self.added_faces.push(tri2);
        self.added_faces.push(tri3);
        self.push_original_face(fid);
        self.push_original_face(fid);
        self.push_original_face(fid);

        self.assert_face_isnt_degenerate(tri1, pos);
        self.assert_face_isnt_degenerate(tri2, pos);
        self.assert_face_isnt_degenerate(tri3, pos);

        self.added_topo.faces.push(TopoFace {
            half_edge: half_edge_ids[0],
        });
        self.added_topo.faces.push(TopoFace {
            half_edge: half_edge_ids[1],
        });
        self.added_topo.faces.push(TopoFace {
            half_edge: half_edge_ids[2],
        });

        // Add 3 half-edges starting at new_point.
        for k in 0u32..3 {
            self.added_topo.half_edges.push(TopoHalfEdge {
                next: half_edge_ids[k as usize],
                twin: half_edge_base_id + k + 3,
                vertex: new_pt_id,
                face: face_base_id + k,
            });
        }

        // Add 3 half-edges ending at new_point.
        for k in 0u32..3 {
            self.added_topo.half_edges.push(TopoHalfEdge {
                next: half_edge_base_id + (k + 2) % 3,
                twin: half_edge_base_id + k,
                vertex: idx[k as usize],
                face: face_base_id + (k + 2) % 3,
            });
        }

        // From existing half-edges:
        // - Modify the `next` index.
        // - Modify the `face` index.
        let nexts = [4, 5, 3];
        for k in 0u32..3 {
            let (hedge, _) = self.substitute_half_edge(half_edge_ids[k as usize]);
            hedge.face = face_base_id + k;
            hedge.next = half_edge_base_id + nexts[k as usize];
        }

        if intersection.first_vertex == u32::MAX {
            intersection.first_vertex = vertex_base_id;
        } else {
            intersection.second_vertex = vertex_base_id;
            // let base = self.added_topo.half_edges.len() - 6;
            for k in 0..6 {
                // if self.added_topo.half_edges[base + k].vertex == intersection.first_vertex {
                if self.half_edge(half_edge_base_id + k).vertex == intersection.first_vertex {
                    assert!(k >= 3);
                    intersection.eid = half_edge_base_id + k as u32;
                }
            }
        }

        (
            [face_base_id, face_base_id + 1, face_base_id + 2],
            [
                half_edge_base_id,
                half_edge_base_id + 1,
                half_edge_base_id + 2,
            ],
        )
    }
}

#[derive(Copy, Clone, Debug, Default)]
struct TriangleTriangleIntersectionPoint {
    pt: Point<Real>,
    f1: FeatureId,
    f2: FeatureId,
}

#[derive(Copy, Clone, Debug, Default)]
struct TriangleTriangleIntersection {
    a: TriangleTriangleIntersectionPoint,
    b: TriangleTriangleIntersectionPoint,
}

fn triangle_triangle_intersection(
    tri1: &Triangle,
    tri2: &Triangle,
) -> Option<TriangleTriangleIntersection> {
    let normal1 = tri1.scaled_normal();
    let normal2 = tri2.scaled_normal();

    let intersection_dir = normal1.cross(&normal2).normalize();

    let mut range1 = [
        (Real::MAX, Point::origin(), FeatureId::Unknown),
        (-Real::MAX, Point::origin(), FeatureId::Unknown),
    ];
    let mut range2 = [
        (Real::MAX, Point::origin(), FeatureId::Unknown),
        (-Real::MAX, Point::origin(), FeatureId::Unknown),
    ];

    let hits1 = [
        segment_plane_intersection(&tri2.a, &normal2, &Segment::new(tri1.a, tri1.b), 0, (0, 1))
            .map(|(p, feat)| (intersection_dir.dot(&p.coords), p, feat)),
        segment_plane_intersection(&tri2.a, &normal2, &Segment::new(tri1.b, tri1.c), 1, (1, 2))
            .map(|(p, feat)| (intersection_dir.dot(&p.coords), p, feat)),
        segment_plane_intersection(&tri2.a, &normal2, &Segment::new(tri1.c, tri1.a), 2, (2, 0))
            .map(|(p, feat)| (intersection_dir.dot(&p.coords), p, feat)),
    ];

    for k in 0..3 {
        if let Some(hit1) = &hits1[k] {
            if hit1.0 < range1[0].0 {
                range1[0] = *hit1;
            }
            if hit1.0 > range1[1].0 {
                range1[1] = *hit1;
            }
        }
    }

    if range1[0].0 >= range1[1].0 {
        // The first triangle doesn’t intersect the second plane.
        return None;
    }

    let hits2 = [
        segment_plane_intersection(&tri1.a, &normal1, &Segment::new(tri2.a, tri2.b), 0, (0, 1))
            .map(|(p, feat)| (intersection_dir.dot(&p.coords), p, feat)),
        segment_plane_intersection(&tri1.a, &normal1, &Segment::new(tri2.b, tri2.c), 1, (1, 2))
            .map(|(p, feat)| (intersection_dir.dot(&p.coords), p, feat)),
        segment_plane_intersection(&tri1.a, &normal1, &Segment::new(tri2.c, tri2.a), 2, (2, 0))
            .map(|(p, feat)| (intersection_dir.dot(&p.coords), p, feat)),
    ];

    for k in 0..3 {
        if let Some(hit2) = &hits2[k] {
            if hit2.0 < range2[0].0 {
                range2[0] = *hit2;
            }
            if hit2.0 > range2[1].0 {
                range2[1] = *hit2;
            }
        }
    }

    if range2[0].0 >= range2[1].0 {
        // The second triangle doesn’t intersect the first plane.
        return None;
    }

    if range1[1].0 <= range2[0].0 + EPS || range2[1].0 <= range1[0].0 + EPS {
        // The two triangles intersect each others’ plane, but these intersections are disjoint.
        return None;
    }

    let edge_between = |a, b| match (a, b) {
        (0, 1) | (1, 0) => FeatureId::Edge(0),
        (1, 2) | (2, 1) => FeatureId::Edge(1),
        (2, 0) | (0, 2) => FeatureId::Edge(2),
        _ => FeatureId::Edge(a),
    };

    let inter_f1 = match (range1[0].2, range1[1].2) {
        (FeatureId::Vertex(a), FeatureId::Vertex(b)) => edge_between(a, b),
        (FeatureId::Vertex(v), FeatureId::Edge(e)) | (FeatureId::Edge(e), FeatureId::Vertex(v)) => {
            if e == (v + 1) % 3 {
                FeatureId::Face(0)
            } else {
                FeatureId::Edge(e)
            }
        }
        _ => FeatureId::Face(0),
    };
    let inter_f2 = match (range2[0].2, range2[1].2) {
        (FeatureId::Vertex(a), FeatureId::Vertex(b)) => edge_between(a, b),
        (FeatureId::Vertex(v), FeatureId::Edge(e)) | (FeatureId::Edge(e), FeatureId::Vertex(v)) => {
            if e == (v + 1) % 3 {
                FeatureId::Face(0)
            } else {
                FeatureId::Edge(e)
            }
        }
        _ => FeatureId::Face(0),
    };

    let a = if range2[0].0 > range1[0].0 + EPS {
        TriangleTriangleIntersectionPoint {
            pt: range2[0].1,
            f1: inter_f1,
            f2: range2[0].2,
        }
    } else if range2[0].0 < range1[0].0 - EPS {
        TriangleTriangleIntersectionPoint {
            pt: range1[0].1,
            f1: range1[0].2,
            f2: inter_f2,
        }
    } else {
        TriangleTriangleIntersectionPoint {
            pt: range1[0].1,
            f1: range1[0].2,
            f2: range2[0].2,
        }
    };

    let b = if range2[1].0 < range1[1].0 - EPS {
        TriangleTriangleIntersectionPoint {
            pt: range2[1].1,
            f1: inter_f1,
            f2: range2[1].2,
        }
    } else if range2[1].0 > range1[1].0 + EPS {
        TriangleTriangleIntersectionPoint {
            pt: range1[1].1,
            f1: range1[1].2,
            f2: inter_f2,
        }
    } else {
        TriangleTriangleIntersectionPoint {
            pt: range1[1].1,
            f1: range1[1].2,
            f2: range2[1].2,
        }
    };

    Some(TriangleTriangleIntersection { a, b })
}

fn segment_plane_intersection(
    plane_center: &Point<Real>,
    plane_normal: &Vector<Real>,
    segment: &Segment,
    eid: u32,
    vids: (u32, u32),
) -> Option<(Point<Real>, FeatureId)> {
    let dir = segment.b - segment.a;
    let dir_norm = dir.norm();

    let toi =
        query::details::line_toi_with_halfspace(plane_center, &plane_normal, &segment.a, &dir)?;
    let scaled_toi = toi * dir_norm;

    if scaled_toi < -EPS || scaled_toi > dir_norm + EPS {
        None
    } else {
        if scaled_toi <= EPS {
            Some((segment.a, FeatureId::Vertex(vids.0)))
        } else if scaled_toi >= dir_norm - EPS {
            Some((segment.b, FeatureId::Vertex(vids.1)))
        } else {
            Some((segment.a + dir * toi, FeatureId::Edge(eid)))
        }
    }
}

#[cfg(test)]
mod test {
    use crate::math::{Isometry, Vector};
    use crate::shape::{Cuboid, TriMesh};

    #[test]
    fn cube_cube_trimesh_intersection() {
        let (vtx, idx) = Cuboid::new(Vector::new(1.0, 1.0, 1.0)).to_trimesh();
        let mut trimesh = TriMesh::new(vtx, idx);
        trimesh.compute_topology().unwrap();
        let pos1 = Isometry::identity();
        let pos2 = Isometry::translation(1.0 - 0.51364326, 1.0 - 0.77000046, 0.73529434);
        let _ = super::intersect_meshes(&pos1, &trimesh, &pos2, &trimesh).unwrap();
    }
}
