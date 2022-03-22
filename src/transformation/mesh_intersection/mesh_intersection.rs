use super::{MeshIntersectionError, EPS};
use crate::math::{Isometry, Point, Real, Vector};
use crate::query::{self, visitors::BoundingVolumeIntersectionsSimultaneousVisitor, PointQuery};
use crate::shape::{FeatureId, TriMesh, Triangle};
use crate::utils::WBasis;
use na::{point, vector, Point2};
use spade::{handles::FixedVertexHandle, ConstrainedDelaunayTriangulation, Triangulation as _};
use std::collections::{HashMap, HashSet};

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

        if super::triangle_triangle_intersection(&tri1, &tri2).is_some() {
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

    let mut new_vertices1 = vec![];
    let mut new_vertices2 = vec![];

    cut_and_triangulate_intersections(
        &pos12,
        &mesh1,
        &mesh2,
        &mut new_vertices1,
        &mut new_vertices2,
        &mut new_indices1,
        &mut new_indices2,
        &mut intersections,
    );

    let old_vertices1 = mesh1.vertices();
    let old_vertices2 = mesh2.vertices();

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
                let vtx = old_vertices1
                    .get(idx1[k] as usize)
                    .copied()
                    .unwrap_or_else(|| new_vertices1[idx1[k] as usize - old_vertices1.len()]);
                new_vertices.push(vtx);
                new_vertices.len() - 1
            });
            idx1[k] = new_id as u32;
        }
    }

    for idx2 in &mut new_indices2 {
        for k in 0..3 {
            let new_id = *index_map2.entry(idx2[k]).or_insert_with(|| {
                let vtx = old_vertices2
                    .get(idx2[k] as usize)
                    .map(|pt| pos12 * pt)
                    .unwrap_or_else(|| new_vertices2[idx2[k] as usize - old_vertices2.len()]);
                new_vertices.push(vtx);
                new_vertices.len() - 1
            });
            idx2[k] = new_id as u32;
        }
    }

    new_indices1.append(&mut new_indices2);

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

                if mesh2.contains_local_point(&pos12.inverse_transform_point(&tri1.center())) {
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

#[derive(Copy, Clone, Debug)]
struct SpadeInfo {
    point2d: spade::Point2<Real>,
    handle: FixedVertexHandle,
}

struct Triangulation {
    delaunay: ConstrainedDelaunayTriangulation<spade::Point2<Real>>,
    basis: [Vector<Real>; 2],
    vtx_handles: [FixedVertexHandle; 3],
    ref_pt: Point<Real>,
    normal: Vector<Real>,
    ref_proj: [Point2<Real>; 3],
}

impl Triangulation {
    fn new(triangle: Triangle) -> Self {
        let mut delaunay = ConstrainedDelaunayTriangulation::<spade::Point2<Real>>::new();
        let normal = triangle.normal().unwrap();
        let basis = normal.orthonormal_basis();

        let ab = triangle.b - triangle.a;
        let ac = triangle.c - triangle.a;

        let ref_proj = [
            Point2::origin(),
            point![ab.dot(&basis[0]), ab.dot(&basis[1])],
            point![ac.dot(&basis[0]), ac.dot(&basis[1])],
        ];

        let vtx_handles = [
            delaunay
                .insert(spade::Point2::new(ref_proj[0].x, ref_proj[0].y))
                .unwrap(),
            delaunay
                .insert(spade::Point2::new(ref_proj[1].x, ref_proj[1].y))
                .unwrap(),
            delaunay
                .insert(spade::Point2::new(ref_proj[2].x, ref_proj[2].y))
                .unwrap(),
        ];

        Self {
            delaunay,
            basis,
            vtx_handles,
            ref_pt: triangle.a,
            normal: normal.into_inner(),
            ref_proj,
        }
    }

    fn project(&self, pt: Point<Real>, orig_fid: FeatureId) -> spade::Point2<Real> {
        let dpt = pt - self.ref_pt;
        let mut proj = point![dpt.dot(&self.basis[0]), dpt.dot(&self.basis[1])];

        match orig_fid {
            FeatureId::Edge(i) => {
                let a = self.ref_proj[i as usize];
                let b = self.ref_proj[(i as usize + 1) % 3];
                let ab = b - a;
                let ap = proj - a;
                let param = ab.dot(&ap) / ab.norm_squared();
                let shift = vector![ab.y, -ab.x];

                // NOTE: if we have intersections exactly on the edge, we nudge
                //       their projection slightly outside of the triangle. That
                //       way, the triangle’s edge gets split automatically by
                //       the triangulation (or, rather, it will be split when we
                //       add the contsraint involving that point).
                // NOTE: this is not ideal though, so we should find a way to simply
                //       delete spurious triangles that are outside of the intersection
                //       curve.
                proj = a + ab * param + shift * EPS;
            }
            _ => {}
        }

        spade::Point2::new(proj.x, proj.y)
    }
}

fn cut_and_triangulate_intersections(
    pos12: &Isometry<Real>,
    mesh1: &TriMesh,
    mesh2: &TriMesh,
    new_vertices1: &mut Vec<Point<Real>>,
    new_vertices2: &mut Vec<Point<Real>>,
    new_indices1: &mut Vec<[u32; 3]>,
    new_indices2: &mut Vec<[u32; 3]>,
    intersections: &mut Vec<(u32, u32)>,
) {
    let mut triangulations1 = HashMap::new();
    let mut triangulations2 = HashMap::new();
    let mut intersection_points = HashMap::new();

    let mut spade_infos = [HashMap::new(), HashMap::new()];
    let mut spade_handle_to_vertex = [HashMap::new(), HashMap::new()];
    let mut spade_handle_to_intersection = [HashMap::new(), HashMap::new()];

    let new_vertices = [&mut *new_vertices1, &mut *new_vertices2];
    let base_vtx_id = [mesh1.vertices().len(), mesh2.vertices().len()];

    for (i1, i2) in intersections.drain(..) {
        let tris = [mesh1.triangle(i1), mesh2.triangle(i2).transformed(pos12)];
        let vids = [mesh1.indices()[i1 as usize], mesh2.indices()[i2 as usize]];

        if let Some(inter) = super::triangle_triangle_intersection(&tris[0], &tris[1]) {
            let tri_ids = [i1, i2];

            let triangulation1 = triangulations1.entry(tri_ids[0]).or_insert_with(|| {
                let triangulation = Triangulation::new(tris[0]);
                for k in 0..3 {
                    let _ = spade_handle_to_vertex[0].insert(
                        (tri_ids[0], triangulation.vtx_handles[k]),
                        vids[0][k] as usize,
                    );
                }
                triangulation
            });

            let triangulation2 = triangulations2.entry(tri_ids[1]).or_insert_with(|| {
                let triangulation = Triangulation::new(tris[1]);
                for k in 0..3 {
                    let _ = spade_handle_to_vertex[1].insert(
                        (tri_ids[1], triangulation.vtx_handles[k]),
                        vids[1][k] as usize,
                    );
                }
                triangulation
            });

            let triangulations = [triangulation1, triangulation2];

            let fa_1 = convert_fid(mesh1, i1, inter.a.f1);
            let fa_2 = convert_fid(mesh2, i2, inter.a.f2);
            let fb_1 = convert_fid(mesh1, i1, inter.b.f1);
            let fb_2 = convert_fid(mesh2, i2, inter.b.f2);

            let orig_fid_a = [inter.a.f1, inter.a.f2];
            let orig_fid_b = [inter.b.f1, inter.b.f2];
            let key_a = (fa_1, fa_2);
            let key_b = (fb_1, fb_2);

            let mut insert_point = |pt: [_; 2],
                                    key: (FeatureId, FeatureId),
                                    orig_fid: [FeatureId; 2],
                                    i: usize| {
                let spade_key = (tri_ids[i], key);

                spade_infos[i]
                    .entry(spade_key)
                    .or_insert_with(|| {
                        let point2d = triangulations[i].project(pt[i], orig_fid[i]);
                        let handle = triangulations[i].delaunay.insert(point2d).unwrap();
                        let _ = spade_handle_to_vertex[i]
                            .insert((tri_ids[i], handle), base_vtx_id[i] + new_vertices[i].len());
                        let _ = spade_handle_to_intersection[i].insert((tri_ids[i], handle), key);
                        new_vertices[i].push(pt[i]);
                        SpadeInfo { point2d, handle }
                    })
                    .handle
            };

            let ins_a = *intersection_points
                .entry(key_a)
                .or_insert([inter.a.p1, inter.a.p2]);
            let ins_b = *intersection_points
                .entry(key_b)
                .or_insert([inter.b.p1, inter.b.p2]);

            let handles_a = [
                insert_point(ins_a, key_a, orig_fid_a, 0),
                insert_point(ins_a, key_a, orig_fid_a, 1),
            ];
            let handles_b = [
                insert_point(ins_b, key_b, orig_fid_b, 0),
                insert_point(ins_b, key_b, orig_fid_b, 1),
            ];

            for i in 0..2 {
                // NOTE: the naming of the `ConstrainedDelaunayTriangulation::can_add_constraint` method is misleading.
                if !triangulations[i]
                    .delaunay
                    .can_add_constraint(handles_a[i], handles_b[i])
                {
                    let _ = triangulations[i]
                        .delaunay
                        .add_constraint(handles_a[i], handles_b[i]);
                }
            }
        }
    }

    extract_result(
        &pos12,
        &mesh1,
        &mesh2,
        &spade_handle_to_vertex,
        &spade_handle_to_intersection,
        &triangulations1,
        &triangulations2,
        new_vertices1,
        new_vertices2,
        new_indices1,
        new_indices2,
    );
}

fn convert_fid(mesh: &TriMesh, tri: u32, fid: FeatureId) -> FeatureId {
    match fid {
        FeatureId::Edge(eid) => {
            let half_edge_id = mesh.topology().face_half_edges_ids(tri)[eid as usize];
            let half_edge = &mesh.topology().half_edges[half_edge_id as usize];
            FeatureId::Edge(half_edge_id.min(half_edge.twin))
        }
        FeatureId::Vertex(vid) => FeatureId::Vertex(mesh.indices()[tri as usize][vid as usize]),
        FeatureId::Face(_) => FeatureId::Face(tri),
        FeatureId::Unknown => FeatureId::Unknown,
    }
}

fn extract_result(
    pos12: &Isometry<Real>,
    mesh1: &TriMesh,
    mesh2: &TriMesh,
    spade_handle_to_vertex: &[HashMap<(u32, FixedVertexHandle), usize>; 2],
    spade_handle_to_intersection: &[HashMap<(u32, FixedVertexHandle), (FeatureId, FeatureId)>; 2],
    triangulations1: &HashMap<u32, Triangulation>,
    triangulations2: &HashMap<u32, Triangulation>,
    new_vertices1: &mut Vec<Point<Real>>,
    new_vertices2: &mut Vec<Point<Real>>,
    new_indices1: &mut Vec<[u32; 3]>,
    new_indices2: &mut Vec<[u32; 3]>,
) {
    for (tri_id, triangulation) in triangulations1.iter() {
        for face in triangulation.delaunay.inner_faces() {
            let vtx = face.vertices();
            let mut tri = [Point::origin(); 3];
            let mut tri_feat = [(FeatureId::Unknown, FeatureId::Unknown); 3];
            let mut idx = [0; 3];
            for k in 0..3 {
                let vid = spade_handle_to_vertex[0][&(*tri_id, vtx[k].fix())];
                idx[k] = vid as u32;

                if vid < mesh1.vertices().len() {
                    tri[k] = mesh1.vertices()[vid];
                    tri_feat[k] = (FeatureId::Vertex(vid as u32), FeatureId::Unknown);
                } else {
                    tri[k] = new_vertices1[vid - mesh1.vertices().len()];
                }

                if let Some(feat) = spade_handle_to_intersection[0].get(&(*tri_id, vtx[k].fix())) {
                    tri_feat[k] = *feat;
                }
            }

            let tri = Triangle::from(tri);

            if !tri.is_affinely_dependent() && mesh2.contains_point(&pos12, &tri.center()) {
                new_indices1.push(idx);
            }
        }
    }

    for (tri_id, triangulation) in triangulations2.iter() {
        for face in triangulation.delaunay.inner_faces() {
            let vtx = face.vertices();
            let mut tri = [Point::origin(); 3];
            let mut tri_feat = [(FeatureId::Unknown, FeatureId::Unknown); 3];
            let mut idx = [0; 3];
            for k in 0..3 {
                let vid = spade_handle_to_vertex[1][&(*tri_id, vtx[k].fix())];
                idx[k] = vid as u32;

                if vid < mesh2.vertices().len() {
                    tri[k] = pos12 * mesh2.vertices()[vid];
                    tri_feat[k] = (FeatureId::Unknown, FeatureId::Vertex(vid as u32));
                } else {
                    tri[k] = new_vertices2[vid - mesh2.vertices().len()];
                }

                if let Some(feat) = spade_handle_to_intersection[1].get(&(*tri_id, vtx[k].fix())) {
                    tri_feat[k] = *feat;
                }
            }

            let tri = Triangle::from(tri);

            if !tri.is_affinely_dependent() && mesh1.contains_local_point(&tri.center()) {
                new_indices2.push(idx);
            }
        }
    }
}
