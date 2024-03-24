use super::{MeshIntersectionError, TriangleTriangleIntersection, EPS};
use crate::math::{Isometry, Point, Real, Vector};
use crate::query::{visitors::BoundingVolumeIntersectionsSimultaneousVisitor, PointQuery};
use crate::shape::{FeatureId, TriMesh, Triangle};
use crate::utils::WBasis;
use na::{Point2, Vector2};
use spade::{handles::FixedVertexHandle, ConstrainedDelaunayTriangulation, Triangulation as _};
use std::collections::{HashMap, HashSet};

/// Computes the intersection of two meshes.
///
/// The meshes must be oriented, have their half-edge topology computed, and must not be self-intersecting.
/// The result mesh vertex coordinates are given in the local-space of `mesh1`.
pub fn intersect_meshes(
    pos1: &Isometry<Real>,
    mesh1: &TriMesh,
    flip1: bool,
    pos2: &Isometry<Real>,
    mesh2: &TriMesh,
    flip2: bool,
) -> Result<Option<TriMesh>, MeshIntersectionError> {
    if mesh1.topology().is_none() || mesh2.topology().is_none() {
        return Err(MeshIntersectionError::MissingTopology);
    }

    if mesh1.pseudo_normals().is_none() || mesh2.pseudo_normals().is_none() {
        return Err(MeshIntersectionError::MissingPseudoNormals);
    }

    // NOTE: remove this, used for debugging only.
    mesh1.assert_half_edge_topology_is_valid();
    mesh2.assert_half_edge_topology_is_valid();

    let pos12 = pos1.inv_mul(pos2);

    // 1: collect all the potential triangle-triangle intersections.
    let mut intersections = vec![];
    let mut visitor = BoundingVolumeIntersectionsSimultaneousVisitor::with_relative_pos(
        pos12,
        |tri1: &u32, tri2: &u32| {
            intersections.push((*tri1, *tri2));
            true
        },
    );

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

    extract_connected_components(
        &pos12,
        mesh1,
        mesh2,
        flip2,
        &deleted_faces1,
        &mut new_indices1,
    );
    extract_connected_components(
        &pos12.inverse(),
        mesh2,
        mesh1,
        flip1,
        &deleted_faces2,
        &mut new_indices2,
    );

    let mut new_vertices12 = vec![];
    let mut new_indices12 = vec![];

    cut_and_triangulate_intersections(
        &pos12,
        mesh1,
        flip1,
        mesh2,
        flip2,
        &mut new_vertices12,
        &mut new_indices12,
        &mut intersections,
    );

    let old_vertices1 = mesh1.vertices();
    let old_vertices2 = mesh2.vertices();

    // At this point, we know what triangles we want from the first mesh,
    // and the ones we want from the second mesh. Now we need to build the
    // vertex buffer and adjust the indices accordingly.
    let mut new_vertices = vec![];

    // Maps from unified index to the final vertex index.
    let mut index_map = HashMap::new();
    let base_id2 = mesh1.vertices().len() as u32;

    // Grab all the triangles from the connected component extracted from the first mesh.
    for idx1 in &mut new_indices1 {
        for k in 0..3 {
            let new_id = *index_map.entry(idx1[k]).or_insert_with(|| {
                let vtx = old_vertices1[idx1[k] as usize];
                new_vertices.push(vtx);
                new_vertices.len() - 1
            });
            idx1[k] = new_id as u32;
        }
    }

    // Grab all the triangles from the connected component extracted from the second mesh.
    for idx2 in &mut new_indices2 {
        for k in 0..3 {
            let new_id = *index_map.entry(base_id2 + idx2[k]).or_insert_with(|| {
                let vtx = pos12 * old_vertices2[idx2[k] as usize];
                new_vertices.push(vtx);
                new_vertices.len() - 1
            });
            idx2[k] = new_id as u32;
        }
    }

    // Grab all the trinangles from the intersections.
    for idx12 in &mut new_indices12 {
        for id12 in idx12 {
            let new_id = *index_map.entry(*id12).or_insert_with(|| {
                let vtx = unified_vertex(mesh1, mesh2, &new_vertices12, &pos12, *id12);
                new_vertices.push(vtx);
                new_vertices.len() - 1
            });
            *id12 = new_id as u32;
        }
    }

    if flip1 {
        new_indices1.iter_mut().for_each(|idx| idx.swap(1, 2));
    }

    if flip2 {
        new_indices2.iter_mut().for_each(|idx| idx.swap(1, 2));
    }

    new_indices1.append(&mut new_indices2);
    new_indices1.append(&mut new_indices12);

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
    flip2: bool,
    deleted_faces1: &HashSet<u32>,
    new_indices1: &mut Vec<[u32; 3]>,
) {
    let topo1 = mesh1.topology().unwrap();
    let mut visited: HashSet<u32> = HashSet::default();
    let mut to_visit = vec![];
    let mut visited_conn_comp = if let Some(cc) = mesh1.connected_components() {
        vec![false; cc.ranges.len()] // TODO: use a Vob instead?
    } else {
        vec![]
    };

    for face in deleted_faces1 {
        if let Some(cc) = mesh1.connected_components() {
            visited_conn_comp[cc.face_colors[*face as usize] as usize] = true;
        }

        let eid = topo1.faces[*face as usize].half_edge;
        let edge_a = &topo1.half_edges[eid as usize];
        let edge_b = &topo1.half_edges[edge_a.next as usize];
        let edge_c = &topo1.half_edges[edge_b.next as usize];
        let edges = [edge_a, edge_b, edge_c];

        for edge in edges {
            if let Some(twin) = topo1.half_edges.get(edge.twin as usize) {
                if !deleted_faces1.contains(&twin.face) {
                    let tri1 = mesh1.triangle(twin.face);

                    if flip2
                        ^ mesh2.contains_local_point(&pos12.inverse_transform_point(&tri1.center()))
                    {
                        to_visit.push(twin.face);
                    }
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
            if let Some(twin) = topo1.half_edges.get(edge.twin as usize) {
                if !deleted_faces1.contains(&twin.face) {
                    to_visit.push(twin.face);
                }
            }
        }
    }

    /*
     * Deal with connected components that don’t intersect the other mesh.
     */
    if let Some(cc) = mesh1.connected_components() {
        for (i, range) in cc.ranges.windows(2).enumerate() {
            if !visited_conn_comp[i] {
                // This connected component doesn’t intersect the second mesh.
                // Classify one of its face (the "representative face", can be any
                // face of the connected copmonent) to determine
                // if the whole thing is inside or outside.
                let repr_face = cc.grouped_faces[range[0]];
                let repr_pt = mesh1.triangle(repr_face).center();
                let indices = mesh1.indices();

                if flip2 ^ mesh2.contains_local_point(&pos12.inverse_transform_point(&repr_pt)) {
                    new_indices1.extend(
                        cc.grouped_faces[range[0]..range[1]]
                            .iter()
                            .map(|fid| indices[*fid as usize]),
                    )
                }
            }
        }
    } else if deleted_faces1.is_empty() {
        // Deal with the case where there is no intersection between the meshes.
        let repr_pt = mesh1.triangle(0).center();

        if flip2 ^ mesh2.contains_local_point(&pos12.inverse_transform_point(&repr_pt)) {
            new_indices1.extend_from_slice(mesh1.indices());
        }
    }
}

#[derive(Copy, Clone, Debug)]
struct SpadeInfo {
    handle: FixedVertexHandle,
}

struct Triangulation {
    delaunay: ConstrainedDelaunayTriangulation<spade::Point2<Real>>,
    basis: [Vector<Real>; 2],
    vtx_handles: [FixedVertexHandle; 3],
    ref_pt: Point<Real>,
    ref_proj: [Point2<Real>; 3],
    normalization: na::Matrix2<Real>,
}

impl Triangulation {
    fn new(triangle: Triangle) -> Self {
        let mut delaunay = ConstrainedDelaunayTriangulation::<spade::Point2<Real>>::new();
        let normal = triangle.normal().unwrap();
        let basis = normal.orthonormal_basis();

        let ab = triangle.b - triangle.a;
        let ac = triangle.c - triangle.a;

        let mut ref_proj = [
            Point2::origin(),
            Point2::new(ab.dot(&basis[0]), ab.dot(&basis[1])),
            Point2::new(ac.dot(&basis[0]), ac.dot(&basis[1])),
        ];

        let normalization_inv =
            na::Matrix2::from_columns(&[ref_proj[1].coords, ref_proj[2].coords]);
        let normalization = normalization_inv
            .try_inverse()
            .unwrap_or(na::Matrix2::identity());

        for ref_proj in &mut ref_proj {
            *ref_proj = normalization * *ref_proj;
        }

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
            ref_proj,
            normalization,
        }
    }

    fn project(&self, pt: Point<Real>, orig_fid: FeatureId) -> spade::Point2<Real> {
        let dpt = pt - self.ref_pt;
        let mut proj =
            self.normalization * Point2::new(dpt.dot(&self.basis[0]), dpt.dot(&self.basis[1]));

        if let FeatureId::Edge(i) = orig_fid {
            let a = self.ref_proj[i as usize];
            let b = self.ref_proj[(i as usize + 1) % 3];
            let ab = b - a;
            let ap = proj - a;
            let param = ab.dot(&ap) / ab.norm_squared();
            let shift = Vector2::new(ab.y, -ab.x);

            // NOTE: if we have intersections exactly on the edge, we nudge
            //       their projection slightly outside of the triangle. That
            //       way, the triangle’s edge gets split automatically by
            //       the triangulation (or, rather, it will be split when we
            //       add the contsraint involving that point).
            // NOTE: this is not ideal though, so we should find a way to simply
            //       delete spurious triangles that are outside of the intersection
            //       curve.
            proj = a + ab * param + shift * EPS * 10.0;
        }

        spade::Point2::new(proj.x, proj.y)
    }
}

fn cut_and_triangulate_intersections(
    pos12: &Isometry<Real>,
    mesh1: &TriMesh,
    flip1: bool,
    mesh2: &TriMesh,
    flip2: bool,
    new_vertices12: &mut Vec<Point<Real>>,
    new_indices12: &mut Vec<[u32; 3]>,
    intersections: &mut Vec<(u32, u32)>,
) {
    let mut triangulations1 = HashMap::new();
    let mut triangulations2 = HashMap::new();
    let mut intersection_points = HashMap::new();

    let mut spade_infos = [HashMap::new(), HashMap::new()];
    let mut spade_handle_to_intersection = [HashMap::new(), HashMap::new()];

    for (i1, i2) in intersections.drain(..) {
        let tris = [mesh1.triangle(i1), mesh2.triangle(i2).transformed(pos12)];
        let vids = [mesh1.indices()[i1 as usize], mesh2.indices()[i2 as usize]];

        if let Some(intersection) = super::triangle_triangle_intersection(&tris[0], &tris[1]) {
            let tri_ids = [i1, i2];

            let triangulation1 = triangulations1.entry(tri_ids[0]).or_insert_with(|| {
                let triangulation = Triangulation::new(tris[0]);
                for k in 0..3 {
                    let _ = spade_handle_to_intersection[0].insert(
                        (tri_ids[0], triangulation.vtx_handles[k]),
                        (FeatureId::Vertex(vids[0][k]), FeatureId::Unknown),
                    );
                }
                triangulation
            });

            let triangulation2 = triangulations2.entry(tri_ids[1]).or_insert_with(|| {
                let triangulation = Triangulation::new(tris[1]);
                for k in 0..3 {
                    let _ = spade_handle_to_intersection[1].insert(
                        (tri_ids[1], triangulation.vtx_handles[k]),
                        (FeatureId::Unknown, FeatureId::Vertex(vids[1][k])),
                    );
                }
                triangulation
            });

            let triangulations = [triangulation1, triangulation2];

            let mut insert_point =
                |pt: [_; 2], key: (FeatureId, FeatureId), orig_fid: [FeatureId; 2], i: usize| {
                    let spade_key = (tri_ids[i], key);

                    spade_infos[i]
                        .entry(spade_key)
                        .or_insert_with(|| {
                            let point2d = triangulations[i].project(pt[i], orig_fid[i]);
                            let handle = triangulations[i].delaunay.insert(point2d).unwrap();
                            let _ =
                                spade_handle_to_intersection[i].insert((tri_ids[i], handle), key);
                            SpadeInfo { handle }
                        })
                        .handle
                };

            match intersection {
                TriangleTriangleIntersection::Segment {
                    a: inter_a,
                    b: inter_b,
                } => {
                    let fa_1 = convert_fid(mesh1, i1, inter_a.f1);
                    let fa_2 = convert_fid(mesh2, i2, inter_a.f2);
                    let fb_1 = convert_fid(mesh1, i1, inter_b.f1);
                    let fb_2 = convert_fid(mesh2, i2, inter_b.f2);

                    let orig_fid_a = [inter_a.f1, inter_a.f2];
                    let orig_fid_b = [inter_b.f1, inter_b.f2];
                    let key_a = (fa_1, fa_2);
                    let key_b = (fb_1, fb_2);

                    let ins_a = *intersection_points
                        .entry(key_a)
                        .or_insert([inter_a.p1, inter_a.p2]);
                    let ins_b = *intersection_points
                        .entry(key_b)
                        .or_insert([inter_b.p1, inter_b.p2]);

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
                TriangleTriangleIntersection::Polygon(intersections) => {
                    for inter in intersections {
                        let f1 = convert_fid(mesh1, i1, inter.f1);
                        let f2 = convert_fid(mesh2, i2, inter.f2);
                        let orig_fid = [inter.f1, inter.f2];
                        let key = (f1, f2);
                        let ins = *intersection_points
                            .entry(key)
                            .or_insert([inter.p1, inter.p2]);

                        let _ = insert_point(ins, key, orig_fid, 0);
                        let _ = insert_point(ins, key, orig_fid, 1);
                    }
                }
            }
        }
    }

    extract_result(
        pos12,
        mesh1,
        flip1,
        mesh2,
        flip2,
        &spade_handle_to_intersection,
        &intersection_points,
        &triangulations1,
        &triangulations2,
        new_vertices12,
        new_indices12,
    );
}

fn convert_fid(mesh: &TriMesh, tri: u32, fid: FeatureId) -> FeatureId {
    match fid {
        FeatureId::Edge(eid) => {
            let topology = mesh.topology().unwrap();
            let half_edge_id = topology.face_half_edges_ids(tri)[eid as usize];
            let half_edge = &topology.half_edges[half_edge_id as usize];
            // NOTE: if the twin doesn’t exist, it’s equal to u32::MAX. So the `min` will
            //       automatically filter it out.
            FeatureId::Edge(half_edge_id.min(half_edge.twin))
        }
        FeatureId::Vertex(vid) => FeatureId::Vertex(mesh.indices()[tri as usize][vid as usize]),
        FeatureId::Face(_) => FeatureId::Face(tri),
        FeatureId::Unknown => FeatureId::Unknown,
    }
}

fn unified_vertex(
    mesh1: &TriMesh,
    mesh2: &TriMesh,
    new_vertices12: &[Point<Real>],
    pos12: &Isometry<Real>,
    vid: u32,
) -> Point<Real> {
    let base_id2 = mesh1.vertices().len() as u32;
    let base_id12 = (mesh1.vertices().len() + mesh2.vertices().len()) as u32;

    if vid < base_id2 {
        mesh1.vertices()[vid as usize]
    } else if vid < base_id12 {
        pos12 * mesh2.vertices()[(vid - base_id2) as usize]
    } else {
        new_vertices12[(vid - base_id12) as usize]
    }
}

fn extract_result(
    pos12: &Isometry<Real>,
    mesh1: &TriMesh,
    flip1: bool,
    mesh2: &TriMesh,
    flip2: bool,
    spade_handle_to_intersection: &[HashMap<(u32, FixedVertexHandle), (FeatureId, FeatureId)>; 2],
    intersection_points: &HashMap<(FeatureId, FeatureId), [Point<Real>; 2]>,
    triangulations1: &HashMap<u32, Triangulation>,
    triangulations2: &HashMap<u32, Triangulation>,
    new_vertices12: &mut Vec<Point<Real>>,
    new_indices12: &mut Vec<[u32; 3]>,
) {
    // Base ids for indexing in the first mesh vertices, second mesh vertices, and new vertices, as if they
    // are part of a single big array.
    let base_id2 = mesh1.vertices().len() as u32;
    let base_id12 = (mesh1.vertices().len() + mesh2.vertices().len()) as u32;

    let mut added_vertices = HashMap::new();
    let mut vertex_remaping = HashMap::new();

    // Generate the new points and setup the mapping between indices from
    // the second mesh, to vertices from the first mash (for cases of vertex/vertex intersections).
    for (fids, pts) in intersection_points.iter() {
        match *fids {
            (FeatureId::Vertex(vid1), FeatureId::Vertex(vid2)) => {
                let _ = vertex_remaping.insert(vid2, vid1);
            }
            (FeatureId::Vertex(_), _) | (_, FeatureId::Vertex(_)) => {}
            _ => {
                let _ = added_vertices.entry(fids).or_insert_with(|| {
                    new_vertices12.push(pts[0]);
                    new_vertices12.len() as u32 - 1
                });
            }
        }
    }

    let fids_to_unified_index = |fids| match fids {
        (FeatureId::Vertex(vid1), _) => vid1,
        (_, FeatureId::Vertex(vid2)) => vertex_remaping
            .get(&vid2)
            .copied()
            .unwrap_or(base_id2 + vid2),
        _ => base_id12 + added_vertices[&fids],
    };

    for (tri_id, triangulation) in triangulations1.iter() {
        for face in triangulation.delaunay.inner_faces() {
            let vtx = face.vertices();
            let mut tri = [Point::origin(); 3];
            let mut idx = [0; 3];
            for k in 0..3 {
                let fids = spade_handle_to_intersection[0][&(*tri_id, vtx[k].fix())];
                let vid = fids_to_unified_index(fids);
                let vertex = unified_vertex(mesh1, mesh2, new_vertices12, pos12, vid);

                idx[k] = vid;
                tri[k] = vertex;
            }

            let tri = Triangle::from(tri);
            let center = tri.center();
            let projection = mesh2.project_point(pos12, &tri.center(), false);

            if !tri.is_affinely_dependent_eps(EPS * 10.0)
                && ((flip2 ^ projection.is_inside)
                    || (projection.point - center).norm() <= EPS * 10.0)
            {
                if flip1 {
                    idx.swap(1, 2);
                }
                new_indices12.push(idx);
            }
        }
    }

    for (tri_id, triangulation) in triangulations2.iter() {
        for face in triangulation.delaunay.inner_faces() {
            let vtx = face.vertices();
            let mut tri = [Point::origin(); 3];
            let mut idx = [0; 3];
            for k in 0..3 {
                let fids = spade_handle_to_intersection[1][&(*tri_id, vtx[k].fix())];
                let vid = fids_to_unified_index(fids);
                let vertex = unified_vertex(mesh1, mesh2, new_vertices12, pos12, vid);

                idx[k] = vid;
                tri[k] = vertex;
            }

            let tri = Triangle::from(tri);
            let center = tri.center();

            // TODO: when two faces are coplanar, they will be present in both `triangulation1` and
            //       `triangulation2`. So we need to only pick one of them. Here we already picked the
            //       face from `triangulation1`. So now we need to ignore the duplicate face. Such face
            //       is detected by looking at the distance from the triangle’s center to the other mesh.
            //       If the center lies on the other mesh, then that we have a duplicate face that was already
            //       added in the previous loop.
            let projection = mesh1.project_local_point(&center, false);
            if !tri.is_affinely_dependent_eps(EPS * 10.0)
                && ((flip1 ^ projection.is_inside)
                    && (projection.point - center).norm() > EPS * 10.0)
            {
                if flip2 {
                    idx.swap(1, 2);
                }
                new_indices12.push(idx);
            }
        }
    }
}
