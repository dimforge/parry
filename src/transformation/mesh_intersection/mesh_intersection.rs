use super::{MeshIntersectionError, TriangleTriangleIntersection, EPS};
use crate::math::{Isometry, Point, Real, Vector};
use crate::query::{visitors::BoundingVolumeIntersectionsSimultaneousVisitor, PointQuery};
use crate::shape::{FeatureId, GenericTriMesh, TriMesh, Triangle};
use crate::utils::{hashmap, DefaultStorage, WBasis};
use core::f64::consts::PI;
use na::{constraint, ComplexField, Point2, Point3, Vector2, Vector3};
use spade::{handles::FixedVertexHandle, ConstrainedDelaunayTriangulation, Triangulation as _};
use std::collections::BTreeMap;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
// Dbg
use obj::{Group, IndexTuple, ObjData, Object, SimplePolygon};
use rstar::RTree;

const EPSILON: f64 = f64::EPSILON * 100.0;
// const EPSILON: f64 = 0.001;

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
    mesh_to_obj(&mesh1, &PathBuf::from("input1.obj"));
    mesh_to_obj(&mesh2, &PathBuf::from("input2.obj"));

    // NOTE: remove this, used for debugging only.
    mesh1.assert_half_edge_topology_is_valid();
    mesh2.assert_half_edge_topology_is_valid();

    if mesh1.topology().is_none() || mesh2.topology().is_none() {
        return Err(MeshIntersectionError::MissingTopology);
    }

    if mesh1.pseudo_normals().is_none() || mesh2.pseudo_normals().is_none() {
        return Err(MeshIntersectionError::MissingPseudoNormals);
    }

    let pos12 = pos1.inv_mul(pos2);

    // 1: collect all the potential triangle-triangle intersections.
    let mut intersection_candidates = vec![];
    let mut visitor = BoundingVolumeIntersectionsSimultaneousVisitor::with_relative_pos(
        pos12,
        |tri1: &u32, tri2: &u32| {
            intersection_candidates.push((*tri1, *tri2));
            true
        },
    );

    mesh1.qbvh().traverse_bvtt(mesh2.qbvh(), &mut visitor);

    let mut deleted_faces1: HashSet<u32> = HashSet::default();
    let mut deleted_faces2: HashSet<u32> = HashSet::default();
    let mut new_indices1 = vec![];
    let mut new_indices2 = vec![];

    let mut dbg_intersections = vec![];
    let mut intersections = vec![];
    for (fid1, fid2) in &intersection_candidates {
        let tri1 = mesh1.triangle(*fid1);
        let tri2 = mesh2.triangle(*fid2).transformed(&pos12);

        if super::triangle_triangle_intersection(&tri1, &tri2).is_some() {
            intersections.push((*fid1, *fid2));
            let _ = deleted_faces1.insert(*fid1);
            let _ = deleted_faces2.insert(*fid2);

            dbg_intersections.push(tri1.a);
            dbg_intersections.push(tri1.b);
            dbg_intersections.push(tri1.c);

            dbg_intersections.push(tri2.a);
            dbg_intersections.push(tri2.b);
            dbg_intersections.push(tri2.c);
        }
    }

    let n = dbg_intersections.len();
    if !intersections.is_empty() {
        mesh_to_obj(
            &TriMesh::new(
                dbg_intersections,
                (0..n)
                    .step_by(3)
                    .map(|i| [i as u32, (i + 1) as u32, (i + 2) as u32])
                    .collect(),
            ),
            &PathBuf::from(format!("intersections_{}.obj", intersections.len())),
        );
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

    let mut point_set = RTree::<TreePoint, _>::new();
    let mut topology_indices = Vec::new();

    {
        let mut insert_point =
            |position: Vector3<f64>| insert_into_set(position, &mut point_set, EPSILON) as u32;
        // Add the inside vertices and triangles from mesh1
        for mut face in new_indices1 {
            if flip1 {
                face.swap(0, 1);
            }
            topology_indices.push([
                insert_point(mesh1.vertices()[face[0] as usize].coords),
                insert_point(mesh1.vertices()[face[1] as usize].coords),
                insert_point(mesh1.vertices()[face[2] as usize].coords),
            ]);
        }

        // Add the inside vertices and triangles from mesh2
        for mut face in new_indices2 {
            if flip2 {
                face.swap(0, 1);
            }

            topology_indices.push([
                insert_point(mesh2.vertices()[face[0] as usize].coords),
                insert_point(mesh2.vertices()[face[1] as usize].coords),
                insert_point(mesh2.vertices()[face[2] as usize].coords),
            ]);
        }
    }

    let mut dbg_vertices: Vec<_> = point_set.iter().copied().collect();
    dbg_vertices.sort_by(|a, b| a.id.cmp(&b.id));

    mesh_to_obj(
        &TriMesh::new(
            dbg_vertices.iter().map(|p| Point3::from(p.point)).collect(),
            topology_indices.clone(),
        ),
        &PathBuf::from("stage1_merging.obj"),
    );

    // For each intersecting triangle, get their intersection points.
    let mut constraints1 = std::collections::BTreeMap::new();
    let mut constraints2 = std::collections::BTreeMap::new();

    for (fid1, fid2) in &intersections {
        let tri1 = mesh1.triangle(*fid1);
        let tri2 = mesh2.triangle(*fid2).transformed(&pos12);

        let list1 = constraints1.entry(fid1).or_insert(vec![]);
        let list2 = constraints2.entry(fid2).or_insert(vec![]);

        let intersection = super::triangle_triangle_intersection(&tri1, &tri2);
        if intersection.is_some() {
            match intersection.unwrap() {
                TriangleTriangleIntersection::Segment { a, b } => {
                    // For both triangles, add the points in the intersection
                    // and their associated edge to the set.
                    // Note this necessarily introduces duplicate points to the
                    // set that need to be filtered out.
                    list1.push([a.p1, b.p1]);
                    list2.push([a.p1, b.p1]);
                }
                TriangleTriangleIntersection::Polygon(polygon) => {
                    panic!()
                }
            }
        }
    }

    merge_triangle_sets(
        mesh1,
        mesh2,
        &constraints1,
        &pos12,
        flip2,
        &mut point_set,
        &mut topology_indices,
    );

    let dbg_vertices: Vec<_> = point_set.iter().copied().collect();
    let pts: Vec<_> = dbg_vertices.iter().map(|p| Point3::from(p.point)).collect();
    let (_, d) = find_closest_distinct_points(&pts);

    let mut dbg_vertices: Vec<_> = point_set.iter().copied().collect();
    dbg_vertices.sort_by(|a, b| a.id.cmp(&b.id));
    mesh_to_obj(
        &TriMesh::new(
            dbg_vertices.iter().map(|p| Point3::from(p.point)).collect(),
            topology_indices.clone(),
        ),
        &PathBuf::from("stage2_merging.obj"),
    );

    merge_triangle_sets(
        mesh2,
        mesh1,
        &constraints2,
        &Isometry::identity(),
        flip1,
        &mut point_set,
        &mut topology_indices,
    );

    let mut dbg_vertices: Vec<_> = point_set.iter().copied().collect();
    dbg_vertices.sort_by(|a, b| a.id.cmp(&b.id));
    mesh_to_obj(
        &TriMesh::new(
            dbg_vertices.iter().map(|p| Point3::from(p.point)).collect(),
            topology_indices,
        ),
        &PathBuf::from("stage3_merging.obj"),
    );
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

            // Why not use `insert_and_split`?
            // NOTE: if we have intersections exactly on the edge, we nudge
            //       their projection slightly outside of the triangle. That
            //       way, the triangle’s edge gets split automatically by
            //       the triangulation (or, rather, it will be split when we
            //       add the constraint involving that point).
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
    merged_vertices: &mut Vec<Point<Real>>,
    merged_indices: &mut Vec<[u32; 3]>,
    intersections: &mut Vec<(u32, u32)>,
) {
    let mut triangulations1 = HashMap::new();
    let mut triangulations2 = HashMap::new();
    let mut intersection_points = HashMap::new();

    let mut spade_infos = [HashMap::new(), HashMap::new()];
    let mut spade_handle_to_intersection = [HashMap::new(), HashMap::new()];

    let mut dbg_points = Vec::new();
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

                    dbg_points.push(inter_a.p1);
                    dbg_points.push(inter_a.p2);
                    dbg_points.push(inter_b.p1);
                    dbg_points.push(inter_b.p2);

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

    if !dbg_points.is_empty() {
        points_to_obj(&dbg_points, &PathBuf::from("inter_points.obj"));
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
        merged_vertices,
        merged_indices,
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
    merged_vertices: &[Point<Real>],
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
        merged_vertices[(vid - base_id12) as usize]
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
    merged_vertices: &mut Vec<Point<Real>>,
    merged_indices: &mut Vec<[u32; 3]>,
) {
    // Base ids for indexing in the first mesh vertices, second mesh vertices, and new vertices, as if they
    // are part of a single big array.
    let base_id2 = mesh1.vertices().len() as u32;
    let base_id12 = (mesh1.vertices().len() + mesh2.vertices().len()) as u32;

    let mut added_vertices = HashMap::new();
    let mut vertex_remaping = HashMap::new();

    // Generate the new points and setup the mapping between indices from
    // the second mesh, to vertices from the first mesh (for cases of vertex/vertex intersections).
    for (fids, pts) in intersection_points.iter() {
        match *fids {
            (FeatureId::Vertex(vid1), FeatureId::Vertex(vid2)) => {
                let _ = vertex_remaping.insert(vid2, vid1);
            }
            (FeatureId::Vertex(_), _) | (_, FeatureId::Vertex(_)) => {}
            _ => {
                let _ = added_vertices.entry(fids).or_insert_with(|| {
                    merged_vertices.push(pts[0]);
                    merged_vertices.len() as u32 - 1
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

    let mut dbg_triangles = Vec::new();
    for (tri_id, triangulation) in triangulations1.iter() {
        for face in triangulation.delaunay.inner_faces() {
            let vtx = face.vertices();
            let mut tri = [Point::origin(); 3];
            let mut idx = [0; 3];

            let mut dbg_pts = Vec::new();
            for k in 0..3 {
                let fids = spade_handle_to_intersection[0][&(*tri_id, vtx[k].fix())];
                let vid = fids_to_unified_index(fids);
                let vertex = unified_vertex(mesh1, mesh2, merged_vertices, pos12, vid);

                idx[k] = vid;
                tri[k] = vertex;

                dbg_pts.push(vertex);
            }

            let tri = Triangle::from(tri);
            let center = tri.center();
            let projection = mesh2.project_point(pos12, &tri.center(), false);

            if !tri.is_affinely_dependent_eps(EPS * 10.0)
                && ((flip2 ^ projection.is_inside)
                    || (projection.point - center).norm() <= EPS * 10.0)
            {
                dbg_triangles.extend(dbg_pts);
                if flip1 {
                    idx.swap(1, 2);
                }
                merged_indices.push(idx);
            }
        }
    }

    let n = dbg_triangles.len();
    mesh_to_obj(
        &TriMesh::new(
            dbg_triangles,
            (0..n)
                .step_by(3)
                .map(|i| [i as u32, (i + 1) as u32, (i + 2) as u32])
                .collect(),
        ),
        &PathBuf::from(format!("sorted_intersections_{}.obj", n)),
    );

    for (tri_id, triangulation) in triangulations2.iter() {
        for face in triangulation.delaunay.inner_faces() {
            let vtx = face.vertices();
            let mut tri = [Point::origin(); 3];
            let mut idx = [0; 3];
            for k in 0..3 {
                let fids = spade_handle_to_intersection[1][&(*tri_id, vtx[k].fix())];
                let vid = fids_to_unified_index(fids);
                let vertex = unified_vertex(mesh1, mesh2, merged_vertices, pos12, vid);

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
                merged_indices.push(idx);
            }
        }
    }
}

fn test_triangle_crossing(
    tri: &Triangle,
    pos12: &Isometry<Real>,
    mesh: &GenericTriMesh<DefaultStorage>,
) -> usize {
    let projections = [
        mesh.project_point(pos12, &tri.a, false),
        mesh.project_point(pos12, &tri.b, false),
        mesh.project_point(pos12, &tri.c, false),
    ];
    let distances = [
        (projections[0].point - tri.a).norm(),
        (projections[1].point - tri.b).norm(),
        (projections[2].point - tri.c).norm(),
    ];
    let inside_tests = [
        projections[0].is_inside || distances[0] <= EPS * 10.0,
        projections[1].is_inside || distances[1] <= EPS * 10.0,
        projections[2].is_inside || distances[2] <= EPS * 10.0,
    ];
    let inside_count: usize = inside_tests.iter().map(|b| *b as usize).sum();

    inside_count
}

fn syncretize_triangulation(
    tri: &Triangle,
    constraints: &[[Point3<f64>; 2]],
) -> (Triangulation, Vec<Point3<f64>>) {
    let mut constraints = constraints.to_vec();

    let epsilon = EPSILON * 10.0;
    // Add the triangle points to the triangulation.
    let mut point_set = RTree::<TreePoint, _>::new();
    let _ = insert_into_set(tri.a.coords, &mut point_set, epsilon);
    let _ = insert_into_set(tri.b.coords, &mut point_set, epsilon);
    let _ = insert_into_set(tri.c.coords, &mut point_set, epsilon);

    // Sometimes, points on the edge of a triangle are slightly off, and this makes
    // spade think that there is a super thin triangle. Project points close to an edge
    // onto the edge to get better performance.
    let triangle = [tri.a.coords, tri.b.coords, tri.c.coords];
    for point_pair in constraints.iter_mut() {
        let p1 = point_pair[0];
        let p2 = point_pair[1];

        for i in 0..3 {
            let q1 = triangle[i];
            let q2 = triangle[(i + 1) % 3];

            let proj1 = project_point_to_segment(&p1.coords, &[q1, q2]);
            if (p1.coords - proj1).norm() < epsilon {
                point_pair[0] = Point3::from(proj1);
            }

            let proj2 = project_point_to_segment(&p2.coords, &[q1, q2]);
            if (p2.coords - proj2).norm() < epsilon {
                point_pair[1] = Point3::from(proj2);
            }
        }
    }

    // Generate edge, taking care to merge duplicate vertices.
    let mut edges = Vec::new();
    for point_pair in constraints {
        let p1_id = insert_into_set(point_pair[0].coords, &mut point_set, EPSILON);
        let p2_id = insert_into_set(point_pair[1].coords, &mut point_set, EPSILON);

        edges.push([p1_id, p2_id]);
    }

    let mut points: Vec<_> = point_set.iter().cloned().collect();
    points.sort_by(|a, b| a.id.cmp(&b.id));

    let tri_points = [tri.a.coords, tri.b.coords, tri.c.coords];
    let best_source = select_angle_closest_to_90(&tri_points);
    let d1 = tri_points[best_source] - tri_points[(best_source + 1) % 3];
    let d2 = tri_points[(best_source + 2) % 3] - tri_points[(best_source + 1) % 3];
    let (e1, e2) = planar_gram_schmidt(d1, d2);
    let project = |p: &Vector3<f64>| spade::Point2::new(e1.dot(p), e2.dot(p));

    // Project points into 2D and triangulate the resulting set.
    let mut triangulation = Triangulation::new(*tri);

    let planar_points: Vec<_> = points
        .iter()
        .copied()
        .map(|point| {
            let point_proj = project(&point.point); //triangulation.project(Point3::from(point.point), FeatureId::Unknown);
            spade::Point2::new(point_proj.x, point_proj.y)
        })
        .collect();
    let cdt_triangulation =
        ConstrainedDelaunayTriangulation::<spade::Point2<f64>>::bulk_load_cdt_stable(
            planar_points,
            edges,
        )
        .unwrap();
    debug_assert!(cdt_triangulation.vertices().len() == points.len());
    triangulation.delaunay = cdt_triangulation;

    let points = points.into_iter().map(|p| Point3::from(p.point)).collect();
    (triangulation, points)
}

fn mesh_to_obj(mesh: &TriMesh, path: &PathBuf) {
    let mut file = std::fs::File::create(path).unwrap();

    ObjData {
        position: mesh
            .vertices()
            .into_iter()
            .map(|v| [v.x as f32, v.y as f32, v.z as f32])
            .collect(),
        objects: vec![Object {
            groups: vec![Group {
                polys: mesh
                    .indices()
                    .into_iter()
                    .map(|tri| {
                        SimplePolygon(vec![
                            IndexTuple(tri[0] as usize, None, None),
                            IndexTuple(tri[1] as usize, None, None),
                            IndexTuple(tri[2] as usize, None, None),
                        ])
                    })
                    .collect(),
                name: "".to_string(),
                index: 0,
                material: None,
            }],
            name: "".to_string(),
        }],
        ..Default::default()
    }
    .write_to_buf(&mut file)
    .unwrap();
}

fn points_to_obj(mesh: &[Point3<f64>], path: &PathBuf) {
    use std::io::Write;
    let mut file = std::fs::File::create(path).unwrap();

    for p in mesh {
        writeln!(file, "v {} {} {}", p.x, p.y, p.z).unwrap();
    }
}

fn points_and_edges_to_obj(mesh: &[Point3<f64>], edges: &[[usize; 2]], path: &PathBuf) {
    use std::io::Write;
    let mut file = std::fs::File::create(path).unwrap();

    for p in mesh {
        writeln!(file, "v {} {} {}", p.x, p.y, p.z).unwrap();
    }

    for e in edges {
        writeln!(file, "l {} {}", e[0] + 1, e[1] + 1).unwrap();
    }
}

#[derive(Copy, Clone, PartialEq, Debug, Default)]
struct TreePoint {
    point: Vector3<f64>,
    id: usize,
}

impl rstar::Point for TreePoint {
    type Scalar = f64;
    const DIMENSIONS: usize = 3;

    fn generate(mut generator: impl FnMut(usize) -> Self::Scalar) -> Self {
        TreePoint {
            point: Vector3::new(generator(0), generator(1), generator(2)),
            id: usize::MAX,
        }
    }

    fn nth(&self, index: usize) -> Self::Scalar {
        match index {
            0 => self.point.x,
            1 => self.point.y,
            2 => self.point.z,
            _ => unreachable!(),
        }
    }

    fn nth_mut(&mut self, index: usize) -> &mut Self::Scalar {
        match index {
            0 => &mut self.point.x,
            1 => &mut self.point.y,
            2 => &mut self.point.z,
            _ => unreachable!(),
        }
    }
}

fn insert_into_set(
    position: Vector3<f64>,
    point_set: &mut RTree<TreePoint>,
    epsilon: f64,
) -> usize {
    let point_count = point_set.size();
    let point_to_insert = TreePoint {
        point: position,
        id: point_count,
    };

    match point_set.nearest_neighbor(&point_to_insert) {
        Some(tree_point) => {
            if (tree_point.point - position).norm_squared() <= epsilon {
                return tree_point.id;
            } else {
                point_set.insert(point_to_insert);
                debug_assert!(point_set.size() == point_count + 1);
                return point_count;
            }
        }
        None => {
            point_set.insert(point_to_insert);
            debug_assert!(point_set.size() == point_count + 1);
            return point_count;
        }
    }
}

fn spade_to_tri_mesh(delaunay: &ConstrainedDelaunayTriangulation<spade::Point2<Real>>) -> TriMesh {
    let pts = delaunay
        .vertices()
        .map(|v| {
            let p = v.position();
            Point3::<f64>::new(p.x, p.y, 0.0)
        })
        .collect::<Vec<_>>();
    let topology = delaunay
        .inner_faces()
        .map(|f| {
            [
                f.vertices()[0].index() as u32,
                f.vertices()[1].index() as u32,
                f.vertices()[2].index() as u32,
            ]
        })
        .collect();

    TriMesh::new(pts, topology)
}

fn find_closest_distinct_points(points: &[Point3<f64>]) -> ([Point3<f64>; 2], f64) {
    let mut distance = f64::MAX;
    let mut pair_points = [points[0], points[1]];
    for i in 0..points.len() {
        for j in 0..points.len() {
            if i == j {
                continue;
            }

            let d = (points[i].coords - points[j].coords).norm();

            if d < distance {
                distance = d;
                pair_points[0] = points[i];
                pair_points[1] = points[j];
            }
        }
    }

    (pair_points, distance)
}

fn select_angle_closest_to_90(points: &[Vector3<f64>]) -> usize {
    let n = points.len();

    let mut best_cos = 2.0;
    let mut selected_i = 0;
    for i in 0..points.len() {
        let d1 = (points[i] - points[(i + 1) % n]).normalize();
        let d2 = (points[(i + 2) % n] - points[(i + 1) % n]).normalize();

        let cos = d1.dot(&d2);

        if cos.abs() < best_cos {
            best_cos = cos.abs();
            selected_i = i;
        }
    }

    selected_i
}

fn smallest_angle(points: &[Vector3<f64>]) -> f64 {
    let n = points.len();

    let mut worst_cos = 2.0;
    for i in 0..points.len() {
        let d1 = (points[i] - points[(i + 1) % n]).normalize();
        let d2 = (points[(i + 2) % n] - points[(i + 1) % n]).normalize();

        let cos = d1.dot(&d2);

        if cos < worst_cos {
            worst_cos = cos.abs();
        }
    }

    worst_cos.acos() * 180. / PI
}

fn planar_gram_schmidt(v1: Vector3<f64>, v2: Vector3<f64>) -> (Vector3<f64>, Vector3<f64>) {
    let u1 = v1;
    let u2 = v2 - (v2.dot(&u1) / u1.norm_squared()) * u1;

    let e1 = u1.normalize();
    let e2 = u2.normalize();

    (e1, e2)
}

fn project_point_to_segment(point: &Vector3<f64>, segment: &[Vector3<f64>; 2]) -> Vector3<f64> {
    let dir = segment[1] - segment[0];
    let local = point - segment[0];

    let norm = dir.norm();
    // restrict the result to the segment portion of the line.
    let coeff = (dir.dot(&local) / norm).clamp(0., norm);

    segment[0] + coeff * dir.normalize()
}

fn project_to_triangle(point: &Vector3<f64>, tri: &Triangle) -> (Vector3<f64>, f64) {
    let points = [tri.a.coords, tri.b.coords, tri.c.coords];

    let mut selected_point = Vector3::default();
    let mut distance = f64::MAX;
    for i in 0..3 {
        let proj = project_point_to_segment(point, &[points[i], points[(i + 1) % 3]]);
        let d = (proj - point).norm();

        if d < distance {
            distance = d;
            selected_point = proj;
        }
    }

    (selected_point, distance)
}

fn largest_side(tri: &Triangle) -> f64 {
    let mut side = (tri.a.coords - tri.b.coords).norm();
    side = side.max((tri.b.coords - tri.c.coords).norm());
    side = side.max((tri.c.coords - tri.a.coords).norm());
    side
}

fn closest_vertex(point: &Vector3<f64>, tri: &Triangle) -> (Vector3<f64>, f64) {
    let points = [tri.a.coords, tri.b.coords, tri.c.coords];

    let mut selected_point = Vector3::default();
    let mut distance = f64::MAX;
    for i in 0..3 {
        let d = (points[i] - point).norm();

        if d < distance {
            distance = d;
            selected_point = points[i];
        }
    }

    (selected_point, distance)
}

/// No matter how smart we are about computing intersections. It is always possible
/// to create ultra thin triangles when a point lies on an edge of a tirangle. These
/// are degenerate and need to be terminated with extreme prejudice.
fn is_triangle_degenerate(
    triangle: &[Vector3<f64>; 3],
    epsilon_degrees: f64,
    epsilon_distance: f64,
) -> bool {
    if smallest_angle(triangle) < epsilon_degrees {
        return true;
    }

    let mut shortest_side = f64::MAX;
    for i in 0..3 {
        let p1 = triangle[i];
        let p2 = triangle[(i + 1) % 3];

        shortest_side = shortest_side.min((p1 - p2).norm());
    }

    let mut worse_projection_distance = f64::MAX;
    for i in 0..3 {
        let dir = triangle[(i + 1) % 3] - triangle[(i + 2) % 3];
        if dir.norm() < epsilon_distance {
            return true;
        }

        let dir = dir.normalize();
        let proj = (triangle[i] - triangle[(i + 2) % 3]).dot(&dir) * dir + triangle[(i + 2) % 3];

        worse_projection_distance = worse_projection_distance.min((proj - triangle[i]).norm());
    }

    if worse_projection_distance < epsilon_distance {
        return true;
    }

    false
}

fn merge_triangle_sets(
    mesh1: &GenericTriMesh<DefaultStorage>,
    mesh2: &GenericTriMesh<DefaultStorage>,
    triangle_constraints: &BTreeMap<&u32, Vec<[Point3<f64>; 2]>>,
    pos12: &Isometry<Real>,
    flip2: bool,
    mut point_set: &mut RTree<TreePoint>,
    topology_indices: &mut Vec<[u32; 3]>,
) {
    // For each triangle, and each constraint edge associated to that triangle,
    // make a triangulation of the face and sort wether or not each generated
    // sub-triangle is part of the intersection.
    // For each sub-triangle that is part of the intersection, add them to the
    // output mesh.
    for (triangle_id, constraints) in triangle_constraints.iter() {
        let tri = mesh1.triangle(**triangle_id);

        let (triangulation, points) = syncretize_triangulation(&tri, &constraints);

        for face in triangulation.delaunay.inner_faces() {
            let verts = face.vertices();
            let p1 = points[verts[0].index()];
            let p2 = points[verts[1].index()];
            let p3 = points[verts[2].index()];

            // Sometimes the triangulation is messed up due to numerical errors. If
            // a triangle does not survive this test. You can bet it should be put out
            // of its misery.
            if is_triangle_degenerate(&[p1.coords, p2.coords, p3.coords], 0.005, EPSILON) {
                continue;
            }

            let center = Triangle {
                a: p1,
                b: p2,
                c: p3,
            }
            .center();

            if flip2 ^ (mesh2.contains_local_point(&pos12.inverse_transform_point(&center))) {
                topology_indices.push([
                    insert_into_set(p1.coords, &mut point_set, EPSILON) as u32,
                    insert_into_set(p2.coords, &mut point_set, EPSILON) as u32,
                    insert_into_set(p3.coords, &mut point_set, EPSILON) as u32,
                ]);

                if flip2 {
                    topology_indices.last_mut().unwrap().swap(0, 1)
                }

                let id1 = topology_indices.last().unwrap()[0];
                let id2 = topology_indices.last().unwrap()[1];
                let id3 = topology_indices.last().unwrap()[2];

                // If this triggers, yell at Camilo because his algorithm is
                // disfunctional.
                if id1 == id2 || id1 == id3 || id2 == id3 {
                    panic!();
                }
            }
        }
    }
}
