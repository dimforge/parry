use super::{MeshIntersectionError, TriangleTriangleIntersection};
use crate::math::{Isometry, Real};
use crate::query::{visitors::BoundingVolumeIntersectionsSimultaneousVisitor, PointQuery};
use crate::shape::{GenericTriMesh, TriMesh, Triangle};
use crate::transformation::mesh_intersection::angle_closest_to_90;
use crate::utils::DefaultStorage;
use core::f64::consts::PI;
use na::{Point3, Vector3};
use spade::{ConstrainedDelaunayTriangulation, Triangulation as _};
use std::collections::BTreeMap;
use std::collections::HashSet;
use std::path::PathBuf;
// Dbg
use obj::{Group, IndexTuple, ObjData, Object, SimplePolygon};
use rstar::RTree;

/// Metadata that specifies thresholds to use when making construction choices
/// in mesh intersections.
pub struct MeshIntersectionMetadata {
    /// The smallest angle (in degrees) that will be tolerated. A triangle with
    /// a smaller angle is considered degenerate and will be deleted.
    pub angle_epsilon: f64,
    /// The maximum distance at which two points are considered to overlap in space
    /// if `||p1 - p2|| < global_insertion_epsilon` then p1 and p2 are considered
    /// to be the same point.
    pub global_insertion_epsilon: f64,
    /// A multiplier coefficient to scale `global_insertion_epsilon` when checking for
    /// point duplicatin within a single triangle. Inside of an individual triangle
    /// the distance at wich two points are considered to be the same is
    /// `global_insertion_epsilon * local_insertion_epsilon_mod`.
    pub local_insertion_epsilon_mod: f64,
}

impl Default for MeshIntersectionMetadata {
    fn default() -> Self {
        Self {
            angle_epsilon: 0.005, // degrees
            global_insertion_epsilon: f64::EPSILON * 100.0,
            local_insertion_epsilon_mod: 10.,
        }
    }
}

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
    intersect_meshes_with_metadata(
        pos1,
        mesh1,
        flip1,
        pos2,
        mesh2,
        flip2,
        MeshIntersectionMetadata::default(),
    )
}

/// Similar to `intersect_meshes`.
///
/// It allows to specify epsilons for how the algorithm will behave.
/// See `MeshIntersectionMetadata` for details.
pub fn intersect_meshes_with_metadata(
    pos1: &Isometry<Real>,
    mesh1: &TriMesh,
    flip1: bool,
    pos2: &Isometry<Real>,
    mesh2: &TriMesh,
    flip2: bool,
    meta_data: MeshIntersectionMetadata,
) -> Result<Option<TriMesh>, MeshIntersectionError> {
    // NOTE: remove this, used for debugging only.
    if cfg!(debug_assertions) {
        mesh1.assert_half_edge_topology_is_valid();
        mesh2.assert_half_edge_topology_is_valid();
    }

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

    // 2: Identify all triangles that do actually intersect.
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

    // 3: Grab all triangles that are inside the other mesh but tdo not intersect it.
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

    // 4: Initialize a new mesh by inserting points into a set. Duplicate points should
    // hash to the same index.
    let mut point_set = RTree::<TreePoint, _>::new();
    let mut topology_indices = Vec::new();
    {
        let mut insert_point = |position: Vector3<f64>| {
            insert_into_set(position, &mut point_set, meta_data.global_insertion_epsilon) as u32
        };
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

    // 5: Associate constraint edges generated by a tringle-triangle intersection
    // to each intersecting triangle where they occur.
    let mut constraints1 = BTreeMap::new();
    let mut constraints2 = BTreeMap::new();
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
                    for i in 0..polygon.len() {
                        let a = polygon[i];
                        let b = polygon[(i + 1) % polygon.len()];

                        list1.push([a.p1, b.p1]);
                        list2.push([a.p1, b.p1]);
                    }
                }
            }
        }
    }

    // 6: Collect all triangles that intersect and their associated constraint edges.
    // For each such triangle, compute a CDT of its constraints. For each face in this CDT,
    // if the face is contained in the opposite mesh, add it to the intersection mesh.
    merge_triangle_sets(
        mesh1,
        mesh2,
        &constraints1,
        &pos12,
        flip2,
        &meta_data,
        &mut point_set,
        &mut topology_indices,
    );

    merge_triangle_sets(
        mesh2,
        mesh1,
        &constraints2,
        &Isometry::identity(),
        flip1,
        &meta_data,
        &mut point_set,
        &mut topology_indices,
    );

    // 7: Sort the ouput points by insertion order.
    let mut vertices: Vec<_> = point_set.iter().copied().collect();
    vertices.sort_by(|a, b| a.id.cmp(&b.id));
    let vertices: Vec<_> = vertices.iter().map(|p| Point3::from(p.point)).collect();

    if !topology_indices.is_empty() {
        Ok(Some(TriMesh::new(vertices, topology_indices)))
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

fn syncretize_triangulation(
    tri: &Triangle,
    constraints: &[[Point3<f64>; 2]],
    epsilon: f64,
) -> (
    ConstrainedDelaunayTriangulation<spade::Point2<Real>>,
    Vec<Point3<f64>>,
) {
    let mut constraints = constraints.to_vec();
    // Add the triangle points to the triangulation.
    let mut point_set = RTree::<TreePoint, _>::new();
    let _ = insert_into_set(tri.a.coords, &mut point_set, epsilon);
    let _ = insert_into_set(tri.b.coords, &mut point_set, epsilon);
    let _ = insert_into_set(tri.c.coords, &mut point_set, epsilon);

    // Sometimes, points on the edge of a triangle are slightly off, and this makes
    // spade think that there is a super thin triangle. Project points close to an edge
    // onto the edge to get better results.
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
        let p1_id = insert_into_set(point_pair[0].coords, &mut point_set, epsilon);
        let p2_id = insert_into_set(point_pair[1].coords, &mut point_set, epsilon);

        edges.push([p1_id, p2_id]);
    }

    let mut points: Vec<_> = point_set.iter().cloned().collect();
    points.sort_by(|a, b| a.id.cmp(&b.id));

    let tri_points = [tri.a.coords, tri.b.coords, tri.c.coords];
    let best_source = angle_closest_to_90(&tri_points);
    let d1 = tri_points[best_source] - tri_points[(best_source + 1) % 3];
    let d2 = tri_points[(best_source + 2) % 3] - tri_points[(best_source + 1) % 3];
    let (e1, e2) = planar_gram_schmidt(d1, d2);
    let project = |p: &Vector3<f64>| spade::Point2::new(e1.dot(p), e2.dot(p));

    // Project points into 2D and triangulate the resulting set.

    let planar_points: Vec<_> = points
        .iter()
        .copied()
        .map(|point| {
            let point_proj = project(&point.point);
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

    let points = points.into_iter().map(|p| Point3::from(p.point)).collect();
    (cdt_triangulation, points)
}

// I heavily recommend that this is left here in case one needs to debug the above code.
fn _mesh_to_obj(mesh: &TriMesh, path: &PathBuf) {
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

// I heavily recommend that this is left here in case one needs to debug the above code.
fn _points_to_obj(mesh: &[Point3<f64>], path: &PathBuf) {
    use std::io::Write;
    let mut file = std::fs::File::create(path).unwrap();

    for p in mesh {
        writeln!(file, "v {} {} {}", p.x, p.y, p.z).unwrap();
    }
}

// I heavily recommend that this is left here in case one needs to debug the above code.
fn _points_and_edges_to_obj(mesh: &[Point3<f64>], edges: &[[usize; 2]], path: &PathBuf) {
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
    metadata: &MeshIntersectionMetadata,
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

        let (delaunay, points) = syncretize_triangulation(
            &tri,
            &constraints,
            metadata.global_insertion_epsilon * metadata.local_insertion_epsilon_mod,
        );

        for face in delaunay.inner_faces() {
            let verts = face.vertices();
            let p1 = points[verts[0].index()];
            let p2 = points[verts[1].index()];
            let p3 = points[verts[2].index()];

            // Sometimes the triangulation is messed up due to numerical errors. If
            // a triangle does not survive this test. You can bet it should be put out
            // of its misery.
            if is_triangle_degenerate(
                &[p1.coords, p2.coords, p3.coords],
                metadata.angle_epsilon,
                metadata.global_insertion_epsilon,
            ) {
                continue;
            }

            let center = Triangle {
                a: p1,
                b: p2,
                c: p3,
            }
            .center();

            let epsilon = metadata.global_insertion_epsilon;
            if flip2 ^ (mesh2.contains_local_point(&pos12.inverse_transform_point(&center))) {
                topology_indices.push([
                    insert_into_set(p1.coords, &mut point_set, epsilon) as u32,
                    insert_into_set(p2.coords, &mut point_set, epsilon) as u32,
                    insert_into_set(p3.coords, &mut point_set, epsilon) as u32,
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
