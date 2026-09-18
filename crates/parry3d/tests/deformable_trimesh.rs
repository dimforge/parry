//! A `TriMesh` flagged `DEFORMABLE` and updated with `set_vertices` must give correct contact
//! manifolds at a fixed relative pose: neither the interfering-triangle cache nor the cached
//! contact points of the rigid path may be reused.

use parry3d::bounding_volume::{Aabb, BoundingVolume};
use parry3d::math::{Pose, Real, Vector};
use parry3d::query::{
    ContactManifold, ContactManifoldsWorkspace, DefaultQueryDispatcher, PersistentQueryDispatcher,
};
use parry3d::shape::{Ball, Cuboid, Shape, TriMesh, TriMeshFlags};

/// A flat `n x n` quad grid at `y = 0` spanning `[-half, half]` on x and z.
fn grid(n: usize, half: Real) -> (Vec<Vector>, Vec<[u32; 3]>) {
    let mut vertices = Vec::new();
    let mut indices = Vec::new();
    let step = 2.0 * half / n as Real;
    for i in 0..=n {
        for j in 0..=n {
            vertices.push(Vector::new(
                -half + i as Real * step,
                0.0,
                -half + j as Real * step,
            ));
        }
    }
    let row = (n + 1) as u32;
    for i in 0..n as u32 {
        for j in 0..n as u32 {
            let a = i * row + j;
            indices.push([a, a + 1, a + row + 1]);
            indices.push([a, a + row + 1, a + row]);
        }
    }
    (vertices, indices)
}

struct Pair {
    manifolds: Vec<ContactManifold<(), ()>>,
    workspace: Option<ContactManifoldsWorkspace>,
}

impl Pair {
    fn new() -> Self {
        Self {
            manifolds: Vec::new(),
            workspace: None,
        }
    }

    fn update(&mut self, mesh: &TriMesh, other: &dyn Shape, pos12: &Pose) {
        DefaultQueryDispatcher
            .contact_manifolds(
                pos12,
                mesh,
                other,
                0.1,
                &mut self.manifolds,
                &mut self.workspace,
            )
            .unwrap();
    }

    fn deepest(&self) -> Option<Real> {
        self.manifolds
            .iter()
            .flat_map(|m| m.points.iter().map(|p| p.dist))
            .fold(None, |acc: Option<Real>, d| {
                Some(acc.map_or(d, |a| a.min(d)))
            })
    }
}

fn translated(vertices: &[Vector], shift: Vector) -> Vec<Vector> {
    vertices.iter().map(|v| *v + shift).collect()
}

/// The ball starts above the middle of the mesh; the mesh then slides sideways so different
/// triangles end up under the ball (same `pos12`): the rigid path keeps its cached interfering
/// triangles (now far away) and finds nothing, the deformable one finds the new triangles.
#[test]
fn ball_finds_the_triangles_moved_under_it() {
    let (vertices, indices) = grid(8, 4.0);
    let ball = Ball::new(0.5);
    let pos12 = Pose::translation(0.0, 0.4, 0.0);

    for deformable in [false, true] {
        let flags = if deformable {
            TriMeshFlags::DEFORMABLE
        } else {
            TriMeshFlags::empty()
        };
        let mut mesh = TriMesh::with_flags(vertices.clone(), indices.clone(), flags).unwrap();
        let mut pair = Pair::new();

        pair.update(&mesh, &ball, &pos12);
        assert!(pair.deepest().unwrap() < 0.0, "initial contact expected");

        // The vertices under the ball move 3 m along x: other triangles are now under it.
        mesh.set_vertices(&translated(&vertices, Vector::new(3.0, 0.0, 0.0)));
        pair.update(&mesh, &ball, &pos12);
        let deepest = pair.deepest();

        if deformable {
            let dist = deepest.expect("the deformable mesh must find the new triangles");
            assert!((dist + 0.1).abs() < 1.0e-4, "dist = {dist}");
        } else {
            assert!(
                deepest.is_none(),
                "the rigid path is expected to keep its stale interfering triangles"
            );
        }
    }
}

/// A cuboid rests on the mesh; the mesh then rises into it at the same `pos12`: the cuboid vs
/// triangle generator reuses its cached points when the pose did not move, so the rigid path
/// reports the stale penetration and the deformable one the new one.
#[test]
fn cuboid_penetration_follows_the_moved_vertices() {
    let (vertices, indices) = grid(4, 2.0);
    let cuboid = Cuboid::new(Vector::splat(0.5));
    let pos12 = Pose::translation(0.0, 0.45, 0.0);

    for deformable in [false, true] {
        let flags = if deformable {
            TriMeshFlags::DEFORMABLE
        } else {
            TriMeshFlags::empty()
        };
        let mut mesh = TriMesh::with_flags(vertices.clone(), indices.clone(), flags).unwrap();
        let mut pair = Pair::new();

        pair.update(&mesh, &cuboid, &pos12);
        let dist0 = pair.deepest().unwrap();
        assert!((dist0 + 0.05).abs() < 1.0e-4, "dist0 = {dist0}");

        mesh.set_vertices(&translated(&vertices, Vector::new(0.0, 0.2, 0.0)));
        pair.update(&mesh, &cuboid, &pos12);
        let dist1 = pair.deepest().unwrap();

        if deformable {
            assert!((dist1 + 0.25).abs() < 1.0e-4, "dist1 = {dist1}");
        } else {
            assert!(
                (dist1 - dist0).abs() < 1.0e-6,
                "the rigid path is expected to keep its stale points (dist1 = {dist1})"
            );
        }
    }
}

/// The ball ends up on the other side of the mesh after the vertices moved past it: the contact
/// normal must flip with the geometry.
#[test]
fn ball_on_the_other_side_after_the_mesh_moved_past_it() {
    let (vertices, indices) = grid(4, 2.0);
    let ball = Ball::new(0.5);
    let pos12 = Pose::translation(0.0, 0.4, 0.0);
    let mut mesh =
        TriMesh::with_flags(vertices.clone(), indices.clone(), TriMeshFlags::DEFORMABLE).unwrap();
    let mut pair = Pair::new();

    pair.update(&mesh, &ball, &pos12);
    let normal0 = pair
        .manifolds
        .iter()
        .find(|m| !m.points.is_empty())
        .unwrap()
        .local_n1;
    assert!(normal0.y > 0.99, "normal0 = {normal0:?}");

    // The mesh moves up through the ball: the ball is now below it.
    mesh.set_vertices(&translated(&vertices, Vector::new(0.0, 0.8, 0.0)));
    pair.update(&mesh, &ball, &pos12);
    let normal1 = pair
        .manifolds
        .iter()
        .find(|m| !m.points.is_empty())
        .unwrap()
        .local_n1;
    assert!(normal1.y < -0.99, "normal1 = {normal1:?}");
    let dist = pair.deepest().unwrap();
    assert!((dist + 0.1).abs() < 1.0e-4, "dist = {dist}");
}

/// After `set_vertices`, the refitted BVH must answer AABB queries exactly like a rebuilt one.
#[test]
fn bvh_refit_after_set_vertices_matches_brute_force() {
    let (vertices, indices) = grid(10, 5.0);
    let mut mesh = TriMesh::new(vertices.clone(), indices.clone()).unwrap();

    // A smooth deformation: a wave along x plus a shear along z.
    let deformed: Vec<Vector> = vertices
        .iter()
        .map(|v| Vector::new(v.x + 0.3 * v.z, (v.x * 1.3).sin() * 1.5, v.z))
        .collect();
    mesh.set_vertices(&deformed);
    assert_eq!(mesh.vertices(), &deformed[..]);

    let rebuilt = TriMesh::new(deformed.clone(), indices.clone()).unwrap();
    assert!(mesh
        .local_aabb()
        .mins
        .abs_diff_eq(rebuilt.local_aabb().mins, 1.0e-6));
    assert!(mesh
        .local_aabb()
        .maxs
        .abs_diff_eq(rebuilt.local_aabb().maxs, 1.0e-6));

    let queries = [
        Aabb::new(Vector::new(-1.0, -0.5, -1.0), Vector::new(1.0, 0.5, 1.0)),
        Aabb::new(Vector::new(2.0, 0.5, -6.0), Vector::new(6.0, 2.0, 6.0)),
        Aabb::new(Vector::new(-6.0, -2.0, 3.0), Vector::new(-3.0, 0.0, 6.0)),
        Aabb::new(Vector::new(-0.2, 1.4, -0.2), Vector::new(0.2, 1.6, 0.2)),
    ];
    for query in &queries {
        let mut found: Vec<u32> = mesh.bvh().intersect_aabb(query).collect();
        found.sort_unstable();
        let mut expected: Vec<u32> = (0..mesh.num_triangles() as u32)
            .filter(|i| mesh.triangle(*i).local_aabb().intersects(query))
            .collect();
        expected.sort_unstable();
        assert_eq!(found, expected, "query {query:?}");
        assert!(!expected.is_empty(), "the query should hit something");
    }
}

/// `update_vertices` is the in-place variant of `set_vertices`.
#[test]
fn update_vertices_refits_too() {
    let (vertices, indices) = grid(3, 1.5);
    let mut mesh = TriMesh::new(vertices, indices).unwrap();
    mesh.update_vertices(|vtx| {
        for v in vtx {
            v.y = 2.0;
        }
    });
    let aabb = mesh.local_aabb();
    assert!((aabb.mins.y - 2.0).abs() < 1.0e-6 && (aabb.maxs.y - 2.0).abs() < 1.0e-6);
    let hits: Vec<u32> = mesh
        .bvh()
        .intersect_aabb(&Aabb::new(
            Vector::new(-2.0, 1.9, -2.0),
            Vector::new(2.0, 2.1, 2.0),
        ))
        .collect();
    assert_eq!(hits.len(), mesh.num_triangles());
}
