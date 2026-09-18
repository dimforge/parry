//! A `Polyline` flagged `DEFORMABLE` and updated with `set_vertices` must give correct contact
//! manifolds at a fixed relative pose: the cached contact points of the rigid path may not be
//! reused after the vertices moved.

use parry2d::bounding_volume::{Aabb, BoundingVolume};
use parry2d::math::{Pose, Real, Vector};
use parry2d::query::{
    ContactManifold, ContactManifoldsWorkspace, DefaultQueryDispatcher, PersistentQueryDispatcher,
};
use parry2d::shape::{Ball, Cuboid, Polyline, PolylineFlags, Shape};

/// A horizontal line strip at `y = 0` from `-half` to `half` with `n` segments.
fn strip(n: usize, half: Real) -> Vec<Vector> {
    let step = 2.0 * half / n as Real;
    (0..=n)
        .map(|i| Vector::new(-half + i as Real * step, 0.0))
        .collect()
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

    fn update(&mut self, polyline: &Polyline, other: &dyn Shape, pos12: &Pose) {
        DefaultQueryDispatcher
            .contact_manifolds(
                pos12,
                polyline,
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

/// A cuboid rests on the polyline; the polyline then rises into it at the same `pos12`: the
/// segment vs cuboid generator reuses its cached points when the pose did not move, so the rigid
/// path reports the stale penetration and the deformable one the new one.
#[test]
fn cuboid_penetration_follows_the_moved_vertices() {
    let vertices = strip(4, 2.0);
    let cuboid = Cuboid::new(Vector::splat(0.5));
    let pos12 = Pose::translation(0.0, 0.45);

    for deformable in [false, true] {
        let flags = if deformable {
            PolylineFlags::DEFORMABLE
        } else {
            PolylineFlags::empty()
        };
        let mut polyline = Polyline::with_flags(vertices.clone(), None, flags);
        let mut pair = Pair::new();

        pair.update(&polyline, &cuboid, &pos12);
        let dist0 = pair.deepest().unwrap();
        assert!((dist0 + 0.05).abs() < 1.0e-4, "dist0 = {dist0}");

        polyline.set_vertices(&translated(&vertices, Vector::new(0.0, 0.2)));
        pair.update(&polyline, &cuboid, &pos12);
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

/// The ball ends up on the other side of the polyline after the vertices moved past it: the
/// contact normal must flip with the geometry.
#[test]
fn ball_on_the_other_side_after_the_polyline_moved_past_it() {
    let vertices = strip(4, 2.0);
    let ball = Ball::new(0.5);
    let pos12 = Pose::translation(0.0, 0.4);
    let mut polyline = Polyline::with_flags(vertices.clone(), None, PolylineFlags::DEFORMABLE);
    let mut pair = Pair::new();

    pair.update(&polyline, &ball, &pos12);
    let normal0 = pair
        .manifolds
        .iter()
        .find(|m| !m.points.is_empty())
        .unwrap()
        .local_n1;
    assert!(normal0.y > 0.99, "normal0 = {normal0:?}");

    polyline.set_vertices(&translated(&vertices, Vector::new(0.0, 0.8)));
    pair.update(&polyline, &ball, &pos12);
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
    let vertices = strip(40, 10.0);
    let mut polyline = Polyline::new(vertices.clone(), None);

    let deformed: Vec<Vector> = vertices
        .iter()
        .map(|v| Vector::new(v.x, (v.x * 0.7).sin() * 2.0))
        .collect();
    polyline.set_vertices(&deformed);
    assert_eq!(polyline.vertices(), &deformed[..]);

    let rebuilt = Polyline::new(deformed.clone(), None);
    assert!(polyline
        .local_aabb()
        .mins
        .abs_diff_eq(rebuilt.local_aabb().mins, 1.0e-6));
    assert!(polyline
        .local_aabb()
        .maxs
        .abs_diff_eq(rebuilt.local_aabb().maxs, 1.0e-6));

    let queries = [
        Aabb::new(Vector::new(-1.0, -0.5), Vector::new(1.0, 0.5)),
        Aabb::new(Vector::new(2.0, 1.0), Vector::new(6.0, 3.0)),
        Aabb::new(Vector::new(-9.0, 0.5), Vector::new(-6.0, 1.5)),
    ];
    for query in &queries {
        let mut found: Vec<u32> = polyline.bvh().intersect_aabb(query).collect();
        found.sort_unstable();
        let mut expected: Vec<u32> = (0..polyline.num_segments() as u32)
            .filter(|i| polyline.segment(*i).local_aabb().intersects(query))
            .collect();
        expected.sort_unstable();
        assert_eq!(found, expected, "query {query:?}");
        assert!(!expected.is_empty(), "the query should hit something");
    }
}

/// An oriented polyline keeps its one-sided behaviour after `set_vertices` (pseudo-normals are
/// recomputed from the moved vertices).
#[test]
fn oriented_polyline_recomputes_pseudo_normals() {
    // Counter-clockwise square: outward is away from the center.
    let vertices = vec![
        Vector::new(-1.0, -1.0),
        Vector::new(1.0, -1.0),
        Vector::new(1.0, 1.0),
        Vector::new(-1.0, 1.0),
    ];
    let indices = Some(vec![[0, 1], [1, 2], [2, 3], [3, 0]]);
    let mut polyline = Polyline::with_flags(
        vertices,
        indices,
        PolylineFlags::ORIENTED | PolylineFlags::DEFORMABLE,
    );
    let bottom0 = polyline.segment_normal_constraints(0).unwrap().face;
    assert!(bottom0.abs_diff_eq(Vector::new(0.0, -1.0), 1.0e-5));

    // Shear the square: the bottom edge now slopes.
    polyline.set_vertices(&[
        Vector::new(-1.0, -1.0),
        Vector::new(1.0, 0.0),
        Vector::new(1.0, 1.0),
        Vector::new(-1.0, 1.0),
    ]);
    let bottom1 = polyline.segment_normal_constraints(0).unwrap().face;
    assert!(
        bottom1.abs_diff_eq(Vector::new(1.0, -2.0).normalize(), 1.0e-5),
        "bottom1 = {bottom1:?}"
    );
}
