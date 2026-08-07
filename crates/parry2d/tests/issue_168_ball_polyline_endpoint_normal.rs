// Regression test for https://github.com/dimforge/rapier/issues/168
//
// A ball contacting a polyline exactly at (or beyond) its endpoint vertex used
// to produce a zero contact normal. The ball-vs-convex manifold computation now
// has a degenerate fallback so the normal is always finite and non-zero.

use parry2d::math::{Pose, Real, Vector};
use parry2d::query::{ContactManifold, DefaultQueryDispatcher, PersistentQueryDispatcher};
use parry2d::shape::{Ball, Polyline};

fn check_manifold_normals(pos12: Pose, prediction: Real) -> usize {
    let polyline = Polyline::new(vec![Vector::new(-10.0, 0.0), Vector::new(10.0, 0.0)], None);
    let ball = Ball::new(1.0);

    let dispatcher = DefaultQueryDispatcher;
    let mut manifolds: Vec<ContactManifold<(), ()>> = Vec::new();
    let mut workspace = None;

    dispatcher
        .contact_manifolds(
            &pos12,
            &polyline,
            &ball,
            prediction,
            &mut manifolds,
            &mut workspace,
        )
        .unwrap();

    let mut num_contacts = 0;

    for manifold in &manifolds {
        if manifold.points.is_empty() {
            continue;
        }

        num_contacts += manifold.points.len();

        for local_n in [manifold.local_n1, manifold.local_n2] {
            assert!(
                local_n.is_finite(),
                "non-finite contact normal: {local_n:?} (pos12: {pos12:?})"
            );
            assert!(
                local_n.length() > 0.9,
                "degenerate contact normal: {local_n:?} (pos12: {pos12:?})"
            );
        }
    }

    num_contacts
}

#[test]
fn ball_beyond_polyline_endpoint() {
    // The geometry from the issue: the polyline ends at x = 10.0, the ball
    // (radius 1) is at x = 11.0, y = 1.0, so it barely misses the endpoint but
    // is within the contact prediction distance.
    let num_contacts = check_manifold_normals(Pose::translation(11.0, 1.0), 1.0);
    assert!(num_contacts > 0);
}

#[test]
fn ball_aligned_with_polyline_beyond_endpoint() {
    // The `______o` case from the issue: ball right beside the line, aligned
    // with it but not touching it.
    let num_contacts = check_manifold_normals(Pose::translation(11.5, 0.0), 1.0);
    assert!(num_contacts > 0);
}

#[test]
fn ball_centered_on_polyline_endpoint() {
    // Fully degenerate case: the ball center lies exactly on the endpoint
    // vertex, so the projection distance is zero.
    let num_contacts = check_manifold_normals(Pose::translation(10.0, 0.0), 0.1);
    assert!(num_contacts > 0);
}

#[test]
fn ball_overlapping_polyline_endpoint() {
    // Overlap case without prediction: the ball overlaps the endpoint vertex.
    let num_contacts = check_manifold_normals(Pose::translation(10.5, 0.5), 0.0);
    assert!(num_contacts > 0);
}
