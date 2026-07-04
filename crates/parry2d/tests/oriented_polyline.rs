//! End-to-end checks that an oriented polyline's [`PolylineFlags`] reach the contact-manifold path
//! (the same path the physics engine uses), not just the standalone pseudo-normal math.

use parry2d::math::{Pose, Real, Rotation, Vector};
use parry2d::query::{
    ContactManifold, ContactManifoldsWorkspace, DefaultQueryDispatcher, PersistentQueryDispatcher,
};
use parry2d::shape::{Ball, Polyline, PolylineFlags};

/// Deepest contact normal (on the polyline, pointing toward the ball) of `ball` at `ball_pos`, or
/// `None` if the manifold is empty.
fn deepest_normal(polyline: &Polyline, ball: &Ball, ball_pos: Vector) -> Option<Vector> {
    let pos12 = Pose::from_parts(ball_pos, Rotation::identity());
    let mut manifolds: Vec<ContactManifold<(), ()>> = Vec::new();
    let mut workspace: Option<ContactManifoldsWorkspace> = None;
    DefaultQueryDispatcher
        .contact_manifolds(&pos12, polyline, ball, 0.0, &mut manifolds, &mut workspace)
        .unwrap();

    let mut deepest: Option<(Vector, Real)> = None;
    for manifold in &manifolds {
        for contact in &manifold.points {
            if deepest.map_or(true, |(_, best)| contact.dist < best) {
                deepest = Some((manifold.local_n1, contact.dist));
            }
        }
    }
    deepest.map(|(normal, _)| normal)
}

#[test]
fn oriented_polyline_keeps_front_contact_outward() {
    // Ball in front (play side): orientation keeps the contact, normal outward (+Y).
    let vertices = vec![Vector::new(1.0, 0.0), Vector::new(-1.0, 0.0)];
    let indices = Some(vec![[0, 1]]);
    let ball = Ball::new(0.7);
    let ball_pos = Vector::new(0.0, 0.5);

    let oriented = Polyline::with_flags(vertices, indices, PolylineFlags::ORIENTED);
    let normal = deepest_normal(&oriented, &ball, ball_pos).expect("front contact kept");

    assert!(normal.dot(Vector::Y) > 0.0);
}

#[test]
fn oriented_polyline_is_one_sided() {
    // Solid below (outward up); a ball behind it collides when double-sided but is ignored once oriented.
    let vertices = vec![Vector::new(1.0, 0.0), Vector::new(-1.0, 0.0)];
    let indices = Some(vec![[0, 1]]);
    let ball = Ball::new(1.0);
    let ball_pos = Vector::new(0.0, -0.5);

    let plain = Polyline::new(vertices.clone(), indices.clone());
    let oriented = Polyline::with_flags(vertices, indices, PolylineFlags::ORIENTED);

    assert!(deepest_normal(&plain, &ball, ball_pos).is_some());
    assert!(deepest_normal(&oriented, &ball, ball_pos).is_none());
}
