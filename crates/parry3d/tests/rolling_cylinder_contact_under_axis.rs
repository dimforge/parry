// A cylinder lying on its side used to get contact points beside the line under its axis, from the
// GJK witness point and from old points kept while it rolled less than a degree, which gave the
// normal force a lever arm so a rolling cylinder sped up by itself.

use parry3d::math::{Pose, Real, Rotation, Vector};
use parry3d::query::{ContactManifold, DefaultQueryDispatcher, PersistentQueryDispatcher};
use parry3d::shape::{Cuboid, Cylinder, RoundCylinder, Shape};

/// Largest lever arm about the shape's axis of a manifold point's normal line, over many poses.
/// With `dphi`, each manifold is first computed at the pose the shape had one step earlier while
/// rolling toward +x by `dphi` radians per step, as a narrow phase does.
fn max_lever_arm(shape: &dyn Shape, radius: Real, pen: Real, dphi: Option<Real>) -> (Real, usize) {
    let ground = Cuboid::new(Vector::new(200.0, 0.5, 200.0));
    let ground_pose = Pose::translation(0.0, -0.5, 0.0);
    let pose = |x: Real, phi: Real| {
        Pose::from_parts(
            Vector::new(x, radius - pen, 0.0),
            Rotation::from_rotation_z(phi)
                * Rotation::from_rotation_x(core::f32::consts::FRAC_PI_2),
        )
    };
    let mut worst: Real = 0.0;
    let mut max_points = 0;
    for i in 0..360 {
        let phi = -(i as Real) * 0.37;
        let x = 3.1 + i as Real * 0.0137;
        let mut manifold = ContactManifold::<(), ()>::new();
        if let Some(dphi) = dphi {
            let prev = ground_pose.inv_mul(&pose(x - radius * dphi, phi + dphi));
            DefaultQueryDispatcher
                .contact_manifold_convex_convex(
                    &prev,
                    &ground,
                    shape,
                    None,
                    None,
                    0.02,
                    &mut manifold,
                )
                .unwrap();
        }
        let cur = pose(x, phi);
        DefaultQueryDispatcher
            .contact_manifold_convex_convex(
                &ground_pose.inv_mul(&cur),
                &ground,
                shape,
                None,
                None,
                0.02,
                &mut manifold,
            )
            .unwrap();
        assert!(!manifold.points.is_empty(), "no contact at pose {i}");
        max_points = max_points.max(manifold.points.len());
        let normal = ground_pose.rotation * manifold.local_n1;
        let axis = cur.rotation * Vector::Y;
        for pt in &manifold.points {
            let arm = (cur.transform_point(pt.local_p2) - cur.translation)
                .cross(normal)
                .dot(axis);
            worst = worst.max(arm.abs());
        }
    }
    (worst, max_points)
}

#[test]
fn cylinder_side_contact_pushes_through_the_axis() {
    let cylinder = Cylinder::new(0.5, 0.5);
    // Penetrating (EPA) and separated within the prediction (GJK closest points).
    for pen in [1.0e-4, -1.0e-3] {
        let (worst, max_points) = max_lever_arm(&cylinder, 0.5, pen, None);
        assert!(
            worst < 1.0e-4,
            "pen {pen}: a normal line passes {worst} m beside the axis"
        );
        assert!(
            max_points >= 2,
            "pen {pen}: the side should touch along its line, got {max_points} point(s)"
        );
    }
}

#[test]
fn rolling_cylinder_contact_pushes_through_the_axis() {
    let cylinder = Cylinder::new(0.5, 0.5);
    let round = RoundCylinder {
        inner_shape: Cylinder::new(0.5, 0.4),
        border_radius: 0.1,
    };
    // 0.48° and 1.9° per step: a 0.5 m cylinder rolling at 1 m/s at 240 Hz and at 60 Hz.
    for dphi in [0.48 as Real, 1.9] {
        let dphi = dphi.to_radians();
        let (worst, _) = max_lever_arm(&cylinder, 0.5, 1.0e-4, Some(dphi));
        assert!(
            worst < 1.0e-4,
            "cylinder, {dphi} rad/step: a normal line passes {worst} m beside the axis"
        );
        let (worst, _) = max_lever_arm(&round, 0.5, -1.0e-3, Some(dphi));
        assert!(
            worst < 1.0e-4,
            "round cylinder, {dphi} rad/step: a normal line passes {worst} m beside the axis"
        );
    }
}
