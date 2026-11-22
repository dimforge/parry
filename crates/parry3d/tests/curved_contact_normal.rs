// A cylinder lying on a box, barely apart from it or barely into it, sometimes got a contact
// normal tens of degrees off the box's face from GJK or EPA, so only one contact point was kept
// and it pushed the cylinder sideways.

use parry3d::math::{Pose, Real, Rotation, Vector};
use parry3d::query::{ContactManifold, DefaultQueryDispatcher, PersistentQueryDispatcher};
use parry3d::shape::{Cuboid, Cylinder, RoundCylinder, Shape};

fn worst_normal_angle(shape: &dyn Shape, radius: Real) -> Real {
    let ground = Cuboid::new(Vector::new(1000.0, 0.5, 1000.0));
    let ground_pose = Pose::translation(0.0, -0.5, 0.0);
    let mut worst: Real = 0.0;
    for gap in [-5.0e-4, -2.0e-4, -5.0e-5, 0.0, 3.0e-5, 1.0e-4] {
        let mut manifold = ContactManifold::<(), ()>::new();
        // Rolling toward +x by 2° a step, each manifold computed from the previous one.
        for i in 0..2000 {
            let phi = -(i as Real) * 0.035;
            let pose = Pose::from_parts(
                Vector::new(i as Real * 0.035 * radius, radius + gap, 0.0),
                Rotation::from_rotation_z(phi)
                    * Rotation::from_rotation_x(core::f32::consts::FRAC_PI_2),
            );
            DefaultQueryDispatcher
                .contact_manifold_convex_convex(
                    &ground_pose.inv_mul(&pose),
                    &ground,
                    shape,
                    None,
                    None,
                    0.02,
                    &mut manifold,
                )
                .unwrap();
            let up = manifold.local_n1.dot(Vector::Y).clamp(-1.0, 1.0);
            worst = worst.max(up.acos().to_degrees());
        }
    }
    worst
}

#[test]
fn rolling_cylinder_normal_stays_on_the_face_normal() {
    let worst = worst_normal_angle(&Cylinder::new(0.5, 0.5), 0.5);
    assert!(worst < 0.1, "cylinder normal {worst}° off the face normal");
    let worst = worst_normal_angle(
        &RoundCylinder {
            inner_shape: Cylinder::new(0.5, 0.4),
            border_radius: 0.1,
        },
        0.5,
    );
    assert!(
        worst < 0.1,
        "round cylinder normal {worst}° off the face normal"
    );
}
