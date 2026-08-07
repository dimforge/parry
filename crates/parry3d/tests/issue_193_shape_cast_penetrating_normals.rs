// Regression test for https://github.com/dimforge/parry/issues/193
// (the four cases from HeartofPhos' 2024-11-29 comment)
//
// A shape slightly penetrating the +Y face of a large cuboid used to get wildly wrong
// shape-cast normals (e.g. [0.707, 0.027, 0.707] instead of ±Y), the GJK-derived normal
// being unreliable for tiny TOIs. Casts with `toi < 1e-4` now fall back on the contact
// query for the impact geometry.

use parry3d::math::{Pose, Vector};
use parry3d::query::{self, ShapeCastOptions};
use parry3d::shape::{Ball, Capsule, Cuboid, Shape};

#[test]
fn slightly_penetrating_casts_report_the_face_normal() {
    let cases: [(Box<dyn Shape>, Pose); 4] = [
        {
            let g = Capsule::new_y(5.0, 1.0);
            let pos = Pose::translation(0.0, g.half_height() - 0.01, 0.0);
            (Box::new(g), pos)
        },
        {
            let g = Capsule::new_y(5.0, 1.0);
            let pos = Pose::translation(0.0, g.half_height() - 0.02, 0.0);
            (Box::new(g), pos)
        },
        {
            let g = Ball::new(1.0);
            let pos = Pose::translation(0.0, g.radius - 0.2, 0.0);
            (Box::new(g), pos)
        },
        {
            let g = Ball::new(1.0);
            let pos = Pose::translation(0.0, g.radius - 0.200001, 0.0);
            (Box::new(g), pos)
        },
    ];

    let g2 = Cuboid::new(Vector::new(20.0, 1.0, 20.0));
    let pos2 = Pose::translation(0.0, -g2.half_extents.y, 0.0);
    let vel2 = Vector::ZERO;

    // Velocity not being straight down played a part in the original failures.
    let vel1 = Vector::new(0.1, 1.0, 0.0);

    let options = ShapeCastOptions {
        compute_impact_geometry_on_penetration: true,
        ..Default::default()
    };

    for (i, (g1, pos1)) in cases.iter().enumerate() {
        let hit = query::cast_shapes(pos1, vel1, &**g1, &pos2, vel2, &g2, options)
            .unwrap()
            .expect("the penetrating cast should return a hit");

        assert!(
            (hit.normal2 - Vector::Y).length() < 1.0e-4,
            "case {i}: normal2 {:?} should be +Y",
            hit.normal2,
        );
        assert!(
            (hit.normal1 + Vector::Y).length() < 1.0e-4,
            "case {i}: normal1 {:?} should be -Y",
            hit.normal1,
        );
    }
}
