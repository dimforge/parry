// Regression test for https://github.com/dimforge/rapier/issues/969
//
// `to_outline` on round shapes with a zero border radius produced index buffers
// referencing vertices that were never pushed, the arc helper skipping the intermediate
// vertices the index buffer assumed.

use parry3d::math::Vector;
use parry3d::shape::{
    Cone, ConvexPolyhedron, Cuboid, Cylinder, RoundCone, RoundCuboid, RoundCylinder, RoundShape,
};

fn assert_indices_in_bounds(shape_name: &str, vtx: &[Vector], idx: &[[u32; 2]]) {
    for (i, seg) in idx.iter().enumerate() {
        for &id in seg {
            assert!(
                (id as usize) < vtx.len(),
                "{shape_name}: index {id} of segment {i} is out of bounds (only {} vertices)",
                vtx.len()
            );
        }
    }

    for pt in vtx {
        assert!(pt.is_finite(), "{shape_name}: non-finite outline vertex");
    }
}

fn check_all_outlines(border_radius: f32) {
    let round_cuboid = RoundCuboid {
        inner_shape: Cuboid::new(Vector::new(1.0, 1.0, 1.0)),
        border_radius,
    };
    let (vtx, idx) = round_cuboid.to_outline(5);
    assert_indices_in_bounds("RoundCuboid", &vtx, &idx);

    let round_cylinder = RoundCylinder {
        inner_shape: Cylinder::new(1.0, 0.5),
        border_radius,
    };
    let (vtx, idx) = round_cylinder.to_outline(10, 5);
    assert_indices_in_bounds("RoundCylinder", &vtx, &idx);

    let round_cone = RoundCone {
        inner_shape: Cone::new(1.0, 0.5),
        border_radius,
    };
    let (vtx, idx) = round_cone.to_outline(10, 5);
    assert_indices_in_bounds("RoundCone", &vtx, &idx);

    let convex = ConvexPolyhedron::from_convex_hull(&[
        Vector::new(-1.0, -1.0, -1.0),
        Vector::new(1.0, -1.0, -1.0),
        Vector::new(-1.0, 1.0, -1.0),
        Vector::new(1.0, 1.0, -1.0),
        Vector::new(-1.0, -1.0, 1.0),
        Vector::new(1.0, -1.0, 1.0),
        Vector::new(-1.0, 1.0, 1.0),
        Vector::new(1.0, 1.0, 1.0),
    ])
    .unwrap();
    let round_convex = RoundShape {
        inner_shape: convex,
        border_radius,
    };
    let (vtx, idx) = round_convex.to_outline(5);
    assert_indices_in_bounds("RoundConvexPolyhedron", &vtx, &idx);
}

#[test]
fn round_shape_outline_zero_border_radius() {
    check_all_outlines(0.0);
}

#[test]
fn round_shape_outline_tiny_border_radius() {
    check_all_outlines(1.0e-8);
}

#[test]
fn round_shape_outline_regular_border_radius() {
    check_all_outlines(0.1);
}

#[test]
fn round_cuboid_outline_zero_radius_matches_issue_repro() {
    // The exact shape from the issue: ColliderBuilder::round_cuboid(1.0, 1.0, 1.0, 0.0).
    let shape = RoundCuboid {
        inner_shape: Cuboid::new(Vector::new(1.0, 1.0, 1.0)),
        border_radius: 0.0,
    };
    // Rapier’s debug-render pipeline uses nsubdivs = 20.
    let (vtx, idx) = shape.to_outline(20);
    assert_indices_in_bounds("RoundCuboid", &vtx, &idx);
}
