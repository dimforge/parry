use parry3d::math::{Pose, Vector};
use parry3d::query::details::closest_points_cuboid_cuboid;
use parry3d::query::{self, ClosestPoints};
use parry3d::shape::Cuboid;

// Two disjoint axis-aligned cuboids separated diagonally along all three axes:
// the closest features are one vertex of each box. Analogous to the 2D diagonal
// case; exercises the separating-normal accumulation across every octant.
#[test]
fn closest_points_cuboid_cuboid_axis_aligned_vertex_vertex() {
    let cuboid = Cuboid::new(Vector::new(1.0, 1.0, 1.0));
    // Gaps of (1, 0.5, 0.25) => vertex-vertex distance.
    let true_dist = (1.0f32 + 0.25 + 0.0625).sqrt();

    for sx in [-1.0f32, 1.0] {
        for sy in [-1.0f32, 1.0] {
            for sz in [-1.0f32, 1.0] {
                let pos12 = Pose::translation(sx * 3.0, sy * 2.5, sz * 2.25);
                let corner1 = Vector::new(sx, sy, sz);
                let corner2 = -corner1;

                match closest_points_cuboid_cuboid(&pos12, &cuboid, &cuboid, 2.0) {
                    ClosestPoints::WithinMargin(p1, p2) => {
                        assert_relative_eq!(p1, corner1, epsilon = 1e-5);
                        assert_relative_eq!(p2, corner2, epsilon = 1e-5);
                    }
                    other => panic!(
                        "octant ({sx}, {sy}, {sz}): expected WithinMargin at \
                         distance {true_dist}, got {other:?}"
                    ),
                }

                let dist = query::distance(&Pose::IDENTITY, &cuboid, &pos12, &cuboid)
                    .unwrap()
                    .distance;
                assert_relative_eq!(dist, true_dist, epsilon = 1e-5);
            }
        }
    }
}

// Diagonal separation along x and y with overlapping extents along z: the closest
// features are two parallel z-aligned edges, so any z inside the overlap realizes
// the same distance. The SAT normal correctly has a zero z component, but the
// support point along that degenerate direction picks an arbitrary *end* of
// cuboid2's edge; when that end lies outside cuboid1's z range, the projection
// inflates the reported distance by the z mismatch.
#[test]
fn closest_points_cuboid_cuboid_axis_aligned_edge_edge() {
    let c1 = Cuboid::new(Vector::new(1.0, 1.0, 1.0));
    let c2 = Cuboid::new(Vector::new(1.0, 1.0, 5.0));
    let pos12 = Pose::translation(3.0, 2.5, 0.0);
    // Gaps of (1, 0.5) in x/y, overlap in z => edge-edge distance.
    let true_dist = (1.0f32 + 0.25).sqrt();

    match closest_points_cuboid_cuboid(&pos12, &c1, &c2, 2.0) {
        ClosestPoints::WithinMargin(p1, p2) => {
            let p2_1 = pos12 * p2;
            assert_relative_eq!((p1 - p2_1).length(), true_dist, epsilon = 1e-5);
            // The realized pair must lie on the facing edges: x/y at the corners,
            // matching z anywhere in the overlap.
            assert_relative_eq!(p1.x, 1.0, epsilon = 1e-5);
            assert_relative_eq!(p1.y, 1.0, epsilon = 1e-5);
            assert_relative_eq!(p2_1.x, 2.0, epsilon = 1e-5);
            assert_relative_eq!(p2_1.y, 1.5, epsilon = 1e-5);
            assert_relative_eq!(p1.z, p2_1.z, epsilon = 1e-5);
        }
        other => panic!("expected WithinMargin at distance {true_dist}, got {other:?}"),
    }

    let dist = query::distance(&Pose::IDENTITY, &c1, &pos12, &c2)
        .unwrap()
        .distance;
    assert_relative_eq!(dist, true_dist, epsilon = 1e-5);
}

// 3D analogue of the corner-touching boundary case: exact zero gap along z,
// small positive gaps along x and y.
#[test]
fn closest_points_cuboid_cuboid_axis_aligned_edge_touching() {
    let cuboid = Cuboid::new(Vector::new(1.0, 1.0, 1.0));
    let (gx, gy) = (1.0e-4f32, 0.5e-4f32);
    let true_dist = (gx * gx + gy * gy).sqrt();
    let pos12 = Pose::translation(2.0 + gx, 2.0 + gy, 2.0);

    match closest_points_cuboid_cuboid(&pos12, &cuboid, &cuboid, 0.5) {
        ClosestPoints::WithinMargin(p1, p2) => {
            assert_relative_eq!(p1, Vector::new(1.0, 1.0, 1.0), epsilon = 1e-5);
            assert_relative_eq!(p2, Vector::new(-1.0, -1.0, -1.0), epsilon = 1e-5);
        }
        other => panic!("expected WithinMargin at distance {true_dist}, got {other:?}"),
    }

    let dist = query::distance(&Pose::IDENTITY, &cuboid, &pos12, &cuboid)
        .unwrap()
        .distance;
    assert_relative_eq!(dist, true_dist, epsilon = 1e-4);
}
