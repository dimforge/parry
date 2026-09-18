use parry2d::math::{Pose, Vector};
use parry2d::query::details::closest_points_cuboid_cuboid;
use parry2d::query::{self, ClosestPoints};
use parry2d::shape::Cuboid;

// Two disjoint axis-aligned cuboids separated diagonally: the closest features are
// one corner of each box. The SAT separating normal is the face normal of the axis
// with the largest gap (here ±X), and the support point of the other cuboid in that
// direction is degenerate: every point of the facing edge is a valid support point.
// The corner actually returned is decided by the sign convention of the zero
// component, not by which corner is closest, so depending on the quadrant the
// projected point can be the corner *away* from the first cuboid, overestimating
// the distance or reporting `Disjoint` within the margin.
#[test]
fn closest_points_cuboid_cuboid_axis_aligned_diagonal() {
    let cuboid = Cuboid::new(Vector::new(1.0, 1.0));
    // Gap of 1.0 along x and 0.5 along y => corner-corner distance.
    let true_dist = (1.0f32 * 1.0 + 0.5 * 0.5).sqrt();

    for (x, y) in [(3.0, 2.5), (3.0, -2.5), (-3.0, 2.5), (-3.0, -2.5)] {
        let pos12 = Pose::translation(x, y);
        // The closest corner of cuboid1 points toward cuboid2, and vice-versa.
        let corner1 = Vector::new(x.signum(), y.signum());
        let corner2 = -corner1;

        match closest_points_cuboid_cuboid(&pos12, &cuboid, &cuboid, 2.0) {
            ClosestPoints::WithinMargin(p1, p2) => {
                assert_relative_eq!(p1, corner1, epsilon = 1e-5);
                assert_relative_eq!(p2, corner2, epsilon = 1e-5);
            }
            other => panic!(
                "pos12 = ({x}, {y}): expected WithinMargin at distance {true_dist}, got {other:?}"
            ),
        }

        // The same failure is reachable through the public API: `query::distance`
        // dispatches cuboid-cuboid pairs to `distance_cuboid_cuboid`, which is
        // implemented on top of `closest_points_cuboid_cuboid`.
        let dist = query::distance(&Pose::IDENTITY, &cuboid, &pos12, &cuboid)
            .unwrap()
            .distance;
        assert_relative_eq!(dist, true_dist, epsilon = 1e-5);
    }
}

// Boundary case of the diagonal configuration: the gap along one axis is exactly
// zero (the boxes are corner-touching along y) while the other axis has a small
// positive gap. Only one axis has a strictly positive separation, so SAT reports
// an axis-aligned normal whose other component is zero, and the support point in
// that degenerate direction again picks an arbitrary corner of the facing edge.
#[test]
fn closest_points_cuboid_cuboid_axis_aligned_corner_touching() {
    let cuboid = Cuboid::new(Vector::new(1.0, 1.0));
    let gap = 1.0e-4;
    let pos12 = Pose::translation(2.0 + gap, 2.0);

    match closest_points_cuboid_cuboid(&pos12, &cuboid, &cuboid, 0.5) {
        ClosestPoints::WithinMargin(p1, p2) => {
            assert_relative_eq!(p1, Vector::new(1.0, 1.0), epsilon = 1e-5);
            assert_relative_eq!(p2, Vector::new(-1.0, -1.0), epsilon = 1e-5);
        }
        other => panic!("expected WithinMargin at distance {gap}, got {other:?}"),
    }

    let dist = query::distance(&Pose::IDENTITY, &cuboid, &pos12, &cuboid)
        .unwrap()
        .distance;
    assert_relative_eq!(dist, gap, epsilon = 1e-5);
}

// Fuzzer-found rounding case: the true gap along x is smaller than one ulp of the
// coordinates, so depending on rounding the computed x separation lands a hair
// below zero and the axis is excluded from the separating-normal accumulation.
// SAT then reports the pure (0, 1) normal, and negating it in
// `closest_points_cuboid_cuboid` produces a *negative zero* x component:
// `copy_sign_to` honors that sign and selects cuboid2's -x corner — the corner
// farthest from cuboid1 (which sits at negative x), yielding a distance of
// exactly 2 * half_extents2.x instead of ~1e-6.
#[test]
fn closest_points_cuboid_cuboid_axis_aligned_subulp_gap() {
    let c1 = Cuboid::new(Vector::new(0.5386359, 0.8874373));
    let c2 = Cuboid::new(Vector::new(1.0999209, 0.17219427));
    let pos12 = Pose::translation(-1.6385567, 1.0596324);

    // Both axis gaps are within a couple of ulps of zero, so the true distance is
    // at most a few microns; only assert a loose upper bound.
    match closest_points_cuboid_cuboid(&pos12, &c1, &c2, 0.5) {
        ClosestPoints::WithinMargin(p1, p2) => {
            let dist = (p1 - pos12 * p2).length();
            assert!(dist < 1.0e-5, "closest points are {dist} apart");
        }
        ClosestPoints::Intersecting => {} // Also acceptable at this scale.
        ClosestPoints::Disjoint => panic!("expected WithinMargin, got Disjoint"),
    }

    let dist = query::distance(&Pose::IDENTITY, &c1, &pos12, &c2)
        .unwrap()
        .distance;
    assert!(dist < 1.0e-5, "distance = {dist}");
}
