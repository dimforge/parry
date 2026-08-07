// Regression test for https://github.com/dimforge/parry/issues/111
//
// `Triangle::area` used Kahan's formula on the rounded f32 side lengths, returning
// wildly wrong values (~1e-3) for degenerate triangles; it now uses half the
// cross-product magnitude, exact (0.0) for bitwise-collinear vertices.

use parry3d::math::{Real, Vector};
use parry3d::shape::Triangle;

/// The previous implementation (Kahan's formula on side lengths), kept here to
/// check that the new formula agrees with it on well-conditioned triangles.
fn kahan_area(tri: &Triangle) -> Real {
    let mut s = [
        tri.b.distance(tri.a),
        tri.c.distance(tri.b),
        tri.a.distance(tri.c),
    ];
    s.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let (c, b, a) = (s[0], s[1], s[2]); // a >= b >= c

    let sqr = (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c));
    sqr.max(0.0).sqrt() * 0.25
}

#[test]
fn degenerate_triangle_area_is_exactly_zero() {
    // Exact values from the issue: three collinear points (the old formula
    // returned ~0.0010679931).
    let tri = Triangle::new(
        Vector::new(1.811, -2.871, 17.464),
        Vector::new(1.811, 1.629, 17.464),
        Vector::new(1.811, -1.521, 17.464),
    );

    assert_eq!(tri.area(), 0.0);

    // Two identical vertices (from issue #76).
    let tri = Triangle::new(
        Vector::new(2.27699995, -7.9000001, 16.3180008),
        Vector::new(-0.569999993, -8.10000038, 16.6070004),
        Vector::new(-0.569999993, -8.10000038, 16.6070004),
    );

    assert_eq!(tri.area(), 0.0);
}

#[test]
fn area_matches_kahan_for_well_conditioned_triangles() {
    let tris = [
        Triangle::new(
            Vector::new(0.0, 0.0, 0.0),
            Vector::new(1.0, 0.0, 0.0),
            Vector::new(0.0, 1.0, 0.0),
        ),
        Triangle::new(
            Vector::new(1.0, 2.0, 3.0),
            Vector::new(4.0, 0.0, -1.0),
            Vector::new(2.0, 5.0, 1.0),
        ),
        Triangle::new(
            Vector::new(-3.0, 1.0, 2.0),
            Vector::new(0.0, 4.0, -2.0),
            Vector::new(5.0, -1.0, 3.0),
        ),
        Triangle::new(
            Vector::new(0.1, 0.2, 0.3),
            Vector::new(-0.4, 0.5, 0.1),
            Vector::new(0.3, -0.2, 0.6),
        ),
    ];

    for tri in &tris {
        let expected = kahan_area(tri);
        let area = tri.area();
        assert!(
            (area - expected).abs() <= expected * 1.0e-6,
            "area {area} != kahan area {expected}"
        );
    }

    // Sanity check on an exactly-known area.
    assert_eq!(
        Triangle::new(
            Vector::new(0.0, 0.0, 0.0),
            Vector::new(1.0, 0.0, 0.0),
            Vector::new(0.0, 1.0, 0.0),
        )
        .area(),
        0.5
    );
}
