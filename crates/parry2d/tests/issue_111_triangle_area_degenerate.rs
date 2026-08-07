// Regression test for https://github.com/dimforge/parry/issues/111 (2D version).
//
// `Triangle::area` now uses half the perp-product magnitude in 2D, which is
// exact (0.0) when the vertices are bitwise collinear.

use parry2d::math::{Real, Vector};
use parry2d::shape::Triangle;

/// The previous implementation (Kahan's formula on side lengths).
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
    // The issue's collinear triangle, projected to 2D (all points share x, so
    // use the (y, z) coordinates).
    let tri = Triangle::new(
        Vector::new(-2.871, 17.464),
        Vector::new(1.629, 17.464),
        Vector::new(-1.521, 17.464),
    );
    assert_eq!(tri.area(), 0.0);

    // Two identical vertices.
    let tri = Triangle::new(
        Vector::new(2.277, -7.9),
        Vector::new(-0.57, -8.1),
        Vector::new(-0.57, -8.1),
    );
    assert_eq!(tri.area(), 0.0);
}

#[test]
fn area_matches_kahan_for_well_conditioned_triangles() {
    let tris = [
        Triangle::new(
            Vector::new(0.0, 0.0),
            Vector::new(1.0, 0.0),
            Vector::new(0.0, 1.0),
        ),
        Triangle::new(
            Vector::new(1.0, 2.0),
            Vector::new(4.0, 0.0),
            Vector::new(2.0, 5.0),
        ),
        Triangle::new(
            Vector::new(-3.0, 1.0),
            Vector::new(0.0, 4.0),
            Vector::new(5.0, -1.0),
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

    assert_eq!(
        Triangle::new(
            Vector::new(0.0, 0.0),
            Vector::new(1.0, 0.0),
            Vector::new(0.0, 1.0),
        )
        .area(),
        0.5
    );
}
