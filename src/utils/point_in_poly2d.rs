use crate::math::Real;
use na::Point2;
use num::Zero;

/// Tests if the given point is inside a convex polygon with arbitrary orientation.
///
/// The polygon is assumed to be closed, i.e., first and last point of the polygon are implicitly
/// assumed to be connected by an edge.
pub fn point_in_convex_poly2d(pt: &Point2<Real>, poly: &[Point2<Real>]) -> bool {
    if poly.is_empty() {
        false
    } else {
        let mut sign = 0.0;

        for i1 in 0..poly.len() {
            let i2 = (i1 + 1) % poly.len();
            let seg_dir = poly[i2] - poly[i1];
            let dpt = pt - poly[i1];
            let perp = dpt.perp(&seg_dir);

            if sign.is_zero() {
                sign = perp;
            } else if sign * perp < 0.0 {
                return false;
            }
        }

        true
    }
}

/// Tests if the given point is inside an arbitrary closed polygon with arbitrary orientation,
/// using a counting winding strategy.
///
/// The polygon is assumed to be closed, i.e., first and last point of the polygon are implicitly
/// assumed to be connected by an edge.
///
/// This handles concave polygons. For a function dedicated to convex polygons, see [`point_in_convex_poly2d`].
pub fn point_in_poly2d(pt: &Point2<Real>, poly: &[Point2<Real>]) -> bool {
    if poly.is_empty() {
        return false;
    }

    let mut winding = 0i32;

    for (i, a) in poly.iter().enumerate() {
        let b = poly[(i + 1) % poly.len()];
        let seg_dir = b - a;
        let dpt = pt - a;
        let perp = dpt.perp(&seg_dir);
        winding += match (dpt.y >= 0.0, b.y > pt.y) {
            (true, true) if perp < 0.0 => 1,
            (false, false) if perp > 0.0 => 1,
            _ => 0,
        };
    }

    winding % 2 == 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_in_poly2d_self_intersecting() {
        let poly = [
            [-1.0, -1.0],
            [0.0, -1.0],
            [0.0, 1.0],
            [-2.0, 1.0],
            [-2.0, -2.0],
            [1.0, -2.0],
            [1.0, 2.0],
            [-1.0, 2.0],
        ]
        .map(Point2::from);
        assert!(!point_in_poly2d(&[-0.5, -0.5].into(), &poly));
        assert!(point_in_poly2d(&[0.5, -0.5].into(), &poly));
    }

    #[test]
    fn point_in_poly2d_concave() {
        let poly = [
            [615.474_2, 279.412_08],
            [617.959_5, 281.897_37],
            [624.172_7, 288.731_93],
            [626.658, 292.459_87],
            [634.735_2, 302.401_06],
            [637.841_86, 306.750_34],
            [642.812_44, 312.963_56],
            [652.753_6, 330.981_93],
            [654.617_6, 334.709_87],
            [661.452_15, 349.000_34],
            [666.422_7, 360.184_14],
            [670.150_7, 367.640_05],
            [675.121_3, 381.309_14],
            [678.227_9, 391.250_34],
            [681.334_5, 402.434_14],
            [683.819_8, 414.239_32],
            [685.062_44, 422.316_53],
            [685.683_8, 431.015],
            [686.305_1, 442.820_2],
            [685.683_8, 454.004_03],
            [683.819_8, 460.838_56],
            [679.470_5, 470.779_72],
            [674.499_94, 480.720_95],
            [670.150_7, 486.934_14],
            [662.073_5, 497.496_64],
            [659.588_2, 499.360_66],
            [653.374_94, 503.709_9],
            [647.783, 506.195_2],
            [642.812_44, 507.437_8],
            [631.628_6, 508.059_14],
            [621.066_1, 508.059_14],
            [605.533, 508.059_14],
            [596.213_2, 508.059_14],
            [586.893_3, 508.059_14],
            [578.816_1, 508.059_14],
            [571.360_2, 506.195_2],
            [559.555_1, 499.360_66],
            [557.069_76, 497.496_64],
            [542.158, 484.448_85],
            [534.702_15, 476.371_64],
            [532.838_2, 473.886_35],
            [527.246_3, 466.430_48],
            [522.275_7, 450.897_4],
            [521.654_36, 444.062_8],
            [521.033, 431.636_35],
            [521.654_36, 422.937_8],
            [523.518_3, 409.268_74],
            [527.246_3, 397.463_56],
            [532.838_2, 385.658_42],
            [540.915_4, 373.231_93],
            [547.128_6, 365.776_06],
            [559.555_1, 354.592_22],
            [573.224_2, 342.165_77],
            [575.709_5, 339.680_48],
            [584.408, 331.603_27],
            [597.455_8, 317.312_84],
            [601.805_1, 311.720_92],
            [607.397, 303.643_7],
            [611.746_3, 296.187_84],
            [614.231_57, 288.110_63],
            [615.474_2, 280.654_72],
            [615.474_2, 279.412_08],
        ]
        .map(Point2::from);
        let pt = Point2::from([596.018_2, 427.916_3]);
        assert!(point_in_poly2d(&pt, &poly));
    }

    #[test]
    #[cfg(all(feature = "dim2", feature = "alloc"))]
    fn point_in_poly2d_concave_exact_vertex_bug() {
        let poly = crate::shape::Ball::new(1.0).to_polyline(10);
        assert!(point_in_poly2d(&Point2::origin(), &poly));
        assert!(point_in_poly2d(&Point2::new(-0.25, 0.0), &poly));
        assert!(point_in_poly2d(&Point2::new(0.25, 0.0), &poly));
        assert!(point_in_poly2d(&Point2::new(0.0, -0.25), &poly));
        assert!(point_in_poly2d(&Point2::new(0.0, 0.25), &poly));
        assert!(!point_in_poly2d(&Point2::new(-2.0, 0.0), &poly));
        assert!(!point_in_poly2d(&Point2::new(2.0, 0.0), &poly));
        assert!(!point_in_poly2d(&Point2::new(0.0, -2.0), &poly));
        assert!(!point_in_poly2d(&Point2::new(0.0, 2.0), &poly));
    }
}
