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
            [615.4741821289063, 279.4120788574219],
            [617.95947265625, 281.8973693847656],
            [624.1727294921875, 288.73193359375],
            [626.6580200195313, 292.4598693847656],
            [634.7352294921875, 302.40106201171875],
            [637.8418579101563, 306.7503356933594],
            [642.8124389648438, 312.96356201171875],
            [652.7536010742188, 330.98193359375],
            [654.6176147460938, 334.7098693847656],
            [661.4521484375, 349.0003356933594],
            [666.4227294921875, 360.18414306640625],
            [670.1506958007813, 367.6400451660156],
            [675.1212768554688, 381.30914306640625],
            [678.2279052734375, 391.2503356933594],
            [681.33447265625, 402.43414306640625],
            [683.81982421875, 414.23931884765625],
            [685.0624389648438, 422.3165283203125],
            [685.6837768554688, 431.0150146484375],
            [686.3051147460938, 442.8201904296875],
            [685.6837768554688, 454.0040283203125],
            [683.81982421875, 460.83856201171875],
            [679.4705200195313, 470.77972412109375],
            [674.4999389648438, 480.720947265625],
            [670.1506958007813, 486.93414306640625],
            [662.073486328125, 497.49664306640625],
            [659.5881958007813, 499.36065673828125],
            [653.3749389648438, 503.70989990234375],
            [647.7830200195313, 506.1951904296875],
            [642.8124389648438, 507.43780517578125],
            [631.6286010742188, 508.05914306640625],
            [621.0661010742188, 508.05914306640625],
            [605.5330200195313, 508.05914306640625],
            [596.2131958007813, 508.05914306640625],
            [586.893310546875, 508.05914306640625],
            [578.8161010742188, 508.05914306640625],
            [571.3602294921875, 506.1951904296875],
            [559.5551147460938, 499.36065673828125],
            [557.0697631835938, 497.49664306640625],
            [542.1580200195313, 484.4488525390625],
            [534.7021484375, 476.37164306640625],
            [532.8381958007813, 473.8863525390625],
            [527.2462768554688, 466.43048095703125],
            [522.2756958007813, 450.89739990234375],
            [521.6543579101563, 444.06280517578125],
            [521.0330200195313, 431.6363525390625],
            [521.6543579101563, 422.93780517578125],
            [523.518310546875, 409.26873779296875],
            [527.2462768554688, 397.46356201171875],
            [532.8381958007813, 385.6584167480469],
            [540.9154052734375, 373.23193359375],
            [547.1286010742188, 365.77606201171875],
            [559.5551147460938, 354.59222412109375],
            [573.2241821289063, 342.165771484375],
            [575.70947265625, 339.68048095703125],
            [584.4080200195313, 331.603271484375],
            [597.455810546875, 317.3128356933594],
            [601.8051147460938, 311.7209167480469],
            [607.39697265625, 303.6437072753906],
            [611.7462768554688, 296.1878356933594],
            [614.2315673828125, 288.1106262207031],
            [615.4741821289063, 280.65472412109375],
            [615.4741821289063, 279.4120788574219],
        ]
        .map(Point2::from);
        let pt = Point2::from([596.0181884765625, 427.9162902832031]);
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
