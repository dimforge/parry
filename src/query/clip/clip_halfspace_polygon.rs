use crate::math::{Point, Real, Vector};
use crate::query::{self, Ray};

/// Cuts a polygon with the given half-space.
///
/// Given the half-space `center` and outward `normal`,
/// this computes the intersecting between the half-space and
/// the polygon. (Note that a point `pt` is considered as inside of
/// the half-space if `normal.dot(&(pt - center)) <= 0.0`.
pub fn clip_halfspace_polygon(
    center: &Point<Real>,
    normal: &Vector<Real>,
    polygon: &[Point<Real>],
    result: &mut Vec<Point<Real>>,
) {
    result.clear();

    if polygon.is_empty() {
        return;
    }

    let keep_point = |pt: &Point<Real>| (pt - center).dot(normal) <= 0.0;
    let last_pt = polygon.last().unwrap();
    let mut last_keep = keep_point(last_pt);

    if last_keep {
        result.push(*last_pt);
    }

    for i in 0..polygon.len() {
        let pt = &polygon[i];
        let keep = keep_point(pt);

        if keep != last_keep {
            // We crossed the plane, so we need
            // to cut the edge.
            let prev_i = if i == 0 { polygon.len() - 1 } else { i - 1 };
            let prev_pt = &polygon[prev_i];
            let ray = Ray::new(*prev_pt, pt - prev_pt);

            if let Some(toi) = query::details::ray_toi_with_halfspace(&center, normal, &ray) {
                if toi > 0.0 && toi < 1.0 {
                    result.push(ray.origin + ray.dir * toi)
                }
            }

            last_keep = keep;
        }

        if keep && i != polygon.len() - 1 {
            result.push(*pt);
        }
    }
}
