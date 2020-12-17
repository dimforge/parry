use crate::math::Real;
use na;

/// Computes the median of a set of values.
#[inline]
pub fn median(vals: &mut [Real]) -> Real {
    assert!(vals.len() > 0, "Cannot compute the median of zero values.");

    vals.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let n = vals.len();

    if n % 2 == 0 {
        (vals[n / 2 - 1] + vals[n / 2]) / na::convert::<f64, Real>(2.0)
    } else {
        vals[n / 2]
    }
}
