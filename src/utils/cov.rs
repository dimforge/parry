use crate::math::{Matrix, Point, Real};

/// Computes the covariance matrix of a set of points.
pub fn cov(pts: &[Point<Real>]) -> Matrix<Real> {
    let center = crate::utils::center(pts);
    let mut cov: Matrix<Real> = na::zero();
    let normalizer: Real = 1.0 / (pts.len() as Real);

    for p in pts.iter() {
        let cp = *p - center;
        // NOTE: this is more numerically stable than using cov.syger.
        cov += cp * (cp * normalizer).transpose();
    }

    cov
}
