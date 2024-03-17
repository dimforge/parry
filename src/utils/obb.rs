use crate::math::*;
use crate::shape::Cuboid;

/// Computes an oriented bounding box for the given set of points.
///
/// The returned OBB is not guaranteed to be the smallest enclosing OBB.
/// Though it should be a pretty good on for most purposes.
pub fn obb(pts: &[Point]) -> (Isometry, Cuboid) {
    let cov = crate::utils::cov(pts);
    let mut eigv = cov.symmetric_eigen().eigenvectors;

    if eigv.determinant() < 0.0 {
        eigv = -eigv;
    }

    let mut mins = Vector::repeat(Real::MAX);
    let mut maxs = Vector::repeat(-Real::MAX);

    for pt in pts {
        for i in 0..DIM {
            let dot = eigv.column(i).dot(pt.as_vector());
            mins[i] = mins[i].min(dot);
            maxs[i] = maxs[i].max(dot);
        }
    }

    #[cfg(all(feature = "dim2", feature = "linalg-glam"))]
    let rot = eigv;
    #[cfg(all(feature = "dim2", feature = "linalg-nalgebra"))]
    let rot = Rotation::from_rotation_matrix(&na::Rotation2::from_matrix_unchecked(eigv));
    #[cfg(feature = "dim3")]
    let rot = Rotation::from_rotation_matrix(&RotationMatrix::from_matrix_unchecked(eigv));

    (
        Isometry::from_parts((rot * (maxs + mins) / 2.0).into(), rot),
        Cuboid::new((maxs - mins) / 2.0),
    )
}
