use na::Unit;

use crate::math::{Isometry, Point, Real, Vector};
use crate::query::{DefaultQueryDispatcher, QueryDispatcher, Unsupported};
use crate::shape::Shape;

/// The status of the time-of-impact computation algorithm.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ShapeCastStatus {
    /// The shape-casting algorithm ran out of iterations before achieving convergence.
    ///
    /// The content of the `ShapeCastHit` will still be a conservative approximation of the actual result so
    /// it is often fine to interpret this case as a success.
    OutOfIterations,
    /// The shape-casting algorithm converged successfully.
    Converged,
    /// Something went wrong during the shape-casting, likely due to numerical instabilities.
    ///
    /// The content of the `ShapeCastHit` will still be a conservative approximation of the actual result so
    /// it is often fine to interpret this case as a success.
    Failed,
    /// The two shape already overlap, or are separated by a distance smaller than
    /// [`ShapeCastOptions::target_distance`] at the time 0.
    ///
    /// The witness points and normals provided by the `ShapeCastHit` will have unreliable values unless
    /// [`ShapeCastOptions::compute_impact_geometry_on_penetration`] was set to `true` when calling
    /// the time-of-impact function.
    PenetratingOrWithinTargetDist,
}

/// The result of a shape casting..
#[derive(Copy, Clone, Debug)]
pub struct ShapeCastHit {
    /// The time at which the objects touch.
    pub time_of_impact: Real,
    /// The local-space closest point on the first shape at the time of impact.
    ///
    /// This value is unreliable if `status` is [`ShapeCastStatus::PenetratingOrWithinTargetDist`]
    /// and [`ShapeCastOptions::compute_impact_geometry_on_penetration`] was set to `false`.
    pub witness1: Point<Real>,
    /// The local-space closest point on the second shape at the time of impact.
    ///
    /// This value is unreliable if `status` is [`ShapeCastStatus::PenetratingOrWithinTargetDist`]
    /// and both [`ShapeCastOptions::compute_impact_geometry_on_penetration`] was set to `false`
    /// when calling the time-of-impact function.
    pub witness2: Point<Real>,
    /// The local-space outward normal on the first shape at the time of impact.
    ///
    /// This value is unreliable if `status` is [`ShapeCastStatus::PenetratingOrWithinTargetDist`]
    /// and both [`ShapeCastOptions::compute_impact_geometry_on_penetration`] was set to `false`
    /// when calling the time-of-impact function.
    pub normal1: Unit<Vector<Real>>,
    /// The local-space outward normal on the second shape at the time of impact.
    ///
    /// This value is unreliable if `status` is [`ShapeCastStatus::PenetratingOrWithinTargetDist`]
    /// and both [`ShapeCastOptions::compute_impact_geometry_on_penetration`] was set to `false`
    /// when calling the time-of-impact function.
    pub normal2: Unit<Vector<Real>>,
    /// The way the shape-casting algorithm terminated.
    pub status: ShapeCastStatus,
}

impl ShapeCastHit {
    /// Swaps every data of this shape-casting result such that the role of both shapes are swapped.
    ///
    /// In practice, this makes it so that `self.witness1` and `self.normal1` are swapped with
    /// `self.witness2` and `self.normal2`.
    pub fn swapped(self) -> Self {
        Self {
            time_of_impact: self.time_of_impact,
            witness1: self.witness2,
            witness2: self.witness1,
            normal1: self.normal2,
            normal2: self.normal1,
            status: self.status,
        }
    }

    /// Transform `self.witness1` and `self.normal1` by `pos`.
    pub fn transform1_by(&self, pos: &Isometry<Real>) -> Self {
        Self {
            time_of_impact: self.time_of_impact,
            witness1: pos * self.witness1,
            witness2: self.witness2,
            normal1: pos * self.normal1,
            normal2: self.normal2,
            status: self.status,
        }
    }
}

/// Configuration for controlling the behavior of time-of-impact (i.e. shape-casting) calculations.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ShapeCastOptions {
    /// The maximum time-of-impacts that can be computed.
    ///
    /// Any impact occurring after this time will be ignored.
    pub max_time_of_impact: Real,
    /// The shapes will be considered as impacting as soon as their distance is smaller or
    /// equal to this target distance. Must be positive or zero.
    ///
    /// If the shapes are separated by a distance smaller than `target_distance` at time 0, the
    /// calculated witness points and normals are only reliable if
    /// [`Self::compute_impact_geometry_on_penetration`] is set to `true`.
    pub target_distance: Real,
    /// If `false`, the time-of-impact algorithm will automatically discard any impact at time
    /// 0 where the velocity is separating (i.e., the relative velocity is such that the distance
    /// between the objects projected on the impact normal is increasing through time).
    pub stop_at_penetration: bool,
    /// If `true`, witness points and normals will be calculated even when the time-of-impact is 0.
    pub compute_impact_geometry_on_penetration: bool,
}

impl ShapeCastOptions {
    // Constructor for the most common use-case.
    /// Crates a [`ShapeCastOptions`] with the default values except for the maximum time of impact.
    pub fn with_max_time_of_impact(max_time_of_impact: Real) -> Self {
        Self {
            max_time_of_impact,
            ..Default::default()
        }
    }
}

impl Default for ShapeCastOptions {
    fn default() -> Self {
        Self {
            max_time_of_impact: Real::MAX,
            target_distance: 0.0,
            stop_at_penetration: true,
            compute_impact_geometry_on_penetration: true,
        }
    }
}

/// Computes the smallest time when two shapes under translational movement are separated by a
/// distance smaller or equal to `distance`.
///
/// Returns `0.0` if the objects are touching or closer than `options.target_distance`,
/// or penetrating.
pub fn cast_shapes(
    pos1: &Isometry<Real>,
    vel1: &Vector<Real>,
    g1: &dyn Shape,
    pos2: &Isometry<Real>,
    vel2: &Vector<Real>,
    g2: &dyn Shape,
    options: ShapeCastOptions,
) -> Result<Option<ShapeCastHit>, Unsupported> {
    let pos12 = pos1.inv_mul(pos2);
    let vel12 = pos1.inverse_transform_vector(&(vel2 - vel1));
    DefaultQueryDispatcher.cast_shapes(&pos12, &vel12, g1, g2, options)
}
