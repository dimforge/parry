use crate::math::{Isometry, Point, Vector};
use crate::query::{
    gjk::VoronoiSimplex, PointProjection, PointQuery, Ray, RayCast, RayIntersection,
};
use crate::shape::{FeatureId, SupportMap};
use na::Unit;

#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[derive(Copy, Clone, Debug)]
/// A shape with rounded borders.
pub struct RoundShape<S> {
    /// The shape being rounded.
    pub base_shape: S,
    /// The radius of the rounded border.
    pub border_radius: f32,
}

impl<S: SupportMap> SupportMap for RoundShape<S> {
    fn local_support_point(&self, dir: &Vector<f32>) -> Point<f32> {
        self.local_support_point_toward(&Unit::new_normalize(*dir))
    }

    fn local_support_point_toward(&self, dir: &Unit<Vector<f32>>) -> Point<f32> {
        self.base_shape.local_support_point_toward(dir) + **dir * self.border_radius
    }
}
