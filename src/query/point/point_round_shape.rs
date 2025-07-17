use log::warn;

use crate::math::{Point, Real};
use crate::query::gjk::{GjkOptions, VoronoiSimplex};
use crate::query::point::point_query::QueryOptions;
use crate::query::{PointProjection, PointQuery};
use crate::shape::{FeatureId, RoundShape, SupportMap};

// TODO: if PointQuery had a `project_point_with_normal` method, we could just
// call this and adjust the projected point accordingly.
impl<S: SupportMap> PointQuery for RoundShape<S> {
    #[inline]
    fn project_local_point(
        &self,
        point: &Point<Real>,
        solid: bool,
        options: &dyn QueryOptions,
    ) -> PointProjection {
        #[cfg(not(feature = "alloc"))]
        return unimplemented!(
            "The projection of points on a round shape isn't supported without alloc yet."
        );

        #[cfg(feature = "alloc")] // TODO: can’t be used without alloc because of EPA
        {
            let options = if let Some(options) = options.as_any().downcast_ref() {
                options
            } else {
                warn!("Incorrect option passed to project_local_point: using default options.");
                &GjkOptions::default()
            };
            crate::query::details::local_point_projection_on_support_map(
                self,
                &mut VoronoiSimplex::new(),
                point,
                solid,
                options,
            )
        }
    }

    #[inline]
    fn project_local_point_and_get_feature(
        &self,
        point: &Point<Real>,
        options: &dyn QueryOptions,
    ) -> (PointProjection, FeatureId) {
        (
            self.project_local_point(point, false, options),
            FeatureId::Unknown,
        )
    }
}
