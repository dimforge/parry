use crate::math::{Isometry, Real};
use crate::shape::{MassProperties, Shape};
use std::sync::Arc;

impl MassProperties {
    pub(crate) fn from_compound(
        density: Real,
        shapes: &[(Isometry<Real>, Arc<dyn Shape>)],
    ) -> Self {
        shapes
            .iter()
            .map(|s| s.1.mass_properties(density).transform_by(&s.0))
            .sum()
    }
}
