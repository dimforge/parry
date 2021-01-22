use crate::mass_properties::MassProperties;
use crate::math::{Isometry, Real};
use crate::shape::SharedShape;
use std::sync::Arc;

impl MassProperties {
    pub fn from_compound(density: Real, shapes: &[(Isometry<Real>, SharedShape)]) -> Self {
        shapes
            .iter()
            .map(|s| s.1.mass_properties(density).transform_by(&s.0))
            .sum()
    }
}
