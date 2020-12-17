use crate::math::{Isometry, Vector};
use crate::shape::SupportMap;
use na::Unit;

#[allow(dead_code)]
pub fn support_map_support_map_compute_separation(
    sm1: &impl SupportMap,
    sm2: &impl SupportMap,
    pos12: &Isometry<f32>,
    dir1: &Unit<Vector<f32>>,
) -> f32 {
    let p1 = sm1.local_support_point_toward(dir1);
    let p2 = sm2.support_point_toward(pos12, &-*dir1);
    (p2 - p1).dot(dir1)
}
