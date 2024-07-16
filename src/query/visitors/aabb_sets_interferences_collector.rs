use crate::math::{Isometry, Matrix, Real};

/// Spatial partitioning data structure visitor collecting interferences with a given bounding volume.
pub struct AabbSetsInterferencesCollector<'a, T: 'a> {
    /// The transform from the local-space of the second bounding volumes to the local space of the first.
    pub ls_m2: &'a Isometry<Real>,
    /// The absolute value of the rotation matrix representing `ls_m2.rotation`.
    ///
    /// Equals to `ls_m2.rotation.to_rotation.matrix().matrix().abs()`.
    pub ls_m2_abs_rot: &'a Matrix<Real>,
    /// A tolerance applied to the interference tests.
    ///
    /// Aabb pairs closer than `tolerance` will be reported as intersecting.
    pub tolerance: Real,
    /// The data contained by the nodes with bounding volumes intersecting `self.bv`.
    pub collector: &'a mut Vec<(T, T)>,
}

impl<'a, T> AabbSetsInterferencesCollector<'a, T> {
    /// Creates a new `AabbSetsInterferencesCollector`.
    #[inline]
    pub fn new(
        tolerance: Real,
        ls_m2: &'a Isometry<Real>,
        ls_m2_abs_rot: &'a Matrix<Real>,
        collector: &'a mut Vec<(T, T)>,
    ) -> AabbSetsInterferencesCollector<'a, T> {
        AabbSetsInterferencesCollector {
            tolerance,
            ls_m2,
            ls_m2_abs_rot,
            collector,
        }
    }
}

// impl<'a, T: Clone> SimultaneousVisitor<T, Aabb> for AabbSetsInterferencesCollector<'a, T> {
//     #[inline]
//     fn visit(
//         &mut self,
//         left_bv: &Aabb,
//         left_data: Option<&T>,
//         right_bv: &Aabb,
//         right_data: Option<&T>,
//     ) -> VisitStatus {
//         let ls_right_bv = Aabb::from_half_extents(
//             self.ls_m2 * right_bv.center(),
//             self.ls_m2_abs_rot * right_bv.half_extents() + Vector::repeat(self.tolerance),
//         );
//
//         if left_bv.intersects(&ls_right_bv) {
//             if let (Some(a), Some(b)) = (left_data, right_data) {
//                 self.collector.push((a.clone(), b.clone()))
//             }
//
//             VisitStatus::Continue
//         } else {
//             VisitStatus::Stop
//         }
//     }
// }
