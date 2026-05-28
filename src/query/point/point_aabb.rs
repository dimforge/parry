use crate::bounding_volume::Aabb;
use crate::math::{Real, Vector, VectorExt, DIM};
use crate::num::Zero;
use crate::query::{PointProjection, PointQuery};
use crate::shape::FeatureId;
use crate::utils::WSign;

impl Aabb {
    fn do_project_local_point(&self, pt: Vector, solid: bool) -> (bool, Vector, Vector) {
        let mins_pt = self.mins - pt;
        let pt_maxs = pt - self.maxs;
        let shift = mins_pt.max(Vector::ZERO) - pt_maxs.max(Vector::ZERO);

        if shift != Vector::ZERO {
            (false, pt + shift, shift)
        } else if solid {
            (true, pt, Vector::ZERO)
        } else {
            // Projection for the case where the point is inside the box.
            let centered_pt = pt - self.center();
            let pt_sgn_with_zero = centered_pt.copy_sign_to(Vector::splat(1.0));
            // This is the sign of pt, or -1 for components that were zero.
            // This bias is arbitrary (we could have picked +1), but we picked it so
            // it matches the previous implementation.
            let pt_sgn = pt_sgn_with_zero + (pt_sgn_with_zero.abs() - Vector::ONE);
            let diff = self.half_extents() - pt_sgn * centered_pt;

            #[cfg(feature = "dim2")]
            let shift = {
                let pick_x = diff.x <= diff.y;
                let shift_x = Vector::new(diff.x * pt_sgn.x, 0.0);
                let shift_y = Vector::new(0.0, diff.y * pt_sgn.y);
                if pick_x {
                    shift_x
                } else {
                    shift_y
                }
            };

            #[cfg(feature = "dim3")]
            let shift = {
                let pick_x = diff.x <= diff.y && diff.x <= diff.z;
                let pick_y = diff.y <= diff.x && diff.y <= diff.z;
                let shift_x = Vector::new(diff.x * pt_sgn.x, 0.0, 0.0);
                let shift_y = Vector::new(0.0, diff.y * pt_sgn.y, 0.0);
                let shift_z = Vector::new(0.0, 0.0, diff.z * pt_sgn.z);
                if pick_x {
                    shift_x
                } else if pick_y {
                    shift_y
                } else {
                    shift_z
                }
            };

            (true, pt + shift, shift)
        }
    }
}

impl PointQuery for Aabb {
    #[inline]
    fn project_local_point(&self, pt: Vector, solid: bool) -> PointProjection {
        let (inside, ls_pt, _) = self.do_project_local_point(pt, solid);
        PointProjection::new(inside, ls_pt)
    }

    #[allow(unused_assignments)] // For last_zero_shift which is used only in 3D.
    #[allow(unused_variables)] // For last_zero_shift which is used only in 3D.
    #[inline]
    fn project_local_point_and_get_feature(&self, pt: Vector) -> (PointProjection, FeatureId) {
        let (inside, ls_pt, shift) = self.do_project_local_point(pt, false);
        let proj = PointProjection::new(inside, ls_pt);
        let mut nzero_shifts = 0;
        let mut last_zero_shift = 0;
        let mut last_not_zero_shift = 0;

        for i in 0..DIM {
            if shift.vget(i).is_zero() {
                nzero_shifts += 1;
                last_zero_shift = i;
            } else {
                last_not_zero_shift = i;
            }
        }

        if nzero_shifts == DIM {
            for i in 0..DIM {
                if ls_pt.vget(i) > self.maxs.vget(i) - crate::math::DEFAULT_EPSILON {
                    return (proj, FeatureId::Face(i as u32));
                }
                if ls_pt.vget(i) <= self.mins.vget(i) + crate::math::DEFAULT_EPSILON {
                    return (proj, FeatureId::Face((i + DIM) as u32));
                }
            }

            (proj, FeatureId::Unknown)
        } else if nzero_shifts == DIM - 1 {
            // On a 3D face.
            if ls_pt.vget(last_not_zero_shift) < self.center().vget(last_not_zero_shift) {
                (proj, FeatureId::Face((last_not_zero_shift + DIM) as u32))
            } else {
                (proj, FeatureId::Face(last_not_zero_shift as u32))
            }
        } else {
            // On a vertex or edge.
            let mut id = 0;
            let center = self.center();

            for i in 0..DIM {
                if ls_pt.vget(i) < center.vget(i) {
                    id |= 1 << i;
                }
            }

            #[cfg(feature = "dim3")]
            {
                if nzero_shifts == 0 {
                    (proj, FeatureId::Vertex(id))
                } else {
                    (proj, FeatureId::Edge((id << 2) | (last_zero_shift as u32)))
                }
            }

            #[cfg(feature = "dim2")]
            {
                (proj, FeatureId::Vertex(id))
            }
        }
    }

    #[inline]
    fn distance_to_local_point(&self, pt: Vector, solid: bool) -> Real {
        let mins_pt = self.mins - pt;
        let pt_maxs = pt - self.maxs;
        let shift = mins_pt.max(pt_maxs).max(Vector::ZERO);

        if solid || shift != Vector::ZERO {
            shift.length()
        } else {
            // TODO: optimize that.
            -pt.distance(self.project_local_point(pt, solid).point)
        }
    }
}
