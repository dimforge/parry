//! Support mapping based Cone shape.

use crate::math::{Point, Real, Vector};
use crate::shape::{FeatureId, SupportMap};
use na::{self, Unit};
use num::Zero;

#[cfg(feature = "alloc")]
use either::Either;

#[cfg(not(feature = "alloc"))]
use na::RealField; // for .copysign()

#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

/// Cone shape with its principal axis aligned with the `y` axis.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "bytemuck", derive(bytemuck::Pod, bytemuck::Zeroable))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[derive(PartialEq, Debug, Copy, Clone)]
#[repr(C)]
pub struct Cone {
    /// The half-height of the cone.
    pub half_height: Real,
    /// The base radius of the cone.
    pub radius: Real,
}

impl Cone {
    /// Creates a new cone.
    ///
    /// # Arguments:
    /// * `half_height` - the half length of the cone along the `y` axis.
    /// * `radius` - the length of the cone along all other axis.
    pub fn new(half_height: Real, radius: Real) -> Cone {
        Cone {
            half_height,
            radius,
        }
    }

    /// Computes a scaled version of this cone.
    ///
    /// If the scaling factor is non-uniform, then it can’t be represented as
    /// cone. Instead, a convex polyhedral approximation (with `nsubdivs`
    /// subdivisions) is returned. Returns `None` if that approximation had degenerate
    /// normals (for example if the scaling factor along one axis is zero).
    #[cfg(feature = "alloc")]
    #[inline]
    pub fn scaled(
        self,
        scale: &Vector<Real>,
        nsubdivs: u32,
    ) -> Option<Either<Self, super::ConvexPolyhedron>> {
        // NOTE: if the y scale is negative, the result cone points downwards,
        //       which can’t be represented with this Cone (without a transform).
        if scale.x != scale.z || scale.y < 0.0 {
            // The scaled shape isn’t a cone.
            let (mut vtx, idx) = self.to_trimesh(nsubdivs);
            vtx.iter_mut()
                .for_each(|pt| pt.coords = pt.coords.component_mul(scale));
            Some(Either::Right(super::ConvexPolyhedron::from_convex_mesh(
                vtx, &idx,
            )?))
        } else {
            Some(Either::Left(Self::new(
                self.half_height * scale.y,
                self.radius * scale.x,
            )))
        }
    }

    /// Computes the normal for given [FeatureId] at a given point.
    ///
    /// Supports only `FeatureId::Face(0)` (base) and `FeatureId::Face(1)` (side).
    pub fn feature_normal_at_point(
        &self,
        feature_id: FeatureId,
        point: &Point<Real>,
    ) -> Option<Unit<Vector<Real>>> {
        match feature_id {
            FeatureId::Vertex(vertex) => {
                // TODO: if top of the cone, return (0,0,1) ?
                None
            }
            // Base
            FeatureId::Face(0) => Some(Unit::new_unchecked(Vector::new(0.0, 0.0, -1.0))),
            // Side
            FeatureId::Face(1) => {
                // TODO: check nan...
                if self.half_height == 0.0 || self.radius == 0.0 {
                    return None;
                }
                fn height_from_right_angle(a: Real, b: Real) -> Real {
                    let hypotenuse = (a.powi(2) + b.powi(2)).sqrt();
                    (a * b) / hypotenuse
                }
                fn opposite_from_adjacent_and_hypotenuse(adjacent: Real, hypotenuse: Real) -> Real {
                    let opp_squared = hypotenuse.powi(2) - adjacent.powi(2);
                    opp_squared.sqrt()
                }
                let height = self.half_height * 2.0;
                let triangle_height = height_from_right_angle(height, self.radius);
                let normal_triangle_opposite_length =
                    opposite_from_adjacent_and_hypotenuse(triangle_height, self.radius);
                let normal_triangle_height =
                    (triangle_height * normal_triangle_opposite_length) / self.radius;
                let normal_triangle_ground_length =
                    opposite_from_adjacent_and_hypotenuse(normal_triangle_height, triangle_height);
                let ground_direction =
                    Vector::new(point.x, 0.0, point.z).try_normalize(Real::EPSILON)?;
                let mut normal = ground_direction * normal_triangle_ground_length;
                normal.y = normal_triangle_height;

                // Normalize
                Some(Unit::try_new(normal, Real::EPSILON)?)
            }
            FeatureId::Face(_) => None,
            FeatureId::Unknown => None,
            FeatureId::Edge(_) => None,
        }
    }
}

impl SupportMap for Cone {
    #[inline]
    fn local_support_point(&self, dir: &Vector<Real>) -> Point<Real> {
        let mut vres = *dir;

        vres[1] = 0.0;

        if vres.normalize_mut().is_zero() {
            vres = na::zero();
            vres[1] = self.half_height.copysign(dir[1]);
        } else {
            vres *= self.radius;
            vres[1] = -self.half_height;

            if dir.dot(&vres) < dir[1] * self.half_height {
                vres = na::zero();
                vres[1] = self.half_height
            }
        }

        Point::from(vres)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn cone_normal() {
        let cone = Cone::new(0.5, 1.0);

        let point = Point::new(1.0, 1.0, 1.0);
        let normal = cone
            .feature_normal_at_point(FeatureId::Face(1), &point)
            .unwrap();
        assert_relative_eq!(
            Vector::new(normal.x, normal.y, normal.z),
            Vector::new(0.5, 0.70710677, 0.5)
        );
        // Higher point, not aligned with normal
        let point = Point::new(1.0, 2.0, 1.0);
        let normal = cone
            .feature_normal_at_point(FeatureId::Face(1), &point)
            .unwrap();
        assert_relative_eq!(
            Vector::new(normal.x, normal.y, normal.z),
            Vector::new(0.5, 0.70710677, 0.5)
        );
        // longer cone
        let cone = Cone::new(1.0, 1.0);
        let point = Point::new(0.5, 2.0, 1.0);
        let normal = cone
            .feature_normal_at_point(FeatureId::Face(1), &point)
            .unwrap();
        assert_relative_eq!(
            Vector::new(normal.x, normal.y, normal.z),
            Vector::new(0.39999998, 0.44721365, 0.79999995)
        );
    }
}
