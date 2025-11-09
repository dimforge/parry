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

/// A 3D cone shape with apex pointing upward along the Y axis.
///
/// A cone is a shape that tapers from a circular base to a point (apex). In Parry,
/// cones are always aligned with the Y axis, with the base at y = -half_height and
/// the apex at y = +half_height.
///
/// # Structure
///
/// - **Axis**: Always aligned with Y axis (apex points up)
/// - **Base**: Circular base at y = -half_height with the given radius
/// - **Apex**: Point at y = +half_height
/// - **Total height**: `2 * half_height`
///
/// # Properties
///
/// - **3D only**: Only available with the `dim3` feature
/// - **Convex**: Yes, cones are convex shapes
/// - **Apex**: Sharp point at the top
/// - **Flat base**: Circular base (not rounded)
/// - **Sharp edge**: Rim where cone surface meets the base
///
/// # Coordinate System
///
/// The cone is centered at the origin with:
/// - Base center at `(0, -half_height, 0)`
/// - Apex at `(0, half_height, 0)`
/// - Circular base in the XZ plane
///
/// # Use Cases
///
/// - Projectile nose cones
/// - Traffic cones and markers
/// - Funnel shapes
/// - Spotlight or vision cone representations
/// - Conical collision bounds
///
/// # Example
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::shape::Cone;
///
/// // Create a cone: base radius 3.0, total height 8.0
/// let cone = Cone::new(4.0, 3.0);
///
/// assert_eq!(cone.half_height, 4.0);
/// assert_eq!(cone.radius, 3.0);
///
/// // The cone:
/// // - Base at y = -4.0 with radius 3.0
/// // - Apex at y = +4.0 (sharp point)
/// // - Total height = 8.0
/// # }
/// ```
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
    /// Half the total height of the cone.
    ///
    /// The cone extends from y = -half_height (base center) to
    /// y = +half_height (apex). Must be positive.
    pub half_height: Real,

    /// The radius of the circular base.
    ///
    /// The base is a circle in the XZ plane at y = -half_height.
    /// Must be positive.
    pub radius: Real,
}

impl Cone {
    /// Creates a new cone with apex pointing upward along the Y axis.
    ///
    /// # Arguments
    ///
    /// * `half_height` - Half the total height (apex to base center distance / 2)
    /// * `radius` - The radius of the circular base
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::shape::Cone;
    ///
    /// // Create a cone with height 6.0 and base radius 2.0
    /// let cone = Cone::new(3.0, 2.0);
    ///
    /// assert_eq!(cone.half_height, 3.0);
    /// assert_eq!(cone.radius, 2.0);
    ///
    /// // The cone structure:
    /// // - Apex at (0, 3, 0)  - the top point
    /// // - Base center at (0, -3, 0)
    /// // - Base radius 2.0 in the XZ plane
    /// // - Total height from apex to base = 6.0
    /// # }
    /// ```
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
    ///
    /// Exposed through [Shape::feature_normal_at_point][crate::shape::Shape::feature_normal_at_point]
    pub fn feature_normal_at(
        &self,
        feature_id: FeatureId,
        point: &Point<Real>,
    ) -> Option<Unit<Vector<Real>>> {
        match feature_id {
            FeatureId::Vertex(_) => {
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
#[cfg(feature = "dim3")]
mod test {
    use approx::RelativeEq;

    use crate::shape::Shape;

    use super::*;

    #[test]
    fn cone_normal() {
        let espilon = 0.00001
        let cone = Cone::new(0.5, 1.0);

        let point = Point::new(1.0, 1.0, 1.0);
        let normal = cone
            .feature_normal_at_point(FeatureId::Face(1), &point)
            .unwrap();
        assert!(Vector::new(normal.x, normal.y, normal.z).relative_eq(
            &Vector::new(0.5, 0.70710677, 0.5),
            espilon,
            espilon
        ),);
        // Higher point, not aligned with normal
        let point = Point::new(1.0, 2.0, 1.0);
        let normal = cone
            .feature_normal_at_point(FeatureId::Face(1), &point)
            .unwrap();
        assert!(Vector::new(normal.x, normal.y, normal.z).relative_eq(
            &Vector::new(0.5, 0.70710677, 0.5),
            espilon,
            espilon
        ),);
        // longer cone
        let cone = Cone::new(1.0, 1.0);
        let point = Point::new(0.5, 2.0, 1.0);
        let normal = cone
            .feature_normal_at_point(FeatureId::Face(1), &point)
            .unwrap();
        assert!(Vector::new(normal.x, normal.y, normal.z).relative_eq(
            &Vector::new(0.39999998, 0.44721365, 0.79999995),
            espilon,
            espilon
        ),);
    }
}
