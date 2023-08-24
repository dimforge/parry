use crate::math::{Isometry, Point, Real, Rotation, Vector};
use crate::shape::{Segment, SupportMap};
use na::Unit;

#[cfg(feature = "std")]
use either::Either;

#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "bytemuck", derive(bytemuck::Pod, bytemuck::Zeroable))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[cfg_attr(feature = "cuda", derive(cust_core::DeviceCopy))]
#[repr(C)]
/// A tapered capsule shape defined as a round segment.
pub struct TaperedCapsule {
    /// The endpoints of the capsule’s principal axis.
    pub segment: Segment,
    /// The radius at the end of the `a` point of the capsule segment.
    pub radius_a: Real,
    /// The radius at the end of the `b` point of the capsule segment.
    pub radius_b: Real,
}

impl TaperedCapsule {
    /// Creates a new tapered capsule aligned with the `x` axis and with the given half-height and radii.
    pub fn new_x(half_height: Real, radius_a: Real, radius_b: Real) -> Self {
        let b = Point::from(Vector::x() * half_height);
        Self::new(-b, b, radius_a, radius_b)
    }

    /// Creates a new tapered capsule aligned with the `y` axis and with the given half-height and radii.
    pub fn new_y(half_height: Real, radius_a: Real, radius_b: Real) -> Self {
        let b = Point::from(Vector::y() * half_height);
        Self::new(-b, b, radius_a, radius_b)
    }

    /// Creates a new tapered capsule aligned with the `z` axis and with the given half-height and radii.
    #[cfg(feature = "dim3")]
    pub fn new_z(half_height: Real, radius_a: Real, radius_b: Real) -> Self {
        let b = Point::from(Vector::z() * half_height);
        Self::new(-b, b, radius_a, radius_b)
    }

    /// Creates a new tapered capsule defined as the segment between `a` and `b` and with the given radii
    /// respective to each point.
    pub fn new(a: Point<Real>, b: Point<Real>, radius_a: Real, radius_b: Real) -> Self {
        let segment = Segment::new(a, b);
        Self {
            segment,
            radius_a,
            radius_b,
        }
    }

    /// The height of this capsule.
    pub fn height(&self) -> Real {
        (self.segment.b - self.segment.a).norm()
    }

    /// The half-height of this capsule.
    pub fn half_height(&self) -> Real {
        self.height() / 2.0
    }

    /// The center of this capsule.
    pub fn center(&self) -> Point<Real> {
        na::center(&self.segment.a, &self.segment.b)
    }

    /// Creates a new capsule equal to `self` with all its endpoints transformed by `pos`.
    pub fn transform_by(&self, pos: &Isometry<Real>) -> Self {
        Self::new(
            pos * self.segment.a,
            pos * self.segment.b,
            self.radius_a,
            self.radius_b,
        )
    }

    /// The transformation such that `t * Y` is collinear with `b - a` and `t * origin` equals
    /// the capsule's center.
    pub fn canonical_transform(&self) -> Isometry<Real> {
        let tra = self.center().coords;
        let rot = self.rotation_wrt_y();
        Isometry::from_parts(tra.into(), rot)
    }

    /// The rotation `r` such that `r * Y` is collinear with `b - a`.
    pub fn rotation_wrt_y(&self) -> Rotation<Real> {
        let mut dir = self.segment.b - self.segment.a;
        if dir.y < 0.0 {
            dir = -dir;
        }

        #[cfg(feature = "dim2")]
        {
            Rotation::rotation_between(&Vector::y(), &dir)
        }

        #[cfg(feature = "dim3")]
        {
            Rotation::rotation_between(&Vector::y(), &dir).unwrap_or(Rotation::identity())
        }
    }

    /// The transform `t` such that `t * Y` is collinear with `b - a` and such that `t * origin = (b + a) / 2.0`.
    pub fn transform_wrt_y(&self) -> Isometry<Real> {
        let rot = self.rotation_wrt_y();
        Isometry::from_parts(self.center().coords.into(), rot)
    }

    /// Computes a scaled version of this tapered capsule.
    ///
    /// If the scaling factor is non-uniform, then it can’t be represented as
    /// capsule. Instead, a convex polygon approximation (with `nsubdivs`
    /// subdivisions) is returned. Returns `None` if that approximation had degenerate
    /// normals (for example if the scaling factor along one axis is zero).
    #[cfg(all(feature = "dim2", feature = "std"))]
    pub fn scaled(
        self,
        scale: &Vector<Real>,
        nsubdivs: u32,
    ) -> Option<Either<Self, super::ConvexPolygon>> {
        if scale.x != scale.y {
            unimplemented!();
            /*
                       // The scaled shape is not a capsule.
                       let mut vtx = self.to_polyline(nsubdivs);
                       vtx.iter_mut()
                           .for_each(|pt| pt.coords = pt.coords.component_mul(&scale));
                       Some(Either::Right(super::ConvexPolygon::from_convex_polyline(
                           vtx,
                       )?))
            */
        } else {
            let uniform_scale = scale.x;
            Some(Either::Left(Self::new(
                self.segment.a * uniform_scale,
                self.segment.b * uniform_scale,
                self.radius_a * uniform_scale.abs(),
                self.radius_b * uniform_scale.abs(),
            )))
        }
    }

    /// Computes a scaled version of this capsule.
    ///
    /// If the scaling factor is non-uniform, then it can’t be represented as
    /// capsule. Instead, a convex polygon approximation (with `nsubdivs`
    /// subdivisions) is returned. Returns `None` if that approximation had degenerate
    /// normals (for example if the scaling factor along one axis is zero).
    #[cfg(all(feature = "dim3", feature = "std"))]
    pub fn scaled(
        self,
        scale: &Vector<Real>,
        nsubdivs: u32,
    ) -> Option<Either<Self, super::ConvexPolyhedron>> {
        if scale.x != scale.y || scale.x != scale.z || scale.y != scale.z {
            unimplemented!();
            /*
            // The scaled shape is not a capsule.
            let (mut vtx, idx) = self.to_trimesh(nsubdivs, nsubdivs);
            vtx.iter_mut()
                .for_each(|pt| pt.coords = pt.coords.component_mul(&scale));
            Some(Either::Right(super::ConvexPolyhedron::from_convex_mesh(
                vtx, &idx,
            )?))
            */
        } else {
            let uniform_scale = scale.x;
            Some(Either::Left(Self::new(
                self.segment.a * uniform_scale,
                self.segment.b * uniform_scale,
                self.radius_a * uniform_scale.abs(),
                self.radius_b * uniform_scale.abs(),
            )))
        }
    }
}

impl SupportMap for TaperedCapsule {
    fn local_support_point(&self, dir: &Vector<Real>) -> Point<Real> {
        let dir = Unit::try_new(*dir, 0.0).unwrap_or(Vector::y_axis());
        self.local_support_point_toward(&dir)
    }

    fn local_support_point_toward(&self, dir: &Unit<Vector<Real>>) -> Point<Real> {
        if dir.dot(&self.segment.a.coords) > dir.dot(&self.segment.b.coords) {
            self.segment.a + **dir * self.radius_a
        } else {
            self.segment.b + **dir * self.radius_b
        }
    }
}
