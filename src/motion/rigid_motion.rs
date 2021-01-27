use crate::math::{Isometry, Point, Real, Translation, Vector};

/// A continuous rigid motion.
///
/// This is a function, assumed to be continuous, that, given a parameter `t` returns a direct isometry.
/// Mathematically speaking this is a one-parameter curve on the space of direct isometries. This curve
/// should have a continuity of at least `C0`.
pub trait RigidMotion {
    /// Get a position at the time `t`.
    fn position_at_time(&self, t: Real) -> Isometry<Real>;
}

impl RigidMotion for Isometry<Real> {
    fn position_at_time(&self, _: Real) -> Isometry<Real> {
        *self
    }
}

/// Interpolation between two isometries using LERP for the translation part and SLERP for the rotation.
pub struct InterpolatedRigidMotion {
    /// The transformation at `t = 0.0`.
    pub start: Isometry<Real>,
    /// The transformation at `t = 1.0`.
    pub end: Isometry<Real>,
}

impl InterpolatedRigidMotion {
    /// Initialize a lerp-slerp motion with the given start and end transformations.
    ///
    /// The `start` is the transformation at the time `t = 0.0` and `end` is the transformation at
    /// the time `t = 1.0`.
    pub fn new(start: Isometry<Real>, end: Isometry<Real>) -> Self {
        InterpolatedRigidMotion { start, end }
    }
}

impl RigidMotion for InterpolatedRigidMotion {
    fn position_at_time(&self, t: Real) -> Isometry<Real> {
        self.start.lerp_slerp(&self.end, t)
    }
}

/// A linear motion from a starting isometry traveling at constant translational velocity.
pub struct ConstantLinearVelocityRigidMotion {
    /// The time at which this parametrization begins. Can be negative.
    pub t0: Real,
    /// The starting isometry at `t = self.t0`.
    pub start: Isometry<Real>,
    /// The translational velocity of this motion.
    pub velocity: Vector<Real>,
}

impl ConstantLinearVelocityRigidMotion {
    /// Initialize a linear motion from a starting isometry and a translational velocity.
    pub fn new(t0: Real, start: Isometry<Real>, velocity: Vector<Real>) -> Self {
        ConstantLinearVelocityRigidMotion {
            t0,
            start,
            velocity,
        }
    }
}

impl RigidMotion for ConstantLinearVelocityRigidMotion {
    fn position_at_time(&self, t: Real) -> Isometry<Real> {
        Isometry::from_parts(
            (self.start.translation.vector + self.velocity * (t - self.t0)).into(),
            self.start.rotation,
        )
    }
}

/// A linear motion from a starting isometry traveling at constant translational and rotational velocity.
#[derive(Debug)]
pub struct ConstantVelocityRigidMotion {
    /// The time at which this parametrization begins. Can be negative.
    pub t0: Real,
    /// The starting isometry at `t = self.t0`.
    pub start: Isometry<Real>,
    /// The local-space point at which the rotational part of this motion is applied.
    pub local_center: Point<Real>,
    /// The translational velocity of this motion.
    pub linvel: Vector<Real>,
    /// The angular velocity of this motion.
    #[cfg(feature = "dim2")]
    pub angvel: Real,
    /// The angular velocity of this motion.
    #[cfg(feature = "dim3")]
    pub angvel: Vector<Real>,
}

impl ConstantVelocityRigidMotion {
    /// Initialize a motion from a starting isometry and linear and angular velocities.
    #[cfg(feature = "dim2")]
    pub fn new(
        t0: Real,
        start: Isometry<Real>,
        local_center: Point<Real>,
        linvel: Vector<Real>,
        angvel: Real,
    ) -> Self {
        ConstantVelocityRigidMotion {
            t0,
            start,
            local_center,
            linvel,
            angvel,
        }
    }

    /// Initialize a motion from a starting isometry and linear and angular velocities.
    #[cfg(feature = "dim3")]
    pub fn new(
        t0: Real,
        start: Isometry<Real>,
        local_center: Point<Real>,
        linvel: Vector<Real>,
        angvel: Vector<Real>,
    ) -> Self {
        ConstantVelocityRigidMotion {
            t0,
            start,
            local_center,
            linvel,
            angvel,
        }
    }
}

impl RigidMotion for ConstantVelocityRigidMotion {
    fn position_at_time(&self, t: Real) -> Isometry<Real> {
        let scaled_linvel = self.linvel * (t - self.t0);
        let scaled_angvel = self.angvel * (t - self.t0);

        let center = self.start.rotation * self.local_center.coords;
        let lhs = self.start.translation * Translation::from(center);
        let rhs = Translation::from(-center) * self.start.rotation;

        lhs * Isometry::new(scaled_linvel, scaled_angvel) * rhs
    }
}

/*
 * For composition.
 */

/// Trait for composing some rigid motions.
pub trait RigidMotionComposition: RigidMotion {
    /// Inverse the transforms returned by the rigid-motion `self`.
    fn inverse(&self) -> Inverse<Self> {
        Inverse { motion: self }
    }

    /// Returns the motion that returns `self.position_at_time(t).inverse() * rhs.position_at_time(t)`
    /// for all `t`.
    fn inv_mul<'a>(&'a self, rhs: &'a dyn RigidMotion) -> InvMul<'a, Self, dyn RigidMotion + 'a> {
        InvMul {
            motion1: self,
            motion2: rhs,
        }
    }

    /// Prepend a translation to the rigid motion `self`.
    fn prepend_translation(&self, translation: Vector<Real>) -> PrependTranslation<Self> {
        PrependTranslation {
            motion: self,
            translation,
        }
    }

    /// Prepend a translation to the rigid motion `self`.
    fn append_translation(&self, translation: Vector<Real>) -> AppendTranslation<Self> {
        AppendTranslation {
            motion: self,
            translation,
        }
    }

    /// Prepend a transformation to the rigid motion `self`.
    fn prepend_transformation(
        &self,
        transformation: Isometry<Real>,
    ) -> PrependTransformation<Self> {
        PrependTransformation {
            motion: self,
            transformation,
        }
    }

    /// Append a transformation to the rigid motion `self`.
    fn append_transformation(&self, transformation: Isometry<Real>) -> AppendTransformation<Self> {
        AppendTransformation {
            motion: self,
            transformation,
        }
    }
}

impl<M: ?Sized + RigidMotion> RigidMotionComposition for M {}

/// The result of prepending a translation to a rigid motion.
pub struct Inverse<'a, M: ?Sized> {
    motion: &'a M,
}

impl<'a, M: ?Sized + RigidMotion> RigidMotion for Inverse<'a, M> {
    fn position_at_time(&self, t: Real) -> Isometry<Real> {
        self.motion.position_at_time(t).inverse()
    }
}

/// The result of inverting `M1.inverse() * M2`.
pub struct InvMul<'a, M1: ?Sized, M2: ?Sized> {
    motion1: &'a M1,
    motion2: &'a M2,
}

impl<'a, M1: ?Sized + RigidMotion, M2: ?Sized + RigidMotion> RigidMotion for InvMul<'a, M1, M2> {
    fn position_at_time(&self, t: Real) -> Isometry<Real> {
        let m1 = self.motion1.position_at_time(t);
        let m2 = self.motion2.position_at_time(t);
        m1.inv_mul(&m2)
    }
}

/// The result of prepending a translation to a rigid motion.
pub struct PrependTranslation<'a, M: ?Sized> {
    motion: &'a M,
    translation: Vector<Real>,
}

impl<'a, M: ?Sized + RigidMotion> RigidMotion for PrependTranslation<'a, M> {
    fn position_at_time(&self, t: Real) -> Isometry<Real> {
        let m = self.motion.position_at_time(t);
        m * Translation::from(self.translation)
    }
}

/// The result of appending a translation to a rigid motion.
pub struct AppendTranslation<'a, M: ?Sized> {
    motion: &'a M,
    translation: Vector<Real>,
}

impl<'a, M: ?Sized + RigidMotion> RigidMotion for AppendTranslation<'a, M> {
    fn position_at_time(&self, t: Real) -> Isometry<Real> {
        let m = self.motion.position_at_time(t);
        Translation::from(self.translation) * m
    }
}

/// The result of prepending an isometric transformation to a rigid motion.
pub struct PrependTransformation<'a, M: ?Sized> {
    motion: &'a M,
    transformation: Isometry<Real>,
}

impl<'a, M: ?Sized + RigidMotion> RigidMotion for PrependTransformation<'a, M> {
    fn position_at_time(&self, t: Real) -> Isometry<Real> {
        let m = self.motion.position_at_time(t);
        m * self.transformation
    }
}

/// The result of prepending an isometric transformation to a rigid motion.
pub struct AppendTransformation<'a, M: ?Sized> {
    motion: &'a M,
    transformation: Isometry<Real>,
}

impl<'a, M: ?Sized + RigidMotion> RigidMotion for AppendTransformation<'a, M> {
    fn position_at_time(&self, t: Real) -> Isometry<Real> {
        let m = self.motion.position_at_time(t);
        self.transformation * m
    }
}
