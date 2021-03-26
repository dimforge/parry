use crate::math::{Isometry, Point, Real, Translation, Vector};

/// A nonlinear motion from a starting isometry traveling at constant translational and rotational velocity.
#[derive(Debug, Copy, Clone)]
pub struct NonlinearRigidMotion {
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

impl NonlinearRigidMotion {
    /// Initialize a motion from a starting isometry and linear and angular velocities.
    #[cfg(feature = "dim2")]
    pub fn new(
        t0: Real,
        start: Isometry<Real>,
        local_center: Point<Real>,
        linvel: Vector<Real>,
        angvel: Real,
    ) -> Self {
        NonlinearRigidMotion {
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
        NonlinearRigidMotion {
            t0,
            start,
            local_center,
            linvel,
            angvel,
        }
    }

    /// Create a `RigidMotion` that always returns the identity matrix.
    pub fn identity() -> Self {
        Self::constant_position(Isometry::identity())
    }

    /// Create a `RigidMotion` that always return `pos`.
    pub fn constant_position(pos: Isometry<Real>) -> Self {
        Self {
            t0: 0.0,
            start: pos,
            linvel: na::zero(),
            angvel: na::zero(),
            local_center: pos * Point::origin(),
        }
    }

    /// Freezes this motion at the time `t`.
    ///
    /// After calling this, any further calls to `self.position_at_time`
    /// will always return `self.position_at_time(t)` (where `t` is the value given
    /// to this method). This sets the linear velocity and angular velocity
    /// of `self` to zero.
    pub fn freeze(&mut self, t: Real) {
        self.start = self.position_at_time(t);
        self.linvel = na::zero();
        self.angvel = na::zero();
    }

    /// Appends a constant translation to this rigid-motion.
    #[must_use]
    pub fn append_translation(&self, tra: Vector<Real>) -> Self {
        let mut result = self.clone();
        result.start = Translation::from(tra) * result.start;
        result
    }

    /// Prepends a constant translation to this rigid-motion.
    #[must_use]
    pub fn prepend_translation(&self, tra: Vector<Real>) -> Self {
        let mut result = self.clone();
        result.start = result.start * Translation::from(tra);
        result
    }

    /// Appends a constant isometry to this rigid-motion.
    #[must_use]
    pub fn append(&self, iso: Isometry<Real>) -> Self {
        let mut result = self.clone();
        result.start = iso * result.start;
        result
    }

    /// Prepends a constant translation to this rigid-motion.
    #[must_use]
    pub fn prepend(&self, iso: Isometry<Real>) -> Self {
        let mut result = self.clone();
        result.start = result.start * iso;
        result
    }

    pub fn position_at_time(&self, t: Real) -> Isometry<Real> {
        let dt = t - self.t0;
        let center = self.start * self.local_center;
        let shift = Translation::from(center.coords);
        (shift * Isometry::new(self.linvel * dt, self.angvel * dt)) * (shift.inverse() * self.start)
    }
}
