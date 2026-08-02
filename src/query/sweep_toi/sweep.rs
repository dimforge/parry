use crate::math::{Pose, Real, Rotation, Vector};

/// Rotates a vector by a rotation.
#[inline]
pub(crate) fn rotate_vec(q: &Rotation, v: Vector) -> Vector {
    #[cfg(feature = "dim2")]
    {
        q.transform_vector(v)
    }
    #[cfg(feature = "dim3")]
    {
        *q * v
    }
}

/// Rotates a vector by the inverse of a rotation.
#[inline]
pub(crate) fn inv_rotate_vec(q: &Rotation, v: Vector) -> Vector {
    #[cfg(feature = "dim2")]
    {
        q.inverse_transform_vector(v)
    }
    #[cfg(feature = "dim3")]
    {
        q.inverse() * v
    }
}

/// Normalized linear interpolation between two rotations (3D: shortest arc).
#[inline]
pub(crate) fn nlerp(q1: &Rotation, q2: &Rotation, t: Real) -> Rotation {
    #[cfg(feature = "dim2")]
    {
        q1.lerp(*q2, t).normalize()
    }
    #[cfg(feature = "dim3")]
    {
        let q1 = if q1.dot(*q2) < 0.0 { -*q1 } else { *q1 };
        (q1 * (1.0 - t) + *q2 * t).normalize()
    }
}

/// Describes the motion of a rigid body over a timestep as linear interpolation between two
/// endpoint poses: the center of mass moves on a straight line while the rotation is
/// interpolated with a normalized lerp (nlerp).
///
/// This is a common motion model for continuous collision detection. It is exact
/// at both endpoints and a good approximation in between as long as the rotation delta stays
/// below ~45°.
#[derive(Copy, Clone, Debug, PartialEq)]
#[cfg_attr(
    feature = "serde-serialize",
    derive(serde::Serialize, serde::Deserialize)
)]
pub struct Sweep {
    /// The center of mass expressed in the shape’s local frame.
    pub local_center: Vector,
    /// The world-space center of mass at the start of the sweep.
    pub c1: Vector,
    /// The world-space center of mass at the end of the sweep.
    pub c2: Vector,
    /// The rotation at the start of the sweep.
    pub q1: Rotation,
    /// The rotation at the end of the sweep.
    pub q2: Rotation,
}

impl Sweep {
    /// Builds a sweep from the start and end poses of a shape’s local frame.
    pub fn from_poses(start: &Pose, end: &Pose, local_center: Vector) -> Self {
        Self {
            local_center,
            c1: start.transform_point(local_center),
            c2: end.transform_point(local_center),
            q1: start.rotation,
            q2: end.rotation,
        }
    }

    /// A degenerate sweep holding the shape stationary at the given pose.
    pub fn constant(pose: &Pose, local_center: Vector) -> Self {
        Self::from_poses(pose, pose, local_center)
    }

    /// The pose of the shape’s local frame at time `t ∈ [0, 1]`.
    ///
    /// The center of mass is lerped, the rotation is nlerped, and the local-frame origin is
    /// recovered by un-shifting the local center.
    pub fn transform_at(&self, t: Real) -> Pose {
        let q = nlerp(&self.q1, &self.q2, t);
        let p = self.c1.lerp(self.c2, t) - rotate_vec(&q, self.local_center);
        Pose::from_parts(p, q)
    }

    /// The pose of the shape’s local frame at the end of the sweep (`t = 1`), computed exactly.
    pub fn final_transform(&self) -> Pose {
        let p = self.c2 - rotate_vec(&self.q2, self.local_center);
        Pose::from_parts(p, self.q2)
    }

    /// Translates the entire sweep by `-origin`.
    ///
    /// Used to re-center the time-of-impact computation for better floating-point accuracy.
    pub fn shifted(&self, origin: Vector) -> Self {
        Self {
            local_center: self.local_center,
            c1: self.c1 - origin,
            c2: self.c2 - origin,
            q1: self.q1,
            q2: self.q2,
        }
    }
}
