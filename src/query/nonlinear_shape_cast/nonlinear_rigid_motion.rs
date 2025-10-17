use crate::math::{Isometry, Point, Real, Translation, Vector};

/// Describes the complete motion of a rigid body with both translation and rotation.
///
/// This type represents the **6 degrees of freedom** motion of a rigid body:
/// - **3 translational** (linear velocity in x, y, z)
/// - **3 rotational** (angular velocity around x, y, z axes)
///
/// The motion is assumed to have **constant velocities** over the time interval, meaning:
/// - Linear velocity doesn't change (no acceleration)
/// - Angular velocity doesn't change (no angular acceleration)
///
/// This is used with [`cast_shapes_nonlinear`](crate::query::cast_shapes_nonlinear) to perform
/// continuous collision detection for rotating objects.
///
/// # Physics Model
///
/// At any time `t`, the object's position is computed as:
/// 1. Rotate around `local_center` by `angvel * t`
/// 2. Translate by `linvel * t`
/// 3. Apply to the starting position `start`
///
/// This creates a **helical trajectory** (螺旋軌跡) when both velocities are non-zero:
/// - Pure translation: straight line
/// - Pure rotation: circular arc
/// - Both: helix/spiral path
///
/// # Fields
///
/// * `start` - The initial position and orientation at time `t = 0`
/// * `local_center` - The point (in the shape's local coordinate system) around which
///   rotation occurs. For most cases, use the center of mass or `Point::origin()`.
/// * `linvel` - Linear velocity vector (units per second). Direction is the direction
///   of motion, magnitude is speed.
/// * `angvel` - Angular velocity:
///   - **2D**: Scalar rotation rate in radians/second (positive = counter-clockwise)
///   - **3D**: Axis-angle vector (direction = rotation axis, magnitude = rotation rate in radians/second)
///
/// # Example: Simple Translation
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::query::NonlinearRigidMotion;
/// use nalgebra::{Isometry3, Point3, Vector3};
///
/// // Object moving right at 5 units/second, no rotation
/// let motion = NonlinearRigidMotion::new(
///     Isometry3::translation(0.0, 0.0, 0.0), // start at origin
///     Point3::origin(),                      // rotation center (irrelevant here)
///     Vector3::new(5.0, 0.0, 0.0),          // moving right
///     Vector3::zeros(),                      // not rotating
/// );
///
/// // At t=2.0 seconds, object has moved 10 units right
/// let pos_at_2 = motion.position_at_time(2.0);
/// assert_eq!(pos_at_2.translation.vector.x, 10.0);
/// # }
/// ```
///
/// # Example: Pure Rotation (3D)
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::query::NonlinearRigidMotion;
/// use nalgebra::{Isometry3, Point3, Vector3};
/// use std::f32::consts::PI;
///
/// // Object spinning around Y axis, no translation
/// let motion = NonlinearRigidMotion::new(
///     Isometry3::translation(0.0, 0.0, 0.0),
///     Point3::origin(),                // rotate around origin
///     Vector3::zeros(),                // not translating
///     Vector3::new(0.0, PI, 0.0),      // rotating around Y at π rad/s (180°/s)
/// );
///
/// // At t=1.0 second, object has rotated 180 degrees
/// let pos_at_1 = motion.position_at_time(1.0);
/// // Object is now rotated 180° around Y axis
/// # }
/// ```
///
/// # Example: Combined Translation and Rotation (Helical Path)
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::query::NonlinearRigidMotion;
/// use nalgebra::{Isometry3, Point3, Vector3};
///
/// // Spinning projectile moving forward
/// let motion = NonlinearRigidMotion::new(
///     Isometry3::translation(0.0, 0.0, 0.0),
///     Point3::origin(),
///     Vector3::new(10.0, 0.0, 0.0),     // moving forward at 10 units/s
///     Vector3::new(20.0, 0.0, 0.0),     // spinning around its movement axis
/// );
///
/// // The object traces a helical path (like a bullet with rifling)
/// let pos_at_half = motion.position_at_time(0.5);
/// // Moved 5 units forward AND rotated 10 radians
/// # }
/// ```
///
/// # Example: Rotation Around Off-Center Point
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::query::NonlinearRigidMotion;
/// use nalgebra::{Isometry3, Point3, Vector3};
///
/// // Object rotating around a point that's NOT its center
/// // Useful for: swinging weapons, rotating around pivot point, etc.
/// let motion = NonlinearRigidMotion::new(
///     Isometry3::translation(5.0, 0.0, 0.0), // object is at x=5
///     Point3::new(-5.0, 0.0, 0.0),           // rotate around x=0 (in local space)
///     Vector3::zeros(),
///     Vector3::new(0.0, 1.0, 0.0),           // rotate around Y axis at 1 rad/s
/// );
///
/// // The object orbits in a circle around the world origin
/// // Like a hammer being swung around
/// # }
/// ```
///
/// # Common Use Cases
///
/// 1. **Spinning Projectiles**: Bullets, thrown objects with spin
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// # use parry3d::query::NonlinearRigidMotion;
/// # use nalgebra::{Isometry3, Point3, Vector3};
/// let bullet = NonlinearRigidMotion::new(
///     Isometry3::translation(0.0, 1.5, 0.0),
///     Point3::origin(),
///     Vector3::new(100.0, -2.0, 0.0),  // fast forward, slight drop
///     Vector3::new(50.0, 0.0, 0.0),    // high spin rate
/// );
/// # }
/// ```
///
/// 2. **Tumbling Debris**: Objects affected by explosion or impact
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// # use parry3d::query::NonlinearRigidMotion;
/// # use nalgebra::{Isometry3, Point3, Vector3};
/// let debris = NonlinearRigidMotion::new(
///     Isometry3::translation(0.0, 2.0, 0.0),
///     Point3::origin(),
///     Vector3::new(3.0, 5.0, -2.0),    // chaotic velocity
///     Vector3::new(2.0, -3.0, 1.5),    // chaotic rotation
/// );
/// # }
/// ```
///
/// 3. **Rotating Machinery**: Blades, gears, rotating parts
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// # use parry3d::query::NonlinearRigidMotion;
/// # use nalgebra::{Isometry3, Point3, Vector3};
/// let blade = NonlinearRigidMotion::new(
///     Isometry3::translation(0.0, 1.0, 0.0),
///     Point3::origin(),           // spin around center
///     Vector3::zeros(),           // blade doesn't translate
///     Vector3::new(0.0, 10.0, 0.0), // fast rotation
/// );
/// # }
/// ```
///
/// 4. **Stationary Objects**: Use `constant_position()` for non-moving obstacles
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// # use parry3d::query::NonlinearRigidMotion;
/// # use nalgebra::Isometry3;
/// let wall = NonlinearRigidMotion::constant_position(
///     Isometry3::translation(10.0, 0.0, 0.0)
/// );
/// # }
/// ```
///
/// # 2D vs 3D Angular Velocity
///
/// The representation of `angvel` differs between 2D and 3D:
///
/// **2D (scalar)**:
/// - Positive value: counter-clockwise rotation
/// - Negative value: clockwise rotation
/// - Magnitude: rotation rate in radians/second
///
/// ```rust
/// # #[cfg(all(feature = "dim2", feature = "f32"))] {
/// use parry2d::query::NonlinearRigidMotion;
/// use nalgebra::{Isometry2, Point2, Vector2};
///
/// let motion = NonlinearRigidMotion::new(
///     Isometry2::translation(0.0, 0.0),
///     Point2::origin(),
///     Vector2::zeros(),
///     3.14,  // rotating counter-clockwise at π rad/s
/// );
/// # }
/// ```
///
/// **3D (axis-angle vector)**:
/// - Direction: axis of rotation (right-hand rule)
/// - Magnitude: rotation rate in radians/second
/// - Zero vector: no rotation
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::query::NonlinearRigidMotion;
/// use nalgebra::{Isometry3, Point3, Vector3};
///
/// let motion = NonlinearRigidMotion::new(
///     Isometry3::translation(0.0, 0.0, 0.0),
///     Point3::origin(),
///     Vector3::zeros(),
///     Vector3::new(0.0, 3.14, 0.0),  // rotate around Y axis at π rad/s
/// );
/// # }
/// ```
///
/// # Important Notes
///
/// 1. **Local vs World Space**: The `local_center` must be in the **shape's local
///    coordinate system**, not world space. For a shape centered at origin,
///    use `Point::origin()`.
///
/// 2. **Angular Velocity Units**: Always use **radians per second**, not degrees!
///    - To convert: `degrees * (PI / 180.0) = radians`
///    - Example: 90°/s = `90.0 * (PI / 180.0)` ≈ 1.571 rad/s
///
/// 3. **Constant Velocities**: This assumes velocities don't change over time.
///    For physics simulations with acceleration/forces, you typically:
///    - Compute motion for a single small timestep (e.g., 1/60 second)
///    - Update velocities after each step
///    - Create new `NonlinearRigidMotion` for next timestep
///
/// 4. **Center of Mass**: For realistic physics, `local_center` should be the
///    shape's center of mass. For simple cases, `Point::origin()` often works.
///
/// # See Also
///
/// - [`cast_shapes_nonlinear`](crate::query::cast_shapes_nonlinear) - Uses this type for collision detection
/// - [`cast_shapes`](crate::query::cast_shapes) - Linear motion (no rotation)
/// - [`position_at_time`](Self::position_at_time) - Compute position at any time
#[derive(Debug, Copy, Clone)]
pub struct NonlinearRigidMotion {
    /// The starting isometry at `t = 0`.
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
    /// Creates a new rigid motion from a starting position and velocities.
    ///
    /// # Arguments
    ///
    /// * `start` - Initial position and orientation at time `t = 0`
    /// * `local_center` - Point (in local coordinates) around which rotation occurs
    /// * `linvel` - Linear velocity vector (units per second)
    /// * `angvel` - Angular velocity (2D: radians/sec scalar, 3D: axis-angle vector)
    ///
    /// # Example (2D)
    ///
    /// ```
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// use parry2d::query::NonlinearRigidMotion;
    /// use nalgebra::{Isometry2, Point2, Vector2};
    /// use std::f32::consts::PI;
    ///
    /// // Object moving right and rotating counter-clockwise
    /// let motion = NonlinearRigidMotion::new(
    ///     Isometry2::translation(0.0, 0.0),
    ///     Point2::origin(),
    ///     Vector2::new(5.0, 0.0),  // 5 units/sec to the right
    ///     PI,                       // π rad/sec counter-clockwise
    /// );
    /// # }
    /// ```
    ///
    /// # Example (3D)
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::query::NonlinearRigidMotion;
    /// use nalgebra::{Isometry3, Point3, Vector3};
    ///
    /// // Object moving forward and spinning around its movement axis
    /// let motion = NonlinearRigidMotion::new(
    ///     Isometry3::translation(0.0, 0.0, 0.0),
    ///     Point3::origin(),
    ///     Vector3::new(10.0, 0.0, 0.0),    // moving forward
    ///     Vector3::new(5.0, 0.0, 0.0),     // spinning around X axis
    /// );
    /// # }
    /// ```
    #[cfg(feature = "dim2")]
    pub fn new(
        start: Isometry<Real>,
        local_center: Point<Real>,
        linvel: Vector<Real>,
        angvel: Real,
    ) -> Self {
        NonlinearRigidMotion {
            start,
            local_center,
            linvel,
            angvel,
        }
    }

    /// Creates a new rigid motion from a starting position and velocities.
    ///
    /// # Arguments
    ///
    /// * `start` - Initial position and orientation at time `t = 0`
    /// * `local_center` - Point (in local coordinates) around which rotation occurs
    /// * `linvel` - Linear velocity vector (units per second)
    /// * `angvel` - Angular velocity as axis-angle vector (radians per second)
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::query::NonlinearRigidMotion;
    /// use nalgebra::{Isometry3, Point3, Vector3};
    ///
    /// // Object moving forward and spinning around its movement axis
    /// let motion = NonlinearRigidMotion::new(
    ///     Isometry3::translation(0.0, 0.0, 0.0),
    ///     Point3::origin(),
    ///     Vector3::new(10.0, 0.0, 0.0),    // moving forward
    ///     Vector3::new(5.0, 0.0, 0.0),     // spinning around X axis
    /// );
    /// # }
    /// ```
    #[cfg(feature = "dim3")]
    pub fn new(
        start: Isometry<Real>,
        local_center: Point<Real>,
        linvel: Vector<Real>,
        angvel: Vector<Real>,
    ) -> Self {
        NonlinearRigidMotion {
            start,
            local_center,
            linvel,
            angvel,
        }
    }

    /// Creates a stationary motion at the origin (identity transformation).
    ///
    /// Equivalent to `constant_position(Isometry::identity())`.
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::query::NonlinearRigidMotion;
    ///
    /// let stationary = NonlinearRigidMotion::identity();
    /// // Object stays at origin with no rotation, forever
    /// # }
    /// ```
    pub fn identity() -> Self {
        Self::constant_position(Isometry::identity())
    }

    /// Creates a motion that stays at a constant position (no translation, no rotation).
    ///
    /// Useful for representing static/stationary objects in collision queries.
    /// Both linear and angular velocities are set to zero.
    ///
    /// # Arguments
    ///
    /// * `pos` - The fixed position and orientation
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::query::NonlinearRigidMotion;
    /// use nalgebra::Isometry3;
    ///
    /// // A wall that never moves
    /// let wall_motion = NonlinearRigidMotion::constant_position(
    ///     Isometry3::translation(10.0, 0.0, 0.0)
    /// );
    ///
    /// // At any time t, position is always (10, 0, 0)
    /// let pos_at_5 = wall_motion.position_at_time(5.0);
    /// assert_eq!(pos_at_5.translation.vector.x, 10.0);
    /// # }
    /// ```
    pub fn constant_position(pos: Isometry<Real>) -> Self {
        Self {
            start: pos,
            linvel: na::zero(),
            angvel: na::zero(),
            local_center: Point::origin(),
        }
    }

    fn set_start(&mut self, new_start: Isometry<Real>) {
        // NOTE: we need to adjust the local_center so that the angular
        // velocity is still expressed wrt. the original center.
        self.local_center = new_start.inverse_transform_point(&(self.start * self.local_center));
        self.start = new_start;
    }

    /// Freezes this motion at a specific time, making it stationary from that point forward.
    ///
    /// After calling this method:
    /// - `self.start` is updated to the position at time `t`
    /// - Both velocities are set to zero
    /// - All future calls to `position_at_time()` return the same frozen position
    ///
    /// This is useful for "stopping" an object mid-motion, or capturing a specific
    /// moment in a motion trajectory.
    ///
    /// # Arguments
    ///
    /// * `t` - The time at which to freeze the motion
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::query::NonlinearRigidMotion;
    /// use nalgebra::{Isometry3, Point3, Vector3};
    ///
    /// let mut motion = NonlinearRigidMotion::new(
    ///     Isometry3::translation(0.0, 0.0, 0.0),
    ///     Point3::origin(),
    ///     Vector3::new(5.0, 0.0, 0.0),  // moving right
    ///     Vector3::zeros(),
    /// );
    ///
    /// // Freeze at t=2.0 (when object is at x=10)
    /// motion.freeze(2.0);
    ///
    /// // Now position is constant at x=10, regardless of time
    /// let pos_at_100 = motion.position_at_time(100.0);
    /// assert_eq!(pos_at_100.translation.vector.x, 10.0);
    /// # }
    /// ```
    pub fn freeze(&mut self, t: Real) {
        self.start = self.position_at_time(t);
        self.linvel = na::zero();
        self.angvel = na::zero();
    }

    /// Appends a constant translation to this rigid-motion.
    #[must_use]
    pub fn append_translation(&self, tra: Vector<Real>) -> Self {
        let mut result = *self;
        result.set_start(Translation::from(tra) * result.start);
        result
    }

    /// Prepends a constant translation to this rigid-motion.
    #[must_use]
    pub fn prepend_translation(&self, tra: Vector<Real>) -> Self {
        let mut result = *self;
        result.set_start(result.start * Translation::from(tra));
        result
    }

    /// Appends a constant isometry to this rigid-motion.
    #[must_use]
    pub fn append(&self, iso: Isometry<Real>) -> Self {
        let mut result = *self;
        result.set_start(iso * result.start);
        result
    }

    /// Prepends a constant translation to this rigid-motion.
    #[must_use]
    pub fn prepend(&self, iso: Isometry<Real>) -> Self {
        let mut result = *self;
        result.set_start(result.start * iso);
        result
    }

    /// Computes the position and orientation at a given time.
    ///
    /// Returns the full isometry (position + orientation) of the rigid body at time `t`,
    /// accounting for both translation and rotation from the starting position.
    ///
    /// The computation follows this sequence:
    /// 1. Rotate around `local_center` by angle `angvel * t`
    /// 2. Translate by displacement `linvel * t`
    /// 3. Apply to the starting position `start`
    ///
    /// # Arguments
    ///
    /// * `t` - Time value (typically in seconds, matching your velocity units)
    ///
    /// # Returns
    ///
    /// The complete transformation (position and orientation) at time `t`.
    ///
    /// # Example: Tracking a Moving Object
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::query::NonlinearRigidMotion;
    /// use nalgebra::{Isometry3, Point3, Vector3};
    ///
    /// let motion = NonlinearRigidMotion::new(
    ///     Isometry3::translation(0.0, 0.0, 0.0),
    ///     Point3::origin(),
    ///     Vector3::new(3.0, 0.0, 0.0),     // 3 units/sec to the right
    ///     Vector3::new(0.0, 1.0, 0.0),     // 1 radian/sec around Y axis
    /// );
    ///
    /// // Position at t=0
    /// let pos_0 = motion.position_at_time(0.0);
    /// assert_eq!(pos_0.translation.vector.x, 0.0);
    ///
    /// // Position at t=2.0 seconds
    /// let pos_2 = motion.position_at_time(2.0);
    /// // Object has moved 6 units to the right
    /// assert!((pos_2.translation.vector.x - 6.0).abs() < 0.01);
    /// // And rotated 2 radians around Y
    /// # }
    /// ```
    ///
    /// # Example: Animation Frame
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::query::NonlinearRigidMotion;
    /// use nalgebra::{Isometry3, Point3, Vector3};
    ///
    /// let motion = NonlinearRigidMotion::new(
    ///     Isometry3::translation(0.0, 5.0, 0.0),
    ///     Point3::origin(),
    ///     Vector3::new(0.0, -9.8, 0.0),    // falling (gravity)
    ///     Vector3::new(1.0, 2.0, 0.5),     // tumbling
    /// );
    ///
    /// // Render at 60 FPS
    /// let dt = 1.0 / 60.0;
    /// for frame in 0..60 {
    ///     let t = frame as f32 * dt;
    ///     let pos = motion.position_at_time(t);
    ///     // Use `pos` to render object at this frame
    /// }
    /// # }
    /// ```
    pub fn position_at_time(&self, t: Real) -> Isometry<Real> {
        let center = self.start * self.local_center;
        let shift = Translation::from(center.coords);
        (shift * Isometry::new(self.linvel * t, self.angvel * t)) * (shift.inverse() * self.start)
    }
}
