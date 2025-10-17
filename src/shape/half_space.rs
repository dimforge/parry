//! Support mapping based HalfSpace shape.
use crate::math::{Real, Vector};
use na::Unit;

#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

/// A half-space delimited by an infinite plane.
///
/// # What is a HalfSpace?
///
/// A half-space represents an infinite region of space on one side of a plane. It divides
/// space into two regions:
/// - The "inside" region (where the normal vector points away from)
/// - The "outside" region (where the normal vector points toward)
///
/// The plane itself passes through the origin of the shape's coordinate system and is defined
/// by its outward normal vector. All points in the direction opposite to the normal are
/// considered "inside" the half-space.
///
/// # When to Use HalfSpace
///
/// Half-spaces are useful for representing:
/// - **Ground planes**: A flat, infinite floor for collision detection
/// - **Walls**: Infinite vertical barriers
/// - **Bounding regions**: Constraining objects to one side of a plane
/// - **Clipping planes**: Cutting off geometry in one direction
///
/// Because half-spaces are infinite, they are very efficient for collision detection and
/// don't require complex shape representations.
///
/// # Coordinate System
///
/// The plane always passes through the origin `(0, 0)` in 2D or `(0, 0, 0)` in 3D of the
/// half-space's local coordinate system. To position the plane elsewhere in your world,
/// use an [`Isometry`](na::Isometry) transformation when performing queries.
///
/// # Examples
///
/// ## Creating a Ground Plane (3D)
///
/// ```
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::shape::HalfSpace;
/// use parry3d::na::{Vector3, Unit};
///
/// // Create a horizontal ground plane with normal pointing up (positive Y-axis)
/// let ground = HalfSpace::new(Unit::new_normalize(Vector3::y()));
///
/// // The ground plane is at Y = 0 in local coordinates
/// // Everything below (negative Y) is "inside" the half-space
/// # }
/// ```
///
/// ## Vertical Wall (2D)
///
/// ```
/// # #[cfg(all(feature = "dim2", feature = "f32"))] {
/// use parry2d::shape::HalfSpace;
/// use parry2d::na::{Vector2, Unit};
///
/// // Create a vertical wall with normal pointing right (positive X-axis)
/// let wall = HalfSpace::new(Unit::new_normalize(Vector2::x()));
///
/// // The wall is at X = 0 in local coordinates
/// // Everything to the left (negative X) is "inside" the half-space
/// # }
/// ```
///
/// ## Collision Detection with a Ball (3D)
///
/// ```
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::shape::{HalfSpace, Ball};
/// use parry3d::query;
/// use parry3d::na::{Isometry3, Vector3, Unit};
///
/// // Create a ground plane at Y = 0, normal pointing up
/// let ground = HalfSpace::new(Unit::new_normalize(Vector3::y()));
/// let ground_pos = Isometry3::identity();
///
/// // Create a ball with radius 1.0 at position (0, 0.5, 0)
/// // The ball is resting on the ground, just touching it
/// let ball = Ball::new(1.0);
/// let ball_pos = Isometry3::translation(0.0, 0.5, 0.0);
///
/// // Check if they're in contact (with a small prediction distance)
/// let contact = query::contact(
///     &ground_pos,
///     &ground,
///     &ball_pos,
///     &ball,
///     0.1
/// );
///
/// assert!(contact.unwrap().is_some());
/// # }
/// ```
///
/// ## Positioned Ground Plane (3D)
///
/// ```
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::shape::HalfSpace;
/// use parry3d::query::{PointQuery};
/// use parry3d::na::{Isometry3, Vector3, Point3, Unit};
///
/// // Create a ground plane with normal pointing up
/// let ground = HalfSpace::new(Unit::new_normalize(Vector3::y()));
///
/// // Position the plane at Y = 5.0 using an isometry
/// let ground_pos = Isometry3::translation(0.0, 5.0, 0.0);
///
/// // Check if a point is below the ground (inside the half-space)
/// let point = Point3::new(0.0, 3.0, 0.0); // Point at Y = 3.0 (below the plane)
///
/// // Project the point onto the ground plane
/// let proj = ground.project_point(&ground_pos, &point, true);
///
/// // The point is below the ground (inside the half-space)
/// assert!(proj.is_inside);
/// # }
/// ```
///
/// ## Tilted Plane (3D)
///
/// ```
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::shape::HalfSpace;
/// use parry3d::na::{Vector3, Unit};
///
/// // Create a plane tilted at 45 degrees
/// // Normal points up and to the right
/// let normal = Vector3::new(1.0, 1.0, 0.0);
/// let tilted_plane = HalfSpace::new(Unit::new_normalize(normal));
///
/// // This plane passes through the origin and divides space diagonally
/// # }
/// ```
#[derive(PartialEq, Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[repr(C)]
pub struct HalfSpace {
    /// The halfspace planar boundary's outward normal.
    ///
    /// This unit vector points in the direction considered "outside" the half-space.
    /// All points in the direction opposite to this normal (when measured from the
    /// plane at the origin) are considered "inside" the half-space.
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::shape::HalfSpace;
    /// use parry3d::na::{Vector3, Unit};
    ///
    /// let ground = HalfSpace::new(Unit::new_normalize(Vector3::y()));
    ///
    /// // The normal points up (positive Y direction)
    /// assert_eq!(*ground.normal, Vector3::y());
    /// # }
    /// ```
    pub normal: Unit<Vector<Real>>,
}

impl HalfSpace {
    /// Builds a new half-space from its outward normal vector.
    ///
    /// The plane defining the half-space passes through the origin of the local coordinate
    /// system and is perpendicular to the given normal vector. The normal points toward
    /// the "outside" region, while the opposite direction is considered "inside."
    ///
    /// # Parameters
    ///
    /// * `normal` - A unit vector defining the plane's outward normal direction. This must
    ///   be a normalized vector (use [`Unit::new_normalize`](na::Unit::new_normalize) to
    ///   create one from any vector).
    ///
    /// # Examples
    ///
    /// ## Creating a Horizontal Ground Plane (3D)
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::shape::HalfSpace;
    /// use parry3d::na::{Vector3, Unit};
    ///
    /// // Ground plane with normal pointing up
    /// let ground = HalfSpace::new(Unit::new_normalize(Vector3::y()));
    /// # }
    /// ```
    ///
    /// ## Creating a Vertical Wall (2D)
    ///
    /// ```
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// use parry2d::shape::HalfSpace;
    /// use parry2d::na::{Vector2, Unit};
    ///
    /// // Wall with normal pointing to the right
    /// let wall = HalfSpace::new(Unit::new_normalize(Vector2::x()));
    /// # }
    /// ```
    ///
    /// ## Custom Normal Direction (3D)
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::shape::HalfSpace;
    /// use parry3d::na::{Vector3, Unit};
    ///
    /// // Plane with normal at 45-degree angle
    /// let custom_normal = Vector3::new(1.0, 1.0, 0.0);
    /// let plane = HalfSpace::new(Unit::new_normalize(custom_normal));
    ///
    /// // Verify the normal is normalized
    /// assert!((plane.normal.norm() - 1.0).abs() < 1e-5);
    /// # }
    /// ```
    #[inline]
    pub fn new(normal: Unit<Vector<Real>>) -> HalfSpace {
        HalfSpace { normal }
    }

    /// Computes a scaled version of this half-space.
    ///
    /// Scaling a half-space applies non-uniform scaling to its normal vector. This is useful
    /// when transforming shapes in a scaled coordinate system. The resulting normal is
    /// re-normalized to maintain the half-space's validity.
    ///
    /// # Parameters
    ///
    /// * `scale` - A vector containing the scaling factors for each axis. For example,
    ///   `Vector3::new(2.0, 1.0, 1.0)` doubles the X-axis scaling.
    ///
    /// # Returns
    ///
    /// * `Some(HalfSpace)` - The scaled half-space with the transformed normal
    /// * `None` - If the scaled normal becomes zero (degenerate case), meaning the
    ///   half-space cannot be represented after scaling
    ///
    /// # When This Returns None
    ///
    /// The method returns `None` when any component of the normal becomes zero after
    /// scaling AND that component was the only non-zero component. For example:
    /// - A horizontal plane (normal = `[0, 1, 0]`) scaled by `[1, 0, 1]` → `None`
    /// - A diagonal plane (normal = `[0.7, 0.7, 0]`) scaled by `[1, 0, 1]` → `Some(...)`
    ///
    /// # Examples
    ///
    /// ## Uniform Scaling (3D)
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::shape::HalfSpace;
    /// use parry3d::na::{Vector3, Unit};
    ///
    /// let ground = HalfSpace::new(Unit::new_normalize(Vector3::y()));
    ///
    /// // Uniform scaling doesn't change the normal direction
    /// let scaled = ground.scaled(&Vector3::new(2.0, 2.0, 2.0)).unwrap();
    /// assert_eq!(*scaled.normal, Vector3::y());
    /// # }
    /// ```
    ///
    /// ## Non-Uniform Scaling (3D)
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::shape::HalfSpace;
    /// use parry3d::na::{Vector3, Unit};
    ///
    /// // Diagonal plane
    /// let plane = HalfSpace::new(
    ///     Unit::new_normalize(Vector3::new(1.0, 1.0, 0.0))
    /// );
    ///
    /// // Scale X-axis by 2.0, Y-axis stays 1.0
    /// let scaled = plane.scaled(&Vector3::new(2.0, 1.0, 1.0)).unwrap();
    ///
    /// // The normal changes direction due to non-uniform scaling
    /// // It's no longer at 45 degrees
    /// assert!(scaled.normal.x != scaled.normal.y);
    /// # }
    /// ```
    ///
    /// ## Degenerate Case (3D)
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::shape::HalfSpace;
    /// use parry3d::na::{Vector3, Unit};
    ///
    /// // Horizontal ground plane
    /// let ground = HalfSpace::new(Unit::new_normalize(Vector3::y()));
    ///
    /// // Scaling Y to zero makes the normal degenerate
    /// let scaled = ground.scaled(&Vector3::new(1.0, 0.0, 1.0));
    /// assert!(scaled.is_none()); // Returns None because normal becomes zero
    /// # }
    /// ```
    ///
    /// ## Practical Use Case (2D)
    ///
    /// ```
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// use parry2d::shape::HalfSpace;
    /// use parry2d::na::{Vector2, Unit};
    ///
    /// // Create a wall in a 2D platformer
    /// let wall = HalfSpace::new(Unit::new_normalize(Vector2::x()));
    ///
    /// // Apply level scaling (e.g., for pixel-perfect rendering)
    /// let pixel_scale = Vector2::new(16.0, 16.0);
    /// if let Some(scaled_wall) = wall.scaled(&pixel_scale) {
    ///     // Use the scaled wall for collision detection
    /// }
    /// # }
    /// ```
    pub fn scaled(self, scale: &Vector<Real>) -> Option<Self> {
        Unit::try_new(self.normal.component_mul(scale), 0.0).map(|normal| Self { normal })
    }
}
