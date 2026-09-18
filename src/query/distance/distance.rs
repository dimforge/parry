use crate::shape::SubShapeId;

/// The distance between two shapes, and the sub-shapes it was measured between.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ShapeDistance {
    /// The separation between the two shapes; `0.0` when they touch or overlap.
    pub distance: Real,
    /// The sub-shape of the first shape the distance was measured to.
    ///
    /// Always `0` for a shape with no sub-shapes.
    pub subshape1: SubShapeId,
    /// The sub-shape of the second shape the distance was measured to.
    ///
    /// Always `0` for a shape with no sub-shapes.
    pub subshape2: SubShapeId,
}

impl ShapeDistance {
    /// A distance measured between shapes with no sub-shapes to distinguish.
    pub fn new(distance: Real) -> Self {
        Self {
            distance,
            subshape1: 0,
            subshape2: 0,
        }
    }

    /// Sets the sub-shapes this distance was measured between.
    pub fn with_subshapes(mut self, subshape1: SubShapeId, subshape2: SubShapeId) -> Self {
        self.subshape1 = subshape1;
        self.subshape2 = subshape2;
        self
    }

    /// Swaps the roles of the two shapes.
    pub fn swapped(mut self) -> Self {
        core::mem::swap(&mut self.subshape1, &mut self.subshape2);
        self
    }
}

#[cfg(feature = "alloc")]
impl crate::partitioning::BvhLeafCost for ShapeDistance {
    #[inline]
    fn cost(&self) -> Real {
        self.distance
    }
}

impl From<Real> for ShapeDistance {
    fn from(distance: Real) -> Self {
        Self::new(distance)
    }
}

use crate::math::{Pose, Real};

use crate::query::{DefaultQueryDispatcher, QueryDispatcher, Unsupported};
use crate::shape::Shape;

/// Computes the minimum distance separating two shapes.
///
/// This is one of the most fundamental geometric queries in collision detection.
/// It calculates the shortest distance between any two points on the surfaces of
/// the two shapes.
///
/// # Behavior
///
/// - Returns the **shortest distance** between the two shape surfaces
/// - Returns `0.0` if the shapes are **touching** (surfaces just make contact)
/// - Returns `0.0` if the shapes are **penetrating** (overlapping)
/// - Always returns a **non-negative** value
///
/// # Arguments
///
/// * `pos1` - Position and orientation of the first shape in world space
/// * `g1` - The first shape (can be any shape implementing the `Shape` trait)
/// * `pos2` - Position and orientation of the second shape in world space
/// * `g2` - The second shape (can be any shape implementing the `Shape` trait)
///
/// # Returns
///
/// * `Ok(distance)` - The minimum distance between the shapes
/// * `Err(Unsupported)` - If this shape pair combination is not supported
///
/// # Performance
///
/// Performance varies by shape type:
/// - **Ball-Ball**: Very fast (analytical solution)
/// - **Cuboid-Cuboid**: Fast (SAT-based)
/// - **Convex-Convex**: Moderate (GJK algorithm)
/// - **Concave shapes**: Slower (requires BVH traversal)
///
/// # Example
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::query::distance;
/// use parry3d::shape::Ball;
/// use parry3d::math::Pose;
///
/// // Create two balls
/// let ball1 = Ball::new(1.0);
/// let ball2 = Ball::new(2.0);
///
/// // Position them 10 units apart along the x-axis
/// let pos1 = Pose::translation(0.0, 0.0, 0.0);
/// let pos2 = Pose::translation(10.0, 0.0, 0.0);
///
/// // Compute distance
/// let dist = distance(&pos1, &ball1, &pos2, &ball2).unwrap().distance;
///
/// // Distance = 10.0 (separation) - 1.0 (radius1) - 2.0 (radius2) = 7.0
/// assert_eq!(dist, 7.0);
/// # }
/// ```
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::query::distance;
/// use parry3d::shape::Cuboid;
/// use parry3d::math::{Pose, Vector};
///
/// // Create two boxes
/// let box1 = Cuboid::new(Vector::new(1.0, 1.0, 1.0));
/// let box2 = Cuboid::new(Vector::new(0.5, 0.5, 0.5));
///
/// // Position them so they're touching
/// let pos1 = Pose::translation(0.0, 0.0, 0.0);
/// let pos2 = Pose::translation(1.5, 0.0, 0.0); // Edge to edge
///
/// let dist = distance(&pos1, &box1, &pos2, &box2).unwrap().distance;
///
/// // They're touching, so distance is 0.0
/// assert_eq!(dist, 0.0);
/// # }
/// ```
///
/// # See Also
///
/// - [`closest_points`](crate::query::closest_points()) - For finding the actual closest points
/// - [`contact`](crate::query::contact()) - For penetration depth when overlapping
/// - [`intersection_test`](crate::query::intersection_test()) - For boolean overlap test
pub fn distance(
    pos1: &Pose,
    g1: &dyn Shape,
    pos2: &Pose,
    g2: &dyn Shape,
) -> Result<ShapeDistance, Unsupported> {
    let pos12 = pos1.inv_mul(pos2);
    DefaultQueryDispatcher.distance(&pos12, g1, g2)
}
