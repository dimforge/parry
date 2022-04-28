use crate::bounding_volume::AABB;
use crate::math::Real;
use crate::query::SplitResult;

impl AABB {
    /// Splits this AABB along the given canonical axis.
    ///
    /// This will split the AABB by a plane with a normal with it’s `axis`-th component set to 1.
    /// The splitting plane is shifted wrt. the origin by the `bias` (i.e. it passes through the point
    /// equal to `normal * bias`).
    ///
    /// # Result
    /// Returns the result of the split. The first AABB returned is the piece lying on the negative
    /// half-space delimited by the splitting plane. The second AABB returned is the piece lying on the
    /// positive half-space delimited by the splitting plane.
    pub fn canonical_split(&self, axis: usize, bias: Real, epsilon: Real) -> SplitResult<Self> {
        if self.mins[axis] >= bias - epsilon {
            SplitResult::Positive
        } else if self.maxs[axis] <= bias + epsilon {
            SplitResult::Negative
        } else {
            let mut left = *self;
            let mut right = *self;
            left.maxs[axis] = bias;
            right.mins[axis] = bias;
            SplitResult::Pair(left, right)
        }
    }
}
