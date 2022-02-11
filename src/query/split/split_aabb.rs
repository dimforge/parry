use crate::bounding_volume::AABB;
use crate::math::Real;
use crate::query::{CanonicalSplit, SplitResult};

impl CanonicalSplit for AABB {
    fn canonical_split(&self, axis: usize, bias: Real, epsilon: Real) -> SplitResult<Self> {
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
