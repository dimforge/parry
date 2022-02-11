use crate::math::{Real, UnitVector, Vector};
use crate::query::{CanonicalSplit, Split, SplitResult};
use crate::shape::Segment;

impl CanonicalSplit for Segment {
    fn canonical_split(&self, axis: usize, bias: Real, epsilon: Real) -> SplitResult<Self> {
        // TODO optimize this.
        self.local_split(&Vector::ith_axis(axis), bias, epsilon)
    }
}

impl Split for Segment {
    fn local_split(
        &self,
        local_axis: &UnitVector<Real>,
        bias: Real,
        epsilon: Real,
    ) -> SplitResult<Self> {
        let dir = self.b - self.a;
        let a = bias - local_axis.dot(&self.a.coords);
        let b = local_axis.dot(&dir);
        let bcoord = a / b;
        let dir_norm = dir.norm();

        if relative_eq!(b, 0.0)
            || bcoord * dir_norm <= epsilon
            || bcoord * dir_norm >= 1.0 - epsilon
        {
            // the plane is parallel to `self`.
            if a >= 0.0 {
                SplitResult::Negative
            } else {
                SplitResult::Positive
            }
        } else {
            let intersection = self.a + dir * bcoord;
            let s1 = Segment::new(self.a, intersection);
            let s2 = Segment::new(intersection, self.b);
            if a >= 0.0 {
                SplitResult::Pair(s1, s2)
            } else {
                SplitResult::Pair(s2, s1)
            }
        }
    }
}
