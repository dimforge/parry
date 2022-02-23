use crate::math::{Point, Real, UnitVector, Vector};
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
        self.local_split_and_get_intersection(local_axis, bias, epsilon)
            .0
    }
}

impl Segment {
    pub fn local_split_and_get_intersection(
        &self,
        local_axis: &UnitVector<Real>,
        bias: Real,
        epsilon: Real,
    ) -> (SplitResult<Self>, Option<(Point<Real>, Real)>) {
        let dir = self.b - self.a;
        let a = bias - local_axis.dot(&self.a.coords);
        let b = local_axis.dot(&dir);
        let bcoord = a / b;
        let dir_norm = dir.norm();

        if relative_eq!(b, 0.0)
            || bcoord * dir_norm <= epsilon
            || bcoord * dir_norm >= dir_norm - epsilon
        {
            if a >= 0.0 {
                (SplitResult::Negative, None)
            } else {
                (SplitResult::Positive, None)
            }
        } else {
            let intersection = self.a + dir * bcoord;
            let s1 = Segment::new(self.a, intersection);
            let s2 = Segment::new(intersection, self.b);
            if a >= 0.0 {
                (SplitResult::Pair(s1, s2), Some((intersection, bcoord)))
            } else {
                (SplitResult::Pair(s2, s1), Some((intersection, bcoord)))
            }
        }
    }
}
