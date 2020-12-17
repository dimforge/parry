//! A reference packed with a cost value.

use std::cmp::Ordering;

/// A reference packed with a cost value.
pub struct RefWithCost<'a, Real, T: 'a> {
    /// The reference to an object.
    pub object: &'a T,
    /// The cost of the object.
    pub cost: Real,
}

impl<'a, Real: PartialEq, T> PartialEq for RefWithCost<'a, Real, T> {
    #[inline]
    fn eq(&self, other: &RefWithCost<'a, Real, T>) -> bool {
        self.cost.eq(&other.cost)
    }
}

impl<'a, Real: PartialEq, T> Eq for RefWithCost<'a, Real, T> {}

impl<'a, Real: PartialOrd, T> PartialOrd for RefWithCost<'a, Real, T> {
    #[inline]
    fn partial_cmp(&self, other: &RefWithCost<'a, Real, T>) -> Option<Ordering> {
        self.cost.partial_cmp(&other.cost)
    }
}

impl<'a, Real: PartialOrd, T> Ord for RefWithCost<'a, Real, T> {
    #[inline]
    fn cmp(&self, other: &RefWithCost<'a, Real, T>) -> Ordering {
        if self.cost < other.cost {
            Ordering::Less
        } else if self.cost > other.cost {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }
}
