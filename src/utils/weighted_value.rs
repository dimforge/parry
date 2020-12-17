use crate::math::Real;
use std::cmp::{Eq, Ord, Ordering, PartialEq, PartialOrd};

pub struct WeightedValue<T> {
    pub value: T,
    pub cost: Real,
}

impl<T> WeightedValue<T> {
    /// Creates a new reference packed with a cost value.
    #[inline]
    pub fn new(value: T, cost: Real) -> WeightedValue<T> {
        WeightedValue { value, cost }
    }
}

impl<T> PartialEq for WeightedValue<T> {
    #[inline]
    fn eq(&self, other: &WeightedValue<T>) -> bool {
        self.cost.eq(&other.cost)
    }
}

impl<T> Eq for WeightedValue<T> {}

impl<T> PartialOrd for WeightedValue<T> {
    #[inline]
    fn partial_cmp(&self, other: &WeightedValue<T>) -> Option<Ordering> {
        self.cost.partial_cmp(&other.cost)
    }
}

impl<T> Ord for WeightedValue<T> {
    #[inline]
    fn cmp(&self, other: &WeightedValue<T>) -> Ordering {
        if self.cost < other.cost {
            Ordering::Less
        } else if self.cost > other.cost {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }
}
