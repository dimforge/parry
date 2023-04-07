use std::cmp::PartialOrd;
use std::mem;
use std::ops::Deref;

/// A pair of elements sorted in increasing order.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
pub struct SortedPair<T: PartialOrd>([T; 2]);

impl<T: PartialOrd> SortedPair<T> {
    /// Sorts two elements in increasing order into a new pair.
    pub fn new(element1: T, element2: T) -> Self {
        if element1 > element2 {
            SortedPair([element2, element1])
        } else {
            SortedPair([element1, element2])
        }
    }
}

impl<T: PartialOrd> Deref for SortedPair<T> {
    type Target = (T, T);

    fn deref(&self) -> &(T, T) {
        unsafe { mem::transmute(self) }
    }
}

// TODO: can we avoid these manual impls of Hash/PartialEq/Eq for the archived types?
#[cfg(feature = "rkyv")]
impl<T: rkyv::Archive + PartialOrd> std::hash::Hash for ArchivedSortedPair<T>
where
    [<T as rkyv::Archive>::Archived; 2]: std::hash::Hash,
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.0.hash(state)
    }
}

#[cfg(feature = "rkyv")]
impl<T: rkyv::Archive + PartialOrd> PartialEq for ArchivedSortedPair<T>
where
    [<T as rkyv::Archive>::Archived; 2]: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

#[cfg(feature = "rkyv")]
impl<T: rkyv::Archive + PartialOrd> Eq for ArchivedSortedPair<T> where
    [<T as rkyv::Archive>::Archived; 2]: Eq
{
}
