use std::mem;
use std::slice;

use na::{Point2, Point3, Vector2, Vector3};
use simba::scalar::RealField;

/// Trait that transforms thing to a slice of u8.
pub trait AsBytes {
    /// Converts `self` to a slice of bytes.
    fn as_bytes(&self) -> &[u8];
}

macro_rules! generic_as_bytes_impl(
    ($T: ident, $dimension: expr) => (
        impl<N: RealField> AsBytes for $T<N> {
            #[inline(always)]
            fn as_bytes(&self) -> &[u8] {
                unsafe {
                    slice::from_raw_parts(self as *const $T<N> as *const u8, mem::size_of::<N>() * $dimension)
                }
            }
        }
    )
);

generic_as_bytes_impl!(Vector2, 2);
generic_as_bytes_impl!(Point2, 2);
generic_as_bytes_impl!(Vector3, 3);
generic_as_bytes_impl!(Point3, 3);

// FIXME: implement for all `T: Copy` instead?
