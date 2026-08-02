/*!
parry
========

**parry** is a 2 and 3-dimensional geometric library written with
the rust programming language.

*/

#![deny(non_camel_case_types)]
#![deny(unused_parens)]
#![deny(non_upper_case_globals)]
#![deny(unused_results)]
#![deny(unused_qualifications)]
#![warn(missing_docs)]
#![warn(unused_imports)]
#![allow(missing_copy_implementations)]
#![allow(clippy::too_many_arguments)] // Maybe revisit this one later.
#![allow(clippy::module_inception)]
#![allow(clippy::manual_range_contains)] // This usually makes it way more verbose that it could be.
#![allow(clippy::type_complexity)] // Complains about closures that are fairly simple.
#![cfg_attr(feature = "dim2", doc(html_root_url = "https://docs.rs/parry2d"))]
#![cfg_attr(feature = "dim3", doc(html_root_url = "https://docs.rs/parry3d"))]
#![no_std]

#[cfg(all(
    feature = "simd-is-enabled",
    not(feature = "simd-stable"),
    not(feature = "simd-nightly")
))]
std::compile_error!("The `simd-is-enabled` feature should not be enabled explicitly. Please enable the `simd-stable` or the `simd-nightly` feature instead.");
#[cfg(all(feature = "simd8", feature = "enhanced-determinism"))]
core::compile_error!(
    "8-lanes SIMD cannot be enabled when the `enhanced-determinism` feature is also enabled because it breaks cross-platform determinism."
);

#[allow(unused_macros)]
macro_rules! array(
    ($callback: expr; SIMD_WIDTH) => {
        {
            #[inline(always)]
            #[allow(dead_code)]
            fn create_arr<T>(callback: impl FnMut(usize) -> T) -> [T; SIMD_WIDTH] {
                // Width-agnostic: `N` is inferred from the `[T; SIMD_WIDTH]` return type,
                // so this covers the 1-, 4-, and 8-lane builds alike.
                core::array::from_fn(callback)
            }

            create_arr($callback)
        }
    }
);

#[cfg(feature = "std")]
extern crate std;

#[cfg(feature = "alloc")]
#[cfg_attr(test, macro_use)]
extern crate alloc;

#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;
#[macro_use]
extern crate approx;
extern crate num_traits as num;

pub extern crate either;
pub extern crate glamx;
pub extern crate simba;

#[macro_use]
mod macros;

pub mod bounding_volume;
pub mod mass_properties;
pub mod math;
pub mod partitioning;
pub mod query;
pub mod shape;
#[cfg(feature = "alloc")]
pub mod transformation;
pub mod utils;

mod simd {
    // 8-lane SIMD (f32 only; simba has no `WideF64x8`). Opt-in via `simd8`
    // on top of `simd-stable`/`simd-nightly`. Requires an AVX-enabled target
    // (`RUSTFLAGS="-C target-feature=+avx2,+fma"` or `-C target-cpu=native`) for
    // the compiler to actually emit 256-bit instructions; otherwise it runs
    // (correctly) as two 128-bit halves.
    #[cfg(all(feature = "simd8", feature = "simd-nightly", feature = "f32"))]
    pub use simba::simd::{f32x8 as SimdReal, mask32x8 as SimdBool};
    #[cfg(all(feature = "simd8", feature = "simd-stable", feature = "f32"))]
    pub use simba::simd::{WideBoolF32x8 as SimdBool, WideF32x8 as SimdReal};

    // 4-lane SIMD (default width).
    #[cfg(all(not(feature = "simd8"), feature = "simd-nightly", feature = "f32"))]
    pub use simba::simd::{f32x4 as SimdReal, mask32x4 as SimdBool};
    #[cfg(all(not(feature = "simd-is-enabled"), feature = "f32"))]
    pub use simba::simd::{AutoBoolx4 as SimdBool, AutoF32x4 as SimdReal};
    #[cfg(all(not(feature = "simd8"), feature = "simd-stable", feature = "f32"))]
    pub use simba::simd::{WideBoolF32x4 as SimdBool, WideF32x4 as SimdReal};

    // f64 stays 4-lane regardless of `simd8` (no 8-lane f64 type in simba).
    #[cfg(all(feature = "simd-nightly", feature = "f64"))]
    pub use simba::simd::{f64x4 as SimdReal, mask64x4 as SimdBool};
    #[cfg(all(not(feature = "simd-is-enabled"), feature = "f64"))]
    pub use simba::simd::{AutoBoolx4 as SimdBool, AutoF64x4 as SimdReal};
    #[cfg(all(feature = "simd-stable", feature = "f64"))]
    pub use simba::simd::{WideBoolF64x4 as SimdBool, WideF64x4 as SimdReal};

    /// The number of lanes of a SIMD number.
    #[cfg(all(feature = "simd8", feature = "f32"))]
    pub const SIMD_WIDTH: usize = 8;
    /// SIMD_WIDTH - 1
    #[cfg(all(feature = "simd8", feature = "f32"))]
    pub const SIMD_LAST_INDEX: usize = 7;

    /// The number of lanes of a SIMD number.
    #[cfg(not(all(feature = "simd8", feature = "f32")))]
    pub const SIMD_WIDTH: usize = 4;
    /// SIMD_WIDTH - 1
    #[cfg(not(all(feature = "simd8", feature = "f32")))]
    pub const SIMD_LAST_INDEX: usize = 3;
}
