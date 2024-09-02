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
#![warn(missing_docs)] // TODO: deny this
#![warn(unused_imports)]
#![allow(missing_copy_implementations)]
#![allow(clippy::too_many_arguments)] // Maybe revisit this one later.
#![allow(clippy::module_inception)]
#![allow(clippy::manual_range_contains)] // This usually makes it way more verbose that it could be.
#![allow(clippy::type_complexity)] // Complains about closures that are fairly simple.
#![doc(html_root_url = "http://docs.rs/parry/0.1.1")]
#![cfg_attr(not(feature = "std"), no_std)]
#![cfg_attr(not(feature = "rkyv"), deny(unused_qualifications))] // TODO: deny that everytime

#[cfg(all(
    feature = "simd-is-enabled",
    not(feature = "simd-stable"),
    not(feature = "simd-nightly")
))]
std::compile_error!("The `simd-is-enabled` feature should not be enabled explicitly. Please enable the `simd-stable` or the `simd-nightly` feature instead.");
#[cfg(all(feature = "simd-is-enabled", feature = "enhanced-determinism"))]
std::compile_error!(
    "SIMD cannot be enabled when the `enhanced-determinism` feature is also enabled."
);

macro_rules! array(
    ($callback: expr; SIMD_WIDTH) => {
        {
            #[inline(always)]
            #[allow(dead_code)]
            fn create_arr<T>(mut callback: impl FnMut(usize) -> T) -> [T; SIMD_WIDTH] {
                [callback(0usize), callback(1usize), callback(2usize), callback(3usize)]
            }

            create_arr($callback)
        }
    }
);

#[cfg(all(feature = "alloc", not(feature = "std")))]
#[cfg_attr(test, macro_use)]
extern crate alloc;

#[cfg(not(feature = "std"))]
extern crate core as std;

#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;
#[macro_use]
extern crate approx;
extern crate num_traits as num;

pub extern crate either;
pub extern crate nalgebra as na;
pub extern crate simba;

pub mod bounding_volume;
pub mod mass_properties;
pub mod partitioning;
pub mod query;
pub mod shape;
#[cfg(feature = "std")]
pub mod transformation;
pub mod utils;

mod real {
    /// The scalar type used throughout this crate.
    #[cfg(feature = "f64")]
    pub use f64 as Real;

    /// The scalar type used throughout this crate.
    #[cfg(feature = "f32")]
    pub use f32 as Real;
}

/// Compilation flags dependent aliases for mathematical types.
#[cfg(feature = "dim3")]
pub mod math {
    pub use super::real::*;
    pub use super::simd::*;
    use na::{
        Isometry3, Matrix3, Point3, Translation3, UnitQuaternion, UnitVector3, Vector3, Vector6,
        U3, U6,
    };

    /// The default tolerance used for geometric operations.
    pub const DEFAULT_EPSILON: Real = Real::EPSILON;

    /// The dimension of the space.
    pub const DIM: usize = 3;

    /// The dimension of the space multiplied by two.
    pub const TWO_DIM: usize = DIM * 2;

    /// The dimension of the ambient space.
    pub type Dim = U3;

    /// The dimension of a spatial vector.
    pub type SpatialDim = U6;

    /// The dimension of the rotations.
    pub type AngDim = U3;

    /// The point type.
    pub type Point<N> = Point3<N>;

    /// The angular vector type.
    pub type AngVector<N> = Vector3<N>;

    /// The vector type.
    pub type Vector<N> = Vector3<N>;

    /// The unit vector type.
    pub type UnitVector<N> = UnitVector3<N>;

    /// The matrix type.
    pub type Matrix<N> = Matrix3<N>;

    /// The vector type with dimension `SpatialDim × 1`.
    pub type SpatialVector<N> = Vector6<N>;

    /// The orientation type.
    pub type Orientation<N> = Vector3<N>;

    /// The transformation matrix type.
    pub type Isometry<N> = Isometry3<N>;

    /// The rotation matrix type.
    pub type Rotation<N> = UnitQuaternion<N>;

    /// The translation type.
    pub type Translation<N> = Translation3<N>;

    /// The angular inertia of a rigid body.
    pub type AngularInertia<N> = crate::utils::SdpMatrix3<N>;

    /// The principal angular inertia of a rigid body.
    pub type PrincipalAngularInertia<N> = Vector3<N>;

    /// A matrix that represent the cross product with a given vector.
    pub type CrossMatrix<N> = Matrix3<N>;

    /// A vector with a dimension equal to the maximum number of degrees of freedom of a rigid body.
    pub type SpacialVector<N> = Vector6<N>;

    /// A 3D symmetric-definite-positive matrix.
    pub type SdpMatrix<N> = crate::utils::SdpMatrix3<N>;
}

/// Compilation flags dependent aliases for mathematical types.
#[cfg(feature = "dim2")]
pub mod math {
    pub use super::real::*;
    pub use super::simd::*;
    use na::{
        Isometry2, Matrix2, Point2, Translation2, UnitComplex, UnitVector2, Vector1, Vector2,
        Vector3, U1, U2,
    };

    /// The default tolerance used for geometric operations.
    pub const DEFAULT_EPSILON: Real = Real::EPSILON;

    /// The dimension of the space.
    pub const DIM: usize = 2;

    /// The dimension of the space multiplied by two.
    pub const TWO_DIM: usize = DIM * 2;

    /// The dimension of the ambient space.
    pub type Dim = U2;

    /// The dimension of the rotations.
    pub type AngDim = U1;

    /// The point type.
    pub type Point<N> = Point2<N>;

    /// The angular vector type.
    pub type AngVector<N> = N;

    /// The vector type.
    pub type Vector<N> = Vector2<N>;

    /// The unit vector type.
    pub type UnitVector<N> = UnitVector2<N>;

    /// The matrix type.
    pub type Matrix<N> = Matrix2<N>;

    /// The orientation type.
    pub type Orientation<N> = Vector1<N>;

    /// The transformation matrix type.
    pub type Isometry<N> = Isometry2<N>;

    /// The rotation matrix type.
    pub type Rotation<N> = UnitComplex<N>;

    /// The translation type.
    pub type Translation<N> = Translation2<N>;

    /// The angular inertia of a rigid body.
    pub type AngularInertia<N> = N;

    /// The principal angular inertia of a rigid body.
    pub type PrincipalAngularInertia<N> = N;

    /// A matrix that represent the cross product with a given vector.
    pub type CrossMatrix<N> = Vector2<N>;

    /// A vector with a dimension equal to the maximum number of degrees of freedom of a rigid body.
    pub type SpacialVector<N> = Vector3<N>;

    /// A 2D symmetric-definite-positive matrix.
    pub type SdpMatrix<N> = crate::utils::SdpMatrix2<N>;
}

#[cfg(not(feature = "simd-is-enabled"))]
mod simd {
    use simba::simd::AutoBoolx4;
    /// The number of lanes of a SIMD number.
    pub const SIMD_WIDTH: usize = 4;
    /// SIMD_WIDTH - 1
    pub const SIMD_LAST_INDEX: usize = 3;

    /// A SIMD float with SIMD_WIDTH lanes.
    #[cfg(feature = "f32")]
    pub type SimdReal = simba::simd::AutoF32x4;

    /// A SIMD float with SIMD_WIDTH lanes.
    #[cfg(feature = "f64")]
    pub type SimdReal = simba::simd::AutoF64x4;

    /// A SIMD bool with SIMD_WIDTH lanes.
    pub type SimdBool = AutoBoolx4;
}

#[cfg(feature = "simd-is-enabled")]
mod simd {
    #[cfg(all(feature = "simd-nightly", feature = "f32"))]
    pub use simba::simd::{f32x4 as SimdReal, mask32x4 as SimdBool};
    #[cfg(all(feature = "simd-stable", feature = "f32"))]
    pub use simba::simd::{WideBoolF32x4 as SimdBool, WideF32x4 as SimdReal};

    #[cfg(all(feature = "simd-nightly", feature = "f64"))]
    pub use simba::simd::{f64x4 as SimdReal, mask64x4 as SimdBool};
    #[cfg(all(feature = "simd-stable", feature = "f64"))]
    pub use simba::simd::{WideBoolF64x4 as SimdBool, WideF64x4 as SimdReal};

    /// The number of lanes of a SIMD number.
    pub const SIMD_WIDTH: usize = 4;
    /// SIMD_WIDTH - 1
    pub const SIMD_LAST_INDEX: usize = 3;
}
