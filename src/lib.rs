/*!
parry
========

**parry** is a 2 and 3-dimensional geometric library written with
the rust programming language.

*/

#![deny(non_camel_case_types)]
#![deny(unused_parens)]
#![deny(non_upper_case_globals)]
#![deny(unused_qualifications)]
#![deny(unused_results)]
#![warn(missing_docs)]
#![warn(unused_imports)]
#![allow(missing_copy_implementations)]
#![doc(html_root_url = "http://docs.rs/parry/0.1.1")]

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
#[cfg(all(feature = "simd-is-enabled", feature = "f64"))]
std::compile_error!(
    "Explicit SIMD optimization are not yet supported when the f64 feature is enabled."
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

#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;
#[macro_use]
extern crate approx;
#[macro_use]
#[cfg(feature = "dim3")]
extern crate bitflags;
extern crate num_traits as num;

pub extern crate nalgebra as na;
pub extern crate simba;

pub mod bounding_volume;
pub mod mass_properties;
pub mod partitioning;
pub mod query;
pub mod shape;
pub mod transformation;
pub mod utils;

mod real {
    /// The scalar type used throughout this crate.
    #[cfg(feature = "f64")]
    pub type Real = f64;

    /// The scalar type used throughout this crate.
    #[cfg(feature = "f32")]
    pub type Real = f32;
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
        Vector3, U2,
    };

    /// The default tolerance used for geometric operations.
    pub const DEFAULT_EPSILON: Real = Real::EPSILON;

    /// The dimension of the space.
    pub const DIM: usize = 2;

    /// The dimension of the ambient space.
    pub type Dim = U2;

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
    #[allow(unused_imports)]
    #[cfg(feature = "simd-nightly")]
    use simba::simd::{f32x16, f32x4, f32x8, m32x16, m32x4, m32x8, u8x16, u8x4, u8x8};
    #[cfg(feature = "simd-stable")]
    use simba::simd::{WideBoolF32x4, WideF32x4};

    /// The number of lanes of a SIMD number.
    pub const SIMD_WIDTH: usize = 4;
    /// SIMD_WIDTH - 1
    pub const SIMD_LAST_INDEX: usize = 3;
    #[cfg(not(feature = "simd-nightly"))]
    /// A SIMD float with SIMD_WIDTH lanes.
    pub type SimdReal = WideF32x4;
    #[cfg(not(feature = "simd-nightly"))]
    /// A SIMD bool with SIMD_WIDTH lanes.
    pub type SimdBool = WideBoolF32x4;
    #[cfg(feature = "simd-nightly")]
    /// A SIMD float with SIMD_WIDTH lanes.
    pub type SimdReal = f32x4;
    #[cfg(feature = "simd-nightly")]
    /// A bool float with SIMD_WIDTH lanes.
    pub type SimdBool = m32x4;

    // pub const SIMD_WIDTH: usize = 8;
    // pub const SIMD_LAST_INDEX: usize = 7;
    // pub type SimdReal = f32x8;
    // pub type SimdBool = m32x8;

    // pub const SIMD_WIDTH: usize = 16;
    // pub const SIMD_LAST_INDEX: usize = 15;
    // pub type SimdReal = f32x16;
    // pub type SimdBool = m32x16;
}
