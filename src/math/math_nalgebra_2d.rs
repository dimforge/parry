pub use super::simd::*;
use super::Real;
use na::{Isometry2, Matrix2, Translation2, UnitComplex, UnitVector2, Vector1, Vector6, U1, U2};
use simba::simd::SimdRealField;

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
pub type Point = na::Point2<Real>;

/// The point type.
pub type Point2 = na::Point2<Real>;

/// An AoSoA SIMD point.
pub type SimdPoint = na::Point2<SimdReal>;

/// The angular vector type.
pub type AngVector = Real;

/// The SIMD angular vector type.
pub type SimdAngVector = SimdReal;

/// The vector type.
pub type Vector = na::Vector2<Real>;

/// The 2D vector type.
pub type Vector2 = na::Vector2<Real>;

/// The u32 vector type.
pub type VectorU32 = na::Vector2<u32>;

/// An AoSoA SIMD vector.
pub type SimdVector = na::Vector2<SimdReal>;

/// The unit vector type.
pub type UnitVector = UnitVector2<Real>;

/// The matrix type.
pub type Matrix = Matrix2<Real>;

/// The SIMD matrix type.
pub type SimdMatrix = Matrix2<SimdReal>;

/// The orientation type.
pub type Orientation = Vector1<Real>;

/// The transformation matrix type.
pub type Isometry = Isometry2<Real>;

/// The AoSoA transformation matrix type.
pub type SimdIsometry = Isometry2<SimdReal>;

/// The rotation matrix type.
pub type Rotation = UnitComplex<Real>;

/// The SIMD rotation type.
pub type SimdRotation = UnitComplex<SimdReal>;

/// The translation type.
pub type Translation = Translation2<Real>;

/// The angular inertia of a rigid body.
pub type AngularInertia<N> = N;

/// The principal angular inertia of a rigid body.
pub type PrincipalAngularInertia = Real;

/// A matrix that represent the cross product with a given vector.
pub type CrossMatrix = na::Vector2<Real>;

/// A 2D symmetric-definite-positive matrix.
pub type SdpMatrix = crate::utils::SdpMatrix2<Real>;

/// A 3D vector, generally used to combine tanslation (2D) + rotation (1D) degrees of freedoms.
pub type SpatialVector<N> = na::Vector3<N>;
