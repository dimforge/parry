pub use super::{Iso2, Real, SimdIso2, SimdMat2, SimdReal, SimdVec2};
use crate::utils::SdpMatrix2;
use glam::{Mat2, UVec2, Vec2};
use na::{Vector3, U1, U2};

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
pub type Point = Vec2;

/// The 2D point type.
pub type Point2 = Vec2;

/// An AoSoA SIMD point type.
pub type SimdPoint = SimdVec2;

/// The angular vector type.
pub type AngVector = Real;

/// The SIMD angular vector type.
pub type SimdAngVector = SimdReal;

/// The vector type.
pub type Vector = Vec2;

/// The 2D vector type.
pub type Vector2 = Vec2;

/// The SIMD 2D vector type.
pub type SimdVector2 = SimdVec2;

/// A 3D vector, generally used to combine tanslation (2D) + rotation (1D) degrees of freedoms.
pub type SpatialVector<N> = Vector3<N>; // TODO: replace by a glam type

/// The u22 vector type.
pub type VectorU32 = UVec2;

/// An AoSoA SIMD vector.
pub type SimdVector = SimdVec2;

/// The unit vector type.
pub type UnitVector = Vec2;

/// The matrix type.
pub type Matrix = Mat2;

/// The SIMD matrix type.
pub type SimdMatrix = SimdMat2;

/// The 2D matrix type.
pub type Matrix2 = Mat2;

/// The orientation type.
pub type Orientation = Real;

/// The transformation matrix type.
pub type Isometry = Iso2;

/// The AoSoA transformation matrix type.
pub type SimdIsometry = SimdIso2;

/// The rotation type.
pub type Rotation = Mat2;

/// The SIMD rotation type.
pub type SimdRotation = SimdMat2;

/// The rotation matrix type.
pub type RotationMatrix = Mat2;

/// The translation type.
pub type Translation = Iso2;

/// The angular inertia of a rigid body.
pub type AngularInertia<N> = N;

/// The principal angular inertia of a rigid body.
pub type PrincipalAngularInertia = Real;

/// A matrix that represent the cross product with a given vector.
pub type CrossMatrix = Mat2;

/// A 2D symmetric-definite-positive matrix.
pub type SdpMatrix<N> = SdpMatrix2<N>;
