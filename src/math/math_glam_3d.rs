pub use super::{Real, SimdIso3, SimdMat3, SimdQuat, SimdVec3};
use crate::math::{Iso3, SimdVec2};
use crate::utils::SdpMatrix3;
use glam::{Mat2, Mat3, Quat, UVec3, Vec2, Vec3};
use na::{Vector6, U3};

/// The default tolerance used for geometric operations.
pub const DEFAULT_EPSILON: Real = Real::EPSILON;

/// The dimension of the space.
pub const DIM: usize = 3;

/// The dimension of the space multiplied by two.
pub const TWO_DIM: usize = DIM * 2;

/// The dimension of the ambient space.
pub type Dim = U3;

/// The dimension of the rotations.
pub type AngDim = U3;

/// The point type.
pub type Point = Vec3;

/// The 2D point type.
pub type Point2 = Vec2;

/// An AoSoA SIMD point type.
pub type SimdPoint = SimdVec3;

/// The angular vector type.
pub type AngVector = Vec3;

/// The SIMD angular vector type.
pub type SimdAngVector = SimdVec3;

/// The vector type.
pub type Vector = Vec3;

/// The 2D vector type.
pub type Vector2 = Vec2;

/// The SIMD 2D vector type.
pub type SimdVector2 = SimdVec2;

/// A 6D vector, generally used to combine tanslation (3D) + rotation (3D) degrees of freedoms.
pub type SpatialVector<N> = Vector6<N>; // TODO: replace by a glam type

/// The u32 vector type.
pub type VectorU32 = UVec3;

/// An AoSoA SIMD vector.
pub type SimdVector = SimdVec3;

/// The unit vector type.
pub type UnitVector = Vec3;

/// The matrix type.
pub type Matrix = Mat3;

/// The SIMD matrix type.
pub type SimdMatrix = SimdMat3;

/// The 2D matrix type.
pub type Matrix2 = Mat2;

/// The orientation type.
pub type Orientation = Vec3;

/// The transformation matrix type.
pub type Isometry = Iso3;

/// The AoSoA transformation matrix type.
pub type SimdIsometry = SimdIso3;

/// The rotation type.
pub type Rotation = Quat;

/// The SIMD rotation type.
pub type SimdRotation = SimdQuat;

/// The rotation matrix type.
pub type RotationMatrix = Mat3;

/// The translation type.
pub type Translation = Iso3;

/// The angular inertia of a rigid body.
pub type AngularInertia<N> = SdpMatrix3<N>;

/// The principal angular inertia of a rigid body.
pub type PrincipalAngularInertia = Vec3;

/// A matrix that represent the cross product with a given vector.
pub type CrossMatrix = Mat3;

/// A 3D symmetric-definite-positive matrix.
pub type SdpMatrix<N> = SdpMatrix3<N>;
