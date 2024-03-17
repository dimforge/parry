pub use super::{Real, SimdReal};
use na::{
    Isometry3, Matrix3, Point3, Rotation3, Translation3, UnitQuaternion, UnitVector3, Vector3,
    Vector6, U3,
};
use simba::simd::SimdRealField;

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
pub type Point = Point3<Real>;

/// The 2D point type.
pub type Point2 = na::Point2<Real>;

/// An AoSoA SIMD point type.
pub type SimdPoint = Point3<SimdReal>;

/// The angular vector type.
pub type AngVector = Vector3<Real>;

/// The SIMD angular vector type.
pub type SimdAngVector = Vector3<SimdReal>;

/// The vector type.
pub type Vector = Vector3<Real>;

/// The 2D vector type.
pub type Vector2 = na::Vector2<Real>;

/// The SIMD 2D vector type.
pub type SimdVector2 = na::Vector2<SimdReal>;

/// The u32 vector type.
pub type VectorU32 = Vector3<u32>;

/// An AoSoA SIMD vector.
pub type SimdVector = Vector3<SimdReal>;

/// The unit vector type.
pub type UnitVector = UnitVector3<Real>;

/// The matrix type.
pub type Matrix = Matrix3<Real>;

/// The SIMD matrix type.
pub type SimdMatrix = Matrix3<SimdReal>;

/// The 2D matrix type.
pub type Matrix2 = na::Matrix2<Real>;

/// The orientation type.
pub type Orientation = Vector3<Real>;

/// The transformation matrix type.
pub type Isometry = Isometry3<Real>;

/// The AoSoA transformation matrix type.
pub type SimdIsometry = Isometry3<SimdReal>;

/// The rotation type.
pub type Rotation = UnitQuaternion<Real>;

/// The SIMD rotation type.
pub type SimdRotation = UnitQuaternion<SimdReal>;

/// The rotation matrix type.
pub type RotationMatrix = Rotation3<Real>;

/// The translation type.
pub type Translation = Translation3<Real>;

/// The angular inertia of a rigid body.
pub type AngularInertia<N> = crate::utils::SdpMatrix3<N>;

/// The principal angular inertia of a rigid body.
pub type PrincipalAngularInertia = Vector3<Real>;

/// A matrix that represent the cross product with a given vector.
pub type CrossMatrix = Matrix3<Real>;

/// A 3D symmetric-definite-positive matrix.
pub type SdpMatrix<N> = crate::utils::SdpMatrix3<N>;

/// A 6D vector, generally used to combine tanslation (3D) + rotation (3D) degrees of freedoms.
pub type SpatialVector<N> = Vector6<N>;
