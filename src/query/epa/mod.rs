//! The EPA algorithm for penetration depth computation.
//!
#[cfg(feature = "dim2")]
pub use self::epa2::EPA;
#[cfg(feature = "dim3")]
pub use self::epa3::EPA;

#[cfg(feature = "dim2")]
pub mod epa2;
#[cfg(feature = "dim3")]
pub mod epa3;
