pub use self::default_gen::generate;
#[cfg(feature = "alloc")]
pub use self::generators::generate_trimesh_around_origin;
pub use self::unref::unref;

mod default_gen;
#[cfg(feature = "alloc")]
mod generators;
mod unref;
