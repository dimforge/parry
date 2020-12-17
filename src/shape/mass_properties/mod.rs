pub use self::mass_properties::MassProperties;

mod mass_properties;
mod mass_properties_ball;
mod mass_properties_capsule;
#[cfg(feature = "dim3")]
mod mass_properties_cone;
mod mass_properties_cuboid;
mod mass_properties_cylinder;
#[cfg(feature = "dim2")]
mod mass_properties_polygon;
