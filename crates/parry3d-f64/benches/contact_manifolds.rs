//! f64 variant of the parry3d contact-manifold benchmarks: the bench body is
//! shared with parry3d (it is written against `parry3d::math::Real`) and the
//! extern-crate alias below rebinds it to the f64 build.
extern crate parry3d_f64 as parry3d;

include!("../../parry3d/benches/contact_manifolds.rs");
