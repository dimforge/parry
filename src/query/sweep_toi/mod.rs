//! Time-of-impact computation on endpoint-interpolated sweeps.
//!
//! Unlike [`cast_shapes_nonlinear`](crate::query::cast_shapes_nonlinear), which models motion
//! as constant velocities, this module models a timestep as a [`Sweep`] between two endpoint
//! poses (linear center-of-mass interpolation + rotation nlerp) and computes the earliest
//! time at which two swept shapes reach a slop-based target separation, using conservative
//! advancement with separation functions.
//!
//! NOTE: this is mostly ported from Box2D which had much better CCD quality than Rapier.
//! TODO: see how we can combine that with `nonlinear_shape_cast` since they serve similar goals.

pub use self::proxy_distance::{proxy_distance, ProxyDistanceOutput, SimplexCache};
pub use self::sweep::Sweep;
pub use self::sweep_toi::{sweep_time_of_impact, SweepToiOutput, SweepToiStatus};
pub use self::toi_proxy::{ToiProxy, TOI_PROXY_INLINE_POINTS};

#[cfg(feature = "alloc")]
pub use self::composite::{sweep_time_of_impact_composite, SweepCompositeFastShape, CORE_FRACTION};

#[cfg(feature = "alloc")]
mod composite;
mod proxy_distance;
mod separation;
mod sweep;
mod sweep_toi;
mod toi_proxy;
