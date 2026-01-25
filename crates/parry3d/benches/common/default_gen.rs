use parry3d::bounding_volume::{Aabb, BoundingSphere};
use parry3d::math::{Pose3, Real, Rot3, Vec3, Vector};
use parry3d::query::Ray;
use parry3d::shape::ConvexPolyhedron;
use parry3d::shape::{Ball, Capsule, Cone, Cuboid, Cylinder, Segment, Triangle};
use rand::Rng;

pub trait DefaultGen {
    fn generate<R: Rng>(rng: &mut R) -> Self;
}

pub fn generate<T: DefaultGen, R: Rng>(rng: &mut R) -> T {
    DefaultGen::generate(rng)
}

macro_rules! impl_rand_default_gen (
    ($t: ty) => {
        impl DefaultGen for $t {
            fn generate<R: Rng>(rng: &mut R) -> $t {
                rng.random::<$t>()
            }
        }
    }
);

impl_rand_default_gen!(f32);
impl_rand_default_gen!(f64);
impl_rand_default_gen!(bool);

// Implement DefaultGen for glam types manually since they don't implement Distribution
impl DefaultGen for Vec3 {
    fn generate<R: Rng>(rng: &mut R) -> Vec3 {
        Vec3::new(rng.random(), rng.random(), rng.random())
    }
}

impl DefaultGen for Pose3 {
    fn generate<R: Rng>(rng: &mut R) -> Pose3 {
        let translation = Vec3::new(rng.random(), rng.random(), rng.random());
        let axis = Vec3::new(
            rng.random::<f32>() - 0.5,
            rng.random::<f32>() - 0.5,
            rng.random::<f32>() - 0.5,
        )
        .normalize_or_zero();
        let angle: f32 = rng.random::<f32>() * core::f32::consts::TAU;
        let rotation = Rot3::from_axis_angle(axis, angle);
        Pose3::from_parts(translation, rotation)
    }
}

impl DefaultGen for Ball {
    fn generate<R: Rng>(rng: &mut R) -> Ball {
        Ball::new(rng.random::<f32>().abs())
    }
}

impl DefaultGen for Cuboid {
    fn generate<R: Rng>(rng: &mut R) -> Cuboid {
        Cuboid::new(Vec3::new(
            rng.random::<f32>().abs(),
            rng.random::<f32>().abs(),
            rng.random::<f32>().abs(),
        ))
    }
}

impl DefaultGen for Capsule {
    fn generate<R: Rng>(rng: &mut R) -> Capsule {
        Capsule::new(
            Vec3::new(rng.random(), rng.random(), rng.random()),
            Vec3::new(rng.random(), rng.random(), rng.random()),
            rng.random::<Real>().abs(),
        )
    }
}

impl DefaultGen for Cone {
    fn generate<R: Rng>(rng: &mut R) -> Cone {
        Cone::new(rng.random::<Real>().abs(), rng.random::<Real>().abs())
    }
}

impl DefaultGen for Cylinder {
    fn generate<R: Rng>(rng: &mut R) -> Cylinder {
        Cylinder::new(rng.random::<Real>().abs(), rng.random::<Real>().abs())
    }
}

impl DefaultGen for Segment {
    fn generate<R: Rng>(rng: &mut R) -> Segment {
        Segment::new(
            Vec3::new(rng.random(), rng.random(), rng.random()),
            Vec3::new(rng.random(), rng.random(), rng.random()),
        )
    }
}

impl DefaultGen for Triangle {
    fn generate<R: Rng>(rng: &mut R) -> Triangle {
        Triangle::new(
            Vec3::new(rng.random(), rng.random(), rng.random()),
            Vec3::new(rng.random(), rng.random(), rng.random()),
            Vec3::new(rng.random(), rng.random(), rng.random()),
        )
    }
}

impl DefaultGen for ConvexPolyhedron {
    fn generate<R: Rng>(rng: &mut R) -> ConvexPolyhedron {
        // It is recommended to have at most 100 points.
        // Otherwise, a smarter structure like the DK hierarchy would be needed.
        let pts: Vec<_> = (0..100)
            .map(|_| Vec3::new(rng.random(), rng.random(), rng.random()))
            .collect();
        ConvexPolyhedron::from_convex_hull(&pts).unwrap()
    }
}

impl DefaultGen for Ray {
    fn generate<R: Rng>(rng: &mut R) -> Ray {
        // The generated ray will always point to the origin.
        let shift = Vec3::new(
            rng.random::<f32>() * 10.0,
            rng.random::<f32>() * 10.0,
            rng.random::<f32>() * 10.0,
        );
        Ray::new(Vector::ZERO + shift, -shift)
    }
}

impl DefaultGen for Aabb {
    fn generate<R: Rng>(rng: &mut R) -> Aabb {
        // an Aabb centered at the origin.
        let half_extents = Vec3::new(
            rng.random::<f32>().abs(),
            rng.random::<f32>().abs(),
            rng.random::<f32>().abs(),
        );
        Aabb::new(Vector::ZERO + (-half_extents), Vector::ZERO + half_extents)
    }
}

impl DefaultGen for BoundingSphere {
    fn generate<R: Rng>(rng: &mut R) -> BoundingSphere {
        // a bounding sphere centered at the origin.
        BoundingSphere::new(Vector::ZERO, rng.random::<Real>().abs())
    }
}
