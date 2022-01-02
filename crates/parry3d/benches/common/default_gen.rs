use na::{
    self, Isometry2, Isometry3, Matrix2, Matrix3, Matrix4, Point2, Point3, Point4, RealField,
    Vector2, Vector3, Vector4,
};
use parry3d::bounding_volume::{BoundingSphere, AABB};
use parry3d::math::{Point, Real, Vector};
use parry3d::query::Ray;
use parry3d::shape::{Ball, Capsule, Cone, ConvexHull, Cuboid, Cylinder, Segment, Triangle};
use rand::distributions::{Distribution, Standard};
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
                rng.gen::<$t>()
            }
        }
    }
);

impl_rand_default_gen!(Vector2<f32>);
impl_rand_default_gen!(Vector3<f32>);
impl_rand_default_gen!(Vector4<f32>);
impl_rand_default_gen!(Point2<f32>);
impl_rand_default_gen!(Point3<f32>);
impl_rand_default_gen!(Point4<f32>);
impl_rand_default_gen!(Matrix2<f32>);
impl_rand_default_gen!(Matrix3<f32>);
impl_rand_default_gen!(Matrix4<f32>);
impl_rand_default_gen!(Isometry2<f32>);
impl_rand_default_gen!(Isometry3<f32>);
impl_rand_default_gen!(Vector2<f64>);
impl_rand_default_gen!(Vector3<f64>);
impl_rand_default_gen!(Vector4<f64>);
impl_rand_default_gen!(Point2<f64>);
impl_rand_default_gen!(Point3<f64>);
impl_rand_default_gen!(Point4<f64>);
impl_rand_default_gen!(Matrix2<f64>);
impl_rand_default_gen!(Matrix3<f64>);
impl_rand_default_gen!(Matrix4<f64>);
impl_rand_default_gen!(Isometry2<f64>);
impl_rand_default_gen!(Isometry3<f64>);
impl_rand_default_gen!(f32);
impl_rand_default_gen!(f64);
impl_rand_default_gen!(bool);

impl DefaultGen for Ball
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> Ball {
        Ball::new(rng.gen::<f32>().abs())
    }
}

impl DefaultGen for Cuboid
where
    Standard: Distribution<Vector<Real>>,
{
    fn generate<R: Rng>(rng: &mut R) -> Cuboid {
        Cuboid::new(rng.gen::<Vector<Real>>().abs())
    }
}

impl DefaultGen for Capsule
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> Capsule {
        Capsule::new(rng.gen::<Real>().abs(), rng.gen::<Real>().abs())
    }
}

impl DefaultGen for Cone
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> Cone {
        Cone::new(rng.gen::<Real>().abs(), rng.gen::<Real>().abs())
    }
}

impl DefaultGen for Cylinder
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> Cylinder {
        Cylinder::new(rng.gen::<Real>().abs(), rng.gen::<Real>().abs())
    }
}

impl DefaultGen for Segment
where
    Standard: Distribution<Point<Real>>,
{
    fn generate<R: Rng>(rng: &mut R) -> Segment {
        Segment::new(rng.gen(), rng.gen())
    }
}

impl DefaultGen for Triangle
where
    Standard: Distribution<Point<Real>>,
{
    fn generate<R: Rng>(rng: &mut R) -> Triangle {
        Triangle::new(rng.gen(), rng.gen(), rng.gen())
    }
}

impl DefaultGen for ConvexHull
where
    Standard: Distribution<Point<Real>>,
{
    fn generate<R: Rng>(rng: &mut R) -> ConvexHull {
        // It is recommended to have at most 100 points.
        // Otherwise, a smarter structure like the DK hierarchy would be needed.
        // let pts: Vec<_> = (0..100).map(|_| rng.gen()).collect();
        // ConvexHull::try_from_points(&pts).unwrap()
        unimplemented!()
    }
}

impl DefaultGen for Ray
where
    Standard: Distribution<Vector<Real>>,
{
    fn generate<R: Rng>(rng: &mut R) -> Ray {
        // The generate ray will always point to the origin.
        let shift = rng.gen::<Vector<Real>>() * na::convert::<_, Real>(10.0f64);
        Ray::new(Point::origin() + shift, -shift)
    }
}

impl DefaultGen for AABB
where
    Standard: Distribution<Vector<Real>>,
{
    fn generate<R: Rng>(rng: &mut R) -> AABB {
        // an AABB centered at the origin.
        let half_extents = rng.gen::<Vector<Real>>().abs();
        AABB::new(
            Point::origin() + (-half_extents),
            Point::origin() + half_extents,
        )
    }
}

impl DefaultGen for BoundingSphere
where
    Standard: Distribution<Real>,
{
    fn generate<R: Rng>(rng: &mut R) -> BoundingSphere {
        // a bounding sphere centered at the origin.
        BoundingSphere::new(Point::origin(), rng.gen::<Real>().abs())
    }
}
