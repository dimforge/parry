// https://github.com/dimforge/parry/issues/242

use na::{Isometry3, Point3, Translation3, UnitQuaternion, Vector3};
use parry3d::query::Ray;
use parry3d::shape::{Ball, Cuboid, Shape};

fn run_test<S>(name: &str, shape: S)
where
    S: Shape,
{
    let mut rng = oorandom::Rand32::new(42);

    for _ in 0..1000 {
        let ray_origin = Point3::from(Vector3::from_fn(|_, _| rng.rand_float()).normalize() * 5.0);
        let ray = Ray::new(ray_origin, Point3::origin() - ray_origin);

        let rotation = if rng.rand_float() < 0.01 {
            UnitQuaternion::identity()
        } else {
            na::Unit::try_new(
                na::Quaternion::new(
                    rng.rand_float(),
                    rng.rand_float(),
                    rng.rand_float(),
                    rng.rand_float(),
                ),
                1.0e-5,
            )
            .unwrap_or(UnitQuaternion::identity())
        };
        let position = Isometry3::from_parts(Translation3::identity(), rotation);

        let intersection = shape
            .cast_ray_and_get_normal(&position, &ray, f32::MAX, true)
            .unwrap_or_else(|| panic!("Ray {ray:?} did not hit Shape {name} rotated with {rotation:?}"));

        let point = ray.origin + ray.dir * intersection.time_of_impact;
        let point_nudged_in = point + intersection.normal * -0.001;
        let point_nudged_out = point + intersection.normal * 0.001;

        assert!(
            shape.contains_point(&position, &point_nudged_in),
            "Shape {} rotated with {:#?} does not contain point nudged in {:#?}",
            name,
            rotation.axis(),
            point_nudged_in,
        );

        assert!(
            !shape.contains_point(&position, &point_nudged_out),
            "Shape {} rotated with {:#?} does contains point nudged out {:#?}",
            name,
            rotation.axis(),
            point_nudged_out,
        );

        let new_ray = Ray::new(point_nudged_out, ray_origin - point_nudged_out);

        assert!(
            shape
                .cast_ray_and_get_normal(&position, &new_ray, f32::MAX, true)
                .is_none(),
            "Ray {:#?} from outside Shape {} rotated with {:#?} did hit at t={}",
            ray,
            name,
            rotation,
            shape
                .cast_ray_and_get_normal(&position, &new_ray, f32::MAX, true)
                .expect("recurring ray cast produced a different answer")
                .time_of_impact
        );
    }
}

#[test]
fn shape_ray_cast_points_to_surface() {
    run_test("ball with radius 1", Ball::new(1.0));
    run_test(
        "cube with half-side 1",
        Cuboid::new(Vector3::new(1.0, 1.0, 1.0)),
    );
    run_test("tall rectangle", Cuboid::new(Vector3::new(1.0, 1.0, 0.5)));
    run_test(
        "tall and slim rectangle",
        Cuboid::new(Vector3::new(0.5, 1.0, 0.5)),
    );
}
