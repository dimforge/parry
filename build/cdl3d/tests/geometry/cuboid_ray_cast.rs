// https://github.com/rustsim/cdl/issues/242

use cdl3d::query::{PointQuery, Ray, RayCast};
use cdl3d::shape::{Ball, Cuboid, Shape};
use na::{Isometry3, Point3, Translation3, UnitQuaternion, Vector3};

fn run_test<S>(name: &str, shape: S)
where
    S: Shape,
{
    for _ in 0..1000 {
        let ray_origin = Point3::from(rand::random::<Vector3<f32>>().normalize() * 5.0);
        let ray = Ray::new(ray_origin, Point3::origin() - ray_origin);

        let rotation = if rand::random::<f32>() < 0.01 {
            UnitQuaternion::identity()
        } else {
            rand::random::<UnitQuaternion<f32>>()
        };
        let position = Isometry3::from_parts(Translation3::identity(), rotation);

        let intersection = shape
            .cast_ray_and_get_normal(&position, &ray, std::f32::MAX, true)
            .expect(&format!(
                "Ray {:?} did not hit Shape {} rotated with {:?}",
                ray, name, rotation
            ));

        let point = ray.origin + ray.dir * intersection.toi;
        let point_nudged_in = point + position * intersection.normal * -0.001;
        let point_nudged_out = point + position * intersection.normal * 0.001;

        assert!(
            shape.contains_point(&position, &point_nudged_in),
            format!(
                "Shape {} rotated with {:#?} does not contain point nudged in {:#?}",
                name,
                rotation.axis(),
                point_nudged_in
            ),
        );

        assert!(
            !shape.contains_point(&position, &point_nudged_out),
            format!(
                "Shape {} rotated with {:#?} does contains point nudged out {:#?}",
                name,
                rotation.axis(),
                point_nudged_out
            ),
        );

        let new_ray = Ray::new(point_nudged_out, ray_origin - point_nudged_out);

        assert!(
            shape
                .cast_ray_and_get_normal(&position, &new_ray, std::f32::MAX, true)
                .is_none(),
            format!(
                "Ray {:#?} from outside Shape {} rotated with {:#?} did hit at t={}",
                ray,
                name,
                rotation,
                shape
                    .cast_ray_and_get_normal(&position, &new_ray, std::f32::MAX, true)
                    .expect("recurring ray cast produced a different answer")
                    .toi
            )
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
