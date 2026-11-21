// Regression test for https://github.com/dimforge/parry/issues/473
//
// Cylinders used the generic GJK ray cast, which converges slowly near the cap rim: it could
// return `None` (a miss) at its iteration cap, or a time of impact off the surface.

use parry3d::math::Vector;
use parry3d::query::{Ray, RayCast};
use parry3d::shape::Cylinder;

#[test]
fn ray_from_far_above_hits_the_cap() {
    let cylinder = Cylinder::new(0.5, 0.5);
    let ray = Ray::new(Vector::new(0.0, 1000.0, 0.0), Vector::new(0.0, -1.0, 0.0));
    let hit = cylinder
        .cast_local_ray_and_get_normal(&ray, f32::MAX, true)
        .expect("the ray goes through the middle of the cap");
    assert!((hit.time_of_impact - 999.5).abs() < 1.0e-3);
    assert_eq!(hit.normal, Vector::new(0.0, 1.0, 0.0));
}

#[test]
fn rays_crossing_the_cylinder_hit_from_both_sides() {
    let (half_height, radius) = (5.0, 5.0);
    let cylinder = Cylinder::new(half_height, radius);
    // Vertical rays through the caps, close to the rim, cast from above and from below.
    for i in 0..200 {
        let x = radius * (i as f32 / 200.0) * 0.999;
        for dir in [-1.0, 1.0] {
            let origin = Vector::new(x, -dir * 20.0, 0.0);
            let ray = Ray::new(origin, Vector::new(0.0, dir, 0.0));
            let hit = cylinder
                .cast_local_ray_and_get_normal(&ray, f32::MAX, true)
                .unwrap_or_else(|| panic!("missed at x = {x}"));
            assert!((hit.time_of_impact - (20.0 - half_height)).abs() < 1.0e-4);
            assert_eq!(hit.normal, Vector::new(0.0, -dir, 0.0));
        }
    }
    // Horizontal rays through the side, cast from both sides.
    for i in 0..200 {
        let z = radius * (i as f32 / 200.0) * 0.999;
        for dir in [-1.0, 1.0] {
            let origin = Vector::new(-dir * 20.0, 4.9, z);
            let ray = Ray::new(origin, Vector::new(dir, 0.0, 0.0));
            let hit = cylinder
                .cast_local_ray_and_get_normal(&ray, f32::MAX, true)
                .unwrap_or_else(|| panic!("missed at z = {z}"));
            let expected = 20.0 - (radius * radius - z * z).sqrt();
            assert!((hit.time_of_impact - expected).abs() < 1.0e-3);
        }
    }
}

#[test]
fn ray_starting_inside_a_hollow_cylinder_hits_on_its_way_out() {
    let cylinder = Cylinder::new(1.0, 2.0);
    let ray = Ray::new(Vector::ZERO, Vector::new(1.0, 0.0, 0.0));
    let solid = cylinder.cast_local_ray(&ray, f32::MAX, true);
    assert_eq!(solid, Some(0.0));
    let hit = cylinder
        .cast_local_ray_and_get_normal(&ray, f32::MAX, false)
        .unwrap();
    assert!((hit.time_of_impact - 2.0).abs() < 1.0e-6);
    assert_eq!(hit.normal, Vector::new(-1.0, 0.0, 0.0));
}

#[test]
fn ray_from_very_far_gets_a_unit_normal() {
    // From 200 km a 5 mm radius is below the float's resolution at the origin: the hit point
    // must still not land on the axis, where the side normal is undefined.
    let cylinder = Cylinder::new(1.0, 0.005);
    let ray = Ray::new(
        Vector::new(200_000.0, 0.0, 0.0),
        Vector::new(-1.0, 0.0, 0.0),
    );
    let hit = cylinder
        .cast_local_ray_and_get_normal(&ray, f32::MAX, true)
        .expect("the ray goes through the cylinder");
    assert_eq!(hit.normal, Vector::new(1.0, 0.0, 0.0));
}

#[test]
fn ray_without_direction_misses() {
    let cylinder = Cylinder::new(1.0, 1.0);
    for origin in [
        Vector::ZERO,
        Vector::new(0.2, 0.1, 0.1),
        Vector::new(5.0, 0.0, 0.0),
    ] {
        for solid in [true, false] {
            let ray = Ray::new(origin, Vector::ZERO);
            assert!(cylinder
                .cast_local_ray_and_get_normal(&ray, f32::MAX, solid)
                .is_none());
        }
    }
}
