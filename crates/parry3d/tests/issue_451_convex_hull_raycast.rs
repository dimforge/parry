use parry3d::math::Vector;
use parry3d::query::{PointQueryWithLocation, Ray};
use parry3d::shape::{SharedShape, Tetrahedron, TetrahedronPointLocation};

#[test]
fn tetrahedron_distant_face_projection() {
    let a = Vector::ZERO;
    let b = Vector::new(1.0, -1.0, 0.0);
    let c = Vector::new(1.0, 0.0, -1.0);
    let d = Vector::splat(-1.0);
    let expected = Vector::new(0.5, -0.25, -0.25);
    for (face, tetrahedron) in [
        Tetrahedron::new(a, b, c, d),
        Tetrahedron::new(a, b, d, c),
        Tetrahedron::new(a, d, b, c),
        Tetrahedron::new(d, a, b, c),
    ]
    .into_iter()
    .enumerate()
    {
        for distance in [0.25, 4096.0] {
            let (projection, location) = tetrahedron
                .project_local_point_and_get_location(expected + Vector::splat(distance), true);
            assert!(!projection.is_inside);
            assert!((projection.point - expected).length() < 1.0e-3);
            let TetrahedronPointLocation::OnFace(index, coordinates) = location else {
                panic!("unexpected location: {location:?}");
            };
            assert_eq!(index as usize, face);
            for (actual, expected) in coordinates.into_iter().zip([0.5, 0.25, 0.25]) {
                assert!((actual - expected).abs() < 1.0e-3);
            }
        }
    }
}

fn hull() -> SharedShape {
    SharedShape::convex_hull(&[
        Vector::new(-34.0, -14.0, 29.0),
        Vector::new(-24.0, -12.0, 29.0),
        Vector::new(-32.0, -12.0, 29.0),
        Vector::new(-32.0, -12.0, 2.0),
        Vector::new(-24.0, -12.0, 2.0),
    ])
    .unwrap()
}

#[test]
fn convex_hull_nonintersecting_ray() {
    let ray = Ray::new(
        Vector::new(-16.0, 8.0, 0.0),
        Vector::new(-0.5956519, -0.67131174, 0.44106603),
    );
    assert!(!hull().intersects_local_ray(&ray, f32::MAX));
}

#[test]
fn convex_hull_intersecting_ray() {
    let ray = Ray::new(Vector::new(-28.0, 0.0, 10.0), -Vector::Y);
    assert!(hull().intersects_local_ray(&ray, f32::MAX));
}
