use parry3d::math::{Pose, Real, Vector};
use parry3d::query::{self, PointQuery, Ray, RayCast, ShapeCastOptions};
use parry3d::shape::{Ball, Compound, Cuboid, SharedShape, TriMesh};

/// Three unit boxes in a row along x, centered at x = 0, 4 and 8, so a query aimed at one of them
/// can only be answered by that part.
fn three_boxes() -> Compound {
    Compound::new(
        (0..3)
            .map(|i| {
                (
                    Pose::from_translation(Vector::new(i as Real * 4.0, 0.0, 0.0)),
                    SharedShape::new(Cuboid::new(Vector::splat(0.5))),
                )
            })
            .collect(),
    )
}

#[test]
fn ray_cast_reports_the_part_it_hit() {
    let compound = three_boxes();

    for part in 0..3u32 {
        let ray = Ray::new(
            Vector::new(part as Real * 4.0, 5.0, 0.0),
            Vector::new(0.0, -1.0, 0.0),
        );
        let hit = compound
            .cast_local_ray_and_get_normal(&ray, Real::MAX, true)
            .expect("hits a box");
        assert_eq!(hit.subshape, part, "ray aimed at part {part}");
    }
}

#[test]
fn point_projection_reports_the_part_it_projects_onto() {
    let compound = three_boxes();

    for part in 0..3u32 {
        let point = Vector::new(part as Real * 4.0, 3.0, 0.0);
        assert_eq!(compound.project_local_point(point, false).subshape, part);
        assert_eq!(
            compound
                .project_local_point_and_get_feature(point)
                .0
                .subshape,
            part
        );
    }
}

#[test]
fn contact_shape_cast_distance_and_intersection_report_the_part() {
    let compound = three_boxes();
    let probe = Ball::new(0.25);

    for part in 0..3u32 {
        // A ball resting just above the top face of one box.
        let pose12 = Pose::from_translation(Vector::new(part as Real * 4.0, 0.7, 0.0));

        let contact = query::contact(&Pose::IDENTITY, &compound, &pose12, &probe, 1.0)
            .unwrap()
            .expect("within prediction");
        assert_eq!(contact.subshape1, part, "contact");
        assert_eq!(contact.subshape2, 0, "the ball has no sub-shapes");

        let dist = query::distance(&Pose::IDENTITY, &compound, &pose12, &probe).unwrap();
        assert_eq!(dist.subshape1, part, "distance");
        assert_eq!(dist.subshape2, 0);

        // Cast the ball straight down into that box.
        let hit = query::cast_shapes(
            &Pose::IDENTITY,
            Vector::ZERO,
            &compound,
            &Pose::from_translation(Vector::new(part as Real * 4.0, 5.0, 0.0)),
            Vector::new(0.0, -1.0, 0.0),
            &probe,
            ShapeCastOptions::default(),
        )
        .unwrap()
        .expect("the ball reaches the box");
        assert_eq!(hit.subshape1, part, "shape cast");

        // Overlapping the box outright.
        let overlapping = Pose::from_translation(Vector::new(part as Real * 4.0, 0.0, 0.0));
        let test =
            query::intersection_test(&Pose::IDENTITY, &compound, &overlapping, &probe).unwrap();
        assert!(test.intersecting);
        assert_eq!(test.subshape1, part, "intersection test");
    }
}

#[test]
fn the_roles_swap_when_the_composite_is_the_second_shape() {
    let compound = three_boxes();
    let probe = Ball::new(0.25);
    let part = 2u32;
    let pose1 = Pose::from_translation(Vector::new(part as Real * 4.0, 0.7, 0.0));

    // Same query with the composite as shape 2: the id must land on `subshape2`.
    let contact = query::contact(&pose1, &probe, &Pose::IDENTITY, &compound, 1.0)
        .unwrap()
        .expect("within prediction");
    assert_eq!(contact.subshape1, 0, "the ball has no sub-shapes");
    assert_eq!(contact.subshape2, part);

    let dist = query::distance(&pose1, &probe, &Pose::IDENTITY, &compound).unwrap();
    assert_eq!(dist.subshape1, 0);
    assert_eq!(dist.subshape2, part);

    let overlapping = Pose::from_translation(Vector::new(part as Real * 4.0, 0.0, 0.0));
    let test = query::intersection_test(&overlapping, &probe, &Pose::IDENTITY, &compound).unwrap();
    assert!(test.intersecting);
    assert_eq!(test.subshape1, 0);
    assert_eq!(test.subshape2, part);
}

#[test]
fn a_trimesh_reports_the_triangle() {
    // Two triangles forming a quad in the xz plane; each query lands on a known one.
    let vertices = vec![
        Vector::new(0.0, 0.0, 0.0),
        Vector::new(1.0, 0.0, 0.0),
        Vector::new(1.0, 0.0, 1.0),
        Vector::new(0.0, 0.0, 1.0),
    ];
    let mesh = TriMesh::new(vertices, vec![[0, 1, 2], [0, 2, 3]]).unwrap();

    // A point over the first triangle's interior, and one over the second's.
    let over_first = Vector::new(0.8, 1.0, 0.4);
    let over_second = Vector::new(0.2, 1.0, 0.6);
    assert_eq!(mesh.project_local_point(over_first, false).subshape, 0);
    assert_eq!(mesh.project_local_point(over_second, false).subshape, 1);

    let ray = Ray::new(over_first, Vector::new(0.0, -1.0, 0.0));
    let hit = mesh
        .cast_local_ray_and_get_normal(&ray, Real::MAX, true)
        .expect("hits the mesh");
    assert_eq!(hit.subshape, 0);
}

/// The feature a query reports is the sub-shape's own; the sub-shape itself is `subshape`.
#[test]
fn features_are_local_to_the_sub_shape() {
    use parry3d::shape::FeatureId;

    // A quad of two triangles in the xz plane, wound so its front faces up.
    let mesh = TriMesh::new(
        vec![
            Vector::new(0.0, 0.0, 0.0),
            Vector::new(1.0, 0.0, 0.0),
            Vector::new(1.0, 0.0, 1.0),
            Vector::new(0.0, 0.0, 1.0),
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
    .unwrap();

    // Hitting triangle 1 from above reports one of the triangle's own two faces, never its index.
    let from_above = Ray::new(Vector::new(0.2, 1.0, 0.6), Vector::new(0.0, -1.0, 0.0));
    let hit = mesh
        .cast_local_ray_and_get_normal(&from_above, Real::MAX, true)
        .expect("hits the mesh");
    assert_eq!(hit.subshape, 1);
    assert!(matches!(
        hit.feature,
        FeatureId::Face(0) | FeatureId::Face(1)
    ));

    // The two sides of a triangle are told apart by the feature alone.
    let from_below = Ray::new(Vector::new(0.2, -1.0, 0.6), Vector::new(0.0, 1.0, 0.0));
    let below = mesh
        .cast_local_ray_and_get_normal(&from_below, Real::MAX, true)
        .expect("hits the mesh");
    assert_eq!(below.subshape, 1);
    assert_ne!(
        mesh.is_backface(hit.feature),
        mesh.is_backface(below.feature),
        "one of the two sides has to be the backface"
    );
}

/// `feature_normal_at_point` needs the sub-shape to pick the triangle: the feature alone is the
/// triangle's own and says nothing about which one it is.
#[test]
fn a_trimesh_normal_follows_the_sub_shape() {
    use parry3d::shape::{FeatureId, Shape};

    // Two triangles with clearly different normals: one flat, one tilted.
    let mesh = TriMesh::new(
        vec![
            Vector::new(0.0, 0.0, 0.0),
            Vector::new(1.0, 0.0, 0.0),
            Vector::new(1.0, 0.0, 1.0),
            Vector::new(0.0, 1.0, 1.0),
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
    .unwrap();

    for triangle in 0..2u32 {
        let expected = mesh.triangle(triangle).normal().unwrap();
        let normal = mesh
            .feature_normal_at_point(triangle, FeatureId::Face(0), Vector::ZERO)
            .expect("the sub-shape names the triangle");
        assert_eq!(mesh.triangle_normal(triangle), Some(normal));
        assert!(
            (normal - expected).length() < 1.0e-5,
            "triangle {triangle}: {normal:?} vs {expected:?}"
        );
    }
}

#[test]
fn shapes_without_sub_shapes_report_zero() {
    let ball = Ball::new(1.0);
    let cuboid = Cuboid::new(Vector::splat(1.0));

    let ray = Ray::new(Vector::new(0.0, 5.0, 0.0), Vector::new(0.0, -1.0, 0.0));
    assert_eq!(
        ball.cast_local_ray_and_get_normal(&ray, Real::MAX, true)
            .unwrap()
            .subshape,
        0
    );
    assert_eq!(
        ball.project_local_point(Vector::new(3.0, 0.0, 0.0), false)
            .subshape,
        0
    );

    let pose12 = Pose::from_translation(Vector::new(1.5, 0.0, 0.0));
    let contact = query::contact(&Pose::IDENTITY, &ball, &pose12, &cuboid, 1.0)
        .unwrap()
        .expect("within prediction");
    assert_eq!((contact.subshape1, contact.subshape2), (0, 0));
}
