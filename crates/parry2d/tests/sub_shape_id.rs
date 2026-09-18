use parry2d::math::{Pose, Real, Vector};
use parry2d::query::{self, PointQuery, Ray, RayCast};
use parry2d::shape::{Ball, Compound, Cuboid, Polyline, SharedShape};

/// Three unit boxes in a row along x, centered at x = 0, 4 and 8.
fn three_boxes() -> Compound {
    Compound::new(
        (0..3)
            .map(|i| {
                (
                    Pose::from_translation(Vector::new(i as Real * 4.0, 0.0)),
                    SharedShape::new(Cuboid::new(Vector::splat(0.5))),
                )
            })
            .collect(),
    )
}

#[test]
fn queries_against_a_compound_report_the_part() {
    let compound = three_boxes();
    let probe = Ball::new(0.25);

    for part in 0..3u32 {
        let x = part as Real * 4.0;

        let ray = Ray::new(Vector::new(x, 5.0), Vector::new(0.0, -1.0));
        assert_eq!(
            compound
                .cast_local_ray_and_get_normal(&ray, Real::MAX, true)
                .expect("hits a box")
                .subshape,
            part
        );

        assert_eq!(
            compound
                .project_local_point(Vector::new(x, 3.0), false)
                .subshape,
            part
        );

        let pose12 = Pose::from_translation(Vector::new(x, 0.7));
        let contact = query::contact(&Pose::IDENTITY, &compound, &pose12, &probe, 1.0)
            .unwrap()
            .expect("within prediction");
        assert_eq!((contact.subshape1, contact.subshape2), (part, 0));

        let dist = query::distance(&Pose::IDENTITY, &compound, &pose12, &probe).unwrap();
        assert_eq!((dist.subshape1, dist.subshape2), (part, 0));

        let overlapping = Pose::from_translation(Vector::new(x, 0.0));
        let test =
            query::intersection_test(&Pose::IDENTITY, &compound, &overlapping, &probe).unwrap();
        assert!(test.intersecting);
        assert_eq!((test.subshape1, test.subshape2), (part, 0));
    }
}

#[test]
fn a_polyline_reports_the_segment() {
    // A staircase: segment 0 spans x in [0,1], segment 1 [1,2], segment 2 [2,3].
    let polyline = Polyline::new(
        vec![
            Vector::new(0.0, 0.0),
            Vector::new(1.0, 0.0),
            Vector::new(2.0, 0.0),
            Vector::new(3.0, 0.0),
        ],
        None,
    );

    for segment in 0..3u32 {
        let x = segment as Real + 0.5;
        assert_eq!(
            polyline
                .project_local_point(Vector::new(x, 2.0), false)
                .subshape,
            segment
        );

        let ray = Ray::new(Vector::new(x, 2.0), Vector::new(0.0, -1.0));
        assert_eq!(
            polyline
                .cast_local_ray_and_get_normal(&ray, Real::MAX, true)
                .expect("hits the polyline")
                .subshape,
            segment
        );
    }
}

/// The conversion resolves a segment endpoint to the polyline vertex it indexes, so the two
/// segments meeting at a corner name the same vertex.
#[test]
fn segment_features_convert_to_polyline_features() {
    use parry2d::shape::FeatureId;

    // Segment 0 spans vertices 0-1, segment 1 spans 1-2: they share vertex 1.
    let polyline = Polyline::new(
        vec![
            Vector::new(0.0, 0.0),
            Vector::new(1.0, 0.0),
            Vector::new(2.0, 0.0),
        ],
        None,
    );

    // Endpoint 1 of segment 0 and endpoint 0 of segment 1 are the same polyline vertex.
    assert_eq!(
        polyline.segment_feature_to_polyline_feature(0, FeatureId::Vertex(1)),
        FeatureId::Vertex(1)
    );
    assert_eq!(
        polyline.segment_feature_to_polyline_feature(1, FeatureId::Vertex(0)),
        FeatureId::Vertex(1)
    );
    assert_eq!(
        polyline.segment_feature_to_polyline_feature(1, FeatureId::Vertex(1)),
        FeatureId::Vertex(2)
    );

    // Each segment has a side facing each way, and every (segment, side) pair is distinct.
    let sides: Vec<_> = (0..2u32)
        .flat_map(|segment| (0..2u32).map(move |side| (segment, side)))
        .map(|(segment, side)| {
            polyline.segment_feature_to_polyline_feature(segment, FeatureId::Face(side))
        })
        .collect();
    for (i, face) in sides.iter().enumerate() {
        assert!(!sides[..i].contains(face), "{face:?} reused");
    }
}

#[test]
fn shapes_without_sub_shapes_report_zero() {
    let ball = Ball::new(1.0);
    let ray = Ray::new(Vector::new(0.0, 5.0), Vector::new(0.0, -1.0));

    assert_eq!(
        ball.cast_local_ray_and_get_normal(&ray, Real::MAX, true)
            .unwrap()
            .subshape,
        0
    );
    assert_eq!(
        ball.project_local_point(Vector::new(3.0, 0.0), false)
            .subshape,
        0
    );

    let dist = query::distance(
        &Pose::IDENTITY,
        &ball,
        &Pose::from_translation(Vector::new(5.0, 0.0)),
        &Ball::new(1.0),
    )
    .unwrap();
    assert_eq!((dist.subshape1, dist.subshape2), (0, 0));
}
