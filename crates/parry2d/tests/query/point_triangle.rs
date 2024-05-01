use parry2d::{math::Point, query::PointQuery, shape::Triangle};

#[test]
fn project_local_point_point_on_ab() {
    let verts = [Point::new(2.0, 1.0), Point::new(0.0, 1.0), Point::new(1.0, 0.0)];
    let tri1 = Triangle::new(verts[0], verts[1], verts[2]);
    let tri2 = Triangle::new(verts[2], verts[0], verts[1]);

    let query_pt = Point::new(1.4, 1.0);

    let proj1 = tri1.project_local_point(&query_pt, false); // Used to fail on 0.14 and earlier
    let proj2 = tri2.project_local_point(&query_pt, false);

    assert_eq!(proj1.point, proj2.point);
    assert_eq!(proj1.point, query_pt);
}
