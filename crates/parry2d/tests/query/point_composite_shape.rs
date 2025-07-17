use na::Point2;
use parry2d::{query::PointQuery, shape::TriMesh};

#[test]
fn project_local_point_and_get_feature_gets_the_enclosing_triangle() {
    let vertices = vec![
        Point2::new(0.0, 1.0),
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
    ];

    let mesh = TriMesh::new(vertices, vec![[0, 1, 2], [3, 0, 2]]).unwrap();
    let options = &();
    let query_pt = Point2::new(0.6, 0.6); // Inside the top-right triangle (index 1)

    let (proj, feat) = mesh.project_local_point_and_get_feature(&query_pt, options);

    let correct_tri_idx = 1;
    let correct_tri = mesh.triangle(correct_tri_idx);

    let is_inside_correct = correct_tri.contains_local_point(&query_pt, options);

    assert!(is_inside_correct);
    assert_eq!(proj.is_inside, is_inside_correct);
    assert_eq!(feat.unwrap_face(), correct_tri_idx);
}

#[test]
fn project_local_point_and_get_feature_projects_correctly_from_outside() {
    let vertices = vec![
        Point2::new(0.0, 1.0),
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
    ];

    let mesh = TriMesh::new(vertices, vec![[0, 1, 2], [3, 0, 2]]).unwrap();
    let options = &();
    {
        let query_pt = Point2::new(-1.0, 0.0); // Left from the bottom-left triangle (index 0)

        let (proj, feat) = mesh.project_local_point_and_get_feature(&query_pt, options);

        let correct_tri_idx = 0;
        let correct_tri = mesh.triangle(correct_tri_idx);

        let is_inside_correct = correct_tri.contains_local_point(&query_pt, options);

        assert_eq!(is_inside_correct, false);
        assert_eq!(proj.is_inside, is_inside_correct);
        assert_eq!(proj.point, Point2::origin());
        assert_eq!(feat.unwrap_face(), correct_tri_idx);
    }
    {
        let query_pt = Point2::new(0.5, 2.0); // Above the top-right triangle (index 1)

        let (proj, feat) = mesh.project_local_point_and_get_feature(&query_pt, options);

        let correct_tri_idx = 1;
        let correct_tri = mesh.triangle(correct_tri_idx);

        let is_inside_correct = correct_tri.contains_local_point(&query_pt, options);

        assert_eq!(is_inside_correct, false);
        assert_eq!(proj.is_inside, is_inside_correct);
        assert_eq!(proj.point, Point2::new(0.5, 1.0));
        assert_eq!(feat.unwrap_face(), correct_tri_idx);
    }
}
