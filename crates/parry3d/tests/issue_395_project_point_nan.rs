// Regression test for https://github.com/dimforge/parry/issues/395
//
// `CompositeShapeRef::project_local_point` (and friends) used to panic with
// `unreachable!()` when `Bvh::find_best` recorded no candidate, which is what a NaN
// query point causes. Point queries now report the query point itself instead.

use parry3d::math::{Pose, Real, Vector};
use parry3d::query::{PointQuery, PointQueryWithLocation};
use parry3d::shape::{Compound, Cuboid, Polyline, SharedShape, TriMesh, TriMeshFlags};

fn cube_mesh() -> (Vec<Vector>, Vec<[u32; 3]>) {
    let vertices = vec![
        Vector::new(-0.5, -0.5, -0.5),
        Vector::new(0.5, -0.5, -0.5),
        Vector::new(0.5, 0.5, -0.5),
        Vector::new(-0.5, 0.5, -0.5),
        Vector::new(-0.5, -0.5, 0.5),
        Vector::new(0.5, -0.5, 0.5),
        Vector::new(0.5, 0.5, 0.5),
        Vector::new(-0.5, 0.5, 0.5),
    ];
    let indices = vec![
        [0, 2, 1],
        [0, 3, 2],
        [4, 5, 6],
        [4, 6, 7],
        [0, 1, 5],
        [0, 5, 4],
        [2, 3, 7],
        [2, 7, 6],
        [1, 2, 6],
        [1, 6, 5],
        [0, 4, 7],
        [0, 7, 3],
    ];
    (vertices, indices)
}

#[test]
fn trimesh_nan_point_queries_do_not_panic() {
    let (vertices, indices) = cube_mesh();
    let mesh = TriMesh::new(vertices, indices).unwrap();
    let nan = Vector::splat(Real::NAN);

    // The projection is the NaN query point itself, not an arbitrary point picked
    // from the mesh.
    for solid in [true, false] {
        let proj = mesh.project_local_point(nan, solid);
        assert!(proj.point.is_nan());
        assert!(!proj.is_inside);

        let (proj, _) = mesh.project_local_point_and_get_location(nan, solid);
        assert!(proj.point.is_nan());
    }

    assert!(mesh
        .project_local_point_and_get_feature(nan)
        .0
        .point
        .is_nan());
    assert!(mesh
        .project_local_point_with_max_dist(nan, true, Real::MAX)
        .is_none());
    assert!(mesh.distance_to_local_point(nan, true).is_nan());
    assert!(!mesh.contains_local_point(nan));
}

#[test]
fn oriented_trimesh_nan_point_queries_do_not_panic() {
    // The pseudo-normals code path is different; make sure it is NaN-safe too.
    let (vertices, indices) = cube_mesh();
    let mesh = TriMesh::with_flags(vertices, indices, TriMeshFlags::ORIENTED).unwrap();
    let nan = Vector::splat(Real::NAN);

    assert!(mesh.project_local_point(nan, true).point.is_nan());
    assert!(mesh
        .project_local_point_and_get_feature(nan)
        .0
        .point
        .is_nan());
    assert!(mesh
        .project_local_point_and_get_location(nan, true)
        .0
        .point
        .is_nan());
    assert!(!mesh.contains_local_point(nan));
}

#[test]
fn polyline_nan_point_queries_do_not_panic() {
    let polyline = Polyline::new(
        vec![
            Vector::new(0.0, 0.0, 0.0),
            Vector::new(1.0, 0.0, 0.0),
            Vector::new(1.0, 1.0, 0.0),
        ],
        None,
    );
    let nan = Vector::splat(Real::NAN);

    assert!(polyline.project_local_point(nan, true).point.is_nan());
    assert!(polyline
        .project_local_point_and_get_feature(nan)
        .0
        .point
        .is_nan());
    assert!(polyline
        .project_local_point_and_get_location(nan, true)
        .0
        .point
        .is_nan());
    assert!(!polyline.contains_local_point(nan));
}

#[test]
fn compound_nan_point_queries_do_not_panic() {
    let compound = Compound::new(vec![
        (
            Pose::IDENTITY,
            SharedShape::new(Cuboid::new(Vector::splat(0.5))),
        ),
        (
            Pose::from_translation(Vector::new(2.0, 0.0, 0.0)),
            SharedShape::new(Cuboid::new(Vector::splat(0.5))),
        ),
    ]);
    let nan = Vector::splat(Real::NAN);

    for solid in [true, false] {
        assert!(compound.project_local_point(nan, solid).point.is_nan());
    }
    assert!(compound
        .project_local_point_and_get_feature(nan)
        .0
        .point
        .is_nan());
    // NOTE: unlike `TriMesh`/`Polyline` above, this currently returns `true`, because
    //       `Aabb::contains_local_point` rejects a point by testing whether it lies
    //       outside the bounds: every such test is false for a NaN, so the point is
    //       reported as contained. That affects `Aabb`/`Cuboid` themselves, not just
    //       composite shapes, so it is left alone here.
    let _ = compound.contains_local_point(nan);
}

#[test]
fn finite_point_projection_still_works() {
    // Sanity check: the non-finite handling must not change regular queries.
    let (vertices, indices) = cube_mesh();
    let mesh = TriMesh::new(vertices, indices).unwrap();

    let proj = mesh.project_local_point(Vector::new(2.0, 0.0, 0.0), true);
    assert!(!proj.is_inside);
    assert!((proj.point - Vector::new(0.5, 0.0, 0.0)).length() < 1.0e-5);
}

#[test]
fn empty_trimesh_is_a_construction_error() {
    assert!(TriMesh::new(Vec::new(), Vec::new()).is_err());
}

#[test]
#[should_panic(expected = "must contain at least one shape")]
fn empty_compound_panics_at_construction() {
    let _ = Compound::new(Vec::new());
}
