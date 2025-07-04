//! The EPA algorithm for penetration depth computation.
//!
#[cfg(feature = "dim2")]
pub use self::epa2::EPA;
#[cfg(feature = "dim3")]
pub use self::epa3::EPA;

#[cfg(feature = "dim2")]
pub mod epa2;
#[cfg(feature = "dim3")]
pub mod epa3;

#[cfg(test)]
mod test {
    use crate::math::Real;
    use crate::{math::Point, shape::Triangle};

    fn point_dim(x: Real, y: Real) -> Point<Real> {
        Point::new(
            x,
            y,
            #[cfg(feature = "dim3")]
            0.0,
        )
    }
    #[test]
    fn same_2_points_triangle() {
        use na::Isometry;

        use crate::query::PointQueryWithLocation;

        let triangle = Triangle::new(
            point_dim(40.0, 0.0),
            point_dim(0.0, 80.0),
            point_dim(0.0, 80.0),
        );
        let res = triangle.project_point_and_get_location(
            &Isometry::identity(),
            &point_dim(10.0, 20.0),
            false,
        );
        res.0.point.iter().for_each(|p| assert!(p.is_finite()));
        let res = triangle.project_point_and_get_location(
            &Isometry::identity(),
            &point_dim(40.0, 0.0),
            false,
        );
        res.0.point.iter().for_each(|p| assert!(p.is_finite()));
    }
    #[test]
    fn collinear_points_triangle() {
        use na::Isometry;

        use crate::query::PointQueryWithLocation;

        let triangle = Triangle::new(
            point_dim(0.0, 0.0),
            point_dim(100.0, 0.0),
            point_dim(160.0, 0.0),
        );
        let res = triangle.project_point_and_get_location(
            &Isometry::identity(),
            &point_dim(10.0, 0.0),
            false,
        );
        assert!(res.0.is_inside);
        res.0.point.iter().for_each(|p| assert!(p.is_finite()));
        let res = triangle.project_point_and_get_location(
            &Isometry::identity(),
            &point_dim(10.0, 10.0),
            false,
        );
        // FIXME: This assert is currently failing !
        // False negative for inside on edge may be understandable due to floating points imprecisions,
        // But false positives for points outside is a problem. (maybe related to triangle orientation?)
        assert!(res.0.is_inside == false);
        res.0.point.iter().for_each(|p| assert!(p.is_finite()));
    }
}
