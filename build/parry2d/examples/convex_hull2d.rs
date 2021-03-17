extern crate nalgebra as na;

use na::Point2;
use parry2d::transformation;

fn main() {
    let points = vec![
        Point2::new(0.77705324, 0.05374551),
        Point2::new(0.35096353, 0.9873069),
        Point2::new(0.09537989, 0.44411153),
        Point2::new(0.108208835, 0.72445065),
        Point2::new(0.7661844, 0.86163324),
        Point2::new(0.5185994, 0.66594696),
        Point2::new(0.768981, 0.23657233),
        Point2::new(0.058607936, 0.09037298),
        Point2::new(0.8818559, 0.3804205),
        Point2::new(0.9571466, 0.17664945),
    ];

    let _ = transformation::convex_hull(&points[..]);
}
