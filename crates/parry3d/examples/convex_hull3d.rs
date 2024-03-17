use parry3d::math::Point;
use parry3d::transformation;

fn main() {
    let points = vec![
        Point::new(0.77705324, 0.05374551, 0.9822232),
        Point::new(0.35096353, 0.9873069, 0.28922123),
        Point::new(0.09537989, 0.44411153, 0.05486667),
        Point::new(0.108208835, 0.72445065, 0.6669141),
        Point::new(0.7661844, 0.86163324, 0.80507314),
        Point::new(0.5185994, 0.66594696, 0.072779536),
        Point::new(0.768981, 0.23657233, 0.44346774),
        Point::new(0.058607936, 0.09037298, 0.017009139),
        Point::new(0.8818559, 0.3804205, 0.25173646),
        Point::new(0.9571466, 0.17664945, 0.6029223),
    ];

    let _ = transformation::convex_hull(&points[..]);
}
