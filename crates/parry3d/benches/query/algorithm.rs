use na::Point3;
use parry3d::query::gjk::{eps_tol, CSOPoint, VoronoiSimplex};
use test::Bencher;

#[bench]
fn bench_johnson_simplex(bh: &mut Bencher) {
    let a = CSOPoint::single_point(Point3::new(-0.5f32, -0.5, -0.5));
    let b = CSOPoint::single_point(Point3::new(0.0, 0.5, 0.0));
    let c = CSOPoint::single_point(Point3::new(0.5, -0.5, -0.5));
    let d = CSOPoint::single_point(Point3::new(0.0, -0.5, -0.5));

    bh.iter(|| {
        let mut spl = VoronoiSimplex::new();

        spl.reset(a);

        spl.add_point(b, eps_tol());
        spl.add_point(c, eps_tol());
        spl.add_point(d, eps_tol());

        test::black_box(spl.project_origin_and_reduce());
    })
}

#[bench]
fn bench_johnson_simplex_tls(bh: &mut Bencher) {
    let a = CSOPoint::single_point(Point3::new(-0.5f32, -0.5, -0.5));
    let b = CSOPoint::single_point(Point3::new(0.0, 0.5, 0.0));
    let c = CSOPoint::single_point(Point3::new(0.5, -0.5, -0.5));
    let d = CSOPoint::single_point(Point3::new(0.0, -0.5, -0.5));

    bh.iter(|| {
        let mut spl = VoronoiSimplex::new();

        spl.reset(a);

        spl.add_point(b, eps_tol());
        spl.add_point(c, eps_tol());
        spl.add_point(d, eps_tol());

        test::black_box(spl.project_origin_and_reduce());
    })
}
