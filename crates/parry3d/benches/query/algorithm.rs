use parry3d::math::Vec3;
use parry3d::query::gjk::{CsoPoint, VoronoiSimplex};
use test::Bencher;

#[bench]
fn bench_johnson_simplex(bh: &mut Bencher) {
    let a = CsoPoint::single_point(Vec3::new(-0.5f32, -0.5, -0.5));
    let b = CsoPoint::single_point(Vec3::new(0.0, 0.5, 0.0));
    let c = CsoPoint::single_point(Vec3::new(0.5, -0.5, -0.5));
    let d = CsoPoint::single_point(Vec3::new(0.0, -0.5, -0.5));

    bh.iter(|| {
        let mut spl = VoronoiSimplex::new();

        spl.reset(a);

        spl.add_point(b);
        spl.add_point(c);
        spl.add_point(d);

        test::black_box(spl.project_origin_and_reduce());
    })
}

#[bench]
fn bench_johnson_simplex_tls(bh: &mut Bencher) {
    let a = CsoPoint::single_point(Vec3::new(-0.5f32, -0.5, -0.5));
    let b = CsoPoint::single_point(Vec3::new(0.0, 0.5, 0.0));
    let c = CsoPoint::single_point(Vec3::new(0.5, -0.5, -0.5));
    let d = CsoPoint::single_point(Vec3::new(0.0, -0.5, -0.5));

    bh.iter(|| {
        let mut spl = VoronoiSimplex::new();

        spl.reset(a);

        spl.add_point(b);
        spl.add_point(c);
        spl.add_point(d);

        test::black_box(spl.project_origin_and_reduce());
    })
}
