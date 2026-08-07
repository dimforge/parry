// Regression test for https://github.com/dimforge/rapier/issues/332
//
// `HeightField::project_local_point` used to iterate on every triangle, making
// ball-vs-heightfield tests O(rows * cols) per query; it now visits only the cells
// overlapping a growing neighborhood. These check it matches brute force and is faster.

use parry3d::math::Vector;
use parry3d::query::PointQuery;
use parry3d::shape::HeightField;
use parry3d::utils::Array2;

/// Tiny deterministic LCG so the tests don't depend on rand seeding.
struct Lcg(u64);

impl Lcg {
    fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 32) as u32
    }

    fn real(&mut self, min: f32, max: f32) -> f32 {
        min + (max - min) * (self.next_u32() as f32 / u32::MAX as f32)
    }
}

fn random_heightfield(rng: &mut Lcg, n: usize, scale: Vector) -> HeightField {
    let heights: Vec<f32> = (0..n * n).map(|_| rng.real(-1.0, 1.0)).collect();
    HeightField::new(Array2::new(n, n, heights), scale)
}

fn brute_force_distance(heightfield: &HeightField, pt: Vector) -> f32 {
    let mut smallest_dist = f32::MAX;
    for tri in heightfield.triangles() {
        let proj = tri.project_local_point(pt, false);
        smallest_dist = smallest_dist.min((pt - proj.point).length());
    }
    smallest_dist
}

#[test]
fn pruned_projection_matches_brute_force() {
    let mut rng = Lcg(0x332);
    let heightfield = random_heightfield(&mut rng, 16, Vector::new(16.0, 2.0, 16.0));

    for i in 0..1000 {
        // Mix of points: near the surface, inside the AABB, far outside, above, below.
        let pt = match i % 4 {
            0 => Vector::new(
                rng.real(-8.0, 8.0),
                rng.real(-2.0, 2.0),
                rng.real(-8.0, 8.0),
            ),
            1 => Vector::new(
                rng.real(-8.0, 8.0),
                rng.real(-100.0, 100.0),
                rng.real(-8.0, 8.0),
            ),
            2 => Vector::new(
                rng.real(-100.0, 100.0),
                rng.real(-10.0, 10.0),
                rng.real(-100.0, 100.0),
            ),
            _ => Vector::new(
                rng.real(-8.5, 8.5),
                rng.real(-3.0, 3.0),
                rng.real(-8.5, 8.5),
            ),
        };

        let proj = heightfield.project_local_point(pt, false);
        let dist = (pt - proj.point).length();
        let brute_dist = brute_force_distance(&heightfield, pt);

        assert!(
            (dist - brute_dist).abs() <= 1.0e-4 * brute_dist.max(1.0),
            "point {pt:?}: pruned distance {dist} != brute-force distance {brute_dist}"
        );
    }
}

#[test]
fn pruned_projection_matches_brute_force_flat_field() {
    // The exact setup from the issue: a constant-height field.
    let n = 16;
    let heights = vec![1.0; n * n];
    let heightfield = HeightField::new(Array2::new(n, n, heights), Vector::new(16.0, 1.0, 16.0));
    let mut rng = Lcg(0x332_332);

    for _ in 0..1000 {
        let pt = Vector::new(
            rng.real(-20.0, 20.0),
            rng.real(-5.0, 5.0),
            rng.real(-20.0, 20.0),
        );

        let proj = heightfield.project_local_point(pt, true);
        let dist = (pt - proj.point).length();
        let brute_dist = brute_force_distance(&heightfield, pt);

        assert!(
            (dist - brute_dist).abs() <= 1.0e-4 * brute_dist.max(1.0),
            "point {pt:?}: pruned distance {dist} != brute-force distance {brute_dist}"
        );
    }
}

#[test]
#[ignore = "benchmark, run manually with --ignored --test-threads=1"]
fn pruned_projection_perf() {
    let mut rng = Lcg(0xbe7c4);
    let n = 512;
    let heightfield = random_heightfield(&mut rng, n, Vector::new(256.0, 2.0, 256.0));

    let pts: Vec<Vector> = (0..10_000)
        .map(|_| {
            Vector::new(
                rng.real(-128.0, 128.0),
                rng.real(-4.0, 4.0),
                rng.real(-128.0, 128.0),
            )
        })
        .collect();

    let t_pruned = std::time::Instant::now();
    let mut acc = 0.0;
    for pt in &pts {
        let proj = heightfield.project_local_point(*pt, false);
        acc += (proj.point - *pt).length();
    }
    let t_pruned = t_pruned.elapsed().as_secs_f64() / pts.len() as f64;

    // Brute force is way too slow for 10k queries; sample it on a few points.
    let t_brute = std::time::Instant::now();
    for pt in &pts[..20] {
        acc += brute_force_distance(&heightfield, *pt);
    }
    let t_brute = t_brute.elapsed().as_secs_f64() / 20.0;

    println!(
        "pruned: {:.3}us/query, brute force: {:.3}us/query, speedup: {:.1}x (acc: {acc})",
        t_pruned * 1.0e6,
        t_brute * 1.0e6,
        t_brute / t_pruned
    );

    assert!(
        t_pruned * 100.0 < t_brute,
        "pruned projection should be at least 100x faster than brute force \
         (pruned: {t_pruned}s, brute: {t_brute}s)"
    );
}
