// Regression test for https://github.com/dimforge/rapier/issues/961
//
// The binned BVH builder computed bin indices without clamping them to the bin count, so
// degenerate leaf AABBs (huge, zero-extent, or non-finite) could panic out-of-bounds.

use parry3d::bounding_volume::Aabb;
use parry3d::math::Vector;
use parry3d::partitioning::{Bvh, BvhBuildStrategy, BvhWorkspace};

/// Tiny deterministic LCG so the fuzz-style tests don't depend on rand seeding.
struct Lcg(u64);

impl Lcg {
    fn next_u32(&mut self) -> u32 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 32) as u32
    }

    fn real(&mut self, min: f32, max: f32) -> f32 {
        min + (max - min) * (self.next_u32() as f32 / u32::MAX as f32)
    }
}

fn aabb(mins: Vector, maxs: Vector) -> Aabb {
    Aabb::new(mins, maxs)
}

fn exercise(leaves: &[Aabb]) {
    // Build (binned) directly from the leaves.
    let mut bvh = Bvh::from_leaves(BvhBuildStrategy::Binned, leaves);
    let mut workspace = BvhWorkspace::default();

    // Explicit full rebuilds.
    bvh.rebuild(&mut workspace, BvhBuildStrategy::Binned);

    // Insert/remove churn + incremental optimization, mimicking the
    // broad-phase update loop from the original report.
    for k in 0..8 {
        for i in 0..leaves.len() {
            if (i + k) % 3 == 0 {
                bvh.remove(i as u32);
            }
        }
        for (i, aabb) in leaves.iter().enumerate() {
            if (i + k) % 3 == 0 {
                bvh.insert(*aabb, i as u32);
            }
        }
        bvh.optimize_incremental(&mut workspace);
    }
}

#[test]
fn coincident_aabbs() {
    let unit = aabb(Vector::new(-1.0, -1.0, -1.0), Vector::new(1.0, 1.0, 1.0));
    exercise(&[unit; 32]);
}

#[test]
fn one_enormous_aabb() {
    let mut leaves = vec![];
    for i in 0..31 {
        let c = Vector::new(i as f32, 0.0, 0.0);
        leaves.push(aabb(c - Vector::splat(0.5), c + Vector::splat(0.5)));
    }
    leaves.push(aabb(Vector::splat(-1.0e30), Vector::splat(1.0e30)));
    exercise(&leaves);
}

#[test]
fn huge_translated_aabbs() {
    // AABBs whose coordinates are large enough for `mins + maxs` to overflow
    // to infinity when computing centers.
    let mut leaves = vec![];
    for i in 0..16 {
        let c = Vector::splat(2.0e38) + Vector::new(i as f32, 0.0, 0.0);
        leaves.push(aabb(c - Vector::splat(0.5), c + Vector::splat(0.5)));
    }
    for i in 0..16 {
        let c = Vector::new(i as f32, 0.0, 0.0);
        leaves.push(aabb(c - Vector::splat(0.5), c + Vector::splat(0.5)));
    }
    exercise(&leaves);
}

#[test]
fn zero_extent_aabbs() {
    let mut leaves = vec![];
    for i in 0..32 {
        let c = Vector::new(i as f32 * 0.25, -(i as f32), 3.0);
        leaves.push(aabb(c, c));
    }
    exercise(&leaves);
}

#[test]
fn non_finite_aabbs() {
    let mut leaves = vec![];
    for i in 0..16 {
        let c = Vector::new(i as f32, 0.0, 0.0);
        leaves.push(aabb(c - Vector::splat(0.5), c + Vector::splat(0.5)));
    }
    leaves.push(aabb(Vector::splat(f32::NAN), Vector::splat(f32::NAN)));
    leaves.push(aabb(Vector::splat(f32::INFINITY), Vector::splat(f32::INFINITY)));
    leaves.push(aabb(Vector::splat(f32::NEG_INFINITY), Vector::splat(f32::INFINITY)));
    leaves.push(Aabb::new_invalid());
    exercise(&leaves);
}

#[test]
fn fuzz_mixed_degenerate_aabbs() {
    let mut rng = Lcg(0x961_961_961);

    for _ in 0..64 {
        let n = 8 + (rng.next_u32() % 64) as usize;
        let mut leaves = vec![];

        for _ in 0..n {
            let kind = rng.next_u32() % 5;
            let leaf = match kind {
                0 => {
                    // Regular AABB.
                    let c = Vector::new(
                        rng.real(-100.0, 100.0),
                        rng.real(-100.0, 100.0),
                        rng.real(-100.0, 100.0),
                    );
                    aabb(c - Vector::splat(0.5), c + Vector::splat(0.5))
                }
                1 => {
                    // Huge AABB.
                    let e = rng.real(1.0e30, 3.0e38);
                    aabb(Vector::splat(-e), Vector::splat(e))
                }
                2 => {
                    // Zero-extent AABB, possibly at a huge coordinate.
                    let c = Vector::splat(rng.real(-3.0e38, 3.0e38));
                    aabb(c, c)
                }
                3 => {
                    // Non-finite AABB.
                    let v = if rng.next_u32() % 2 == 0 {
                        f32::NAN
                    } else {
                        f32::INFINITY
                    };
                    aabb(Vector::splat(v), Vector::splat(v))
                }
                _ => {
                    // Tiny cluster: nearly-coincident centroids.
                    let c = Vector::new(1.0, 2.0, 3.0);
                    let e = rng.real(0.0, 1.0e-30);
                    aabb(c - Vector::splat(e), c + Vector::splat(e))
                }
            };
            leaves.push(leaf);
        }

        exercise(&leaves);
    }
}
