// Regression test for https://github.com/dimforge/parry/issues/399
//
// The request was a reusable BVH: rebuilding from moved leaves every tick without
// reallocating. This locks in the recipe — `insert_or_update_partially` + `rebuild` with
// a reused `BvhWorkspace` — and checks queries stay correct and the leaf set is stable.

use parry3d::bounding_volume::{Aabb, BoundingVolume};
use parry3d::math::Vector;
use parry3d::partitioning::{Bvh, BvhBuildStrategy, BvhWorkspace};

fn leaf_aabb(i: usize, tick: usize) -> Aabb {
    // Leaves on a moving diagonal so every tick changes every AABB.
    let offset = tick as f32 * 3.5;
    let x = i as f32 * 2.0 + offset;
    let y = (i % 7) as f32 - offset;
    let mins = Vector::new(x, y, -1.0);
    Aabb::new(mins, mins + Vector::new(1.0, 1.0, 2.0))
}

#[test]
fn bvh_tick_loop_rebuild_with_reused_workspace() {
    const NUM_LEAVES: usize = 100;
    const NUM_TICKS: usize = 4;

    let mut bvh = Bvh::new();
    let mut workspace = BvhWorkspace::default();

    for tick in 0..NUM_TICKS {
        // Update (or insert, on the first tick) every leaf in place.
        for i in 0..NUM_LEAVES {
            bvh.insert_or_update_partially(leaf_aabb(i, tick), i as u32, 0.0);
        }

        // Rebuild in place, reusing both the tree's own buffers and the
        // workspace.
        bvh.rebuild(&mut workspace, BvhBuildStrategy::Binned);

        // The leaf set must not grow tick over tick (updates must not
        // duplicate leaves).
        assert_eq!(bvh.leaf_count() as usize, NUM_LEAVES);

        // The queries must reflect the new leaf positions.
        for i in [0, 1, NUM_LEAVES / 2, NUM_LEAVES - 1] {
            let query = leaf_aabb(i, tick);
            let hits: Vec<u32> = bvh.intersect_aabb(&query).collect();
            assert!(
                hits.contains(&(i as u32)),
                "tick {tick}: leaf {i} not found at its updated position"
            );
            // And no leaf may be reported at its previous (tick - 1) position
            // unless it actually overlaps the new one.
            for &hit in &hits {
                let hit_aabb = leaf_aabb(hit as usize, tick);
                assert!(
                    hit_aabb.intersects(&query),
                    "tick {tick}: stale leaf {hit} reported"
                );
            }
        }
    }
}
