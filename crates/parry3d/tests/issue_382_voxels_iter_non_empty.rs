// Regression test for https://github.com/dimforge/parry/issues/382
//
// The sparse chunk storage makes `Voxels::voxels()` and `Voxels::voxels_in_range()`
// yield only non-empty voxels, which the doc-comments used to contradict.

use parry3d::math::{IVector, Vector};
use parry3d::shape::Voxels;

#[test]
fn voxels_iterators_only_yield_non_empty_voxels() {
    // An L-shaped set of voxels: the domain bounding box contains empty cells
    // (e.g. (1, 1, 0), (2, 1, 0), ...).
    let coords = [
        IVector::new(0, 0, 0),
        IVector::new(1, 0, 0),
        IVector::new(2, 0, 0),
        IVector::new(0, 1, 0),
        IVector::new(0, 2, 0),
    ];
    let voxels = Voxels::new(Vector::new(1.0, 1.0, 1.0), &coords);

    // `voxels()` yields exactly the filled voxels, all non-empty.
    let all: Vec<_> = voxels.voxels().collect();
    assert_eq!(all.len(), coords.len());
    assert!(all.iter().all(|v| !v.state.is_empty()));
    for c in &coords {
        assert!(all.iter().any(|v| v.grid_coords == *c));
    }

    // `voxels_in_range()` over a range strictly larger than the domain still
    // only yields the non-empty voxels.
    let in_range: Vec<_> = voxels
        .voxels_in_range(IVector::new(-10, -10, -10), IVector::new(10, 10, 10))
        .collect();
    assert_eq!(in_range.len(), coords.len());
    assert!(in_range.iter().all(|v| !v.state.is_empty()));

    // A range covering only empty cells of the domain's bounding box yields
    // nothing.
    let empty_range: Vec<_> = voxels
        .voxels_in_range(IVector::new(1, 1, 0), IVector::new(3, 3, 1))
        .collect();
    assert!(empty_range.is_empty());
}
