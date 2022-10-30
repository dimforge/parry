use crate::bounding_volume::Aabb;
use crate::math::{Point, Real};

pub fn split_indices_wrt_dim<'a>(
    indices: &'a mut [usize],
    aabbs: &[Aabb],
    split_point: &Point<Real>,
    dim: usize,
    enable_fallback_split: bool,
) -> (&'a mut [usize], &'a mut [usize]) {
    let mut icurr = 0;
    let mut ilast = indices.len();

    // The loop condition we can just do 0..indices.len()
    // instead of the test icurr < ilast because we know
    // we will iterate exactly once per index.
    for _ in 0..indices.len() {
        let i = indices[icurr];
        let center = aabbs[i].center();

        if center[dim] > split_point[dim] {
            ilast -= 1;
            indices.swap(icurr, ilast);
        } else {
            icurr += 1;
        }
    }

    if enable_fallback_split && (icurr == 0 || icurr == indices.len()) {
        // We don't want to return one empty set. But
        // this can happen if all the coordinates along the
        // given dimension are equal.
        // In this is the case, we just split in the middle.
        let half = indices.len() / 2;
        indices.split_at_mut(half)
    } else {
        indices.split_at_mut(icurr)
    }
}
