use crate::math::{Point, Real};
use std::iter;

/// Given an index buffer, remove from `points` every point that is not indexed.
pub fn remove_unused_points(points: &mut Vec<Point<Real>>, idx: &mut [[u32; 3]]) {
    let mut used: Vec<bool> = iter::repeat(false).take(points.len()).collect();
    let mut remap: Vec<usize> = (0..points.len()).map(|i| i).collect();
    let used = &mut used[..];
    let remap = &mut remap[..];

    for i in idx.iter() {
        used[i[0] as usize] = true;
        used[i[1] as usize] = true;
        used[i[2] as usize] = true;
    }

    let mut i = 0;
    while i != points.len() {
        if !used[i] {
            let _ = points.swap_remove(i);
            remap[points.len()] = i;
            used[i] = used[points.len()];
        } else {
            i = i + 1;
        }
    }

    for id in idx.iter_mut() {
        id[0] = remap[id[0] as usize] as u32;
        id[1] = remap[id[1] as usize] as u32;
        id[2] = remap[id[2] as usize] as u32;
    }
}
