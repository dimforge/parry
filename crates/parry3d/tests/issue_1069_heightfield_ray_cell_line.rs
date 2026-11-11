// Regression test for https://github.com/dimforge/rapier/issues/1069
//
// A ray parallel to the y axis that lands exactly on a cell line of a heightfield used to miss it
// for some lines when the cell count was not a power of two: `cell_at_point` could round to the
// cell on the other side of the line from the triangles containing the point, and a vertical ray
// never steps into the neighboring cell.

use parry3d::math::{Pose, Vector};
use parry3d::query::{Ray, RayCast};
use parry3d::shape::HeightField;
use parry3d::utils::Array2;

fn heightfield(rows: usize, cols: usize, cell_size: f32) -> HeightField {
    let h = |x: f32, z: f32| 0.03 * x + 0.4 * (x / 9.0).sin() + 0.3 * (z / 7.0).cos();
    let mut heights = Vec::with_capacity((rows + 1) * (cols + 1));
    for j in 0..=cols {
        for i in 0..=rows {
            heights.push(h(j as f32 * cell_size, i as f32 * cell_size));
        }
    }
    HeightField::new(
        Array2::new(rows + 1, cols + 1, heights),
        Vector::new(cols as f32 * cell_size, 1.0, rows as f32 * cell_size),
    )
}

/// Casts straight down at every quarter-cell point (cell lines included) and returns the
/// points that missed.
fn vertical_ray_misses(rows: usize, cols: usize, cell_size: f32) -> Vec<(f32, f32)> {
    let hf = heightfield(rows, cols, cell_size);
    let width = cols as f32 * cell_size;
    let depth = rows as f32 * cell_size;
    let pos = Pose::from_translation(Vector::new(width / 2.0, 0.0, depth / 2.0));
    let step = cell_size / 4.0;
    let mut misses = Vec::new();

    for ix in 0..=cols * 4 {
        for iz in 0..=rows * 4 {
            let (x, z) = (ix as f32 * step, iz as f32 * step);
            let ray = Ray::new(Vector::new(x, 1000.0, z), Vector::new(0.0, -1.0, 0.0));
            if hf
                .cast_ray_and_get_normal(&pos, &ray, 4000.0, true)
                .is_none()
            {
                misses.push((x, z));
            }
        }
    }

    misses
}

#[test]
fn vertical_ray_hits_heightfield_on_cell_lines() {
    for (rows, cols, cell_size) in [
        (64, 64, 2.0),
        (48, 48, 2.0),
        (24, 24, 2.0),
        (16, 48, 2.0),
        (48, 16, 2.0),
        (16, 17, 2.0),
        (24, 24, 1.0),
        (16, 17, 0.5),
        (48, 48, 1.0),
    ] {
        let misses = vertical_ray_misses(rows, cols, cell_size);
        assert!(
            misses.is_empty(),
            "{rows} x {cols} (cell size {cell_size}): {} rays missed, first at {:?}",
            misses.len(),
            &misses[..misses.len().min(3)]
        );
    }
}
