// A ray cast straight down onto a heightfield, or tilted along a grid line only, used to miss
// within a few rounding steps of that line: the ray's start cell comes from `cell_at_point`, the
// cell's triangles from `x_at`/`z_at`, which round differently, and a ray that does not move
// across the line never visits the cell on its other side.

use approx::{assert_relative_eq, relative_eq};
use parry3d::math::{Real, Vector};
use parry3d::query::{Ray, RayCast};
use parry3d::shape::{HeightField, HeightFieldCellStatus};
use parry3d::utils::Array2;

/// The sizes, in samples per side and length per side, of the heightfields cast at. On the last
/// two, `-0.5 + unit_cell_width * cells` rounds short of `0.5`.
const HEIGHTFIELDS: [(usize, Real); 6] = [
    (41, 40.0),
    (1001, 1000.0),
    (513, 1000.0),
    (101, 1000.0),
    (42, 1000.0),
    (210, 40.0),
];

/// How many rounding steps away from a grid line the rays are cast, on each side.
const ULPS: i32 = 16;

fn flat_heightfield(samples: usize, size: Real) -> HeightField {
    HeightField::new(
        Array2::zeros(samples, samples),
        Vector::new(size, 1.0, size),
    )
}

/// The value `steps` rounding steps away from `x`.
fn ulps_away(x: Real, steps: i32) -> Real {
    let mut x = x;
    for _ in 0..steps.abs() {
        x = if steps > 0 {
            x.next_up()
        } else {
            x.next_down()
        };
    }
    x
}

/// The coordinates at and around grid line `j` of `heightfield`, `size` long:
/// the line as `x_at` places it, as a caller would compute it, and as the nearest `Real` to the
/// exact value, each with its neighbours up to [`ULPS`] rounding steps away, and only those
/// within the heightfield.
fn coordinates_around_line(heightfield: &HeightField, j: usize, size: Real) -> Vec<Real> {
    let cells = heightfield.ncols();
    let half = size / 2.0;
    let exact = (j as f64 * size as f64 / cells as f64 - half as f64) as Real;
    let by_caller = j as Real * (size / cells as Real) - half;
    let mut coordinates = Vec::new();
    for base in [heightfield.x_at(j), by_caller, exact] {
        for steps in -ULPS..=ULPS {
            let x = ulps_away(base, steps);
            if (-half..=half).contains(&x) {
                coordinates.push(x);
            }
        }
    }
    coordinates
}

fn cast_down(heightfield: &HeightField, x: Real, z: Real) -> Option<Real> {
    let ray = Ray::new(Vector::new(x, 1.0, z), -Vector::Y);
    heightfield.cast_local_ray(&ray, Real::MAX, true)
}

#[test]
fn a_vertical_ray_on_a_grid_line_hits_a_flat_heightfield() {
    for (samples, size) in HEIGHTFIELDS {
        let heightfield = flat_heightfield(samples, size);
        let cell = size / heightfield.ncols() as Real;
        let mut misses = Vec::new();

        for j in 0..=heightfield.ncols() {
            // Mid-cell along the other axis, so that only the one grid line is near.
            let across = heightfield.x_at(j.min(heightfield.ncols() - 1)) + cell / 2.0;
            for along in coordinates_around_line(&heightfield, j, size) {
                for (x, z) in [(along, across), (across, along), (along, along)] {
                    match cast_down(&heightfield, x, z) {
                        Some(time_of_impact) => assert_relative_eq!(time_of_impact, 1.0),
                        None => misses.push((x, z)),
                    }
                }
            }
        }

        assert!(
            misses.is_empty(),
            "{samples} samples, {size} long: {} vertical rays missed, e.g. at {:?}",
            misses.len(),
            &misses[..misses.len().min(5)]
        );
    }
}

/// Tiny deterministic LCG so the tests don't depend on rand seeding.
struct Lcg(u64);

impl Lcg {
    fn real(&mut self, min: Real, max: Real) -> Real {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        min + (max - min) * ((self.0 >> 32) as u32 as Real / u32::MAX as Real)
    }
}

/// The nearest hit of `ray` among the triangles of `heightfield`, each cast on its own.
fn nearest_triangle_hit(heightfield: &HeightField, ray: &Ray) -> Option<Real> {
    heightfield
        .triangles()
        .filter_map(|triangle| triangle.cast_local_ray(ray, Real::MAX, true))
        .min_by(Real::total_cmp)
}

#[test]
fn a_ray_near_a_grid_line_gives_the_nearest_triangle_hit() {
    let (samples, size) = (11, 10.0);
    let mut rng = Lcg(0x5eed);
    let heights = (0..samples * samples)
        .map(|_| rng.real(-1.0, 1.0))
        .collect();
    let heights = Array2::new(samples, samples, heights);
    let mut heightfield = HeightField::new(heights, Vector::new(size, 2.0, size));
    // Both ways of splitting a cell into triangles, in a checkerboard.
    for i in 0..heightfield.nrows() {
        for j in (i % 2..heightfield.ncols()).step_by(2) {
            heightfield.set_cell_status(i, j, HeightFieldCellStatus::ZIGZAG_SUBDIVISION);
        }
    }

    // Straight down; tilted by 1e-4, so that the ray stays within half a millimetre of where it
    // started; and tilted across several cells.
    let directions = [
        Vector::new(0.0, -1.0, 0.0),
        Vector::new(1e-4, -1.0, 0.0),
        Vector::new(-1e-4, -1.0, 0.0),
        Vector::new(0.0, -1.0, 1e-4),
        Vector::new(0.0, -1.0, -1e-4),
        Vector::new(0.3, -1.0, -0.2),
        Vector::new(-0.2, -1.0, 0.3),
    ];
    let mut casts = 0;

    for j in 0..=heightfield.ncols() {
        for along in coordinates_around_line(&heightfield, j, size) {
            let across = rng.real(-size / 2.0, size / 2.0);
            for (x, z) in [(along, across), (across, along), (along, along)] {
                for dir in directions {
                    let ray = Ray::new(Vector::new(x, 3.0, z), dir.normalize());
                    let Some(expected) = nearest_triangle_hit(&heightfield, &ray) else {
                        continue;
                    };
                    casts += 1;
                    let hit = heightfield.cast_local_ray(&ray, Real::MAX, true);
                    assert!(
                        hit.is_some_and(|time_of_impact| {
                            relative_eq!(time_of_impact, expected, max_relative = 1e-5)
                        }),
                        "{ray:?}: the heightfield gave {hit:?}, the nearest triangle {expected}"
                    );
                }
            }
        }
    }

    // The comparison ran: most rays hit a triangle, the others leave through a side.
    assert!(casts > 1000, "only {casts} rays hit a triangle");
}

#[test]
fn a_vertical_ray_on_a_grid_line_hits_beside_a_removed_cell() {
    // A point on a grid line lies in the cells on both of its sides, and a point where two lines
    // cross in the four cells around it: with all but one of them removed, a ray there, straight
    // down or tilted along the line, still hits.
    let (samples, size) = (41, 40.0);
    let mut misses = Vec::new();

    for line in 1..samples - 1 {
        for removed in [line - 1, line] {
            let mut heightfield = flat_heightfield(samples, size);
            for k in 0..heightfield.ncols() {
                heightfield.set_cell_status(k, removed, HeightFieldCellStatus::CELL_REMOVED);
                heightfield.set_cell_status(removed, k, HeightFieldCellStatus::CELL_REMOVED);
            }

            let on_line = heightfield.x_at(line);
            assert_eq!(heightfield.z_at(line), on_line);
            // Mid-cell, away from the removed row and column.
            let far = if line < samples / 2 { samples - 2 } else { 0 };
            let across = (heightfield.x_at(far) + heightfield.x_at(far + 1)) / 2.0;
            let rays = [
                (on_line, across, Vector::new(0.0, -1.0, 0.0)),
                (on_line, across, Vector::new(0.0, -1.0, 1e-4)),
                (across, on_line, Vector::new(0.0, -1.0, 0.0)),
                (across, on_line, Vector::new(1e-4, -1.0, 0.0)),
                (on_line, on_line, Vector::new(0.0, -1.0, 0.0)),
            ];
            for (x, z, dir) in rays {
                let ray = Ray::new(Vector::new(x, 1.0, z), dir.normalize());
                if heightfield.cast_local_ray(&ray, Real::MAX, true).is_none() {
                    misses.push((line, removed, ray));
                }
            }
        }
    }

    assert!(
        misses.is_empty(),
        "{} rays missed (line, removed row and column, ray), e.g. {:?}",
        misses.len(),
        &misses[..misses.len().min(5)]
    );
}

#[test]
fn a_ray_along_a_grid_line_hits_beside_removed_cells_after_walking() {
    // A ray that runs along a grid line, not straight down, walks the cells on one side of it,
    // and the cells beside them stay the ones to test as it walks. On a slope rising towards the
    // far side, the ray meets the surface several rows (columns) in, with the column (row) it
    // starts in removed.
    let (samples, size) = (11, 10.0);
    let heights = (0..samples * samples)
        .map(|k| 0.1 * ((k / samples) + (k % samples)) as Real)
        .collect();
    let slope = HeightField::new(
        Array2::new(samples, samples, heights),
        Vector::new(size, 1.0, size),
    );
    let mut misses = Vec::new();
    let mut casts = 0;

    for line in 1..samples - 1 {
        let on_line = slope.x_at(line);
        assert_eq!(slope.z_at(line), on_line);
        // Along the column line from the near side in z, then along the row line in x.
        let rays = [
            Ray::new(
                Vector::new(on_line, 2.1, -6.0),
                Vector::new(0.0, -0.1, 1.0).normalize(),
            ),
            Ray::new(
                Vector::new(-6.0, 2.1, on_line),
                Vector::new(1.0, -0.1, 0.0).normalize(),
            ),
        ];
        for (along_column, ray) in [true, false].into_iter().zip(rays) {
            let mut heightfield = slope.clone();
            // The cell the walk starts in, mid-cell along the walk.
            let start = if along_column {
                Vector::new(on_line, 0.0, -4.5)
            } else {
                Vector::new(-4.5, 0.0, on_line)
            };
            let (row, column) = heightfield.cell_at_point(start).unwrap();
            if along_column {
                for i in 0..heightfield.nrows() {
                    heightfield.set_cell_status(i, column, HeightFieldCellStatus::CELL_REMOVED);
                }
            } else {
                for j in 0..heightfield.ncols() {
                    heightfield.set_cell_status(row, j, HeightFieldCellStatus::CELL_REMOVED);
                }
            }

            // A ray on a triangle's edge may miss it by rounding, as it may miss the heightfield.
            let Some(expected) = nearest_triangle_hit(&heightfield, &ray) else {
                continue;
            };
            casts += 1;
            let hit = heightfield.cast_local_ray_and_get_normal(&ray, Real::MAX, true);
            let walked = hit.is_some_and(|hit| {
                let point = ray.point_at(hit.time_of_impact);
                let (walked, first_line) = if along_column {
                    (point.z, heightfield.z_at(1))
                } else {
                    (point.x, heightfield.x_at(1))
                };
                // The hit lies past the first row (column), on a triangle the ray hits there.
                walked > first_line
                    && relative_eq!(hit.time_of_impact, expected, max_relative = 1e-5)
                    && heightfield
                        .triangle_at_id(hit.subshape)
                        .and_then(|triangle| triangle.cast_local_ray(&ray, Real::MAX, true))
                        .is_some_and(|time_of_impact| {
                            relative_eq!(time_of_impact, hit.time_of_impact, max_relative = 1e-5)
                        })
            });
            if !walked {
                misses.push((line, along_column, hit, expected));
            }
        }
    }

    assert!(
        misses.is_empty(),
        "{} rays gave no hit, or the wrong one (line, along the column, hit, nearest triangle), \
         e.g. {:?}",
        misses.len(),
        &misses[..misses.len().min(5)]
    );
    // The comparison ran: most rays hit a triangle of the cells left.
    assert!(casts > samples, "only {casts} rays hit a triangle");
}

#[test]
fn the_cell_at_a_point_of_a_mirrored_heightfield_lies_under_the_point() {
    let heightfield = HeightField::new(Array2::zeros(41, 41), Vector::new(-40.0, 1.0, -40.0));

    for j in 0..heightfield.ncols() {
        let mid = (heightfield.x_at(j) + heightfield.x_at(j + 1)) / 2.0;
        let point = Vector::new(mid, 0.0, mid);
        assert_eq!(
            heightfield.cell_at_point(point),
            Some((j, j)),
            "at {point:?}"
        );
    }
}
