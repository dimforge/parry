// A ray cast straight down onto a heightfield used to miss within a few rounding steps of a
// grid point: the ray's start cell comes from `cell_at_point` and the cell's segment from
// `segment_at`, which round differently; a segment did not end where the next one starts; and
// the line-line parameters put a ray on a sloped segment's end outside of it. A vertical ray
// never visits a neighbouring cell.

use approx::{assert_relative_eq, relative_eq};
use parry2d::math::{Real, Vector};
use parry2d::query::{Ray, RayCast};
use parry2d::shape::HeightField;

/// The sizes, in samples and length, of the heightfields cast at. On the last two,
/// `-0.5 + unit_cell_width * cells` rounds short of `0.5`.
const HEIGHTFIELDS: [(usize, Real); 6] = [
    (41, 40.0),
    (1001, 1000.0),
    (513, 1000.0),
    (101, 1000.0),
    (42, 1000.0),
    (210, 40.0),
];

/// How many rounding steps away from a grid point the rays are cast, on each side.
const ULPS: i32 = 16;

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

/// The coordinates at and around grid point `i` of a heightfield `size` long: the point as the
/// segments on either side of it end, as a caller would compute it, and as the nearest `Real` to
/// the exact value, each with its neighbours up to [`ULPS`] rounding steps away, and only those
/// within the heightfield.
fn coordinates_around_point(heightfield: &HeightField, i: usize, size: Real) -> Vec<Real> {
    let cells = heightfield.num_cells();
    let half = size / 2.0;
    let mut bases = vec![
        i as Real * (size / cells as Real) - half,
        (i as f64 * size as f64 / cells as f64 - half as f64) as Real,
    ];
    if let Some(right) = heightfield.segment_at(i) {
        bases.push(right.a.x);
    }
    if let Some(left) = i.checked_sub(1).and_then(|i| heightfield.segment_at(i)) {
        bases.push(left.b.x);
    }

    let mut coordinates = Vec::new();
    for base in bases {
        for steps in -ULPS..=ULPS {
            let x = ulps_away(base, steps);
            if (-half..=half).contains(&x) {
                coordinates.push(x);
            }
        }
    }
    coordinates
}

#[test]
fn a_vertical_ray_on_a_grid_point_hits_a_flat_heightfield() {
    for (samples, size) in HEIGHTFIELDS {
        let heightfield = HeightField::new(vec![0.0; samples], Vector::new(size, 1.0));
        let mut misses = Vec::new();

        for i in 0..=heightfield.num_cells() {
            for x in coordinates_around_point(&heightfield, i, size) {
                let ray = Ray::new(Vector::new(x, 1.0), -Vector::Y);
                match heightfield.cast_local_ray(&ray, Real::MAX, true) {
                    Some(time_of_impact) => assert_relative_eq!(time_of_impact, 1.0),
                    None => misses.push(x),
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

/// The nearest hit of `ray` among the segments of `heightfield`, each cast on its own.
fn nearest_segment_hit(heightfield: &HeightField, ray: &Ray) -> Option<Real> {
    heightfield
        .segments()
        .filter_map(|segment| segment.cast_local_ray(ray, Real::MAX, true))
        .min_by(Real::total_cmp)
}

#[test]
fn a_ray_near_a_grid_point_gives_the_nearest_segment_hit() {
    let (samples, size) = (11, 10.0);
    let mut rng = Lcg(0x5eed);
    let heights = (0..samples).map(|_| rng.real(-1.0, 1.0)).collect();
    let heightfield = HeightField::new(heights, Vector::new(size, 2.0));

    // Straight down; tilted by 1e-4, so that the ray stays within half a millimetre of where it
    // started; and tilted across several cells.
    let directions = [
        Vector::new(0.0, -1.0),
        Vector::new(1e-4, -1.0),
        Vector::new(-1e-4, -1.0),
        Vector::new(0.3, -1.0),
        Vector::new(-0.3, -1.0),
    ];
    let mut casts = 0;

    for i in 0..=heightfield.num_cells() {
        for x in coordinates_around_point(&heightfield, i, size) {
            for dir in directions {
                let ray = Ray::new(Vector::new(x, 3.0), dir.normalize());
                let Some(expected) = nearest_segment_hit(&heightfield, &ray) else {
                    continue;
                };
                casts += 1;
                let hit = heightfield.cast_local_ray(&ray, Real::MAX, true);
                assert!(
                    hit.is_some_and(|time_of_impact| {
                        relative_eq!(time_of_impact, expected, max_relative = 1e-5)
                    }),
                    "{ray:?}: the heightfield gave {hit:?}, the nearest segment {expected}"
                );
            }
        }
    }

    // The comparison ran: most rays hit a segment, the others leave through a side.
    assert!(casts > 300, "only {casts} rays hit a segment");
}

#[test]
fn a_vertical_ray_on_a_grid_point_hits_beside_a_removed_segment() {
    // A grid point is the end of the segments on both of its sides: with either removed, a ray
    // there still hits the other.
    let samples = 41;
    let mut misses = Vec::new();

    for point in 1..samples - 1 {
        for removed in [point - 1, point] {
            let mut heightfield = HeightField::new(vec![0.0; samples], Vector::new(40.0, 1.0));
            let x = heightfield.segment_at(point).unwrap().a.x;
            heightfield.set_segment_removed(removed, true);

            let ray = Ray::new(Vector::new(x, 1.0), -Vector::Y);
            if heightfield.cast_local_ray(&ray, Real::MAX, true).is_none() {
                misses.push((point, removed, x));
            }
        }
    }

    assert!(
        misses.is_empty(),
        "{} vertical rays missed (grid point, removed segment, x), e.g. {:?}",
        misses.len(),
        &misses[..misses.len().min(5)]
    );
}

#[test]
fn the_cell_at_a_point_of_a_mirrored_heightfield_lies_under_the_point() {
    let heightfield = HeightField::new(vec![0.0; 41], Vector::new(-40.0, 1.0));

    for i in 0..heightfield.num_cells() {
        let segment = heightfield.segment_at(i).unwrap();
        let point = Vector::new((segment.a.x + segment.b.x) / 2.0, 0.0);
        assert_eq!(heightfield.cell_at_point(point), Some(i), "at {point:?}");
    }
}
