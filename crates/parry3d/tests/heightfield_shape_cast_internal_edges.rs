//! A shape cast onto a flat face of a heightfield with `HeightFieldFlags::FIX_INTERNAL_EDGES`
//! returns the face's normal, not the normal of an internal edge between two of its coplanar
//! triangles.
//!
//! The cast tests each triangle on its own, and its time of impact can come out early by up to
//! about a millimetre, so the contact with the edge of a neighbouring triangle can come out first,
//! tilted off the face by a few degrees.

use parry3d::math::{Pose, Real, Rotation, Vector};
use parry3d::query::{cast_shapes, ShapeCastHit, ShapeCastOptions};
use parry3d::shape::{
    Capsule, Cuboid, HeightField, HeightFieldCellStatus, HeightFieldFlags, Shape,
};
use parry3d::utils::Array2;

/// How far a normal may lie from the face's: the noise of a cast onto a face, well below the
/// normal of an internal edge, which lies a few degrees off.
const TOLERANCE: Real = 1.0e-3;

/// A heightfield 40 m wide along `x` on 1 m columns and `length` long along `z` on 40 rows, centred
/// on the origin, with `height` at `x`.
fn heightfield(
    height: impl Fn(Real) -> Real,
    flags: HeightFieldFlags,
    length: Real,
) -> HeightField {
    let samples = 41;
    // Rows run along `z` and columns along `x`, stored column by column.
    let heights = (0..samples)
        .flat_map(|j| (0..samples).map(move |_| -20.0 + j as Real))
        .map(height)
        .collect();
    HeightField::with_flags(
        Array2::new(samples, samples, heights),
        Vector::new(40.0, 1.0, length),
        flags,
    )
}

/// The normal of a face of slope `degrees` rising along `+x`, and its gradient.
fn face(degrees: Real) -> (Vector, Real) {
    let angle = degrees.to_radians();
    (Vector::new(-angle.sin(), angle.cos(), 0.0), angle.tan())
}

/// The upright capsule of a character: radius 0.4 m, height 1.8 m.
fn capsule() -> Capsule {
    Capsule::new_y(0.5, 0.4)
}

/// How far down the capsule is cast.
const MAX_TIME_OF_IMPACT: Real = 0.32;

/// The start points of the casts: a grid of 40 × 100 points 1 cm apart starting at `x0`, near and
/// off the heightfield's grid lines and the diagonals of its cells, at `start_height(x)`.
fn starts(x0: Real, start_height: impl Fn(Real) -> Real) -> impl Iterator<Item = Vector> {
    (0..40).flat_map(move |ix| {
        let x = x0 + 0.0013 + ix as Real * 0.01;
        let y = start_height(x);
        (0..100).map(move |iz| Vector::new(x, y, -1.0 + 0.00071 + iz as Real * 0.01))
    })
}

/// The hits of the capsule cast straight down onto `surface` at `surface_pose` from the
/// [`starts`]. Each hit comes with its start point, and its normals in world space.
fn casts(
    surface: &dyn Shape,
    surface_pose: &Pose,
    x0: Real,
    start_height: impl Fn(Real) -> Real,
) -> Vec<(Vector, ShapeCastHit)> {
    let options = ShapeCastOptions::with_max_time_of_impact(MAX_TIME_OF_IMPACT);
    starts(x0, start_height)
        .map(|origin| {
            let hit = cast_shapes(
                surface_pose,
                Vector::ZERO,
                surface,
                &Pose::from_translation(origin),
                -Vector::Y,
                &capsule(),
                options,
            )
            .expect("the cast is supported")
            .unwrap_or_else(|| panic!("the cast from {origin:?} hits the surface"));
            (origin, hit.transform1_by(surface_pose))
        })
        .collect()
}

/// The height at `x` of the capsule's centre above the face of slope `degrees` rising along `+x`
/// through the origin: 0.2 m above where its bottom touches the face.
fn above_face(degrees: Real) -> impl Fn(Real) -> Real {
    let (normal, gradient) = face(degrees);
    let capsule = capsule();
    let touching = capsule.radius / normal.y + capsule.half_height();
    move |x| x * gradient + touching + 0.2
}

/// The hits of the capsule cast down onto the face of slope `degrees` rising along `+x` through
/// the origin, `surface` at `surface_pose`.
fn casts_onto_face(
    surface: &dyn Shape,
    surface_pose: &Pose,
    degrees: Real,
) -> Vec<(Vector, ShapeCastHit)> {
    casts(surface, surface_pose, 1.6, above_face(degrees))
}

/// The hits whose normals lie more than [`TOLERANCE`] off `face`.
fn off_face(hits: &[(Vector, ShapeCastHit)], face: Vector) -> Vec<(Vector, Vector, Vector)> {
    hits.iter()
        .filter(|(_, hit)| {
            (hit.normal1 - face).length() > TOLERANCE || (hit.normal2 + face).length() > TOLERANCE
        })
        .map(|(origin, hit)| (*origin, hit.normal1, hit.normal2))
        .collect()
}

/// The heights of a ridge: flat up to `x = 1`, a grid line, then falling along `+x` at `degrees`.
fn ridge(degrees: Real) -> impl Fn(Real) -> Real + Copy {
    let (_, gradient) = face(degrees);
    move |x| -(x - 1.0).max(0.0) * gradient
}

/// The normal of the ridge of slope `degrees` at its point nearest `centre`, a point above it.
fn ridge_normal(degrees: Real, centre: Vector) -> Vector {
    let angle = degrees.to_radians();
    let down_slope = Vector::new(angle.cos(), -angle.sin(), 0.0);
    let on_flat = Vector::new(centre.x.min(1.0), 0.0, centre.z);
    let on_ridge = Vector::new(1.0, 0.0, centre.z);
    let on_slope = on_ridge + down_slope * (centre - on_ridge).dot(down_slope).max(0.0);
    let nearest = if centre.distance(on_flat) <= centre.distance(on_slope) {
        on_flat
    } else {
        on_slope
    };
    (centre - nearest).normalize()
}

#[test]
fn cast_onto_a_tilted_heightfield_face_returns_the_face_normal() {
    for degrees in [45.16, 46.0, 48.0] {
        let (normal, gradient) = face(degrees);
        let hf = heightfield(|x| x * gradient, HeightFieldFlags::FIX_INTERNAL_EDGES, 40.0);

        let off = off_face(&casts_onto_face(&hf, &Pose::identity(), degrees), normal);
        assert!(
            off.is_empty(),
            "{} casts onto the {degrees}° face return a normal off the face {normal:?}, e.g. {:?}",
            off.len(),
            &off[..off.len().min(3)],
        );
    }
}

#[test]
fn cast_onto_a_cuboid_in_the_face_plane_returns_the_face_normal() {
    for degrees in [45.16, 46.0, 48.0] {
        let (normal, _) = face(degrees);
        // A slab whose top face lies in the face's plane through the origin.
        let half_thickness = 0.5;
        let slab = Cuboid::new(Vector::new(20.0, half_thickness, 20.0));
        let slab_pose = Pose::from_parts(
            -normal * half_thickness,
            Rotation::from_rotation_z(degrees.to_radians()),
        );

        let off = off_face(&casts_onto_face(&slab, &slab_pose, degrees), normal);
        assert!(
            off.is_empty(),
            "{} casts onto the {degrees}° cuboid return a normal off the face {normal:?}, e.g. {:?}",
            off.len(),
            &off[..off.len().min(3)],
        );
    }
}

#[test]
fn the_flag_corrects_only_the_normals_of_a_cast() {
    let degrees = 46.0;
    let (normal, gradient) = face(degrees);
    let fixed = heightfield(|x| x * gradient, HeightFieldFlags::FIX_INTERNAL_EDGES, 40.0);
    let plain = heightfield(|x| x * gradient, HeightFieldFlags::empty(), 40.0);
    let fixed_hits = casts_onto_face(&fixed, &Pose::identity(), degrees);
    let plain_hits = casts_onto_face(&plain, &Pose::identity(), degrees);

    for ((origin, fixed), (_, plain)) in fixed_hits.iter().zip(&plain_hits) {
        assert_eq!(
            (fixed.time_of_impact, fixed.witness1, fixed.witness2),
            (plain.time_of_impact, plain.witness1, plain.witness2),
            "the cast from {origin:?} stops elsewhere with the flag",
        );
    }
    // Without the flag, the normals stay those of the triangles hit, internal edges included.
    assert!(
        !off_face(&plain_hits, normal).is_empty(),
        "no cast stops on an internal edge, so the test no longer shows what the flag corrects",
    );
}

#[test]
fn cast_onto_a_ridge_keeps_the_normal_of_its_edge() {
    // Flat up to the ridge at `x = 1`, a grid line, then falling along `+x`. A capsule whose bottom
    // sphere's centre lies past the ridge by less than `radius * sin(slope)` touches the ridge
    // itself, with a normal between those of its two faces that the flag must keep. At 60°, one
    // of the casts touches the ridge at a vertex of a slope triangle that meets the flat part
    // only there.
    for degrees in [46.0, 60.0] {
        let (rising, _) = face(degrees);
        let falling = Vector::new(-rising.x, rising.y, 0.0);
        let fixed = heightfield(ridge(degrees), HeightFieldFlags::FIX_INTERNAL_EDGES, 40.0);
        let plain = heightfield(ridge(degrees), HeightFieldFlags::empty(), 40.0);
        let start_height = |_| capsule().radius + capsule().half_height() + 0.1;
        let fixed_hits = casts(&fixed, &Pose::identity(), 0.8, start_height);
        let plain_hits = casts(&plain, &Pose::identity(), 0.8, start_height);

        let mut on_ridge = 0;
        for ((origin, fixed), (_, plain)) in fixed_hits.iter().zip(&plain_hits) {
            assert!(
                (fixed.normal1 - plain.normal1).length() <= TOLERANCE,
                "the cast from {origin:?} bends the {degrees}° ridge's normal {:?} to {:?}",
                plain.normal1,
                fixed.normal1,
            );
            if (plain.normal1 - Vector::Y).length() > TOLERANCE
                && (plain.normal1 - falling).length() > TOLERANCE
            {
                on_ridge += 1;
            }
        }
        assert!(
            on_ridge > 0,
            "no cast stops on the {degrees}° ridge, so the test no longer shows that its normal \
             is kept",
        );
    }
}

#[test]
fn cast_onto_a_ridge_of_long_cells_keeps_the_normal_of_its_edge() {
    // Cells 1 km long along the ridge make slivers of the triangles beside it, whose edges run
    // nearly parallel to it. The casts there are noisier, so each normal is held against the
    // ridge's own rather than the normal without the flag: the flag may bring it closer, never
    // move it away.
    let degrees = 46.0;
    let fixed = heightfield(
        ridge(degrees),
        HeightFieldFlags::FIX_INTERNAL_EDGES,
        40_000.0,
    );
    let plain = heightfield(ridge(degrees), HeightFieldFlags::empty(), 40_000.0);
    let start_height = |_| capsule().radius + capsule().half_height() + 0.1;
    let fixed_hits = casts(&fixed, &Pose::identity(), 0.8, start_height);
    let plain_hits = casts(&plain, &Pose::identity(), 0.8, start_height);

    // The angle between a hit's normal and the ridge's normal at the capsule's bottom sphere.
    let error = |origin: Vector, hit: &ShapeCastHit| {
        let centre = origin - Vector::Y * (hit.time_of_impact + capsule().half_height());
        hit.normal1.angle_between(ridge_normal(degrees, centre))
    };
    for ((origin, fixed), (_, plain)) in fixed_hits.iter().zip(&plain_hits) {
        assert!(
            error(*origin, fixed) <= error(*origin, plain) + TOLERANCE,
            "the cast from {origin:?} moves the ridge's normal {:?} away to {:?}",
            plain.normal1,
            fixed.normal1,
        );
    }
}

#[test]
fn the_flag_corrects_the_normals_in_the_frame_of_each_shape() {
    // The heightfield is turned about `y` and comes second in the cast, so each normal is in a
    // frame of its own: the capsule's, the world's, and the heightfield's turned one.
    let degrees = 46.0;
    let (normal, gradient) = face(degrees);
    let hf = heightfield(|x| x * gradient, HeightFieldFlags::FIX_INTERNAL_EDGES, 40.0);
    let hf_pose = Pose::from_rotation(Rotation::from_rotation_y(0.5));
    let options = ShapeCastOptions::with_max_time_of_impact(MAX_TIME_OF_IMPACT);

    let mut off = Vec::new();
    for start in starts(1.6, above_face(degrees)) {
        let origin = hf_pose * start;
        let hit = cast_shapes(
            &Pose::from_translation(origin),
            -Vector::Y,
            &capsule(),
            &hf_pose,
            Vector::ZERO,
            &hf,
            options,
        )
        .expect("the cast is supported")
        .unwrap_or_else(|| panic!("the cast from {origin:?} hits the heightfield"));
        if (hit.normal1 + hf_pose.rotation * normal).length() > TOLERANCE
            || (hit.normal2 - normal).length() > TOLERANCE
        {
            off.push((origin, hit.normal1, hit.normal2));
        }
    }
    assert!(
        off.is_empty(),
        "{} casts onto the turned face return a normal off the face, e.g. {:?}",
        off.len(),
        &off[..off.len().min(3)],
    );
}

/// The face of slope `degrees` rising along `+x` with a hole: the right triangle of the cell
/// between `x = 2` and `x = 3`, `z = -1` and `z = 0` removed, which meets the grid line `x = 2`
/// only at `z = 0`.
fn face_with_hole(degrees: Real, flags: HeightFieldFlags) -> HeightField {
    let (_, gradient) = face(degrees);
    let mut hf = heightfield(|x| x * gradient, flags, 40.0);
    hf.set_cell_status(19, 22, HeightFieldCellStatus::RIGHT_TRIANGLE_REMOVED);
    hf
}

/// How far the point at `x`, `z` of the face's plane lies outside the hole of [`face_with_hole`],
/// whose points satisfy `x + z >= 2` within its cell.
fn outside_hole(x: Real, z: Real) -> Real {
    (2.0 - x - z) / Real::sqrt(2.0)
}

#[test]
fn cast_beside_a_removed_triangle_returns_the_face_normal() {
    // The casts touch the face 5 cm or more from the hole, some of them on the grid line
    // `x = 2`, an internal edge between two intact triangles, one of them in the cell with the
    // removed triangle.
    let degrees = 46.0;
    let (normal, _) = face(degrees);
    let hf = face_with_hole(degrees, HeightFieldFlags::FIX_INTERNAL_EDGES);
    // How far along `x` the capsule's bottom touches the face from its axis.
    let reach = capsule().radius * -normal.x;
    let options = ShapeCastOptions::with_max_time_of_impact(MAX_TIME_OF_IMPACT);

    let mut hits = Vec::new();
    for origin in starts(1.6, above_face(degrees)) {
        if outside_hole(origin.x + reach, origin.z) < 0.05 {
            continue;
        }
        let hit = cast_shapes(
            &Pose::identity(),
            Vector::ZERO,
            &hf,
            &Pose::from_translation(origin),
            -Vector::Y,
            &capsule(),
            options,
        )
        .expect("the cast is supported")
        .unwrap_or_else(|| panic!("the cast from {origin:?} hits the face"));
        hits.push((origin, hit));
    }

    let off = off_face(&hits, normal);
    assert!(
        off.is_empty(),
        "{} casts beside the hole return a normal off the face {normal:?}, e.g. {:?}",
        off.len(),
        &off[..off.len().min(3)],
    );
}

#[test]
fn cast_into_a_hole_keeps_the_normal_of_its_border() {
    // The casts whose bottom would touch the face's plane inside the hole drop into it until they
    // touch its border, with a normal off the face that the flag must keep.
    let degrees = 46.0;
    let (normal, _) = face(degrees);
    let fixed = face_with_hole(degrees, HeightFieldFlags::FIX_INTERNAL_EDGES);
    let plain = face_with_hole(degrees, HeightFieldFlags::empty());
    let reach = capsule().radius * -normal.x;
    // Deep enough for the capsule to reach the border of the hole.
    let options = ShapeCastOptions::with_max_time_of_impact(1.0);

    let mut on_border = 0;
    for origin in starts(1.6, above_face(degrees)) {
        if outside_hole(origin.x + reach, origin.z) > 0.0 {
            continue;
        }
        let cast = |hf: &HeightField| {
            cast_shapes(
                &Pose::identity(),
                Vector::ZERO,
                hf,
                &Pose::from_translation(origin),
                -Vector::Y,
                &capsule(),
                options,
            )
            .expect("the cast is supported")
            .unwrap_or_else(|| panic!("the cast from {origin:?} hits the border of the hole"))
        };
        let (fixed, plain) = (cast(&fixed), cast(&plain));
        assert!(
            (fixed.normal1 - plain.normal1).length() <= TOLERANCE,
            "the cast from {origin:?} bends the normal {:?} of the hole's border to {:?}",
            plain.normal1,
            fixed.normal1,
        );
        if (plain.normal1 - normal).length() > TOLERANCE {
            on_border += 1;
        }
    }
    assert!(
        on_border > 0,
        "no cast stops on the border of the hole, so the test no longer shows that its normal is \
         kept",
    );
}

#[test]
fn cast_into_a_hole_with_a_far_removed_corner_keeps_the_normal_of_its_border() {
    // The free corner of the removed triangle, at `x = 3`, `z = 0`, lies 1e20 m up, in removed
    // triangles only. Its distance to a point on the hole's border then rounds to NaN, which must
    // still leave the border's normal as it is.
    let degrees = 46.0;
    let (normal, gradient) = face(degrees);
    let spiked = |flags| {
        let samples = 41;
        let mut heights: Vec<Real> = (0..samples)
            .flat_map(|j| (0..samples).map(move |_| (-20.0 + j as Real) * gradient))
            .collect();
        // Column 23 (`x = 3`), row 20 (`z = 0`).
        heights[23 * samples + 20] = 1.0e20;
        let mut hf = HeightField::with_flags(
            Array2::new(samples, samples, heights),
            Vector::new(40.0, 1.0, 40.0),
            flags,
        );
        hf.set_cell_status(19, 22, HeightFieldCellStatus::RIGHT_TRIANGLE_REMOVED);
        for (i, j) in [(19, 23), (20, 22), (20, 23)] {
            hf.set_cell_status(i, j, HeightFieldCellStatus::CELL_REMOVED);
        }
        hf
    };
    let fixed = spiked(HeightFieldFlags::FIX_INTERNAL_EDGES);
    let plain = spiked(HeightFieldFlags::empty());
    let reach = capsule().radius * -normal.x;
    let options = ShapeCastOptions::with_max_time_of_impact(1.0);

    let mut on_border = 0;
    for origin in starts(1.6, above_face(degrees)) {
        if outside_hole(origin.x + reach, origin.z) > 0.0 {
            continue;
        }
        let cast = |hf: &HeightField| {
            cast_shapes(
                &Pose::identity(),
                Vector::ZERO,
                hf,
                &Pose::from_translation(origin),
                -Vector::Y,
                &capsule(),
                options,
            )
            .expect("the cast is supported")
        };
        let (Some(fixed), Some(plain)) = (cast(&fixed), cast(&plain)) else {
            continue;
        };
        assert!(
            (fixed.normal1 - plain.normal1).length() <= TOLERANCE,
            "the cast from {origin:?} bends the normal {:?} of the hole's border to {:?}",
            plain.normal1,
            fixed.normal1,
        );
        if (plain.normal1 - normal).length() > TOLERANCE {
            on_border += 1;
        }
    }
    assert!(
        on_border > 0,
        "no cast stops on the border of the hole, so the test no longer shows that its normal is \
         kept",
    );
}
