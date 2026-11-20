// Contacts between a small shape and a face of a very large cuboid, with the small shape
// placed from 40 micrometers off the face to 40 micrometers into it.
//
// Near touching, GJK can stall (the distance stops shrinking or a support point repeats) with
// a lower bound on the distance close to zero. The contact normal must then still be the
// face's normal: a sideways normal makes a resting or landing body stop dead or jump.

use parry3d::math::{Pose, Real, Rotation, Vector};
use parry3d::query::details::contact_support_map_support_map_with_params;
use parry3d::query::gjk::{GJKResult, VoronoiSimplex};
use parry3d::shape::{ConvexPolyhedron, Cuboid, SupportMap};

/// The largest accepted angle between the contact normal and the face's normal.
const MAX_ANGLE_DEG: Real = 1.0;

/// A wedge-shaped convex hull about 2 m wide, 1.15 m high and 4.8 m long.
fn hull_points() -> Vec<Vector> {
    let mut points = Vec::new();
    for sx in [-1.0, 1.0] {
        for sz in [-1.0, 1.0] {
            points.push(Vector::new(0.95 * sx, -0.5, 2.3 * sz));
            points.push(Vector::new(0.98 * sx, 0.0, 2.4 * sz));
            points.push(Vector::new(
                0.7 * sx,
                0.6,
                if sz > 0.0 { 1.2 } else { -1.0 },
            ));
            points.push(Vector::new(0.9 * sx, -0.55, 1.6 * sz));
        }
    }
    points
}

/// The corners of the small cuboid.
fn cuboid_points(half: Vector) -> Vec<Vector> {
    let mut points = Vec::new();
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            for sz in [-1.0, 1.0] {
                points.push(Vector::new(half.x * sx, half.y * sy, half.z * sz));
            }
        }
    }
    points
}

/// Places `small` (whose vertices are `points`) against the top face of the cuboid
/// (`half`, centered at `center`) at the point `at` of that face, at 45 attitudes (tilted
/// 8 to 16 degrees, with some roll and yaw), from 40 micrometers above the face to 40 below
/// it, and checks every contact, starting GJK from `init_dir`.
fn check_top_face_contacts<S: SupportMap>(
    small: &S,
    points: &[Vector],
    half: Vector,
    center: Vector,
    at: Vector,
    init_dir: Option<Vector>,
) {
    let normal = Vector::Y;
    let big = Cuboid::new(half);
    let big_pose = Pose::from_translation(center);
    let mut asked = 0;
    let mut lost = Vec::new();
    let mut off = Vec::new();

    let pitches: [Real; 5] = [-8.0, -10.0, -12.0, -14.0, -16.0];
    let rolls: [Real; 3] = [-3.0, 0.0, 3.0];
    let yaws: [Real; 3] = [-4.0, 0.0, 4.0];
    for pitch in pitches {
        for roll in rolls {
            for yaw in yaws {
                let rotation = Rotation::from_rotation_y(yaw.to_radians())
                    * Rotation::from_rotation_x(pitch.to_radians())
                    * Rotation::from_rotation_z(roll.to_radians());
                // The lowest vertex along the face's normal.
                let lowest = points
                    .iter()
                    .map(|p| (rotation * *p).dot(normal))
                    .fold(Real::MAX, Real::min);
                for step in -40..=40 {
                    let gap = step as Real * 1.0e-6;
                    let pose = Pose::from_parts(at - normal * (lowest - gap), rotation);
                    let pos12 = big_pose.inv_mul(&pose);
                    asked += 1;
                    let result = contact_support_map_support_map_with_params(
                        &pos12,
                        &big,
                        small,
                        0.01,
                        &mut VoronoiSimplex::new(),
                        init_dir,
                    );
                    let GJKResult::ClosestPoints(_, _, n) = result else {
                        lost.push((pitch, roll, yaw, gap));
                        continue;
                    };
                    let angle = n.dot(normal).clamp(-1.0, 1.0).acos().to_degrees();
                    if angle > MAX_ANGLE_DEG {
                        off.push((pitch, roll, yaw, gap, angle));
                    }
                }
            }
        }
    }

    assert!(
        lost.is_empty(),
        "{} of {asked} contacts lost, e.g. (pitch, roll, yaw, gap) {:?}",
        lost.len(),
        &lost[..lost.len().min(5)]
    );
    let worst = off.iter().map(|o| o.4).fold(0.0, Real::max);
    assert!(
        off.is_empty(),
        "{} of {asked} normals more than {MAX_ANGLE_DEG} degree off the face's normal \
         (worst {worst:.1} degrees), e.g. (pitch, roll, yaw, gap, angle) {:?}",
        off.len(),
        &off[..off.len().min(5)]
    );
}

/// A 300 m by 760 m slab, 1 m thick, and a point of its top face 340 m from its center.
const SLAB_HALF: Vector = Vector::new(150.0, 0.5, 380.0);
const SLAB_CENTER: Vector = Vector::new(0.0, -0.5, 0.0);
const SLAB_AT: Vector = Vector::new(-30.0, 0.0, -340.0);

#[test]
fn hull_touching_a_large_slab_gets_the_face_normal() {
    let points = hull_points();
    let hull = ConvexPolyhedron::from_convex_hull(&points).unwrap();
    for init_dir in [Some(Vector::Y), None] {
        check_top_face_contacts(&hull, &points, SLAB_HALF, SLAB_CENTER, SLAB_AT, init_dir);
    }
}

#[test]
fn cuboid_touching_a_large_slab_gets_the_face_normal() {
    let half = Vector::new(0.9, 0.6, 2.2);
    let cuboid = Cuboid::new(half);
    let points = cuboid_points(half);
    for init_dir in [Some(Vector::Y), None] {
        check_top_face_contacts(&cuboid, &points, SLAB_HALF, SLAB_CENTER, SLAB_AT, init_dir);
    }
}

/// A 2 km by 2 km slab, 2 m thick, and three points of its top face: its center, a point
/// 500 m from it, and a point 20 m in from a corner.
const BIG_SLAB_HALF: Vector = Vector::new(1000.0, 1.0, 1000.0);
const BIG_SLAB_CENTER: Vector = Vector::new(0.0, -1.0, 0.0);
const BIG_SLAB_AT: [Vector; 3] = [
    Vector::ZERO,
    Vector::new(300.0, 0.0, -400.0),
    Vector::new(980.0, 0.0, 980.0),
];

// On the 2 km slab the support points' coordinates are rounded to about 1e-4 m, more than
// the gap, so when the search stalls its last direction is noise. The direction that gave
// the best lower bound on the distance is the face's normal.
#[test]
fn hull_touching_a_2km_slab_gets_the_face_normal() {
    let points = hull_points();
    let hull = ConvexPolyhedron::from_convex_hull(&points).unwrap();
    for at in BIG_SLAB_AT {
        for init_dir in [Some(Vector::Y), None] {
            check_top_face_contacts(&hull, &points, BIG_SLAB_HALF, BIG_SLAB_CENTER, at, init_dir);
        }
    }
}

#[test]
fn cuboid_touching_a_2km_slab_gets_the_face_normal() {
    let half = Vector::new(0.9, 0.6, 2.2);
    let cuboid = Cuboid::new(half);
    let points = cuboid_points(half);
    for at in BIG_SLAB_AT {
        for init_dir in [Some(Vector::Y), None] {
            check_top_face_contacts(
                &cuboid,
                &points,
                BIG_SLAB_HALF,
                BIG_SLAB_CENTER,
                at,
                init_dir,
            );
        }
    }
}
