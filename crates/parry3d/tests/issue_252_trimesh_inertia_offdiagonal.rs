// Regression test for https://github.com/dimforge/parry/issues/252
//
// The issue reported swapped `xy`/`xz` off-diagonal terms in a trimesh's inertia tensor.
// This checks the full tensor against the analytic cuboid one, both axis-aligned and
// with a rotation baked into the vertices.

use parry3d::mass_properties::MassProperties;
use parry3d::math::{Mat3, Rot3, Vector};

fn box_mesh(he: Vector) -> (Vec<Vector>, Vec<[u32; 3]>) {
    // Outward-oriented box triangulation.
    let vertices = vec![
        Vector::new(-he.x, -he.y, -he.z),
        Vector::new(he.x, -he.y, -he.z),
        Vector::new(he.x, he.y, -he.z),
        Vector::new(-he.x, he.y, -he.z),
        Vector::new(-he.x, -he.y, he.z),
        Vector::new(he.x, -he.y, he.z),
        Vector::new(he.x, he.y, he.z),
        Vector::new(-he.x, he.y, he.z),
    ];
    let indices = vec![
        [0, 2, 1],
        [0, 3, 2], // -z
        [4, 5, 6],
        [4, 6, 7], // +z
        [0, 1, 5],
        [0, 5, 4], // -y
        [2, 3, 7],
        [2, 7, 6], // +y
        [1, 2, 6],
        [1, 6, 5], // +x
        [0, 4, 7],
        [0, 7, 3], // -x
    ];
    (vertices, indices)
}

fn assert_mat_relative_eq(actual: Mat3, expected: Mat3, scale: f32) {
    let a = actual.to_cols_array();
    let e = expected.to_cols_array();
    for i in 0..9 {
        assert!(
            (a[i] - e[i]).abs() <= scale * 1.0e-3,
            "inertia mismatch at {i}: {} vs {} (actual {a:?}, expected {e:?})",
            a[i],
            e[i]
        );
    }
}

#[test]
fn aligned_box_mesh_inertia_matches_cuboid() {
    let he = Vector::new(0.4, 0.6, 0.8);
    let density = 2.5;

    let (vertices, indices) = box_mesh(he);
    let mp_mesh = MassProperties::from_trimesh(density, &vertices, &indices);
    let mp_cuboid = MassProperties::from_cuboid(density, he);

    let mass_mesh = 1.0 / mp_mesh.inv_mass;
    let mass_cuboid = 1.0 / mp_cuboid.inv_mass;
    assert!((mass_mesh - mass_cuboid).abs() <= mass_cuboid * 1.0e-4);
    assert!(mp_mesh.local_com.length() <= 1.0e-5);

    let i_mesh = mp_mesh.reconstruct_inertia_matrix();
    let i_cuboid = mp_cuboid.reconstruct_inertia_matrix();

    // Scale for the relative tolerance: the largest diagonal term.
    let scale = i_cuboid.to_cols_array().iter().fold(0.0f32, |m, x| m.max(x.abs()));
    assert_mat_relative_eq(i_mesh, i_cuboid, scale);

    // The aligned case must have (near-)zero off-diagonals.
    let m = i_mesh.to_cols_array_2d();
    for (i, col) in m.iter().enumerate() {
        for (j, v) in col.iter().enumerate() {
            if i != j {
                assert!(v.abs() <= scale * 1.0e-4, "off-diagonal [{i}][{j}] = {v}");
            }
        }
    }
}

#[test]
fn rotated_box_mesh_inertia_matches_rotated_cuboid_tensor() {
    let he = Vector::new(0.4, 0.6, 0.8);
    let density = 2.5;
    let rot = Rot3::from_axis_angle(Vector::new(1.0, 2.0, 3.0).normalize(), 0.7);

    let (vertices, indices) = box_mesh(he);
    let rotated: Vec<_> = vertices.iter().map(|v| rot * *v).collect();
    let mp_mesh = MassProperties::from_trimesh(density, &rotated, &indices);
    let mp_cuboid = MassProperties::from_cuboid(density, he);

    let i_mesh = mp_mesh.reconstruct_inertia_matrix();

    // Expected: R * I * R^T. This has non-zero off-diagonal terms, so an
    // xy/xz swap in `reconstruct_inertia_matrix` (the original report) or in
    // the trimesh accumulation would be caught here.
    let r = Mat3::from_quat(rot);
    let i_expected = r * mp_cuboid.reconstruct_inertia_matrix() * r.transpose();

    let scale = i_expected.to_cols_array().iter().fold(0.0f32, |m, x| m.max(x.abs()));
    assert_mat_relative_eq(i_mesh, i_expected, scale);

    // Sanity: the rotated tensor really has significant off-diagonal terms.
    let e = i_expected.to_cols_array_2d();
    assert!(e[0][1].abs() > 1.0e-3 || e[0][2].abs() > 1.0e-3 || e[1][2].abs() > 1.0e-3);
}
