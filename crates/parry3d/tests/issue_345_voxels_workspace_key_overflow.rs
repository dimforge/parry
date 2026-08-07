// Regression test for https://github.com/dimforge/parry/issues/345
// (originally reported as https://github.com/dimforge/rapier/issues/845)
//
// The contact-manifold workspace key was computed over the undilated voxel domain, so a
// voxel with a non-free face toward a domain minimum produced a relative key of -1,
// overflowing on the cast to u32. It is now computed over a domain dilated by 1.

use parry3d::math::{IVector, Pose, Vector};
use parry3d::query::{ContactManifold, DefaultQueryDispatcher, PersistentQueryDispatcher};
use parry3d::shape::{Cuboid, Voxels};

#[test]
fn voxels_contact_manifold_at_domain_minimum_does_not_overflow() {
    // A 5x2x5 slab of unit voxels. Every top-layer voxel has a neighbor below
    // it, so its canonical shape extends past the domain minimum along -y
    // (and along -x/-z for the interior ones), which is exactly the situation
    // that overflowed the workspace-key computation.
    let mut coords = Vec::new();
    for i in 0..5 {
        for j in 0..2 {
            for k in 0..5 {
                coords.push(IVector::new(i, j, k));
            }
        }
    }
    let voxels = Voxels::new(Vector::new(1.0, 1.0, 1.0), &coords);

    // A small cuboid resting on the center of the slab's top face (y = 2),
    // penetrating by 0.05. The dispatcher routes cuboid-vs-voxels through
    // `contact_manifolds_voxels_shape`, which calls
    // `CanonicalVoxelShape::from_voxel` for each candidate voxel.
    let cuboid = Cuboid::new(Vector::new(0.4, 0.4, 0.4));
    let pos12 = Pose::translation(2.5, 2.35, 2.5);

    let mut manifolds: Vec<ContactManifold<(), ()>> = Vec::new();
    let mut workspace = None;

    // Run two frames to also exercise the workspace-reuse (occupied entry) path.
    for _ in 0..2 {
        DefaultQueryDispatcher
            .contact_manifolds(
                &pos12,
                &voxels,
                &cuboid,
                0.05,
                &mut manifolds,
                &mut workspace,
            )
            .expect("the voxels/cuboid pair must be supported");
    }

    // No overflow panic reached this point; also check the result is sane.
    let deepest = manifolds
        .iter()
        .flat_map(|m| m.points.iter())
        .map(|pt| pt.dist)
        .fold(f32::MAX, f32::min);
    assert!(
        (deepest + 0.05).abs() < 1.0e-3,
        "expected a contact at dist ~ -0.05, got {deepest}"
    );
}
