//! `TriMeshFlags::FIX_INTERNAL_EDGES_TWO_SIDED` keeps the internal-edge fix on both faces of the
//! mesh instead of discarding the contacts coming from the back.

use parry3d::math::{Pose, Real, Vector};
use parry3d::query::{
    ContactManifold, ContactManifoldsWorkspace, DefaultQueryDispatcher, PersistentQueryDispatcher,
};
use parry3d::shape::{Ball, TriMesh, TriMeshFlags};

/// A flat two-quad strip at `y = 0`, normals up.
fn strip() -> (Vec<Vector>, Vec<[u32; 3]>) {
    let vertices = vec![
        Vector::new(-2.0, 0.0, -1.0),
        Vector::new(0.0, 0.0, -1.0),
        Vector::new(2.0, 0.0, -1.0),
        Vector::new(-2.0, 0.0, 1.0),
        Vector::new(0.0, 0.0, 1.0),
        Vector::new(2.0, 0.0, 1.0),
    ];
    let indices = vec![[0, 4, 1], [0, 3, 4], [1, 5, 2], [1, 4, 5]];
    (vertices, indices)
}

/// Deepest contact (normal on the mesh, distance) of a ball at `ball_pos`.
fn deepest(mesh: &TriMesh, ball: &Ball, ball_pos: Vector) -> Option<(Vector, Real)> {
    let pos12 = Pose::from_translation(ball_pos);
    let mut manifolds: Vec<ContactManifold<(), ()>> = Vec::new();
    let mut workspace: Option<ContactManifoldsWorkspace> = None;
    DefaultQueryDispatcher
        .contact_manifolds(&pos12, mesh, ball, 0.0, &mut manifolds, &mut workspace)
        .unwrap();
    let mut best: Option<(Vector, Real)> = None;
    for m in &manifolds {
        for pt in &m.points {
            if best.is_none_or(|(_, d)| pt.dist < d) {
                best = Some((m.local_n1, pt.dist));
            }
        }
    }
    best
}

#[test]
fn two_sided_internal_edges_keep_back_contacts() {
    let (vertices, indices) = strip();
    let ball = Ball::new(0.5);
    let above = Vector::new(0.0, 0.4, 0.0);
    let below = Vector::new(0.0, -0.4, 0.0);

    let one_sided = TriMesh::with_flags(
        vertices.clone(),
        indices.clone(),
        TriMeshFlags::FIX_INTERNAL_EDGES,
    )
    .unwrap();
    let two_sided = TriMesh::with_flags(
        vertices,
        indices,
        TriMeshFlags::FIX_INTERNAL_EDGES_TWO_SIDED,
    )
    .unwrap();
    assert!(two_sided.flags().contains(TriMeshFlags::FIX_INTERNAL_EDGES));

    // Front contacts are identical.
    let (n1, d1) = deepest(&one_sided, &ball, above).unwrap();
    let (n2, d2) = deepest(&two_sided, &ball, above).unwrap();
    assert!(n1.y > 0.99 && n2.y > 0.99);
    assert!((d1 - d2).abs() < 1.0e-6 && (d1 + 0.1).abs() < 1.0e-4);

    // Back contacts are discarded by the one-sided fix and kept (normal flipped) by the two-sided one.
    assert!(deepest(&one_sided, &ball, below).is_none());
    let (n3, d3) = deepest(&two_sided, &ball, below).unwrap();
    assert!(n3.y < -0.99, "n3 = {n3:?}");
    assert!((d3 + 0.1).abs() < 1.0e-4, "d3 = {d3}");
}
