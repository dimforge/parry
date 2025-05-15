use crate::bounding_volume::{Aabb, BoundingVolume};
use crate::math::{Isometry, Real, Vector};
use crate::query::ContactManifold;
use crate::shape::{Shape, VoxelType, Voxels};
use crate::query::contact_manifolds::contact_manifolds_voxels_ball::detect_hit_voxel_ball;
use alloc::vec::Vec;

/// Computes the contact manifold between a convex shape and a ball, both represented as a `Shape` trait-object.
pub fn contact_manifolds_voxels_voxels_shapes<ManifoldData, ContactData>(
    pos12: &Isometry<Real>,
    shape1: &dyn Shape,
    shape2: &dyn Shape,
    prediction: Real,
    manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
) where
    ManifoldData: Default,
    ContactData: Default + Copy,
{
    if let (Some(voxels1), Some(voxels2)) = (shape1.as_voxels(), shape2.as_voxels()) {
        contact_manifolds_voxels_voxels(pos12, voxels1, voxels2, prediction, manifolds);
    }
}

/// Computes the contact manifold between a convex shape and a ball.
pub fn contact_manifolds_voxels_voxels<'a, ManifoldData, ContactData>(
    pos12: &Isometry<Real>,
    voxels1: &'a Voxels,
    voxels2: &'a Voxels,
    prediction: Real,
    manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
) where
    ManifoldData: Default,
    ContactData: Default + Copy,
{
    // TODO: don’t generate one manifold per voxel.
    manifolds.clear();

    let radius1 = voxels1.voxel_size() / 2.0;
    let radius2 = voxels2.voxel_size() / 2.0;

    // FIXME: optimize this.
    let aabb1 = voxels1.local_aabb().loosened(prediction / 2.0);
    let aabb2 = voxels2.local_aabb().loosened(prediction / 2.0);
    let geometry1 = voxels1.primitive_geometry();
    let geometry2 = voxels2.primitive_geometry();

    if let Some((intersection_aabb1, intersection_aabb2)) =
        aabb1.aligned_intersections(pos12, &aabb2)
    {
        let pos21 = pos12.inverse();

        for vox1 in voxels1.voxels_intersecting_local_aabb(&intersection_aabb1) {
            let type1 = vox1.state.voxel_type();
            match type1 {
                #[cfg(feature = "dim2")]
                VoxelType::Vertex => { /* Ok */ }
                #[cfg(feature = "dim3")]
                VoxelType::Vertex | VoxelType::Edge => { /* Ok */ }
                _ => continue,
            }

            let centered_aabb1_2 =
                Aabb::from_half_extents(pos21 * vox1.center, radius1 + Vector::repeat(prediction));

            for vox2 in voxels2.voxels_intersecting_local_aabb(&centered_aabb1_2) {
                let type2 = vox2.state.voxel_type();
                #[cfg(feature = "dim2")]
                match (type1, type2) {
                    (VoxelType::Vertex, VoxelType::Vertex)
                    | (VoxelType::Vertex, VoxelType::Face) => {
                        detect_hit_voxel_ball(
                            pos21, vox2.center, radius2, vox2.state, geometry2, vox1.center, radius1.x,
                            prediction, true, manifolds,
                        );
                    }
                    (VoxelType::Face, VoxelType::Vertex) => {
                        detect_hit_voxel_ball(
                            *pos12, vox1.center, radius1, vox1.state, geometry1, vox2.center, radius2.x,
                            prediction, false, manifolds,
                        );
                    }
                    _ => continue, /* Ignore */
                }

                #[cfg(feature = "dim3")]
                match (type1, type2) {
                    (VoxelType::Vertex, VoxelType::Vertex)
                    | (VoxelType::Vertex, VoxelType::Edge)
                    | (VoxelType::Vertex, VoxelType::Face) => {
                        detect_hit_voxel_ball(
                            pos21, vox2.center, radius2, vox2.state, geometry2, vox1.center, radius1.x,
                            prediction, true, manifolds,
                        );
                    }
                    (VoxelType::Edge, VoxelType::Vertex)
                    | (VoxelType::Face, VoxelType::Vertex)
                    | (VoxelType::Edge, VoxelType::Edge) => {
                        detect_hit_voxel_ball(
                            *pos12, vox1.center, radius1, vox1.state, geometry1, vox2.center, radius2.x,
                            prediction, false, manifolds,
                        );
                    }
                    _ => continue, /* Ignore */
                }
            }
        }

        for vox2 in voxels2.voxels_intersecting_local_aabb(&intersection_aabb2) {
            if vox2.state.voxel_type() != VoxelType::Vertex {
                continue;
            }

            let centered_aabb2_1 =
                Aabb::from_half_extents(pos12 * vox2.center, radius2 + Vector::repeat(prediction));

            for vox1 in voxels1.voxels_intersecting_local_aabb(&centered_aabb2_1) {
                if vox1.state.voxel_type() != VoxelType::Face {
                    continue;
                }

                detect_hit_voxel_ball(
                    *pos12, vox1.center, radius1, vox1.state, geometry1, vox2.center, radius2.x, prediction,
                    false, manifolds,
                );
            }
        }
    }
}
