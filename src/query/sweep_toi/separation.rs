//! Separation functions used by the conservative-advancement time-of-impact loop.
//!
//! A separation function measures the signed separation of two swept proxies along a
//! progressively chosen axis, either re-picking the deepest support points
//! ([`SeparationFunction::find_min_separation`]) or holding witness features fixed for root
//! finding ([`SeparationFunction::evaluate`]).

use super::proxy_distance::SimplexCache;
use super::sweep::{inv_rotate_vec, rotate_vec, Sweep};
use super::toi_proxy::ToiProxy;
use crate::math::{Real, Vector};

/// Sentinel index meaning "no witness index" (face types keep a fixed local point instead).
pub(crate) const INVALID_INDEX: u32 = u32::MAX;

#[cfg(feature = "dim2")]
#[derive(Copy, Clone, PartialEq)]
enum SeparationType {
    Points,
    FaceA,
    FaceB,
}

#[cfg(feature = "dim2")]
pub(crate) struct SeparationFunction<'a, 'b> {
    proxy_a: &'a ToiProxy<'b>,
    proxy_b: &'a ToiProxy<'b>,
    sweep_a: Sweep,
    sweep_b: Sweep,
    local_point: Vector,
    axis: Vector,
    ty: SeparationType,
}

#[cfg(feature = "dim2")]
impl<'a, 'b> SeparationFunction<'a, 'b> {
    pub fn new(
        cache: &SimplexCache,
        proxy_a: &'a ToiProxy<'b>,
        sweep_a: &Sweep,
        proxy_b: &'a ToiProxy<'b>,
        sweep_b: &Sweep,
        _world_normal: Vector,
        t1: Real,
    ) -> Self {
        let count = cache.count as usize;
        debug_assert!(0 < count && count < 3);

        let xf_a = sweep_a.transform_at(t1);
        let xf_b = sweep_b.transform_at(t1);

        if count == 1 {
            let local_point_a = proxy_a.points()[cache.index_a[0] as usize];
            let local_point_b = proxy_b.points()[cache.index_b[0] as usize];
            let point_a = xf_a.transform_point(local_point_a);
            let point_b = xf_b.transform_point(local_point_b);
            return Self {
                proxy_a,
                proxy_b,
                sweep_a: *sweep_a,
                sweep_b: *sweep_b,
                local_point: Vector::ZERO,
                axis: (point_b - point_a).normalize_or_zero(),
                ty: SeparationType::Points,
            };
        }

        if cache.index_a[0] == cache.index_a[1] {
            // Two points on B and one on A.
            let local_point_b1 = proxy_b.points()[cache.index_b[0] as usize];
            let local_point_b2 = proxy_b.points()[cache.index_b[1] as usize];

            // Perpendicular of the edge: cross(edge, 1.0).
            let edge = local_point_b2 - local_point_b1;
            let mut axis = Vector::new(edge.y, -edge.x).normalize_or_zero();
            let normal = rotate_vec(&xf_b.rotation, axis);

            let local_point = 0.5 * (local_point_b1 + local_point_b2);
            let point_b = xf_b.transform_point(local_point);

            let local_point_a = proxy_a.points()[cache.index_a[0] as usize];
            let point_a = xf_a.transform_point(local_point_a);

            if (point_a - point_b).dot(normal) < 0.0 {
                axis = -axis;
            }

            Self {
                proxy_a,
                proxy_b,
                sweep_a: *sweep_a,
                sweep_b: *sweep_b,
                local_point,
                axis,
                ty: SeparationType::FaceB,
            }
        } else {
            // Two points on A and one or two points on B.
            let local_point_a1 = proxy_a.points()[cache.index_a[0] as usize];
            let local_point_a2 = proxy_a.points()[cache.index_a[1] as usize];

            let edge = local_point_a2 - local_point_a1;
            let mut axis = Vector::new(edge.y, -edge.x).normalize_or_zero();
            let normal = rotate_vec(&xf_a.rotation, axis);

            let local_point = 0.5 * (local_point_a1 + local_point_a2);
            let point_a = xf_a.transform_point(local_point);

            let local_point_b = proxy_b.points()[cache.index_b[0] as usize];
            let point_b = xf_b.transform_point(local_point_b);

            if (point_b - point_a).dot(normal) < 0.0 {
                axis = -axis;
            }

            Self {
                proxy_a,
                proxy_b,
                sweep_a: *sweep_a,
                sweep_b: *sweep_b,
                local_point,
                axis,
                ty: SeparationType::FaceA,
            }
        }
    }

    /// Finds the deepest support points at time `t` and returns their signed separation.
    pub fn find_min_separation(&self, t: Real) -> (Real, u32, u32) {
        let xf_a = self.sweep_a.transform_at(t);
        let xf_b = self.sweep_b.transform_at(t);

        match self.ty {
            SeparationType::Points => {
                let axis_a = inv_rotate_vec(&xf_a.rotation, self.axis);
                let axis_b = inv_rotate_vec(&xf_b.rotation, -self.axis);

                let index_a = self.proxy_a.support(axis_a);
                let index_b = self.proxy_b.support(axis_b);

                let point_a = xf_a.transform_point(self.proxy_a.points()[index_a as usize]);
                let point_b = xf_b.transform_point(self.proxy_b.points()[index_b as usize]);

                ((point_b - point_a).dot(self.axis), index_a, index_b)
            }
            SeparationType::FaceA => {
                let normal = rotate_vec(&xf_a.rotation, self.axis);
                let point_a = xf_a.transform_point(self.local_point);

                let axis_b = inv_rotate_vec(&xf_b.rotation, -normal);
                let index_b = self.proxy_b.support(axis_b);
                let point_b = xf_b.transform_point(self.proxy_b.points()[index_b as usize]);

                ((point_b - point_a).dot(normal), INVALID_INDEX, index_b)
            }
            SeparationType::FaceB => {
                let normal = rotate_vec(&xf_b.rotation, self.axis);
                let point_b = xf_b.transform_point(self.local_point);

                let axis_a = inv_rotate_vec(&xf_a.rotation, -normal);
                let index_a = self.proxy_a.support(axis_a);
                let point_a = xf_a.transform_point(self.proxy_a.points()[index_a as usize]);

                ((point_a - point_b).dot(normal), index_a, INVALID_INDEX)
            }
        }
    }

    /// Evaluates the separation of fixed witness features at time `t`.
    pub fn evaluate(&self, index_a: u32, index_b: u32, t: Real) -> Real {
        let xf_a = self.sweep_a.transform_at(t);
        let xf_b = self.sweep_b.transform_at(t);

        match self.ty {
            SeparationType::Points => {
                let point_a = xf_a.transform_point(self.proxy_a.points()[index_a as usize]);
                let point_b = xf_b.transform_point(self.proxy_b.points()[index_b as usize]);
                (point_b - point_a).dot(self.axis)
            }
            SeparationType::FaceA => {
                let normal = rotate_vec(&xf_a.rotation, self.axis);
                let point_a = xf_a.transform_point(self.local_point);
                let point_b = xf_b.transform_point(self.proxy_b.points()[index_b as usize]);
                (point_b - point_a).dot(normal)
            }
            SeparationType::FaceB => {
                let normal = rotate_vec(&xf_b.rotation, self.axis);
                let point_b = xf_b.transform_point(self.local_point);
                let point_a = xf_a.transform_point(self.proxy_a.points()[index_a as usize]);
                (point_a - point_b).dot(normal)
            }
        }
    }

    /// Whether this function uses the edge/edge cross-product axis (3D only; always false in
    /// 2D).
    pub fn uses_edge_axis(&self) -> bool {
        false
    }

    /// Converts an edge-axis function to a fixed world axis (3D only; no-op in 2D).
    pub fn force_fixed_axis(&mut self, _t: Real) {}
}

// ==================== 3D ====================

#[cfg(feature = "dim3")]
#[derive(Copy, Clone, PartialEq)]
enum SeparationType {
    Vertices,
    Edges,
    FaceA,
    FaceB,
}

#[cfg(feature = "dim3")]
pub(crate) struct SeparationFunction<'a, 'b> {
    proxy_a: &'a ToiProxy<'b>,
    proxy_b: &'a ToiProxy<'b>,
    sweep_a: Sweep,
    sweep_b: Sweep,
    // These are associated with different bodies depending on the separation function type.
    // It could be two local vectors/points on the same body (for example, both on bodyA).
    witness1: Vector,
    witness2: Vector,
    ty: SeparationType,
}

#[cfg(feature = "dim3")]
fn unique_count(count: usize, indices: &[u32; 3]) -> usize {
    match count {
        1 => 1,
        2 => {
            if indices[0] != indices[1] {
                2
            } else {
                1
            }
        }
        3 => {
            if indices[0] != indices[1] && indices[0] != indices[2] && indices[1] != indices[2] {
                3
            } else if indices[0] == indices[1] && indices[0] == indices[2] {
                1
            } else {
                2
            }
        }
        _ => unreachable!(),
    }
}

// This checks if the cross product of two edges switches direction over the sweep.
#[cfg(feature = "dim3")]
fn check_fast_edges(
    sweep_a: &Sweep,
    local_edge_a: Vector,
    sweep_b: &Sweep,
    local_edge_b: Vector,
    axis0: Vector,
) -> bool {
    // By taking the local witness axes we make sure that we get the correct orientations
    // (e.g. if one axis was flipped)!
    let xf_a2 = sweep_a.final_transform();
    let xf_b2 = sweep_b.final_transform();
    let edge_a = rotate_vec(&xf_a2.rotation, local_edge_a);
    let edge_b = rotate_vec(&xf_b2.rotation, local_edge_b);
    let axis = edge_a.cross(edge_b);
    axis.dot(axis0) < 0.0
}

#[cfg(feature = "dim3")]
impl<'a, 'b> SeparationFunction<'a, 'b> {
    pub fn new(
        cache: &SimplexCache,
        proxy_a: &'a ToiProxy<'b>,
        sweep_a: &Sweep,
        proxy_b: &'a ToiProxy<'b>,
        sweep_b: &Sweep,
        world_normal: Vector,
        t1: Real,
    ) -> Self {
        let count = cache.count as usize;
        debug_assert!(1 <= count && count <= 3);

        let mut index_a = cache.index_a;
        let mut index_b = cache.index_b;

        let unique_count_a = unique_count(count, &index_a);
        let unique_count_b = unique_count(count, &index_b);

        let xf_a1 = sweep_a.transform_at(t1);
        let xf_b1 = sweep_b.transform_at(t1);

        let qa = xf_a1.rotation;
        let qb = xf_b1.rotation;

        // Minimize round-off
        let delta_p = xf_b1.translation - xf_a1.translation;

        let mut result = Self {
            proxy_a,
            proxy_b,
            sweep_a: *sweep_a,
            sweep_b: *sweep_b,
            witness1: world_normal,
            witness2: Vector::ZERO,
            ty: SeparationType::Vertices,
        };

        match count {
            1 => {
                // Witness is the world space direction
                result.ty = SeparationType::Vertices;
                result.witness1 = world_normal;
            }
            2 => {
                if unique_count_a == 2 && unique_count_b == 2 {
                    Self::init_edges(
                        &mut result,
                        proxy_a,
                        proxy_b,
                        &index_a,
                        &index_b,
                        &qa,
                        &qb,
                        delta_p,
                        world_normal,
                        0.05,
                    );
                } else {
                    // Vertex versus edge, use world axis witness
                    result.ty = SeparationType::Vertices;
                    result.witness1 = world_normal;
                }
            }
            3 => {
                if unique_count_a == 3 {
                    let va1 = proxy_a.points()[index_a[0] as usize];
                    let va2 = proxy_a.points()[index_a[1] as usize];
                    let va3 = proxy_a.points()[index_a[2] as usize];
                    let mut local_axis_a = (va2 - va1).cross(va3 - va1).normalize_or_zero();
                    let axis_a = rotate_vec(&qa, local_axis_a);

                    let local_point_a = (va1 + va2 + va3) / 3.0;
                    let local_point_b = proxy_b.points()[index_b[0] as usize];
                    let delta =
                        rotate_vec(&qb, local_point_b) - rotate_vec(&qa, local_point_a) + delta_p;

                    if delta.dot(axis_a) < 0.0 {
                        // Make axis point from A to B
                        local_axis_a = -local_axis_a;
                    }

                    // Witness is the local plane of faceA
                    result.ty = SeparationType::FaceA;
                    result.witness1 = local_axis_a;
                    result.witness2 = local_point_a;
                } else if unique_count_b == 3 {
                    let vb1 = proxy_b.points()[index_b[0] as usize];
                    let vb2 = proxy_b.points()[index_b[1] as usize];
                    let vb3 = proxy_b.points()[index_b[2] as usize];
                    let mut local_axis_b = (vb2 - vb1).cross(vb3 - vb1).normalize_or_zero();
                    let axis_b = rotate_vec(&qb, local_axis_b);

                    let local_point_a = proxy_a.points()[index_a[0] as usize];
                    let local_point_b = (vb1 + vb2 + vb3) / 3.0;
                    let delta =
                        rotate_vec(&qa, local_point_a) - rotate_vec(&qb, local_point_b) - delta_p;

                    if delta.dot(axis_b) < 0.0 {
                        // Make axis point from B to A
                        local_axis_b = -local_axis_b;
                    }

                    // Witness is the local plane of faceB
                    result.ty = SeparationType::FaceB;
                    result.witness1 = local_axis_b;
                    result.witness2 = local_point_b;
                } else {
                    debug_assert!(unique_count_a == 2 && unique_count_b == 2);

                    if index_a[0] == index_a[1] {
                        // Make first two indices unique
                        index_a[1] = index_a[2];
                    }
                    if index_b[0] == index_b[1] {
                        // Make first two indices unique
                        index_b[1] = index_b[2];
                    }

                    Self::init_edges(
                        &mut result,
                        proxy_a,
                        proxy_b,
                        &index_a,
                        &index_b,
                        &qa,
                        &qb,
                        delta_p,
                        world_normal,
                        0.005,
                    );
                }
            }
            _ => unreachable!(),
        }

        result
    }

    #[allow(clippy::too_many_arguments)]
    fn init_edges(
        result: &mut Self,
        proxy_a: &'a ToiProxy<'b>,
        proxy_b: &'a ToiProxy<'b>,
        index_a: &[u32; 3],
        index_b: &[u32; 3],
        qa: &crate::math::Rotation,
        qb: &crate::math::Rotation,
        delta_p: Vector,
        world_normal: Vector,
        parallel_tolerance: Real,
    ) {
        // Edge/Edge
        let va1 = proxy_a.points()[index_a[0] as usize];
        let local_edge_a = (proxy_a.points()[index_a[1] as usize] - va1).normalize_or_zero();

        let vb1 = proxy_b.points()[index_b[0] as usize];
        let mut local_edge_b = (proxy_b.points()[index_b[1] as usize] - vb1).normalize_or_zero();

        let edge_a = rotate_vec(qa, local_edge_a);
        let edge_b = rotate_vec(qb, local_edge_b);

        let mut axis = edge_a.cross(edge_b);
        let length_squared = axis.length_squared();

        // Skip near parallel edges: |e1 x e2| = sin(alpha) * |e1| * |e2|
        let tolerance_squared = parallel_tolerance * parallel_tolerance;
        if length_squared < tolerance_squared {
            // The axis is not safe to normalize so we use a world axis instead!
            result.ty = SeparationType::Vertices;
            result.witness1 = world_normal;
            return;
        }

        let delta = rotate_vec(qb, vb1) - rotate_vec(qa, va1) + delta_p;
        if delta.dot(axis) < 0.0 {
            // Make axis point from A to B
            axis = -axis;
            local_edge_b = -local_edge_b;
        }

        // Check for possible sign flip in edge/edge cross product
        if check_fast_edges(
            &result.sweep_a,
            local_edge_a,
            &result.sweep_b,
            local_edge_b,
            axis,
        ) {
            // Not safe to use local edges, fall back to initial world space axis instead
            result.ty = SeparationType::Vertices;
            result.witness1 = axis.normalize_or_zero();
        } else {
            // Edge cross product is safe. This converges faster than a fixed axis.
            result.ty = SeparationType::Edges;
            result.witness1 = local_edge_a;
            result.witness2 = local_edge_b;
        }
    }

    /// Finds the deepest support points at time `t` and returns their signed separation.
    pub fn find_min_separation(&self, t: Real) -> (Real, u32, u32) {
        let xf_a = self.sweep_a.transform_at(t);
        let xf_b = self.sweep_b.transform_at(t);

        match self.ty {
            SeparationType::Vertices => {
                let axis = self.witness1;

                let local_axis_a = inv_rotate_vec(&xf_a.rotation, axis);
                let local_axis_b = inv_rotate_vec(&xf_b.rotation, -axis);

                let index_a = self.proxy_a.support(local_axis_a);
                let index_b = self.proxy_b.support(local_axis_b);

                let delta_p = xf_b.translation - xf_a.translation;
                let local_point_a = self.proxy_a.points()[index_a as usize];
                let local_point_b = self.proxy_b.points()[index_b as usize];
                let delta = rotate_vec(&xf_b.rotation, local_point_b)
                    - rotate_vec(&xf_a.rotation, local_point_a)
                    + delta_p;

                (delta.dot(axis), index_a, index_b)
            }
            SeparationType::Edges => {
                let edge_a = rotate_vec(&xf_a.rotation, self.witness1);
                let edge_b = rotate_vec(&xf_b.rotation, self.witness2);
                let axis = edge_a.cross(edge_b).normalize_or_zero();

                let axis_a = inv_rotate_vec(&xf_a.rotation, axis);
                let index_a = self.proxy_a.support(axis_a);

                let axis_b = inv_rotate_vec(&xf_b.rotation, axis);
                let index_b = self.proxy_b.support(-axis_b);

                let delta_p = xf_b.translation - xf_a.translation;
                let local_point_a = self.proxy_a.points()[index_a as usize];
                let local_point_b = self.proxy_b.points()[index_b as usize];
                let delta = rotate_vec(&xf_b.rotation, local_point_b)
                    - rotate_vec(&xf_a.rotation, local_point_a)
                    + delta_p;

                (delta.dot(axis), index_a, index_b)
            }
            SeparationType::FaceA => {
                let normal = rotate_vec(&xf_a.rotation, self.witness1);
                let point_a = xf_a.transform_point(self.witness2);

                let axis_b = inv_rotate_vec(&xf_b.rotation, normal);
                let index_b = self.proxy_b.support(-axis_b);
                let point_b = xf_b.transform_point(self.proxy_b.points()[index_b as usize]);

                ((point_b - point_a).dot(normal), INVALID_INDEX, index_b)
            }
            SeparationType::FaceB => {
                let normal = rotate_vec(&xf_b.rotation, self.witness1);

                let axis_a = inv_rotate_vec(&xf_a.rotation, normal);
                let index_a = self.proxy_a.support(-axis_a);
                let point_a = xf_a.transform_point(self.proxy_a.points()[index_a as usize]);

                let point_b = xf_b.transform_point(self.witness2);

                ((point_a - point_b).dot(normal), index_a, INVALID_INDEX)
            }
        }
    }

    /// Evaluates the separation of fixed witness features at time `t`.
    pub fn evaluate(&self, index_a: u32, index_b: u32, t: Real) -> Real {
        let xf_a = self.sweep_a.transform_at(t);
        let xf_b = self.sweep_b.transform_at(t);

        match self.ty {
            SeparationType::Vertices => {
                let axis = self.witness1;
                let point_a = xf_a.transform_point(self.proxy_a.points()[index_a as usize]);
                let point_b = xf_b.transform_point(self.proxy_b.points()[index_b as usize]);
                (point_b - point_a).dot(axis)
            }
            SeparationType::Edges => {
                let edge_a = rotate_vec(&xf_a.rotation, self.witness1);
                let edge_b = rotate_vec(&xf_b.rotation, self.witness2);
                let axis = edge_a.cross(edge_b).normalize_or_zero();

                let point_a = xf_a.transform_point(self.proxy_a.points()[index_a as usize]);
                let point_b = xf_b.transform_point(self.proxy_b.points()[index_b as usize]);
                (point_b - point_a).dot(axis)
            }
            SeparationType::FaceA => {
                let axis = rotate_vec(&xf_a.rotation, self.witness1);
                let point_a = xf_a.transform_point(self.witness2);
                let point_b = xf_b.transform_point(self.proxy_b.points()[index_b as usize]);
                (point_b - point_a).dot(axis)
            }
            SeparationType::FaceB => {
                let axis = rotate_vec(&xf_b.rotation, self.witness1);
                let point_a = xf_a.transform_point(self.proxy_a.points()[index_a as usize]);
                let point_b = xf_b.transform_point(self.witness2);
                (point_a - point_b).dot(axis)
            }
        }
    }

    /// Whether this function uses the edge/edge cross-product axis.
    pub fn uses_edge_axis(&self) -> bool {
        self.ty == SeparationType::Edges
    }

    /// Converts an edge-axis function to a frozen world axis evaluated at time `t`.
    pub fn force_fixed_axis(&mut self, t: Real) {
        debug_assert!(self.ty == SeparationType::Edges);

        let xf_a = self.sweep_a.transform_at(t);
        let xf_b = self.sweep_b.transform_at(t);

        let edge_a = rotate_vec(&xf_a.rotation, self.witness1);
        let edge_b = rotate_vec(&xf_b.rotation, self.witness2);
        let axis = edge_a.cross(edge_b).normalize_or_zero();

        self.ty = SeparationType::Vertices;
        self.witness1 = axis;
        self.witness2 = Vector::ZERO;
    }
}
