//! GJK distance between two point-cloud proxies with an index-based simplex cache.
//!
//! The cache warm-starts successive distance queries on the same pair, which is the
//! backbone of the conservative-advancement time-of-impact loop.

use super::toi_proxy::ToiProxy;
use crate::math::{Pose, Real, Vector};

/// Warm-starting simplex cache for [`proxy_distance`].
#[derive(Copy, Clone, Debug, Default)]
pub struct SimplexCache {
    /// Simplex size measure used to detect a stale cache (3D only).
    #[cfg(feature = "dim3")]
    pub metric: Real,
    /// Number of cached simplex vertices (0 to 3).
    pub count: u8,
    /// Cached support indices on the first proxy.
    pub index_a: [u32; 3],
    /// Cached support indices on the second proxy.
    pub index_b: [u32; 3],
}

/// Result of [`proxy_distance`]. All geometric quantities are expressed in the local frame of
/// the first proxy.
#[derive(Copy, Clone, Debug, Default)]
pub struct ProxyDistanceOutput {
    /// Closest point on the first proxy (frame A).
    pub point_a: Vector,
    /// Closest point on the second proxy (frame A).
    pub point_b: Vector,
    /// Separation direction pointing from A to B (frame A). Zero if overlapped.
    pub normal: Vector,
    /// Distance between the two proxies (0 if overlapped).
    pub distance: Real,
    /// Number of GJK iterations used.
    pub iterations: u32,
}

#[derive(Copy, Clone, Default)]
struct SimplexVertex {
    wa: Vector,
    wb: Vector,
    w: Vector,
    a: Real,
    index_a: u32,
    index_b: u32,
}

fn cache_is_valid(cache: &SimplexCache, proxy_a: &ToiProxy, proxy_b: &ToiProxy) -> bool {
    let (na, nb) = (proxy_a.points().len() as u32, proxy_b.points().len() as u32);
    cache.count <= 3
        && cache.index_a[..cache.count as usize].iter().all(|i| *i < na)
        && cache.index_b[..cache.count as usize].iter().all(|i| *i < nb)
}

/// Computes the distance between two point-cloud proxies, warm-started by `cache`.
///
/// `pos12` is the pose of the second proxy relative to the first; the query runs entirely in
/// the first proxy’s local frame. When `use_radii` is `false` the proxy radii are ignored
/// (core-shape distance), which is what the time-of-impact loop uses.
#[cfg(feature = "dim2")]
pub fn proxy_distance(
    pos12: &Pose,
    proxy_a: &ToiProxy,
    proxy_b: &ToiProxy,
    use_radii: bool,
    cache: &mut SimplexCache,
) -> ProxyDistanceOutput {
    let points_a = proxy_a.points();
    let points_b = proxy_b.points();
    let point_b_in_a = |i: u32| pos12.transform_point(points_b[i as usize]);

    let mut output = ProxyDistanceOutput::default();

    // Initialize the simplex from the cache.
    let mut simplex = [SimplexVertex::default(); 3];
    let mut count = if cache_is_valid(cache, proxy_a, proxy_b) {
        cache.count as usize
    } else {
        0
    };
    for i in 0..count {
        let v = &mut simplex[i];
        v.index_a = cache.index_a[i];
        v.index_b = cache.index_b[i];
        v.wa = points_a[v.index_a as usize];
        v.wb = point_b_in_a(v.index_b);
        v.w = v.wa - v.wb;
        // Invalid coefficient; set by the simplex solvers.
        v.a = -1.0;
    }

    if count == 0 {
        let v = &mut simplex[0];
        v.index_a = 0;
        v.index_b = 0;
        v.wa = points_a[0];
        v.wb = point_b_in_a(0);
        v.w = v.wa - v.wb;
        v.a = 1.0;
        count = 1;
    }

    let mut non_unit_normal = Vector::ZERO;
    let mut save_a = [0u32; 3];
    let mut save_b = [0u32; 3];

    // Main iteration loop. All computations are done in frame A.
    const MAX_ITERATIONS: u32 = 20;
    let mut iteration = 0;
    while iteration < MAX_ITERATIONS {
        // Copy simplex indices so we can identify duplicates.
        let save_count = count;
        for i in 0..save_count {
            save_a[i] = simplex[i].index_a;
            save_b[i] = simplex[i].index_b;
        }

        let d = match count {
            1 => -simplex[0].w,
            2 => solve_simplex2(&mut simplex, &mut count),
            3 => solve_simplex3(&mut simplex, &mut count),
            _ => unreachable!(),
        };

        // If we have 3 points, then the origin is in the corresponding triangle.
        if count == 3 {
            let (pa, pb) = witness_points(&simplex, count);
            output.point_a = pa;
            output.point_b = pb;
            output.iterations = iteration;
            return output;
        }

        // Ensure the search direction is numerically fit; a degenerate direction means the
        // origin is contained by a segment, i.e. the shapes are overlapped.
        if d.dot(d) < Real::EPSILON * Real::EPSILON {
            let (pa, pb) = witness_points(&simplex, count);
            output.point_a = pa;
            output.point_b = pb;
            output.iterations = iteration;
            return output;
        }

        non_unit_normal = d;

        // Compute a tentative new simplex vertex using support points:
        // support = support(a, d) - support(b, -d).
        let index_a = proxy_a.support(d);
        let index_b = proxy_b.support(pos12.rotation.inverse_transform_vector(-d));
        let vertex = &mut simplex[count];
        vertex.index_a = index_a;
        vertex.wa = points_a[index_a as usize];
        vertex.index_b = index_b;
        vertex.wb = point_b_in_a(index_b);
        vertex.w = vertex.wa - vertex.wb;

        // Iteration count is equated to the number of support point calls.
        iteration += 1;

        // Check for duplicate support points. This is the main termination criteria.
        let duplicate = (0..save_count).any(|i| index_a == save_a[i] && index_b == save_b[i]);
        if duplicate {
            break;
        }

        // New vertex is valid and needed.
        count += 1;
    }

    let normal = non_unit_normal.normalize_or_zero();
    let (pa, pb) = witness_points(&simplex, count);
    output.normal = normal;
    output.distance = pa.distance(pb);
    output.point_a = pa;
    output.point_b = pb;
    output.iterations = iteration;

    // Cache the simplex.
    cache.count = count as u8;
    for i in 0..count {
        cache.index_a[i] = simplex[i].index_a;
        cache.index_b[i] = simplex[i].index_b;
    }

    // Apply radii if requested.
    if use_radii {
        let radius_a = proxy_a.radius;
        let radius_b = proxy_b.radius;
        output.distance = (output.distance - radius_a - radius_b).max(0.0);

        // Keep closest points on perimeter even if overlapped, this way the points move
        // smoothly.
        output.point_a += radius_a * normal;
        output.point_b -= radius_b * normal;
    }

    output
}

#[cfg(feature = "dim2")]
fn witness_points(simplex: &[SimplexVertex; 3], count: usize) -> (Vector, Vector) {
    match count {
        1 => (simplex[0].wa, simplex[0].wb),
        2 => (
            simplex[0].a * simplex[0].wa + simplex[1].a * simplex[1].wa,
            simplex[0].a * simplex[0].wb + simplex[1].a * simplex[1].wb,
        ),
        3 => {
            let pa = simplex[0].a * simplex[0].wa
                + simplex[1].a * simplex[1].wa
                + simplex[2].a * simplex[2].wa;
            (pa, pa)
        }
        _ => unreachable!(),
    }
}

// Returns a vector pointing towards the origin, reducing the simplex if needed.
#[cfg(feature = "dim2")]
fn solve_simplex2(simplex: &mut [SimplexVertex; 3], count: &mut usize) -> Vector {
    let w1 = simplex[0].w;
    let w2 = simplex[1].w;
    let e12 = w2 - w1;

    // w1 region
    let d12_2 = -w1.dot(e12);
    if d12_2 <= 0.0 {
        // a2 <= 0, so we clamp it to 0
        simplex[0].a = 1.0;
        *count = 1;
        return -w1;
    }

    // w2 region
    let d12_1 = w2.dot(e12);
    if d12_1 <= 0.0 {
        // a1 <= 0, so we clamp it to 0
        simplex[1].a = 1.0;
        *count = 1;
        simplex[0] = simplex[1];
        return -w2;
    }

    // Must be in e12 region.
    let inv_d12 = 1.0 / (d12_1 + d12_2);
    simplex[0].a = d12_1 * inv_d12;
    simplex[1].a = d12_2 * inv_d12;
    *count = 2;
    // cross(cross(w1 + w2, e12), e12)
    let s = (w1 + w2).perp_dot(e12);
    Vector::new(-s * e12.y, s * e12.x)
}

#[cfg(feature = "dim2")]
fn solve_simplex3(simplex: &mut [SimplexVertex; 3], count: &mut usize) -> Vector {
    let w1 = simplex[0].w;
    let w2 = simplex[1].w;
    let w3 = simplex[2].w;

    // Edge12: a3 = 0
    let e12 = w2 - w1;
    let w1e12 = w1.dot(e12);
    let w2e12 = w2.dot(e12);
    let d12_1 = w2e12;
    let d12_2 = -w1e12;

    // Edge13: a2 = 0
    let e13 = w3 - w1;
    let w1e13 = w1.dot(e13);
    let w3e13 = w3.dot(e13);
    let d13_1 = w3e13;
    let d13_2 = -w1e13;

    // Edge23: a1 = 0
    let e23 = w3 - w2;
    let w2e23 = w2.dot(e23);
    let w3e23 = w3.dot(e23);
    let d23_1 = w3e23;
    let d23_2 = -w2e23;

    // Triangle123
    let n123 = e12.perp_dot(e13);
    let d123_1 = n123 * w2.perp_dot(w3);
    let d123_2 = n123 * w3.perp_dot(w1);
    let d123_3 = n123 * w1.perp_dot(w2);

    // w1 region
    if d12_2 <= 0.0 && d13_2 <= 0.0 {
        simplex[0].a = 1.0;
        *count = 1;
        return -w1;
    }

    // e12
    if d12_1 > 0.0 && d12_2 > 0.0 && d123_3 <= 0.0 {
        let inv_d12 = 1.0 / (d12_1 + d12_2);
        simplex[0].a = d12_1 * inv_d12;
        simplex[1].a = d12_2 * inv_d12;
        *count = 2;
        let s = (w1 + w2).perp_dot(e12);
        return Vector::new(-s * e12.y, s * e12.x);
    }

    // e13
    if d13_1 > 0.0 && d13_2 > 0.0 && d123_2 <= 0.0 {
        let inv_d13 = 1.0 / (d13_1 + d13_2);
        simplex[0].a = d13_1 * inv_d13;
        simplex[2].a = d13_2 * inv_d13;
        *count = 2;
        simplex[1] = simplex[2];
        let s = (w1 + w3).perp_dot(e13);
        return Vector::new(-s * e13.y, s * e13.x);
    }

    // w2 region
    if d12_1 <= 0.0 && d23_2 <= 0.0 {
        simplex[1].a = 1.0;
        *count = 1;
        simplex[0] = simplex[1];
        return -w2;
    }

    // w3 region
    if d13_1 <= 0.0 && d23_1 <= 0.0 {
        simplex[2].a = 1.0;
        *count = 1;
        simplex[0] = simplex[2];
        return -w3;
    }

    // e23
    if d23_1 > 0.0 && d23_2 > 0.0 && d123_1 <= 0.0 {
        let inv_d23 = 1.0 / (d23_1 + d23_2);
        simplex[1].a = d23_1 * inv_d23;
        simplex[2].a = d23_2 * inv_d23;
        *count = 2;
        simplex[0] = simplex[2];
        let s = (w2 + w3).perp_dot(e23);
        return Vector::new(-s * e23.y, s * e23.x);
    }

    // Must be in triangle123
    let inv_d123 = 1.0 / (d123_1 + d123_2 + d123_3);
    simplex[0].a = d123_1 * inv_d123;
    simplex[1].a = d123_2 * inv_d123;
    simplex[2].a = d123_3 * inv_d123;
    *count = 3;

    // No search direction
    Vector::ZERO
}

// ==================== 3D ====================

#[cfg(feature = "dim3")]
const MAX_GJK_ITERATIONS: u32 = 32;

#[cfg(feature = "dim3")]
fn barycentric_coords_edge(a: Vector, b: Vector) -> [Real; 3] {
    let ab = b - a;
    // Last element is divisor
    [b.dot(ab), -a.dot(ab), ab.dot(ab)]
}

#[cfg(feature = "dim3")]
fn barycentric_coords_tri(a: Vector, b: Vector, c: Vector) -> [Real; 4] {
    let ab = b - a;
    let ac = c - a;

    let b_x_c = b.cross(c);
    let c_x_a = c.cross(a);
    let a_x_b = a.cross(b);

    let ab_x_ac = ab.cross(ac);

    // Last element is divisor
    [
        b_x_c.dot(ab_x_ac),
        c_x_a.dot(ab_x_ac),
        a_x_b.dot(ab_x_ac),
        ab_x_ac.dot(ab_x_ac),
    ]
}

#[cfg(feature = "dim3")]
fn scalar_triple_product(a: Vector, b: Vector, c: Vector) -> Real {
    a.cross(b).dot(c)
}

#[cfg(feature = "dim3")]
fn barycentric_coords_tet(a: Vector, b: Vector, c: Vector, d: Vector) -> [Real; 5] {
    let ab = b - a;
    let ac = c - a;
    let ad = d - a;

    // Last element is divisor (forced to be positive)
    let divisor = scalar_triple_product(ab, ac, ad);
    let sign = if divisor < 0.0 { -1.0 } else { 1.0 };

    [
        sign * scalar_triple_product(b, c, d),
        sign * scalar_triple_product(a, d, c),
        sign * scalar_triple_product(a, b, d),
        sign * scalar_triple_product(a, c, b),
        sign * divisor,
    ]
}

#[cfg(feature = "dim3")]
fn simplex_metric(simplex: &[SimplexVertex; 4], count: usize) -> Real {
    match count {
        1 => 0.0,
        2 => simplex[0].w.distance(simplex[1].w),
        3 => {
            let a = simplex[0].w;
            let b = simplex[1].w;
            let c = simplex[2].w;
            (b - a).cross(c - a).length() / 2.0
        }
        4 => {
            let a = simplex[0].w;
            let b = simplex[1].w;
            let c = simplex[2].w;
            let d = simplex[3].w;
            scalar_triple_product(b - a, c - a, d - a) / 6.0
        }
        _ => unreachable!(),
    }
}

#[cfg(feature = "dim3")]
fn witness_points3(simplex: &[SimplexVertex; 4], count: usize) -> (Vector, Vector) {
    let vs = simplex;
    match count {
        1 => (vs[0].wa, vs[0].wb),
        2 => (
            vs[0].a * vs[0].wa + vs[1].a * vs[1].wa,
            vs[0].a * vs[0].wb + vs[1].a * vs[1].wb,
        ),
        3 => (
            vs[0].a * vs[0].wa + vs[1].a * vs[1].wa + vs[2].a * vs[2].wa,
            vs[0].a * vs[0].wb + vs[1].a * vs[1].wb + vs[2].a * vs[2].wb,
        ),
        4 => {
            // Force identical points and *zero* distance
            let sum = vs[0].a * vs[0].wa
                + vs[1].a * vs[1].wa
                + vs[2].a * vs[2].wa
                + vs[3].a * vs[3].wa;
            (sum, sum)
        }
        _ => unreachable!(),
    }
}

/// Solves the 2-simplex. Returns `false` when the barycentric divisor degenerates.
#[cfg(feature = "dim3")]
fn solve_simplex2_3d(simplex: &mut [SimplexVertex; 4], count: &mut usize) -> bool {
    let a = simplex[0].w;
    let b = simplex[1].w;
    let ab = b - a;

    let divisor = ab.dot(ab);
    let u = b.dot(ab);
    let v = -a.dot(ab);

    // V(A)
    if v <= 0.0 {
        *count = 1;
        simplex[0].a = 1.0;
        return true;
    }

    // V(B)
    if u <= 0.0 {
        *count = 1;
        simplex[0] = simplex[1];
        simplex[0].a = 1.0;
        return true;
    }

    // Edge region
    if divisor <= 0.0 {
        return false;
    }

    let denominator = 1.0 / divisor;
    simplex[0].a = denominator * u;
    simplex[1].a = denominator * v;
    true
}

#[cfg(feature = "dim3")]
fn solve_simplex3_3d(simplex: &mut [SimplexVertex; 4], count: &mut usize) -> bool {
    let v1 = simplex[0];
    let v2 = simplex[1];
    let v3 = simplex[2];

    let w_ab = barycentric_coords_edge(v1.w, v2.w);
    let w_bc = barycentric_coords_edge(v2.w, v3.w);
    let w_ca = barycentric_coords_edge(v3.w, v1.w);

    // VR(A)
    if w_ab[1] <= 0.0 && w_ca[0] <= 0.0 {
        *count = 1;
        simplex[0] = v1;
        simplex[0].a = 1.0;
        return true;
    }

    // VR(B)
    if w_bc[1] <= 0.0 && w_ab[0] <= 0.0 {
        *count = 1;
        simplex[0] = v2;
        simplex[0].a = 1.0;
        return true;
    }

    // VR(C)
    if w_ca[1] <= 0.0 && w_bc[0] <= 0.0 {
        *count = 1;
        simplex[0] = v3;
        simplex[0].a = 1.0;
        return true;
    }

    let w_abc = barycentric_coords_tri(v1.w, v2.w, v3.w);

    // VR(AB)
    if w_abc[2] <= 0.0 && w_ab[0] > 0.0 && w_ab[1] > 0.0 {
        *count = 2;
        simplex[0] = v1;
        simplex[1] = v2;
        let divisor = w_ab[2];
        if divisor <= 0.0 {
            return false;
        }
        simplex[0].a = w_ab[0] / divisor;
        simplex[1].a = w_ab[1] / divisor;
        return true;
    }

    // VR(BC)
    if w_abc[0] <= 0.0 && w_bc[0] > 0.0 && w_bc[1] > 0.0 {
        *count = 2;
        simplex[0] = v2;
        simplex[1] = v3;
        let divisor = w_bc[2];
        if divisor <= 0.0 {
            return false;
        }
        simplex[0].a = w_bc[0] / divisor;
        simplex[1].a = w_bc[1] / divisor;
        return true;
    }

    // VR(CA)
    if w_abc[1] <= 0.0 && w_ca[0] > 0.0 && w_ca[1] > 0.0 {
        *count = 2;
        simplex[0] = v3;
        simplex[1] = v1;
        let divisor = w_ca[2];
        if divisor <= 0.0 {
            return false;
        }
        simplex[0].a = w_ca[0] / divisor;
        simplex[1].a = w_ca[1] / divisor;
        return true;
    }

    // Face region
    let divisor = w_abc[3];
    if divisor <= 0.0 {
        return false;
    }

    // VR(ABC)
    simplex[0].a = w_abc[0] / divisor;
    simplex[1].a = w_abc[1] / divisor;
    simplex[2].a = w_abc[2] / divisor;
    true
}

#[cfg(feature = "dim3")]
fn solve_simplex4_3d(simplex: &mut [SimplexVertex; 4], count: &mut usize) -> bool {
    let vertex_a = simplex[0];
    let vertex_b = simplex[1];
    let vertex_c = simplex[2];
    let vertex_d = simplex[3];

    let w_ab = barycentric_coords_edge(vertex_a.w, vertex_b.w);
    let w_ac = barycentric_coords_edge(vertex_a.w, vertex_c.w);
    let w_ad = barycentric_coords_edge(vertex_a.w, vertex_d.w);
    let w_bc = barycentric_coords_edge(vertex_b.w, vertex_c.w);
    let w_cd = barycentric_coords_edge(vertex_c.w, vertex_d.w);
    let w_db = barycentric_coords_edge(vertex_d.w, vertex_b.w);

    // VR(A)
    if w_ab[1] <= 0.0 && w_ac[1] <= 0.0 && w_ad[1] <= 0.0 {
        *count = 1;
        simplex[0] = vertex_a;
        simplex[0].a = 1.0;
        return true;
    }

    // VR(B)
    if w_ab[0] <= 0.0 && w_db[0] <= 0.0 && w_bc[1] <= 0.0 {
        *count = 1;
        simplex[0] = vertex_b;
        simplex[0].a = 1.0;
        return true;
    }

    // VR(C)
    if w_ac[0] <= 0.0 && w_bc[0] <= 0.0 && w_cd[1] <= 0.0 {
        *count = 1;
        simplex[0] = vertex_c;
        simplex[0].a = 1.0;
        return true;
    }

    // VR(D)
    if w_ad[0] <= 0.0 && w_cd[0] <= 0.0 && w_db[1] <= 0.0 {
        *count = 1;
        simplex[0] = vertex_d;
        simplex[0].a = 1.0;
        return true;
    }

    let w_acb = barycentric_coords_tri(vertex_a.w, vertex_c.w, vertex_b.w);
    let w_abd = barycentric_coords_tri(vertex_a.w, vertex_b.w, vertex_d.w);
    let w_adc = barycentric_coords_tri(vertex_a.w, vertex_d.w, vertex_c.w);
    let w_bcd = barycentric_coords_tri(vertex_b.w, vertex_c.w, vertex_d.w);

    // VR(AB)
    if w_abd[2] <= 0.0 && w_acb[1] <= 0.0 && w_ab[0] > 0.0 && w_ab[1] > 0.0 {
        *count = 2;
        simplex[0] = vertex_a;
        simplex[1] = vertex_b;
        let divisor = w_ab[2];
        if divisor <= 0.0 {
            return false;
        }
        simplex[0].a = w_ab[0] / divisor;
        simplex[1].a = w_ab[1] / divisor;
        return true;
    }

    // VR(AC)
    if w_acb[2] <= 0.0 && w_adc[1] <= 0.0 && w_ac[0] > 0.0 && w_ac[1] > 0.0 {
        *count = 2;
        simplex[0] = vertex_a;
        simplex[1] = vertex_c;
        let divisor = w_ac[2];
        if divisor <= 0.0 {
            return false;
        }
        simplex[0].a = w_ac[0] / divisor;
        simplex[1].a = w_ac[1] / divisor;
        return true;
    }

    // VR(AD)
    if w_adc[2] <= 0.0 && w_abd[1] <= 0.0 && w_ad[0] > 0.0 && w_ad[1] > 0.0 {
        *count = 2;
        simplex[0] = vertex_a;
        simplex[1] = vertex_d;
        let divisor = w_ad[2];
        if divisor <= 0.0 {
            return false;
        }
        simplex[0].a = w_ad[0] / divisor;
        simplex[1].a = w_ad[1] / divisor;
        return true;
    }

    // VR(BC)
    if w_acb[0] <= 0.0 && w_bcd[2] <= 0.0 && w_bc[0] > 0.0 && w_bc[1] > 0.0 {
        *count = 2;
        simplex[0] = vertex_b;
        simplex[1] = vertex_c;
        let divisor = w_bc[2];
        if divisor <= 0.0 {
            return false;
        }
        simplex[0].a = w_bc[0] / divisor;
        simplex[1].a = w_bc[1] / divisor;
        return true;
    }

    // VR(CD)
    if w_adc[0] <= 0.0 && w_bcd[0] <= 0.0 && w_cd[0] > 0.0 && w_cd[1] > 0.0 {
        *count = 2;
        simplex[0] = vertex_c;
        simplex[1] = vertex_d;
        let divisor = w_cd[2];
        if divisor <= 0.0 {
            return false;
        }
        simplex[0].a = w_cd[0] / divisor;
        simplex[1].a = w_cd[1] / divisor;
        return true;
    }

    // VR(DB)
    if w_abd[0] <= 0.0 && w_bcd[1] <= 0.0 && w_db[0] > 0.0 && w_db[1] > 0.0 {
        *count = 2;
        simplex[0] = vertex_d;
        simplex[1] = vertex_b;
        let divisor = w_db[2];
        if divisor <= 0.0 {
            return false;
        }
        simplex[0].a = w_db[0] / divisor;
        simplex[1].a = w_db[1] / divisor;
        return true;
    }

    let w_abcd = barycentric_coords_tet(vertex_a.w, vertex_b.w, vertex_c.w, vertex_d.w);

    // VR(ACB)
    if w_abcd[3] < 0.0 && w_acb[0] > 0.0 && w_acb[1] > 0.0 && w_acb[2] > 0.0 {
        *count = 3;
        simplex[0] = vertex_a;
        simplex[1] = vertex_c;
        simplex[2] = vertex_b;
        let divisor = w_acb[3];
        if divisor <= 0.0 {
            return false;
        }
        simplex[0].a = w_acb[0] / divisor;
        simplex[1].a = w_acb[1] / divisor;
        simplex[2].a = w_acb[2] / divisor;
        return true;
    }

    // VR(ABD)
    if w_abcd[2] < 0.0 && w_abd[0] > 0.0 && w_abd[1] > 0.0 && w_abd[2] > 0.0 {
        *count = 3;
        simplex[0] = vertex_a;
        simplex[1] = vertex_b;
        simplex[2] = vertex_d;
        let divisor = w_abd[3];
        if divisor <= 0.0 {
            return false;
        }
        simplex[0].a = w_abd[0] / divisor;
        simplex[1].a = w_abd[1] / divisor;
        simplex[2].a = w_abd[2] / divisor;
        return true;
    }

    // VR(ADC)
    if w_abcd[1] < 0.0 && w_adc[0] > 0.0 && w_adc[1] > 0.0 && w_adc[2] > 0.0 {
        *count = 3;
        simplex[0] = vertex_a;
        simplex[1] = vertex_d;
        simplex[2] = vertex_c;
        let divisor = w_adc[3];
        if divisor <= 0.0 {
            return false;
        }
        simplex[0].a = w_adc[0] / divisor;
        simplex[1].a = w_adc[1] / divisor;
        simplex[2].a = w_adc[2] / divisor;
        return true;
    }

    // VR(BCD)
    if w_abcd[0] < 0.0 && w_bcd[0] > 0.0 && w_bcd[1] > 0.0 && w_bcd[2] > 0.0 {
        *count = 3;
        simplex[0] = vertex_b;
        simplex[1] = vertex_c;
        simplex[2] = vertex_d;
        let divisor = w_bcd[3];
        if divisor <= 0.0 {
            return false;
        }
        simplex[0].a = w_bcd[0] / divisor;
        simplex[1].a = w_bcd[1] / divisor;
        simplex[2].a = w_bcd[2] / divisor;
        return true;
    }

    // *** Inside tetrahedron ***
    let divisor = w_abcd[4];
    if divisor <= 0.0 {
        return false;
    }

    // VR(ABCD)
    simplex[0].a = w_abcd[0] / divisor;
    simplex[1].a = w_abcd[1] / divisor;
    simplex[2].a = w_abcd[2] / divisor;
    simplex[3].a = w_abcd[3] / divisor;
    true
}

/// Computes the distance between two point-cloud proxies, warm-started by `cache`.
///
/// `pos12` is the pose of the second proxy relative to the first; the query runs entirely in
/// the first proxy’s local frame. When `use_radii` is `false` the proxy radii are ignored
/// (core-shape distance), which is what the time-of-impact loop uses.
#[cfg(feature = "dim3")]
pub fn proxy_distance(
    pos12: &Pose,
    proxy_a: &ToiProxy,
    proxy_b: &ToiProxy,
    use_radii: bool,
    cache: &mut SimplexCache,
) -> ProxyDistanceOutput {
    let points_a = proxy_a.points();
    let points_b = proxy_b.points();
    let point_b_in_a = |i: u32| pos12.transform_point(points_b[i as usize]);

    let mut output = ProxyDistanceOutput::default();

    // Compute initial simplex from cache. Note that in 3D the CSO points are w = wB - wA
    // (the opposite of the 2D convention).
    let mut simplex = [SimplexVertex::default(); 4];
    let mut count = if cache_is_valid(cache, proxy_a, proxy_b) {
        cache.count as usize
    } else {
        0
    };
    for i in 0..count {
        let v = &mut simplex[i];
        v.index_a = cache.index_a[i];
        v.index_b = cache.index_b[i];
        v.wa = points_a[v.index_a as usize];
        v.wb = point_b_in_a(v.index_b);
        v.w = v.wb - v.wa;
        v.a = 0.0;
    }

    // Compute the new simplex metric; if it is substantially different than the old metric
    // flush the simplex.
    if count > 0 {
        let metric1 = cache.metric;
        let metric2 = simplex_metric(&simplex, count);
        if 2.0 * metric1 < metric2 || metric2 < 0.5 * metric1 || metric2 < Real::EPSILON {
            count = 0;
        }
    }

    // If the cache is invalid or empty.
    if count == 0 {
        let v = &mut simplex[0];
        v.index_a = 0;
        v.index_b = 0;
        v.wa = points_a[0];
        v.wb = point_b_in_a(0);
        v.w = v.wb - v.wa;
        v.a = 0.0;
        count = 1;
    }

    let mut backup = simplex;
    let mut backup_count = 0usize;

    // Keep track of squared distance.
    let mut distance_sq = Real::MAX;
    let mut normal = Vector::ZERO;

    // Run GJK.
    let mut iteration = 0;
    while iteration < MAX_GJK_ITERATIONS {
        // Solve simplex.
        let solved = match count {
            1 => {
                simplex[0].a = 1.0;
                true
            }
            2 => solve_simplex2_3d(&mut simplex, &mut count),
            3 => solve_simplex3_3d(&mut simplex, &mut count),
            4 => solve_simplex4_3d(&mut simplex, &mut count),
            _ => unreachable!(),
        };

        if !solved {
            // No progress - reconstruct last simplex.
            if backup_count == 0 {
                break;
            }
            simplex = backup;
            count = backup_count;
            break;
        }

        if count == 4 {
            // Overlap
            let (pa, pb) = witness_points3(&simplex, count);
            output.point_a = pa;
            output.point_b = pb;
            output.iterations = iteration;
            return output;
        }

        // Assure distance progression.
        let old_distance_sq = distance_sq;

        // Compute closest point.
        let closest_point = match count {
            1 => simplex[0].w,
            2 => simplex[0].a * simplex[0].w + simplex[1].a * simplex[1].w,
            3 => {
                simplex[0].a * simplex[0].w
                    + simplex[1].a * simplex[1].w
                    + simplex[2].a * simplex[2].w
            }
            _ => unreachable!(),
        };

        distance_sq = closest_point.dot(closest_point);

        if distance_sq >= old_distance_sq {
            // No progress - reconstruct last simplex.
            if backup_count == 0 {
                break;
            }
            simplex = backup;
            count = backup_count;
            break;
        }

        // Build new tentative support point.
        let search_direction = match count {
            1 => -simplex[0].w,
            2 => {
                // v = (AB x AO) x AB
                let a = simplex[0].w;
                let b = simplex[1].w;
                let ab = b - a;
                ab.cross(-a).cross(ab)
            }
            3 => {
                // v = AB x AC or v = AC x AB
                let a = simplex[0].w;
                let b = simplex[1].w;
                let c = simplex[2].w;
                let n = (b - a).cross(c - a);
                if n.dot(a) < 0.0 {
                    n
                } else {
                    -n
                }
            }
            _ => unreachable!(),
        };

        if search_direction.length_squared() < 1000.0 * Real::MIN_POSITIVE {
            // The origin is probably contained by a line segment or triangle.
            // Thus the shapes are overlapped.
            let (pa, pb) = witness_points3(&simplex, count);
            output.point_a = pa;
            output.point_b = pb;
            output.iterations = iteration;
            return output;
        }

        normal = -search_direction;

        // Get new support points.
        let index_a = proxy_a.support(-search_direction);
        let support_a = points_a[index_a as usize];
        let index_b = proxy_b.support(pos12.rotation.inverse() * search_direction);
        let support_b = point_b_in_a(index_b);

        // Save current simplex and add new vertex - this can fail if we detect cycling.
        backup = simplex;
        backup_count = count;

        // Check for duplicate support points. This is the main termination criteria.
        let duplicate =
            (0..count).any(|i| simplex[i].index_a == index_a && simplex[i].index_b == index_b);
        if duplicate {
            break;
        }

        simplex[count].index_a = index_a;
        simplex[count].index_b = index_b;
        simplex[count].wa = support_a;
        simplex[count].wb = support_b;
        simplex[count].w = support_b - support_a;
        count += 1;

        iteration += 1;
    }

    let normal = normal.normalize_or_zero();
    if normal == Vector::ZERO {
        // Treat as overlap.
        output.iterations = iteration;
        return output;
    }

    // Build witness points and save cache.
    let (pa, pb) = witness_points3(&simplex, count);
    cache.metric = simplex_metric(&simplex, count);
    cache.count = count.min(3) as u8;
    for i in 0..count.min(3) {
        cache.index_a[i] = simplex[i].index_a;
        cache.index_b[i] = simplex[i].index_b;
    }

    // Results stay in frame A.
    output.point_a = pa;
    output.point_b = pb;
    output.distance = pa.distance(pb);
    output.normal = normal;
    output.iterations = iteration;

    // Apply radii if requested.
    if use_radii {
        let ra = proxy_a.radius;
        let rb = proxy_b.radius;
        output.distance = (output.distance - ra - rb).max(0.0);

        // Keep closest points on perimeter even if overlapped, this way the points move
        // smoothly.
        output.point_a += ra * normal;
        output.point_b -= rb * normal;
    }

    output
}
