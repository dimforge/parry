use num::{Bounded, Zero};
use std::cmp::Ordering;

use crate::math::Real;
use crate::shape::Triangle;
use crate::transformation::{
    self,
    convex_hull_utils::{indexed_support_point_id, normalize, support_point_id},
};
use crate::utils;
use na::{self, Matrix3, Point2, Point3, Vector3};

/// Computes the convex hull of a set of 3d points.
pub fn convex_hull3(points: &[Point3<Real>]) -> (Vec<Point3<Real>>, Vec<Point3<u32>>) {
    if points.is_empty() {
        return (Vec::new(), Vec::new());
    }

    let mut normalized_points = points.to_vec();
    let _ = normalize(&mut normalized_points[..]);

    let mut undecidable_points = Vec::new();
    let mut horizon_loop_facets = Vec::new();
    let mut horizon_loop_ids = Vec::new();
    let mut removed_facets = Vec::new();

    let mut triangles;

    match get_initial_mesh(points, &mut normalized_points[..], &mut undecidable_points) {
        InitialMesh::Facets(facets) => {
            triangles = facets;
        }
        InitialMesh::ResultMesh(vertices, indices) => {
            return (vertices, indices);
        }
    }

    let mut i = 0;
    while i != triangles.len() {
        horizon_loop_facets.clear();
        horizon_loop_ids.clear();

        if !triangles[i].valid {
            i = i + 1;
            continue;
        }

        // FIXME: use triangles[i].furthest_point instead.
        let pt_id = indexed_support_point_id(
            &triangles[i].normal,
            &normalized_points[..],
            &triangles[i].visible_points[..],
        );

        if let Some(point) = pt_id {
            removed_facets.clear();

            triangles[i].valid = false;
            removed_facets.push(i);

            for j in 0usize..3 {
                compute_silhouette(
                    triangles[i].adj[j],
                    triangles[i].indirect_adj_id[j],
                    point,
                    &mut horizon_loop_facets,
                    &mut horizon_loop_ids,
                    &normalized_points[..],
                    &mut removed_facets,
                    &mut triangles[..],
                );
            }

            if horizon_loop_facets.is_empty() {
                // Due to inaccuracies, the silhouette could not be computed
                // (the point seems to be visible from… every triangle).
                let mut any_valid = false;
                for j in i + 1..triangles.len() {
                    if triangles[j].valid {
                        any_valid = true;
                    }
                }

                if any_valid {
                    //                        println!("Warning: exitting an unfinished work.");
                }

                // FIXME: this is verry harsh.
                triangles[i].valid = true;
                break;
            }

            attach_and_push_facets3(
                &horizon_loop_facets[..],
                &horizon_loop_ids[..],
                point,
                &normalized_points[..],
                &mut triangles,
                &removed_facets[..],
                &mut undecidable_points,
            );
        }

        i = i + 1;
    }

    let mut idx = Vec::new();

    for facet in triangles.iter() {
        if facet.valid {
            idx.push(Point3::new(
                facet.pts[0] as u32,
                facet.pts[1] as u32,
                facet.pts[2] as u32,
            ));
        }
    }

    let mut points = points.to_vec();
    utils::remove_unused_points(&mut points, &mut idx[..]);

    assert!(points.len() != 0, "Internal error: empty output mesh.");

    (points, idx)
}

enum InitialMesh {
    Facets(Vec<TriangleFacet>),
    ResultMesh(Vec<Point3<Real>>, Vec<Point3<u32>>),
}

fn build_degenerate_mesh_point(point: Point3<Real>) -> (Vec<Point3<Real>>, Vec<Point3<u32>>) {
    let ta = Point3::new(0u32, 0, 0);
    let tb = Point3::new(0u32, 0, 0);

    (vec![point], vec![ta, tb])
}

fn build_degenerate_mesh_segment(
    dir: &Vector3<Real>,
    points: &[Point3<Real>],
) -> (Vec<Point3<Real>>, Vec<Point3<u32>>) {
    let a = utils::point_cloud_support_point(dir, points);
    let b = utils::point_cloud_support_point(&-*dir, points);

    let ta = Point3::new(0u32, 1, 0);
    let tb = Point3::new(1u32, 0, 0);

    (vec![a, b], vec![ta, tb])
}

fn get_initial_mesh(
    original_points: &[Point3<Real>],
    normalized_points: &mut [Point3<Real>],
    undecidable: &mut Vec<usize>,
) -> InitialMesh {
    /*
     * Compute the eigenvectors to see if the input data live on a subspace.
     */
    let cov_mat;
    let eigvec;
    let eigval;

    #[cfg(not(feature = "improved_fixed_point_support"))]
    {
        cov_mat = crate::utils::cov(normalized_points);
        let eig = cov_mat.symmetric_eigen();
        eigvec = eig.eigenvectors;
        eigval = eig.eigenvalues;
    }

    #[cfg(feature = "improved_fixed_point_support")]
    {
        cov_mat = Matrix3::identity();
        eigvec = Matrix3::identity();
        eigval = Vector3::repeat(1.0);
    }

    let mut eigpairs = [
        (eigvec.column(0).into_owned(), eigval[0]),
        (eigvec.column(1).into_owned(), eigval[1]),
        (eigvec.column(2).into_owned(), eigval[2]),
    ];

    /*
     * Sort in decreasing order wrt. eigenvalues.
     */
    eigpairs.sort_by(|a, b| {
        if a.1 > b.1 {
            Ordering::Less // `Less` and `Greater` are reversed.
        } else if a.1 < b.1 {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    });

    /*
     * Count the dimension the data lives in.
     */
    let mut dimension = 0;
    while dimension < 3 {
        if relative_eq!(eigpairs[dimension].1, 0.0, epsilon = 1.0e-7) {
            break;
        }

        dimension = dimension + 1;
    }

    match dimension {
        0 => {
            // The hull is a point.
            let (vtx, idx) = build_degenerate_mesh_point(original_points[0].clone());
            InitialMesh::ResultMesh(vtx, idx)
        }
        1 => {
            // The hull is a segment.
            let (vtx, idx) = build_degenerate_mesh_segment(&eigpairs[0].0, original_points);
            InitialMesh::ResultMesh(vtx, idx)
        }
        2 => {
            // The hull is a triangle.
            // Project into the principal halfspace…
            let axis1 = &eigpairs[0].0;
            let axis2 = &eigpairs[1].0;

            let mut subspace_points = Vec::with_capacity(normalized_points.len());

            for point in normalized_points.iter() {
                subspace_points.push(Point2::new(
                    point.coords.dot(axis1),
                    point.coords.dot(axis2),
                ))
            }

            // … and compute the 2d convex hull.
            let idx = transformation::convex_hull2_idx(&subspace_points[..]);

            // Finalize the result, triangulating the polyline.
            let npoints = idx.len();
            let coords = idx
                .into_iter()
                .map(|i| original_points[i].clone())
                .collect();
            let mut triangles = Vec::with_capacity(npoints + npoints - 4);

            let a = 0u32;

            for id in 1u32..npoints as u32 - 1 {
                triangles.push(Point3::new(a, id, id + 1));
                triangles.push(Point3::new(id, a, id + 1));
            }

            InitialMesh::ResultMesh(coords, triangles)
        }
        3 => {
            // The hull is a polyhedron.
            // Find a initial triangle lying on the principal halfspace…
            let diag = eigval.map(|e| 1.0 / e);
            let diag = Matrix3::from_diagonal(&diag);
            let icov = eigvec * diag * eigvec.transpose();

            for point in normalized_points.iter_mut() {
                *point = Point3::origin() + icov * point.coords;
            }

            let p1 = support_point_id(&eigpairs[0].0, normalized_points).unwrap();
            let p2 = support_point_id(&-eigpairs[0].0, normalized_points).unwrap();

            let mut max_area = 0.0;
            let mut p3 = usize::max_value();

            for (i, point) in normalized_points.iter().enumerate() {
                let area =
                    Triangle::new(normalized_points[p1], normalized_points[p2], *point).area();

                if area > max_area {
                    max_area = area;
                    p3 = i;
                }
            }

            assert!(
                p3 != usize::max_value(),
                "Internal convex hull error: no triangle found."
            );

            // Build two facets with opposite normals
            let mut f1 = TriangleFacet::new(p1, p2, p3, normalized_points);
            let mut f2 = TriangleFacet::new(p2, p1, p3, normalized_points);

            // Link the facets together
            f1.set_facets_adjascency(1, 1, 1, 0, 2, 1);
            f2.set_facets_adjascency(0, 0, 0, 0, 2, 1);

            let mut facets = vec![f1, f2];

            // … and attribute visible points to each one of them.
            // FIXME: refactor this with the two others.
            let mut ignored = 0usize;
            for point in 0..normalized_points.len() {
                if point == p1 || point == p2 || point == p3 {
                    continue;
                }

                let mut furthest = usize::max_value();
                let mut furthest_dist = 0.0;

                for (i, curr_facet) in facets.iter().enumerate() {
                    if curr_facet.can_be_seen_by(point, normalized_points) {
                        let distance = curr_facet.distance_to_point(point, normalized_points);

                        if distance > furthest_dist {
                            furthest = i;
                            furthest_dist = distance;
                        }
                    }
                }

                if furthest != usize::max_value() {
                    facets[furthest].add_visible_point(point, normalized_points);
                } else {
                    undecidable.push(point);
                    ignored = ignored + 1;
                }

                // If none of the facet can be seen from the point, it is naturally deleted.
            }

            verify_facet_links(0, &facets[..]);
            verify_facet_links(1, &facets[..]);

            InitialMesh::Facets(facets)
        }
        _ => unreachable!(),
    }
}

fn compute_silhouette(
    facet: usize,
    indirect_id: usize,
    point: usize,
    out_facets: &mut Vec<usize>,
    out_adj_idx: &mut Vec<usize>,
    points: &[Point3<Real>],
    removed_facets: &mut Vec<usize>,
    triangles: &mut [TriangleFacet],
) {
    if triangles[facet].valid {
        if !triangles[facet].can_be_seen_by_or_is_affinely_dependent_with_contour(
            point,
            points,
            indirect_id,
        ) {
            out_facets.push(facet);
            out_adj_idx.push(indirect_id);
        } else {
            triangles[facet].valid = false; // The facet must be removed from the convex hull.
            removed_facets.push(facet);

            compute_silhouette(
                triangles[facet].adj[(indirect_id + 1) % 3],
                triangles[facet].indirect_adj_id[(indirect_id + 1) % 3],
                point,
                out_facets,
                out_adj_idx,
                points,
                removed_facets,
                triangles,
            );
            compute_silhouette(
                triangles[facet].adj[(indirect_id + 2) % 3],
                triangles[facet].indirect_adj_id[(indirect_id + 2) % 3],
                point,
                out_facets,
                out_adj_idx,
                points,
                removed_facets,
                triangles,
            );
        }
    }
}

fn verify_facet_links(ifacet: usize, facets: &[TriangleFacet]) {
    let facet = &facets[ifacet];

    for i in 0usize..3 {
        let adji = &facets[facet.adj[i]];

        assert!(
            adji.adj[facet.indirect_adj_id[i]] == ifacet
                && adji.first_point_from_edge(facet.indirect_adj_id[i])
                    == facet.second_point_from_edge(adji.indirect_adj_id[facet.indirect_adj_id[i]])
                && adji.second_point_from_edge(facet.indirect_adj_id[i])
                    == facet.first_point_from_edge(adji.indirect_adj_id[facet.indirect_adj_id[i]])
        )
    }
}

fn attach_and_push_facets3(
    horizon_loop_facets: &[usize],
    horizon_loop_ids: &[usize],
    point: usize,
    points: &[Point3<Real>],
    triangles: &mut Vec<TriangleFacet>,
    removed_facets: &[usize],
    undecidable: &mut Vec<usize>,
) {
    // The horizon is built to be in CCW order.
    let mut new_facets = Vec::with_capacity(horizon_loop_facets.len());

    // Create new facets.
    let mut adj_facet: usize;
    let mut indirect_id: usize;

    for i in 0..horizon_loop_facets.len() {
        adj_facet = horizon_loop_facets[i];
        indirect_id = horizon_loop_ids[i];

        let facet = TriangleFacet::new(
            point,
            triangles[adj_facet].second_point_from_edge(indirect_id),
            triangles[adj_facet].first_point_from_edge(indirect_id),
            points,
        );
        new_facets.push(facet);
    }

    // Link the facets together.
    for i in 0..horizon_loop_facets.len() {
        let prev_facet;

        if i == 0 {
            prev_facet = triangles.len() + horizon_loop_facets.len() - 1;
        } else {
            prev_facet = triangles.len() + i - 1;
        }

        let middle_facet = horizon_loop_facets[i];
        let next_facet = triangles.len() + (i + 1) % horizon_loop_facets.len();
        let middle_id = horizon_loop_ids[i];

        new_facets[i].set_facets_adjascency(prev_facet, middle_facet, next_facet, 2, middle_id, 0);
        triangles[middle_facet].adj[middle_id] = triangles.len() + i; // The future id of curr_facet.
        triangles[middle_facet].indirect_adj_id[middle_id] = 1;
    }

    // Assign to each facets some of the points which can see it.
    // FIXME: refactor this with the others.
    for curr_facet in removed_facets.iter() {
        for visible_point in triangles[*curr_facet].visible_points.iter() {
            if *visible_point == point {
                continue;
            }

            let mut furthest = usize::max_value();
            let mut furthest_dist = 0.0;

            for (i, curr_facet) in new_facets.iter_mut().enumerate() {
                if curr_facet.can_be_seen_by(*visible_point, points) {
                    let distance = curr_facet.distance_to_point(*visible_point, points);

                    if distance > furthest_dist {
                        furthest = i;
                        furthest_dist = distance;
                    }
                }
            }

            if furthest != usize::max_value() {
                new_facets[furthest].add_visible_point(*visible_point, points);
            }

            // If none of the facet can be seen from the point, it is naturally deleted.
        }
    }

    // Try to assign collinear points to one of the new facets.
    let mut i = 0;

    while i != undecidable.len() {
        let mut furthest = usize::max_value();
        let mut furthest_dist = 0.0;
        let undecidable_point = undecidable[i];

        for (j, curr_facet) in new_facets.iter_mut().enumerate() {
            if curr_facet.can_be_seen_by(undecidable_point, points) {
                let distance = curr_facet.distance_to_point(undecidable_point, points);

                if distance > furthest_dist {
                    furthest = j;
                    furthest_dist = distance;
                }
            }
        }

        if furthest != usize::max_value() {
            new_facets[furthest].add_visible_point(undecidable_point, points);
            let _ = undecidable.swap_remove(i);
        } else {
            i = i + 1;
        }
    }

    // Push facets.
    // FIXME: can we avoid the tmp vector `new_facets` ?
    for curr_facet in new_facets.into_iter() {
        triangles.push(curr_facet);
    }
}

struct TriangleFacet {
    valid: bool,
    normal: Vector3<Real>,
    adj: [usize; 3],
    indirect_adj_id: [usize; 3],
    pts: [usize; 3],
    visible_points: Vec<usize>,
    furthest_point: usize,
    furthest_distance: Real,
}

impl TriangleFacet {
    pub fn new(p1: usize, p2: usize, p3: usize, points: &[Point3<Real>]) -> TriangleFacet {
        let p1p2 = points[p2] - points[p1];
        let p1p3 = points[p3] - points[p1];

        let mut normal = p1p2.cross(&p1p3);
        if normal.normalize_mut().is_zero() {
            panic!("ConvexHull hull failure: a facet must not be affinely dependent.");
        }

        TriangleFacet {
            valid: true,
            normal,
            adj: [0, 0, 0],
            indirect_adj_id: [0, 0, 0],
            pts: [p1, p2, p3],
            visible_points: Vec::new(),
            furthest_point: Bounded::max_value(),
            furthest_distance: 0.0,
        }
    }

    pub fn add_visible_point(&mut self, pid: usize, points: &[Point3<Real>]) {
        let distance = self.distance_to_point(pid, points);

        if distance > self.furthest_distance {
            self.furthest_distance = distance;
            self.furthest_point = pid;
        }

        self.visible_points.push(pid);
    }

    pub fn distance_to_point(&self, point: usize, points: &[Point3<Real>]) -> Real {
        self.normal.dot(&(points[point] - points[self.pts[0]]))
    }

    pub fn set_facets_adjascency(
        &mut self,
        adj1: usize,
        adj2: usize,
        adj3: usize,
        id_adj1: usize,
        id_adj2: usize,
        id_adj3: usize,
    ) {
        self.indirect_adj_id[0] = id_adj1;
        self.indirect_adj_id[1] = id_adj2;
        self.indirect_adj_id[2] = id_adj3;

        self.adj[0] = adj1;
        self.adj[1] = adj2;
        self.adj[2] = adj3;
    }

    pub fn first_point_from_edge(&self, id: usize) -> usize {
        self.pts[id]
    }

    pub fn second_point_from_edge(&self, id: usize) -> usize {
        self.pts[(id + 1) % 3]
    }

    pub fn can_be_seen_by(&self, point: usize, points: &[Point3<Real>]) -> bool {
        let p0 = points[self.pts[0]];
        let p1 = points[self.pts[1]];
        let p2 = points[self.pts[2]];
        let pt = points[point];

        let _eps = crate::math::DEFAULT_EPSILON;

        (pt - p0).dot(&self.normal) > _eps * 100.0
            && !Triangle::new(p0, p1, pt).is_affinely_dependent()
            && !Triangle::new(p0, p2, pt).is_affinely_dependent()
            && !Triangle::new(p1, p2, pt).is_affinely_dependent()
    }

    pub fn can_be_seen_by_or_is_affinely_dependent_with_contour(
        &self,
        point: usize,
        points: &[Point3<Real>],
        edge: usize,
    ) -> bool {
        let p0 = points[self.first_point_from_edge(edge)];
        let p1 = points[self.second_point_from_edge(edge)];
        let pt = points[point];

        let aff_dep = Triangle::new(p0, p1, pt).is_affinely_dependent()
            || Triangle::new(p0, pt, p1).is_affinely_dependent()
            || Triangle::new(p1, p0, pt).is_affinely_dependent()
            || Triangle::new(p1, pt, p0).is_affinely_dependent()
            || Triangle::new(pt, p0, p1).is_affinely_dependent()
            || Triangle::new(pt, p1, p0).is_affinely_dependent();

        (pt - p0).dot(&self.normal) >= 0.0 || aff_dep
    }
}

#[cfg(test)]
mod test {
    use crate::transformation;
    #[cfg(feature = "dim2")]
    use na::Point2;

    #[cfg(feature = "dim2")]
    #[test]
    fn test_simple_convex_hull() {
        let points = [
            Point2::new(4.723881f32, 3.597233),
            Point2::new(3.333363, 3.429991),
            Point2::new(3.137215, 2.812263),
        ];

        let chull = transformation::convex_hull(points.as_slice());

        assert!(chull.coords.len() == 3);
    }

    #[cfg(feature = "dim3")]
    #[test]
    fn test_ball_convex_hull() {
        use crate::shape::Ball;

        // This triggered a failure to an affinely dependent facet.
        let (points, _) = Ball::new(0.4).to_trimesh(20, 20);
        let (vertices, _) = transformation::convex_hull(points.as_slice());

        // dummy test, we are just checking that the construction did not fail.
        assert!(vertices.len() == vertices.len());
    }
}
