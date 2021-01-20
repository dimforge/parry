use super::TriangleFacet;
use crate::math::Real;
use crate::shape::Triangle;
use crate::transformation;
use crate::transformation::convex_hull_utils::support_point_id;
use crate::utils;
use na::{Point2, Point3, Vector3};
use std::cmp::Ordering;

pub enum InitialMesh {
    Facets(Vec<TriangleFacet>),
    ResultMesh(Vec<Point3<Real>>, Vec<[u32; 3]>),
}

fn build_degenerate_mesh_point(point: Point3<Real>) -> (Vec<Point3<Real>>, Vec<[u32; 3]>) {
    let ta = [0u32; 3];
    let tb = [0u32; 3];

    (vec![point], vec![ta, tb])
}

fn build_degenerate_mesh_segment(
    dir: &Vector3<Real>,
    points: &[Point3<Real>],
) -> (Vec<Point3<Real>>, Vec<[u32; 3]>) {
    let a = utils::point_cloud_support_point(dir, points);
    let b = utils::point_cloud_support_point(&-*dir, points);

    let ta = [0u32, 1, 0];
    let tb = [1u32, 0, 0];

    (vec![a, b], vec![ta, tb])
}

pub fn get_initial_mesh(
    original_points: &[Point3<Real>],
    normalized_points: &mut [Point3<Real>],
    undecidable: &mut Vec<usize>,
) -> InitialMesh {
    /*
     * Compute the eigenvectors to see if the input data live on a subspace.
     */
    let eigvec;
    let eigval;

    #[cfg(not(feature = "improved_fixed_point_support"))]
    {
        let cov_mat = crate::utils::cov(normalized_points);
        let eig = cov_mat.symmetric_eigen();
        eigvec = eig.eigenvectors;
        eigval = eig.eigenvalues;
    }

    #[cfg(feature = "improved_fixed_point_support")]
    {
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

            for id in 1u32..npoints as u32 - 1 {
                triangles.push([0, id, id + 1]);
            }

            // NOTE: We use a different starting point for the triangulation
            // of the bottom faces in order to avoid bad topology where
            // and edge would end be being shared by more than two triangles.
            for id in 0u32..npoints as u32 - 2 {
                let a = npoints as u32 - 1;
                triangles.push([a, id + 1, id]);
            }

            InitialMesh::ResultMesh(coords, triangles)
        }
        3 => {
            // The hull is a polyhedron.
            // Find a initial triangle lying on the principal halfspace…
            let center = crate::utils::center(&normalized_points);

            for point in normalized_points.iter_mut() {
                *point = Point3::from((*point - center) / eigval.amax());
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
                if normalized_points[point] == normalized_points[p1]
                    || normalized_points[point] == normalized_points[p2]
                    || normalized_points[point] == normalized_points[p3]
                {
                    continue;
                }

                let mut furthest = usize::max_value();
                let mut furthest_dist = 0.0;

                for (i, curr_facet) in facets.iter().enumerate() {
                    if curr_facet.can_see_point(point, normalized_points) {
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

            super::check_facet_links(0, &facets[..]);
            super::check_facet_links(1, &facets[..]);

            InitialMesh::Facets(facets)
        }
        _ => unreachable!(),
    }
}
