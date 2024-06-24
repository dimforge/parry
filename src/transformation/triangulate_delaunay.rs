use crate::math::Real;
use crate::utils::point_in_poly2d;
use na::Point2;
use spade::{ConstrainedDelaunayTriangulation, Point2 as Pt2, Triangulation};

pub fn triangulate_delaunay(
    poly1: &[Point2<Real>],
    holes: &[Vec<Point2<Real>>],
) -> (Vec<Point2<Real>>, Vec<[u32; 3]>) {
    let mut cdt = ConstrainedDelaunayTriangulation::<Pt2<_>>::new();
    let mut handles = vec![];

    for poly in std::iter::once(poly1).chain(holes.iter().map(|h| &h[..])) {
        handles.clear();
        handles.extend(poly.iter().map(|pt| cdt.insert(Pt2::new(pt.x, pt.y))));

        for ia in 0..poly.len() {
            let ib = (ia + 1) % poly.len();

            if let (Ok(handle_a), Ok(handle_b)) = (handles[ia], handles[ib]) {
                let _ = cdt.add_constraint_and_split(handle_a, handle_b, |v| v);
            }
        }
    }

    // Eliminate unwanted triangles.
    let mut result_idx = vec![];
    let mut result_pts = vec![];
    let mut handle_to_pt_id = vec![None; cdt.num_vertices()];

    for face in cdt.inner_faces() {
        let tri_handles = face.vertices();

        let tri_idx = tri_handles.map(|v| {
            if let Some(id) = handle_to_pt_id[v.fix().index()] {
                id
            } else {
                let pt = v.data();
                let id = result_pts.len() as u32;
                result_pts.push(Point2::new(pt.x, pt.y));
                handle_to_pt_id[v.fix().index()] = Some(id);
                id
            }
        });
        let tri_pts = tri_idx.map(|i| result_pts[i as usize]);
        let tri_center = crate::utils::center(&tri_pts);

        if point_in_poly2d(&tri_center, poly1)
            && holes.iter().all(|hole| !point_in_poly2d(&tri_center, hole))
        {
            // Keep the triangle, its center is inside the polygon, but outside all the holes.
            result_idx.push(tri_idx);
        }
    }

    (result_pts, result_idx)
}
