use super::TriangleFacet;
use crate::math::Real;
use na::Point3;

pub fn check_facet_links(ifacet: usize, facets: &[TriangleFacet]) {
    let facet = &facets[ifacet];

    for i in 0..3 {
        assert!(facets[facet.adj[i]].valid);
    }

    for i in 0..3 {
        let adj_facet = &facets[facet.adj[i]];

        assert_eq!(adj_facet.adj[facet.indirect_adj_id[i]], ifacet);
        assert_eq!(adj_facet.indirect_adj_id[facet.indirect_adj_id[i]], i);
        assert_eq!(
            adj_facet.first_point_from_edge(facet.indirect_adj_id[i]),
            facet.second_point_from_edge(i)
        );
        assert_eq!(
            adj_facet.second_point_from_edge(facet.indirect_adj_id[i]),
            facet.first_point_from_edge(i)
        );
    }
}

/// Checks if a convex-hull is properly formed.
pub fn check_convex_hull(points: &[Point3<Real>], triangles: &[[u32; 3]]) {
    use crate::utils::hashmap::{Entry, HashMap};
    use crate::utils::SortedPair;
    let mut edges = HashMap::default();

    struct EdgeData {
        adjacent_triangles: [usize; 2],
    }

    // println!(
    //     "Investigating with {} triangles, and {} points.",
    //     triangles.len(),
    //     points.len()
    // );
    // print_buildable_vec("vertices", points);
    // print_buildable_vec("indices", triangles);

    for i in 0..points.len() {
        for j in i + 1..points.len() {
            if points[i] == points[j] {
                println!("Duplicate: {}", points[i]);
                panic!("Found duplicate points.")
            }
        }
    }

    for (itri, tri) in triangles.iter().enumerate() {
        assert!(tri[0] != tri[1]);
        assert!(tri[0] != tri[2]);
        assert!(tri[2] != tri[1]);

        for i in 0..3 {
            let ivtx1 = tri[i as usize];
            let ivtx2 = tri[(i as usize + 1) % 3];
            let edge_key = SortedPair::new(ivtx1, ivtx2);

            match edges.entry(edge_key) {
                Entry::Vacant(e) => {
                    let _ = e.insert(EdgeData {
                        adjacent_triangles: [itri, usize::MAX],
                    });
                }
                Entry::Occupied(mut e) => {
                    if e.get().adjacent_triangles[1] != usize::MAX {
                        panic!(
                            "Detected t-junction for triangle {}, edge: {:?}.",
                            itri,
                            (ivtx1, ivtx2)
                        );
                    }

                    e.get_mut().adjacent_triangles[1] = itri;
                }
            }
        }
    }

    for edge in &edges {
        if edge.1.adjacent_triangles[1] == usize::MAX {
            panic!("Detected unfinished triangle.");
        }
    }

    // Check Euler characteristic.
    assert_eq!(points.len() + triangles.len() - edges.len(), 2);
}

// fn print_buildable_vec<T: std::fmt::Display + na::Scalar>(desc: &str, elts: &[Point3<T>]) {
//     print!("let {} = vec![", desc);
//     for elt in elts {
//         print!("Point3::new({},{},{}),", elt.x, elt.y, elt.z);
//     }
//     println!("];")
// }
