use crate::math::Vector;

#[cfg(feature = "alloc")]
use crate::{math::Vector3, query::details::NormalConstraints};

// NOTE: ideally, the normal cone should take into account the point where the normal cone is
//       considered. But as long as we assume that the triangles are one-way we can get away with
//       just relying on the normal directions.
//       Taking the point into account would be technically doable (and desirable if we wanted
//       to define, e.g., a one-way mesh) but requires:
//       1. To make sure the edge pseudo-normals are given in the correct edge order.
//       2. To have access to the contact feature.
//       We can have access to both during the narrow-phase, but leave that as a future
//       potential improvements.
// NOTE: this isn’t equal to the "true" normal cones since concave edges will have pseudo-normals
//       still pointing outward (instead of inward or being empty).
/// The pseudo-normals of a triangle providing approximations of its feature’s normal cones.
#[derive(Clone, Debug)]
pub struct TrianglePseudoNormals {
    /// The triangle’s face normal.
    pub face: Vector,
    // TODO: if we switch this to providing pseudo-normals in a specific order
    //       (e.g. in the same order as the triangle’s edges), then we should
    //       think of fixing that order in the heightfield
    //       triangle_pseudo_normals code.
    /// The edges pseudo-normals, in no particular order.
    pub edges: [Vector; 3],
    /// If `true`, a direction pointing to the back of the triangle is projected into the
    /// mirrored cone (`-face`, `-edges`) and kept, instead of being discarded.
    pub two_sided: bool,
}

#[cfg(feature = "alloc")]
impl NormalConstraints for TrianglePseudoNormals {
    /// Projects the given direction to it is contained in the polygonal
    /// cone defined `self`.
    fn project_local_normal_mut(&self, dir: &mut Vector) -> bool {
        if self.two_sided && dir.dot(self.face) < 0.0 {
            // The back cone is the front cone mirrored through the triangle's plane.
            let mut mirrored = -*dir;
            let _ = self.project_front(&mut mirrored);
            *dir = -mirrored;
            return true;
        }

        self.project_front(dir)
    }
}

#[cfg(feature = "alloc")]
impl TrianglePseudoNormals {
    fn project_front(&self, dir: &mut Vector) -> bool {
        // Find the closest pseudo-normal.
        let dots = Vector3::new(
            dir.dot(self.edges[0]),
            dir.dot(self.edges[1]),
            dir.dot(self.edges[2]),
        );
        let closest_edge = self.edges[dots.max_position()];
        crate::shape::pseudo_normals::project_into_cone(self.face, closest_edge, dir)
    }
}

#[cfg(test)]
#[cfg(all(feature = "dim3", feature = "alloc"))]
mod test {
    use super::NormalConstraints;
    use crate::math::{Real, Vector};
    use crate::shape::TrianglePseudoNormals;

    fn bisector(v1: Vector, v2: Vector) -> Vector {
        (v1 + v2).normalize()
    }

    fn bisector_y(v: Vector) -> Vector {
        bisector(v, Vector::Y)
    }

    #[test]
    fn trivial_pseudo_normals_projection() {
        let pn = TrianglePseudoNormals {
            face: Vector::Y,
            edges: [Vector::Y; 3],
            two_sided: false,
        };

        assert_eq!(
            pn.project_local_normal(Vector::new(1.0, 1.0, 1.0)),
            Some(Vector::Y)
        );
        assert!(pn.project_local_normal(-Vector::Y).is_none());
    }

    #[test]
    fn edge_pseudo_normals_projection_strictly_positive() {
        let bisector = |v1: Vector, v2: Vector| (v1 + v2).normalize();
        let bisector_y = |v: Vector| bisector(v, Vector::Y);

        // The normal cones for this test will be fully contained in the +Y half-space.
        let cones_ref_dir = [
            -Vector::Z,
            -Vector::X,
            Vector::new(1.0, 0.0, 1.0).normalize(),
        ];
        let cones_ends = cones_ref_dir.map(bisector_y);
        let cones_axes = cones_ends.map(bisector_y);

        let pn = TrianglePseudoNormals {
            face: Vector::Y,
            edges: cones_axes.map(|v| v.normalize()),
            two_sided: false,
        };

        for i in 0..3 {
            assert!(pn
                .project_local_normal(cones_ends[i])
                .unwrap()
                .abs_diff_eq(cones_ends[i], 1.0e-5));
            assert_eq!(pn.project_local_normal(cones_axes[i]), Some(cones_axes[i]));

            // Guaranteed to be inside the normal cone of edge i.
            let subdivs = 100;

            for k in 1..100 {
                let v = Vector::Y
                    .lerp(cones_ends[i], k as Real / (subdivs as Real))
                    .normalize();
                assert_eq!(pn.project_local_normal(v).unwrap(), v);
            }

            // Guaranteed to be outside the normal cone of edge i.
            for k in 1..subdivs {
                let v = cones_ref_dir[i]
                    .lerp(cones_ends[i], k as Real / (subdivs as Real))
                    .normalize();
                assert!(pn
                    .project_local_normal(v)
                    .unwrap()
                    .abs_diff_eq(cones_ends[i], 1.0e-5));
            }

            // Guaranteed to be outside the normal cone, and in the -Y half-space.
            for k in 1..subdivs {
                let v = cones_ref_dir[i]
                    .lerp(-Vector::Y, k as Real / (subdivs as Real))
                    .normalize();
                assert!(pn.project_local_normal(v).is_none(),);
            }
        }
    }

    #[test]
    fn edge_pseudo_normals_projection_negative() {
        // The normal cones for this test will be fully contained in the +Y half-space.
        let cones_ref_dir = [
            -Vector::Z,
            -Vector::X,
            Vector::new(1.0, 0.0, 1.0).normalize(),
        ];
        let cones_ends = cones_ref_dir.map(|v| bisector(v, -Vector::Y));
        let cones_axes = [
            bisector(bisector_y(cones_ref_dir[0]), cones_ref_dir[0]),
            bisector(bisector_y(cones_ref_dir[1]), cones_ref_dir[1]),
            bisector(bisector_y(cones_ref_dir[2]), cones_ref_dir[2]),
        ];

        let pn = TrianglePseudoNormals {
            face: Vector::Y,
            edges: cones_axes.map(|v| v.normalize()),
            two_sided: false,
        };

        for i in 0..3 {
            assert_eq!(pn.project_local_normal(cones_axes[i]), Some(cones_axes[i]));

            // Guaranteed to be inside the normal cone of edge i.
            let subdivs = 100;

            for k in 1..subdivs {
                let v = Vector::Y
                    .lerp(cones_ends[i], k as Real / (subdivs as Real))
                    .normalize();
                assert_eq!(pn.project_local_normal(v).unwrap(), v);
            }

            // Guaranteed to be outside the normal cone of edge i.
            // Since it is additionally guaranteed to be in the -Y half-space, we should get None.
            for k in 1..subdivs {
                let v = (-Vector::Y)
                    .lerp(cones_ends[i], k as Real / (subdivs as Real))
                    .normalize();
                assert!(pn.project_local_normal(v).is_none());
            }
        }
    }

    #[test]
    fn two_sided_pseudo_normals_mirror_the_cone() {
        let one_sided = TrianglePseudoNormals {
            face: Vector::Y,
            edges: [Vector::Y; 3],
            two_sided: false,
        };
        let two_sided = TrianglePseudoNormals {
            two_sided: true,
            ..one_sided.clone()
        };

        // Front directions behave the same.
        let front = Vector::new(1.0, 1.0, 1.0);
        assert_eq!(two_sided.project_local_normal(front), Some(Vector::Y));

        // Back directions are kept and constrained to the mirrored cone.
        let back = Vector::new(1.0, -1.0, 1.0);
        assert!(one_sided.project_local_normal(back).is_none());
        assert_eq!(two_sided.project_local_normal(back), Some(-Vector::Y));

        // The mirrored projection is the negation of the front projection.
        let cone_axes = [
            bisector_y(-Vector::Z),
            bisector_y(-Vector::X),
            bisector_y(Vector::new(1.0, 0.0, 1.0).normalize()),
        ];
        let pn = TrianglePseudoNormals {
            face: Vector::Y,
            edges: cone_axes,
            two_sided: true,
        };
        for dir in [
            Vector::new(0.3, 1.0, -0.8).normalize(),
            Vector::new(-1.0, 0.2, 0.1).normalize(),
            Vector::new(0.9, 0.5, 0.9).normalize(),
        ] {
            let front = pn.project_local_normal(dir).unwrap();
            let back = pn.project_local_normal(-dir).unwrap();
            assert!(back.abs_diff_eq(-front, 1.0e-6));
        }
    }
}
