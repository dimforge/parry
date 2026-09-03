use crate::math::Vector;

#[cfg(feature = "alloc")]
use crate::query::details::NormalConstraints;

/// The pseudo-normals of a polyline segment, approximating the outward normal cones of its
/// features.
///
/// This is the 2D segment analog of [`TrianglePseudoNormals`](crate::shape::TrianglePseudoNormals):
/// `face` is the segment's outward normal and `edges` are the outward pseudo-normals at its two
/// endpoints. An oriented polyline uses them to clamp contact normals to one side, so it acts as a
/// one-sided surface.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct SegmentPseudoNormals {
    /// The segment's outward normal.
    pub face: Vector,
    /// The outward pseudo-normals at the segment's two endpoints.
    pub edges: [Vector; 2],
}

#[cfg(feature = "alloc")]
impl NormalConstraints for SegmentPseudoNormals {
    /// Projects `dir` so it lies within the outward cone defined by `self`.
    fn project_local_normal_mut(&self, dir: &mut Vector) -> bool {
        // Find the closest pseudo-normal.
        let closest_edge = if dir.dot(self.edges[0]) >= dir.dot(self.edges[1]) {
            self.edges[0]
        } else {
            self.edges[1]
        };
        crate::shape::pseudo_normals::project_into_cone(self.face, closest_edge, dir)
    }
}

#[cfg(test)]
#[cfg(all(feature = "dim2", feature = "alloc"))]
mod test {
    use super::{NormalConstraints, SegmentPseudoNormals};
    use crate::math::Vector;

    fn bisector(v1: Vector, v2: Vector) -> Vector {
        (v1 + v2).normalize()
    }

    #[test]
    fn degenerate_cone_collapses_to_face() {
        let pn = SegmentPseudoNormals {
            face: Vector::Y,
            edges: [Vector::Y, Vector::Y],
        };

        assert_eq!(
            pn.project_local_normal(Vector::new(1.0, 1.0)),
            Some(Vector::Y)
        );
        assert!(pn.project_local_normal(-Vector::Y).is_none());
    }

    #[test]
    fn clamps_into_the_outward_cone() {
        // A cone centered on +Y, bounded at +-45 degrees (endpoint pseudo-normals at +-22.5 degrees).
        let ends = [
            bisector(Vector::Y, Vector::X),
            bisector(Vector::Y, -Vector::X),
        ];
        let edges = [bisector(Vector::Y, ends[0]), bisector(Vector::Y, ends[1])];
        let pn = SegmentPseudoNormals {
            face: Vector::Y,
            edges,
        };

        // Inside the cone: returned unchanged.
        assert_eq!(pn.project_local_normal(Vector::Y), Some(Vector::Y));
        for edge in edges {
            assert_eq!(pn.project_local_normal(edge), Some(edge));
        }
        let inside = Vector::new(0.2, 1.0).normalize();
        assert!(pn
            .project_local_normal(inside)
            .unwrap()
            .abs_diff_eq(inside, 1.0e-5));

        // Outside the cone but still outward: pulled onto the boundary, not left crossing the surface.
        let sideways = Vector::new(1.0, 0.2).normalize();
        let clamped = pn.project_local_normal(sideways).unwrap();
        assert!(clamped.dot(Vector::Y) > sideways.dot(Vector::Y));
        assert!(clamped.abs_diff_eq(ends[0], 1.0e-5));

        // Into the solid (-face half-space): rejected.
        assert!(pn.project_local_normal(-Vector::Y).is_none());
        assert!(pn
            .project_local_normal(Vector::new(0.3, -1.0).normalize())
            .is_none());
    }
}
