use crate::math::{Real, Vector};

#[cfg(feature = "alloc")]
use crate::query::details::NormalConstraints;
#[cfg(feature = "alloc")]
use alloc::vec::Vec;

/// Whether the turn from `from` to `to` is clockwise.
#[cfg(all(feature = "alloc", feature = "dim2"))]
fn turns_clockwise(from: Vector, to: Vector) -> bool {
    from.x * to.y - from.y * to.x < 0.0
}

/// The directions one boundary edge of a [`Compound`](crate::shape::Compound) part may push along.
///
/// Each limit reaches halfway towards the neighbouring boundary edge, so the two edges meeting at a
/// corner split its range between them. A limit facing a cut stops at `face`: the surface ends
/// there.
#[cfg(all(feature = "alloc", feature = "dim2"))]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Copy, Debug)]
pub struct CompoundEdgeCone {
    /// The edge's outward normal.
    pub face: Vector,
    /// How far the cone reaches at the edge's clockwise end.
    pub clockwise_limit: Vector,
    /// How far the cone reaches at the edge's counter-clockwise end.
    pub counter_clockwise_limit: Vector,
}

#[cfg(all(feature = "alloc", feature = "dim2"))]
impl CompoundEdgeCone {
    fn contains(&self, dir: Vector) -> bool {
        dir.dot(self.face) > 0.0
            && !turns_clockwise(self.clockwise_limit, dir)
            && !turns_clockwise(dir, self.counter_clockwise_limit)
    }

    /// The direction in this cone closest to `dir`.
    fn clamped(&self, dir: Vector) -> Vector {
        if self.contains(dir) {
            dir
        } else if dir.dot(self.clockwise_limit) >= dir.dot(self.counter_clockwise_limit) {
            self.clockwise_limit
        } else {
            self.counter_clockwise_limit
        }
    }
}

/// The outward normal cones of one [`Compound`](crate::shape::Compound) part.
///
/// Parts of a decomposed outline meet along cut edges that lie inside the union, where nothing can
/// reach them -- but each part is convex and reports contacts against its own faces regardless, so
/// a body sliding across a join catches on a face that is not a surface. Dropping the edges a
/// sibling part covers leaves every part facing the way the union does.
///
/// The [`Compound`](crate::shape::Compound) counterpart of
/// [`SegmentPseudoNormals`](crate::shape::SegmentPseudoNormals) and
/// [`TrianglePseudoNormals`](crate::shape::TrianglePseudoNormals), which describe a part that is a
/// single face; a compound part is a whole polygon, so it carries one cone per boundary edge.
#[cfg(all(feature = "alloc", feature = "dim2"))]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug, Default)]
pub struct CompoundPseudoNormals {
    /// One cone per boundary edge, in the part's local frame. Cut edges are absent.
    pub boundary_edges: Vec<CompoundEdgeCone>,
}

#[cfg(all(feature = "alloc", feature = "dim2"))]
impl NormalConstraints for CompoundPseudoNormals {
    /// Projects `dir` onto the closest direction the part's outline can face.
    ///
    /// Returns `false` when nothing faces that way, so a contact reached from inside the union is
    /// dropped rather than turned around.
    fn project_local_normal_mut(&self, dir: &mut Vector) -> bool {
        let mut nearest: Option<(Vector, Real)> = None;

        for edge in &self.boundary_edges {
            let candidate = edge.clamped(*dir);
            let alignment = dir.dot(candidate);

            match nearest {
                Some((_, best)) if alignment <= best => {}
                _ => nearest = Some((candidate, alignment)),
            }
        }

        match nearest {
            Some((direction, alignment)) if alignment > 0.0 => {
                *dir = direction;
                true
            }
            _ => false,
        }
    }
}

/// The directions one boundary face of a [`Compound`](crate::shape::Compound) part may push
/// along.
///
/// Each edge pseudo-normal reaches halfway towards the face across that edge; at a cut it keeps
/// the face's own normal, since the surface ends there.
#[cfg(all(feature = "alloc", feature = "dim3"))]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct CompoundFaceCone {
    /// The face's outward normal.
    pub face: Vector,
    /// The outward pseudo-normal at each edge of the face.
    pub edge_pseudo_normals: Vec<Vector>,
}

/// The outward normal cones of one [`Compound`](crate::shape::Compound) part.
///
/// The 3D counterpart of the 2D edge cones: faces two parts share are cuts made by a convex
/// decomposition, interior to the union where nothing can reach them, so they are absent and no
/// contact normal can be clamped onto one.
#[cfg(all(feature = "alloc", feature = "dim3"))]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug, Default)]
pub struct CompoundPseudoNormals {
    /// One cone per boundary face of the part, in the part's local frame.
    pub boundary_faces: Vec<CompoundFaceCone>,
}

#[cfg(all(feature = "alloc", feature = "dim3"))]
impl NormalConstraints for CompoundPseudoNormals {
    /// Projects `dir` into the normal cone of the boundary face it points most directly out of,
    /// the way [`TrianglePseudoNormals`](crate::shape::TrianglePseudoNormals) does for one
    /// triangle.
    ///
    /// Returns `false` for a part with no boundary face at all: it is buried inside the union,
    /// where nothing can touch it.
    fn project_local_normal_mut(&self, dir: &mut Vector) -> bool {
        let mut best: Option<(&CompoundFaceCone, Real)> = None;

        for cone in &self.boundary_faces {
            let alignment = dir.dot(cone.face);
            match best {
                Some((_, best_alignment)) if alignment <= best_alignment => {}
                _ => best = Some((cone, alignment)),
            }
        }

        let Some((cone, _)) = best else {
            return false;
        };

        let closest_edge = cone
            .edge_pseudo_normals
            .iter()
            .copied()
            .max_by(|a, b| dir.dot(*a).total_cmp(&dir.dot(*b)))
            .unwrap_or(cone.face);
        crate::shape::pseudo_normals::project_into_cone(cone.face, closest_edge, dir)
    }
}

#[cfg(test)]
#[cfg(all(feature = "dim2", feature = "alloc"))]
mod test {
    use super::{CompoundEdgeCone, CompoundPseudoNormals, NormalConstraints};
    use crate::math::Vector;
    use alloc::vec;

    fn cone(face: Vector, clockwise: Vector, counter_clockwise: Vector) -> CompoundEdgeCone {
        CompoundEdgeCone {
            face,
            clockwise_limit: clockwise,
            counter_clockwise_limit: counter_clockwise,
        }
    }

    /// A cone with nowhere to open, as happens where a surface runs into a cut at both ends.
    fn pinned(face: Vector) -> CompoundEdgeCone {
        cone(face, face, face)
    }

    #[test]
    fn buried_part_rejects_everything() {
        let normals = CompoundPseudoNormals::default();

        assert!(normals.project_local_normal(Vector::Y).is_none());
        assert!(normals.project_local_normal(-Vector::Y).is_none());
    }

    #[test]
    fn a_direction_already_on_the_outline_is_left_alone() {
        let normals = CompoundPseudoNormals {
            boundary_edges: vec![cone(
                Vector::Y,
                Vector::new(1.0, 1.0).normalize(),
                Vector::new(-1.0, 1.0).normalize(),
            )],
        };

        assert_eq!(normals.project_local_normal(Vector::Y), Some(Vector::Y));

        let inside_the_cone = Vector::new(0.3, 1.0).normalize();
        assert_eq!(
            normals.project_local_normal(inside_the_cone),
            Some(inside_the_cone)
        );
    }

    #[test]
    fn a_direction_past_the_corner_is_pulled_back_to_it() {
        let corner = Vector::new(1.0, 1.0).normalize();
        let normals = CompoundPseudoNormals {
            boundary_edges: vec![cone(Vector::Y, corner, Vector::new(-1.0, 1.0).normalize())],
        };

        assert_eq!(normals.project_local_normal(Vector::X), Some(corner));
    }

    #[test]
    fn picks_the_edge_the_normal_points_out_of() {
        let normals = CompoundPseudoNormals {
            boundary_edges: vec![pinned(Vector::Y), pinned(Vector::X)],
        };

        assert_eq!(normals.project_local_normal(Vector::Y), Some(Vector::Y));
        assert_eq!(normals.project_local_normal(Vector::X), Some(Vector::X));
    }

    /// The case a convex decomposition creates: a wedge whose tip is a boundary edge running into a
    /// cut. A body driving into the tip reports a normal along the wedge axis, which is not a
    /// surface direction, and has to be pulled onto the one boundary edge that is.
    #[test]
    fn a_direction_reaching_the_part_from_behind_is_dropped() {
        let normals = CompoundPseudoNormals {
            boundary_edges: vec![pinned(Vector::Y)],
        };

        assert!(normals.project_local_normal(-Vector::Y).is_none());
    }

    #[test]
    fn clamps_a_wedge_tip_normal_onto_the_boundary_edge() {
        let ramp = Vector::new(-0.6, 0.8).normalize();
        let normals = CompoundPseudoNormals {
            boundary_edges: vec![pinned(ramp)],
        };

        let clamped = normals.project_local_normal(-Vector::X).unwrap();

        assert!((clamped - ramp).length() < 1.0e-5);
    }
}
