use crate::math::Vector;

/// Projects `dir` into the outward cone formed by a feature's `face` normal and its `closest_edge`
/// pseudo-normal, shared by [`TrianglePseudoNormals`](crate::shape::TrianglePseudoNormals) and
/// [`SegmentPseudoNormals`](crate::shape::SegmentPseudoNormals).
///
/// Returns whether `dir` lies in the front half-space; a contact whose normal points into the back
/// is meant to be discarded by the caller.
pub(crate) fn project_into_cone(face: Vector, closest_edge: Vector, dir: &mut Vector) -> bool {
    let dot_face = dir.dot(face);

    // Apply the projection. Note that this isn’t 100% accurate since this approximates the
    // vertex normal cone using the closest edge’s normal cone instead of the
    // true polygonal cone on S² (but taking into account this polygonal cone exactly
    // would be quite computationally expensive).

    if closest_edge == face {
        // The normal cone is degenerate, there is only one possible direction.
        *dir = face;
        return dot_face >= 0.0;
    }

    // TODO: take into account the two closest edges instead for better continuity
    //       of vertex normals?
    let dot_edge_face = face.dot(closest_edge);
    let dot_dir_face = face.dot(*dir);
    let dot_corrected_dir_face = 2.0 * dot_edge_face * dot_edge_face - 1.0; // cos(2 * angle(closest_edge, face))

    if dot_dir_face >= dot_corrected_dir_face {
        // The direction is in the pseudo-normal cone. No correction to apply.
        return true;
    }

    // We need to correct.
    let edge_on_normal = face * dot_edge_face;
    let edge_orthogonal_to_normal = closest_edge - edge_on_normal;

    let dir_on_normal = face * dot_dir_face;
    let dir_orthogonal_to_normal = *dir - dir_on_normal;
    let Some(unit_dir_orthogonal_to_normal) = dir_orthogonal_to_normal.try_normalize() else {
        return dot_face >= 0.0;
    };

    // NOTE: the normalization might be redundant as the result vector is guaranteed to be
    //       unit sized. Though some rounding errors could throw it off.
    let adjusted_pseudo_normal =
        edge_on_normal + unit_dir_orthogonal_to_normal * edge_orthogonal_to_normal.length();
    let (adjusted_pseudo_normal, length) = adjusted_pseudo_normal.normalize_and_length();
    if length <= 1.0e-6 {
        return dot_face >= 0.0;
    };

    // The reflection of the face normal wrt. the adjusted pseudo-normal gives us the
    // second end of the pseudo-normal cone the direction is projected on.
    *dir = adjusted_pseudo_normal * (2.0 * face.dot(adjusted_pseudo_normal)) - face;
    dot_face >= 0.0
}
