use crate::math::{Isometry, Point, Real, Vector};
use crate::shape::{ConvexPolygonalFeature, SupportMap};
use na::Unit;

/// An identifier of a feature of a convex polyhedron.
///
/// This identifier is shape-dependent and is seach that it
/// allows an efficient retrieval of the geometric information of the
/// feature.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Copy, Clone, Debug, Hash, PartialEq, Eq)]
pub enum FeatureId {
    /// Shape-dependent identifier of a vertex.
    Vertex(u32),
    #[cfg(feature = "dim3")]
    /// Shape-dependent identifier of an edge.
    Edge(u32),
    /// Shape-dependent identifier of a face.
    Face(u32),
    // XXX: remove this variant.
    /// Unknown identifier.
    Unknown,
}

impl FeatureId {
    /// Revries the value of the identifier if `self` is a vertex.
    pub fn unwrap_vertex(self) -> u32 {
        match self {
            FeatureId::Vertex(id) => id,
            _ => panic!("The feature id does not identify a vertex."),
        }
    }

    /// Revries the value of the identifier if `self` is an edge.
    #[cfg(feature = "dim3")]
    pub fn unwrap_edge(self) -> u32 {
        match self {
            FeatureId::Edge(id) => id,
            _ => panic!("The feature id does not identify an edge."),
        }
    }

    /// Retrieves the value of the identifier if `self` is a face.
    pub fn unwrap_face(self) -> u32 {
        match self {
            FeatureId::Face(id) => id,
            _ => panic!("The feature id does not identify a face."),
        }
    }
}

/// Trait implemented by all convex polyhedron.
pub trait ConvexPolyhedron: SupportMap {
    /// Gets the specified vertex in the shape local-space.
    fn vertex(&self, id: FeatureId) -> Point<Real>;
    /// Fill `face` with the geometric description of the specified face, in the shape's local-space.
    fn face(&self, id: FeatureId, face: &mut ConvexPolygonalFeature);
    #[cfg(feature = "dim3")]
    /// Get the specified edge's vertices (in the shape local-space) and the vertices' identifiers.
    fn edge(&self, id: FeatureId) -> (Point<Real>, Point<Real>, FeatureId, FeatureId);

    /// Returns any normal from the normal cone of the given feature.
    fn feature_normal(&self, feature: FeatureId) -> Unit<Vector<Real>>;

    /// Retrieve the face (in world-space) with a normal that maximizes the scalar product with `dir`.
    fn support_face_toward(
        &self,
        transform: &Isometry<Real>,
        dir: &Unit<Vector<Real>>,
        out: &mut ConvexPolygonalFeature,
    );

    /// Retrieve the feature (in world-space) which normal cone contains `dir`.
    fn support_feature_toward(
        &self,
        transform: &Isometry<Real>,
        dir: &Unit<Vector<Real>>,
        _angle: Real,
        out: &mut ConvexPolygonalFeature,
    );

    /// Retrieve the identifier of the feature which normal cone contains `dir`.
    fn support_feature_id_toward(&self, local_dir: &Unit<Vector<Real>>) -> FeatureId;
}
