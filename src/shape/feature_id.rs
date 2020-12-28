/// An identifier of a feature of a convex polyhedron.
///
/// This identifier is shape-dependent and is such that it
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
