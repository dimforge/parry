#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

/// An identifier of a feature of a convex polyhedron.
///
/// This identifier is shape-dependent and is such that it
/// allows an efficient retrieval of the geometric information of the
/// feature.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(as = "Self")
)]
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

impl Default for FeatureId {
    fn default() -> Self {
        FeatureId::Unknown
    }
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

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
#[derive(Copy, Clone, Debug, Hash, PartialEq, Eq)]
/// A feature id where the feature type is packed into the same value as the feature index.
pub struct PackedFeatureId(pub u32);

impl PackedFeatureId {
    /// Packed feature id identifying an unknown feature.
    pub const UNKNOWN: Self = Self(0);

    const CODE_MASK: u32 = 0x3fff_ffff;
    const HEADER_MASK: u32 = !Self::CODE_MASK;
    const HEADER_VERTEX: u32 = 0b01 << 30;
    #[cfg(feature = "dim3")]
    const HEADER_EDGE: u32 = 0b10 << 30;
    const HEADER_FACE: u32 = 0b11 << 30;

    /// Converts a vertex feature id into a packed feature id.
    pub fn vertex(code: u32) -> Self {
        assert_eq!(code & Self::HEADER_MASK, 0);
        Self(Self::HEADER_VERTEX | code)
    }

    /// Converts a edge feature id into a packed feature id.
    #[cfg(feature = "dim3")]
    pub fn edge(code: u32) -> Self {
        assert_eq!(code & Self::HEADER_MASK, 0);
        Self(Self::HEADER_EDGE | code)
    }

    /// Converts a face feature id into a packed feature id.
    pub fn face(code: u32) -> Self {
        assert_eq!(code & Self::HEADER_MASK, 0);
        Self(Self::HEADER_FACE | code)
    }

    #[cfg(feature = "dim2")]
    /// Converts an array of vertex feature ids into an array of packed feature ids.
    pub(crate) fn vertices(code: [u32; 2]) -> [Self; 2] {
        [Self::vertex(code[0]), Self::vertex(code[1])]
    }

    #[cfg(feature = "dim3")]
    /// Converts an array of vertex feature ids into an array of packed feature ids.
    pub(crate) fn vertices(code: [u32; 4]) -> [Self; 4] {
        [
            Self::vertex(code[0]),
            Self::vertex(code[1]),
            Self::vertex(code[2]),
            Self::vertex(code[3]),
        ]
    }

    #[cfg(feature = "dim3")]
    /// Converts an array of edge feature ids into an array of packed feature ids.
    pub(crate) fn edges(code: [u32; 4]) -> [Self; 4] {
        [
            Self::edge(code[0]),
            Self::edge(code[1]),
            Self::edge(code[2]),
            Self::edge(code[3]),
        ]
    }

    /// Unpacks this feature id into an explicit enum.
    pub fn unpack(self) -> FeatureId {
        let header = self.0 & Self::HEADER_MASK;
        let code = self.0 & Self::CODE_MASK;
        match header {
            Self::HEADER_VERTEX => FeatureId::Vertex(code),
            #[cfg(feature = "dim3")]
            Self::HEADER_EDGE => FeatureId::Edge(code),
            Self::HEADER_FACE => FeatureId::Face(code),
            _ => FeatureId::Unknown,
        }
    }

    /// Is the identified feature a face?
    pub fn is_face(self) -> bool {
        self.0 & Self::HEADER_MASK == Self::HEADER_FACE
    }

    /// Is the identified feature a vertex?
    pub fn is_vertex(self) -> bool {
        self.0 & Self::HEADER_MASK == Self::HEADER_VERTEX
    }

    /// Is the identified feature an edge?
    #[cfg(feature = "dim3")]
    pub fn is_edge(self) -> bool {
        self.0 & Self::HEADER_MASK == Self::HEADER_EDGE
    }

    /// Is the identified feature unknown?
    pub fn is_unknown(self) -> bool {
        self == Self::UNKNOWN
    }
}

impl From<FeatureId> for PackedFeatureId {
    fn from(value: FeatureId) -> Self {
        match value {
            FeatureId::Face(fid) => Self::face(fid),
            #[cfg(feature = "dim3")]
            FeatureId::Edge(fid) => Self::edge(fid),
            FeatureId::Vertex(fid) => Self::vertex(fid),
            FeatureId::Unknown => Self::UNKNOWN,
        }
    }
}
