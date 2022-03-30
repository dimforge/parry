use core::fmt;

/// Error indicating that a query is not supported between certain shapes
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum MeshIntersectionError {
    MissingTopology,
    MissingPseudoNormals,
    TriTriError,
}

impl fmt::Display for MeshIntersectionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingTopology => {
                f.pad("at least on of the meshes is missing its topology information. Call `mesh.compute_topology` on the mesh")
            }
            Self::MissingPseudoNormals => {
                f.pad("at least on of the meshes is missing its pseudo-normals. Call `mesh.compute_pseudo_normals` on the mesh")
            }
            Self::TriTriError => f.pad("internal failure while intersecting two triangles"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for MeshIntersectionError {}
