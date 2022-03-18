use crate::math::{Isometry, Point, Real, Vector};
use crate::query::{self, visitors::BoundingVolumeIntersectionsSimultaneousVisitor, PointQuery};
use crate::shape::trimesh::{TopoFace, TopoHalfEdge, TopoVertex, TriMeshTopology};
use crate::shape::{FeatureId, Segment, TriMesh, Triangle};
use core::fmt;
use std::collections::{HashMap, HashSet, VecDeque};

/// Error indicating that a query is not supported between certain shapes
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum MeshIntersectionError {
    MissingTopology,
    TriTriError,
}

impl fmt::Display for MeshIntersectionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingTopology => {
                f.pad("at least on of the meshes is missing its topology information")
            }
            Self::TriTriError => f.pad("internal failure while intersecting two triangles"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for MeshIntersectionError {}
