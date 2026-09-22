#![allow(clippy::multiple_bound_locations)] // for impl_downcast

use alloc::boxed::Box;
use downcast_rs::{impl_downcast, DowncastSync};

use crate::query::contact_manifolds::{
    CompositeShapeCompositeShapeContactManifoldsWorkspace,
    CompositeShapeShapeContactManifoldsWorkspace,
    HeightFieldCompositeShapeContactManifoldsWorkspace, HeightFieldShapeContactManifoldsWorkspace,
    TriMeshShapeContactManifoldsWorkspace, VoxelsShapeContactManifoldsWorkspace,
};

#[derive(Copy, Clone)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize))]
/// Enum representing workspace data of a specific type.
pub enum TypedWorkspaceData<'a> {
    /// A trimesh workspace.
    TriMeshShapeContactManifoldsWorkspace(&'a TriMeshShapeContactManifoldsWorkspace),
    /// A heightfield vs. shape workspace.
    HeightfieldShapeContactManifoldsWorkspace(&'a HeightFieldShapeContactManifoldsWorkspace),
    /// A heightfield vs. composite shape workspace.
    HeightfieldCompositeShapeContactManifoldsWorkspace(
        &'a HeightFieldCompositeShapeContactManifoldsWorkspace,
    ),
    /// A composite shape vs. composite shape workspace.
    CompositeShapeCompositeShapeContactManifoldsWorkspace(
        &'a CompositeShapeCompositeShapeContactManifoldsWorkspace,
    ),
    /// A composite shape vs. shape workspace.
    CompositeShapeShapeContactManifoldsWorkspace(&'a CompositeShapeShapeContactManifoldsWorkspace),
    /// A voxels vs. shape workspace.
    VoxelsShapeContactManifoldsWorkspace(&'a VoxelsShapeContactManifoldsWorkspace<2>),
    /// A voxels vs. composite shape workspace.
    VoxelsCompositeShapeContactManifoldsWorkspace(&'a VoxelsShapeContactManifoldsWorkspace<3>),
    /// A voxels vs. voxels workspace.
    VoxelsVoxelsContactManifoldsWorkspace(&'a VoxelsShapeContactManifoldsWorkspace<4>),
    /// A custom workspace.
    Custom,
}

// NOTE: must match the TypedWorkspaceData enum.
#[cfg(feature = "serde-serialize")]
#[derive(Deserialize)]
enum DeserializableWorkspaceData {
    TriMeshShapeContactManifoldsWorkspace(TriMeshShapeContactManifoldsWorkspace),
    HeightfieldShapeContactManifoldsWorkspace(HeightFieldShapeContactManifoldsWorkspace),
    HeightfieldCompositeShapeContactManifoldsWorkspace(
        HeightFieldCompositeShapeContactManifoldsWorkspace,
    ),
    CompositeShapeCompositeShapeContactManifoldsWorkspace(
        CompositeShapeCompositeShapeContactManifoldsWorkspace,
    ),
    CompositeShapeShapeContactManifoldsWorkspace(CompositeShapeShapeContactManifoldsWorkspace),
    VoxelsShapeContactManifoldsWorkspace(VoxelsShapeContactManifoldsWorkspace<2>),
    VoxelsCompositeShapeContactManifoldsWorkspace(VoxelsShapeContactManifoldsWorkspace<3>),
    VoxelsVoxelsContactManifoldsWorkspace(VoxelsShapeContactManifoldsWorkspace<4>),
    #[allow(dead_code)]
    Custom,
}

#[cfg(feature = "serde-serialize")]
impl DeserializableWorkspaceData {
    pub fn into_contact_manifold_workspace(self) -> ContactManifoldsWorkspace {
        match self {
            DeserializableWorkspaceData::TriMeshShapeContactManifoldsWorkspace(w) => {
                ContactManifoldsWorkspace(Box::new(w))
            }
            DeserializableWorkspaceData::HeightfieldShapeContactManifoldsWorkspace(w) => {
                ContactManifoldsWorkspace(Box::new(w))
            }
            DeserializableWorkspaceData::HeightfieldCompositeShapeContactManifoldsWorkspace(w) => {
                ContactManifoldsWorkspace(Box::new(w))
            }
            DeserializableWorkspaceData::CompositeShapeCompositeShapeContactManifoldsWorkspace(
                w,
            ) => ContactManifoldsWorkspace(Box::new(w)),
            DeserializableWorkspaceData::CompositeShapeShapeContactManifoldsWorkspace(w) => {
                ContactManifoldsWorkspace(Box::new(w))
            }
            DeserializableWorkspaceData::VoxelsShapeContactManifoldsWorkspace(w) => {
                ContactManifoldsWorkspace(Box::new(w))
            }
            DeserializableWorkspaceData::VoxelsCompositeShapeContactManifoldsWorkspace(w) => {
                ContactManifoldsWorkspace(Box::new(w))
            }
            DeserializableWorkspaceData::VoxelsVoxelsContactManifoldsWorkspace(w) => {
                ContactManifoldsWorkspace(Box::new(w))
            }
            // The data of a custom workspace is unknown to parry, so we fall back to a placeholder.
            // This can later be overwritten by algorithms that need the workspace and know its type.
            DeserializableWorkspaceData::Custom => {
                ContactManifoldsWorkspace(Box::new(UnknownContactManifoldsWorkspace))
            }
        }
    }
}

/// A placeholder workspace used for custom workspace types unknown to parry.
///
/// This is yielded when deserializing a custom [`ContactManifoldsWorkspace`],
/// as its data is lost during serialization. Algorithms that access the workspace
/// can replace the placeholder with the correct workspace type on their next call.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct UnknownContactManifoldsWorkspace;

impl WorkspaceData for UnknownContactManifoldsWorkspace {
    fn as_typed_workspace_data(&self) -> TypedWorkspaceData<'_> {
        TypedWorkspaceData::Custom
    }

    fn clone_dyn(&self) -> Box<dyn WorkspaceData> {
        Box::new(*self)
    }
}

/// Data from a [`ContactManifoldsWorkspace`].
pub trait WorkspaceData: DowncastSync {
    /// Gets the underlying workspace as an enum.
    fn as_typed_workspace_data(&self) -> TypedWorkspaceData<'_>;

    /// Clones `self`.
    fn clone_dyn(&self) -> Box<dyn WorkspaceData>;
}

impl_downcast!(sync WorkspaceData);

// Note we have this newtype because it simplifies the serialization/deserialization code.
/// A serializable workspace used by some contact-manifolds computation algorithms.
pub struct ContactManifoldsWorkspace(pub Box<dyn WorkspaceData>);

impl Clone for ContactManifoldsWorkspace {
    fn clone(&self) -> Self {
        ContactManifoldsWorkspace(self.0.clone_dyn())
    }
}

impl<T: WorkspaceData> From<T> for ContactManifoldsWorkspace {
    fn from(data: T) -> Self {
        Self(Box::new(data) as Box<dyn WorkspaceData>)
    }
}

#[cfg(feature = "serde-serialize")]
impl serde::Serialize for ContactManifoldsWorkspace {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.0.as_typed_workspace_data().serialize(serializer)
    }
}

#[cfg(feature = "serde-serialize")]
impl<'de> serde::Deserialize<'de> for ContactManifoldsWorkspace {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        Ok(DeserializableWorkspaceData::deserialize(deserializer)?
            .into_contact_manifold_workspace())
    }
}
