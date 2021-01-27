use crate::utils::MaybeSerializableData;

#[cfg(feature = "serde-serialize")]
use num_derive::FromPrimitive;

#[derive(Copy, Clone, Debug, FromPrimitive)]
#[cfg(feature = "serde-serialize")]
pub(super) enum WorkspaceSerializationTag {
    TriMeshShapeContactManifoldsWorkspace = 0,
    HeightfieldShapeContactManifoldsWorkspace,
    HeightfieldCompositeShapeContactManifoldsWorkspace,
    CompositeShapeCompositeShapeContactManifoldsWorkspace,
    CompositeShapeShapeContactManifoldsWorkspace,
}

// Note we have this newtype because it simplifies the serialization/deserialization code.
/// A serializable workspace used by some contact-manifolds computation algorithms.
pub struct ContactManifoldsWorkspace(pub Box<dyn MaybeSerializableData>);

impl Clone for ContactManifoldsWorkspace {
    fn clone(&self) -> Self {
        ContactManifoldsWorkspace(self.0.clone_dyn())
    }
}

impl<T: MaybeSerializableData> From<T> for ContactManifoldsWorkspace {
    fn from(data: T) -> Self {
        Self(Box::new(data) as Box<dyn MaybeSerializableData>)
    }
}

#[cfg(feature = "serde-serialize")]
impl serde::Serialize for ContactManifoldsWorkspace {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use crate::serde::ser::SerializeStruct;

        if let Some((tag, ser)) = self.0.as_serialize() {
            let mut state = serializer.serialize_struct("ContactManifoldsWorkspace", 2)?;
            state.serialize_field("tag", &tag)?;
            state.serialize_field("inner", ser)?;
            state.end()
        } else {
            Err(serde::ser::Error::custom(
                "Found a non-serializable contact generator workspace.",
            ))
        }
    }
}

#[cfg(feature = "serde-serialize")]
impl<'de> serde::Deserialize<'de> for ContactManifoldsWorkspace {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use super::{
            CompositeShapeCompositeShapeContactManifoldsWorkspace,
            CompositeShapeShapeContactManifoldsWorkspace,
            HeightFieldCompositeShapeContactManifoldsWorkspace,
            HeightFieldShapeContactManifoldsWorkspace, TriMeshShapeContactManifoldsWorkspace,
        };

        struct Visitor {};
        impl<'de> serde::de::Visitor<'de> for Visitor {
            type Value = ContactManifoldsWorkspace;
            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(formatter, "one shape type tag and the inner shape data")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: serde::de::SeqAccess<'de>,
            {
                use num::cast::FromPrimitive;

                let tag: u32 = seq
                    .next_element()?
                    .ok_or_else(|| serde::de::Error::invalid_length(0, &self))?;

                fn deser<'de, A, S: MaybeSerializableData + serde::Deserialize<'de>>(
                    seq: &mut A,
                ) -> Result<Box<dyn MaybeSerializableData>, A::Error>
                where
                    A: serde::de::SeqAccess<'de>,
                {
                    let workspace: S = seq.next_element()?.ok_or_else(|| {
                        serde::de::Error::custom("Failed to deserialize builtin workspace.")
                    })?;
                    Ok(Box::new(workspace) as Box<dyn MaybeSerializableData>)
                }

                let workspace = match WorkspaceSerializationTag::from_u32(tag) {
                    Some(WorkspaceSerializationTag::HeightfieldShapeContactManifoldsWorkspace) => {
                        deser::<A, HeightFieldShapeContactManifoldsWorkspace>(&mut seq)?
                    }
                    Some(WorkspaceSerializationTag::HeightfieldCompositeShapeContactManifoldsWorkspace) => {
                        deser::<A, HeightFieldCompositeShapeContactManifoldsWorkspace>(&mut seq)?
                    }
                    Some(WorkspaceSerializationTag::TriMeshShapeContactManifoldsWorkspace) => {
                        deser::<A, TriMeshShapeContactManifoldsWorkspace>(&mut seq)?
                    }
                    Some(WorkspaceSerializationTag::CompositeShapeCompositeShapeContactManifoldsWorkspace) => {
                        deser::<A, CompositeShapeCompositeShapeContactManifoldsWorkspace>(&mut seq)?
                    }
                    Some(WorkspaceSerializationTag::CompositeShapeShapeContactManifoldsWorkspace) => {
                        deser::<A, CompositeShapeShapeContactManifoldsWorkspace>(&mut seq)?
                    }
                    None => {
                        return Err(serde::de::Error::custom(
                            "found invalid contact generator workspace type to deserialize",
                        ))
                    }
                };

                Ok(ContactManifoldsWorkspace(workspace))
            }
        }

        deserializer.deserialize_struct("ContactManifoldsWorkspace", &["tag", "inner"], Visitor {})
    }
}
