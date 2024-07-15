use crate::shape::TriMesh;
use obj::{Group, IndexTuple, ObjData, ObjError, Object, SimplePolygon};
use std::path::PathBuf;

impl TriMesh {
    pub fn to_obj_file(&self, path: &PathBuf) -> Result<(), ObjError> {
        let mut file = std::fs::File::create(path).unwrap();

        ObjData {
            position: self
                .vertices()
                .into_iter()
                .map(|v| [v.x as f32, v.y as f32, v.z as f32])
                .collect(),
            objects: vec![Object {
                groups: vec![Group {
                    polys: self
                        .indices()
                        .into_iter()
                        .map(|tri| {
                            SimplePolygon(vec![
                                IndexTuple(tri[0] as usize, None, None),
                                IndexTuple(tri[1] as usize, None, None),
                                IndexTuple(tri[2] as usize, None, None),
                            ])
                        })
                        .collect(),
                    name: "".to_string(),
                    index: 0,
                    material: None,
                }],
                name: "".to_string(),
            }],
            ..Default::default()
        }
        .write_to_buf(&mut file)
    }
}
