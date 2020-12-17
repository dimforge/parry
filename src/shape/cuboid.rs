//! Support mapping based Cuboid shape.

use crate::math::{Isometry, Point, Real, Vector, DIM};
use crate::shape::{ConvexPolygonalFeature, ConvexPolyhedron, FeatureId, SupportMap};
use crate::utils::WSign;
use na::{self, Unit};
use std::f64;

/// Shape of a box.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(PartialEq, Debug, Copy, Clone)]
pub struct Cuboid {
    /// The half-extents of the cuboid.
    pub half_extents: Vector<Real>,
}

// NOTE: format of the cuboid feature id:
//
// FeatureId::Vertex(id): the i-th bit of `id` is set to 1 iff. the i-th component of the vertex is negative.
// FeatureId::Edge(id): the part `id & 0b11` contains a number in [0,1] (or [0,2] in 3D) to indicate the axis (x, y, z).
//                      the part `id >> 2` follow the same rule as the vertex id.
// FeatureId::Face(id): if `id` lies in [0,1] (or [0 2] in 3D) indicates the axis (x, y, z) corresponding to the face normal.
//                      If `id` is greater than 1 (or 2 in 3D), then the negative axis (-x, -y, -z) is given by `id - 2` (or `id - 3` in 3D).
impl Cuboid {
    /// Creates a new box from its half-extents. Half-extents are the box half-width along each
    /// axis. Each half-extent must be positive.
    #[inline]
    pub fn new(half_extents: Vector<Real>) -> Cuboid {
        Cuboid { half_extents }
    }
}

impl SupportMap for Cuboid {
    #[inline]
    fn local_support_point(&self, dir: &Vector<Real>) -> Point<Real> {
        dir.copy_sign_to(self.half_extents).into()
    }
}

impl ConvexPolyhedron for Cuboid {
    fn vertex(&self, id: FeatureId) -> Point<Real> {
        let vid = id.unwrap_vertex();
        let mut res = self.half_extents;

        for i in 0..DIM {
            if vid & (1 << i) != 0 {
                res[i] = -res[i]
            }
        }

        Point::from(res)
    }

    #[cfg(feature = "dim3")]
    fn edge(&self, id: FeatureId) -> (Point<Real>, Point<Real>, FeatureId, FeatureId) {
        let eid = id.unwrap_edge();
        let mut res = self.half_extents;

        let edge_i = eid & 0b11;
        let vertex_i = eid >> 2;

        for i in 0..DIM {
            if i as u32 != edge_i && (vertex_i & (1 << i) != 0) {
                res[i] = -res[i]
            }
        }

        let p1 = Point::from(res);
        res[edge_i as usize] = -res[edge_i as usize];
        let p2 = Point::from(res);
        let vid1 = FeatureId::Vertex(vertex_i & !(1 << edge_i));
        let vid2 = FeatureId::Vertex(vertex_i | (1 << edge_i));

        (p1, p2, vid1, vid2)
    }

    fn face(&self, id: FeatureId, out: &mut ConvexPolygonalFeature) {
        out.clear();

        let i = id.unwrap_face() as usize;
        let i1;
        let sign;

        if i < DIM {
            i1 = i;
            sign = na::one::<Real>();
        } else {
            i1 = i - DIM;
            sign = -na::one::<Real>();
        }

        #[cfg(feature = "dim2")]
        {
            let i2 = (i1 + 1) % 2;

            let mut vertex = self.half_extents;
            vertex[i1] *= sign;
            vertex[i2] *= if i1 == 0 { -sign } else { sign };

            let p1 = Point::from(vertex);
            vertex[i2] = -vertex[i2];
            let p2 = Point::from(vertex);

            let mut vertex_id1 = if sign < na::zero::<Real>() {
                1 << i1
            } else {
                0
            };
            let mut vertex_id2 = vertex_id1;
            if p1[i2] < na::zero::<Real>() {
                vertex_id1 |= 1 << i2;
            } else {
                vertex_id2 |= 1 << i2;
            }

            out.push(p1, FeatureId::Vertex(vertex_id1));
            out.push(p2, FeatureId::Vertex(vertex_id2));

            let mut normal: Vector<Real> = na::zero();
            normal[i1] = sign;
            out.set_normal(Unit::new_unchecked(normal));
            out.set_feature_id(FeatureId::Face(i as u32));
        }
        #[cfg(feature = "dim3")]
        {
            let i2 = (i1 + 1) % 3;
            let i3 = (i1 + 2) % 3;
            let (edge_i2, edge_i3) = if sign > na::zero::<Real>() {
                (i2, i3)
            } else {
                (i3, i2)
            };
            let mask_i2 = !(1 << edge_i2); // The masks are for ensuring each edge has a unique ID.
            let mask_i3 = !(1 << edge_i3);
            let mut vertex = self.half_extents;
            vertex[i1] *= sign;

            let (sbit, msbit) = if sign < na::zero::<Real>() {
                (1, 0)
            } else {
                (0, 1)
            };
            let mut vertex_id = sbit << i1;
            out.push(Point::from(vertex), FeatureId::Vertex(vertex_id));
            out.push_edge_feature_id(FeatureId::Edge(
                edge_i2 as u32 | ((vertex_id & mask_i2) << 2),
            ));

            vertex[i2] = -sign * self.half_extents[i2];
            vertex[i3] = sign * self.half_extents[i3];
            vertex_id |= msbit << i2 | sbit << i3;
            out.push(Point::from(vertex), FeatureId::Vertex(vertex_id));
            out.push_edge_feature_id(FeatureId::Edge(
                edge_i3 as u32 | ((vertex_id & mask_i3) << 2),
            ));

            vertex[i2] = -self.half_extents[i2];
            vertex[i3] = -self.half_extents[i3];
            vertex_id |= 1 << i2 | 1 << i3;
            out.push(Point::from(vertex), FeatureId::Vertex(vertex_id));
            out.push_edge_feature_id(FeatureId::Edge(
                edge_i2 as u32 | ((vertex_id & mask_i2) << 2),
            ));

            vertex[i2] = sign * self.half_extents[i2];
            vertex[i3] = -sign * self.half_extents[i3];
            vertex_id = sbit << i1 | sbit << i2 | msbit << i3;
            out.push(Point::from(vertex), FeatureId::Vertex(vertex_id));
            out.push_edge_feature_id(FeatureId::Edge(
                edge_i3 as u32 | ((vertex_id & mask_i3) << 2),
            ));

            let mut normal: Vector<Real> = na::zero();
            normal[i1] = sign;
            out.set_normal(Unit::new_unchecked(normal));

            if sign > na::zero::<Real>() {
                out.set_feature_id(FeatureId::Face(i1 as u32));
            } else {
                out.set_feature_id(FeatureId::Face(i1 as u32 + 3));
            }

            out.recompute_edge_normals();
        }
    }

    fn support_face_toward(
        &self,
        m: &Isometry<Real>,
        dir: &Unit<Vector<Real>>,
        out: &mut ConvexPolygonalFeature,
    ) {
        out.clear();
        let local_dir = m.inverse_transform_vector(dir);

        let mut iamax = 0;
        let mut amax = local_dir[0].abs();

        // FIXME: we should use nalgebra's iamax method.
        for i in 1..DIM {
            let candidate = local_dir[i].abs();
            if candidate > amax {
                amax = candidate;
                iamax = i;
            }
        }

        if local_dir[iamax] > na::zero::<Real>() {
            self.face(FeatureId::Face(iamax as u32), out);
            out.transform_by(m);
        } else {
            self.face(FeatureId::Face((iamax + DIM) as u32), out);
            out.transform_by(m);
        }
    }

    fn support_feature_toward(
        &self,
        m: &Isometry<Real>,
        dir: &Unit<Vector<Real>>,
        angle: Real,
        out: &mut ConvexPolygonalFeature,
    ) {
        let local_dir = m.inverse_transform_vector(dir);
        let cang = angle.cos();
        let mut support_point = self.half_extents;

        out.clear();

        #[cfg(feature = "dim2")]
        {
            let mut support_point_id = 0;
            for i1 in 0..2 {
                let sign = local_dir[i1].signum();
                if sign * local_dir[i1] >= cang {
                    if sign > na::zero::<Real>() {
                        self.face(FeatureId::Face(i1 as u32), out);
                        out.transform_by(m);
                    } else {
                        self.face(FeatureId::Face(i1 as u32 + 2), out);
                        out.transform_by(m);
                    }
                    return;
                } else {
                    if sign < na::zero::<Real>() {
                        support_point_id |= 1 << i1;
                    }
                    support_point[i1] *= sign;
                }
            }

            // We are not on a face, return the support vertex.
            out.push(
                m * Point::from(support_point),
                FeatureId::Vertex(support_point_id),
            );
            out.set_feature_id(FeatureId::Vertex(support_point_id));
        }

        #[cfg(feature = "dim3")]
        {
            let sang = angle.sin();
            let mut support_point_id = 0;

            // Check faces.
            for i1 in 0..3 {
                let sign = local_dir[i1].signum();
                if sign * local_dir[i1] >= cang {
                    if sign > na::zero::<Real>() {
                        self.face(FeatureId::Face(i1 as u32), out);
                        out.transform_by(m);
                    } else {
                        self.face(FeatureId::Face(i1 as u32 + 3), out);
                        out.transform_by(m);
                    }
                    return;
                } else {
                    if sign < na::zero::<Real>() {
                        support_point[i1] *= sign;
                        support_point_id |= 1 << i1;
                    }
                }
            }

            // Check edges.
            for i in 0..3 {
                let sign = local_dir[i].signum();

                // sign * local_dir[i] <= cos(pi / 2 - angle)
                if sign * local_dir[i] <= sang {
                    support_point[i] = -self.half_extents[i];
                    let p1 = Point::from(support_point);
                    support_point[i] = self.half_extents[i];
                    let p2 = Point::from(support_point);
                    let p2_id = support_point_id & !(1 << i);
                    out.push(m * p1, FeatureId::Vertex(support_point_id | (1 << i)));
                    out.push(m * p2, FeatureId::Vertex(p2_id));

                    let edge_id = FeatureId::Edge(i as u32 | (p2_id << 2));
                    out.push_edge_feature_id(edge_id);
                    out.set_feature_id(edge_id);
                    return;
                }
            }

            // We are not on a face or edge, return the support vertex.
            out.push(
                m * Point::from(support_point),
                FeatureId::Vertex(support_point_id),
            );
            out.set_feature_id(FeatureId::Vertex(support_point_id));
        }
    }

    fn support_feature_id_toward(&self, local_dir: &Unit<Vector<Real>>) -> FeatureId {
        let one_degree: Real = na::convert::<f64, Real>(f64::consts::PI / 180.0);
        let cang = one_degree.cos();

        #[cfg(feature = "dim2")]
        {
            let mut support_point_id = 0;
            for i1 in 0..2 {
                let sign = local_dir[i1].signum();
                if sign * local_dir[i1] >= cang {
                    if sign > na::zero::<Real>() {
                        return FeatureId::Face(i1 as u32);
                    } else {
                        return FeatureId::Face(i1 as u32 + 2);
                    }
                } else {
                    if sign < na::zero::<Real>() {
                        support_point_id |= 1 << i1;
                    }
                }
            }

            // We are not on a face, return the support vertex.
            FeatureId::Vertex(support_point_id)
        }

        #[cfg(feature = "dim3")]
        {
            let sang = one_degree.sin();
            let mut support_point_id = 0;

            // Check faces.
            for i1 in 0..3 {
                let sign = local_dir[i1].signum();
                if sign * local_dir[i1] >= cang {
                    if sign > na::zero::<Real>() {
                        return FeatureId::Face(i1 as u32);
                    } else {
                        return FeatureId::Face(i1 as u32 + 3);
                    }
                } else {
                    if sign < na::zero::<Real>() {
                        support_point_id |= 1 << i1;
                    }
                }
            }

            // Check edges.
            for i in 0..3 {
                let sign = local_dir[i].signum();

                // sign * local_dir[i] <= cos(pi / 2 - angle)
                if sign * local_dir[i] <= sang {
                    let mask_i = !(1 << i); // To ensure each edge has a unique id.
                    return FeatureId::Edge(i as u32 | ((support_point_id & mask_i) << 2));
                }
            }

            FeatureId::Vertex(support_point_id)
        }
    }

    #[cfg(feature = "dim2")]
    fn feature_normal(&self, feature: FeatureId) -> Unit<Vector<Real>> {
        match feature {
            FeatureId::Face(id) => {
                let mut dir: Vector<Real> = na::zero();

                if id < 2 {
                    dir[id as usize] = na::one::<Real>();
                } else {
                    dir[id as usize - 2] = -na::one::<Real>();
                }
                Unit::new_unchecked(dir)
            }
            FeatureId::Vertex(id) => {
                let mut dir: Vector<Real> = na::zero();

                match id {
                    0b00 => {
                        dir[0] = na::one::<Real>();
                        dir[1] = na::one::<Real>();
                    }
                    0b01 => {
                        dir[1] = na::one::<Real>();
                        dir[0] = -na::one::<Real>();
                    }
                    0b11 => {
                        dir[0] = -na::one::<Real>();
                        dir[1] = -na::one::<Real>();
                    }
                    0b10 => {
                        dir[1] = -na::one::<Real>();
                        dir[0] = na::one::<Real>();
                    }
                    _ => panic!("Invalid feature ID: {:?}", feature),
                }

                Unit::new_normalize(dir)
            }
            _ => panic!("Invalid feature ID {:?}.", feature),
        }
    }

    #[cfg(feature = "dim3")]
    fn feature_normal(&self, feature: FeatureId) -> Unit<Vector<Real>> {
        match feature {
            FeatureId::Face(id) => {
                let mut dir: Vector<Real> = na::zero();

                if id < 3 {
                    dir[id as usize] = na::one::<Real>();
                } else {
                    dir[id as usize - 3] = -na::one::<Real>();
                }
                Unit::new_unchecked(dir)
            }
            FeatureId::Edge(id) => {
                let edge = id & 0b011;
                let face1 = (edge + 1) % 3;
                let face2 = (edge + 2) % 3;
                let signs = id >> 2;

                let mut dir: Vector<Real> = na::zero();
                let _1: Real = na::one();

                if signs & (1 << face1) != 0 {
                    dir[face1 as usize] = -_1
                } else {
                    dir[face1 as usize] = _1
                }

                if signs & (1 << face2) != 0 {
                    dir[face2 as usize] = -_1
                } else {
                    dir[face2 as usize] = _1;
                }

                Unit::new_normalize(dir)
            }
            FeatureId::Vertex(id) => {
                let mut dir: Vector<Real> = na::zero();
                for i in 0..3 {
                    let _1: Real = na::one();

                    if id & (1 << i) != 0 {
                        dir[i] = -_1;
                    } else {
                        dir[i] = _1
                    }
                }

                Unit::new_normalize(dir)
            }
            _ => panic!("Invalid feature ID: {:?}", feature),
        }
    }
}
