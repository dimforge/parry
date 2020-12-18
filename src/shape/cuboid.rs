//! Support mapping based Cuboid shape.

use crate::math::{Isometry, Point, Real, Vector, DIM};
use crate::shape::{FeatureId, SupportMap, PolygonalFeature};
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

impl Cuboid {
    /// Creates a new box from its half-extents. Half-extents are the box half-width along each
    /// axis. Each half-extent must be positive.
    #[inline]
    pub fn new(half_extents: Vector<Real>) -> Cuboid {
        Cuboid { half_extents }
    }

    #[cfg(feature = "dim2")]
    pub fn vertex_feature_id(vertex: Point<Real>) -> u8 {
        ((vertex.x.to_bits() >> 31) & 0b001 | (vertex.y.to_bits() >> 30) & 0b010) as u8
    }

    #[cfg(feature = "dim2")]
    pub fn support_feature(&self, local_dir: Vector<Real>) -> PolygonalFeature {
        // In 2D, it is best for stability to always return a face.
        // It won't have any notable impact on performances anyway.
        self.support_face(local_dir)
    }

    #[cfg(feature = "dim2")]
    pub fn support_face(&self, local_dir: Vector<Real>) -> PolygonalFeature {
        let he = self.half_extents;
        let i = local_dir.iamin();
        let j = (i + 1) % 2;
        let mut a = Point::origin();
        a[i] = he[i];
        a[j] = he[j].copysign(local_dir[j]);

        let mut b = a;
        b[i] = -he[i];

        let vid1 = Self::vertex_feature_id(a);
        let vid2 = Self::vertex_feature_id(b);
        let fid = (vid1.max(vid2) << 2) | vid1.min(vid2) | 0b11_00_00;

        PolygonalFeature {
            vertices: [a, b],
            vids: [vid1, vid2],
            fid,
            num_vertices: 2,
        }
    }

    #[cfg(feature = "dim3")]
    pub fn support_feature(&self, local_dir: Vector<Real>) -> PolygonalFeature {
        // FIXME: this should actually return the feature.
        // And we should change all the callers of this method to use
        // `.support_face` instead of this method to preserve their old behavior.
        self.support_face(local_dir)
        /*
        const MAX_DOT_THRESHOLD: Real = crate::utils::COS_10_DEGREES;
        const MIN_DOT_THRESHOLD: Real = 1.0 - MAX_DOT_THRESHOLD;

        let amax = local_dir.amax();
        let amin = local_dir.amin();

        if amax > MAX_DOT_THRESHOLD {
            // Support face.
            CuboidFeature::Face(support_face(self, local_dir))
        } else if amin < MIN_DOT_THRESHOLD {
            // Support edge.
            CuboidFeature::Edge(support_edge(self, local_dir))
        } else {
            // Support vertex.
            CuboidFeature::Vertex(support_vertex(self, local_dir))
        }
        */
    }

    // #[cfg(feature = "dim3")
    // pub(crate) fn support_vertex(&self, local_dir: Vector<Real>) -> CuboidFeatureVertex {
    //     let vertex = local_support_point(self, local_dir);
    //     let vid = vertex_feature_id(vertex);
    //
    //     CuboidFeatureVertex { vertex, vid }
    // }

    // #[cfg(feature = "dim3")
    // pub(crate) fn support_edge(&self, local_dir: Vector<Real>) -> CuboidFeatureEdge {
    //     let he = self.half_extents;
    //     let i = local_dir.iamin();
    //     let j = (i + 1) % 3;
    //     let k = (i + 2) % 3;
    //     let mut a = Point::origin();
    //     a[i] = he[i];
    //     a[j] = he[j].copysign(local_dir[j]);
    //     a[k] = he[k].copysign(local_dir[k]);
    //
    //     let mut b = a;
    //     b[i] = -he[i];
    //
    //     let vid1 = vertex_feature_id(a);
    //     let vid2 = vertex_feature_id(b);
    //     let eid = (vid1.max(vid2) << 3) | vid1.min(vid2) | 0b11_000_000;
    //
    //     CuboidFeatureEdge {
    //         vertices: [a, b],
    //         vids: [vid1, vid2],
    //         eid,
    //     }
    // }

    #[cfg(feature = "dim3")]
    pub fn support_face(&self, local_dir: Vector<Real>) -> PolygonalFeature {
        // NOTE: can we use the orthonormal basis of local_dir
        // to make this AoSoA friendly?
        let he = self.half_extents;
        let iamax = local_dir.iamax();
        let sign = (1.0 as Real).copysign(local_dir[iamax]);

        let vertices = match iamax {
            0 => [
                Point::new(he.x * sign, he.y, he.z),
                Point::new(he.x * sign, -he.y, he.z),
                Point::new(he.x * sign, -he.y, -he.z),
                Point::new(he.x * sign, he.y, -he.z),
            ],
            1 => [
                Point::new(he.x, he.y * sign, he.z),
                Point::new(-he.x, he.y * sign, he.z),
                Point::new(-he.x, he.y * sign, -he.z),
                Point::new(he.x, he.y * sign, -he.z),
            ],
            2 => [
                Point::new(he.x, he.y, he.z * sign),
                Point::new(he.x, -he.y, he.z * sign),
                Point::new(-he.x, -he.y, he.z * sign),
                Point::new(-he.x, he.y, he.z * sign),
            ],
            _ => unreachable!(),
        };

        pub fn vid(i: u8) -> u8 {
            // Each vertex has an even feature id.
            i * 2
        }

        let sign_index = ((sign as i8 + 1) / 2) as usize;
        // The vertex id as numbered depending on the sign of the vertex
        // component. A + sign means the corresponding bit is 0 while a -
        // sign means the corresponding bit is 1.
        // For exampl the vertex [2.0, -1.0, -3.0] has the id 0b011
        let vids = match iamax {
            0 => [
                [vid(0b000), vid(0b010), vid(0b011), vid(0b001)],
                [vid(0b100), vid(0b110), vid(0b111), vid(0b101)],
            ][sign_index],
            1 => [
                [vid(0b000), vid(0b100), vid(0b101), vid(0b001)],
                [vid(0b010), vid(0b110), vid(0b111), vid(0b011)],
            ][sign_index],
            2 => [
                [vid(0b000), vid(0b010), vid(0b110), vid(0b100)],
                [vid(0b001), vid(0b011), vid(0b111), vid(0b101)],
            ][sign_index],
            _ => unreachable!(),
        };

        // The feature ids of edges is obtained from the vertex ids
        // of their endpoints.
        // Assuming vid1 > vid2, we do:   (vid1 << 3) | vid2 | 0b11000000
        //
        let eids = match iamax {
            0 => [
                [0b11_010_000, 0b11_011_010, 0b11_011_001, 0b11_001_000],
                [0b11_110_100, 0b11_111_110, 0b11_111_101, 0b11_101_100],
            ][sign_index],
            1 => [
                [0b11_100_000, 0b11_101_100, 0b11_101_001, 0b11_001_000],
                [0b11_110_010, 0b11_111_110, 0b11_111_011, 0b11_011_010],
            ][sign_index],
            2 => [
                [0b11_010_000, 0b11_110_010, 0b11_110_100, 0b11_100_000],
                [0b11_011_001, 0b11_111_011, 0b11_111_101, 0b11_101_001],
            ][sign_index],
            _ => unreachable!(),
        };

        // The face with normals [x, y, z] are numbered [10, 11, 12].
        // The face with negated normals are numbered [13, 14, 15].
        let fid = iamax + sign_index * 3 + 10;

        PolygonalFeature {
            vertices,
            vids,
            eids,
            fid: fid as u8,
            num_vertices: 4,
        }
    }
}

impl SupportMap for Cuboid {
    #[inline]
    fn local_support_point(&self, dir: &Vector<Real>) -> Point<Real> {
        dir.copy_sign_to(self.half_extents).into()
    }
}

/*
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
*/