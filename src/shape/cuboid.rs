//! Support mapping based Cuboid shape.

#[cfg(feature = "dim3")]
use crate::math::Real;
use crate::math::{Vector, VectorExt};
#[cfg(feature = "dim3")]
use crate::shape::Segment;
use crate::shape::{FeatureId, PackedFeatureId, PolygonalFeature, SupportMap};
use crate::utils::WSign;

/// A cuboid shape, also known as a box or rectangle.
///
/// A cuboid is defined by its **half-extents**, which are half the width, height
/// (and depth in 3D) along each axis. The cuboid is always axis-aligned in its
/// local coordinate system and centered at the origin.
///
/// # Properties
///
/// - **In 2D**: Represents a rectangle with dimensions `2 * half_extents.x` by `2 * half_extents.y`
/// - **In 3D**: Represents a box with dimensions `2 * half_extents.x/y/z`
/// - **Convex**: Yes, cuboids are always convex shapes
/// - **Axis-aligned**: In local space, yes (but can be rotated via transformation)
///
/// # Why Half-Extents?
///
/// Using half-extents instead of full dimensions makes many calculations simpler
/// and more efficient. For example, checking if a point is inside a cuboid becomes:
/// `abs(point.x) <= half_extents.x && abs(point.y) <= half_extents.y`
///
/// # Use Cases
///
/// Cuboids are ideal for:
/// - Boxes, crates, and containers
/// - Walls, floors, and platforms
/// - Simple collision bounds for complex objects
/// - AABB (Axis-Aligned Bounding Box) representations
///
/// # Example
///
/// ```rust
/// # #[cfg(all(feature = "dim3", feature = "f32"))] {
/// use parry3d::shape::Cuboid;
/// use parry3d::math::Vector;
///
/// // Create a box that is 4 units wide, 2 units tall, and 6 units deep
/// // (half-extents are half of each dimension)
/// let cuboid = Cuboid::new(Vector::new(2.0, 1.0, 3.0));
///
/// assert_eq!(cuboid.half_extents.x, 2.0);
/// assert_eq!(cuboid.half_extents.y, 1.0);
/// assert_eq!(cuboid.half_extents.z, 3.0);
///
/// // Full dimensions would be:
/// // width = 4.0, height = 2.0, depth = 6.0
/// # }
/// ```
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "bytemuck", derive(bytemuck::Pod, bytemuck::Zeroable))]
#[cfg_attr(feature = "encase", derive(encase::ShaderType))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize)
)]
#[derive(PartialEq, Debug, Copy, Clone)]
#[repr(C)]
pub struct Cuboid {
    /// The half-extents of the cuboid along each axis.
    ///
    /// Each component represents half the dimension along that axis:
    /// - `half_extents.x`: Half the width
    /// - `half_extents.y`: Half the height
    /// - `half_extents.z`: Half the depth (3D only)
    ///
    /// All components should be positive.
    pub half_extents: Vector,
}

impl Cuboid {
    /// Creates a new cuboid from its half-extents.
    ///
    /// Half-extents represent half the width along each axis. To create a cuboid
    /// with full dimensions (width, height, depth), divide each by 2.
    ///
    /// # Arguments
    ///
    /// * `half_extents` - Half the dimensions along each axis. All components should be positive.
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::shape::Cuboid;
    /// use parry3d::math::Vector;
    ///
    /// // Create a 10x6x4 box (full dimensions)
    /// let cuboid = Cuboid::new(Vector::new(5.0, 3.0, 2.0));
    ///
    /// // Verify the half-extents
    /// assert_eq!(cuboid.half_extents.x, 5.0);
    /// assert_eq!(cuboid.half_extents.y, 3.0);
    /// assert_eq!(cuboid.half_extents.z, 2.0);
    /// # }
    /// ```
    ///
    /// ```
    /// # #[cfg(all(feature = "dim2", feature = "f32"))] {
    /// // In 2D:
    /// use parry2d::shape::Cuboid;
    /// use parry2d::math::Vector;
    ///
    /// // Create a 20x10 rectangle
    /// let rect = Cuboid::new(Vector::new(10.0, 5.0));
    /// assert_eq!(rect.half_extents.x, 10.0);
    /// assert_eq!(rect.half_extents.y, 5.0);
    /// # }
    /// ```
    #[inline]
    pub fn new(half_extents: Vector) -> Cuboid {
        Cuboid { half_extents }
    }

    /// Computes a scaled version of this cuboid.
    ///
    /// Each dimension is multiplied by the corresponding component of the `scale` vector.
    /// Unlike balls, cuboids can be scaled non-uniformly (different scale factors per axis)
    /// and still remain valid cuboids.
    ///
    /// # Arguments
    ///
    /// * `scale` - The scaling factors for each axis
    ///
    /// # Returns
    ///
    /// A new cuboid with scaled dimensions
    ///
    /// # Example
    ///
    /// ```
    /// # #[cfg(all(feature = "dim3", feature = "f32"))] {
    /// use parry3d::shape::Cuboid;
    /// use parry3d::math::Vector;
    ///
    /// let cuboid = Cuboid::new(Vector::new(1.0, 2.0, 3.0));
    ///
    /// // Uniform scaling: double all dimensions
    /// let scaled_uniform = cuboid.scaled(Vector::new(2.0, 2.0, 2.0));
    /// assert_eq!(scaled_uniform.half_extents, Vector::new(2.0, 4.0, 6.0));
    ///
    /// // Non-uniform scaling: different scale per axis
    /// let scaled_non_uniform = cuboid.scaled(Vector::new(2.0, 1.0, 0.5));
    /// assert_eq!(scaled_non_uniform.half_extents, Vector::new(2.0, 2.0, 1.5));
    /// # }
    /// ```
    pub fn scaled(self, scale: Vector) -> Self {
        let new_hext = self.half_extents * scale;
        Self {
            half_extents: new_hext,
        }
    }

    /// Return the id of the vertex of this cuboid with a normal that maximizes
    /// the dot product with `dir`.
    #[cfg(feature = "dim2")]
    pub fn vertex_feature_id(vertex: Vector) -> u32 {
        // TODO: is this still correct with the f64 version?
        #[allow(clippy::unnecessary_cast)] // Unnecessary for f32 but necessary for f64.
        {
            ((vertex.x.to_bits() >> 31) & 0b001 | (vertex.y.to_bits() >> 30) & 0b010) as u32
        }
    }

    /// Return the feature of this cuboid with a normal that maximizes
    /// the dot product with `dir`.
    #[cfg(feature = "dim2")]
    pub fn support_feature(&self, local_dir: Vector) -> PolygonalFeature {
        // In 2D, it is best for stability to always return a face.
        // It won't have any notable impact on performances anyway.
        self.support_face(local_dir)
    }

    /// Return the face of this cuboid with a normal that maximizes
    /// the dot product with `local_dir`.
    #[cfg(feature = "dim2")]
    pub fn support_face(&self, local_dir: Vector) -> PolygonalFeature {
        let he = self.half_extents;
        let i = local_dir.abs().min_position();

        let vertices = match i {
            0 => [
                Vector::new(he.x, local_dir.y.copy_sign_to(he.y)),
                Vector::new(-he.x, local_dir.y.copy_sign_to(he.y)),
            ],
            _ => [
                Vector::new(local_dir.x.copy_sign_to(he.x), he.y),
                Vector::new(local_dir.x.copy_sign_to(he.x), -he.y),
            ],
        };

        let vid1 = Self::vertex_feature_id(vertices[0]);
        let vid2 = Self::vertex_feature_id(vertices[1]);
        let fid = (vid1.max(vid2) << 2) | vid1.min(vid2) | 0b11_00_00;

        PolygonalFeature {
            vertices,
            vids: PackedFeatureId::vertices([vid1, vid2]),
            fid: PackedFeatureId::face(fid),
            num_vertices: 2,
        }
    }

    /// Return the face of this cuboid with a normal that maximizes
    /// the dot product with `local_dir`.
    #[cfg(feature = "dim3")]
    pub fn support_feature(&self, local_dir: Vector) -> PolygonalFeature {
        // TODO: this should actually return the feature.
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
    // pub(crate) fn support_vertex(&self, local_dir: Vector) -> CuboidFeatureVertex {
    //     let vertex = local_support_point(self, local_dir);
    //     let vid = vertex_feature_id(vertex);
    //
    //     CuboidFeatureVertex { vertex, vid }
    // }

    /// Return the edge segment of this cuboid with a normal cone containing
    /// a direction that that maximizes the dot product with `local_dir`.
    #[cfg(feature = "dim3")]
    pub fn local_support_edge_segment(&self, local_dir: Vector) -> Segment {
        let he = self.half_extents;
        let i = local_dir.abs().min_position();
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        let mut a = Vector::ZERO;
        a.vset(i, he.vget(i));
        a.vset(j, local_dir.vget(j).copy_sign_to(he.vget(j)));
        a.vset(k, local_dir.vget(k).copy_sign_to(he.vget(k)));

        let mut b = a;
        b.vset(i, -he.vget(i));

        Segment::new(a, b)
    }

    /// Computes the face with a normal that maximizes the dot-product with `local_dir`.
    #[cfg(feature = "dim3")]
    // The identity ors/shifts below are kept so the bit patterns line up with the comments
    // documenting the vertex/edge numbering.
    #[allow(clippy::identity_op)]
    pub fn support_face(&self, local_dir: Vector) -> PolygonalFeature {
        // NOTE: can we use the orthonormal basis of local_dir
        // to make this AoSoA friendly?
        let he = self.half_extents;
        let imax = local_dir.abs().max_position();
        #[expect(clippy::unnecessary_cast)]
        let sign = match imax {
            0 => local_dir.x.copy_sign_to(1.0 as Real),
            1 => local_dir.y.copy_sign_to(1.0 as Real),
            _ => local_dir.z.copy_sign_to(1.0 as Real),
        };

        let vertices = match imax {
            0 => [
                Vector::new(he.x * sign, he.y, he.z),
                Vector::new(he.x * sign, -he.y, he.z),
                Vector::new(he.x * sign, -he.y, -he.z),
                Vector::new(he.x * sign, he.y, -he.z),
            ],
            1 => [
                Vector::new(he.x, he.y * sign, he.z),
                Vector::new(-he.x, he.y * sign, he.z),
                Vector::new(-he.x, he.y * sign, -he.z),
                Vector::new(he.x, he.y * sign, -he.z),
            ],
            _ => [
                Vector::new(he.x, he.y, he.z * sign),
                Vector::new(he.x, -he.y, he.z * sign),
                Vector::new(-he.x, -he.y, he.z * sign),
                Vector::new(-he.x, he.y, he.z * sign),
            ],
        };

        pub fn vid(i: u32) -> u32 {
            // Each vertex has an even feature id.
            i * 2
        }

        let sign_index = ((sign as isize + 1) / 2) as u32;
        // The vertex id as numbered depending on the sign of the vertex
        // component. A + sign means the corresponding bit is 0 while a -
        // sign means the corresponding bit is 1.
        // For example the vertex [2.0, -1.0, -3.0] has the id 0b011
        let vids = match imax {
            0 => {
                let sbit = sign_index << 2;
                [
                    vid(0b000 | sbit),
                    vid(0b010 | sbit),
                    vid(0b011 | sbit),
                    vid(0b001 | sbit),
                ]
            }
            1 => {
                let sbit = sign_index << 1;
                [
                    vid(0b000 | sbit),
                    vid(0b100 | sbit),
                    vid(0b101 | sbit),
                    vid(0b001 | sbit),
                ]
            }
            _ => {
                let sbit = sign_index;
                [
                    vid(0b000 | sbit),
                    vid(0b010 | sbit),
                    vid(0b110 | sbit),
                    vid(0b100 | sbit),
                ]
            }
        };

        // The feature ids of edges is obtained from the vertex ids
        // of their endpoints.
        // Assuming vid1 > vid2, we do:   (vid1 << 3) | vid2 | 0b11000000
        //
        let eids = match imax {
            0 => {
                let sbits = (sign_index << 2) | (sign_index << 5); // 0b00_100_100
                [
                    0b11_010_000 | sbits,
                    0b11_011_010 | sbits,
                    0b11_011_001 | sbits,
                    0b11_001_000 | sbits,
                ]
            }
            1 => {
                let sbits = (sign_index << 1) | (sign_index << 4); // 0b00_010_010
                [
                    0b11_100_000 | sbits,
                    0b11_101_100 | sbits,
                    0b11_101_001 | sbits,
                    0b11_001_000 | sbits,
                ]
            }
            _ => {
                let sbits = (sign_index << 0) | (sign_index << 3); // 0b00_001_001
                [
                    0b11_010_000 | sbits,
                    0b11_110_010 | sbits,
                    0b11_110_100 | sbits,
                    0b11_100_000 | sbits,
                ]
            }
        };

        // The face with normals [x, y, z] are numbered [10, 11, 12].
        // The face with negated normals are numbered [13, 14, 15].
        let fid = imax as u32 + sign_index * 3 + 10;

        PolygonalFeature {
            vertices,
            vids: PackedFeatureId::vertices(vids),
            eids: PackedFeatureId::edges(eids),
            fid: PackedFeatureId::face(fid),
            num_vertices: 4,
        }
    }

    /// The normal of the given feature of this shape.
    #[cfg(feature = "dim2")]
    pub fn feature_normal(&self, feature: FeatureId) -> Option<Vector> {
        match feature {
            FeatureId::Face(id) => {
                let mut dir: Vector = Vector::ZERO;

                if id < 2 {
                    dir.vset(id as usize, 1.0);
                } else {
                    dir.vset(id as usize - 2, -1.0);
                }
                Some(dir)
            }
            FeatureId::Vertex(id) => {
                let mut dir: Vector = Vector::ZERO;

                match id {
                    0b00 => {
                        dir.x = 1.0;
                        dir.y = 1.0;
                    }
                    0b01 => {
                        dir.y = 1.0;
                        dir.x = -1.0;
                    }
                    0b11 => {
                        dir.x = -1.0;
                        dir.y = -1.0;
                    }
                    0b10 => {
                        dir.y = -1.0;
                        dir.x = 1.0;
                    }
                    _ => return None,
                }

                Some(dir.normalize())
            }
            _ => None,
        }
    }

    /// The normal of the given feature of this shape.
    #[cfg(feature = "dim3")]
    pub fn feature_normal(&self, feature: FeatureId) -> Option<Vector> {
        match feature {
            FeatureId::Face(id) => {
                let mut dir: Vector = Vector::ZERO;

                if id < 3 {
                    dir.vset(id as usize, 1.0);
                } else {
                    dir.vset(id as usize - 3, -1.0);
                }
                Some(dir)
            }
            FeatureId::Edge(id) => {
                let edge = id & 0b011;
                let face1 = (edge + 1) % 3;
                let face2 = (edge + 2) % 3;
                let signs = id >> 2;

                let mut dir: Vector = Vector::ZERO;

                if signs & (1 << face1) != 0 {
                    dir.vset(face1 as usize, -1.0)
                } else {
                    dir.vset(face1 as usize, 1.0)
                }

                if signs & (1 << face2) != 0 {
                    dir.vset(face2 as usize, -1.0)
                } else {
                    dir.vset(face2 as usize, 1.0);
                }

                Some(dir.normalize())
            }
            FeatureId::Vertex(id) => {
                let mut dir: Vector = Vector::ZERO;
                for i in 0..3 {
                    if id & (1 << i) != 0 {
                        dir.vset(i, -1.0);
                    } else {
                        dir.vset(i, 1.0)
                    }
                }

                Some(dir.normalize())
            }
            _ => None,
        }
    }
}

impl SupportMap for Cuboid {
    #[inline]
    fn local_support_point(&self, dir: Vector) -> Vector {
        dir.copy_sign_to(self.half_extents)
    }
}

/*
impl ConvexPolyhedron for Cuboid {
    fn vertex(&self, id: FeatureId) -> Vector {
        let vid = id.unwrap_vertex();
        let mut res = self.half_extents;

        for i in 0..DIM {
            if vid & (1 << i) != 0 {
                res.vset(i, -res.vget(i))
            }
        }

        res
    }

    #[cfg(feature = "dim3")]
    fn edge(&self, id: FeatureId) -> (Vector, Vector, FeatureId, FeatureId) {
        let eid = id.unwrap_edge();
        let mut res = self.half_extents;

        let edge_i = eid & 0b11;
        let vertex_i = eid >> 2;

        for i in 0..DIM {
            if i as u32 != edge_i && (vertex_i & (1 << i) != 0) {
                res.vset(i, -res.vget(i))
            }
        }

        let p1 = res;
        res.vset(edge_i as usize, -res.vget(edge_i as usize));
        let p2 = res;
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
            sign = 1.0;
        } else {
            i1 = i - DIM;
            sign = -1.0;
        }

        #[cfg(feature = "dim2")]
        {
            let i2 = (i1 + 1) % 2;

            let mut vertex = self.half_extents;
            vertex.vset(i1, vertex.vget(i1) * sign);
            vertex.vset(i2, vertex.vget(i2) * if i1 == 0 { -sign } else { sign });

            let p1 = vertex;
            vertex.vset(i2, -vertex.vget(i2));
            let p2 = vertex;

            let mut vertex_id1 = if sign < 0.0 {
                1 << i1
            } else {
                0
            };
            let mut vertex_id2 = vertex_id1;
            if p1.vget(i2) < 0.0 {
                vertex_id1 |= 1 << i2;
            } else {
                vertex_id2 |= 1 << i2;
            }

            out.push(p1, FeatureId::Vertex(vertex_id1));
            out.push(p2, FeatureId::Vertex(vertex_id2));

            let mut normal: Vector = Vector::ZERO;
            normal.vset(i1, sign);
            out.set_normal(normal);
            out.set_feature_id(FeatureId::Face(i as u32));
        }
        #[cfg(feature = "dim3")]
        {
            let i2 = (i1 + 1) % 3;
            let i3 = (i1 + 2) % 3;
            let (edge_i2, edge_i3) = if sign > 0.0 {
                (i2, i3)
            } else {
                (i3, i2)
            };
            let mask_i2 = !(1 << edge_i2); // The masks are for ensuring each edge has a unique ID.
            let mask_i3 = !(1 << edge_i3);
            let mut vertex = self.half_extents;
            vertex.vset(i1, vertex.vget(i1) * sign);

            let (sbit, msbit) = if sign < 0.0 {
                (1, 0)
            } else {
                (0, 1)
            };
            let mut vertex_id = sbit << i1;
            out.push(vertex, FeatureId::Vertex(vertex_id));
            out.push_edge_feature_id(FeatureId::Edge(
                edge_i2 as u32 | ((vertex_id & mask_i2) << 2),
            ));

            vertex.vset(i2, -sign * self.half_extents.vget(i2));
            vertex.vset(i3, sign * self.half_extents.vget(i3));
            vertex_id |= msbit << i2 | sbit << i3;
            out.push(vertex, FeatureId::Vertex(vertex_id));
            out.push_edge_feature_id(FeatureId::Edge(
                edge_i3 as u32 | ((vertex_id & mask_i3) << 2),
            ));

            vertex.vset(i2, -self.half_extents.vget(i2));
            vertex.vset(i3, -self.half_extents.vget(i3));
            vertex_id |= 1 << i2 | 1 << i3;
            out.push(vertex, FeatureId::Vertex(vertex_id));
            out.push_edge_feature_id(FeatureId::Edge(
                edge_i2 as u32 | ((vertex_id & mask_i2) << 2),
            ));

            vertex.vset(i2, sign * self.half_extents.vget(i2));
            vertex.vset(i3, -sign * self.half_extents.vget(i3));
            vertex_id = sbit << i1 | sbit << i2 | msbit << i3;
            out.push(vertex, FeatureId::Vertex(vertex_id));
            out.push_edge_feature_id(FeatureId::Edge(
                edge_i3 as u32 | ((vertex_id & mask_i3) << 2),
            ));

            let mut normal: Vector = Vector::ZERO;
            normal.vset(i1, sign);
            out.set_normal(normal);

            if sign > 0.0 {
                out.set_feature_id(FeatureId::Face(i1 as u32));
            } else {
                out.set_feature_id(FeatureId::Face(i1 as u32 + 3));
            }

            out.recompute_edge_normals();
        }
    }

    fn support_face_toward(
        &self,
        m: &Pose,
        dir: Vector,
        out: &mut ConvexPolygonalFeature,
    ) {
        out.clear();
        let local_dir = m.inverse_transform_vector(dir);
        let imax = iamax(local_dir);

        if local_dir.vget(imax) > 0.0 {
            self.face(FeatureId::Face(imax as u32), out);
            out.transform_by(m);
        } else {
            self.face(FeatureId::Face((imax + DIM) as u32), out);
            out.transform_by(m);
        }
    }

    fn support_feature_toward(
        &self,
        m: &Pose,
        dir: Vector,
        angle: Real,
        out: &mut ConvexPolygonalFeature,
    ) {
        let local_dir = m.inverse_transform_vector(dir);
        let cang = <Real as ComplexField>::cos(angle);
        let mut support_point = self.half_extents;

        out.clear();

        #[cfg(feature = "dim2")]
        {
            let mut support_point_id = 0;
            for i1 in 0..2 {
                let sign = local_dir.vget(i1).signum();
                if sign * local_dir.vget(i1) >= cang {
                    if sign > 0.0 {
                        self.face(FeatureId::Face(i1 as u32), out);
                        out.transform_by(m);
                    } else {
                        self.face(FeatureId::Face(i1 as u32 + 2), out);
                        out.transform_by(m);
                    }
                    return;
                } else {
                    if sign < 0.0 {
                        support_point_id |= 1 << i1;
                    }
                    support_point.vset(i1, support_point.vget(i1) * sign);
                }
            }

            // We are not on a face, return the support vertex.
            out.push(
                m * support_point,
                FeatureId::Vertex(support_point_id),
            );
            out.set_feature_id(FeatureId::Vertex(support_point_id));
        }

        #[cfg(feature = "dim3")]
        {
            let sang = <Real as ComplexField>::sin(angle);
            let mut support_point_id = 0;

            // Check faces.
            for i1 in 0..3 {
                let sign = local_dir.vget(i1).signum();
                if sign * local_dir.vget(i1) >= cang {
                    if sign > 0.0 {
                        self.face(FeatureId::Face(i1 as u32), out);
                        out.transform_by(m);
                    } else {
                        self.face(FeatureId::Face(i1 as u32 + 3), out);
                        out.transform_by(m);
                    }
                    return;
                } else {
                    if sign < 0.0 {
                        support_point.vset(i1, support_point.vget(i1) * sign);
                        support_point_id |= 1 << i1;
                    }
                }
            }

            // Check edges.
            for i in 0..3 {
                let sign = local_dir.vget(i).signum();

                // sign * local_dir.vget(i) <= cos(pi / 2 - angle)
                if sign * local_dir.vget(i) <= sang {
                    support_point.vset(i, -self.half_extents.vget(i));
                    let p1 = support_point;
                    support_point.vset(i, self.half_extents.vget(i));
                    let p2 = support_point;
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
                m * support_point,
                FeatureId::Vertex(support_point_id),
            );
            out.set_feature_id(FeatureId::Vertex(support_point_id));
        }
    }

    fn support_feature_id_toward(&self, local_dir: Vector) -> FeatureId {
        let one_degree: Real = (f64::consts::PI / 180.0) as Real;
        let cang = <Real as ComplexField>::cos(one_degree);

        #[cfg(feature = "dim2")]
        {
            let mut support_point_id = 0;
            for i1 in 0..2 {
                let sign = local_dir.vget(i1).signum();
                if sign * local_dir.vget(i1) >= cang {
                    if sign > 0.0 {
                        return FeatureId::Face(i1 as u32);
                    } else {
                        return FeatureId::Face(i1 as u32 + 2);
                    }
                } else {
                    if sign < 0.0 {
                        support_point_id |= 1 << i1;
                    }
                }
            }

            // We are not on a face, return the support vertex.
            FeatureId::Vertex(support_point_id)
        }

        #[cfg(feature = "dim3")]
        {
            let sang = <Real as ComplexField>::sin(one_degree);
            let mut support_point_id = 0;

            // Check faces.
            for i1 in 0..3 {
                let sign = local_dir.vget(i1).signum();
                if sign * local_dir.vget(i1) >= cang {
                    if sign > 0.0 {
                        return FeatureId::Face(i1 as u32);
                    } else {
                        return FeatureId::Face(i1 as u32 + 3);
                    }
                } else {
                    if sign < 0.0 {
                        support_point_id |= 1 << i1;
                    }
                }
            }

            // Check edges.
            for i in 0..3 {
                let sign = local_dir.vget(i).signum();

                // sign * local_dir.vget(i) <= cos(pi / 2 - angle)
                if sign * local_dir.vget(i) <= sang {
                    let mask_i = !(1 << i); // To ensure each edge has a unique id.
                    return FeatureId::Edge(i as u32 | ((support_point_id & mask_i) << 2));
                }
            }

            FeatureId::Vertex(support_point_id)
        }
    }
}
*/
