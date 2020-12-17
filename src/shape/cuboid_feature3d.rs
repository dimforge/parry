use crate::approx::AbsDiffEq;
use crate::math::{Isometry, Point, Real, Vector};
use crate::query::{ContactManifold, TrackedContact};
use crate::shape::{Cuboid, PolygonalFeature};
use crate::utils::WBasis;
use na::Point2;

#[derive(Debug)]
#[allow(dead_code)]
pub enum CuboidFeature {
    Face(CuboidFeatureFace),
    Edge(CuboidFeatureEdge),
    Vertex(CuboidFeatureVertex),
}

#[derive(Debug)]
pub struct CuboidFeatureVertex {
    pub vertex: Point<Real>,
    pub vid: u8,
}

impl CuboidFeatureVertex {
    pub fn transform_by(&mut self, iso: &Isometry<Real>) {
        self.vertex = iso * self.vertex;
    }
}

#[derive(Debug)]
pub struct CuboidFeatureEdge {
    pub vertices: [Point<Real>; 2],
    pub vids: [u8; 2],
    pub eid: u8,
}

impl CuboidFeatureEdge {
    pub fn transform_by(&mut self, iso: &Isometry<Real>) {
        self.vertices[0] = iso * self.vertices[0];
        self.vertices[1] = iso * self.vertices[1];
    }
}

#[derive(Debug)]
pub struct CuboidFeatureFace {
    pub vertices: [Point<Real>; 4],
    pub vids: [u8; 4], // Feature ID of the vertices.
    pub eids: [u8; 4], // Feature ID of the edges.
    pub fid: u8,       // Feature ID of the face.
}

impl CuboidFeatureFace {
    pub fn transform_by(&mut self, iso: &Isometry<Real>) {
        self.vertices[0] = iso * self.vertices[0];
        self.vertices[1] = iso * self.vertices[1];
        self.vertices[2] = iso * self.vertices[2];
        self.vertices[3] = iso * self.vertices[3];
    }
}

impl CuboidFeature {
    pub fn transform_by(&mut self, iso: &Isometry<Real>) {
        match self {
            CuboidFeature::Face(face) => face.transform_by(iso),
            CuboidFeature::Edge(edge) => edge.transform_by(iso),
            CuboidFeature::Vertex(vertex) => vertex.transform_by(iso),
        }
    }

    /// Compute contacts points between a face and a vertex.
    ///
    /// This method assume we already know that at least one contact exists.
    pub fn face_vertex_contacts<ManifoldData, ContactData: Default + Copy>(
        pos12: &Isometry<Real>,
        face1: &CuboidFeatureFace,
        sep_axis1: &Vector<Real>,
        vertex2: &CuboidFeatureVertex,
        manifold: &mut ContactManifold<ManifoldData, ContactData>,
        flipped: bool,
    ) {
        let v2_1 = pos12 * vertex2.vertex;
        let normal1 =
            (face1.vertices[0] - face1.vertices[1]).cross(&(face1.vertices[2] - face1.vertices[1]));
        let denom = -normal1.dot(&sep_axis1);
        let dist = (face1.vertices[0] - v2_1).dot(&normal1) / denom;
        let local_p2_1 = v2_1;
        let local_p1 = v2_1 - dist * sep_axis1;

        let contact = TrackedContact::flipped(
            local_p1,
            pos12.inverse_transform_point(&local_p2_1),
            face1.fid,
            vertex2.vid,
            dist,
            flipped,
        );
        manifold.points.push(contact);
    }

    /// Compute contacts points between a face and an edge.
    ///
    /// This method assume we already know that at least one contact exists.
    pub fn face_edge_contacts<ManifoldData, ContactData: Default + Copy>(
        pos12: &Isometry<Real>,
        face1: &CuboidFeatureFace,
        sep_axis1: &Vector<Real>,
        edge2: &CuboidFeatureEdge,
        prediction_distance: Real,
        manifold: &mut ContactManifold<ManifoldData, ContactData>,
        flipped: bool,
    ) {
        // Project the faces to a 2D plane for contact clipping.
        // The plane they are projected onto has normal sep_axis1
        // and contains the origin (this is numerically OK because
        // we are not working in world-space here).
        let basis = sep_axis1.orthonormal_basis();
        let projected_face1 = [
            Point2::new(
                face1.vertices[0].coords.dot(&basis[0]),
                face1.vertices[0].coords.dot(&basis[1]),
            ),
            Point2::new(
                face1.vertices[1].coords.dot(&basis[0]),
                face1.vertices[1].coords.dot(&basis[1]),
            ),
            Point2::new(
                face1.vertices[2].coords.dot(&basis[0]),
                face1.vertices[2].coords.dot(&basis[1]),
            ),
            Point2::new(
                face1.vertices[3].coords.dot(&basis[0]),
                face1.vertices[3].coords.dot(&basis[1]),
            ),
        ];

        let vertices2_1 = [pos12 * edge2.vertices[0], pos12 * edge2.vertices[1]];

        let projected_edge2 = [
            Point2::new(
                vertices2_1[0].coords.dot(&basis[0]),
                vertices2_1[0].coords.dot(&basis[1]),
            ),
            Point2::new(
                vertices2_1[1].coords.dot(&basis[0]),
                vertices2_1[1].coords.dot(&basis[1]),
            ),
        ];

        // Now we have to compute the intersection between all pairs of
        // edges from the face 1 with the edge 2
        for i in 0..4 {
            let projected_edge1 = [projected_face1[i], projected_face1[(i + 1) % 4]];

            if let Some(bcoords) = closest_points_line2d(projected_edge1, projected_edge2) {
                if bcoords.0 > 0.0 && bcoords.0 < 1.0 && bcoords.1 > 0.0 && bcoords.1 < 1.0 {
                    // Found a contact between the two edges.
                    let edge1 = [face1.vertices[i], face1.vertices[(i + 1) % 4]];
                    let local_p1 = edge1[0] * (1.0 - bcoords.0) + edge1[1].coords * bcoords.0;
                    let local_p2_1 =
                        vertices2_1[0] * (1.0 - bcoords.1) + vertices2_1[1].coords * bcoords.1;
                    let dist = (local_p2_1 - local_p1).dot(&sep_axis1);

                    if dist < prediction_distance {
                        let contact = TrackedContact::flipped(
                            local_p1,
                            pos12.inverse_transform_point(&local_p2_1),
                            face1.eids[i],
                            edge2.eid,
                            dist,
                            flipped,
                        );
                        manifold.points.push(contact);
                    }
                }
            }
        }

        // Project the two points from the edge into the face.
        let normal1 =
            (face1.vertices[2] - face1.vertices[1]).cross(&(face1.vertices[0] - face1.vertices[1]));

        'point_loop2: for i in 0..2 {
            let p2 = projected_edge2[i];

            let sign = (projected_face1[0] - projected_face1[3]).perp(&(p2 - projected_face1[3]));
            for j in 0..3 {
                let new_sign =
                    (projected_face1[j + 1] - projected_face1[j]).perp(&(p2 - projected_face1[j]));
                if new_sign * sign < 0.0 {
                    // The point lies outside.
                    continue 'point_loop2;
                }
            }

            // All the perp had the same sign: the point is inside of the other shapes projection.
            // Output the contact.
            let denom = -normal1.dot(&sep_axis1);
            let dist = (face1.vertices[0] - vertices2_1[i]).dot(&normal1) / denom;
            let local_p2_1 = vertices2_1[i];
            let local_p1 = vertices2_1[i] - dist * sep_axis1;

            if dist < prediction_distance {
                let contact = TrackedContact::flipped(
                    local_p1,
                    pos12.inverse_transform_point(&local_p2_1),
                    face1.fid,
                    edge2.vids[i],
                    dist,
                    flipped,
                );
                manifold.points.push(contact);
            }
        }
    }

    /// Compute contacts points between two edges.
    ///
    /// This method assume we already know that at least one contact exists.
    pub fn edge_edge_contacts<ManifoldData, ContactData: Default + Copy>(
        pos12: &Isometry<Real>,
        edge1: &CuboidFeatureEdge,
        sep_axis1: &Vector<Real>,
        edge2: &CuboidFeatureEdge,
        manifold: &mut ContactManifold<ManifoldData, ContactData>,
        flipped: bool,
    ) {
        let basis = sep_axis1.orthonormal_basis();
        let projected_edge1 = [
            Point2::new(
                edge1.vertices[0].coords.dot(&basis[0]),
                edge1.vertices[0].coords.dot(&basis[1]),
            ),
            Point2::new(
                edge1.vertices[1].coords.dot(&basis[0]),
                edge1.vertices[1].coords.dot(&basis[1]),
            ),
        ];

        let vertices2_1 = [pos12 * edge2.vertices[0], pos12 * edge2.vertices[1]];
        let projected_edge2 = [
            Point2::new(
                vertices2_1[0].coords.dot(&basis[0]),
                vertices2_1[0].coords.dot(&basis[1]),
            ),
            Point2::new(
                vertices2_1[1].coords.dot(&basis[0]),
                vertices2_1[1].coords.dot(&basis[1]),
            ),
        ];

        if let Some(bcoords) = closest_points_line2d(projected_edge1, projected_edge2) {
            let local_p1 =
                edge1.vertices[0] * (1.0 - bcoords.0) + edge1.vertices[1].coords * bcoords.0;
            let local_p2_1 = vertices2_1[0] * (1.0 - bcoords.1) + vertices2_1[1].coords * bcoords.1;
            let dist = (local_p2_1 - local_p1).dot(&sep_axis1);

            manifold.points.push(TrackedContact::flipped(
                local_p1,
                pos12.inverse_transform_point(&local_p2_1),
                edge1.eid,
                edge2.eid,
                dist,
                flipped,
            ));
        }
    }

    pub fn face_face_contacts<ManifoldData, ContactData: Default + Copy>(
        pos12: &Isometry<Real>,
        face1: &CuboidFeatureFace,
        sep_axis1: &Vector<Real>,
        face2: &CuboidFeatureFace,
        _prediction_distance: Real,
        manifold: &mut ContactManifold<ManifoldData, ContactData>,
        flipped: bool,
    ) {
        // Project the faces to a 2D plane for contact clipping.
        // The plane they are projected onto has normal sep_axis1
        // and contains the origin (this is numerically OK because
        // we are not working in world-space here).
        let basis = sep_axis1.orthonormal_basis();
        let projected_face1 = [
            Point2::new(
                face1.vertices[0].coords.dot(&basis[0]),
                face1.vertices[0].coords.dot(&basis[1]),
            ),
            Point2::new(
                face1.vertices[1].coords.dot(&basis[0]),
                face1.vertices[1].coords.dot(&basis[1]),
            ),
            Point2::new(
                face1.vertices[2].coords.dot(&basis[0]),
                face1.vertices[2].coords.dot(&basis[1]),
            ),
            Point2::new(
                face1.vertices[3].coords.dot(&basis[0]),
                face1.vertices[3].coords.dot(&basis[1]),
            ),
        ];

        let vertices2_1 = [
            pos12 * face2.vertices[0],
            pos12 * face2.vertices[1],
            pos12 * face2.vertices[2],
            pos12 * face2.vertices[3],
        ];
        let projected_face2 = [
            Point2::new(
                vertices2_1[0].coords.dot(&basis[0]),
                vertices2_1[0].coords.dot(&basis[1]),
            ),
            Point2::new(
                vertices2_1[1].coords.dot(&basis[0]),
                vertices2_1[1].coords.dot(&basis[1]),
            ),
            Point2::new(
                vertices2_1[2].coords.dot(&basis[0]),
                vertices2_1[2].coords.dot(&basis[1]),
            ),
            Point2::new(
                vertices2_1[3].coords.dot(&basis[0]),
                vertices2_1[3].coords.dot(&basis[1]),
            ),
        ];

        // Also find all the vertices located inside of the other projected face.
        let normal1 =
            (face1.vertices[2] - face1.vertices[1]).cross(&(face1.vertices[0] - face1.vertices[1]));
        let normal2_1 = (vertices2_1[2] - vertices2_1[1]).cross(&(vertices2_1[0] - vertices2_1[1]));

        // NOTE: The loop iterating on all the vertices for face1 is different from
        // the one iterating on all the vertices of face2. In the second loop, we
        // classify every point wrt. every edge on the other face. This will give
        // us bit masks to filter out several edge-edge intersections.
        'point_loop1: for i in 0..4 {
            let p1 = projected_face1[i];

            let sign = (projected_face2[0] - projected_face2[3]).perp(&(p1 - projected_face2[3]));
            for j in 0..3 {
                let new_sign =
                    (projected_face2[j + 1] - projected_face2[j]).perp(&(p1 - projected_face2[j]));
                if new_sign * sign < 0.0 {
                    // The point lies outside.
                    continue 'point_loop1;
                }
            }

            // All the perp had the same sign: the point is inside of the other shapes projection.
            // Output the contact.
            let denom = normal2_1.dot(&sep_axis1);
            let dist = (vertices2_1[0] - face1.vertices[i]).dot(&normal2_1) / denom;
            let local_p1 = face1.vertices[i];
            let local_p2_1 = face1.vertices[i] + dist * sep_axis1;

            if true {
                // dist <= prediction_distance {
                manifold.points.push(TrackedContact::flipped(
                    local_p1,
                    pos12.inverse_transform_point(&local_p2_1),
                    face1.vids[i],
                    face2.fid,
                    dist,
                    flipped,
                ));
            }
        }

        let is_clockwise1 = (projected_face1[0] - projected_face1[1])
            .perp(&(projected_face1[2] - projected_face1[1]))
            >= 0.0;
        let mut vertex_class2 = [0u8; 4];

        for i in 0..4 {
            let p2 = projected_face2[i];

            let sign = (projected_face1[0] - projected_face1[3]).perp(&(p2 - projected_face1[3]));
            vertex_class2[i] |= ((sign >= 0.0) as u8) << 3;

            for j in 0..3 {
                let sign =
                    (projected_face1[j + 1] - projected_face1[j]).perp(&(p2 - projected_face1[j]));
                vertex_class2[i] |= ((sign >= 0.0) as u8) << j;
            }

            if !is_clockwise1 {
                vertex_class2[i] = (!vertex_class2[i]) & 0b01111;
            }

            if vertex_class2[i] == 0 {
                // All the perp had the same sign: the point is inside of the other shapes projection.
                // Output the contact.
                let denom = -normal1.dot(&sep_axis1);
                let dist = (face1.vertices[0] - vertices2_1[i]).dot(&normal1) / denom;
                let local_p2_1 = vertices2_1[i];
                let local_p1 = vertices2_1[i] - dist * sep_axis1;

                if true {
                    // dist < prediction_distance {
                    manifold.points.push(TrackedContact::flipped(
                        local_p1,
                        pos12.inverse_transform_point(&local_p2_1),
                        face1.fid,
                        face2.vids[i],
                        dist,
                        flipped,
                    ));
                }
            }
        }

        // Now we have to compute the intersection between all pairs of
        // edges from the face 1 and from the face2.
        for j in 0..4 {
            let projected_edge2 = [projected_face2[j], projected_face2[(j + 1) % 4]];

            if (vertex_class2[j] & vertex_class2[(j + 1) % 4]) != 0 {
                continue;
            }

            let edge_class2 = vertex_class2[j] | vertex_class2[(j + 1) % 4];

            for i in 0..4 {
                if (edge_class2 & (1 << i)) != 0 {
                    let projected_edge1 = [projected_face1[i], projected_face1[(i + 1) % 4]];
                    if let Some(bcoords) = closest_points_line2d(projected_edge1, projected_edge2) {
                        if bcoords.0 > 0.0 && bcoords.0 < 1.0 && bcoords.1 > 0.0 && bcoords.1 < 1.0
                        {
                            // Found a contact between the two edges.
                            let edge1 = (face1.vertices[i], face1.vertices[(i + 1) % 4]);
                            let edge2 = (vertices2_1[j], vertices2_1[(j + 1) % 4]);
                            let local_p1 = edge1.0 * (1.0 - bcoords.0) + edge1.1.coords * bcoords.0;
                            let local_p2_1 =
                                edge2.0 * (1.0 - bcoords.1) + edge2.1.coords * bcoords.1;
                            let dist = (local_p2_1 - local_p1).dot(&sep_axis1);

                            if true {
                                // dist <= prediction_distance {
                                manifold.points.push(TrackedContact::flipped(
                                    local_p1,
                                    pos12.inverse_transform_point(&local_p2_1),
                                    face1.eids[i],
                                    face2.eids[j],
                                    dist,
                                    flipped,
                                ));
                            }
                        }
                    }
                }
            }
        }
    }
}

/// Compute the barycentric coordinates of the intersection between the two given lines.
/// Returns `None` if the lines are parallel.
fn closest_points_line2d(
    edge1: [Point2<Real>; 2],
    edge2: [Point2<Real>; 2],
) -> Option<(Real, Real)> {
    let inter = crate::query::details::closest_points_line_line_parameters_eps(
        &edge1[0],
        &(edge1[1] - edge1[0]),
        &edge2[0],
        &(edge2[1] - edge2[0]),
        Real::default_epsilon(),
    );

    if !inter.2 {
        Some((inter.0, inter.1))
    } else {
        None
    }
}

impl Cuboid {
    // pub fn vertex_feature_id(vertex: Point<Real>) -> u8 {
    //     ((vertex.x.to_bits() >> 31) & 0b001
    //         | (vertex.y.to_bits() >> 30) & 0b010
    //         | (vertex.z.to_bits() >> 29) & 0b100) as u8
    // }

    pub fn polyhedron_support_face(&self, local_dir: Vector<Real>) -> PolygonalFeature {
        self.support_face(local_dir).into()
    }

    pub fn support_feature(&self, local_dir: Vector<Real>) -> CuboidFeature {
        // FIXME: this should actually return the feature.
        // And we should change all the callers of this method to use
        // `.support_face` instead of this method to preserve their old behavior.
        CuboidFeature::Face(self.support_face(local_dir))
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

    // pub(crate) fn support_vertex(&self, local_dir: Vector<Real>) -> CuboidFeatureVertex {
    //     let vertex = local_support_point(self, local_dir);
    //     let vid = vertex_feature_id(vertex);
    //
    //     CuboidFeatureVertex { vertex, vid }
    // }

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

    pub fn support_face(&self, local_dir: Vector<Real>) -> CuboidFeatureFace {
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

        CuboidFeatureFace {
            vertices,
            vids,
            eids,
            fid: fid as u8,
        }
    }
}
