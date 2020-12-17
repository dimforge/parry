use crate::math::{Isometry, Point, Real, Vector};
use crate::query::{self, ContactManifold, TrackedContact};
use crate::shape::{Cuboid, Segment};

#[derive(Debug)]
#[allow(dead_code)]
pub enum CuboidFeature {
    Face(CuboidFeatureFace),
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
pub struct CuboidFeatureFace {
    pub vertices: [Point<Real>; 2],
    pub vids: [u8; 2],
    pub fid: u8,
}

impl From<Segment> for CuboidFeatureFace {
    fn from(seg: Segment) -> Self {
        CuboidFeatureFace {
            vertices: [seg.a, seg.b],
            vids: [0, 2],
            fid: 1,
        }
    }
}

impl CuboidFeatureFace {
    pub fn transform_by(&mut self, iso: &Isometry<Real>) {
        self.vertices[0] = iso * self.vertices[0];
        self.vertices[1] = iso * self.vertices[1];
    }
}

impl CuboidFeature {
    pub fn transform_by(&mut self, iso: &Isometry<Real>) {
        match self {
            CuboidFeature::Face(face) => face.transform_by(iso),
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
        let tangent1 = face1.vertices[1] - face1.vertices[0];
        let normal1 = Vector::new(-tangent1.y, tangent1.x);
        let denom = -normal1.dot(&sep_axis1);
        let dist = (face1.vertices[0] - v2_1).dot(&normal1) / denom;
        let local_p2 = v2_1;
        let local_p1 = v2_1 - dist * normal1;

        let contact = TrackedContact::flipped(
            local_p1,
            pos12.inverse_transform_point(&local_p2),
            face1.fid,
            vertex2.vid,
            dist,
            flipped,
        );

        manifold.points.push(contact);
    }

    pub fn face_face_contacts<ManifoldData, ContactData: Default + Copy>(
        pos12: &Isometry<Real>,
        face1: &CuboidFeatureFace,
        normal1: &Vector<Real>,
        face2: &CuboidFeatureFace,
        _prediction_distance: Real,
        manifold: &mut ContactManifold<ManifoldData, ContactData>,
        flipped: bool,
    ) {
        if let Some((clip_a, clip_b)) = query::details::clip_segment_segment(
            (face1.vertices[0], face1.vertices[1]),
            (pos12 * face2.vertices[0], pos12 * face2.vertices[1]),
        ) {
            let fids1 = [face1.vids[0], face1.fid, face1.vids[1]];
            let fids2 = [face2.vids[0], face2.fid, face2.vids[1]];

            let dist = (clip_a.1 - clip_a.0).dot(normal1);
            if true {
                // dist < prediction_distance {
                let contact = TrackedContact::flipped(
                    clip_a.0,
                    pos12.inverse_transform_point(&clip_a.1),
                    fids1[clip_a.2],
                    fids2[clip_a.3],
                    dist,
                    flipped,
                );
                manifold.points.push(contact);
            }

            let dist = (clip_b.1 - clip_b.0).dot(normal1);
            if true {
                // dist < prediction_distance {
                let contact = TrackedContact::flipped(
                    clip_b.0,
                    pos12.inverse_transform_point(&clip_b.1),
                    fids1[clip_b.2],
                    fids2[clip_b.3],
                    dist,
                    flipped,
                );
                manifold.points.push(contact);
            }
        }
    }
}

impl Cuboid {
    // pub fn polygon_ref(
    //     cuboid: Cuboid,
    //     out_vertices: &mut [Point<Real>; 4],
    //     out_normals: &mut [Vector<Real>; 4],
    // ) -> PolygonRef {
    //     *out_vertices = [
    //         Point::new(cuboid.half_extents.x, -cuboid.half_extents.y),
    //         Point::new(cuboid.half_extents.x, cuboid.half_extents.y),
    //         Point::new(-cuboid.half_extents.x, cuboid.half_extents.y),
    //         Point::new(-cuboid.half_extents.x, -cuboid.half_extents.y),
    //     ];
    //     *out_normals = [Vector::x(), Vector::y(), -Vector::x(), -Vector::y()];
    //
    //     PolygonRef {
    //         vertices: &out_vertices[..],
    //         normals: &out_normals[..],
    //     }
    // }

    pub fn vertex_feature_id(vertex: Point<Real>) -> u8 {
        ((vertex.x.to_bits() >> 31) & 0b001 | (vertex.y.to_bits() >> 30) & 0b010) as u8
    }

    pub fn support_feature(&self, local_dir: Vector<Real>) -> CuboidFeature {
        // In 2D, it is best for stability to always return a face.
        // It won't have any notable impact on performances anyway.
        CuboidFeature::Face(self.support_face(local_dir))

        /*
        let amax = local_dir.amax();

        const MAX_DOT_THRESHOLD: Real = 0.98480775301; // 10 degrees.

        if amax > MAX_DOT_THRESHOLD {
            // Support face.
            CuboidFeature::Face(self.support_face(local_dir))
        } else {
            // Support vertex
            CuboidFeature::Vertex(self.support_vertex(local_dir))
        }
        */
    }

    pub fn support_face(&self, local_dir: Vector<Real>) -> CuboidFeatureFace {
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

        CuboidFeatureFace {
            vertices: [a, b],
            vids: [vid1, vid2],
            fid,
        }
    }
}
