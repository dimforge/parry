use crate::bounding_volume::AABB;
use crate::math::{Isometry, Point, Real};
use crate::partitioning::QBVH;
use crate::query::{PointProjection, PointQueryWithLocation};
use crate::shape::composite_shape::SimdCompositeShape;
use crate::shape::{FeatureId, Segment, SegmentPointLocation, Shape, TypedSimdCompositeShape};

#[cfg(not(feature = "std"))]
use na::ComplexField; // for .abs()

#[derive(Clone)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
/// A polyline.
pub struct Polyline {
    qbvh: QBVH<u32>,
    vertices: Vec<Point<Real>>,
    indices: Vec<[u32; 2]>,
}

impl Polyline {
    /// Creates a new polyline from a vertex buffer and an index buffer.
    pub fn new(vertices: Vec<Point<Real>>, indices: Option<Vec<[u32; 2]>>) -> Self {
        let indices =
            indices.unwrap_or_else(|| (0..vertices.len() as u32 - 1).map(|i| [i, i + 1]).collect());
        let data = indices.iter().enumerate().map(|(i, idx)| {
            let aabb =
                Segment::new(vertices[idx[0] as usize], vertices[idx[1] as usize]).local_aabb();
            (i as u32, aabb)
        });

        let mut qbvh = QBVH::new();
        // NOTE: we apply no dilation factor because we won't
        // update this tree dynamically.
        qbvh.clear_and_rebuild(data, 0.0);

        Self {
            qbvh,
            vertices,
            indices,
        }
    }

    /// Compute the axis-aligned bounding box of this polyline.
    pub fn aabb(&self, pos: &Isometry<Real>) -> AABB {
        self.qbvh.root_aabb().transform_by(pos)
    }

    /// Gets the local axis-aligned bounding box of this polyline.
    pub fn local_aabb(&self) -> &AABB {
        &self.qbvh.root_aabb()
    }

    pub(crate) fn qbvh(&self) -> &QBVH<u32> {
        &self.qbvh
    }

    /// The number of segments forming this polyline.
    pub fn num_segments(&self) -> usize {
        self.indices.len()
    }

    /// An iterator through all the segments of this mesh.
    pub fn segments(&self) -> impl Iterator<Item = Segment> + '_ {
        self.indices.iter().map(move |ids| {
            Segment::new(
                self.vertices[ids[0] as usize],
                self.vertices[ids[1] as usize],
            )
        })
    }

    /// Get the `i`-th segment of this mesh.
    pub fn segment(&self, i: u32) -> Segment {
        let idx = self.indices[i as usize];
        Segment::new(
            self.vertices[idx[0] as usize],
            self.vertices[idx[1] as usize],
        )
    }

    /// Transforms  the feature-id of a segment to the feature-id of this polyline.
    pub fn segment_feature_to_polyline_feature(
        &self,
        segment: u32,
        _feature: FeatureId,
    ) -> FeatureId {
        // TODO: return a vertex feature when it makes sense.
        #[cfg(feature = "dim2")]
        return FeatureId::Face(segment);
        #[cfg(feature = "dim3")]
        return FeatureId::Edge(segment);
    }

    /// The vertex buffer of this mesh.
    pub fn vertices(&self) -> &[Point<Real>] {
        &self.vertices[..]
    }

    /// The index buffer of this mesh.
    pub fn indices(&self) -> &[[u32; 2]] {
        &self.indices
    }

    /// A flat view of the index buffer of this mesh.
    pub fn flat_indices(&self) -> &[u32] {
        unsafe {
            let len = self.indices.len() * 3;
            let data = self.indices.as_ptr() as *const u32;
            std::slice::from_raw_parts(data, len)
        }
    }

    /// Reverse the orientation of this polyline by swapping the indices of all
    /// its segments and reverting its index buffer.
    pub fn reverse(&mut self) {
        for idx in &mut self.indices {
            idx.swap(0, 1);
        }

        self.indices.reverse();

        // Because we reversed the indices, we need to
        // adjust the segment indices stored in the QBVH.
        for (_, seg_id) in self.qbvh.iter_data_mut() {
            *seg_id = self.indices.len() as u32 - *seg_id - 1;
        }
    }

    /// Perform a point projection assuming a solid interior based on a counter-clock-wise orientation.
    ///
    /// This is similar to `self.project_local_point_and_get_location` except that the resulting
    /// `PointProjection::is_inside` will be set to true if the point is inside of the area delimited
    /// by this polyline, assuming that:
    /// - This polyline isn’t self-crossing.
    /// - This polyline is closed with `self.indices[i][1] == self.indices[(i + 1) % num_indices][0]` where
    ///   `num_indices == self.indices.len()`.
    /// - This polyline is oriented counter-clockwise.
    /// - In 3D, the polyline is assumed to be fully coplanar, on a plane with normal given by
    /// These properties are not checked.
    pub fn project_local_point_assuming_solid_interior_ccw(
        &self,
        point: Point<Real>,
        #[cfg(feature = "dim3")] axis: u8,
    ) -> (PointProjection, (u32, SegmentPointLocation)) {
        let mut proj = self.project_local_point_and_get_location(&point, false);
        let segment1 = self.segment((proj.1).0);

        #[cfg(feature = "dim2")]
        let normal1 = segment1.normal();
        #[cfg(feature = "dim3")]
        let normal1 = segment1.planar_normal(axis);

        if let Some(normal1) = normal1 {
            proj.0.is_inside = match proj.1 .1 {
                SegmentPointLocation::OnVertex(i) => {
                    let dir2 = if i == 0 {
                        let adj_seg = if proj.1 .0 == 0 {
                            self.indices().len() as u32 - 1
                        } else {
                            proj.1 .0 - 1
                        };

                        assert_eq!(segment1.a, self.segment(adj_seg).b);
                        -self.segment(adj_seg).scaled_direction()
                    } else {
                        assert_eq!(i, 1);
                        let adj_seg = (proj.1 .0 + 1) % self.indices().len() as u32;
                        assert_eq!(segment1.b, self.segment(adj_seg).a);

                        self.segment(adj_seg).scaled_direction()
                    };

                    let dot = normal1.dot(&dir2);
                    // TODO: is this threshold too big? This corresponds to an angle equal to
                    //       abs(acos(1.0e-3)) = (90 - 0.057) degrees.
                    //       We did encounter some cases where this was needed, but perhaps the
                    //       actual problem was an issue with the SegmentPointLocation (which should
                    //       perhaps have been Edge instead of Vertex)?
                    let threshold = 1.0e-3 * dir2.norm();
                    if dot.abs() > threshold {
                        // If the vertex is a reentrant vertex, then the point is
                        // inside. Otherwise, it is outside.
                        dot >= 0.0
                    } else {
                        // If the two edges are collinear, we can’t classify the vertex.
                        // So check against the edge’s normal instead.
                        (point - proj.0.point).dot(&normal1) <= 0.0
                    }
                }
                SegmentPointLocation::OnEdge(_) => (point - proj.0.point).dot(&normal1) <= 0.0,
            };
        }

        proj
    }
}

impl SimdCompositeShape for Polyline {
    fn map_part_at(&self, i: u32, f: &mut dyn FnMut(Option<&Isometry<Real>>, &dyn Shape)) {
        let tri = self.segment(i);
        f(None, &tri)
    }

    fn qbvh(&self) -> &QBVH<u32> {
        &self.qbvh
    }
}

impl TypedSimdCompositeShape for Polyline {
    type PartShape = Segment;
    type PartId = u32;

    #[inline(always)]
    fn map_typed_part_at(
        &self,
        i: u32,
        mut f: impl FnMut(Option<&Isometry<Real>>, &Self::PartShape),
    ) {
        let seg = self.segment(i);
        f(None, &seg)
    }

    #[inline(always)]
    fn map_untyped_part_at(&self, i: u32, mut f: impl FnMut(Option<&Isometry<Real>>, &dyn Shape)) {
        let seg = self.segment(i);
        f(None, &seg)
    }

    fn typed_qbvh(&self) -> &QBVH<u32> {
        &self.qbvh
    }
}
