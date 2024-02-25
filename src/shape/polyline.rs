use crate::bounding_volume::Aabb;
use crate::math::{Isometry, Point, Real, Vector};
use crate::partitioning::Qbvh;
use crate::query::{PointProjection, PointQueryWithLocation};
use crate::shape::composite_shape::SimdCompositeShape;
use crate::shape::{FeatureId, Segment, SegmentPointLocation, Shape, TypedSimdCompositeShape};

use crate::utils::DefaultStorage;
#[cfg(not(feature = "std"))]
use na::ComplexField; // for .abs()

#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
    archive(check_bytes)
)]
/// A polyline.
pub struct Polyline {
    qbvh: Qbvh<u32>,
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

        let mut qbvh = Qbvh::new();
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
    pub fn aabb(&self, pos: &Isometry<Real>) -> Aabb {
        self.qbvh.root_aabb().transform_by(pos)
    }

    /// Gets the local axis-aligned bounding box of this polyline.
    pub fn local_aabb(&self) -> &Aabb {
        self.qbvh.root_aabb()
    }

    pub(crate) fn qbvh(&self) -> &Qbvh<u32> {
        &self.qbvh
    }

    /// The number of segments forming this polyline.
    pub fn num_segments(&self) -> usize {
        self.indices.len()
    }

    /// An iterator through all the segments of this mesh.
    pub fn segments(&self) -> impl ExactSizeIterator<Item = Segment> + '_ {
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
            let len = self.indices.len() * 2;
            let data = self.indices.as_ptr() as *const u32;
            std::slice::from_raw_parts(data, len)
        }
    }

    /// Computes a scaled version of this polyline.
    pub fn scaled(mut self, scale: &Vector<Real>) -> Self {
        self.vertices
            .iter_mut()
            .for_each(|pt| pt.coords.component_mul_assign(scale));
        Self {
            qbvh: self.qbvh.scaled(scale),
            vertices: self.vertices,
            indices: self.indices,
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
        // adjust the segment indices stored in the Qbvh.
        for (_, seg_id) in self.qbvh.iter_data_mut() {
            *seg_id = self.indices.len() as u32 - *seg_id - 1;
        }
    }

    /// Extracts the connected components of this polyline, consuming `self`.
    ///
    /// This method is currently quite restrictive on the kind of allowed input. The polyline
    /// represented by `self` must already have an index buffer sorted such that:
    /// - Each connected component appears in the index buffer one after the other, i.e., a
    ///   connected component of this polyline must be a contiguous range of this polyline’s
    ///   index buffer.
    /// - Each connected component is closed, i.e., each range of this polyline index buffer
    ///   `self.indices[i_start..=i_end]` forming a complete connected component, we must have
    ///   `self.indices[i_start][0] == self.indices[i_end][1]`.
    /// - The indices for each component must already be in order, i.e., if the segments
    ///   `self.indices[i]` and `self.indices[i + 1]` are part of the same connected component then
    ///   we must have `self.indices[i][1] == self.indices[i + 1][0]`.
    ///
    /// # Output
    /// Returns the set of polylines. If the inputs fulfill the constraints mentioned above, each
    /// polyline will be a closed loop with consistent edge orientations, i.e., for all indices `i`,
    /// we have `polyline.indices[i][1] == polyline.indices[i + 1][0]`.
    ///
    /// The orientation of each closed loop (clockwise or counterclockwise) are identical to their
    /// original orientation in `self`.
    pub fn extract_connected_components(&self) -> Vec<Polyline> {
        let vertices = self.vertices();
        let indices = self.indices();

        if indices.is_empty() {
            // Polyline is empty, return empty Vec
            Vec::new()
        } else {
            let mut components = Vec::new();

            let mut start_i = 0; // Start position of component
            let mut start_node = indices[0][0]; // Start vertex index of component

            let mut component_vertices = Vec::new();
            let mut component_indices: Vec<[u32; 2]> = Vec::new();

            // Iterate over indices, building polylines as we go
            for (i, idx) in indices.iter().enumerate() {
                component_vertices.push(vertices[idx[0] as usize]);

                if idx[1] != start_node {
                    // Keep scanning and adding data
                    component_indices.push([(i - start_i) as u32, (i - start_i + 1) as u32]);
                } else {
                    // Start node reached: build polyline and start next component
                    component_indices.push([(i - start_i) as u32, 0]);
                    components.push(Polyline::new(
                        std::mem::take(&mut component_vertices),
                        Some(std::mem::take(&mut component_indices)),
                    ));

                    if i + 1 < indices.len() {
                        // More components to find
                        start_node = indices[i + 1][0];
                        start_i = i + 1;
                    }
                }
            }

            components
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

    fn qbvh(&self) -> &Qbvh<u32> {
        &self.qbvh
    }
}

impl TypedSimdCompositeShape for Polyline {
    type PartShape = Segment;
    type PartId = u32;
    type QbvhStorage = DefaultStorage;

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

    fn typed_qbvh(&self) -> &Qbvh<u32> {
        &self.qbvh
    }
}
