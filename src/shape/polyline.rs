use crate::bounding_volume::AABB;
use crate::math::{Isometry, Point, Real};
use crate::partitioning::SimdQuadTree;
use crate::shape::composite_shape::SimdCompositeShape;
use crate::shape::{FeatureId, Segment, Shape, TypedSimdCompositeShape};
use na::Point2;

#[derive(Clone)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
/// A polyline.
pub struct Polyline {
    quadtree: SimdQuadTree<u32>,
    aabb: AABB,
    vertices: Vec<Point<Real>>,
    indices: Vec<Point2<u32>>,
}

impl Polyline {
    /// Creates a new polyline from a vertex buffer and an index buffer.
    pub fn new(vertices: Vec<Point<Real>>, indices: Option<Vec<Point2<u32>>>) -> Self {
        let indices = indices.unwrap_or_else(|| {
            (0..vertices.len() as u32 - 2)
                .map(|i| Point2::new(i, i + 1))
                .collect()
        });
        let aabb = AABB::from_points(&vertices);
        let data = indices.iter().enumerate().map(|(i, idx)| {
            let aabb =
                Segment::new(vertices[idx[0] as usize], vertices[idx[1] as usize]).local_aabb();
            (i as u32, aabb)
        });

        let mut quadtree = SimdQuadTree::new();
        // NOTE: we apply no dilation factor because we won't
        // update this tree dynamically.
        quadtree.clear_and_rebuild(data, 0.0);

        Self {
            quadtree,
            aabb,
            vertices,
            indices,
        }
    }

    /// Compute the axis-aligned bounding box of this polyline.
    pub fn aabb(&self, pos: &Isometry<Real>) -> AABB {
        self.aabb.transform_by(pos)
    }

    /// Gets the local axis-aligned bounding box of this polyline.
    pub fn local_aabb(&self) -> &AABB {
        &self.aabb
    }

    pub(crate) fn quadtree(&self) -> &SimdQuadTree<u32> {
        &self.quadtree
    }

    /// The number of segments forming this polyline.
    pub fn num_segments(&self) -> usize {
        self.indices.len()
    }

    /// An iterator through all the segments of this mesh.
    pub fn segments(&self) -> impl Iterator<Item = Segment> + '_ {
        self.indices.iter().map(move |ids| {
            Segment::new(self.vertices[ids.x as usize], self.vertices[ids.y as usize])
        })
    }

    /// Get the `i`-th segment of this mesh.
    pub fn segment(&self, i: u32) -> Segment {
        let idx = self.indices[i as usize];
        Segment::new(self.vertices[idx.x as usize], self.vertices[idx.y as usize])
    }

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
    pub fn indices(&self) -> &[Point2<u32>] {
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
}

impl SimdCompositeShape for Polyline {
    fn nparts(&self) -> usize {
        self.num_segments()
    }

    fn map_part_at(&self, i: u32, f: &mut dyn FnMut(Option<&Isometry<Real>>, &dyn Shape)) {
        let tri = self.segment(i);
        f(None, &tri)
    }

    fn quadtree(&self) -> &SimdQuadTree<u32> {
        &self.quadtree
    }
}

impl TypedSimdCompositeShape for Polyline {
    type PartShape = Segment;

    #[inline(always)]
    fn map_typed_part_at(
        &self,
        i: u32,
        mut f: impl FnMut(Option<&Isometry<Real>>, &Self::PartShape),
    ) {
        let seg = self.segment(i);
        f(None, &seg)
    }
}
