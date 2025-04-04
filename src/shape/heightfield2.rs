use core::ops::Range;
#[cfg(not(feature = "std"))]
use na::ComplexField;
#[cfg(feature = "alloc")]
use na::DVector;

use na::Point2;

use crate::bounding_volume::Aabb;
use crate::math::{Real, Vector};

use crate::shape::Segment;

/// Indicates if a cell of a heightfield is removed or not. Set this to `false` for
/// a removed cell.
pub type HeightFieldCellStatus = bool;

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
// TODO: Archive isn’t implemented for VecStorage yet.
// #[cfg_attr(
//     feature = "rkyv",
//     derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize),
//     archive(check_bytes)
// )]
#[derive(Debug, Clone)]
#[repr(C)]
/// A 2D heightfield with a generic storage buffer for its heights.
pub struct HeightField {
    heights: DVector<Real>,
    status: DVector<HeightFieldCellStatus>,

    scale: Vector<Real>,
    aabb: Aabb,
}

#[cfg(feature = "alloc")]
impl HeightField {
    /// Creates a new 2D heightfield with the given heights and scale factor.
    pub fn new(heights: DVector<Real>, scale: Vector<Real>) -> Self {
        assert!(
            heights.len() > 1,
            "A heightfield heights must have at least 2 elements."
        );

        let max = heights.max();
        let min = heights.min();
        let hscale = scale * 0.5;
        let aabb = Aabb::new(
            Point2::new(-hscale.x, min * scale.y),
            Point2::new(hscale.x, max * scale.y),
        );
        let num_segments = heights.len() - 1;

        HeightField {
            heights,
            status: DVector::repeat(num_segments, true),
            scale,
            aabb,
        }
    }
}

impl HeightField {
    /// The number of cells of this heightfield.
    pub fn num_cells(&self) -> usize {
        self.heights.len() - 1
    }

    /// The height at each cell endpoint.
    pub fn heights(&self) -> &DVector<Real> {
        &self.heights
    }

    /// The scale factor applied to this heightfield.
    pub fn scale(&self) -> &Vector<Real> {
        &self.scale
    }

    /// Sets the scale factor applied to this heightfield.
    pub fn set_scale(&mut self, new_scale: Vector<Real>) {
        let ratio = new_scale.component_div(&self.scale);
        self.aabb.mins.coords.component_mul_assign(&ratio);
        self.aabb.maxs.coords.component_mul_assign(&ratio);
        self.scale = new_scale;
    }

    /// Returns a scaled version of this heightfield.
    pub fn scaled(mut self, scale: &Vector<Real>) -> Self {
        self.set_scale(self.scale.component_mul(scale));
        self
    }

    /// The [`Aabb`] of this heightfield.
    pub fn root_aabb(&self) -> &Aabb {
        &self.aabb
    }

    /// The width of a single cell of this heightfield.
    pub fn cell_width(&self) -> Real {
        self.unit_cell_width() * self.scale.x
    }

    /// The width of a single cell of this heightfield, without taking the scale factor into account.
    pub fn unit_cell_width(&self) -> Real {
        1.0 / (self.heights.len() as Real - 1.0)
    }

    /// The left-most x-coordinate of this heightfield.
    pub fn start_x(&self) -> Real {
        self.scale.x * -0.5
    }

    fn quantize_floor_unclamped(&self, val: Real, seg_length: Real) -> isize {
        ((val + 0.5) / seg_length).floor() as isize
    }

    fn quantize_ceil_unclamped(&self, val: Real, seg_length: Real) -> isize {
        ((val + 0.5) / seg_length).ceil() as isize
    }

    fn quantize_floor(&self, val: Real, seg_length: Real) -> usize {
        na::clamp(
            ((val + 0.5) / seg_length).floor(),
            0.0,
            (self.num_cells() - 1) as Real,
        ) as usize
    }

    fn quantize_ceil(&self, val: Real, seg_length: Real) -> usize {
        na::clamp(
            ((val + 0.5) / seg_length).ceil(),
            0.0,
            self.num_cells() as Real,
        ) as usize
    }

    /// Index of the cell a point is on after vertical projection.
    pub fn cell_at_point(&self, pt: &Point2<Real>) -> Option<usize> {
        let scaled_pt = pt.coords.component_div(&self.scale);
        let seg_length = self.unit_cell_width();

        if scaled_pt.x < -0.5 || scaled_pt.x > 0.5 {
            // Outside of the heightfield bounds.
            None
        } else {
            Some(self.quantize_floor(scaled_pt.x, seg_length))
        }
    }

    /// Height of the heightfield at the given point after vertical projection on the heightfield surface.
    pub fn height_at_point(&self, pt: &Point2<Real>) -> Option<Real> {
        let cell = self.cell_at_point(pt)?;
        let seg = self.segment_at(cell)?;
        let inter = crate::query::details::closest_points_line_line_parameters(
            &seg.a,
            &seg.scaled_direction(),
            pt,
            &Vector::y(),
        );
        Some(seg.a.y + inter.1)
    }

    /// Iterator through all the segments of this heightfield.
    pub fn segments(&self) -> impl Iterator<Item = Segment> + '_ {
        // TODO: this is not very efficient since this will
        // recompute shared points twice.
        (0..self.num_cells()).filter_map(move |i| self.segment_at(i))
    }

    /// The i-th segment of the heightfield if it has not been removed.
    pub fn segment_at(&self, i: usize) -> Option<Segment> {
        if i >= self.num_cells() || self.is_segment_removed(i) {
            return None;
        }

        let seg_length = 1.0 / (self.heights.len() as Real - 1.0);

        let x0 = -0.5 + seg_length * (i as Real);
        let x1 = x0 + seg_length;

        let y0 = self.heights[i];
        let y1 = self.heights[i + 1];

        let mut p0 = Point2::new(x0, y0);
        let mut p1 = Point2::new(x1, y1);

        // Apply scales:
        p0.coords.component_mul_assign(&self.scale);
        p1.coords.component_mul_assign(&self.scale);

        Some(Segment::new(p0, p1))
    }

    /// Mark the i-th segment of this heightfield as removed or not.
    pub fn set_segment_removed(&mut self, i: usize, removed: bool) {
        self.status[i] = !removed
    }

    /// Checks if the i-th segment has been removed.
    pub fn is_segment_removed(&self, i: usize) -> bool {
        !self.status[i]
    }

    /// The range of segment ids that may intersect the given local Aabb.
    pub fn unclamped_elements_range_in_local_aabb(&self, aabb: &Aabb) -> Range<isize> {
        let ref_mins = aabb.mins.coords.component_div(&self.scale);
        let ref_maxs = aabb.maxs.coords.component_div(&self.scale);
        let seg_length = 1.0 / (self.heights.len() as Real - 1.0);

        let min_x = self.quantize_floor_unclamped(ref_mins.x, seg_length);
        let max_x = self.quantize_ceil_unclamped(ref_maxs.x, seg_length);
        min_x..max_x
    }

    /// Applies `f` to each segment of this heightfield that intersects the given `aabb`.
    pub fn map_elements_in_local_aabb(&self, aabb: &Aabb, f: &mut impl FnMut(u32, &Segment)) {
        let ref_mins = aabb.mins.coords.component_div(&self.scale);
        let ref_maxs = aabb.maxs.coords.component_div(&self.scale);
        let seg_length = 1.0 / (self.heights.len() as Real - 1.0);

        if ref_maxs.x < -0.5 || ref_mins.x > 0.5 {
            // Outside of the heightfield bounds.
            return;
        }

        let min_x = self.quantize_floor(ref_mins.x, seg_length);
        let max_x = self.quantize_ceil(ref_maxs.x, seg_length);

        // TODO: find a way to avoid recomputing the same vertices
        // multiple times.
        for i in min_x..max_x {
            if self.is_segment_removed(i) {
                continue;
            }

            let x0 = -0.5 + seg_length * (i as Real);
            let x1 = x0 + seg_length;

            let y0 = self.heights[i];
            let y1 = self.heights[i + 1];

            if (y0 > ref_maxs.y && y1 > ref_maxs.y) || (y0 < ref_mins.y && y1 < ref_mins.y) {
                continue;
            }

            let mut p0 = Point2::new(x0, y0);
            let mut p1 = Point2::new(x1, y1);

            // Apply scales:
            p0.coords.component_mul_assign(&self.scale);
            p1.coords.component_mul_assign(&self.scale);

            // Build the segment.
            let seg = Segment::new(p0, p1);

            // Call the callback.
            f(i as u32, &seg);
        }
    }
}
