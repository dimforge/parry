#[cfg(not(feature = "std"))]
use na::ComplexField;
#[cfg(feature = "std")]
use na::DVector;
#[cfg(all(feature = "std", feature = "cuda"))]
use {crate::utils::CudaArray1, cust::error::CudaResult};

use na::{Point2, Scalar};

use crate::bounding_volume::AABB;
use crate::math::{Real, Vector};

use crate::shape::Segment;

/// Indicates if a cell of an heightfield is removed or not. Set this to `false` for
/// a removed cell.
pub type HeightFieldCellStatus = bool;

/// Abstraction over the storage type of an heightfield’s heights.
pub trait HeightFieldStorage {
    /// The type of heights.
    type Item;
    /// The number of heights on this storage.
    fn len(&self) -> usize;
    /// Gets the `i`-th height of the heightfield.
    fn get(&self, i: usize) -> Self::Item;
    /// Sets the `i`-th height of the heightfield.
    fn set(&mut self, i: usize, val: Self::Item);
}

#[cfg(feature = "std")]
impl<T: Scalar> HeightFieldStorage for DVector<T> {
    type Item = T;

    #[inline]
    fn len(&self) -> usize {
        self.len()
    }

    #[inline]
    fn get(&self, i: usize) -> Self::Item {
        self[i].clone()
    }

    #[inline]
    fn set(&mut self, i: usize, val: Self::Item) {
        self[i] = val
    }
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "cuda", derive(cust_core::DeviceCopy))]
#[derive(Copy, Clone, Debug)]
#[repr(C)] // Needed for Cuda.
/// A 2D heightfield with a generic storage buffer for its heights.
pub struct GenericHeightField<Heights, Status> {
    heights: Heights,
    status: Status,

    scale: Vector<Real>,
    aabb: AABB,
}

/// A 2D heightfield.
#[cfg(feature = "std")]
pub type HeightField = GenericHeightField<DVector<Real>, DVector<HeightFieldCellStatus>>;

/// A 2D heightfield stored in the CUDA memory, initializable from the host.
#[cfg(all(feature = "std", feature = "cuda"))]
pub type CudaHeightField = GenericHeightField<CudaArray1<Real>, CudaArray1<HeightFieldCellStatus>>;

/// A 3D heightfield stored in the CUDA memory, accessible from within a Cuda kernel.
#[cfg(feature = "cuda")]
pub type CudaHeightFieldPointer = GenericHeightField<
    crate::utils::CudaArrayPointer1<Real>,
    crate::utils::CudaArrayPointer1<HeightFieldCellStatus>,
>;

#[cfg(feature = "std")]
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
        let aabb = AABB::new(
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

    /// Converts this RAM-based heightfield to an heightfield based on CUDA memory.
    #[cfg(all(feature = "cuda"))]
    pub fn to_cuda(&self) -> CudaResult<CudaHeightField> {
        Ok(CudaHeightField {
            heights: CudaArray1::from_vector(&self.heights)?,
            status: CudaArray1::from_vector(&self.status)?,
            scale: self.scale,
            aabb: self.aabb,
        })
    }
}

#[cfg(all(feature = "std", feature = "cuda"))]
impl CudaHeightField {
    /// Returns the heightfield usable from within a CUDA kernel.
    pub fn as_device_ptr(&self) -> CudaHeightFieldPointer {
        CudaHeightFieldPointer {
            heights: self.heights.as_device_ptr(),
            status: self.status.as_device_ptr(),
            aabb: self.aabb,
            scale: self.scale,
        }
    }
}

impl<Heights, Status> GenericHeightField<Heights, Status>
where
    Heights: HeightFieldStorage<Item = Real>,
    Status: HeightFieldStorage<Item = HeightFieldCellStatus>,
{
    /// The number of cells of this heightfield.
    pub fn num_cells(&self) -> usize {
        self.heights.len() - 1
    }

    /// The height at each cell endpoint.
    pub fn heights(&self) -> &Heights {
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
        self.set_scale(self.scale.component_mul(&scale));
        self
    }

    /// The AABB of this heightfield.
    pub fn root_aabb(&self) -> &AABB {
        &self.aabb
    }

    /// The width of a single cell of this heightfield.
    pub fn cell_width(&self) -> Real {
        self.unit_cell_width() * self.scale.x
    }

    /// The width of a single cell of this heightfield, without taking the scale factor into account.
    pub fn unit_cell_width(&self) -> Real {
        1.0 / na::convert::<f64, Real>(self.heights.len() as f64 - 1.0)
    }

    /// The left-most x-coordinate of this heightfield.
    pub fn start_x(&self) -> Real {
        self.scale.x * na::convert::<f64, Real>(-0.5)
    }

    fn quantize_floor(&self, val: Real, seg_length: Real) -> usize {
        let _0_5: Real = na::convert::<f64, Real>(0.5);
        let i = na::clamp(
            ((val + _0_5) / seg_length).floor(),
            0.0,
            na::convert::<f64, Real>((self.num_cells() - 1) as f64),
        );
        na::convert_unchecked::<Real, f64>(i) as usize
    }

    fn quantize_ceil(&self, val: Real, seg_length: Real) -> usize {
        let _0_5: Real = na::convert::<f64, Real>(0.5);
        let i = na::clamp(
            ((val + _0_5) / seg_length).ceil(),
            0.0,
            na::convert::<f64, Real>(self.num_cells() as f64),
        );
        na::convert_unchecked::<Real, f64>(i) as usize
    }

    /// Index of the cell a point is on after vertical projection.
    pub fn cell_at_point(&self, pt: &Point2<Real>) -> Option<usize> {
        let _0_5: Real = na::convert::<f64, Real>(0.5);
        let scaled_pt = pt.coords.component_div(&self.scale);
        let seg_length = self.unit_cell_width();

        if scaled_pt.x < -_0_5 || scaled_pt.x > _0_5 {
            // Outside of the heightfield bounds.
            None
        } else {
            Some(self.quantize_floor(scaled_pt.x, seg_length))
        }
    }

    /// Height of the heightfield a the given point after vertical projection on the heightfield surface.
    pub fn height_at_point(&self, pt: &Point2<Real>) -> Option<Real> {
        let cell = self.cell_at_point(pt)?;
        let seg = self.segment_at(cell)?;
        let inter = crate::query::details::closest_points_line_line_parameters(
            &seg.a,
            &seg.scaled_direction(),
            &pt,
            &Vector::y(),
        );
        Some(seg.a.y + inter.1)
    }

    /// Iterator through all the segments of this heightfield.
    pub fn segments<'a>(&'a self) -> impl Iterator<Item = Segment> + 'a {
        // FIXME: this is not very efficient since this wil
        // recompute shared points twice.
        (0..self.num_cells()).filter_map(move |i| self.segment_at(i))
    }

    /// The i-th segment of the heightfield if it has not been removed.
    pub fn segment_at(&self, i: usize) -> Option<Segment> {
        if i >= self.num_cells() || self.is_segment_removed(i) {
            return None;
        }

        let _0_5: Real = na::convert::<f64, Real>(0.5);
        let seg_length = 1.0 / na::convert::<f64, Real>(self.heights.len() as f64 - 1.0);

        let x0 = -_0_5 + seg_length * na::convert::<f64, Real>(i as f64);
        let x1 = x0 + seg_length;

        let y0 = self.heights.get(i + 0);
        let y1 = self.heights.get(i + 1);

        let mut p0 = Point2::new(x0, y0);
        let mut p1 = Point2::new(x1, y1);

        // Apply scales:
        p0.coords.component_mul_assign(&self.scale);
        p1.coords.component_mul_assign(&self.scale);

        Some(Segment::new(p0, p1))
    }

    /// Mark the i-th segment of this heightfield as removed or not.
    pub fn set_segment_removed(&mut self, i: usize, removed: bool) {
        self.status.set(i, !removed)
    }

    /// Checks if the i-th segment has been removed.
    pub fn is_segment_removed(&self, i: usize) -> bool {
        !self.status.get(i)
    }

    /// Applies `f` to each segment of this heightfield that intersects the given `aabb`.
    pub fn map_elements_in_local_aabb(&self, aabb: &AABB, f: &mut impl FnMut(u32, &Segment)) {
        let _0_5: Real = na::convert::<f64, Real>(0.5);
        let ref_mins = aabb.mins.coords.component_div(&self.scale);
        let ref_maxs = aabb.maxs.coords.component_div(&self.scale);
        let seg_length = 1.0 / na::convert::<f64, Real>(self.heights.len() as f64 - 1.0);

        if ref_maxs.x < -_0_5 || ref_mins.x > _0_5 {
            // Outside of the heightfield bounds.
            return;
        }

        let min_x = self.quantize_floor(ref_mins.x, seg_length);
        let max_x = self.quantize_ceil(ref_maxs.x, seg_length);

        // FIXME: find a way to avoid recomputing the same vertices
        // multiple times.
        for i in min_x..max_x {
            if self.is_segment_removed(i) {
                continue;
            }

            let x0 = -_0_5 + seg_length * na::convert::<f64, Real>(i as f64);
            let x1 = x0 + seg_length;

            let y0 = self.heights.get(i + 0);
            let y1 = self.heights.get(i + 1);

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
