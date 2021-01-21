use na::{DMatrix, Point3};

use crate::bounding_volume::AABB;
use crate::math::{Real, Vector};
use crate::shape::{FeatureId, Triangle};

bitflags! {
    #[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
    #[derive(Default)]
    /// The status of the cell of an heightfield.
    pub struct HeightFieldCellStatus: u8 {
        /// If this bit is set, the concerned heightfield cell is subdivided using a Z pattern.
        const ZIGZAG_SUBDIVISION = 0b00000001;
        /// If this bit is set, the leftmost triangle of the concerned heightfield cell is removed.
        const LEFT_TRIANGLE_REMOVED = 0b00000010;
        /// If this bit is set, the rightmost triangle of the concerned heightfield cell is removed.
        const RIGHT_TRIANGLE_REMOVED = 0b00000100;
        /// If this bit is set, both triangles of the concerned heightfield cell are removed.
        const CELL_REMOVED = Self::LEFT_TRIANGLE_REMOVED.bits | Self::RIGHT_TRIANGLE_REMOVED.bits;
    }
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
/// An heightfield implicitly discretized with triangles.
pub struct HeightField {
    heights: DMatrix<Real>,
    scale: Vector<Real>,
    aabb: AABB,
    num_triangles: usize,
    status: DMatrix<HeightFieldCellStatus>,
}

impl HeightField {
    /// Initializes a new heightfield with the given heights and a scaling factor.
    pub fn new(heights: DMatrix<Real>, scale: Vector<Real>) -> Self {
        assert!(
            heights.nrows() > 1 && heights.ncols() > 1,
            "A heightfield heights must have at least 2 rows and columns."
        );
        let max = heights.max();
        let min = heights.min();
        let hscale = scale * na::convert::<_, Real>(0.5);
        let aabb = AABB::new(
            Point3::new(-hscale.x, min * scale.y, -hscale.z),
            Point3::new(hscale.x, max * scale.y, hscale.z),
        );
        let num_triangles = (heights.nrows() - 1) * (heights.ncols() - 1) * 2;
        let status = DMatrix::repeat(
            heights.nrows() - 1,
            heights.ncols() - 1,
            HeightFieldCellStatus::default(),
        );

        HeightField {
            heights,
            scale,
            aabb,
            num_triangles,
            status,
        }
    }

    /// The number of rows of this heightfield.
    pub fn nrows(&self) -> usize {
        self.heights.nrows() - 1
    }

    /// The number of columns of this heightfield.
    pub fn ncols(&self) -> usize {
        self.heights.ncols() - 1
    }

    fn triangle_id(&self, i: usize, j: usize, left: bool) -> u32 {
        let tid = j * (self.heights.nrows() - 1) + i;
        if left {
            tid as u32
        } else {
            (tid + self.num_triangles / 2) as u32
        }
    }

    fn face_id(&self, i: usize, j: usize, left: bool, front: bool) -> u32 {
        let tid = self.triangle_id(i, j, left);
        if front {
            tid
        } else {
            tid + self.num_triangles as u32
        }
    }

    fn quantize_floor(&self, val: Real, cell_size: Real, num_cells: usize) -> usize {
        let _0_5: Real = na::convert::<f64, Real>(0.5);
        let i = na::clamp(
            ((val + _0_5) / cell_size).floor(),
            0.0,
            na::convert::<f64, Real>((num_cells - 1) as f64),
        );
        na::convert_unchecked::<Real, f64>(i) as usize
    }

    fn quantize_ceil(&self, val: Real, cell_size: Real, num_cells: usize) -> usize {
        let _0_5: Real = na::convert::<f64, Real>(0.5);
        let i = na::clamp(
            ((val + _0_5) / cell_size).ceil(),
            0.0,
            na::convert::<f64, Real>(num_cells as f64),
        );
        na::convert_unchecked::<Real, f64>(i) as usize
    }

    /// The pair of index of the cell containing the vertical projection of the given point.
    pub fn cell_at_point(&self, pt: &Point3<Real>) -> Option<(usize, usize)> {
        let _0_5: Real = na::convert::<f64, Real>(0.5);
        let scaled_pt = pt.coords.component_div(&self.scale);
        let cell_width = self.unit_cell_width();
        let cell_height = self.unit_cell_height();
        let ncells_x = self.ncols();
        let ncells_z = self.nrows();

        if scaled_pt.x < -_0_5 || scaled_pt.x > _0_5 || scaled_pt.z < -_0_5 || scaled_pt.z > _0_5 {
            // Outside of the heightfield bounds.
            None
        } else {
            let j = self.quantize_floor(scaled_pt.x, cell_width, ncells_x);
            let i = self.quantize_floor(scaled_pt.z, cell_height, ncells_z);
            Some((i, j))
        }
    }

    /// The smallest x coordinate of the `j`-th column of this heightfield.
    pub fn x_at(&self, j: usize) -> Real {
        let _0_5: Real = na::convert::<f64, Real>(0.5);
        (-_0_5 + self.unit_cell_width() * na::convert::<f64, Real>(j as f64)) * self.scale.x
    }

    /// The smallest z coordinate of the start of the `i`-th row of this heightfield.
    pub fn z_at(&self, i: usize) -> Real {
        let _0_5: Real = na::convert::<f64, Real>(0.5);
        (-_0_5 + self.unit_cell_height() * na::convert::<f64, Real>(i as f64)) * self.scale.z
    }

    /// An iterator through all the triangles of this heightfield.
    pub fn triangles<'a>(&'a self) -> impl Iterator<Item = Triangle> + 'a {
        HeightfieldTriangles {
            heightfield: self,
            curr: (0, 0),
            tris: self.triangles_at(0, 0),
        }
    }

    /// The two triangles at the cell (i, j) of this heightfield.
    ///
    /// Returns `None` fore triangles that have been removed because of their user-defined status
    /// flags (described by the `HeightFieldCellStatus` bitfield).
    pub fn triangles_at(&self, i: usize, j: usize) -> (Option<Triangle>, Option<Triangle>) {
        let status = self.status[(i, j)];

        if status.contains(
            HeightFieldCellStatus::LEFT_TRIANGLE_REMOVED
                | HeightFieldCellStatus::RIGHT_TRIANGLE_REMOVED,
        ) {
            return (None, None);
        }

        let cell_width = self.unit_cell_width();
        let cell_height = self.unit_cell_height();

        let _0_5: Real = na::convert::<f64, Real>(0.5);
        let z0 = -_0_5 + cell_height * na::convert::<f64, Real>(i as f64);
        let z1 = z0 + cell_height;

        let x0 = -_0_5 + cell_width * na::convert::<f64, Real>(j as f64);
        let x1 = x0 + cell_width;

        let y00 = self.heights[(i + 0, j + 0)];
        let y10 = self.heights[(i + 1, j + 0)];
        let y01 = self.heights[(i + 0, j + 1)];
        let y11 = self.heights[(i + 1, j + 1)];

        let mut p00 = Point3::new(x0, y00, z0);
        let mut p10 = Point3::new(x0, y10, z1);
        let mut p01 = Point3::new(x1, y01, z0);
        let mut p11 = Point3::new(x1, y11, z1);

        // Apply scales:
        p00.coords.component_mul_assign(&self.scale);
        p10.coords.component_mul_assign(&self.scale);
        p01.coords.component_mul_assign(&self.scale);
        p11.coords.component_mul_assign(&self.scale);

        if status.contains(HeightFieldCellStatus::ZIGZAG_SUBDIVISION) {
            let tri1 = if status.contains(HeightFieldCellStatus::LEFT_TRIANGLE_REMOVED) {
                None
            } else {
                Some(Triangle::new(p00, p10, p11))
            };

            let tri2 = if status.contains(HeightFieldCellStatus::RIGHT_TRIANGLE_REMOVED) {
                None
            } else {
                Some(Triangle::new(p00, p11, p01))
            };

            (tri1, tri2)
        } else {
            let tri1 = if status.contains(HeightFieldCellStatus::LEFT_TRIANGLE_REMOVED) {
                None
            } else {
                Some(Triangle::new(p00, p10, p01))
            };

            let tri2 = if status.contains(HeightFieldCellStatus::RIGHT_TRIANGLE_REMOVED) {
                None
            } else {
                Some(Triangle::new(p10, p11, p01))
            };

            (tri1, tri2)
        }
    }

    /// The status of the `(i, j)`-th cell.
    pub fn cell_status(&self, i: usize, j: usize) -> HeightFieldCellStatus {
        self.status[(i, j)]
    }

    /// Set the status of the `(i, j)`-th cell.
    pub fn set_cell_status(&mut self, i: usize, j: usize, status: HeightFieldCellStatus) {
        self.status[(i, j)] = status
    }

    /// The statuses of all the cells of this heightfield.
    pub fn cells_statuses(&self) -> &DMatrix<HeightFieldCellStatus> {
        &self.status
    }

    /// The mutable statuses of all the cells of this heightfield.
    pub fn cells_statuses_mut(&mut self) -> &mut DMatrix<HeightFieldCellStatus> {
        &mut self.status
    }

    /// The heights of this heightfield.
    pub fn heights(&self) -> &DMatrix<Real> {
        &self.heights
    }

    /// The scale factor applied to this heightfield.
    pub fn scale(&self) -> &Vector<Real> {
        &self.scale
    }

    /// The width (extent along its local `x` axis) of each cell of this heightmap, including the scale factor.
    pub fn cell_width(&self) -> Real {
        self.unit_cell_width() * self.scale.x
    }

    /// The height (extent along its local `z` axis) of each cell of this heightmap, including the scale factor.
    pub fn cell_height(&self) -> Real {
        self.unit_cell_height() * self.scale.z
    }

    /// The width (extent along its local `x` axis) of each cell of this heightmap, excluding the scale factor.
    pub fn unit_cell_width(&self) -> Real {
        1.0 / na::convert::<f64, Real>(self.heights.ncols() as f64 - 1.0)
    }

    /// The height (extent along its local `z` axis) of each cell of this heightmap, excluding the scale factor.
    pub fn unit_cell_height(&self) -> Real {
        1.0 / na::convert::<f64, Real>(self.heights.nrows() as f64 - 1.0)
    }

    /// The AABB of this heightmap.
    pub fn root_aabb(&self) -> &AABB {
        &self.aabb
    }

    /// Converts the FeatureID of the left or right triangle at the cell `(i, j)` into a FeatureId
    /// of the whole heightfield.
    pub fn convert_triangle_feature_id(
        &self,
        i: usize,
        j: usize,
        left: bool,
        fid: FeatureId,
    ) -> FeatureId {
        match fid {
            FeatureId::Vertex(ivertex) => {
                let nrows = self.heights.nrows();
                let ij = i + j * nrows;

                if self.status[(i, j)].contains(HeightFieldCellStatus::ZIGZAG_SUBDIVISION) {
                    if left {
                        FeatureId::Vertex([ij, ij + 1, ij + 1 + nrows][ivertex as usize] as u32)
                    } else {
                        FeatureId::Vertex([ij, ij + 1 + nrows, ij + nrows][ivertex as usize] as u32)
                    }
                } else {
                    if left {
                        FeatureId::Vertex([ij, ij + 1, ij + nrows][ivertex as usize] as u32)
                    } else {
                        FeatureId::Vertex(
                            [ij + 1, ij + 1 + nrows, ij + nrows][ivertex as usize] as u32,
                        )
                    }
                }
            }
            FeatureId::Edge(iedge) => {
                let (nrows, ncols) = self.heights.shape();
                let vshift = 0; // First vertical line index.
                let hshift = (nrows - 1) * ncols; // First horizontal line index.
                let dshift = hshift + nrows * (ncols - 1); // First diagonal line index.
                let idiag = dshift + i + j * (nrows - 1);
                let itop = hshift + i + j * nrows;
                let ibottom = itop + 1;
                let ileft = vshift + i + j * (nrows - 1);
                let iright = ileft + nrows - 1;

                if self.status[(i, j)].contains(HeightFieldCellStatus::ZIGZAG_SUBDIVISION) {
                    if left {
                        // Triangle:
                        //
                        // |\
                        // |_\
                        //
                        FeatureId::Edge([ileft, ibottom, idiag][iedge as usize] as u32)
                    } else {
                        // Triangle:
                        // ___
                        // \ |
                        //  \|
                        //
                        FeatureId::Edge([idiag, iright, itop][iedge as usize] as u32)
                    }
                } else {
                    if left {
                        // Triangle:
                        // ___
                        // | /
                        // |/
                        //
                        FeatureId::Edge([ileft, idiag, itop][iedge as usize] as u32)
                    } else {
                        // Triangle:
                        //
                        //  /|
                        // /_|
                        //
                        FeatureId::Edge([ibottom, iright, idiag][iedge as usize] as u32)
                    }
                }
            }
            FeatureId::Face(iface) => {
                if iface == 0 {
                    FeatureId::Face(self.face_id(i, j, left, true))
                } else {
                    FeatureId::Face(self.face_id(i, j, left, false))
                }
            }
            FeatureId::Unknown => FeatureId::Unknown,
        }
    }

    /// Applies the function `f` to all the triangles of this heightfield intersecting the given AABB.
    pub fn map_elements_in_local_aabb(&self, aabb: &AABB, f: &mut impl FnMut(u32, &Triangle)) {
        let _0_5: Real = na::convert::<f64, Real>(0.5);
        let ncells_x = self.ncols();
        let ncells_z = self.nrows();

        let ref_mins = aabb.mins.coords.component_div(&self.scale);
        let ref_maxs = aabb.maxs.coords.component_div(&self.scale);
        let cell_width = self.unit_cell_width();
        let cell_height = self.unit_cell_height();

        if ref_maxs.x <= -_0_5 || ref_maxs.z <= -_0_5 || ref_mins.x >= _0_5 || ref_mins.z >= _0_5 {
            // Outside of the heightfield bounds.
            return;
        }

        let min_x = self.quantize_floor(ref_mins.x, cell_width, ncells_x);
        let min_z = self.quantize_floor(ref_mins.z, cell_height, ncells_z);

        let max_x = self.quantize_ceil(ref_maxs.x, cell_width, ncells_x);
        let max_z = self.quantize_ceil(ref_maxs.z, cell_height, ncells_z);

        // FIXME: find a way to avoid recomputing the same vertices
        // multiple times.
        for j in min_x..max_x {
            for i in min_z..max_z {
                let status = self.status[(i, j)];

                if status.contains(HeightFieldCellStatus::CELL_REMOVED) {
                    continue;
                }

                let z0 = -_0_5 + cell_height * na::convert::<f64, Real>(i as f64);
                let z1 = z0 + cell_height;

                let x0 = -_0_5 + cell_width * na::convert::<f64, Real>(j as f64);
                let x1 = x0 + cell_width;

                let y00 = self.heights[(i + 0, j + 0)];
                let y10 = self.heights[(i + 1, j + 0)];
                let y01 = self.heights[(i + 0, j + 1)];
                let y11 = self.heights[(i + 1, j + 1)];

                if (y00 > ref_maxs.y && y10 > ref_maxs.y && y01 > ref_maxs.y && y11 > ref_maxs.y)
                    || (y00 < ref_mins.y
                        && y10 < ref_mins.y
                        && y01 < ref_mins.y
                        && y11 < ref_mins.y)
                {
                    continue;
                }

                let mut p00 = Point3::new(x0, y00, z0);
                let mut p10 = Point3::new(x0, y10, z1);
                let mut p01 = Point3::new(x1, y01, z0);
                let mut p11 = Point3::new(x1, y11, z1);

                // Apply scales:
                p00.coords.component_mul_assign(&self.scale);
                p10.coords.component_mul_assign(&self.scale);
                p01.coords.component_mul_assign(&self.scale);
                p11.coords.component_mul_assign(&self.scale);

                // Build the two triangles, contact processors and call f.
                if !status.contains(HeightFieldCellStatus::LEFT_TRIANGLE_REMOVED) {
                    let tri1 = if status.contains(HeightFieldCellStatus::ZIGZAG_SUBDIVISION) {
                        Triangle::new(p00, p10, p11)
                    } else {
                        Triangle::new(p00, p10, p01)
                    };

                    let tid = self.triangle_id(i, j, true);
                    f(tid as u32, &tri1);
                }

                if !status.contains(HeightFieldCellStatus::RIGHT_TRIANGLE_REMOVED) {
                    let tri2 = if status.contains(HeightFieldCellStatus::ZIGZAG_SUBDIVISION) {
                        Triangle::new(p00, p11, p01)
                    } else {
                        Triangle::new(p10, p11, p01)
                    };
                    let tid = self.triangle_id(i, j, false);
                    f(tid as u32, &tri2);
                }
            }
        }
    }
}

struct HeightfieldTriangles<'a> {
    heightfield: &'a HeightField,
    curr: (usize, usize),
    tris: (Option<Triangle>, Option<Triangle>),
}

impl<'a> Iterator for HeightfieldTriangles<'a> {
    type Item = Triangle;

    fn next(&mut self) -> Option<Triangle> {
        loop {
            if let Some(tri1) = self.tris.0.take() {
                return Some(tri1);
            } else if let Some(tri2) = self.tris.1.take() {
                return Some(tri2);
            } else {
                self.curr.0 += 1;

                if self.curr.0 >= self.heightfield.nrows() {
                    if self.curr.1 >= self.heightfield.ncols() - 1 {
                        return None;
                    }

                    self.curr.0 = 0;
                    self.curr.1 += 1;
                }

                // tri1 and tri2 are None
                self.tris = self.heightfield.triangles_at(self.curr.0, self.curr.1);
            }
        }
    }
}
