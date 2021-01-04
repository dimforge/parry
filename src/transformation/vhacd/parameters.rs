#[derive(Copy, Clone, PartialEq, Eq)]
pub enum FillMode {
    SurfaceOnly,
    FloodFill,
    // RaycastFill
}

pub struct Parameters {
    pub concavity: Real,
    pub alpha: Real,
    pub beta: Real,
    pub min_volume_per_ch: Real,
    pub resolution: u32,
    pub max_num_vertices_per_ch: u32,
    pub plane_downsampling: u32,
    pub convex_hull_downsampling: u32,
    pub pca: u32,
    // 0: voxel-based (recommended), 1: tetrahedron-based.
    pub mode: u32,
    pub fill_mode: FillMode,
    pub convex_hull_approximation: bool,
    pub max_convex_hulls: u32,
    // This will project the output convex hull vertices onto the original source mesh to increase the floating point accuracy of the results.
    pub project_hull_vertices: bool,
}

impl Default for Parameters {
    fn default() -> Self {
        Self {
            resolution: 100000,
            concavity: 0.001,
            plane_downsampling: 4,
            convex_hull_downsampling: 4,
            alpha: 0.05,
            beta: 0.05,
            pca: 0,
            mode: 0,
            max_num_vertices_per_ch: 64,
            min_volume_per_ch: 0.0001,
            convex_hull_approximation: true,
            max_convex_hulls: 1024,
            project_hull_vertices: true,
            fill_mode: FillMode::FloodFill,
        }
    }
}
