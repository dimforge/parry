use crate::math::Real;
use crate::transformation::voxelization::FillMode;

pub struct VHACDParameters {
    pub concavity: Real,
    pub alpha: Real,
    pub beta: Real,
    pub resolution: u32,
    pub plane_downsampling: u32,
    pub convex_hull_downsampling: u32,
    pub fill_mode: FillMode,
    pub convex_hull_approximation: bool,
    pub max_convex_hulls: u32,
}

impl Default for VHACDParameters {
    fn default() -> Self {
        Self {
            #[cfg(feature = "dim3")]
            resolution: 64,
            #[cfg(feature = "dim3")]
            concavity: 0.001,
            #[cfg(feature = "dim2")]
            resolution: 512,
            #[cfg(feature = "dim2")]
            concavity: 0.01,
            plane_downsampling: 4,
            convex_hull_downsampling: 4,
            alpha: 0.05,
            beta: 0.05,
            convex_hull_approximation: true,
            max_convex_hulls: 1024,
            fill_mode: FillMode::FloodFill,
        }
    }
}
