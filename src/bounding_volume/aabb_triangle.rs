use crate::math::*;
use crate::{bounding_volume::Aabb, shape::Triangle};

impl Triangle {
    /// Computes the world-space [`Aabb`] of this triangle, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry) -> Aabb {
        self.transformed(pos).local_aabb()
    }

    /// Computes the local-space [`Aabb`] of this triangle.
    #[inline]
    pub fn local_aabb(&self) -> Aabb {
        let a = self.a.as_vector();
        let b = self.b.as_vector();
        let c = self.c.as_vector();

        let mut min = Point::origin();
        let mut max = Point::origin();

        for d in 0..DIM {
            min.as_vector_mut()[d] = a[d].min(b[d]).min(c[d]);
            max.as_vector_mut()[d] = a[d].max(b[d]).max(c[d]);
        }

        Aabb::new(min, max)
    }
}

#[cfg(test)]
#[cfg(feature = "dim3")]
mod test {
    use crate::{
        bounding_volume::details::support_map_aabb,
        math::{Isometry, Point, Real, Vector},
        shape::Triangle,
    };
    use na::RealField;

    #[test]
    fn triangle_aabb_matches_support_map_aabb() {
        let t = Triangle::new(
            Point::new(0.3, -0.1, 0.2),
            Point::new(-0.7, 1.0, 0.0),
            Point::new(-0.7, 1.5, 0.0),
        );

        let m = Isometry::new(
            Vector::new(-0.2, 5.0, 0.2),
            Vector::new(0.0, Real::frac_pi_2(), 0.0),
        );

        assert_eq!(t.aabb(&m), support_map_aabb(&m, &t));

        // TODO: also test local Aabb once support maps have a local Aabb
        // function too
    }
}
