use crate::{
    bounding_volume::Aabb,
    math::{Isometry, Point, Real, DIM},
    shape::Triangle,
};

impl Triangle {
    /// Computes the world-space Aabb of this triangle, transformed by `pos`.
    #[inline]
    pub fn aabb(&self, pos: &Isometry<Real>) -> Aabb {
        self.transformed(pos).local_aabb()
    }

    /// Computes the local-space Aabb of this triangle.
    #[inline]
    pub fn local_aabb(&self) -> Aabb {
        let a = self.a.coords;
        let b = self.b.coords;
        let c = self.c.coords;

        let mut min = Point::origin();
        let mut max = Point::origin();

        for d in 0..DIM {
            min.coords[d] = a[d].min(b[d]).min(c[d]);
            max.coords[d] = a[d].max(b[d]).max(c[d]);
        }

        Aabb::new(min, max)
    }
}

#[cfg(test)]
#[cfg(feature = "dim3")]
mod test {
    use crate::{
        bounding_volume::details::support_map_aabb,
        math::{Isometry, Point, Real, Translation, real},
        shape::Triangle,
    };
    use na::{RealField, UnitQuaternion};

    #[test]
    fn triangle_aabb_matches_support_map_aabb() {
        let t = Triangle::new(
            Point::new(real!(0.3), real!(-0.1), real!(0.2)),
            Point::new(real!(-0.7), real!(1.0), real!(0.0)),
            Point::new(real!(-0.7), real!(1.5), real!(0.0)),
        );

        let m = Isometry::from_parts(
            Translation::new(real!(-0.2), real!(5.0), real!(0.2)),
            UnitQuaternion::from_euler_angles(real!(0.0), Real::frac_pi_2(), real!(0.0)),
        );

        assert_eq!(t.aabb(&m), support_map_aabb(&m, &t));

        // TODO: also test local Aabb once support maps have a local Aabb
        // function too
    }
}
