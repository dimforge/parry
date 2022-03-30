use crate::{
    bounding_volume::AABB,
    math::{Isometry, Real, UnitVector},
    shape::Cuboid,
};

pub enum SplitResult<T> {
    Pair(T, T),
    Negative,
    Positive,
}

pub trait CanonicalSplit: Sized {
    /// Splits this shape along the given canonical axis.
    ///
    /// This will split the shape by a plane with a normal with it’s `axis`-th component set to 1.
    /// The splitting plane is shifted wrt. the origin by the `bias` (i.e. it passes through the point
    /// equal to `normal * bias`).
    ///
    /// # Result
    /// Returns the result of the split. The first sahpe returned is the piece lying on the negative
    /// half-space delimited by the splitting plane. The second sahpe returned is the piece lying on the
    /// positive half-space delimited by the splitting plane.
    fn canonical_split(&self, axis: usize, bias: Real, eps: Real) -> SplitResult<Self>;

    /// Compute the intersection volume between `self` and the given AABB.
    fn intersection_with_aabb(&self, _aabb: AABB, _eps: Real) -> Option<Self> {
        todo!()
    }
}

pub trait Split: Sized + CanonicalSplit {
    fn split(
        &self,
        position: &Isometry<Real>,
        axis: &UnitVector<Real>,
        bias: Real,
        epsilon: Real,
    ) -> SplitResult<Self> {
        let local_axis = position.inverse_transform_unit_vector(axis);
        let added_bias = -position.translation.vector.dot(&axis);
        self.local_split(&local_axis, bias + added_bias, epsilon)
    }

    fn local_split(
        &self,
        local_axis: &UnitVector<Real>,
        bias: Real,
        epsilon: Real,
    ) -> SplitResult<Self>;

    fn intersection_with_cuboid(
        &self,
        position: &Isometry<Real>,
        cuboid: &Cuboid,
        cuboid_position: &Isometry<Real>,
        epsilon: Real,
    ) -> Option<Self> {
        self.intersection_with_local_cuboid(cuboid, &position.inv_mul(cuboid_position), epsilon)
    }

    fn intersection_with_local_cuboid(
        &self,
        _cuboid: &Cuboid,
        _cuboid_position: &Isometry<Real>,
        _epsilon: Real,
    ) -> Option<Self> {
        todo!()
    }

    fn intersection_with_aabb(
        &self,
        position: &Isometry<Real>,
        aabb: &AABB,
        epsilon: Real,
    ) -> Option<Self> {
        let cuboid = Cuboid::new(aabb.half_extents());
        let cuboid_pos = Isometry::from(aabb.center());
        self.intersection_with_cuboid(position, &cuboid, &cuboid_pos, epsilon)
    }

    fn intersection_with_local_aabb(&self, aabb: &AABB, epsilon: Real) -> Option<Self> {
        let cuboid = Cuboid::new(aabb.half_extents());
        let cuboid_pos = Isometry::from(aabb.center());
        self.intersection_with_local_cuboid(&cuboid, &cuboid_pos, epsilon)
    }
}

// TODO: impl CanonicalSplit for TriMesh, Polyline, ConvexPolygon, ConvexPolyhedron, Cuboid.
// TODO: impl Split for TriMesh, Polyline, ConvexPolygon, ConvexPolyhedron.
