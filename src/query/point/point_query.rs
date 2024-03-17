use crate::math::*;
use crate::shape::FeatureId;

#[cfg(feature = "rkyv")]
use rkyv::{bytecheck, CheckBytes};

/// Description of the projection of a point on a shape.
#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize, CheckBytes),
    archive(as = "Self")
)]
pub struct PointProjection {
    /// Whether or not the point to project was inside of the shape.
    pub is_inside: bool,
    /// The projection result.
    pub point: Point,
}

impl PointProjection {
    /// Initializes a new `PointProjection`.
    pub fn new(is_inside: bool, point: Point) -> Self {
        PointProjection { is_inside, point }
    }

    /// Transforms `self.point` by `pos`.
    pub fn transform_by(&self, pos: &Isometry) -> Self {
        PointProjection {
            is_inside: self.is_inside,
            point: pos.transform_point(&self.point),
        }
    }
}

/// Trait of objects that can be tested for point inclusion and projection.
pub trait PointQuery {
    /// Projects a point on `self`, unless the projection lies further than the given max distance.
    ///
    /// The point is assumed to be expressed in the local-space of `self`.
    fn project_local_point_with_max_dist(
        &self,
        pt: &Point,
        solid: bool,
        max_dist: Real,
    ) -> Option<PointProjection> {
        let proj = self.project_local_point(pt, solid);
        if distance(proj.point, *pt) > max_dist {
            None
        } else {
            Some(proj)
        }
    }

    /// Projects a point on `self` transformed by `m`, unless the projection lies further than the given max distance.
    fn project_point_with_max_dist(
        &self,
        m: &Isometry,
        pt: &Point,
        solid: bool,
        max_dist: Real,
    ) -> Option<PointProjection> {
        self.project_local_point_with_max_dist(&m.inverse_transform_point(pt), solid, max_dist)
            .map(|proj| proj.transform_by(m))
    }

    /// Projects a point on `self`.
    ///
    /// The point is assumed to be expressed in the local-space of `self`.
    fn project_local_point(&self, pt: &Point, solid: bool) -> PointProjection;

    /// Projects a point on the boundary of `self` and returns the id of the
    /// feature the point was projected on.
    fn project_local_point_and_get_feature(&self, pt: &Point) -> (PointProjection, FeatureId);

    /// Computes the minimal distance between a point and `self`.
    fn distance_to_local_point(&self, pt: &Point, solid: bool) -> Real {
        let proj = self.project_local_point(pt, solid);
        let dist = distance(*pt, proj.point);

        if solid || !proj.is_inside {
            dist
        } else {
            -dist
        }
    }

    /// Tests if the given point is inside of `self`.
    fn contains_local_point(&self, pt: &Point) -> bool {
        self.project_local_point(pt, true).is_inside
    }

    /// Projects a point on `self` transformed by `m`.
    fn project_point(&self, m: &Isometry, pt: &Point, solid: bool) -> PointProjection {
        self.project_local_point(&m.inverse_transform_point(pt), solid)
            .transform_by(m)
    }

    /// Computes the minimal distance between a point and `self` transformed by `m`.
    #[inline]
    fn distance_to_point(&self, m: &Isometry, pt: &Point, solid: bool) -> Real {
        self.distance_to_local_point(&m.inverse_transform_point(pt), solid)
    }

    /// Projects a point on the boundary of `self` transformed by `m` and returns the id of the
    /// feature the point was projected on.
    fn project_point_and_get_feature(
        &self,
        m: &Isometry,
        pt: &Point,
    ) -> (PointProjection, FeatureId) {
        let res = self.project_local_point_and_get_feature(&m.inverse_transform_point(pt));
        (res.0.transform_by(m), res.1)
    }

    /// Tests if the given point is inside of `self` transformed by `m`.
    #[inline]
    fn contains_point(&self, m: &Isometry, pt: &Point) -> bool {
        self.contains_local_point(&m.inverse_transform_point(pt))
    }
}

/// Returns shape-specific info in addition to generic projection information
///
/// One requirement for the `PointQuery` trait is to be usable as a trait
/// object. Unfortunately this precludes us from adding an associated type to it
/// that might allow us to return shape-specific information in addition to the
/// general information provided in `PointProjection`. This is where
/// `PointQueryWithLocation` comes in. It forgoes the ability to be used as a trait
/// object in exchange for being able to provide shape-specific projection
/// information.
///
/// Any shapes that implement `PointQuery` but are able to provide extra
/// information, can implement `PointQueryWithLocation` in addition and have their
/// `PointQuery::project_point` implementation just call out to
/// `PointQueryWithLocation::project_point_and_get_location`.
pub trait PointQueryWithLocation {
    /// Additional shape-specific projection information
    ///
    /// In addition to the generic projection information returned in
    /// `PointProjection`, implementations might provide shape-specific
    /// projection info. The type of this shape-specific information is defined
    /// by this associated type.
    type Location;

    /// Projects a point on `self`.
    fn project_local_point_and_get_location(
        &self,
        pt: &Point,
        solid: bool,
    ) -> (PointProjection, Self::Location);

    /// Projects a point on `self` transformed by `m`.
    fn project_point_and_get_location(
        &self,
        m: &Isometry,
        pt: &Point,
        solid: bool,
    ) -> (PointProjection, Self::Location) {
        let res = self.project_local_point_and_get_location(&m.inverse_transform_point(pt), solid);
        (res.0.transform_by(m), res.1)
    }

    /// Projects a point on `self`, with a maximum projection distance.
    fn project_local_point_and_get_location_with_max_dist(
        &self,
        pt: &Point,
        solid: bool,
        max_dist: Real,
    ) -> Option<(PointProjection, Self::Location)> {
        let (proj, location) = self.project_local_point_and_get_location(pt, solid);
        if distance(proj.point, *pt) > max_dist {
            None
        } else {
            Some((proj, location))
        }
    }

    /// Projects a point on `self` transformed by `m`, with a maximum projection distance.
    fn project_point_and_get_location_with_max_dist(
        &self,
        m: &Isometry,
        pt: &Point,
        solid: bool,
        max_dist: Real,
    ) -> Option<(PointProjection, Self::Location)> {
        self.project_local_point_and_get_location_with_max_dist(
            &m.inverse_transform_point(pt),
            solid,
            max_dist,
        )
        .map(|res| (res.0.transform_by(m), res.1))
    }
}
