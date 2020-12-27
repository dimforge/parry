use crate::math::{Isometry, Point, Real, Vector};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
/// The type local linear approximation of the neighborhood of a pair contact points on two shapes
pub enum KinematicsCategory {
    /// Both neighborhoods are assimilated to a single point.
    PointPoint,
    /// The first shape's neighborhood at the contact point is assimilated to a plane while
    /// the second is assimilated to a point.
    PlanePoint,
}

#[derive(Copy, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
/// Local contact geometry at the neighborhood of a pair of contact points.
pub struct ContactKinematics {
    /// The local contact geometry.
    pub category: KinematicsCategory,
    /// The dilation applied to the first contact geometry.
    pub radius1: Real,
    /// The dilation applied to the second contact geometry.
    pub radius2: Real,
}

impl Default for ContactKinematics {
    fn default() -> Self {
        ContactKinematics {
            category: KinematicsCategory::PointPoint,
            radius1: 0.0,
            radius2: 0.0,
        }
    }
}

#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
/// A single contact between two shape.
pub struct TrackedContact<Data> {
    /// The contact point in the local-space of the first shape.
    pub local_p1: Point<Real>,
    /// The contact point in the local-space of the second shape.
    pub local_p2: Point<Real>,
    /// The identifier of the subshape of the first shape involved in this contact.
    ///
    /// For primitive shapes like cuboid, ball, etc., this is 0.
    /// For shapes like trimesh and heightfield this identifies the specific triangle
    /// involved in the contact.
    pub fid1: u8,
    /// The identifier of the subshape of the second shape involved in this contact.
    ///
    /// For primitive shapes like cuboid, ball, etc., this is 0.
    /// For shapes like trimesh and heightfield this identifies the specific triangle
    /// involved in the contact.
    pub fid2: u8,
    /// The distance between the two shapes along the contact normal.
    ///
    /// If this is negative, the shapes are penetrating.
    pub dist: Real,
    /// Data associated to this contact.
    pub data: Data,
}

impl<Data: Default + Copy> TrackedContact<Data> {
    pub fn new(
        local_p1: Point<Real>,
        local_p2: Point<Real>,
        fid1: u8,
        fid2: u8,
        dist: Real,
    ) -> Self {
        Self {
            local_p1,
            local_p2,
            fid1,
            fid2,
            dist,
            data: Data::default(),
        }
    }

    pub fn flipped(
        local_p1: Point<Real>,
        local_p2: Point<Real>,
        fid1: u8,
        fid2: u8,
        dist: Real,
        flipped: bool,
    ) -> Self {
        if !flipped {
            Self::new(local_p1, local_p2, fid1, fid2, dist)
        } else {
            Self::new(local_p2, local_p1, fid2, fid1, dist)
        }
    }

    pub fn copy_geometry_from(&mut self, contact: Self) {
        self.local_p1 = contact.local_p1;
        self.local_p2 = contact.local_p2;
        self.fid1 = contact.fid1;
        self.fid2 = contact.fid2;
        self.dist = contact.dist;
    }
}

#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
/// A contact manifold between two shapes.
///
/// A contact manifold describes a set of contacts between two shapes. All the contact
/// part of the same contact manifold share the same contact normal and contact kinematics.
pub struct ContactManifold<ManifoldData, ContactData> {
    // NOTE: use a SmallVec instead?
    // And for 2D use an ArrayVec since there will never be more than 2 contacts anyways.
    #[cfg(feature = "dim2")]
    pub points: arrayvec::ArrayVec<[TrackedContact<ContactData>; 2]>,
    #[cfg(feature = "dim3")]
    pub points: Vec<TrackedContact<ContactData>>,
    /// The number of active contacts on this contact manifold.
    ///
    /// Active contacts are these that may result in contact forces.
    pub num_active_contacts: usize,
    /// The contact normal of all the contacts of this manifold, expressed in the local space of the first shape.
    pub local_n1: Vector<Real>,
    /// The contact normal of all the contacts of this manifold, expressed in the local space of the second shape.
    pub local_n2: Vector<Real>,
    /// The contact kinematics of all the contacts of this manifold.
    pub kinematics: ContactKinematics,
    /// The pair of subshapes involved in this contact manifold.
    pub subshape_index_pair: (usize, usize),
    /// Additional tracked data associated to this contact manifold.
    pub data: ManifoldData,
}

impl<ManifoldData, ContactData: Default + Copy> ContactManifold<ManifoldData, ContactData> {
    pub fn new() -> Self
    where
        ManifoldData: Default,
    {
        Self::with_data((0, 0), ManifoldData::default())
    }

    pub fn with_data(subshape_index_pair: (usize, usize), data: ManifoldData) -> Self {
        Self {
            #[cfg(feature = "dim2")]
            points: arrayvec::ArrayVec::new(),
            #[cfg(feature = "dim3")]
            points: Vec::new(),
            num_active_contacts: 0,
            local_n1: Vector::zeros(),
            local_n2: Vector::zeros(),
            subshape_index_pair,
            kinematics: ContactKinematics::default(),
            data,
        }
    }

    pub fn take(&mut self) -> Self
    where
        ManifoldData: Clone,
    {
        ContactManifold {
            #[cfg(feature = "dim2")]
            points: self.points.clone(),
            #[cfg(feature = "dim3")]
            points: std::mem::replace(&mut self.points, Vec::new()),
            num_active_contacts: self.num_active_contacts,
            local_n1: self.local_n1,
            local_n2: self.local_n2,
            kinematics: self.kinematics,
            subshape_index_pair: self.subshape_index_pair,
            data: self.data.clone(),
        }
    }

    /*
    pub(crate) fn single_manifold<'a, 'b>(
        manifolds: &mut Vec<Self>,
        data: &dyn Fn() -> ManifoldData,
    ) -> &'a mut Self {
        if manifolds.len() == 0 {
            let manifold_data = data();
            manifolds.push(ContactManifold::with_data((0, 0), manifold_data));
        }

        &mut manifolds[0]
    }
    */

    /// Number of active contacts on this contact manifold.
    #[inline]
    pub fn num_active_contacts(&self) -> usize {
        self.num_active_contacts
    }

    /// The slice of all the active contacts on this contact manifold.
    ///
    /// Active contacts are contacts that may end up generating contact forces.
    #[inline]
    pub fn active_contacts(&self) -> &[TrackedContact<ContactData>] {
        &self.points[..self.num_active_contacts]
    }

    #[inline]
    pub fn active_contacts_mut(&mut self) -> &mut [TrackedContact<ContactData>] {
        &mut self.points[..self.num_active_contacts]
    }

    /// The slice of all the contacts, active or not, on this contact manifold.
    #[inline]
    pub fn all_contacts(&self) -> &[TrackedContact<ContactData>] {
        &self.points
    }

    #[inline]
    pub fn try_update_contacts(&mut self, pos12: &Isometry<Real>) -> bool {
        //        const DOT_THRESHOLD: Real = 0.crate::COS_10_DEGREES;
        const DOT_THRESHOLD: Real = crate::utils::COS_5_DEGREES;
        const DIST_SQ_THRESHOLD: Real = 0.001; // FIXME: this should not be hard-coded.
        self.try_update_contacts_eps(pos12, DOT_THRESHOLD, DIST_SQ_THRESHOLD)
    }

    #[inline]
    pub fn try_update_contacts_eps(
        &mut self,
        pos12: &Isometry<Real>,
        angle_dot_threshold: Real,
        dist_sq_threshold: Real,
    ) -> bool {
        if self.points.len() == 0 {
            return false;
        }

        let local_n2 = pos12 * self.local_n2;

        if -self.local_n1.dot(&local_n2) < angle_dot_threshold {
            return false;
        }

        for pt in &mut self.points {
            let local_p2 = pos12 * pt.local_p2;
            let dpt = local_p2 - pt.local_p1;
            let dist = dpt.dot(&self.local_n1);

            if dist * pt.dist < 0.0 {
                // We switched between penetrating/non-penetrating.
                // The may result in other contacts to appear.
                return false;
            }
            let new_local_p1 = local_p2 - self.local_n1 * dist;

            if na::distance_squared(&pt.local_p1, &new_local_p1) > dist_sq_threshold {
                return false;
            }

            pt.dist = dist;
            pt.local_p1 = new_local_p1;
        }

        true
    }

    pub fn match_contacts(&mut self, old_contacts: &[TrackedContact<ContactData>]) {
        for contact in &mut self.points {
            for old_contact in old_contacts {
                if contact.fid1 == old_contact.fid1 && contact.fid2 == old_contact.fid2 {
                    // Transfer the tracked data.
                    contact.data = old_contact.data;
                }
            }
        }
    }

    pub fn match_contacts_using_positions(
        &mut self,
        old_contacts: &[TrackedContact<ContactData>],
        dist_threshold: Real,
    ) {
        let sq_threshold = dist_threshold * dist_threshold;
        for contact in &mut self.points {
            for old_contact in old_contacts {
                if na::distance_squared(&contact.local_p1, &old_contact.local_p1) < sq_threshold
                    && na::distance_squared(&contact.local_p2, &old_contact.local_p2) < sq_threshold
                {
                    // Transfer the tracked data.
                    contact.data = old_contact.data;
                }
            }
        }
    }

    pub fn clear(&mut self) {
        self.points.clear();
        self.num_active_contacts = 0;
    }

    pub fn find_deepest_contact(&self) -> Option<&TrackedContact<ContactData>> {
        let mut deepest = self.points.get(0)?;

        for pt in &self.points {
            if pt.dist < deepest.dist {
                deepest = pt;
            }
        }

        Some(deepest)
    }

    /// Sort the contacts of this contact manifold such that the active contacts are in the first
    /// positions of the array.
    #[inline]
    pub fn sort_contacts(&mut self, prediction_distance: Real) {
        let num_contacts = self.points.len();
        match num_contacts {
            0 => {
                self.num_active_contacts = 0;
            }
            1 => {
                self.num_active_contacts = (self.points[0].dist < prediction_distance) as usize;
            }
            _ => {
                let mut first_inactive_index = num_contacts;

                self.num_active_contacts = 0;
                while self.num_active_contacts != first_inactive_index {
                    if self.points[self.num_active_contacts].dist >= prediction_distance {
                        // Swap with the last contact.
                        self.points
                            .swap(self.num_active_contacts, first_inactive_index - 1);
                        first_inactive_index -= 1;
                    } else {
                        self.num_active_contacts += 1;
                    }
                }
            }
        }
    }
}
