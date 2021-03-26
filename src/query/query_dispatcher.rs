use crate::math::{Isometry, Real, Vector};
use crate::query::contact_manifolds::ContactManifoldsWorkspace;
use crate::query::{
    ClosestPoints, Contact, ContactManifold, NonlinearRigidMotion, Unsupported, TOI,
};
use crate::shape::Shape;

/// A query dispatcher for queries relying on spatial coherence, including contact-manifold computation.
pub trait PersistentQueryDispatcher<ManifoldData = (), ContactData = ()>: QueryDispatcher {
    /// Compute all the contacts between two shapes.
    ///
    /// The output is written into `manifolds` and `context`. Both can persist
    /// between multiple calls to `contacts` by re-using the result of the previous
    /// call to `contacts`. This persistence can significantly improve collision
    /// detection performances by allowing the underlying algorithms to exploid
    /// spatial and temporal coherence.
    fn contact_manifolds(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        prediction: Real,
        manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
        workspace: &mut Option<ContactManifoldsWorkspace>,
    ) -> Result<(), Unsupported>;

    /// Computes the contact-manifold between two convex shapes.
    fn contact_manifold_convex_convex(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        prediction: Real,
        manifold: &mut ContactManifold<ManifoldData, ContactData>,
    ) -> Result<(), Unsupported>;
}

/// Dispatcher for pairwise queries.
///
/// Custom implementations allow crates that support an abstract `QueryDispatcher` to handle custom
/// shapes.
pub trait QueryDispatcher: Send + Sync {
    /// Tests whether two shapes are intersecting.
    fn intersection_test(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
    ) -> Result<bool, Unsupported>;

    /// Computes the minimum distance separating two shapes.
    ///
    /// Returns `0.0` if the objects are touching or penetrating.
    fn distance(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
    ) -> Result<Real, Unsupported>;

    /// Computes one pair of contact points point between two shapes.
    ///
    /// Returns `None` if the objects are separated by a distance greater than `prediction`.
    fn contact(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        prediction: Real,
    ) -> Result<Option<Contact>, Unsupported>;

    /// Computes the pair of closest points between two shapes.
    ///
    /// Returns `ClosestPoints::Disjoint` if the objects are separated by a distance greater than `max_dist`.
    fn closest_points(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        max_dist: Real,
    ) -> Result<ClosestPoints, Unsupported>;

    /// Computes the smallest time when two shapes under translational movement are separated by a
    /// distance smaller or equal to `distance`.
    ///
    /// Returns `0.0` if the objects are touching or penetrating.
    fn time_of_impact(
        &self,
        pos12: &Isometry<Real>,
        vel12: &Vector<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        max_toi: Real,
        target_distance: Real,
    ) -> Result<Option<TOI>, Unsupported>;

    /// Construct a `QueryDispatcher` that falls back on `other` for cases not handled by `self`
    fn chain<U: QueryDispatcher>(self, other: U) -> QueryDispatcherChain<Self, U>
    where
        Self: Sized,
    {
        QueryDispatcherChain(self, other)
    }

    /// Computes the smallest time of impact of two shapes under translational movement.
    fn nonlinear_time_of_impact(
        &self,
        motion1: &NonlinearRigidMotion,
        g1: &dyn Shape,
        motion2: &NonlinearRigidMotion,
        g2: &dyn Shape,
        start_time: Real,
        end_time: Real,
        stop_at_penetration: bool,
    ) -> Result<Option<TOI>, Unsupported>;
}

/// The composition of two dispatchers
pub struct QueryDispatcherChain<T, U>(T, U);

macro_rules! chain_method {
    ($name:ident ( $( $arg:ident : $ty:ty,)*) -> $result:ty) => {
        fn $name(&self, $($arg : $ty,)*
        ) -> Result<$result, Unsupported> {
            (self.0).$name($($arg,)*)
                .or_else(|Unsupported| (self.1).$name($($arg,)*))
        }
    }
}

impl<T, U> QueryDispatcher for QueryDispatcherChain<T, U>
where
    T: QueryDispatcher,
    U: QueryDispatcher,
{
    chain_method!(intersection_test(
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
    ) -> bool);

    chain_method!(distance(pos12: &Isometry<Real>, g1: &dyn Shape, g2: &dyn Shape,) -> Real);

    chain_method!(contact(
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        prediction: Real,
    ) -> Option<Contact>);

    chain_method!(closest_points(
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        max_dist: Real,
    ) -> ClosestPoints);

    chain_method!(time_of_impact(
        pos12: &Isometry<Real>,
        vel12: &Vector<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        max_toi: Real,
        target_distance: Real,
    ) -> Option<TOI>);

    chain_method!(nonlinear_time_of_impact(
        motion1: &NonlinearRigidMotion,
        g1: &dyn Shape,
        motion2: &NonlinearRigidMotion,
        g2: &dyn Shape,
        start_time: Real,
        end_time: Real,
        stop_at_penetration: bool,
    ) -> Option<TOI>);
}
