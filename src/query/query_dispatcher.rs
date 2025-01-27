use crate::math::{Isometry, Real, Vector};
use crate::query::details::ShapeCastOptions;
#[cfg(feature = "std")]
use crate::query::{
    contact_manifolds::{ContactManifoldsWorkspace, NormalConstraints},
    ContactManifold,
};
use crate::query::{ClosestPoints, Contact, NonlinearRigidMotion, ShapeCastHit, Unsupported};
use crate::shape::Shape;

#[cfg(feature = "std")]
/// A query dispatcher for queries relying on spatial coherence, including contact-manifold computation.
///
/// Third-party crates adding custom shapes should implement
/// this trait to allow pair-wise queries between the added shape and other shapes.
///
/// The `root_dispatcher` argument is a reference to the root dispatcher of the composite dispatcher.
/// This is necessary to support recursive dispatching for composite shapes.
pub trait PersistentQueryDispatcherComposite<ManifoldData = (), ContactData = ()>:
    QueryDispatcherComposite
{
    /// Compute all the contacts between two shapes.
    ///
    /// The output is written into `manifolds` and `context`. Both can persist
    /// between multiple calls to `contacts` by re-using the result of the previous
    /// call to `contacts`. This persistence can significantly improve collision
    /// detection performances by allowing the underlying algorithms to exploit
    /// spatial and temporal coherence.
    fn contact_manifolds(
        &self,
        root_dispatcher: &dyn PersistentQueryDispatcher<ManifoldData, ContactData>,
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
        normal_constraints1: Option<&dyn NormalConstraints>,
        normal_constraints2: Option<&dyn NormalConstraints>,
        prediction: Real,
        manifold: &mut ContactManifold<ManifoldData, ContactData>,
    ) -> Result<(), Unsupported>;
}

#[cfg(feature = "std")]
/// A query dispatcher for queries relying on spatial coherence, including contact-manifold computation.
///
/// This trait should not be implemented. Instead, implement the [`PersistentQueryDispatcherComposite`] trait.
pub trait PersistentQueryDispatcher<ManifoldData = (), ContactData = ()>: QueryDispatcher {
    /// Compute all the contacts between two shapes.
    ///
    /// The output is written into `manifolds` and `context`. Both can persist
    /// between multiple calls to `contacts` by re-using the result of the previous
    /// call to `contacts`. This persistence can significantly improve collision
    /// detection performances by allowing the underlying algorithms to exploit
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
        normal_constraints1: Option<&dyn NormalConstraints>,
        normal_constraints2: Option<&dyn NormalConstraints>,
        prediction: Real,
        manifold: &mut ContactManifold<ManifoldData, ContactData>,
    ) -> Result<(), Unsupported>;
}

#[cfg(feature = "std")]
impl<T, ManifoldData, ContactData> PersistentQueryDispatcher<ManifoldData, ContactData> for T
where
    T: PersistentQueryDispatcherComposite<ManifoldData, ContactData>,
{
    fn contact_manifolds(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        prediction: Real,
        manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
        workspace: &mut Option<ContactManifoldsWorkspace>,
    ) -> Result<(), Unsupported> {
        <Self as PersistentQueryDispatcherComposite<ManifoldData, ContactData>>::contact_manifolds(
            self, self, pos12, g1, g2, prediction, manifolds, workspace,
        )
    }

    fn contact_manifold_convex_convex(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        normal_constraints1: Option<&dyn NormalConstraints>,
        normal_constraints2: Option<&dyn NormalConstraints>,
        prediction: Real,
        manifold: &mut ContactManifold<ManifoldData, ContactData>,
    ) -> Result<(), Unsupported> {
        <Self as PersistentQueryDispatcherComposite<ManifoldData, ContactData>>::contact_manifold_convex_convex(self, pos12, g1, g2, normal_constraints1, normal_constraints2, prediction, manifold)
    }
}

/// Dispatcher for pairwise queries.
///
/// Third-party crates adding custom shapes should implement
/// this trait to allow pair-wise queries between the added shape and other shapes.
///
/// The `pos12` argument to most queries is the transform from the local space of `g2` to that of
/// `g1`.
///
/// The `root_dispatcher` argument is a reference to the root dispatcher of the composite dispatcher.
/// This is necessary to support recursive dispatching for composite shapes.
pub trait QueryDispatcherComposite: Send + Sync {
    /// Tests whether two shapes are intersecting.
    fn intersection_test(
        &self,
        root_dispatcher: &dyn QueryDispatcher,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
    ) -> Result<bool, Unsupported>;

    /// Computes the minimum distance separating two shapes.
    ///
    /// Returns `0.0` if the objects are touching or penetrating.
    fn distance(
        &self,
        root_dispatcher: &dyn QueryDispatcher,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
    ) -> Result<Real, Unsupported>;

    /// Computes one pair of contact points point between two shapes.
    ///
    /// Returns `None` if the objects are separated by a distance greater than `prediction`.
    fn contact(
        &self,
        root_dispatcher: &dyn QueryDispatcher,
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
        root_dispatcher: &dyn QueryDispatcher,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        max_dist: Real,
    ) -> Result<ClosestPoints, Unsupported>;

    /// Computes the smallest time when two shapes under translational movement are separated by a
    /// distance smaller or equal to `distance`.
    ///
    /// Returns `0.0` if the objects are touching or penetrating.
    ///
    /// # Parameters
    /// - `pos12`: the position of the second shape relative to the first shape.
    /// - `local_vel12`: the relative velocity between the two shapes, expressed in the local-space
    ///                  of the first shape. In other world: `pos1.inverse() * (vel2 - vel1)`.
    /// - `g1`: the first shape involved in the shape-cast.
    /// - `g2`: the second shape involved in the shape-cast.
    /// - `target_dist`: a hit will be returned as soon as the two shapes get closer than `target_dist`.
    /// - `max_time_of_impact`: the maximum allowed travel time. This method returns `None` if the time-of-impact
    ///              detected is theater than this value.
    fn cast_shapes(
        &self,
        root_dispatcher: &dyn QueryDispatcher,
        pos12: &Isometry<Real>,
        local_vel12: &Vector<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        options: ShapeCastOptions,
    ) -> Result<Option<ShapeCastHit>, Unsupported>;

    /// Computes the smallest time of impact of two shapes under translational and rotational movement.
    ///
    /// # Parameters
    /// * `motion1` - The motion of the first shape.
    /// * `g1` - The first shape involved in the query.
    /// * `motion2` - The motion of the second shape.
    /// * `g2` - The second shape involved in the query.
    /// * `start_time` - The starting time of the interval where the motion takes place.
    /// * `end_time` - The end time of the interval where the motion takes place.
    /// * `stop_at_penetration` - If the casted shape starts in a penetration state with any
    ///    collider, two results are possible. If `stop_at_penetration` is `true` then, the
    ///    result will have a `time_of_impact` equal to `start_time`. If `stop_at_penetration` is `false`
    ///    then the nonlinear shape-casting will see if further motion wrt. the penetration normal
    ///    would result in tunnelling. If it does not (i.e. we have a separating velocity along
    ///    that normal) then the nonlinear shape-casting will attempt to find another impact,
    ///    at a time `> start_time` that could result in tunnelling.
    fn cast_shapes_nonlinear(
        &self,
        root_dispatcher: &dyn QueryDispatcher,
        motion1: &NonlinearRigidMotion,
        g1: &dyn Shape,
        motion2: &NonlinearRigidMotion,
        g2: &dyn Shape,
        start_time: Real,
        end_time: Real,
        stop_at_penetration: bool,
    ) -> Result<Option<ShapeCastHit>, Unsupported>;
}

/// Dispatcher for pairwise queries.
///
/// This trait should not be implemented. Instead, implement the [`QueryDispatcherComposite`] trait.
///
/// The `pos12` argument to most queries is the transform from the local space of `g2` to that of
/// `g1`.
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
    ///
    /// # Parameters
    /// - `pos12`: the position of the second shape relative to the first shape.
    /// - `local_vel12`: the relative velocity between the two shapes, expressed in the local-space
    ///                  of the first shape. In other world: `pos1.inverse() * (vel2 - vel1)`.
    /// - `g1`: the first shape involved in the shape-cast.
    /// - `g2`: the second shape involved in the shape-cast.
    /// - `target_dist`: a hit will be returned as soon as the two shapes get closer than `target_dist`.
    /// - `max_time_of_impact`: the maximum allowed travel time. This method returns `None` if the time-of-impact
    ///              detected is theater than this value.
    fn cast_shapes(
        &self,
        pos12: &Isometry<Real>,
        local_vel12: &Vector<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        options: ShapeCastOptions,
    ) -> Result<Option<ShapeCastHit>, Unsupported>;

    /// Construct a `QueryDispatcher` that falls back on `other` for cases not handled by `self`
    fn chain<U: QueryDispatcherComposite>(self, other: U) -> QueryDispatcherChain<Self, U>
    where
        Self: Sized,
    {
        QueryDispatcherChain(self, other)
    }

    /// Computes the smallest time of impact of two shapes under translational and rotational movement.
    ///
    /// # Parameters
    /// * `motion1` - The motion of the first shape.
    /// * `g1` - The first shape involved in the query.
    /// * `motion2` - The motion of the second shape.
    /// * `g2` - The second shape involved in the query.
    /// * `start_time` - The starting time of the interval where the motion takes place.
    /// * `end_time` - The end time of the interval where the motion takes place.
    /// * `stop_at_penetration` - If the casted shape starts in a penetration state with any
    ///    collider, two results are possible. If `stop_at_penetration` is `true` then, the
    ///    result will have a `time_of_impact` equal to `start_time`. If `stop_at_penetration` is `false`
    ///    then the nonlinear shape-casting will see if further motion wrt. the penetration normal
    ///    would result in tunnelling. If it does not (i.e. we have a separating velocity along
    ///    that normal) then the nonlinear shape-casting will attempt to find another impact,
    ///    at a time `> start_time` that could result in tunnelling.
    fn cast_shapes_nonlinear(
        &self,
        motion1: &NonlinearRigidMotion,
        g1: &dyn Shape,
        motion2: &NonlinearRigidMotion,
        g2: &dyn Shape,
        start_time: Real,
        end_time: Real,
        stop_at_penetration: bool,
    ) -> Result<Option<ShapeCastHit>, Unsupported>;
}

impl<T: QueryDispatcherComposite> QueryDispatcher for T {
    fn intersection_test(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
    ) -> Result<bool, Unsupported> {
        <Self as QueryDispatcherComposite>::intersection_test(self, self, pos12, g1, g2)
    }

    fn distance(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
    ) -> Result<Real, Unsupported> {
        <Self as QueryDispatcherComposite>::distance(self, self, pos12, g1, g2)
    }

    fn contact(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        prediction: Real,
    ) -> Result<Option<Contact>, Unsupported> {
        <Self as QueryDispatcherComposite>::contact(self, self, pos12, g1, g2, prediction)
    }

    fn closest_points(
        &self,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        max_dist: Real,
    ) -> Result<ClosestPoints, Unsupported> {
        <Self as QueryDispatcherComposite>::closest_points(self, self, pos12, g1, g2, max_dist)
    }

    fn cast_shapes(
        &self,
        pos12: &Isometry<Real>,
        local_vel12: &Vector<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        options: ShapeCastOptions,
    ) -> Result<Option<ShapeCastHit>, Unsupported> {
        <Self as QueryDispatcherComposite>::cast_shapes(
            self,
            self,
            pos12,
            local_vel12,
            g1,
            g2,
            options,
        )
    }

    fn cast_shapes_nonlinear(
        &self,
        motion1: &NonlinearRigidMotion,
        g1: &dyn Shape,
        motion2: &NonlinearRigidMotion,
        g2: &dyn Shape,
        start_time: Real,
        end_time: Real,
        stop_at_penetration: bool,
    ) -> Result<Option<ShapeCastHit>, Unsupported> {
        <Self as QueryDispatcherComposite>::cast_shapes_nonlinear(
            self,
            self,
            motion1,
            g1,
            motion2,
            g2,
            start_time,
            end_time,
            stop_at_penetration,
        )
    }
}

/// The composition of two dispatchers
#[derive(Clone, Copy, Debug)]
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

impl<T, U> QueryDispatcherComposite for QueryDispatcherChain<T, U>
where
    T: QueryDispatcherComposite,
    U: QueryDispatcherComposite,
{
    chain_method!(intersection_test(
        root_dispatcher: &dyn QueryDispatcher,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
    ) -> bool);

    chain_method!(distance(
        root_dispatcher: &dyn QueryDispatcher, pos12: &Isometry<Real>, g1: &dyn Shape, g2: &dyn Shape,) -> Real);

    chain_method!(contact(
        root_dispatcher: &dyn QueryDispatcher,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        prediction: Real,
    ) -> Option<Contact>);

    chain_method!(closest_points(
        root_dispatcher: &dyn QueryDispatcher,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        max_dist: Real,
    ) -> ClosestPoints);

    chain_method!(cast_shapes(
        root_dispatcher: &dyn QueryDispatcher,
        pos12: &Isometry<Real>,
        vel12: &Vector<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        options: ShapeCastOptions,
    ) -> Option<ShapeCastHit>);

    chain_method!(cast_shapes_nonlinear(
        root_dispatcher: &dyn QueryDispatcher,
        motion1: &NonlinearRigidMotion,
        g1: &dyn Shape,
        motion2: &NonlinearRigidMotion,
        g2: &dyn Shape,
        start_time: Real,
        end_time: Real,
        stop_at_penetration: bool,
    ) -> Option<ShapeCastHit>);
}

#[cfg(feature = "std")]
impl<ManifoldData, ContactData, T, U> PersistentQueryDispatcherComposite<ManifoldData, ContactData>
    for QueryDispatcherChain<T, U>
where
    T: PersistentQueryDispatcherComposite<ManifoldData, ContactData>,
    U: PersistentQueryDispatcherComposite<ManifoldData, ContactData>,
{
    chain_method!(contact_manifolds(
        root_dispatcher: &dyn PersistentQueryDispatcher<ManifoldData, ContactData>,
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        prediction: Real,
        manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
        workspace: &mut Option<ContactManifoldsWorkspace>,
    ) -> ());

    chain_method!(contact_manifold_convex_convex(
        pos12: &Isometry<Real>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        normal_constraints1: Option<&dyn NormalConstraints>,
        normal_constraints2: Option<&dyn NormalConstraints>,
        prediction: Real,
        manifold: &mut ContactManifold<ManifoldData, ContactData>,
    ) -> ());
}
