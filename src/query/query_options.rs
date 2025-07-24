use core::any::Any;

/// Options to pass to [QueryOptions]. The type to be passed depends on each shape implementation.
/// If you're not sure, use `&()` which will evaluate to default options for included shapes in
/// [QueryOptions] implementations.
///
/// # See Also
/// - [GjkOptions][crate::query::gjk::GjkOptions]
/// - [QueryOptionsDispatcher][crate::query::point::QueryOptionsDispatcher]
///   - [QueryOptionsDispatcherMap][crate::query::point::QueryOptionsDispatcherMap]
pub trait QueryOptions {
    /// Downcast to [Any] to be compatible with [PointQuery][crate::query::PointQuery] or [RayCast][crate::query::RayCast] query options.
    fn as_any(&self) -> &dyn Any;
    /// Downcast to [Any] mutably.
    fn as_any_mut(&mut self) -> &mut dyn Any;
}

/// Unit type used to communicate that an option is not used downstream.
///
/// Pass this to a method awaiting `&dyn QueryOptions` that we know doesn't use any algorithm options,
/// so the flow is easier to understand.
///
/// Its presence should be challenged any time a new algorithm is modified or a new option is introduced.
pub struct QueryOptionsNotUsed;

impl QueryOptions for QueryOptionsNotUsed {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}
/// Unit type to help with passing default options, it can be challenging to know which [QueryOptions]
/// implementation should be passed to a function so this is an easy default.
///
/// Note that you can also pass `&()` in such cases.
pub struct DefaultQueryOptions;

impl QueryOptions for DefaultQueryOptions {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

impl QueryOptions for () {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}
