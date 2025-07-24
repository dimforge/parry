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
    fn as_any(&self) -> &dyn Any;
    fn as_any_mut(&mut self) -> &mut dyn Any;
}

/// Type alias used to communicate that an option is not used downstream.
///
/// Pass this to a method awaiting `&dyn QueryOptions` that we know doesn't use any algorithm options,
/// so the flow is easier to understand.
///
/// Its presence should be challenged any time a new algorithm is modified or a new option is introduced.
pub use DefaultQueryOptions as QueryOptionsNotUsed;

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
