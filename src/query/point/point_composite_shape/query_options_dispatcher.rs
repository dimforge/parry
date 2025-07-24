use alloc::boxed::Box;
use core::any::{Any, TypeId};

use crate::query::gjk::GjkOptions;
use crate::query::{DefaultQueryOptions, QueryOptions, QueryOptionsNotUsed};
use crate::shape::Compound;
use hashbrown::HashMap;

/// Retrieves a stored option which should be used for given shape.
pub trait QueryOptionsDispatcher {
    /// Retrieves a stored option which should be used for given shape.
    fn get_option_for_shape(&self, shape_type_id: &TypeId) -> &dyn QueryOptions;
}

impl QueryOptionsDispatcher for () {
    fn get_option_for_shape(&self, _shape_type_id: &TypeId) -> &dyn QueryOptions {
        self
    }
}
impl QueryOptionsDispatcher for DefaultQueryOptions {
    fn get_option_for_shape(&self, _shape_type_id: &TypeId) -> &dyn QueryOptions {
        self
    }
}
impl QueryOptionsDispatcher for QueryOptionsNotUsed {
    fn get_option_for_shape(&self, _shape_type_id: &TypeId) -> &dyn QueryOptions {
        self
    }
}
impl QueryOptionsDispatcher for GjkOptions {
    fn get_option_for_shape(&self, _shape_type_id: &TypeId) -> &dyn QueryOptions {
        self
    }
}

/// Implements [QueryOptionsDispatcher] with a generic dispatch of options.
/// This supports Compound shapes, user-defined shapes and user-defined algorithms.
pub struct QueryOptionsDispatcherMap {
    /// Contains options to apply for shapes query algorithms.
    ///
    /// - Because algorithms can be recursive, Rc is used.
    /// - Because algorithms can contain user-defined shapes, TypeId is used.
    /// - Because algorithms can contain user-defined algorithms (and options), `Box<dyn>` is used.
    pub options_for_shape: HashMap<TypeId, Box<fn(&Self) -> &dyn QueryOptions>>,

    /// Store options inside here.
    /// It's used to capture variables, as [Self::options_for_shape] should be cloneable.
    pub options: HashMap<TypeId, Box<dyn QueryOptions>>,
}

impl Default for QueryOptionsDispatcherMap {
    fn default() -> QueryOptionsDispatcherMap {
        let self_ref: Box<fn(&Self) -> &dyn QueryOptions> =
            Box::new(|s: &QueryOptionsDispatcherMap| -> &dyn QueryOptions { s });
        let gjk_options: Box<fn(&Self) -> &dyn QueryOptions> =
            Box::new(|s: &QueryOptionsDispatcherMap| -> &dyn QueryOptions {
                s.options[&TypeId::of::<GjkOptions>()].as_ref()
            });
        QueryOptionsDispatcherMap {
            options_for_shape: [
                (TypeId::of::<Compound>(), self_ref),
                #[cfg(feature = "dim2")]
                (
                    TypeId::of::<crate::shape::ConvexPolygon>(),
                    gjk_options.clone(),
                ),
                #[cfg(feature = "dim3")]
                (TypeId::of::<crate::shape::ConvexPolyhedron>(), gjk_options),
            ]
            .into(),
            options: [(
                TypeId::of::<GjkOptions>(),
                Box::new(GjkOptions::default()) as Box<dyn QueryOptions>,
            )]
            .into(),
        }
    }
}

impl QueryOptions for QueryOptionsDispatcherMap {
    fn as_any(&self) -> &dyn Any {
        self
    }
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

impl QueryOptionsDispatcherMap {
    /// Retrieves the stored [QueryOptions] of type T.
    pub fn get_option<T: 'static + QueryOptions>(&self) -> Option<&T> {
        let option = self.options.get(&TypeId::of::<T>())?;
        option.as_ref().as_any().downcast_ref::<T>()
    }
    /// Retrieves mutably the stored [QueryOptions] of type T.
    pub fn get_option_mut<T: 'static + QueryOptions>(&mut self) -> Option<&mut T> {
        let option = self.options.get_mut(&TypeId::of::<T>())?;
        option.as_mut().as_any_mut().downcast_mut::<T>()
    }

    /// Downcasts a [QueryOptions] into a [QueryOptionsDispatcher], if it's one, otherwise returns `()`,
    /// which will result in default options being used.
    pub fn from_dyn_or_default(options: &dyn QueryOptions) -> &dyn QueryOptionsDispatcher {
        let options = options
            .as_any()
            .downcast_ref::<QueryOptionsDispatcherMap>()
            .map_or(&() as &dyn QueryOptionsDispatcher, |m| {
                m as &dyn QueryOptionsDispatcher
            });
        options
    }
}

impl QueryOptionsDispatcher for QueryOptionsDispatcherMap {
    fn get_option_for_shape(&self, shape_type_id: &TypeId) -> &dyn QueryOptions {
        let option = self.options_for_shape.get(shape_type_id);
        if let Some(option) = option {
            (*option)(self)
        } else {
            &()
        }
    }
}
