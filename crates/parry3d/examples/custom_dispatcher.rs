extern crate nalgebra as na;

use parry3d::{
    bounding_volume::{Aabb, BoundingSphere},
    mass_properties::MassProperties,
    math::{Isometry, Point, Vector},
    query::{
        ClosestPoints, Contact, DefaultQueryDispatcher, NonlinearRigidMotion, PointProjection,
        PointQuery, QueryDispatcher, Ray, RayCast, RayIntersection, ShapeCastHit, ShapeCastOptions,
        Unsupported,
    },
    shape::{Ball, Cuboid, FeatureId, Shape, ShapeType, TypedShape},
};

pub struct CustomShape(pub Ball);

fn main() {
    let cube = Cuboid::new(Vector::new(1.0, 1.0, 1.0));
    let custom_shape = CustomShape(Ball::new(1.0));

    let pos12 = Isometry::identity();
    // You can chain your custom dispatcher with the default one.
    let dispatcher = CustomShapeDispatcher.chain(DefaultQueryDispatcher);

    let contact = dispatcher.contact(&pos12, &cube, &custom_shape.0, 0.0);

    _ = dbg!(contact);
}

impl PointQuery for CustomShape {
    fn project_local_point(&self, pt: &Point<f32>, solid: bool) -> PointProjection {
        self.0.project_local_point(pt, solid)
    }

    fn project_local_point_and_get_feature(&self, pt: &Point<f32>) -> (PointProjection, FeatureId) {
        self.0.project_local_point_and_get_feature(pt)
    }
}

impl RayCast for CustomShape {
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_time_of_impact: f32,
        solid: bool,
    ) -> Option<RayIntersection> {
        self.0
            .cast_local_ray_and_get_normal(ray, max_time_of_impact, solid)
    }
}

impl Shape for CustomShape {
    fn compute_local_aabb(&self) -> Aabb {
        self.0.compute_local_aabb()
    }

    fn compute_local_bounding_sphere(&self) -> BoundingSphere {
        self.0.compute_local_bounding_sphere()
    }

    fn clone_dyn(&self) -> Box<dyn Shape> {
        Box::new(Self(self.0))
    }

    fn scale_dyn(&self, scale: &Vector<f32>, num_subdivisions: u32) -> Option<Box<dyn Shape>> {
        Some(self.0.scale_dyn(scale, num_subdivisions)?)
    }

    fn mass_properties(&self, density: f32) -> MassProperties {
        self.0.mass_properties(density)
    }

    fn shape_type(&self) -> ShapeType {
        self.0.shape_type()
    }

    fn as_typed_shape(&self) -> TypedShape {
        self.0.as_typed_shape()
    }

    fn ccd_thickness(&self) -> f32 {
        self.0.ccd_thickness()
    }

    fn ccd_angular_thickness(&self) -> f32 {
        self.0.ccd_angular_thickness()
    }
}

pub struct CustomShapeDispatcher;

impl QueryDispatcher for CustomShapeDispatcher {
    fn intersection_test(
        &self,
        pos12: &Isometry<f32>,
        g1: &dyn Shape,
        g2: &dyn Shape,
    ) -> Result<bool, Unsupported> {
        let (maybe_custom1, maybe_custom2) = (
            g1.downcast_ref::<CustomShape>(),
            g2.downcast_ref::<CustomShape>(),
        );

        match (maybe_custom1, maybe_custom2) {
            (Some(custom1), Some(custom2)) => {
                let p12 = Point::from(pos12.translation.vector);
                return Ok(parry3d::query::details::intersection_test_ball_ball(
                    &p12, &custom1.0, &custom2.0,
                ));
            }
            // TODO: implement the algorithm for each shape type combination you want to support.
            _ => Err(Unsupported),
        }
    }

    fn distance(
        &self,
        pos12: &Isometry<f32>,
        g1: &dyn Shape,
        g2: &dyn Shape,
    ) -> Result<f32, Unsupported> {
        let (maybe_custom1, maybe_custom2) = (
            g1.downcast_ref::<CustomShape>(),
            g2.downcast_ref::<CustomShape>(),
        );

        match (maybe_custom1, maybe_custom2) {
            (Some(custom1), Some(custom2)) => {
                let p2 = Point::from(pos12.translation.vector);
                return Ok(parry3d::query::details::distance_ball_ball(
                    &custom1.0, &p2, &custom2.0,
                ));
            }
            // TODO: implement the algorithm for each shape type combination you want to support.
            _ => Err(Unsupported),
        }
    }

    fn contact(
        &self,
        pos12: &Isometry<f32>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        prediction: f32,
    ) -> Result<Option<Contact>, Unsupported> {
        let (maybe_custom1, maybe_custom2) = (
            g1.downcast_ref::<CustomShape>(),
            g2.downcast_ref::<CustomShape>(),
        );

        match (maybe_custom1, maybe_custom2) {
            (Some(custom1), Some(custom2)) => {
                return Ok(parry3d::query::details::contact_ball_ball(
                    pos12, &custom1.0, &custom2.0, prediction,
                ));
            }
            // TODO: implement the algorithm for each shape type combination you want to support.
            _ => Err(Unsupported),
        }
    }

    fn closest_points(
        &self,
        pos12: &Isometry<f32>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        max_dist: f32,
    ) -> Result<ClosestPoints, Unsupported> {
        let (maybe_custom1, maybe_custom2) = (
            g1.downcast_ref::<CustomShape>(),
            g2.downcast_ref::<CustomShape>(),
        );

        match (maybe_custom1, maybe_custom2) {
            (Some(custom1), Some(custom2)) => {
                return Ok(parry3d::query::details::closest_points_ball_ball(
                    pos12, &custom1.0, &custom2.0, max_dist,
                ));
            }
            // TODO: implement the algorithm for each shape type combination you want to support.
            _ => Err(Unsupported),
        }
    }

    fn cast_shapes(
        &self,
        pos12: &Isometry<f32>,
        local_vel12: &Vector<f32>,
        g1: &dyn Shape,
        g2: &dyn Shape,
        options: ShapeCastOptions,
    ) -> Result<Option<ShapeCastHit>, Unsupported> {
        let (maybe_custom1, maybe_custom2) = (
            g1.downcast_ref::<CustomShape>(),
            g2.downcast_ref::<CustomShape>(),
        );

        match (maybe_custom1, maybe_custom2) {
            (Some(custom1), Some(custom2)) => {
                return Ok(parry3d::query::details::cast_shapes_ball_ball(
                    pos12,
                    local_vel12,
                    &custom1.0,
                    &custom2.0,
                    options,
                ));
            }
            // TODO: implement the algorithm for each shape type combination you want to support.
            _ => Err(Unsupported),
        }
    }

    fn cast_shapes_nonlinear(
        &self,
        motion1: &NonlinearRigidMotion,
        g1: &dyn Shape,
        motion2: &NonlinearRigidMotion,
        g2: &dyn Shape,
        start_time: f32,
        end_time: f32,
        stop_at_penetration: bool,
    ) -> Result<Option<ShapeCastHit>, Unsupported> {
        let (maybe_custom1, maybe_custom2) = (
            g1.downcast_ref::<CustomShape>(),
            g2.downcast_ref::<CustomShape>(),
        );

        match (maybe_custom1, maybe_custom2) {
            (Some(custom1), Some(custom2)) => {
                return parry3d::query::details::cast_shapes_nonlinear(
                    motion1,
                    &custom1.0,
                    motion2,
                    &custom2.0,
                    start_time,
                    end_time,
                    stop_at_penetration,
                );
            }
            // TODO: implement the algorithm for each shape type combination you want to support.
            _ => Err(Unsupported),
        }
    }
}
