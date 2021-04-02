use crate::math::{Isometry, Point, Real, Vector};
use crate::query::contact_manifolds::ContactManifoldsWorkspace;
use crate::query::query_dispatcher::PersistentQueryDispatcher;
use crate::query::{
    self, details::NonlinearTOIMode, ClosestPoints, Contact, ContactManifold, NonlinearRigidMotion,
    QueryDispatcher, Unsupported, TOI,
};
use crate::shape::{HalfSpace, Segment, Shape, ShapeType};

/// A dispatcher that exposes built-in queries
#[derive(Debug, Clone)]
pub struct DefaultQueryDispatcher;

impl QueryDispatcher for DefaultQueryDispatcher {
    fn intersection_test(
        &self,
        pos12: &Isometry<Real>,
        shape1: &dyn Shape,
        shape2: &dyn Shape,
    ) -> Result<bool, Unsupported> {
        if let (Some(b1), Some(b2)) = (shape1.as_ball(), shape2.as_ball()) {
            let p12 = Point::from(pos12.translation.vector);
            Ok(query::details::intersection_test_ball_ball(&p12, b1, b2))
        } else if let (Some(c1), Some(c2)) = (shape1.as_cuboid(), shape2.as_cuboid()) {
            Ok(query::details::intersection_test_cuboid_cuboid(
                &pos12, c1, c2,
            ))
        } else if let (Some(t1), Some(c2)) = (shape1.as_triangle(), shape2.as_cuboid()) {
            Ok(query::details::intersection_test_triangle_cuboid(
                &pos12, t1, c2,
            ))
        } else if let (Some(c1), Some(t2)) = (shape1.as_cuboid(), shape2.as_triangle()) {
            Ok(query::details::intersection_test_cuboid_triangle(
                &pos12, c1, t2,
            ))
        } else if let Some(b1) = shape1.as_ball() {
            Ok(query::details::intersection_test_ball_point_query(
                &pos12, b1, shape2,
            ))
        } else if let Some(b2) = shape2.as_ball() {
            Ok(query::details::intersection_test_point_query_ball(
                &pos12, shape1, b2,
            ))
        } else if let (Some(p1), Some(s2)) =
            (shape1.as_shape::<HalfSpace>(), shape2.as_support_map())
        {
            Ok(query::details::intersection_test_halfspace_support_map(
                pos12, p1, s2,
            ))
        } else if let (Some(s1), Some(p2)) =
            (shape1.as_support_map(), shape2.as_shape::<HalfSpace>())
        {
            Ok(query::details::intersection_test_support_map_halfspace(
                pos12, s1, p2,
            ))
        } else if let (Some(s1), Some(s2)) = (shape1.as_support_map(), shape2.as_support_map()) {
            Ok(query::details::intersection_test_support_map_support_map(
                pos12, s1, s2,
            ))
        } else if let Some(c1) = shape1.as_composite_shape() {
            Ok(query::details::intersection_test_composite_shape_shape(
                self, pos12, c1, shape2,
            ))
        } else if let Some(c2) = shape2.as_composite_shape() {
            Ok(query::details::intersection_test_shape_composite_shape(
                self, pos12, shape1, c2,
            ))
        } else {
            Err(Unsupported)
        }
    }

    /// Computes the minimum distance separating two shapes.
    ///
    /// Returns `0.0` if the objects are touching or penetrating.
    fn distance(
        &self,
        pos12: &Isometry<Real>,
        shape1: &dyn Shape,
        shape2: &dyn Shape,
    ) -> Result<Real, Unsupported> {
        let ball1 = shape1.as_ball();
        let ball2 = shape2.as_ball();

        if let (Some(b1), Some(b2)) = (ball1, ball2) {
            let p2 = Point::from(pos12.translation.vector);
            Ok(query::details::distance_ball_ball(b1, &p2, b2))
        } else if let (Some(b1), true) = (ball1, shape2.is_convex()) {
            Ok(query::details::distance_ball_convex_polyhedron(
                pos12, b1, shape2,
            ))
        } else if let (true, Some(b2)) = (shape1.is_convex(), ball2) {
            Ok(query::details::distance_convex_polyhedron_ball(
                pos12, shape1, b2,
            ))
        } else if let (Some(c1), Some(c2)) = (shape1.as_cuboid(), shape2.as_cuboid()) {
            Ok(query::details::distance_cuboid_cuboid(pos12, c1, c2))
        } else if let (Some(s1), Some(s2)) = (shape1.as_segment(), shape2.as_segment()) {
            Ok(query::details::distance_segment_segment(pos12, s1, s2))
        } else if let (Some(p1), Some(s2)) =
            (shape1.as_shape::<HalfSpace>(), shape2.as_support_map())
        {
            Ok(query::details::distance_halfspace_support_map(
                pos12, p1, s2,
            ))
        } else if let (Some(s1), Some(p2)) =
            (shape1.as_support_map(), shape2.as_shape::<HalfSpace>())
        {
            Ok(query::details::distance_support_map_halfspace(
                pos12, s1, p2,
            ))
        } else if let (Some(s1), Some(s2)) = (shape1.as_support_map(), shape2.as_support_map()) {
            Ok(query::details::distance_support_map_support_map(
                pos12, s1, s2,
            ))
        } else if let Some(c1) = shape1.as_composite_shape() {
            Ok(query::details::distance_composite_shape_shape(
                self, pos12, c1, shape2,
            ))
        } else if let Some(c2) = shape2.as_composite_shape() {
            Ok(query::details::distance_shape_composite_shape(
                self, pos12, shape1, c2,
            ))
        } else {
            Err(Unsupported)
        }
    }

    fn contact(
        &self,
        pos12: &Isometry<Real>,
        shape1: &dyn Shape,
        shape2: &dyn Shape,
        prediction: Real,
    ) -> Result<Option<Contact>, Unsupported> {
        let ball1 = shape1.as_ball();
        let ball2 = shape2.as_ball();

        if let (Some(b1), Some(b2)) = (ball1, ball2) {
            Ok(query::details::contact_ball_ball(pos12, b1, b2, prediction))
        // } else if let (Some(c1), Some(c2)) = (shape1.as_cuboid(), shape2.as_cuboid()) {
        //     Ok(query::details::contact_cuboid_cuboid(
        //         pos12, c1, c2, prediction,
        //     ))
        } else if let (Some(p1), Some(s2)) =
            (shape1.as_shape::<HalfSpace>(), shape2.as_support_map())
        {
            Ok(query::details::contact_halfspace_support_map(
                pos12, p1, s2, prediction,
            ))
        } else if let (Some(s1), Some(p2)) =
            (shape1.as_support_map(), shape2.as_shape::<HalfSpace>())
        {
            Ok(query::details::contact_support_map_halfspace(
                pos12, s1, p2, prediction,
            ))
        } else if let (Some(b1), true) = (ball1, shape2.is_convex()) {
            Ok(query::details::contact_ball_convex_polyhedron(
                pos12, b1, shape2, prediction,
            ))
        } else if let (true, Some(b2)) = (shape1.is_convex(), ball2) {
            Ok(query::details::contact_convex_polyhedron_ball(
                pos12, shape1, b2, prediction,
            ))
        } else if let (Some(s1), Some(s2)) = (shape1.as_support_map(), shape2.as_support_map()) {
            Ok(query::details::contact_support_map_support_map(
                pos12, s1, s2, prediction,
            ))
        } else if let Some(c1) = shape1.as_composite_shape() {
            Ok(query::details::contact_composite_shape_shape(
                self, pos12, c1, shape2, prediction,
            ))
        } else if let Some(c2) = shape2.as_composite_shape() {
            Ok(query::details::contact_shape_composite_shape(
                self, pos12, shape1, c2, prediction,
            ))
        } else {
            Err(Unsupported)
        }
    }

    fn closest_points(
        &self,
        pos12: &Isometry<Real>,
        shape1: &dyn Shape,
        shape2: &dyn Shape,
        max_dist: Real,
    ) -> Result<ClosestPoints, Unsupported> {
        let ball1 = shape1.as_ball();
        let ball2 = shape2.as_ball();

        if let (Some(b1), Some(b2)) = (ball1, ball2) {
            Ok(query::details::closest_points_ball_ball(
                &pos12, b1, b2, max_dist,
            ))
        } else if let (Some(b1), true) = (ball1, shape2.is_convex()) {
            Ok(query::details::closest_points_ball_convex_polyhedron(
                pos12, b1, shape2, max_dist,
            ))
        } else if let (true, Some(b2)) = (shape1.is_convex(), ball2) {
            Ok(query::details::closest_points_convex_polyhedron_ball(
                pos12, shape1, b2, max_dist,
            ))
        } else if let (Some(s1), Some(s2)) =
            (shape1.as_shape::<Segment>(), shape2.as_shape::<Segment>())
        {
            Ok(query::details::closest_points_segment_segment(
                &pos12, s1, s2, max_dist,
            ))
        // } else if let (Some(c1), Some(c2)) = (shape1.as_cuboid(), shape2.as_cuboid()) {
        //     Ok(query::details::closest_points_cuboid_cuboid(
        //         pos12, c1, c2, max_dist,
        //     ))
        } else if let (Some(s1), Some(s2)) = (shape1.as_segment(), shape2.as_segment()) {
            Ok(query::details::closest_points_segment_segment(
                pos12, s1, s2, max_dist,
            ))
        // } else if let (Some(c1), Some(t2)) = (shape1.as_cuboid(), shape2.as_triangle()) {
        //     Ok(query::details::closest_points_cuboid_triangle(
        //         pos12, c1, t2, max_dist,
        //     ))
        } else if let (Some(t1), Some(c2)) = (shape1.as_triangle(), shape2.as_cuboid()) {
            Ok(query::details::closest_points_triangle_cuboid(
                pos12, t1, c2, max_dist,
            ))
        } else if let (Some(p1), Some(s2)) =
            (shape1.as_shape::<HalfSpace>(), shape2.as_support_map())
        {
            Ok(query::details::closest_points_halfspace_support_map(
                &pos12, p1, s2, max_dist,
            ))
        } else if let (Some(s1), Some(p2)) =
            (shape1.as_support_map(), shape2.as_shape::<HalfSpace>())
        {
            Ok(query::details::closest_points_support_map_halfspace(
                &pos12, s1, p2, max_dist,
            ))
        } else if let (Some(s1), Some(s2)) = (shape1.as_support_map(), shape2.as_support_map()) {
            Ok(query::details::closest_points_support_map_support_map(
                &pos12, s1, s2, max_dist,
            ))
        } else if let Some(c1) = shape1.as_composite_shape() {
            Ok(query::details::closest_points_composite_shape_shape(
                self, &pos12, c1, shape2, max_dist,
            ))
        } else if let Some(c2) = shape2.as_composite_shape() {
            Ok(query::details::closest_points_shape_composite_shape(
                self, &pos12, shape1, c2, max_dist,
            ))
        } else {
            Err(Unsupported)
        }
    }

    fn time_of_impact(
        &self,
        pos12: &Isometry<Real>,
        local_vel12: &Vector<Real>,
        shape1: &dyn Shape,
        shape2: &dyn Shape,
        max_toi: Real,
    ) -> Result<Option<TOI>, Unsupported> {
        if let (Some(b1), Some(b2)) = (shape1.as_ball(), shape2.as_ball()) {
            Ok(query::details::time_of_impact_ball_ball(
                pos12,
                local_vel12,
                b1,
                b2,
                max_toi,
            ))
        } else if let (Some(p1), Some(s2)) =
            (shape1.as_shape::<HalfSpace>(), shape2.as_support_map())
        {
            Ok(query::details::time_of_impact_halfspace_support_map(
                pos12,
                local_vel12,
                p1,
                s2,
                max_toi,
            ))
        } else if let (Some(s1), Some(p2)) =
            (shape1.as_support_map(), shape2.as_shape::<HalfSpace>())
        {
            Ok(query::details::time_of_impact_support_map_halfspace(
                pos12,
                local_vel12,
                s1,
                p2,
                max_toi,
            ))
        } else if let (Some(s1), Some(s2)) = (shape1.as_support_map(), shape2.as_support_map()) {
            Ok(query::details::time_of_impact_support_map_support_map(
                pos12,
                local_vel12,
                s1,
                s2,
                max_toi,
            ))
        } else if let Some(c1) = shape1.as_composite_shape() {
            Ok(query::details::time_of_impact_composite_shape_shape(
                self,
                pos12,
                local_vel12,
                c1,
                shape2,
                max_toi,
            ))
        } else if let Some(c2) = shape2.as_composite_shape() {
            Ok(query::details::time_of_impact_shape_composite_shape(
                self,
                pos12,
                local_vel12,
                shape1,
                c2,
                max_toi,
            ))
        } else {
            Err(Unsupported)
        }
    }

    fn nonlinear_time_of_impact(
        &self,
        motion1: &NonlinearRigidMotion,
        shape1: &dyn Shape,
        motion2: &NonlinearRigidMotion,
        shape2: &dyn Shape,
        start_time: Real,
        end_time: Real,
        stop_at_penetration: bool,
    ) -> Result<Option<TOI>, Unsupported> {
        if let (Some(sm1), Some(sm2)) = (shape1.as_support_map(), shape2.as_support_map()) {
            let mode = if stop_at_penetration {
                NonlinearTOIMode::StopAtPenetration
            } else {
                NonlinearTOIMode::directional_toi(shape1, shape2)
            };

            Ok(
                query::details::nonlinear_time_of_impact_support_map_support_map(
                    self, motion1, sm1, shape1, motion2, sm2, shape2, start_time, end_time, mode,
                ),
            )
        } else if let Some(c1) = shape1.as_composite_shape() {
            Ok(
                query::details::nonlinear_time_of_impact_composite_shape_shape(
                    self,
                    motion1,
                    c1,
                    motion2,
                    shape2,
                    start_time,
                    end_time,
                    stop_at_penetration,
                ),
            )
        } else if let Some(c2) = shape2.as_composite_shape() {
            Ok(
                query::details::nonlinear_time_of_impact_shape_composite_shape(
                    self,
                    motion1,
                    shape1,
                    motion2,
                    c2,
                    start_time,
                    end_time,
                    stop_at_penetration,
                ),
            )
        /* } else if let (Some(p1), Some(s2)) = (shape1.as_shape::<HalfSpace>(), shape2.as_support_map()) {
        //        query::details::nonlinear_time_of_impact_halfspace_support_map(m1, vel1, p1, m2, vel2, s2)
                unimplemented!()
            } else if let (Some(s1), Some(p2)) = (shape1.as_support_map(), shape2.as_shape::<HalfSpace>()) {
        //        query::details::nonlinear_time_of_impact_support_map_halfspace(m1, vel1, s1, m2, vel2, p2)
                unimplemented!() */
        } else {
            Err(Unsupported)
        }
    }
}

impl<ManifoldData, ContactData> PersistentQueryDispatcher<ManifoldData, ContactData>
    for DefaultQueryDispatcher
where
    ManifoldData: Default + Clone,
    ContactData: Default + Copy,
{
    fn contact_manifolds(
        &self,
        pos12: &Isometry<Real>,
        shape1: &dyn Shape,
        shape2: &dyn Shape,
        prediction: Real,
        manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
        workspace: &mut Option<ContactManifoldsWorkspace>,
    ) -> Result<(), Unsupported> {
        use crate::query::contact_manifolds::*;

        let composite1 = shape1.as_composite_shape();
        let composite2 = shape2.as_composite_shape();

        if let (Some(composite1), Some(composite2)) = (composite1, composite2) {
            contact_manifolds_composite_shape_composite_shape(
                self, pos12, composite1, composite2, prediction, manifolds, workspace,
            );

            return Ok(());
        }

        match (shape1.shape_type(), shape2.shape_type()) {
            (ShapeType::TriMesh, _) | (_, ShapeType::TriMesh) => {
                contact_manifolds_trimesh_shape_shapes(
                    self, pos12, shape1, shape2, prediction, manifolds, workspace,
                );
            }
            (ShapeType::HeightField, _) => {
                if let Some(composite2) = composite2 {
                    contact_manifolds_heightfield_composite_shape(
                        self,
                        pos12,
                        &pos12.inverse(),
                        shape1.as_heightfield().unwrap(),
                        composite2,
                        prediction,
                        manifolds,
                        workspace,
                        false,
                    )
                } else {
                    contact_manifolds_heightfield_shape_shapes(
                        self, pos12, shape1, shape2, prediction, manifolds, workspace,
                    );
                }
            }
            (_, ShapeType::HeightField) => {
                if let Some(composite1) = composite1 {
                    contact_manifolds_heightfield_composite_shape(
                        self,
                        &pos12.inverse(),
                        pos12,
                        shape2.as_heightfield().unwrap(),
                        composite1,
                        prediction,
                        manifolds,
                        workspace,
                        true,
                    )
                } else {
                    contact_manifolds_heightfield_shape_shapes(
                        self, pos12, shape1, shape2, prediction, manifolds, workspace,
                    );
                }
            }
            _ => {
                if let Some(composite1) = composite1 {
                    contact_manifolds_composite_shape_shape(
                        self, pos12, composite1, shape2, prediction, manifolds, workspace, false,
                    );
                } else if let Some(composite2) = composite2 {
                    contact_manifolds_composite_shape_shape(
                        self,
                        &pos12.inverse(),
                        composite2,
                        shape1,
                        prediction,
                        manifolds,
                        workspace,
                        true,
                    );
                } else {
                    if manifolds.is_empty() {
                        manifolds.push(ContactManifold::new());
                    }

                    return self.contact_manifold_convex_convex(
                        pos12,
                        shape1,
                        shape2,
                        prediction,
                        &mut manifolds[0],
                    );
                }
            }
        }

        Ok(())
    }

    fn contact_manifold_convex_convex(
        &self,
        pos12: &Isometry<Real>,
        shape1: &dyn Shape,
        shape2: &dyn Shape,
        prediction: Real,
        manifold: &mut ContactManifold<ManifoldData, ContactData>,
    ) -> Result<(), Unsupported> {
        use crate::query::contact_manifolds::*;

        match (shape1.shape_type(), shape2.shape_type()) {
            (ShapeType::Ball, ShapeType::Ball) => {
                contact_manifold_ball_ball_shapes(pos12, shape1, shape2, prediction, manifold)
            }
            (ShapeType::Cuboid, ShapeType::Cuboid) =>
                contact_manifold_cuboid_cuboid_shapes(pos12, shape1, shape2, prediction, manifold)
            ,
            // (ShapeType::Polygon, ShapeType::Polygon) => (
            //     PrimitiveContactGenerator {
            //         generate_contacts: super::generate_contacts_polygon_polygon,
            //         ..PrimitiveContactGenerator::default()
            //     },
            //     None,
            // ),
            (ShapeType::Capsule, ShapeType::Capsule) => {
                contact_manifold_capsule_capsule_shapes(pos12, shape1, shape2, prediction, manifold)
            }
            (_, ShapeType::Ball) | (ShapeType::Ball, _) => {
                contact_manifold_convex_ball_shapes(pos12, shape1, shape2, prediction, manifold)
            }
            // (ShapeType::Capsule, ShapeType::Cuboid) | (ShapeType::Cuboid, ShapeType::Capsule) =>
            //     contact_manifold_cuboid_capsule_shapes(pos12, shape1, shape2, prediction, manifold),
            (ShapeType::Triangle, ShapeType::Cuboid) | (ShapeType::Cuboid, ShapeType::Triangle) => {
                contact_manifold_cuboid_triangle_shapes(pos12, shape1, shape2, prediction, manifold)
            }
            (ShapeType::HalfSpace, _) => {
                if let Some((pfm2, border_radius2)) = shape2.as_polygonal_feature_map() {
                    contact_manifold_halfspace_pfm(
                        pos12,
                        shape1.as_halfspace().unwrap(),
                        pfm2,
                        border_radius2,
                        prediction,
                        manifold,
                        false
                    )
                } else {
                    return Err(Unsupported)
                }
            }
            (_, ShapeType::HalfSpace) => {
                if let Some((pfm1, border_radius1)) = shape1.as_polygonal_feature_map() {
                    contact_manifold_halfspace_pfm(
                        &pos12.inverse(),
                        shape2.as_halfspace().unwrap(),
                        pfm1,
                        border_radius1,
                        prediction,
                        manifold,
                        true
                    )
                } else {
                    return Err(Unsupported)
                }
            }
            _ => {
                if let (Some(pfm1), Some(pfm2)) = (
                    shape1.as_polygonal_feature_map(),
                    shape2.as_polygonal_feature_map(),
                ) {
                    contact_manifold_pfm_pfm(
                        pos12, pfm1.0, pfm1.1, pfm2.0, pfm2.1, prediction, manifold,
                    )
                } else {
                    return Err(Unsupported);
                }
            }
        }

        Ok(())
    }
}
