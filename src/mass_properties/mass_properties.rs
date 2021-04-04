use crate::math::{AngVector, AngularInertia, Isometry, Point, Real, Rotation, Vector};
use crate::utils;
use na::ComplexField;
use num::Zero;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Sub, SubAssign};
#[cfg(feature = "dim3")]
use {na::Matrix3, std::ops::MulAssign};

const EPSILON: Real = f32::EPSILON as Real;

#[derive(Copy, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
/// The local mass properties of a rigid-body.
pub struct MassProperties {
    /// The center of mass of a rigid-body expressed in its local-space.
    pub local_com: Point<Real>,
    /// The inverse of the mass of a rigid-body.
    ///
    /// If this is zero, the rigid-body is assumed to have infinite mass.
    pub inv_mass: Real,
    /// The inverse of the principal angular inertia of the rigid-body.
    ///
    /// Components set to zero are assumed to be infinite along the corresponding principal axis.
    pub inv_principal_inertia_sqrt: AngVector<Real>,
    #[cfg(feature = "dim3")]
    /// The principal vectors of the local angular inertia tensor of the rigid-body.
    pub principal_inertia_local_frame: Rotation<Real>,
}

impl MassProperties {
    /// Initializes the mass properties with the given center-of-mass, mass, and angular inertia.
    ///
    /// The center-of-mass is specified in the local-space of the rigid-body.
    #[cfg(feature = "dim2")]
    pub fn new(local_com: Point<Real>, mass: Real, principal_inertia: Real) -> Self {
        let inv_mass = utils::inv(mass);
        let inv_principal_inertia_sqrt = utils::inv(ComplexField::sqrt(principal_inertia));
        Self {
            local_com,
            inv_mass,
            inv_principal_inertia_sqrt,
        }
    }

    /// Initializes the mass properties from the given center-of-mass, mass, and principal angular inertia.
    ///
    /// The center-of-mass is specified in the local-space of the rigid-body.
    /// The principal angular inertia are the angular inertia along the coordinate axes in the local-space
    /// of the rigid-body.
    #[cfg(feature = "dim3")]
    pub fn new(local_com: Point<Real>, mass: Real, principal_inertia: AngVector<Real>) -> Self {
        Self::with_principal_inertia_frame(local_com, mass, principal_inertia, Rotation::identity())
    }

    /// Initializes the mass properties from the given center-of-mass, mass, and principal angular inertia.
    ///
    /// The center-of-mass is specified in the local-space of the rigid-body.
    /// The principal angular inertia are the angular inertia along the coordinate axes defined by
    /// the `principal_inertia_local_frame` expressed in the local-space of the rigid-body.
    #[cfg(feature = "dim3")]
    pub fn with_principal_inertia_frame(
        local_com: Point<Real>,
        mass: Real,
        principal_inertia: AngVector<Real>,
        principal_inertia_local_frame: Rotation<Real>,
    ) -> Self {
        let inv_mass = utils::inv(mass);
        let inv_principal_inertia_sqrt =
            principal_inertia.map(|e| utils::inv(ComplexField::sqrt(e)));
        Self {
            local_com,
            inv_mass,
            inv_principal_inertia_sqrt,
            principal_inertia_local_frame,
        }
    }

    /// Initialize a new `MassProperties` from a given center-of-mass, mass, and angular inertia matrix.
    ///
    /// The angular inertia matrix will be diagonalized in order to extract the principal inertia
    /// values and principal inertia frame.
    #[cfg(feature = "dim3")]
    pub fn with_inertia_matrix(local_com: Point<Real>, mass: Real, inertia: Matrix3<Real>) -> Self {
        let eigen = inertia.symmetric_eigen();
        let principal_inertia_local_frame =
            Rotation::from_matrix_eps(&eigen.eigenvectors, 1.0e-6, 10, na::one());
        // Drop negative eigenvalues.
        let principal_inertia = eigen.eigenvalues.map(|e| if e < EPSILON { 0.0 } else { e });

        Self::with_principal_inertia_frame(
            local_com,
            mass,
            principal_inertia,
            principal_inertia_local_frame,
        )
    }

    /// The world-space center of mass of the rigid-body.
    pub fn world_com(&self, pos: &Isometry<Real>) -> Point<Real> {
        pos * self.local_com
    }

    #[cfg(feature = "dim2")]
    /// The world-space inverse angular inertia tensor of the rigid-body.
    pub fn world_inv_inertia_sqrt(&self, _rot: &Rotation<Real>) -> AngularInertia<Real> {
        self.inv_principal_inertia_sqrt
    }

    #[cfg(feature = "dim3")]
    /// The world-space inverse angular inertia tensor of the rigid-body.
    pub fn world_inv_inertia_sqrt(&self, rot: &Rotation<Real>) -> AngularInertia<Real> {
        if !self.inv_principal_inertia_sqrt.is_zero() {
            let mut lhs = (rot * self.principal_inertia_local_frame)
                .to_rotation_matrix()
                .into_inner();
            let rhs = lhs.transpose();
            lhs.column_mut(0)
                .mul_assign(self.inv_principal_inertia_sqrt.x);
            lhs.column_mut(1)
                .mul_assign(self.inv_principal_inertia_sqrt.y);
            lhs.column_mut(2)
                .mul_assign(self.inv_principal_inertia_sqrt.z);
            let inertia = lhs * rhs;
            AngularInertia::from_sdp_matrix(inertia)
        } else {
            AngularInertia::zero()
        }
    }

    #[cfg(feature = "dim3")]
    /// Reconstructs the inverse angular inertia tensor of the rigid body from its principal inertia values and axes.
    pub fn reconstruct_inverse_inertia_matrix(&self) -> Matrix3<Real> {
        let inv_principal_inertia = self.inv_principal_inertia_sqrt.map(|e| e * e);
        self.principal_inertia_local_frame.to_rotation_matrix()
            * Matrix3::from_diagonal(&inv_principal_inertia)
            * self
                .principal_inertia_local_frame
                .inverse()
                .to_rotation_matrix()
    }

    #[cfg(feature = "dim3")]
    /// Reconstructs the angular inertia tensor of the rigid body from its principal inertia values and axes.
    pub fn reconstruct_inertia_matrix(&self) -> Matrix3<Real> {
        let principal_inertia = self.inv_principal_inertia_sqrt.map(|e| utils::inv(e * e));
        self.principal_inertia_local_frame.to_rotation_matrix()
            * Matrix3::from_diagonal(&principal_inertia)
            * self
                .principal_inertia_local_frame
                .inverse()
                .to_rotation_matrix()
    }

    #[cfg(feature = "dim2")]
    pub(crate) fn construct_shifted_inertia_matrix(&self, shift: Vector<Real>) -> Real {
        let i = utils::inv(self.inv_principal_inertia_sqrt * self.inv_principal_inertia_sqrt);

        if self.inv_mass != 0.0 {
            let mass = 1.0 / self.inv_mass;
            i + shift.norm_squared() * mass
        } else {
            i
        }
    }

    #[cfg(feature = "dim3")]
    pub(crate) fn construct_shifted_inertia_matrix(&self, shift: Vector<Real>) -> Matrix3<Real> {
        let matrix = self.reconstruct_inertia_matrix();

        if self.inv_mass != 0.0 {
            let mass = 1.0 / self.inv_mass;
            let diag = shift.norm_squared();
            let diagm = Matrix3::from_diagonal_element(diag);
            matrix + (diagm + shift * shift.transpose()) * mass
        } else {
            matrix
        }
    }

    /// Transform each element of the mass properties.
    pub fn transform_by(&self, m: &Isometry<Real>) -> Self {
        // NOTE: we don't apply the parallel axis theorem here
        // because the center of mass is also transformed.
        Self {
            local_com: m * self.local_com,
            inv_mass: self.inv_mass,
            inv_principal_inertia_sqrt: self.inv_principal_inertia_sqrt,
            #[cfg(feature = "dim3")]
            principal_inertia_local_frame: m.rotation * self.principal_inertia_local_frame,
        }
    }

    /// Changes the mass on these mass-properties.
    ///
    /// The `adjust_angular_inertia` argument should always be `true`, unless
    /// there are some specific reasons not to do so. Setting this to `true`
    /// will automatically adjust the angular inertia of `self` to account
    /// for the mass change (i.e. it will multiply the angular inertia by
    /// `new_mass / prev_mass`). Setting it to `false` will not change the
    /// current angular inertia.
    pub fn set_mass(&mut self, new_mass: Real, adjust_angular_inertia: bool) {
        let new_inv_mass = utils::inv(new_mass);

        if adjust_angular_inertia {
            let curr_mass = utils::inv(self.inv_mass);
            self.inv_principal_inertia_sqrt *= new_inv_mass.sqrt() * curr_mass.sqrt();
        }

        self.inv_mass = new_inv_mass;
    }
}

impl Zero for MassProperties {
    fn zero() -> Self {
        Self {
            inv_mass: 0.0,
            inv_principal_inertia_sqrt: na::zero(),
            #[cfg(feature = "dim3")]
            principal_inertia_local_frame: Rotation::identity(),
            local_com: Point::origin(),
        }
    }

    fn is_zero(&self) -> bool {
        *self == Self::zero()
    }
}

impl Sub<MassProperties> for MassProperties {
    type Output = Self;

    #[cfg(feature = "dim2")]
    fn sub(self, other: MassProperties) -> Self {
        if self.is_zero() || other.is_zero() {
            return self;
        }

        let m1 = utils::inv(self.inv_mass);
        let m2 = utils::inv(other.inv_mass);

        let mut new_mass = m1 - m2;

        if new_mass < EPSILON {
            // Account for small numerical errors.
            new_mass = 0.0;
        }

        let inv_mass = utils::inv(new_mass);

        let local_com = (self.local_com * m1 - other.local_com.coords * m2) * inv_mass;
        let i1 = self.construct_shifted_inertia_matrix(local_com - self.local_com);
        let i2 = other.construct_shifted_inertia_matrix(local_com - other.local_com);
        let mut inertia = i1 - i2;

        if inertia < EPSILON {
            // Account for small numerical errors.
            inertia = 0.0;
        }

        // NOTE: we drop the negative eigenvalues that may result from subtraction rounding errors.
        let inv_principal_inertia_sqrt = utils::inv(ComplexField::sqrt(inertia));

        Self {
            local_com,
            inv_mass,
            inv_principal_inertia_sqrt,
        }
    }

    #[cfg(feature = "dim3")]
    fn sub(self, other: MassProperties) -> Self {
        if self.is_zero() || other.is_zero() {
            return self;
        }

        let m1 = utils::inv(self.inv_mass);
        let m2 = utils::inv(other.inv_mass);
        let mut new_mass = m1 - m2;

        if new_mass < EPSILON {
            new_mass = 0.0;
        }

        let inv_mass = utils::inv(new_mass);
        let local_com = (self.local_com * m1 - other.local_com.coords * m2) * inv_mass;
        let i1 = self.construct_shifted_inertia_matrix(local_com - self.local_com);
        let i2 = other.construct_shifted_inertia_matrix(local_com - other.local_com);
        let inertia = i1 - i2;
        Self::with_inertia_matrix(local_com, new_mass, inertia)
    }
}

impl SubAssign<MassProperties> for MassProperties {
    fn sub_assign(&mut self, rhs: MassProperties) {
        *self = *self - rhs
    }
}

impl Add<MassProperties> for MassProperties {
    type Output = Self;

    #[cfg(feature = "dim2")]
    fn add(self, other: MassProperties) -> Self {
        if self.is_zero() {
            return other;
        } else if other.is_zero() {
            return self;
        }

        let m1 = utils::inv(self.inv_mass);
        let m2 = utils::inv(other.inv_mass);
        let inv_mass = utils::inv(m1 + m2);
        let local_com = (self.local_com * m1 + other.local_com.coords * m2) * inv_mass;
        let i1 = self.construct_shifted_inertia_matrix(local_com - self.local_com);
        let i2 = other.construct_shifted_inertia_matrix(local_com - other.local_com);
        let inertia = i1 + i2;
        let inv_principal_inertia_sqrt = utils::inv(ComplexField::sqrt(inertia));

        Self {
            local_com,
            inv_mass,
            inv_principal_inertia_sqrt,
        }
    }

    #[cfg(feature = "dim3")]
    fn add(self, other: MassProperties) -> Self {
        if self.is_zero() {
            return other;
        } else if other.is_zero() {
            return self;
        }

        let m1 = utils::inv(self.inv_mass);
        let m2 = utils::inv(other.inv_mass);
        let inv_mass = utils::inv(m1 + m2);
        let local_com = (self.local_com * m1 + other.local_com.coords * m2) * inv_mass;
        let i1 = self.construct_shifted_inertia_matrix(local_com - self.local_com);
        let i2 = other.construct_shifted_inertia_matrix(local_com - other.local_com);
        let inertia = i1 + i2;

        Self::with_inertia_matrix(local_com, m1 + m2, inertia)
    }
}

impl AddAssign<MassProperties> for MassProperties {
    fn add_assign(&mut self, rhs: MassProperties) {
        *self = *self + rhs
    }
}

impl Sum<MassProperties> for MassProperties {
    #[cfg(feature = "dim2")]
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let mut total_mass = 0.0;
        let mut total_com = Point::origin();
        let mut total_inertia = 0.0;
        // TODO: avoid this allocation.
        // This is needed because we iterate twice.
        let mut all_props = Vec::new();

        for props in iter {
            let mass = utils::inv(props.inv_mass);
            total_mass += mass;
            total_com += props.local_com.coords * mass;
            all_props.push(props);
        }

        if total_mass > 0.0 {
            total_com /= total_mass;
        }

        for props in all_props {
            total_inertia += props.construct_shifted_inertia_matrix(total_com - props.local_com);
        }

        Self {
            local_com: total_com,
            inv_mass: utils::inv(total_mass),
            inv_principal_inertia_sqrt: utils::inv(ComplexField::sqrt(total_inertia)),
        }
    }

    #[cfg(feature = "dim3")]
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        let mut total_mass = 0.0;
        let mut total_com = Point::origin();
        let mut total_inertia = Matrix3::zeros();
        // TODO: avoid this allocation.
        // This is needed because we iterate twice.
        let mut all_props = Vec::new();

        for props in iter {
            let mass = utils::inv(props.inv_mass);
            total_mass += mass;
            total_com += props.local_com.coords * mass;
            all_props.push(props);
        }

        if total_mass > 0.0 {
            total_com /= total_mass;
        }

        for props in all_props {
            total_inertia += props.construct_shifted_inertia_matrix(total_com - props.local_com);
        }

        Self::with_inertia_matrix(total_com, total_mass, total_inertia)
    }
}

impl approx::AbsDiffEq for MassProperties {
    type Epsilon = Real;
    fn default_epsilon() -> Self::Epsilon {
        Real::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        #[cfg(feature = "dim2")]
        let inertia_is_ok = self
            .inv_principal_inertia_sqrt
            .abs_diff_eq(&other.inv_principal_inertia_sqrt, epsilon);

        #[cfg(feature = "dim3")]
        let inertia_is_ok = self
            .reconstruct_inverse_inertia_matrix()
            .abs_diff_eq(&other.reconstruct_inverse_inertia_matrix(), epsilon);

        inertia_is_ok
            && self.local_com.abs_diff_eq(&other.local_com, epsilon)
            && self.inv_mass.abs_diff_eq(&other.inv_mass, epsilon)
            && self
                .inv_principal_inertia_sqrt
                .abs_diff_eq(&other.inv_principal_inertia_sqrt, epsilon)
    }
}

impl approx::RelativeEq for MassProperties {
    fn default_max_relative() -> Self::Epsilon {
        Real::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        #[cfg(feature = "dim2")]
        let inertia_is_ok = self.inv_principal_inertia_sqrt.relative_eq(
            &other.inv_principal_inertia_sqrt,
            epsilon,
            max_relative,
        );

        #[cfg(feature = "dim3")]
        let inertia_is_ok = self.reconstruct_inverse_inertia_matrix().relative_eq(
            &other.reconstruct_inverse_inertia_matrix(),
            epsilon,
            max_relative,
        );

        inertia_is_ok
            && self
                .local_com
                .relative_eq(&other.local_com, epsilon, max_relative)
            && self
                .inv_mass
                .relative_eq(&other.inv_mass, epsilon, max_relative)
    }
}

#[cfg(test)]
mod test {
    use super::MassProperties;
    use crate::math::Point;
    #[cfg(feature = "dim3")]
    use crate::math::{Rotation, Vector};
    use crate::shape::{Ball, Capsule, Shape};
    use approx::assert_relative_eq;
    use num::Zero;

    #[test]
    fn mass_properties_add_partial_zero() {
        let m1 = MassProperties {
            local_com: Point::origin(),
            inv_mass: 2.0,
            inv_principal_inertia_sqrt: na::zero(),
            #[cfg(feature = "dim3")]
            principal_inertia_local_frame: Rotation::identity(),
        };
        let m2 = MassProperties {
            local_com: Point::origin(),
            inv_mass: 0.0,
            #[cfg(feature = "dim2")]
            inv_principal_inertia_sqrt: 1.0,
            #[cfg(feature = "dim3")]
            inv_principal_inertia_sqrt: Vector::new(1.0, 2.0, 3.0),
            #[cfg(feature = "dim3")]
            principal_inertia_local_frame: Rotation::identity(),
        };
        let result = MassProperties {
            local_com: Point::origin(),
            inv_mass: 2.0,
            #[cfg(feature = "dim2")]
            inv_principal_inertia_sqrt: 1.0,
            #[cfg(feature = "dim3")]
            inv_principal_inertia_sqrt: Vector::new(1.0, 2.0, 3.0),
            #[cfg(feature = "dim3")]
            principal_inertia_local_frame: Rotation::identity(),
        };

        assert_eq!(m1 + m2, result);
        assert_eq!(m2 + m1, result);
    }

    #[test]
    fn mass_properties_add_sub() {
        // Check that addition and subtraction of mass properties behave as expected.
        let c1 = Capsule::new_x(1.0, 2.0);
        let c2 = Capsule::new_y(3.0, 4.0);
        let c3 = Ball::new(5.0);

        let m1 = c1.mass_properties(1.0);
        let m2 = c2.mass_properties(1.0);
        let m3 = c3.mass_properties(1.0);
        let m1m2m3 = m1 + m2 + m3;

        assert_relative_eq!(m1 + m2, m2 + m1, epsilon = 1.0e-6);
        assert_relative_eq!(m1m2m3 - m1, m2 + m3, epsilon = 1.0e-6);
        assert_relative_eq!(m1m2m3 - m2, m1 + m3, epsilon = 1.0e-6);
        assert_relative_eq!(m1m2m3 - m3, m1 + m2, epsilon = 1.0e-6);
        assert_relative_eq!(m1m2m3 - (m1 + m2), m3, epsilon = 1.0e-6);
        assert_relative_eq!(m1m2m3 - (m1 + m3), m2, epsilon = 1.0e-6);
        assert_relative_eq!(m1m2m3 - (m2 + m3), m1, epsilon = 1.0e-6);
        assert_relative_eq!(m1m2m3 - m1 - m2, m3, epsilon = 1.0e-6);
        assert_relative_eq!(m1m2m3 - m1 - m3, m2, epsilon = 1.0e-6);
        assert_relative_eq!(m1m2m3 - m2 - m3, m1, epsilon = 1.0e-6);
        assert_relative_eq!(
            m1m2m3 - m1 - m2 - m3,
            MassProperties::zero(),
            epsilon = 1.0e-6
        );
    }

    #[test]
    fn mass_properties_sum_no_nan() {
        let mp: MassProperties = [MassProperties::zero()].iter().map(|v| *v).sum();
        assert!(!mp.local_com.x.is_nan() && !mp.local_com.y.is_nan());
        #[cfg(feature = "dim3")]
        assert!(!mp.local_com.z.is_nan());
    }
}
