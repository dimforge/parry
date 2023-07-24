use crate::math::{Real, real};

pub fn inv(val: Real) -> Real {
    if val == real!(0.0) {
        real!(0.0)
    } else {
        real!(1.0) / val
    }
}
