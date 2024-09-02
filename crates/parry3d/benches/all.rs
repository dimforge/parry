#![feature(test)]
#![allow(unused_macros)]

use parry3d::math::{Point, Real};

extern crate nalgebra as na;
extern crate parry3d;
extern crate rand;
extern crate test;

mod bounding_volume;
mod common;
mod query;
mod support_map;

#[cfg(feature = "dim2")]
type ConvexHull = Vec<Point<Real>>;
#[cfg(feature = "dim3")]
type ConvexHull = (Vec<Point<Real>>, Vec<[u32; 3]>);
