//! # csx
//!
//! Pure Rust continuous-structure geometry and material modeling, compatible
//! with CSXCAD and openEMS. See `PLAN.md` for the roadmap.

pub mod coords;
pub mod error;
pub mod expr;
pub mod grid;
pub mod io;
pub mod math;
pub mod param;
pub mod primitives;
pub mod properties;
pub mod structure;
pub mod transform;

pub use coords::CoordSystem;
pub use error::{Error, Result};
pub use grid::RectilinearGrid;
pub use math::{BoundingBox, Mat4, Vec3};
pub use param::{Parameter, ParameterSet, Scalar, ScalarPoint};
pub use primitives::{Primitive, Shape};
pub use properties::{Material, Property, PropertyKind, PropertyType, TypeMask};
pub use structure::{ContinuousStructure, PropertyId};
pub use transform::{Transform, TransformOp};
