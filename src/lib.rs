//! # csx
//!
//! Pure Rust continuous-structure geometry and material modeling, compatible
//! with CSXCAD and openEMS. See `PLAN.md` for the roadmap.

pub mod coords;
pub mod math;

pub use coords::CoordSystem;
pub use math::{BoundingBox, Mat4, Vec3};
