//! # csx
//!
//! Pure Rust continuous-structure geometry and material modeling, compatible
//! with CSXCAD and openEMS.
//!
//! A [`ContinuousStructure`] holds named [`Property`]s (materials, metals,
//! excitations, probes, ...). Each property owns [`Primitive`]s describing
//! where it applies; overlapping primitives are resolved by priority. The
//! structure also carries a [`RectilinearGrid`] (the FDTD mesh), a
//! [`ParameterSet`] for parametric models and the background material.
//!
//! ```
//! use csx::primitives::{BoxShape, Cylinder};
//! use csx::properties::Material;
//! use csx::{ContinuousStructure, CoordSystem, Primitive, Property, TypeMask, Vec3};
//!
//! let mut model = ContinuousStructure::new();
//! model.parameters.set("h", 1.6);
//!
//! // FR4 substrate, its thickness given by the parameter `h`.
//! let fr4 = model.add_property(Property::material("FR4", Material::new(4.3, 0.0)));
//! model.add_primitive(fr4, Primitive::new(BoxShape {
//!     start: [-10.0, -10.0, 0.0].into(),
//!     stop: csx::ScalarPoint(["10".into(), "10".into(), "h".into()]),
//! }))?;
//!
//! // A metal via with higher priority.
//! let pec = model.add_property(Property::metal("via"));
//! model.add_primitive(pec, Primitive::new(Cylinder::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.6], 0.3))
//!     .with_priority(10))?;
//!
//! model.update()?; // evaluate parametric values
//!
//! let hit = model.property_at(Vec3::new(5.0, 0.0, 0.8), CoordSystem::Cartesian, TypeMask::ANY);
//! assert_eq!(hit.unwrap().property.name, "FR4");
//! let hit = model.property_at(Vec3::new(0.1, 0.0, 0.8), CoordSystem::Cartesian, TypeMask::ANY);
//! assert_eq!(hit.unwrap().property.name, "via");
//! # Ok::<(), csx::Error>(())
//! ```
//!
//! This is a clean-room implementation written from public documentation and
//! the observable file format only.

pub mod coords;
pub mod error;
pub mod expr;
pub mod grid;
pub mod io;
pub mod math;
pub mod mesh;
pub mod param;
pub mod primitives;
pub mod properties;
pub mod structure;
pub mod transform;

pub use coords::CoordSystem;
pub use error::{Error, Result};
pub use grid::RectilinearGrid;
#[cfg(feature = "xml")]
pub use io::xml::OpenEmsDocument;
pub use math::{BoundingBox, Mat4, Vec3};
pub use param::{Parameter, ParameterSet, Scalar, ScalarPoint};
pub use primitives::{Primitive, Shape};
pub use properties::{Material, Property, PropertyKind, PropertyType, TypeMask};
pub use structure::{ContinuousStructure, PropertyId};
pub use transform::{Transform, TransformOp};
