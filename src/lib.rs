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
//! ## Module overview
//!
//! * [`structure`], [`properties`], [`primitives`]: the model.
//! * [`param`], [`expr`]: parametric values and expressions.
//! * [`grid`], [`mesh`]: mesh lines and automatic mesh generation.
//! * [`discretize`]: sampling the model onto the FDTD grid.
//! * [`tessellate`], [`io`]: triangle meshes, file import/export and (with
//!   the `xml` feature) CSXCAD / openEMS XML.
//! * [`fdtd`]: openEMS simulation settings.
//! * [`prelude`]: glob-importable common types.
//!
//! ## Conventions
//!
//! * Axis indices are `usize`: `0` = x (ρ), `1` = y (α), `2` = z.
//! * Lengths are in drawing units (see
//!   [`RectilinearGrid::delta_unit`]); material parameters are SI.
//! * Angles are radians.
//! * Point queries always take the [`CoordSystem`] of the point.
//! * Parametric values are re-evaluated by [`ContinuousStructure::update`].
//! * NaN marks unset values ([`Scalar::is_unset`]) and is false in
//!   expressions.
//!
//! Details: [`coords`](coords#crate-wide-conventions) and [`expr`].
//!
//! ## Stability
//!
//! Enums and option/result structs that are expected to grow (errors, shape
//! and property kinds, type ids, `*Options`, ...) are `#[non_exhaustive]`:
//! match them with a wildcard arm and build options from `Default` plus
//! their `with_*` builders.
//!
//! This is a clean-room implementation written from public documentation and
//! the observable file format only.

pub mod coords;
pub mod discretize;
pub mod error;
pub mod expr;
pub mod fdtd;
pub mod grid;
pub mod io;
pub mod math;
pub mod mesh;
pub mod param;
pub mod primitives;
pub mod properties;
pub mod structure;
pub mod tessellate;
pub mod transform;

pub use coords::CoordSystem;
pub use error::{Error, Result};
pub use grid::RectilinearGrid;
#[cfg(feature = "xml")]
pub use io::xml::OpenEmsDocument;
pub use math::{BoundingBox, Mat4, Vec3};
pub use param::{Parameter, ParameterSet, Scalar, ScalarPoint};
pub use primitives::{Primitive, Shape};
pub use properties::{Axis, Material, Property, PropertyKind, PropertyType, TypeMask};
pub use structure::{BackgroundMaterial, ContinuousStructure, Hit, Medium, PropertyId};
pub use transform::{Transform, TransformOp};

/// The most commonly used types, for glob import:
///
/// ```
/// use csx::prelude::*;
///
/// let mut s = ContinuousStructure::new();
/// let id = s.add_property(Property::metal("pec"));
/// s.add_primitive(id, BoxShape::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]))?;
/// assert!(s.property_at(Vec3::new(0.5, 0.5, 0.5), CoordSystem::Cartesian, TypeMask::ANY).is_some());
/// # Ok::<(), csx::Error>(())
/// ```
pub mod prelude {
    pub use crate::fdtd::{Boundary, BoundaryConditions, ExcitationSignal, FdtdSettings};
    #[cfg(feature = "xml")]
    pub use crate::io::xml::OpenEmsDocument;
    pub use crate::mesh::{MeshBuilder, MeshOptions};
    pub use crate::primitives::{
        BoxShape, Curve, Cylinder, CylindricalShell, LinPoly, MultiBox, PointShape, Polygon,
        Polyhedron, PolyhedronReader, RotPoly, Sphere, SphericalShell, UserDefined, Wire,
    };
    pub use crate::properties::{
        ConductingSheet, DebyeMaterial, DumpBox, DumpType, Excitation, ExcitationType,
        LorentzMaterial, LumpedElement, Material, ProbeBox, ProbeType, Rgba,
    };
    pub use crate::{
        Axis, BoundingBox, ContinuousStructure, CoordSystem, Error, ParameterSet, Primitive,
        Property, PropertyKind, PropertyType, RectilinearGrid, Result, Scalar, ScalarPoint, Shape,
        Transform, TypeMask, Vec3,
    };
}
