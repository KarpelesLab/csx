//! The continuous structure: the complete geometry and material model.
//!
//! A [`ContinuousStructure`] is the root object of a model. It owns
//!
//! | part | type | purpose |
//! |---|---|---|
//! | properties | [`Property`] (via [`PropertyId`] handles) | *what* a region is: material, metal, port, probe, ... |
//! | primitives | [`Primitive`], owned by each property | *where* the property applies |
//! | mesh | [`RectilinearGrid`] | FDTD mesh lines |
//! | parameters | [`ParameterSet`] | named values used by expressions |
//! | background | [`BackgroundMaterial`] | the medium outside every material |
//!
//! # Queries
//!
//! * [`ContinuousStructure::property_at`] — the winning property at a point
//!   (highest primitive priority; ties go to the property added first).
//! * [`ContinuousStructure::properties_at`] — every property at a point.
//! * [`ContinuousStructure::medium_at`] — ε, μ, κ, σ at a point.
//! * [`ContinuousStructure::spatial_index`] — an accelerated
//!   [`SpatialIndex`] for many point queries.
//! * [`ContinuousStructure::bounding_box`],
//!   [`ContinuousStructure::primitives_in`] — spatial extent.
//! * [`ContinuousStructure::validate`] — consistency warnings
//!   ([`ValidationIssue`]).
//!
//! # Example
//!
//! ```
//! use csx::primitives::BoxShape;
//! use csx::properties::Material;
//! use csx::{ContinuousStructure, CoordSystem, Primitive, Property, TypeMask, Vec3};
//!
//! let mut s = ContinuousStructure::new();
//! let sub = s.add_property(Property::material("substrate", Material::new(3.38, 0.0)));
//! s.add_primitive(sub, BoxShape::new([0.0, 0.0, 0.0], [10.0, 10.0, 0.8]))?;
//! let gnd = s.add_property(Property::metal("ground"));
//! s.add_primitive(gnd, Primitive::new(BoxShape::new([0.0, 0.0, 0.0], [10.0, 10.0, 0.0]))
//!     .with_priority(10))?;
//!
//! let p = Vec3::new(5.0, 5.0, 0.4);
//! let hit = s.property_at(p, CoordSystem::Cartesian, TypeMask::ANY).unwrap();
//! assert_eq!(hit.property.name, "substrate");
//! assert_eq!(s.medium_at(p, CoordSystem::Cartesian).unwrap().epsilon, [3.38; 3]);
//! // metal is reported as `None`
//! assert!(s.medium_at(Vec3::new(5.0, 5.0, 0.0), CoordSystem::Cartesian).is_none());
//! # Ok::<(), csx::Error>(())
//! ```

use std::path::PathBuf;

use crate::coords::CoordSystem;
use crate::error::{Error, Result};
use crate::grid::RectilinearGrid;
use crate::math::{BoundingBox, Vec3};
use crate::param::{ParameterSet, Scalar};
use crate::primitives::bvh::Bvh;
use crate::primitives::{Primitive, Shape};
use crate::properties::{MaterialQuantity, Property, PropertyKind, TypeMask};

/// Material filling all space not covered by any material primitive
/// (`<BackgroundMaterial>`; default: vacuum).
///
/// ```
/// use csx::ContinuousStructure;
/// use csx::properties::MaterialQuantity;
///
/// let mut s = ContinuousStructure::new();
/// assert_eq!(s.background.value(MaterialQuantity::Epsilon), 1.0);
/// s.background.epsilon.set(1.0006); // air
/// assert_eq!(s.background.value(MaterialQuantity::Epsilon), 1.0006);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct BackgroundMaterial {
    /// Relative permittivity.
    pub epsilon: Scalar,
    /// Relative permeability.
    pub mue: Scalar,
    /// Electric conductivity.
    pub kappa: Scalar,
    /// Magnetic conductivity.
    pub sigma: Scalar,
    /// Unrecognised attributes preserved from files.
    pub extra_attributes: Vec<(String, String)>,
}

impl Default for BackgroundMaterial {
    fn default() -> Self {
        BackgroundMaterial {
            epsilon: Scalar::new(1.0),
            mue: Scalar::new(1.0),
            kappa: Scalar::new(0.0),
            sigma: Scalar::new(0.0),
            extra_attributes: Vec::new(),
        }
    }
}

impl BackgroundMaterial {
    fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.epsilon.update(p)?;
        self.mue.update(p)?;
        self.kappa.update(p)?;
        self.sigma.update(p)
    }

    /// Cached value of a quantity (see the update contract on
    /// [`ContinuousStructure`]).
    pub fn value(&self, q: MaterialQuantity) -> f64 {
        match q {
            MaterialQuantity::Epsilon => self.epsilon.value(),
            MaterialQuantity::Mue => self.mue.value(),
            MaterialQuantity::Kappa => self.kappa.value(),
            MaterialQuantity::Sigma => self.sigma.value(),
        }
    }
}

/// Opaque handle of a property inside a [`ContinuousStructure`].
///
/// Handles are returned by [`ContinuousStructure::add_property`] and
/// [`ContinuousStructure::property_id`]; they stay valid until the property
/// is removed and are never reused within one structure. Handles order by
/// insertion.
///
/// ```
/// use csx::{ContinuousStructure, Property};
///
/// let mut s = ContinuousStructure::new();
/// let a = s.add_property(Property::metal("a"));
/// let b = s.add_property(Property::metal("b"));
/// assert!(a < b);
/// assert_eq!(s.property_id("b"), Some(b));
/// assert_eq!(s.property(a).unwrap().id(), a.get());
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct PropertyId(u32);

impl PropertyId {
    /// The numeric id, equal to [`Property::id`] of the referenced
    /// property.
    pub const fn get(self) -> u32 {
        self.0
    }
}

impl core::fmt::Display for PropertyId {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "#{}", self.0)
    }
}

/// A consistency problem reported by [`ContinuousStructure::validate`].
///
/// ```
/// use csx::structure::ValidationIssue;
/// use csx::{ContinuousStructure, Property};
///
/// let mut s = ContinuousStructure::new();
/// s.add_property(Property::metal("m"));
/// let issues = s.validate();
/// assert!(matches!(&issues[0], ValidationIssue::EmptyProperty { name } if name == "m"));
/// assert_eq!(issues[0].to_string(), "property \"m\" has no primitives");
/// ```
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ValidationIssue {
    /// Several properties share a name (reported for each later one);
    /// name-based lookups return the first.
    DuplicateName {
        /// The repeated name.
        name: String,
    },
    /// A property without primitives has no effect.
    EmptyProperty {
        /// Property name.
        name: String,
    },
    /// A primitive whose bounding box is empty (e.g. a polygon without
    /// vertices, a mesh file that is not loaded).
    EmptyPrimitive {
        /// Property name.
        property: String,
        /// Index of the primitive within the property.
        index: usize,
        /// Element name of the shape (`"Box"`, `"Polygon"`, ...).
        shape: &'static str,
    },
}

impl core::fmt::Display for ValidationIssue {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ValidationIssue::DuplicateName { name } => {
                write!(f, "duplicate property name {name:?}")
            }
            ValidationIssue::EmptyProperty { name } => {
                write!(f, "property {name:?} has no primitives")
            }
            ValidationIssue::EmptyPrimitive {
                property,
                index,
                shape,
            } => write!(
                f,
                "property {property:?}: primitive #{index} ({shape}) is empty"
            ),
        }
    }
}

/// The result of a point query: a property and the primitive that contains
/// the point.
#[non_exhaustive]
#[derive(Debug, Clone, Copy)]
pub struct Hit<'a> {
    /// The property.
    pub property: &'a Property,
    /// The primitive that contains the point.
    pub primitive: &'a Primitive,
    /// Index of the primitive within the property.
    pub primitive_index: usize,
}

/// Electromagnetic medium at a point, as returned by
/// [`ContinuousStructure::medium_at`]: per-axis values with weighting
/// functions applied.
///
/// | field | quantity | unit |
/// |---|---|---|
/// | `epsilon` | relative permittivity εᵣ | – |
/// | `mue` | relative permeability μᵣ | – |
/// | `kappa` | electric conductivity κ | S/m |
/// | `sigma` | magnetic conductivity σ | Ω/m |
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Medium {
    /// Relative permittivity per axis.
    pub epsilon: [f64; 3],
    /// Relative permeability per axis.
    pub mue: [f64; 3],
    /// Electric conductivity per axis.
    pub kappa: [f64; 3],
    /// Magnetic conductivity per axis.
    pub sigma: [f64; 3],
}

/// The complete model: properties with their primitives, the mesh, the
/// parameters and the background material.
///
/// # The update contract
///
/// Values that may be expressions ([`Scalar`], polygon vertices, weighting
/// functions, transformations, user-defined functions) cache their evaluated
/// form. Plain numbers are usable immediately; anything referring to
/// parameters is evaluated by [`ContinuousStructure::update`]. Call it after
/// changing [`ContinuousStructure::parameters`] or adding parametric
/// content, before querying. Caches that could silently go stale
/// (compiled user functions, loaded mesh files) are only reachable through
/// setters that invalidate them.
///
/// ```
/// use csx::primitives::BoxShape;
/// use csx::{ContinuousStructure, CoordSystem, Property, ScalarPoint, TypeMask, Vec3};
///
/// let mut s = ContinuousStructure::new();
/// s.parameters.set("w", 2.0);
/// let id = s.add_property(Property::metal("strip"));
/// s.add_primitive(id, BoxShape::new([0.0, 0.0, 0.0], ScalarPoint::from(["w", "1", "0"])))?;
/// s.update()?; // evaluate `w`
///
/// let probe = Vec3::new(1.5, 0.5, 0.0);
/// assert!(s.property_at(probe, CoordSystem::Cartesian, TypeMask::ANY).is_some());
///
/// s.parameters.set("w", 1.0);
/// s.update()?; // geometry follows the parameter
/// assert!(s.property_at(probe, CoordSystem::Cartesian, TypeMask::ANY).is_none());
/// # Ok::<(), csx::Error>(())
/// ```
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ContinuousStructure {
    properties: Vec<Property>,
    next_id: u32,
    /// Mesh.
    pub grid: RectilinearGrid,
    /// Named parameters used by expressions.
    pub parameters: ParameterSet,
    /// Background material.
    pub background: BackgroundMaterial,
    /// Default coordinate system of primitives.
    pub coord_system: CoordSystem,
    /// Directory against which relative file names (mesh files) are
    /// resolved.
    pub base_dir: Option<PathBuf>,
    /// Unrecognised attributes preserved from files.
    pub extra_attributes: Vec<(String, String)>,
    /// Unrecognised child elements preserved from files.
    pub extra_elements: Vec<crate::io::dom::Element>,
}

impl ContinuousStructure {
    /// Empty structure: Cartesian, no properties, empty grid with
    /// `delta_unit = 1`, vacuum background.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a property and returns its handle. Parametric values that do not
    /// depend on parameters are usable immediately; call
    /// [`ContinuousStructure::update`] after adding parametric content.
    pub fn add_property(&mut self, mut p: Property) -> PropertyId {
        self.next_id += 1;
        p.id = self.next_id;
        self.properties.push(p);
        PropertyId(self.next_id)
    }

    /// Adds a primitive to the property `id` and returns it for further
    /// configuration.
    ///
    /// ```
    /// use csx::primitives::Sphere;
    /// use csx::{ContinuousStructure, Property};
    ///
    /// let mut s = ContinuousStructure::new();
    /// let id = s.add_property(Property::metal("ball"));
    /// s.add_primitive(id, Sphere::new([0.0, 0.0, 0.0], 1.0))?.priority = 5;
    /// assert_eq!(s.property(id).unwrap().primitives[0].priority, 5);
    /// # Ok::<(), csx::Error>(())
    /// ```
    ///
    /// # Errors
    ///
    /// [`Error::Invalid`] if `id` does not refer to a property of this
    /// structure (e.g. it was removed).
    pub fn add_primitive(
        &mut self,
        id: PropertyId,
        prim: impl Into<Primitive>,
    ) -> Result<&mut Primitive> {
        let p = self
            .property_mut(id)
            .ok_or_else(|| Error::Invalid(format!("no property {id:?}")))?;
        Ok(p.add_primitive(prim))
    }

    /// Removes a property, returning it (detached, with [`Property::id`]
    /// reset to 0). Returns `None` if `id` is not in this structure.
    pub fn remove_property(&mut self, id: PropertyId) -> Option<Property> {
        let i = self.properties.iter().position(|p| p.id == id.get())?;
        let mut p = self.properties.remove(i);
        p.id = 0;
        Some(p)
    }

    /// All properties in insertion order.
    pub fn properties(&self) -> &[Property] {
        &self.properties
    }

    /// Mutable iterator over all properties.
    pub fn properties_mut(&mut self) -> impl Iterator<Item = &mut Property> {
        self.properties.iter_mut()
    }

    /// Properties matching `mask`, in insertion order.
    ///
    /// ```
    /// use csx::properties::Material;
    /// use csx::{ContinuousStructure, Property, PropertyType, TypeMask};
    ///
    /// let mut s = ContinuousStructure::new();
    /// s.add_property(Property::material("fr4", Material::new(4.3, 0.0)));
    /// s.add_property(Property::metal("pec"));
    /// assert_eq!(s.properties_of(TypeMask::MATERIALS).count(), 1);
    /// assert_eq!(s.properties_of(PropertyType::Metal).next().unwrap().name, "pec");
    /// ```
    pub fn properties_of(&self, mask: impl Into<TypeMask>) -> impl Iterator<Item = &Property> {
        let mask = mask.into();
        self.properties
            .iter()
            .filter(move |p| mask.contains(p.property_type()))
    }

    /// Property by handle.
    pub fn property(&self, id: PropertyId) -> Option<&Property> {
        self.properties.iter().find(|p| p.id == id.get())
    }

    /// Mutable property by handle.
    pub fn property_mut(&mut self, id: PropertyId) -> Option<&mut Property> {
        self.properties.iter_mut().find(|p| p.id == id.get())
    }

    /// First property with the given name.
    pub fn property_by_name(&self, name: &str) -> Option<&Property> {
        self.properties.iter().find(|p| p.name == name)
    }

    /// Mutable first property with the given name.
    pub fn property_by_name_mut(&mut self, name: &str) -> Option<&mut Property> {
        self.properties.iter_mut().find(|p| p.name == name)
    }

    /// Handle of the first property with the given name.
    pub fn property_id(&self, name: &str) -> Option<PropertyId> {
        self.property_by_name(name).map(|p| PropertyId(p.id))
    }

    /// Iterates over all `(property, primitive)` pairs.
    pub fn primitives(&self) -> impl Iterator<Item = (&Property, &Primitive)> {
        self.properties
            .iter()
            .flat_map(|p| p.primitives.iter().map(move |q| (p, q)))
    }

    /// Re-evaluates every parametric value against [`Self::parameters`] and
    /// loads external mesh files of `PolyhedronReader` primitives that are
    /// not loaded yet (relative paths are resolved against
    /// [`Self::base_dir`]).
    ///
    /// # Errors
    ///
    /// * [`Error::UnknownVariable`] / [`Error::Parse`] if an expression
    ///   refers to an undefined parameter or is malformed.
    /// * [`Error::Invalid`] for values outside their domain (negative radius,
    ///   singular transformation, NaN result, ...), prefixed with the
    ///   property name.
    /// * [`Error::Io`] / [`Error::Format`] if a mesh file cannot be read.
    ///
    /// Values updated before the failing one keep their new value.
    pub fn update(&mut self) -> Result<()> {
        self.update_parameters()?;
        for p in &mut self.properties {
            for prim in &mut p.primitives {
                if let Shape::PolyhedronReader(r) = &mut prim.shape
                    && !r.is_loaded()
                {
                    r.load(self.base_dir.as_deref())?;
                }
            }
        }
        Ok(())
    }

    /// Re-evaluates every parametric value against [`Self::parameters`]
    /// without touching external files.
    ///
    /// # Errors
    ///
    /// Same as [`Self::update`], except for the mesh file errors.
    pub fn update_parameters(&mut self) -> Result<()> {
        let params = &self.parameters;
        self.background.update(params)?;
        for p in &mut self.properties {
            p.update(params)?;
        }
        Ok(())
    }

    /// Finds the highest-priority primitive containing `p` (expressed in
    /// `cs`) among properties matching `mask`. On equal priority the
    /// primitive of the earlier property (then the earlier primitive) wins,
    /// as in CSXCAD.
    ///
    /// Primitives without their own coordinate system use
    /// [`Self::coord_system`]. Use [`Self::spatial_index`] for many queries.
    ///
    /// ```
    /// use csx::primitives::BoxShape;
    /// use csx::properties::Material;
    /// use csx::{ContinuousStructure, CoordSystem, Primitive, Property, TypeMask, Vec3};
    ///
    /// let mut s = ContinuousStructure::new();
    /// let low = s.add_property(Property::material("low", Material::new(2.0, 0.0)));
    /// let high = s.add_property(Property::material("high", Material::new(9.0, 0.0)));
    /// s.add_primitive(low, Primitive::new(BoxShape::new([0.0; 3], [2.0; 3])).with_priority(1))?;
    /// s.add_primitive(high, Primitive::new(BoxShape::new([1.0; 3], [3.0; 3])).with_priority(2))?;
    ///
    /// let p = Vec3::new(1.5, 1.5, 1.5);
    /// let hit = s.property_at(p, CoordSystem::Cartesian, TypeMask::ANY).unwrap();
    /// assert_eq!(hit.property.name, "high");
    /// assert_eq!(hit.primitive_index, 0);
    /// assert_eq!(s.properties_at(p, CoordSystem::Cartesian, TypeMask::ANY).len(), 2);
    /// // queries may be given in cylindrical coordinates (ρ, α, z)
    /// let cyl = Vec3::new(0.5, 0.0, 0.5);
    /// let hit = s.property_at(cyl, CoordSystem::Cylindrical, TypeMask::ANY).unwrap();
    /// assert_eq!(hit.property.name, "low");
    /// # Ok::<(), csx::Error>(())
    /// ```
    pub fn property_at(
        &self,
        p: Vec3,
        cs: CoordSystem,
        mask: impl Into<TypeMask>,
    ) -> Option<Hit<'_>> {
        let mask = mask.into();
        let mut best: Option<Hit<'_>> = None;
        for prop in self
            .properties
            .iter()
            .filter(|q| mask.contains(q.property_type()))
        {
            for (i, prim) in prop.primitives.iter().enumerate() {
                if best.is_some_and(|b| b.primitive.priority >= prim.priority) {
                    continue;
                }
                if prim.contains_in(p, cs, self.coord_system) {
                    best = Some(Hit {
                        property: prop,
                        primitive: prim,
                        primitive_index: i,
                    });
                }
            }
        }
        best
    }

    /// Every primitive (of properties matching `mask`) that contains `p`,
    /// highest priority first; ties keep insertion order. The first entry is
    /// the result of [`Self::property_at`].
    pub fn properties_at(
        &self,
        p: Vec3,
        cs: CoordSystem,
        mask: impl Into<TypeMask>,
    ) -> Vec<Hit<'_>> {
        let mask = mask.into();
        let mut hits: Vec<Hit<'_>> = Vec::new();
        for prop in self
            .properties
            .iter()
            .filter(|q| mask.contains(q.property_type()))
        {
            for (i, prim) in prop.primitives.iter().enumerate() {
                if prim.contains_in(p, cs, self.coord_system) {
                    hits.push(Hit {
                        property: prop,
                        primitive: prim,
                        primitive_index: i,
                    });
                }
            }
        }
        // stable sort keeps insertion order, so earlier ones win ties
        hits.sort_by_key(|h| core::cmp::Reverse(h.primitive.priority));
        hits
    }

    /// The physical medium at `p`: the highest-priority material/metal
    /// primitive, or the background material. Metals are reported as
    /// `None`. Weighting functions are evaluated at the Cartesian position
    /// of `p`.
    ///
    /// ```
    /// use csx::primitives::BoxShape;
    /// use csx::properties::Material;
    /// use csx::{ContinuousStructure, CoordSystem, Property, Vec3};
    ///
    /// let mut m = Material::new(2.0, 0.01);
    /// m.epsilon.set_weight("1 + x"); // ε grows along x
    /// let mut s = ContinuousStructure::new();
    /// let id = s.add_property(Property::material("graded", m));
    /// s.add_primitive(id, BoxShape::new([0.0; 3], [1.0; 3]))?;
    ///
    /// let med = s.medium_at(Vec3::new(0.5, 0.5, 0.5), CoordSystem::Cartesian).unwrap();
    /// assert_eq!(med.epsilon, [3.0; 3]);
    /// assert_eq!(med.kappa, [0.01; 3]);
    /// // outside: background (vacuum)
    /// let air = s.medium_at(Vec3::new(5.0, 0.0, 0.0), CoordSystem::Cartesian).unwrap();
    /// assert_eq!(air.epsilon, [1.0; 3]);
    /// # Ok::<(), csx::Error>(())
    /// ```
    pub fn medium_at(&self, p: Vec3, cs: CoordSystem) -> Option<Medium> {
        let hit = self.property_at(
            p,
            cs,
            TypeMask::MATERIALS | TypeMask::of(crate::PropertyType::Metal),
        );
        let cart = cs.to_cartesian(p);
        match hit.map(|h| &h.property.kind) {
            Some(PropertyKind::Metal) => None,
            Some(k) if let Some(m) = k.material() => Some(Medium {
                epsilon: m.values_at(MaterialQuantity::Epsilon, cart),
                mue: m.values_at(MaterialQuantity::Mue, cart),
                kappa: m.values_at(MaterialQuantity::Kappa, cart),
                sigma: m.values_at(MaterialQuantity::Sigma, cart),
            }),
            // No hit (the mask only matches materials and metals).
            _ => {
                let b = &self.background;
                Some(Medium {
                    epsilon: [b.epsilon.value(); 3],
                    mue: [b.mue.value(); 3],
                    kappa: [b.kappa.value(); 3],
                    sigma: [b.sigma.value(); 3],
                })
            }
        }
    }

    /// Union of the bounding boxes of all primitives (in `cs`). Unbounded
    /// primitives (user-defined regions) are ignored; an empty structure
    /// gives [`BoundingBox::EMPTY`].
    pub fn bounding_box(&self, cs: CoordSystem) -> BoundingBox {
        self.primitives()
            .map(|(_, q)| q.bounding_box_in(cs, self.coord_system))
            .filter(|b| b.min.is_finite() && b.max.is_finite())
            .fold(BoundingBox::EMPTY, |a, b| a.union(&b))
    }

    /// Primitives (of properties matching `mask`) whose bounding box in `cs`
    /// intersects `region`.
    pub fn primitives_in<'a>(
        &'a self,
        region: &'a BoundingBox,
        cs: CoordSystem,
        mask: impl Into<TypeMask>,
    ) -> impl Iterator<Item = (&'a Property, &'a Primitive)> + 'a {
        let mask = mask.into();
        let dcs = self.coord_system;
        self.primitives()
            .filter(move |(p, _)| mask.contains(p.property_type()))
            .filter(move |(_, q)| q.bounding_box_in(cs, dcs).intersects(region))
    }

    /// Builds an index for fast repeated point queries in coordinate system
    /// `cs`, restricted to properties matching `mask`.
    ///
    /// The index borrows the structure; rebuild it after modifying the
    /// model.
    ///
    /// ```
    /// use csx::primitives::BoxShape;
    /// use csx::{ContinuousStructure, CoordSystem, Property, TypeMask, Vec3};
    ///
    /// let mut s = ContinuousStructure::new();
    /// for i in 0..100 {
    ///     let id = s.add_property(Property::metal(format!("m{i}")));
    ///     let x = i as f64;
    ///     s.add_primitive(id, BoxShape::new([x, 0.0, 0.0], [x + 0.5, 1.0, 1.0]))?;
    /// }
    /// let index = s.spatial_index(CoordSystem::Cartesian, TypeMask::ANY);
    /// assert_eq!(index.query(Vec3::new(42.25, 0.5, 0.5)).unwrap().property.name, "m42");
    /// assert!(index.query(Vec3::new(42.75, 0.5, 0.5)).is_none());
    /// # Ok::<(), csx::Error>(())
    /// ```
    #[must_use]
    pub fn spatial_index(&self, cs: CoordSystem, mask: impl Into<TypeMask>) -> SpatialIndex<'_> {
        SpatialIndex::new(self, cs, mask.into())
    }

    /// Consistency checks: duplicate names, properties without primitives,
    /// primitives with empty bounds. Returns the issues found (empty if the
    /// model is consistent); they are warnings, not errors.
    ///
    /// ```
    /// use csx::primitives::BoxShape;
    /// use csx::{ContinuousStructure, Property};
    ///
    /// let mut s = ContinuousStructure::new();
    /// let id = s.add_property(Property::metal("pec"));
    /// s.add_primitive(id, BoxShape::new([0.0; 3], [1.0; 3]))?;
    /// assert!(s.validate().is_empty());
    ///
    /// s.add_property(Property::metal("pec"));
    /// for issue in s.validate() {
    ///     println!("warning: {issue}");
    /// }
    /// assert_eq!(s.validate().len(), 2); // duplicate name + no primitives
    /// # Ok::<(), csx::Error>(())
    /// ```
    #[must_use]
    pub fn validate(&self) -> Vec<ValidationIssue> {
        let mut w = Vec::new();
        for (i, p) in self.properties.iter().enumerate() {
            if self.properties[..i].iter().any(|q| q.name == p.name) {
                w.push(ValidationIssue::DuplicateName {
                    name: p.name.clone(),
                });
            }
            if p.primitives.is_empty() {
                w.push(ValidationIssue::EmptyProperty {
                    name: p.name.clone(),
                });
            }
            for (k, q) in p.primitives.iter().enumerate() {
                let b = q.bounding_box(self.coord_system);
                if b.is_empty() {
                    w.push(ValidationIssue::EmptyPrimitive {
                        property: p.name.clone(),
                        index: k,
                        shape: q.shape.type_name(),
                    });
                }
            }
        }
        w
    }
}

struct Entry<'a> {
    bbox: BoundingBox,
    hit: Hit<'a>,
}

/// Point-query accelerator over a borrowed structure.
///
/// Primitives are ranked by priority (highest first, ties: first added
/// first) and their bounding boxes are stored in a bounding volume
/// hierarchy, so a query only runs exact inside tests on the few primitives
/// whose boxes contain the point, in rank order.
pub struct SpatialIndex<'a> {
    entries: Vec<Entry<'a>>,
    bvh: Bvh,
    cs: CoordSystem,
    default_cs: CoordSystem,
}

impl core::fmt::Debug for SpatialIndex<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("SpatialIndex")
            .field("primitives", &self.entries.len())
            .field("cs", &self.cs)
            .field("default_cs", &self.default_cs)
            .finish_non_exhaustive()
    }
}

impl<'a> SpatialIndex<'a> {
    fn new(s: &'a ContinuousStructure, cs: CoordSystem, mask: TypeMask) -> Self {
        let mut entries: Vec<Entry<'a>> = Vec::new();
        for prop in s
            .properties
            .iter()
            .filter(|q| mask.contains(q.property_type()))
        {
            for (i, prim) in prop.primitives.iter().enumerate() {
                let b = prim.bounding_box(s.coord_system);
                // Pad by a tolerance comparable to the one used by the
                // shapes' own inside tests (relative to size and position).
                let scale = b.size().norm().max(b.min.norm()).max(b.max.norm());
                let bbox = if b.is_empty() {
                    b
                } else {
                    b.expanded(1e-9 * scale)
                };
                entries.push(Entry {
                    bbox,
                    hit: Hit {
                        property: prop,
                        primitive: prim,
                        primitive_index: i,
                    },
                });
            }
        }
        entries.sort_by_key(|e| core::cmp::Reverse(e.hit.primitive.priority));
        let boxes: Vec<BoundingBox> = entries.iter().map(|e| e.bbox).collect();
        SpatialIndex {
            bvh: Bvh::build(&boxes),
            entries,
            cs,
            default_cs: s.coord_system,
        }
    }

    /// The winning property at `p` (expressed in the index's coordinate
    /// system); same semantics as [`ContinuousStructure::property_at`].
    pub fn query(&self, p: Vec3) -> Option<Hit<'a>> {
        let c = self.cs.to_cartesian(p);
        // Candidate ranks (indices into `entries`), kept on the stack when
        // there are few of them.
        let mut buf = [0u32; 32];
        let mut n = 0usize;
        let mut spill: Vec<u32> = Vec::new();
        self.bvh.visit(
            |b| b.contains(c),
            |i| {
                if self.entries[i as usize].bbox.contains(c) {
                    if n < buf.len() {
                        buf[n] = i;
                        n += 1;
                    } else {
                        spill.push(i);
                    }
                }
                true
            },
        );
        let cand: &mut [u32] = if spill.is_empty() {
            &mut buf[..n]
        } else {
            spill.extend_from_slice(&buf[..n]);
            &mut spill
        };
        cand.sort_unstable();
        cand.iter()
            .map(|&i| &self.entries[i as usize])
            .find(|e| e.hit.primitive.contains_in(p, self.cs, self.default_cs))
            .map(|e| e.hit)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{BoxShape, Sphere};
    use crate::properties::Material;

    fn model() -> ContinuousStructure {
        let mut s = ContinuousStructure::new();
        let sub = s.add_property(Property::material("sub", Material::new(4.0, 0.0)));
        s.add_primitive(
            sub,
            Primitive::new(BoxShape::new([0.0, 0.0, 0.0], [10.0, 10.0, 1.0])).with_priority(1),
        )
        .unwrap();
        let m = s.add_property(Property::metal("pec"));
        s.add_primitive(
            m,
            Primitive::new(Sphere::new([5.0, 5.0, 1.0], 0.5)).with_priority(5),
        )
        .unwrap();
        s
    }

    /// On equal priority the earlier property wins (verified in CSXCAD).
    #[test]
    fn priority_ties() {
        let mut s = ContinuousStructure::new();
        let a = s.add_property(Property::material("A", Material::new(2.0, 0.0)));
        let b = s.add_property(Property::material("B", Material::new(3.0, 0.0)));
        let unit =
            || Primitive::new(BoxShape::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])).with_priority(3);
        s.add_primitive(b, unit()).unwrap();
        s.add_primitive(a, unit()).unwrap();
        let p = Vec3::new(0.5, 0.5, 0.5);
        let hit = s
            .property_at(p, CoordSystem::Cartesian, TypeMask::ANY)
            .unwrap();
        assert_eq!(hit.property.name, "A");
        let all = s.properties_at(p, CoordSystem::Cartesian, TypeMask::ANY);
        assert_eq!(all[0].property.name, "A");
        let idx = s.spatial_index(CoordSystem::Cartesian, TypeMask::ANY);
        assert_eq!(idx.query(p).unwrap().property.name, "A");
    }

    #[test]
    fn priority() {
        let s = model();
        let h = s
            .property_at(
                Vec3::new(5.0, 5.0, 0.8),
                CoordSystem::Cartesian,
                TypeMask::ANY,
            )
            .unwrap();
        assert_eq!(h.property.name, "pec");
        let h = s
            .property_at(
                Vec3::new(1.0, 1.0, 0.5),
                CoordSystem::Cartesian,
                TypeMask::ANY,
            )
            .unwrap();
        assert_eq!(h.property.name, "sub");
        assert!(
            s.property_at(
                Vec3::new(1.0, 1.0, 2.0),
                CoordSystem::Cartesian,
                TypeMask::ANY
            )
            .is_none()
        );
        assert_eq!(
            s.properties_at(
                Vec3::new(5.0, 5.0, 0.8),
                CoordSystem::Cartesian,
                TypeMask::ANY
            )
            .len(),
            2
        );
        let idx = s.spatial_index(CoordSystem::Cartesian, TypeMask::ANY);
        for p in [
            Vec3::new(5.0, 5.0, 0.8),
            Vec3::new(1.0, 1.0, 0.5),
            Vec3::new(1.0, 1.0, 2.0),
        ] {
            assert_eq!(
                idx.query(p).map(|h| &h.property.name),
                s.property_at(p, CoordSystem::Cartesian, TypeMask::ANY)
                    .map(|h| &h.property.name)
            );
        }
    }

    #[test]
    fn media() {
        let s = model();
        assert_eq!(
            s.medium_at(Vec3::new(1.0, 1.0, 0.5), CoordSystem::Cartesian)
                .unwrap()
                .epsilon,
            [4.0; 3]
        );
        assert_eq!(
            s.medium_at(Vec3::new(1.0, 1.0, 5.0), CoordSystem::Cartesian)
                .unwrap()
                .epsilon,
            [1.0; 3]
        );
        assert!(
            s.medium_at(Vec3::new(5.0, 5.0, 1.0), CoordSystem::Cartesian)
                .is_none()
        );
        let b = s.bounding_box(CoordSystem::Cartesian);
        assert_eq!(b.max, Vec3::new(10.0, 10.0, 1.5));
    }

    #[test]
    fn parametric() {
        let mut s = ContinuousStructure::new();
        s.parameters.set("w", 2.0);
        let id = s.add_property(Property::metal("m"));
        s.add_primitive(
            id,
            BoxShape {
                start: [0.0, 0.0, 0.0].into(),
                stop: crate::ScalarPoint(["w".into(), "w*2".into(), 1.0.into()]),
            },
        )
        .unwrap();
        s.update().unwrap();
        assert!(
            s.property_at(
                Vec3::new(1.9, 3.9, 0.5),
                CoordSystem::Cartesian,
                TypeMask::ANY
            )
            .is_some()
        );
        s.parameters.set("w", 1.0);
        s.update().unwrap();
        assert!(
            s.property_at(
                Vec3::new(1.9, 3.9, 0.5),
                CoordSystem::Cartesian,
                TypeMask::ANY
            )
            .is_none()
        );
        assert!(s.validate().is_empty());
        assert!(s.remove_property(id).is_some());
    }

    #[test]
    fn index_matches_linear_scan() {
        use crate::primitives::UserDefined;

        let mut seed = 11u64;
        let mut rnd = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut s = ContinuousStructure::new();
        for k in 0..60 {
            let id = s.add_property(Property::material(format!("m{k}"), Material::new(1.0, 0.0)));
            // Several primitives per property, few distinct priorities so
            // that ties are common.
            for _ in 0..1 + k % 3 {
                let c = Vec3::new(rnd() * 20.0, rnd() * 20.0, rnd() * 20.0);
                let d = Vec3::new(rnd() * 4.0, rnd() * 4.0, rnd() * 4.0);
                let prim = if rnd() < 0.7 {
                    Primitive::new(BoxShape::new(c - d, c + d))
                } else {
                    Primitive::new(Sphere::new(c, d[0]))
                };
                s.add_primitive(id, prim.with_priority((rnd() * 4.0) as i32))
                    .unwrap();
            }
        }
        // An unbounded primitive must be found everywhere it applies.
        let u = s.add_property(Property::metal("half"));
        s.add_primitive(
            u,
            Primitive::new(UserDefined::new("z > 15")).with_priority(2),
        )
        .unwrap();
        let idx = s.spatial_index(CoordSystem::Cartesian, TypeMask::ANY);
        for _ in 0..5000 {
            let p = Vec3::new(rnd() * 22.0 - 1.0, rnd() * 22.0 - 1.0, rnd() * 22.0 - 1.0);
            let a = idx.query(p);
            let b = s.property_at(p, CoordSystem::Cartesian, TypeMask::ANY);
            match (a, b) {
                (None, None) => {}
                (Some(a), Some(b)) => {
                    assert!(core::ptr::eq(a.primitive, b.primitive), "point {p:?}");
                    assert_eq!(a.property.id(), b.property.id());
                    assert_eq!(a.primitive_index, b.primitive_index);
                }
                _ => panic!("index and linear scan disagree at {p:?}"),
            }
        }
    }
}
