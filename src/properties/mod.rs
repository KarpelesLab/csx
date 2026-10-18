//! Properties: what a region of space *is* (material, metal, excitation,
//! probe, ...). Every property owns the primitives describing *where* it is.
//!
//! A [`Property`] combines a name, visualisation settings, a kind-specific
//! payload ([`PropertyKind`]) and its [`Primitive`]s. Each kind maps to one
//! element of the CSXCAD file format:
//!
//! | kind | XML element | payload | meaning |
//! |---|---|---|---|
//! | [`PropertyKind::Material`] | `<Material>` | [`Material`] | dielectric / magnetic / lossy material |
//! | [`PropertyKind::LorentzMaterial`] | `<LorentzMaterial>` | [`LorentzMaterial`] | Drude / Lorentz dispersive material |
//! | [`PropertyKind::DebyeMaterial`] | `<DebyeMaterial>` | [`DebyeMaterial`] | Debye dispersive material |
//! | [`PropertyKind::DiscMaterial`] | `<DiscMaterial>` | [`DiscMaterial`] | voxel material from a data file |
//! | [`PropertyKind::Metal`] | `<Metal>` | – | perfect electric conductor (PEC) |
//! | [`PropertyKind::ConductingSheet`] | `<ConductingSheet>` | [`ConductingSheet`] | thin sheet of finite conductivity |
//! | [`PropertyKind::LumpedElement`] | `<LumpedElement>` | [`LumpedElement`] | lumped R, L, C |
//! | [`PropertyKind::Excitation`] | `<Excitation>` | [`Excitation`] | field source |
//! | [`PropertyKind::ProbeBox`] | `<ProbeBox>` | [`ProbeBox`] | voltage / current / field probe |
//! | [`PropertyKind::DumpBox`] | `<DumpBox>` | [`DumpBox`] | field dump |
//! | [`PropertyKind::AbsorbingBC`] | `<AbsorbingBC>` | [`AbsorbingBC`] | absorbing boundary sheet |
//! | [`PropertyKind::Unknown`] | any other | – | preserved for round-trips |
//!
//! [`PropertyType`] flags identify kinds; a [`TypeMask`] combines them to
//! filter structure queries.
//!
//! # Material quantities
//!
//! | quantity | field | file attribute | unit | default |
//! |---|---|---|---|---|
//! | relative permittivity εᵣ | [`Material::epsilon`] | `Epsilon` | – | 1 |
//! | relative permeability μᵣ | [`Material::mue`] | `Mue` | – | 1 |
//! | electric conductivity κ | [`Material::kappa`] | `Kappa` | S/m | 0 |
//! | magnetic conductivity σ | [`Material::sigma`] | `Sigma` | Ω/m | 0 |
//! | mass density ρₘ | [`Material::density`] | `Density` | kg/m³ | 0 |
//!
//! Each quantity has a value and a weighting function per axis
//! ([`Aniso`]); with [`Material::isotropic`] set (the default) only the
//! first axis is used. In a cylindrical model the three axes are ρ, α, z.
//!
//! # Weighting functions
//!
//! A [`WeightFn`] multiplies a base value by an expression of the position
//! (in drawing units):
//!
//! | variable | meaning |
//! |---|---|
//! | `x`, `y`, `z` | Cartesian coordinates |
//! | `rho` | distance to the z axis, √(x² + y²) |
//! | `a` | azimuth atan2(y, x) |
//! | `r` | distance to the origin |
//! | `t` | polar angle from +z |
//!
//! # Dispersive materials
//!
//! With angular frequency ω = 2πf (all model frequencies are given in Hz):
//!
//! * Drude / Lorentz pole ([`LorentzPole`]):
//!   `ε(ω) = ε∞ − ω_p² / (ω² − ω_L² − iω/τ)` with plasma frequency
//!   `ω_p = 2π·f_p`, pole frequency `ω_L = 2π·f_L` (0 for a Drude pole) and
//!   relaxation time τ. The magnetic part is analogous.
//! * Debye pole ([`DebyePole`]): `ε(ω) = ε∞ + Δε / (1 + iωτ)`.
//!
//! Several poles add up; ε∞ is the static [`Material::epsilon`].
//!
//! # Type identifiers
//!
//! Settings stored as numbers in files are typed enums with an `Other(i32)`
//! fallback so unknown ids round-trip losslessly. Build them from numbers
//! with `from_id`: an `Other` holding a known id compares unequal to the
//! named variant.
//!
//! All of them provide `from_id` / `id`, `name`, `ALL`, `Display` and
//! `FromStr` (symbolic name, case-insensitive with `_`/`-` ignored, or a
//! number, also float-formatted as some writers produce):
//!
//! ```
//! use csx::properties::{DumpType, ExcitationType};
//!
//! assert_eq!(ExcitationType::from_id(1), ExcitationType::EFieldHard);
//! assert_eq!(ExcitationType::from_id(42), ExcitationType::Other(42));
//! assert_eq!("plane_wave".parse::<ExcitationType>()?, ExcitationType::PlaneWave);
//! assert_eq!("1.000000e+01".parse::<DumpType>()?, DumpType::EFieldFd);
//! assert_eq!(DumpType::Sar1g.id(), 21);
//! assert_eq!(DumpType::Sar1g.to_string(), "SAR1g");
//! # Ok::<(), csx::Error>(())
//! ```
//!
//! | enum | file attribute | ids |
//! |---|---|---|
//! | [`ExcitationType`] | `<Excitation Type>` | 0 E soft, 1 E hard, 2 H soft, 3 H hard, 10 plane wave |
//! | [`ProbeType`] | `<ProbeBox Type>` | 0 voltage, 1 current, 2 E field, 3 H field, 10/11 waveguide voltage/current |
//! | [`DumpType`] | `<DumpBox DumpType>` | 0–5 E, H, J, rot H, D, B (time domain); 10–15 same in frequency domain; 20, 21, 22, 29 SAR |
//! | [`DumpMode`] | `<DumpBox DumpMode>` | 0 none, 1 node, 2 cell |
//! | [`DumpFileType`] | `<DumpBox FileType>` | 0 VTK, 1 HDF5 |
//! | [`LumpedElementType`] | `<LumpedElement LEtype>` | 0 parallel, 1 series |
//! | [`AbcType`] | `<AbsorbingBC AbsorbingBoundaryType>` | 0 Mur 1st, 1 Mur 1st SA |
//! | [`Axis`] | `Direction`, `NormDir` | 0 x/ρ, 1 y/α, 2 z |
//!
//! # Example
//!
//! ```
//! use csx::primitives::BoxShape;
//! use csx::properties::{
//!     Axis, DumpBox, DumpType, Excitation, LumpedElement, Material, ProbeBox, ProbeType,
//! };
//! use csx::{ContinuousStructure, Property, PropertyType, TypeMask};
//!
//! let mut s = ContinuousStructure::new();
//! s.add_property(
//!     Property::material("FR4", Material::new(4.3, 0.0))
//!         .with_primitive(BoxShape::new([-10.0, -10.0, 0.0], [10.0, 10.0, 1.6])),
//! );
//! s.add_property(Property::metal("patch"));
//! s.add_property(Property::lumped_element("port_R", LumpedElement::resistor(Axis::Z, 50.0)));
//! s.add_property(Property::excitation("port_E", Excitation::new([0.0, 0.0, 1.0])));
//! s.add_property(Property::probe("port_U", ProbeBox::new(ProbeType::Voltage)));
//! s.add_property(Property::dump("Et", DumpBox::new(DumpType::EField)));
//!
//! assert_eq!(s.properties_of(TypeMask::PHYSICAL).count(), 3);
//! assert_eq!(s.properties_of(PropertyType::Excitation).count(), 1);
//! ```

mod kinds;
mod types;
mod weight;

pub use kinds::*;
pub use types::{
    AbcType, Axis, DumpFileType, DumpMode, DumpType, ExcitationType, LumpedElementType, ProbeType,
};
pub use weight::{Aniso, WEIGHT_VARS, WeightFn, weight_vars};

use crate::error::Result;
use crate::param::{ParameterSet, Scalar};
use crate::primitives::Primitive;

/// An RGBA colour (`<FillColor>` / `<EdgeColor>` in files).
///
/// ```
/// use csx::properties::Rgba;
///
/// let red = Rgba::new(255, 0, 0, 255);
/// assert_eq!(Rgba::default(), Rgba::new(0, 0, 0, 255)); // opaque black
/// assert_eq!(Rgba::from_seed(7), Rgba::from_seed(7)); // deterministic
/// assert_ne!(Rgba::from_seed(7), red);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Rgba {
    /// Red.
    pub r: u8,
    /// Green.
    pub g: u8,
    /// Blue.
    pub b: u8,
    /// Alpha (255 = opaque).
    pub a: u8,
}

impl Default for Rgba {
    fn default() -> Self {
        Rgba {
            r: 0,
            g: 0,
            b: 0,
            a: 255,
        }
    }
}

impl Rgba {
    /// Colour from components.
    pub const fn new(r: u8, g: u8, b: u8, a: u8) -> Self {
        Rgba { r, g, b, a }
    }

    /// A deterministic, reasonably distinct opaque colour derived from `seed`.
    pub fn from_seed(seed: u32) -> Self {
        let mut x = seed.wrapping_mul(0x9E37_79B9).wrapping_add(0x7F4A_7C15);
        x ^= x >> 16;
        x = x.wrapping_mul(0x85EB_CA6B);
        x ^= x >> 13;
        Rgba {
            r: x as u8,
            g: (x >> 8) as u8,
            b: (x >> 16) as u8,
            a: 255,
        }
    }
}

/// The kind of a property, as a bit flag (see [`TypeMask`]).
///
/// ```
/// use csx::{Property, PropertyType, TypeMask};
///
/// let p = Property::metal("pec");
/// assert_eq!(p.property_type(), PropertyType::Metal);
/// let mask = PropertyType::Metal | PropertyType::ConductingSheet;
/// assert!(mask.contains(p.property_type()));
/// ```
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u32)]
pub enum PropertyType {
    /// Unknown / user-defined.
    Unknown = 0x0001,
    /// Plain material.
    Material = 0x0002,
    /// Perfect electric conductor.
    Metal = 0x0004,
    /// Excitation.
    Excitation = 0x0008,
    /// Probe.
    ProbeBox = 0x0010,
    /// Field dump.
    DumpBox = 0x0020,
    /// Drude/Lorentz dispersive material.
    LorentzMaterial = 0x0040,
    /// Debye dispersive material.
    DebyeMaterial = 0x0080,
    /// Discrete (voxel) material.
    DiscMaterial = 0x0100,
    /// Conducting sheet.
    ConductingSheet = 0x0200,
    /// Lumped element.
    LumpedElement = 0x0400,
    /// Absorbing boundary condition.
    AbsorbingBC = 0x0800,
}

/// A set of [`PropertyType`]s used to filter queries.
///
/// Masks are built from the constants, from single types
/// ([`TypeMask::of`], `From<PropertyType>`) and with `|`:
///
/// ```
/// use csx::{PropertyType, TypeMask};
///
/// let mut m = TypeMask::MATERIALS | PropertyType::Metal;
/// assert!(m.contains(PropertyType::DebyeMaterial));
/// assert!(m.contains(PropertyType::Metal));
/// assert!(!m.contains(PropertyType::ProbeBox));
/// m |= PropertyType::ProbeBox;
/// assert!(m.contains(PropertyType::ProbeBox));
/// assert!(!TypeMask::NONE.contains(PropertyType::Metal));
/// assert!(TypeMask::ANY.contains(PropertyType::Unknown));
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TypeMask(pub u32);

impl TypeMask {
    /// Matches every property.
    pub const ANY: TypeMask = TypeMask(u32::MAX);
    /// Matches no property.
    pub const NONE: TypeMask = TypeMask(0);
    /// Every kind of material (plain, dispersive, discrete).
    pub const MATERIALS: TypeMask = TypeMask(
        PropertyType::Material as u32
            | PropertyType::LorentzMaterial as u32
            | PropertyType::DebyeMaterial as u32
            | PropertyType::DiscMaterial as u32,
    );
    /// Every property that affects the physical structure (materials,
    /// metals, sheets, lumped elements).
    pub const PHYSICAL: TypeMask = TypeMask(
        TypeMask::MATERIALS.0
            | PropertyType::Metal as u32
            | PropertyType::ConductingSheet as u32
            | PropertyType::LumpedElement as u32,
    );

    /// Mask for a single type.
    pub const fn of(t: PropertyType) -> TypeMask {
        TypeMask(t as u32)
    }

    /// True if `t` is in the mask.
    pub const fn contains(self, t: PropertyType) -> bool {
        self.0 & t as u32 != 0
    }
}

impl From<PropertyType> for TypeMask {
    fn from(t: PropertyType) -> Self {
        TypeMask::of(t)
    }
}

impl core::ops::BitOr for TypeMask {
    type Output = TypeMask;
    fn bitor(self, o: TypeMask) -> TypeMask {
        TypeMask(self.0 | o.0)
    }
}

impl core::ops::BitOr<PropertyType> for TypeMask {
    type Output = TypeMask;
    fn bitor(self, o: PropertyType) -> TypeMask {
        TypeMask(self.0 | o as u32)
    }
}

impl core::ops::BitOrAssign for TypeMask {
    fn bitor_assign(&mut self, o: TypeMask) {
        self.0 |= o.0;
    }
}

impl core::ops::BitOrAssign<PropertyType> for TypeMask {
    fn bitor_assign(&mut self, o: PropertyType) {
        self.0 |= o as u32;
    }
}

impl core::ops::BitOr for PropertyType {
    type Output = TypeMask;
    fn bitor(self, o: PropertyType) -> TypeMask {
        TypeMask(self as u32 | o as u32)
    }
}

/// Kind-specific data of a property; see the [module docs](self) for the
/// mapping to file elements.
///
/// ```
/// use csx::properties::{Material, PropertyKind};
///
/// let k = PropertyKind::Material(Material::new(2.2, 0.0));
/// assert_eq!(k.type_name(), "Material");
/// assert_eq!(k.material().unwrap().epsilon.base(0), 2.2);
/// assert!(PropertyKind::Metal.material().is_none());
/// ```
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum PropertyKind {
    /// Unknown property; `property` holds the free-form type text.
    Unknown {
        /// Free-form description.
        property: String,
        /// Element name the property was read from when it was not
        /// `<Unknown>` (e.g. a property type this crate does not model);
        /// written back under the same name.
        element: Option<String>,
    },
    /// Plain material.
    Material(Material),
    /// Perfect electric conductor.
    Metal,
    /// Excitation.
    Excitation(Excitation),
    /// Probe.
    ProbeBox(ProbeBox),
    /// Field dump.
    DumpBox(DumpBox),
    /// Drude/Lorentz dispersive material.
    LorentzMaterial(LorentzMaterial),
    /// Debye dispersive material.
    DebyeMaterial(DebyeMaterial),
    /// Discrete material.
    DiscMaterial(DiscMaterial),
    /// Conducting sheet.
    ConductingSheet(ConductingSheet),
    /// Lumped element.
    LumpedElement(LumpedElement),
    /// Absorbing boundary condition.
    AbsorbingBC(AbsorbingBC),
}

impl PropertyKind {
    /// The type flag.
    pub fn property_type(&self) -> PropertyType {
        match self {
            PropertyKind::Unknown { .. } => PropertyType::Unknown,
            PropertyKind::Material(_) => PropertyType::Material,
            PropertyKind::Metal => PropertyType::Metal,
            PropertyKind::Excitation(_) => PropertyType::Excitation,
            PropertyKind::ProbeBox(_) => PropertyType::ProbeBox,
            PropertyKind::DumpBox(_) => PropertyType::DumpBox,
            PropertyKind::LorentzMaterial(_) => PropertyType::LorentzMaterial,
            PropertyKind::DebyeMaterial(_) => PropertyType::DebyeMaterial,
            PropertyKind::DiscMaterial(_) => PropertyType::DiscMaterial,
            PropertyKind::ConductingSheet(_) => PropertyType::ConductingSheet,
            PropertyKind::LumpedElement(_) => PropertyType::LumpedElement,
            PropertyKind::AbsorbingBC(_) => PropertyType::AbsorbingBC,
        }
    }

    /// Element name used in CSXCAD files.
    pub fn type_name(&self) -> &'static str {
        match self {
            PropertyKind::Unknown { .. } => "Unknown",
            PropertyKind::Material(_) => "Material",
            PropertyKind::Metal => "Metal",
            PropertyKind::Excitation(_) => "Excitation",
            PropertyKind::ProbeBox(_) => "ProbeBox",
            PropertyKind::DumpBox(_) => "DumpBox",
            PropertyKind::LorentzMaterial(_) => "LorentzMaterial",
            PropertyKind::DebyeMaterial(_) => "DebyeMaterial",
            PropertyKind::DiscMaterial(_) => "DiscMaterial",
            PropertyKind::ConductingSheet(_) => "ConductingSheet",
            PropertyKind::LumpedElement(_) => "LumpedElement",
            PropertyKind::AbsorbingBC(_) => "AbsorbingBC",
        }
    }

    /// The (non-dispersive part of the) material, for material kinds
    /// (plain, Lorentz, Debye, discrete); `None` otherwise.
    pub fn material(&self) -> Option<&Material> {
        match self {
            PropertyKind::Material(m) => Some(m),
            PropertyKind::LorentzMaterial(m) => Some(&m.material),
            PropertyKind::DebyeMaterial(m) => Some(&m.material),
            PropertyKind::DiscMaterial(m) => Some(&m.material),
            _ => None,
        }
    }

    /// Mutable access to the material, for material kinds.
    pub fn material_mut(&mut self) -> Option<&mut Material> {
        match self {
            PropertyKind::Material(m) => Some(m),
            PropertyKind::LorentzMaterial(m) => Some(&mut m.material),
            PropertyKind::DebyeMaterial(m) => Some(&mut m.material),
            PropertyKind::DiscMaterial(m) => Some(&mut m.material),
            _ => None,
        }
    }

    fn update(&mut self, p: &ParameterSet) -> Result<()> {
        match self {
            PropertyKind::Material(m) => m.update(p),
            PropertyKind::LorentzMaterial(m) => m.update(p),
            PropertyKind::DebyeMaterial(m) => m.update(p),
            PropertyKind::DiscMaterial(m) => {
                m.material.update(p)?;
                m.scale.update(p)
            }
            PropertyKind::ConductingSheet(s) => {
                s.conductivity.update(p)?;
                s.thickness.update(p)
            }
            PropertyKind::Excitation(e) => e.update(p),
            PropertyKind::ProbeBox(b) => update_probe(b, p),
            PropertyKind::DumpBox(d) => {
                update_probe(&mut d.probe, p)?;
                if let Some(r) = &mut d.opt_resolution {
                    r.iter_mut().try_for_each(|s| s.update(p))?;
                }
                Ok(())
            }
            PropertyKind::LumpedElement(l) => l.update(p),
            PropertyKind::AbsorbingBC(a) => a.phase_velocity.update(p),
            PropertyKind::Unknown { .. } | PropertyKind::Metal => Ok(()),
        }
    }
}

fn update_probe(b: &mut ProbeBox, p: &ParameterSet) -> Result<()> {
    b.weight.update(p)?;
    b.start_time.update(p)?;
    b.stop_time.update(p)
}

/// A named property together with its primitives.
///
/// Create properties with the kind-specific constructors
/// ([`Property::material`], [`Property::metal`], [`Property::excitation`],
/// ...) or [`Property::new`], attach primitives, then add them to a
/// [`ContinuousStructure`](crate::ContinuousStructure).
///
/// ```
/// use csx::primitives::{BoxShape, Cylinder};
/// use csx::properties::{Material, Rgba};
/// use csx::Property;
///
/// let mut p = Property::material("teflon", Material::new(2.1, 0.0))
///     .with_color(Rgba::new(200, 200, 255, 128))
///     .with_primitive(BoxShape::new([0.0; 3], [1.0; 3]));
/// p.add_primitive(Cylinder::new([0.0; 3], [0.0, 0.0, 2.0], 0.5)).priority = 3;
/// p.set_attribute("vendor", "acme");
///
/// assert_eq!(p.primitives.len(), 2);
/// assert_eq!(p.attribute("vendor"), Some("acme"));
/// assert_eq!(p.id(), 0); // not part of a structure yet
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Property {
    pub(crate) id: u32,
    /// Name (should be unique within a structure).
    pub name: String,
    /// Fill colour for visualisation.
    pub fill_color: Rgba,
    /// Edge colour for visualisation.
    pub edge_color: Rgba,
    /// Visibility flag for visualisation.
    pub visible: bool,
    /// Kind-specific data.
    pub kind: PropertyKind,
    /// The primitives defining where this property applies.
    pub primitives: Vec<Primitive>,
    /// Free-form user attributes (`name → value`).
    pub attributes: Vec<(String, String)>,
    /// Unrecognised attributes preserved from files.
    pub extra_attributes: Vec<(String, String)>,
    /// Unrecognised child elements preserved from files.
    pub extra_elements: Vec<crate::io::dom::Element>,
}

impl Property {
    /// New property with default colours (opaque black), visible, and no
    /// primitives.
    ///
    /// ```
    /// use csx::{Property, PropertyKind};
    ///
    /// let p = Property::new("pec", PropertyKind::Metal);
    /// assert_eq!(p, Property::metal("pec"));
    /// ```
    pub fn new(name: impl Into<String>, kind: PropertyKind) -> Self {
        Property {
            id: 0,
            name: name.into(),
            fill_color: Rgba::default(),
            edge_color: Rgba::default(),
            visible: true,
            kind,
            primitives: Vec::new(),
            attributes: Vec::new(),
            extra_attributes: Vec::new(),
            extra_elements: Vec::new(),
        }
    }

    /// New plain material.
    ///
    /// ```
    /// use csx::properties::Material;
    /// use csx::{Property, PropertyType};
    ///
    /// let p = Property::material("Rogers 4003C", Material::new(3.38, 0.0027));
    /// assert_eq!(p.property_type(), PropertyType::Material);
    /// ```
    pub fn material(name: impl Into<String>, m: Material) -> Self {
        Property::new(name, PropertyKind::Material(m))
    }

    /// New perfect electric conductor.
    pub fn metal(name: impl Into<String>) -> Self {
        Property::new(name, PropertyKind::Metal)
    }

    /// New excitation.
    ///
    /// ```
    /// use csx::properties::{Excitation, ExcitationType};
    /// use csx::Property;
    ///
    /// // hard E-field source polarised along z, port number 1
    /// let e = Excitation::new([0.0, 0.0, 1.0])
    ///     .with_type(ExcitationType::EFieldHard)
    ///     .with_number(1);
    /// let p = Property::excitation("port1", e);
    /// assert_eq!(p.kind.type_name(), "Excitation");
    /// ```
    pub fn excitation(name: impl Into<String>, e: Excitation) -> Self {
        Property::new(name, PropertyKind::Excitation(e))
    }

    /// New probe.
    ///
    /// ```
    /// use csx::properties::{Axis, ProbeBox, ProbeType};
    /// use csx::Property;
    ///
    /// let i = Property::probe("port1_I", ProbeBox::new(ProbeType::Current).with_norm_dir(Axis::Z));
    /// assert_eq!(i.kind.type_name(), "ProbeBox");
    /// ```
    pub fn probe(name: impl Into<String>, p: ProbeBox) -> Self {
        Property::new(name, PropertyKind::ProbeBox(p))
    }

    /// New field dump.
    ///
    /// ```
    /// use csx::properties::{DumpBox, DumpFileType, DumpType};
    /// use csx::Property;
    ///
    /// let d = Property::dump("Et", DumpBox::new(DumpType::EField).with_file_type(DumpFileType::Hdf5));
    /// assert_eq!(d.kind.type_name(), "DumpBox");
    /// ```
    pub fn dump(name: impl Into<String>, d: DumpBox) -> Self {
        Property::new(name, PropertyKind::DumpBox(d))
    }

    /// New lumped element.
    ///
    /// ```
    /// use csx::properties::{Axis, LumpedElement};
    /// use csx::Property;
    ///
    /// let r = Property::lumped_element("R50", LumpedElement::resistor(Axis::Z, 50.0));
    /// assert_eq!(r.kind.type_name(), "LumpedElement");
    /// ```
    pub fn lumped_element(name: impl Into<String>, l: LumpedElement) -> Self {
        Property::new(name, PropertyKind::LumpedElement(l))
    }

    /// New conducting sheet with conductivity (S/m) and thickness (m).
    ///
    /// ```
    /// use csx::Property;
    ///
    /// let copper = Property::conducting_sheet("copper", 5.8e7, 35e-6);
    /// assert_eq!(copper.kind.type_name(), "ConductingSheet");
    /// ```
    pub fn conducting_sheet(
        name: impl Into<String>,
        conductivity: impl Into<Scalar>,
        thickness: impl Into<Scalar>,
    ) -> Self {
        Property::new(
            name,
            PropertyKind::ConductingSheet(ConductingSheet::new(conductivity, thickness)),
        )
    }

    /// New Drude/Lorentz dispersive material; see [`LorentzMaterial`].
    pub fn lorentz_material(name: impl Into<String>, m: LorentzMaterial) -> Self {
        Property::new(name, PropertyKind::LorentzMaterial(m))
    }

    /// New Debye dispersive material; see [`DebyeMaterial`].
    pub fn debye_material(name: impl Into<String>, m: DebyeMaterial) -> Self {
        Property::new(name, PropertyKind::DebyeMaterial(m))
    }

    /// New discrete (voxel file) material.
    pub fn disc_material(name: impl Into<String>, m: DiscMaterial) -> Self {
        Property::new(name, PropertyKind::DiscMaterial(m))
    }

    /// New absorbing boundary condition sheet.
    pub fn absorbing_bc(name: impl Into<String>, a: AbsorbingBC) -> Self {
        Property::new(name, PropertyKind::AbsorbingBC(a))
    }

    /// Builder: sets fill and edge colour.
    #[must_use]
    pub fn with_color(mut self, c: Rgba) -> Self {
        self.fill_color = c;
        self.edge_color = c;
        self
    }

    /// Unique id within the owning structure (0 for detached properties);
    /// equal to [`PropertyId::get`](crate::PropertyId::get) of its handle.
    pub fn id(&self) -> u32 {
        self.id
    }

    /// The type flag.
    pub fn property_type(&self) -> PropertyType {
        self.kind.property_type()
    }

    /// Adds a primitive, returning `self` for chaining.
    #[must_use]
    pub fn with_primitive(mut self, p: impl Into<Primitive>) -> Self {
        self.primitives.push(p.into());
        self
    }

    /// Adds a primitive and returns it for further configuration.
    pub fn add_primitive(&mut self, p: impl Into<Primitive>) -> &mut Primitive {
        let i = self.primitives.len();
        self.primitives.push(p.into());
        &mut self.primitives[i]
    }

    /// Value of a user attribute.
    pub fn attribute(&self, name: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(k, _)| k == name)
            .map(|(_, v)| v.as_str())
    }

    /// Sets a user attribute (written to the `<Attributes>` element),
    /// replacing an existing value of the same name.
    pub fn set_attribute(&mut self, name: impl Into<String>, value: impl Into<String>) {
        let (name, value) = (name.into(), value.into());
        match self.attributes.iter_mut().find(|(k, _)| *k == name) {
            Some(e) => e.1 = value,
            None => self.attributes.push((name, value)),
        }
    }

    /// Re-evaluates parametric values of the property and its primitives.
    /// Usually called through
    /// [`ContinuousStructure::update`](crate::ContinuousStructure::update).
    ///
    /// # Errors
    ///
    /// [`Error::UnknownVariable`](crate::Error::UnknownVariable) or
    /// [`Error::Parse`](crate::Error::Parse) for unresolvable expressions,
    /// [`Error::Invalid`](crate::Error::Invalid) (prefixed with the property
    /// name) for out-of-domain values such as a negative radius.
    pub fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.kind.update(p).map_err(|e| context(e, &self.name))?;
        for prim in &mut self.primitives {
            prim.update(p).map_err(|e| context(e, &self.name))?;
        }
        Ok(())
    }
}

fn context(e: crate::Error, name: &str) -> crate::Error {
    match e {
        crate::Error::Invalid(m) => crate::Error::Invalid(format!("property {name:?}: {m}")),
        other => other,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Vec3;

    #[test]
    fn material_values() {
        let mut m = Material::new(4.0, 0.1);
        assert_eq!(m.value_at(MaterialQuantity::Epsilon, 2, Vec3::ZERO), 4.0);
        m.epsilon.set_weight("1+x");
        assert_eq!(
            m.value_at(MaterialQuantity::Epsilon, 1, Vec3::new(1.0, 0.0, 0.0)),
            8.0
        );
        m.isotropic = false;
        m.epsilon.value[2].set(2.0);
        assert_eq!(
            m.values_at(MaterialQuantity::Epsilon, Vec3::ZERO),
            [4.0, 4.0, 2.0]
        );
        m.isotropic = true;
        assert_eq!(
            m.values_at(MaterialQuantity::Epsilon, Vec3::ZERO),
            [4.0, 4.0, 4.0]
        );
    }

    #[test]
    fn constructors() {
        let e = Property::excitation(
            "port",
            Excitation::new([0.0, 0.0, 1.0])
                .with_type(ExcitationType::EFieldHard)
                .with_number(2),
        );
        let PropertyKind::Excitation(x) = &e.kind else {
            panic!()
        };
        assert_eq!(
            (x.excitation_type, x.number, x.excite[2].value()),
            (ExcitationType::EFieldHard, 2, 1.0)
        );
        assert!(x.enabled);
        let l = LumpedElement::resistor(Axis::Z, 50.0).with_c(1e-12);
        assert!(l.l.value().is_nan() && l.c.value() == 1e-12);
        assert_eq!(
            Property::lumped_element("r", l.clone()).kind,
            PropertyKind::LumpedElement(l)
        );
        let p = Property::probe(
            "u",
            ProbeBox::new(ProbeType::Voltage)
                .with_weight(-1.0)
                .with_frequencies([1e9]),
        );
        assert_eq!(p.property_type(), PropertyType::ProbeBox);
        assert_eq!(
            Property::dump("et", DumpBox::new(DumpType::EFieldFd)).property_type(),
            PropertyType::DumpBox
        );
        let c = Property::conducting_sheet("cu", 5.8e7, 35e-6).with_color(Rgba::new(1, 2, 3, 255));
        assert_eq!(c.fill_color, c.edge_color);
        let pt: crate::ScalarPoint = ["w/2", "1", "h"].into();
        assert!(pt.0[0].is_expr() && !pt.0[1].is_expr());
    }

    #[test]
    fn masks() {
        assert!(TypeMask::MATERIALS.contains(PropertyType::DebyeMaterial));
        assert!(!TypeMask::MATERIALS.contains(PropertyType::Metal));
        assert!((PropertyType::Metal | PropertyType::Excitation).contains(PropertyType::Metal));
        let p = Property::metal("pec");
        assert_eq!(p.kind.type_name(), "Metal");
        assert!(TypeMask::PHYSICAL.contains(p.property_type()));
    }
}
