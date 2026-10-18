//! Properties: what a region of space *is* (material, metal, excitation,
//! probe, ...). Every property owns the primitives describing *where* it is.

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

/// An RGBA colour.
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TypeMask(pub u32);

impl TypeMask {
    /// Matches every property.
    pub const ANY: TypeMask = TypeMask(u32::MAX);
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

impl core::ops::BitOr for PropertyType {
    type Output = TypeMask;
    fn bitor(self, o: PropertyType) -> TypeMask {
        TypeMask(self as u32 | o as u32)
    }
}

/// Kind-specific data of a property.
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

    /// The (non-dispersive part of the) material, for material kinds.
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
    /// New property with default colours and no primitives.
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
    pub fn material(name: impl Into<String>, m: Material) -> Self {
        Property::new(name, PropertyKind::Material(m))
    }

    /// New perfect electric conductor.
    pub fn metal(name: impl Into<String>) -> Self {
        Property::new(name, PropertyKind::Metal)
    }

    /// New excitation.
    pub fn excitation(name: impl Into<String>, e: Excitation) -> Self {
        Property::new(name, PropertyKind::Excitation(e))
    }

    /// New probe.
    pub fn probe(name: impl Into<String>, p: ProbeBox) -> Self {
        Property::new(name, PropertyKind::ProbeBox(p))
    }

    /// New field dump.
    pub fn dump(name: impl Into<String>, d: DumpBox) -> Self {
        Property::new(name, PropertyKind::DumpBox(d))
    }

    /// New lumped element.
    pub fn lumped_element(name: impl Into<String>, l: LumpedElement) -> Self {
        Property::new(name, PropertyKind::LumpedElement(l))
    }

    /// New conducting sheet with conductivity (S/m) and thickness (m).
    pub fn conducting_sheet(
        name: impl Into<String>,
        conductivity: impl Into<Scalar>,
        thickness: impl Into<Scalar>,
    ) -> Self {
        Property::new(
            name,
            PropertyKind::ConductingSheet(ConductingSheet {
                conductivity: conductivity.into(),
                thickness: thickness.into(),
            }),
        )
    }

    /// New Drude/Lorentz dispersive material.
    pub fn lorentz_material(name: impl Into<String>, m: LorentzMaterial) -> Self {
        Property::new(name, PropertyKind::LorentzMaterial(m))
    }

    /// New Debye dispersive material.
    pub fn debye_material(name: impl Into<String>, m: DebyeMaterial) -> Self {
        Property::new(name, PropertyKind::DebyeMaterial(m))
    }

    /// Builder: sets fill and edge colour.
    #[must_use]
    pub fn with_color(mut self, c: Rgba) -> Self {
        self.fill_color = c;
        self.edge_color = c;
        self
    }

    /// Unique id within the owning structure (0 for detached properties).
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

    /// Adds a primitive.
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

    /// Sets a user attribute.
    pub fn set_attribute(&mut self, name: impl Into<String>, value: impl Into<String>) {
        let (name, value) = (name.into(), value.into());
        match self.attributes.iter_mut().find(|(k, _)| *k == name) {
            Some(e) => e.1 = value,
            None => self.attributes.push((name, value)),
        }
    }

    /// Re-evaluates parametric values of the property and its primitives.
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
