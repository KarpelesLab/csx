//! CSXCAD / openEMS XML files (feature `xml`).
//!
//! The reader accepts both geometry-only files (root `<ContinuousStructure>`)
//! and complete openEMS simulation files (root `<openEMS>` containing
//! `<FDTD>` and `<ContinuousStructure>`). The writer follows the layout of
//! the files written by CSXCAD itself: grid, background material, parameter
//! set and properties; every property carries an `ID`, a `Name`, fill and
//! edge colours and its primitives.
//!
//! Unknown attributes of properties, primitives and the structure are kept
//! in `extra_attributes`, and unknown child elements in `extra_elements`
//! (as [`dom::Element`] trees); both are written back. Unknown children are
//! re-emitted after the known ones.
//!
//! # Format notes
//!
//! The format was verified black-box against files written and read by
//! CSXCAD 0.7 itself (via its documented Python API):
//!
//! * Properties: `FillColor`, `EdgeColor`, `Attributes`, `Primitives`, then
//!   kind specific children (`Property`/`Weight` for materials, `Weight` for
//!   excitations, `FD_Samples` for probes, `DiscFile` for discrete
//!   materials).
//! * `<ParameterSet>` entries: `<Parameter Type="Const" name=".." Sweep="1"
//!   value=".."/>` and `Type="Linear"` with `min`, `max`, `step`.
//! * Transformations: `Translate`, `Rotate_X/Y/Z` (radians),
//!   `Rotate_Origin Argument="ax,ay,az,angle"`, `Scale` (one value),
//!   `Scale3` (three values), `Matrix` (16 row-major values acting on column
//!   vectors, translation in the last column).
//! * `<RotPoly>`: `NormDir`, `Elevation`, `RotAxisDir`, `<Vertex X1 X2/>` and
//!   `<Angles Start=".." Stop=".."/>` in radians.
//! * `<Polyhedron>`: `<Vertex>x,y,z</Vertex>` and `<Face>i,j,k</Face>`.
//! * `<MultiBox QtyBox="n">` with `<StartP/>`/`<EndP/>` pairs; `<Point X Y Z/>`.
//! * Dispersive materials: pole `n` (1-based) uses a `_n` suffix in both
//!   `Property` and `Weight` (`EpsilonPlasmaFrequency_1`,
//!   `MueLorPoleFrequency_2`, `EpsilonDelta_1`, `EpsilonRelaxTime_1`, ...);
//!   a name without suffix is read as the first pole.
//! * `<AbsorbingBC NormalSignPositive PhaseVelocity AbsorbingBoundaryType>`;
//!   `<DiscMaterial>` with `<DiscFile Type File UseDBBackground Scale/>`.
//! * Unset lumped element values are omitted.
//!
//! `<UserDefined>` is not understood by CSXCAD 0.7; this crate writes it
//! with an `<Equation>` text and a `<CoordShift X Y Z/>` origin.
//!
//! Property elements with an unknown name are read as
//! [`PropertyKind::Unknown`](crate::PropertyKind::Unknown) remembering the
//! element name, and written back under that name with their attributes
//! and child elements.

pub use super::dom;
mod fdtd;
mod read;
mod write;

use std::path::Path;

pub use dom::{Element, Node};

use crate::error::{Error, Result};
use crate::fdtd::FdtdSettings;
use crate::structure::ContinuousStructure;

/// A complete openEMS simulation description: the FDTD settings (kept as raw
/// XML, with typed access through [`OpenEmsDocument::fdtd_settings`]) and
/// the structure.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct OpenEmsDocument {
    /// The `<FDTD>` element, if present.
    pub fdtd: Option<Element>,
    /// The geometry and material model.
    pub structure: ContinuousStructure,
    /// Other top-level elements, preserved verbatim.
    pub extra: Vec<Element>,
}

impl OpenEmsDocument {
    /// Parses a document. A bare `<ContinuousStructure>` root yields a
    /// document without FDTD settings.
    pub fn from_xml_str(xml: &str) -> Result<Self> {
        let root = dom::parse(xml)?;
        match root.name.as_str() {
            "ContinuousStructure" => Ok(OpenEmsDocument {
                structure: load_structure(&root)?,
                ..Default::default()
            }),
            "openEMS" => {
                let mut doc = OpenEmsDocument::default();
                let mut found = false;
                for c in root.elements() {
                    match c.name.as_str() {
                        "FDTD" if doc.fdtd.is_none() => doc.fdtd = Some(c.clone()),
                        "ContinuousStructure" if !found => {
                            doc.structure = load_structure(c)?;
                            found = true;
                        }
                        _ => doc.extra.push(c.clone()),
                    }
                }
                Ok(doc)
            }
            other => Err(Error::Format(format!("unexpected root element <{other}>"))),
        }
    }

    /// Typed view of the `<FDTD>` element; `None` when the document has no
    /// FDTD settings.
    pub fn fdtd_settings(&self) -> Option<Result<FdtdSettings>> {
        self.fdtd.as_ref().map(FdtdSettings::from_xml_element)
    }

    /// Replaces the `<FDTD>` element by `settings`. Child elements other
    /// than `<Excitation>` and `<BoundaryCond>` of an existing element are
    /// kept.
    pub fn set_fdtd_settings(&mut self, settings: &FdtdSettings) {
        self.fdtd = Some(match &self.fdtd {
            Some(old) => settings.merge_into_element(old),
            None => settings.to_xml_element(),
        });
    }

    /// Serialises to an `<openEMS>` document.
    pub fn to_xml_string(&self) -> String {
        let mut root = Element::new("openEMS");
        if let Some(f) = &self.fdtd {
            root.push(f.clone());
        }
        root.push(write::structure(&self.structure));
        for e in &self.extra {
            root.push(e.clone());
        }
        root.to_document()
    }

    /// Reads a file; the structure's `base_dir` is set to the file's
    /// directory.
    pub fn read_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let mut doc = Self::from_xml_str(&read_text(path)?)?;
        doc.structure.base_dir = path.parent().map(Path::to_path_buf);
        Ok(doc)
    }

    /// Writes the document to a file.
    pub fn write_file(&self, path: impl AsRef<Path>) -> Result<()> {
        std::fs::write(path, self.to_xml_string())?;
        Ok(())
    }
}

fn read_text(path: &Path) -> Result<String> {
    let bytes = std::fs::read(path)?;
    // Legacy files may carry Latin-1 names; decode them byte by byte.
    Ok(String::from_utf8(bytes)
        .unwrap_or_else(|e| e.as_bytes().iter().map(|&c| c as char).collect()))
}

fn load_structure(el: &Element) -> Result<ContinuousStructure> {
    let mut s = read::structure(el)?;
    s.update_parameters()?;
    Ok(s)
}

impl ContinuousStructure {
    /// Parses a structure from XML text. The root may be
    /// `<ContinuousStructure>` or `<openEMS>` (whose structure is used).
    /// Parametric values are evaluated; external mesh files are not loaded
    /// until [`ContinuousStructure::update`] is called.
    pub fn from_xml_str(xml: &str) -> Result<Self> {
        Ok(OpenEmsDocument::from_xml_str(xml)?.structure)
    }

    /// Serialises the structure with a `<ContinuousStructure>` root.
    pub fn to_xml_string(&self) -> String {
        write::structure(self).to_document()
    }

    /// The `<ContinuousStructure>` element.
    pub fn to_xml_element(&self) -> Element {
        write::structure(self)
    }

    /// Builds a structure from a `<ContinuousStructure>` element.
    pub fn from_xml_element(el: &Element) -> Result<Self> {
        load_structure(el)
    }

    /// Reads a structure from a file; `base_dir` is set to the file's
    /// directory.
    pub fn read_file(path: impl AsRef<Path>) -> Result<Self> {
        Ok(OpenEmsDocument::read_file(path)?.structure)
    }

    /// Writes the structure to a file.
    pub fn write_file(&self, path: impl AsRef<Path>) -> Result<()> {
        std::fs::write(path, self.to_xml_string())?;
        Ok(())
    }
}

#[cfg(test)]
mod tests;
