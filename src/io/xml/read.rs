//! CSXCAD XML → data model.

use core::cell::Cell;

use super::dom::Element;
use crate::coords::CoordSystem;
use crate::error::{Error, Result};
use crate::expr::parse_number_lenient;
use crate::math::Vec3;
use crate::param::{Parameter, ParameterKind, ParameterSet, Scalar, ScalarPoint};
use crate::primitives::{
    BoxShape, Curve, Cylinder, CylindricalShell, LinPoly, MeshFileType, MultiBox, PointShape,
    PolyVertices, Polygon, Polyhedron, PolyhedronReader, Primitive, RotPoly, Shape, Sphere,
    SphericalShell, UserDefined, Wire,
};
use crate::properties::{
    AbsorbingBC, Aniso, ConductingSheet, DebyeMaterial, DebyePole, DiscMaterial, DumpBox,
    Excitation, LorentzMaterial, LorentzPole, LumpedElement, Material, MaterialQuantity, ProbeBox,
    Property, PropertyKind, Rgba, WeightFn,
};
use crate::properties::{
    Axis, DumpFileType, DumpMode, DumpType, ExcitationType, LumpedElementType, ProbeType,
};
use crate::structure::ContinuousStructure;
use crate::transform::{Transform, TransformOp};

/// Attribute accessor that remembers which attributes were consumed so the
/// rest can be preserved verbatim.
pub(super) struct Attrs<'a> {
    el: &'a Element,
    used: Vec<Cell<bool>>,
}

impl<'a> Attrs<'a> {
    pub(super) fn new(el: &'a Element) -> Self {
        Attrs {
            el,
            used: vec![Cell::new(false); el.attrs.len()],
        }
    }

    /// Value of `name` (first occurrence), marking it as consumed.
    pub(super) fn get(&self, name: &str) -> Option<&'a str> {
        let i = self.el.attrs.iter().position(|(k, _)| k == name)?;
        self.used[i].set(true);
        Some(self.el.attrs[i].1.as_str())
    }

    /// First of several alternative names.
    fn get_any(&self, names: &[&str]) -> Option<&'a str> {
        names.iter().find_map(|n| self.get(n))
    }

    fn scalar(&self, name: &str) -> Option<Scalar> {
        self.get(name).map(scalar)
    }

    fn int(&self, name: &str) -> Result<Option<i32>> {
        self.get(name).map(|v| int(v, name)).transpose()
    }

    fn bool(&self, name: &str) -> Option<bool> {
        self.get(name).map(boolean)
    }

    /// Attributes that were not consumed.
    pub(super) fn rest(&self) -> Vec<(String, String)> {
        self.el
            .attrs
            .iter()
            .zip(&self.used)
            .filter(|(_, u)| !u.get())
            .map(|(a, _)| a.clone())
            .collect()
    }
}

/// Child elements of `el` whose names are not in `known`, preserved
/// verbatim for round-tripping.
fn unknown_children(el: &Element, known: &[&str]) -> Vec<Element> {
    el.elements()
        .filter(|c| !known.contains(&c.name.as_str()))
        .cloned()
        .collect()
}

/// Splits a comma separated list at commas outside parentheses.
pub(super) fn split_list(s: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let mut depth = 0i32;
    let mut start = 0;
    for (i, c) in s.char_indices() {
        match c {
            '(' => depth += 1,
            ')' => depth -= 1,
            ',' if depth <= 0 => {
                out.push(s[start..i].trim());
                start = i + 1;
            }
            _ => {}
        }
    }
    let last = s[start..].trim();
    if !last.is_empty() || !out.is_empty() {
        out.push(last);
    }
    out
}

/// Strips the `term:` prefix used by old Octave front-ends to mark
/// expressions.
fn term(s: &str) -> &str {
    let t = s.trim();
    t.strip_prefix("term:").map_or(t, str::trim)
}

fn scalar(s: &str) -> Scalar {
    Scalar::expr(term(s))
}

fn weight(s: &str) -> WeightFn {
    WeightFn::new(term(s))
}

fn number(s: &str, what: &str) -> Result<f64> {
    parse_number_lenient(s).ok_or_else(|| Error::Format(format!("{what}: {s:?} is not a number")))
}

fn int(s: &str, what: &str) -> Result<i32> {
    let v = number(s, what)?;
    if v.is_finite() && v.abs() <= i32::MAX as f64 {
        Ok(v.round() as i32)
    } else {
        Err(Error::Format(format!("{what}: {s:?} is not an integer")))
    }
}

fn boolean(s: &str) -> bool {
    match s.trim().to_ascii_lowercase().as_str() {
        "true" | "yes" => true,
        "false" | "no" | "" => false,
        t => parse_number_lenient(t).is_some_and(|v| v != 0.0),
    }
}

/// Three scalars from a list: one value applies to every axis.
fn scalar3(s: &str) -> [Scalar; 3] {
    let v = split_list(s);
    match v.len() {
        0 => Default::default(),
        1 => {
            let a = scalar(v[0]);
            [a.clone(), a.clone(), a]
        }
        _ => core::array::from_fn(|i| scalar(v.get(i).copied().unwrap_or("0"))),
    }
}

fn weight3(s: &str) -> [WeightFn; 3] {
    let v = split_list(s);
    match v.len() {
        0 => Default::default(),
        1 => core::array::from_fn(|_| weight(v[0])),
        _ => core::array::from_fn(|i| weight(v.get(i).copied().unwrap_or("1"))),
    }
}

fn point(el: &Element) -> ScalarPoint {
    let c = |n: &str| el.attr(n).map(scalar).unwrap_or_default();
    ScalarPoint([c("X"), c("Y"), c("Z")])
}

fn child_point(el: &Element, name: &str) -> Result<ScalarPoint> {
    el.child(name)
        .map(point)
        .ok_or_else(|| Error::Format(format!("<{}> without <{name}>", el.name)))
}

fn color(el: Option<&Element>) -> Result<Rgba> {
    let Some(el) = el else {
        return Ok(Rgba::default());
    };
    let c = |n: &str, d: u8| -> Result<u8> {
        match el.attr(n) {
            Some(v) => Ok(int(v, n)?.clamp(0, 255) as u8),
            None => Ok(d),
        }
    };
    Ok(Rgba {
        r: c("R", 0)?,
        g: c("G", 0)?,
        b: c("B", 0)?,
        a: c("a", 255)?,
    })
}

fn coord_system(s: &str) -> Option<CoordSystem> {
    let v = parse_number_lenient(s)?;
    // `as` would map NaN to 0 (Cartesian); only accept exact ids.
    if v.fract() != 0.0 || !v.is_finite() {
        return None;
    }
    CoordSystem::from_id(v as i32)
}

/// Builds a structure from a `<ContinuousStructure>` element.
pub(super) fn structure(el: &Element) -> Result<ContinuousStructure> {
    let mut s = ContinuousStructure::new();
    let a = Attrs::new(el);
    if let Some(cs) = a.get("CoordSystem") {
        s.coord_system =
            coord_system(cs).ok_or_else(|| Error::Format(format!("CoordSystem {cs:?}")))?;
    }
    s.extra_attributes = a.rest();

    if let Some(ps) = el.child("ParameterSet") {
        s.parameters = parameters(ps)?;
    }
    if let Some(g) = el.child("RectilinearGrid") {
        grid(g, &mut s)?;
    }
    if let Some(b) = el.child("BackgroundMaterial") {
        let a = Attrs::new(b);
        let bg = &mut s.background;
        for (name, slot) in [
            ("Epsilon", &mut bg.epsilon),
            ("Mue", &mut bg.mue),
            ("Kappa", &mut bg.kappa),
            ("Sigma", &mut bg.sigma),
        ] {
            if let Some(v) = a.scalar(name) {
                *slot = v;
            }
        }
        bg.extra_attributes = a.rest();
    }
    if let Some(props) = el.child("Properties") {
        for p in props.elements() {
            s.add_property(property(p, &s.parameters)?);
        }
    }
    s.extra_elements = unknown_children(
        el,
        &[
            "ParameterSet",
            "RectilinearGrid",
            "BackgroundMaterial",
            "Properties",
        ],
    );
    Ok(s)
}

/// Parameter entries, as written by CSXCAD:
/// `<Parameter Type="Const" name=".." Sweep="1" value=".."/>` and
/// `<Parameter Type="Linear" ... min=".." max=".." step=".."/>`
/// (`Start`/`Stop`/`Step` are accepted too).
fn parameters(el: &Element) -> Result<ParameterSet> {
    let mut ps = ParameterSet::new();
    for p in el.elements() {
        let a = Attrs::new(p);
        let Some(name) = a.get_any(&["name", "Name"]) else {
            continue;
        };
        let value = match a.get_any(&["value", "Value"]) {
            Some(v) => ps.eval(v)?,
            None => 0.0,
        };
        let mut param = Parameter::new(name, value);
        param.sweep = a.bool("Sweep").unwrap_or(true);
        let ty = a.get_any(&["Type", "type"]).unwrap_or("Const");
        if ty.eq_ignore_ascii_case("linear") || p.name == "LinearParameter" {
            let f = |names: &[&str]| -> Result<f64> {
                a.get_any(names).map_or(Ok(0.0), |v| number(v, names[0]))
            };
            param.kind = ParameterKind::Linear {
                start: f(&["min", "Min", "Start"])?,
                stop: f(&["max", "Max", "Stop"])?,
                step: f(&["step", "Step"])?,
            };
        }
        ps.insert(param);
    }
    Ok(ps)
}

fn grid(el: &Element, s: &mut ContinuousStructure) -> Result<()> {
    let a = Attrs::new(el);
    let g = &mut s.grid;
    if let Some(du) = a.get("DeltaUnit") {
        g.delta_unit = s.parameters.eval(du)?;
    }
    if let Some(cs) = a.get("CoordSystem") {
        g.coord_system =
            coord_system(cs).ok_or_else(|| Error::Format(format!("CoordSystem {cs:?}")))?;
    }
    for (dir, name) in ["XLines", "YLines", "ZLines"].into_iter().enumerate() {
        if let Some(l) = el.child(name) {
            let text = l.text();
            let mut v = Vec::new();
            for t in split_list(&text) {
                if !t.is_empty() {
                    v.push(s.parameters.eval(t)?);
                }
            }
            g.set_lines(dir, v)?;
        }
    }
    g.extra_attributes = a.rest();
    Ok(())
}

fn property(el: &Element, ps: &ParameterSet) -> Result<Property> {
    let a = Attrs::new(el);
    a.get("ID");
    let name = a.get("Name").unwrap_or("").to_owned();
    let kind = match el.name.as_str() {
        "Material" => PropertyKind::Material(material(el, &a)?),
        "LorentzMaterial" => {
            let material = material(el, &a)?;
            PropertyKind::LorentzMaterial(LorentzMaterial {
                material,
                poles: lorentz_poles(el),
            })
        }
        "DebyeMaterial" => {
            let material = material(el, &a)?;
            PropertyKind::DebyeMaterial(DebyeMaterial {
                material,
                poles: debye_poles(el),
            })
        }
        "DiscMaterial" => {
            let material = material(el, &a)?;
            // CSXCAD stores the data file in a `<DiscFile Type File
            // UseDBBackground Scale/>` child; attributes on the property
            // itself are accepted as well.
            let mut file = a.get("File").unwrap_or("").to_owned();
            let mut scale = a.scalar("Scale").unwrap_or(Scalar::new(1.0));
            let mut attributes = Vec::new();
            if let Some(f) = el.child("DiscFile") {
                for (k, v) in &f.attrs {
                    match k.as_str() {
                        "File" => file = v.clone(),
                        "Scale" => scale = scalar(v),
                        // Default file type, written back implicitly.
                        "Type" if v.trim() == "0" => {}
                        _ => attributes.push((k.clone(), v.clone())),
                    }
                }
            }
            PropertyKind::DiscMaterial(DiscMaterial {
                material,
                file,
                scale,
                attributes,
            })
        }
        "Metal" => PropertyKind::Metal,
        "ConductingSheet" => PropertyKind::ConductingSheet(ConductingSheet {
            conductivity: a.scalar("Conductivity").unwrap_or_default(),
            thickness: a.scalar("Thickness").unwrap_or_default(),
        }),
        "Excitation" => {
            let mut e = Excitation {
                number: a.int("Number")?.unwrap_or(0).max(0) as u32,
                enabled: a.bool("Enabled").unwrap_or(true),
                excitation_type: ExcitationType::from_id(a.int("Type")?.unwrap_or(0)),
                ..Default::default()
            };
            if let Some(v) = a.get("Excite") {
                e.excite = scalar3(v);
            }
            if let Some(v) = a.get("PropDir") {
                e.prop_dir = scalar3(v);
            }
            if let Some(v) = a.scalar("Delay") {
                e.delay = v;
            }
            if let Some(v) = a.scalar("Frequency") {
                e.frequency = v;
            }
            if let Some(w) = el.child("Weight") {
                for (i, n) in ["X", "Y", "Z"].into_iter().enumerate() {
                    if let Some(f) = w.attr(n) {
                        e.weight[i] = weight(f);
                    }
                }
            }
            PropertyKind::Excitation(e)
        }
        "ProbeBox" => PropertyKind::ProbeBox(probe(el, &a)?),
        "DumpBox" => {
            let probe = probe(el, &a)?;
            let mut d = DumpBox {
                probe,
                ..Default::default()
            };
            d.dump_type = DumpType::from_id(a.int("DumpType")?.unwrap_or(0));
            d.dump_mode = DumpMode::from_id(a.int("DumpMode")?.unwrap_or(1));
            d.file_type = DumpFileType::from_id(a.int("FileType")?.unwrap_or(0));
            d.multi_grid_level = a.int("MultiGridLevel")?.unwrap_or(0);
            if let Some(v) = a.get("SubSampling") {
                let v = split_list(v)
                    .into_iter()
                    .map(|x| int(x, "SubSampling").map(|i| i.max(1) as u32))
                    .collect::<Result<Vec<_>>>()?;
                d.sub_sampling = match v.len() {
                    0 => None,
                    1 => Some([v[0]; 3]),
                    _ => Some(core::array::from_fn(|i| v.get(i).copied().unwrap_or(1))),
                };
            }
            if let Some(v) = a.get("OptResolution") {
                d.opt_resolution = Some(scalar3(v));
            }
            PropertyKind::DumpBox(d)
        }
        "LumpedElement" => PropertyKind::LumpedElement(LumpedElement {
            direction: axis_attr(a.int("Direction")?.unwrap_or(0), "Direction")?,
            caps: a.bool("Caps").unwrap_or(true),
            r: a.scalar("R").unwrap_or(Scalar::new(f64::NAN)),
            c: a.scalar("C").unwrap_or(Scalar::new(f64::NAN)),
            l: a.scalar("L").unwrap_or(Scalar::new(f64::NAN)),
            le_type: a.int("LEtype")?.map(LumpedElementType::from_id),
        }),
        "AbsorbingBC" => PropertyKind::AbsorbingBC(AbsorbingBC {
            normal_sign_positive: a.bool("NormalSignPositive").unwrap_or(true),
            phase_velocity: a.scalar("PhaseVelocity").unwrap_or_default(),
            abc_type: a
                .get("AbsorbingBoundaryType")
                .map(str::parse)
                .transpose()?
                .unwrap_or_default(),
        }),
        "Unknown" => PropertyKind::Unknown {
            property: a.get("Property").unwrap_or("").to_owned(),
            element: None,
        },
        other => PropertyKind::Unknown {
            property: a.get("Property").unwrap_or(other).to_owned(),
            element: Some(other.to_owned()),
        },
    };
    let mut p = Property::new(name, kind);
    p.visible = a.bool("Visible").unwrap_or(true);
    p.fill_color = color(el.child("FillColor"))?;
    p.edge_color = color(el.child("EdgeColor"))?;
    if let Some(attrs) = el.child("Attributes") {
        p.attributes = attrs.attrs.clone();
    }
    if let Some(prims) = el.child("Primitives") {
        for q in prims.elements() {
            p.primitives.push(primitive(q, ps)?);
        }
    }
    let specific: &[&str] = match el.name.as_str() {
        "Material" | "LorentzMaterial" | "DebyeMaterial" => &["Property", "Weight"],
        "DiscMaterial" => &["Property", "Weight", "DiscFile"],
        "Excitation" => &["Weight"],
        "ProbeBox" | "DumpBox" => &["FD_Samples"],
        _ => &[],
    };
    let mut known = vec!["FillColor", "EdgeColor", "Attributes", "Primitives"];
    known.extend_from_slice(specific);
    p.extra_elements = unknown_children(el, &known);
    p.extra_attributes = a.rest();
    Ok(p)
}

fn material(el: &Element, a: &Attrs<'_>) -> Result<Material> {
    let mut m = Material {
        isotropic: a.bool("Isotropy").unwrap_or(true),
        ..Default::default()
    };
    if let Some(p) = el.child("Property") {
        for q in MaterialQuantity::ALL {
            if let Some(v) = p.attr(q.name()) {
                m.quantity_mut(q).value = scalar3(v);
            }
        }
        if let Some(v) = p.attr("Density") {
            m.density = scalar(v);
        }
    }
    if let Some(w) = el.child("Weight") {
        for q in MaterialQuantity::ALL {
            if let Some(v) = w.attr(q.name()) {
                m.quantity_mut(q).weight = weight3(v);
            }
        }
        if let Some(v) = w.attr("Density") {
            m.density_weight = weight(v);
        }
    }
    Ok(m)
}

/// Order (1-based) encoded in an attribute name `base` / `base_<n>`.
pub(super) fn order_of(name: &str, base: &str) -> Option<usize> {
    let rest = name.strip_prefix(base)?;
    if rest.is_empty() {
        return Some(1);
    }
    rest.strip_prefix('_')?
        .parse::<usize>()
        .ok()
        .filter(|&n| (1..=1024).contains(&n))
}

fn max_order(el: &Element, bases: &[&str]) -> usize {
    ["Property", "Weight"]
        .iter()
        .filter_map(|c| el.child(c))
        .flat_map(|c| c.attrs.iter())
        .filter_map(|(k, _)| bases.iter().find_map(|b| order_of(k, b)))
        .max()
        .unwrap_or(0)
}

/// Reads the dispersive quantity `base` of order `n` from `<Property>` and
/// `<Weight>` into `q` (values default to 0, weights to 1).
fn dispersive(el: &Element, base: &str, n: usize, q: &mut Aniso) {
    let find = |c: &Element| -> Option<String> {
        let suffixed = format!("{base}_{n}");
        c.attr(&suffixed)
            .or_else(|| if n == 1 { c.attr(base) } else { None })
            .map(str::to_owned)
    };
    if let Some(v) = el.child("Property").and_then(find) {
        q.value = scalar3(&v);
    }
    if let Some(v) = el.child("Weight").and_then(find) {
        q.weight = weight3(&v);
    }
}

const LORENTZ_NAMES: [&str; 6] = [
    "EpsilonPlasmaFrequency",
    "EpsilonLorPoleFrequency",
    "EpsilonRelaxTime",
    "MuePlasmaFrequency",
    "MueLorPoleFrequency",
    "MueRelaxTime",
];

const DEBYE_NAMES: [&str; 2] = ["EpsilonDelta", "EpsilonRelaxTime"];

fn lorentz_poles(el: &Element) -> Vec<LorentzPole> {
    (1..=max_order(el, &LORENTZ_NAMES))
        .map(|n| {
            let mut p = LorentzPole::zero();
            for (base, q) in LORENTZ_NAMES.iter().zip([
                &mut p.eps_plasma_frequency,
                &mut p.eps_lorentz_pole_frequency,
                &mut p.eps_relax_time,
                &mut p.mue_plasma_frequency,
                &mut p.mue_lorentz_pole_frequency,
                &mut p.mue_relax_time,
            ]) {
                dispersive(el, base, n, q);
            }
            p
        })
        .collect()
}

fn debye_poles(el: &Element) -> Vec<DebyePole> {
    (1..=max_order(el, &DEBYE_NAMES))
        .map(|n| {
            let mut p = DebyePole::default();
            dispersive(el, DEBYE_NAMES[0], n, &mut p.eps_delta);
            dispersive(el, DEBYE_NAMES[1], n, &mut p.eps_relax_time);
            p
        })
        .collect()
}

/// Converts an axis index read from attribute `what`.
fn axis_attr(d: i32, what: &str) -> Result<Axis> {
    Axis::try_from(d).map_err(|_| Error::Format(format!("{what}: {d} is not an axis (0, 1 or 2)")))
}

fn probe(el: &Element, a: &Attrs<'_>) -> Result<ProbeBox> {
    let mut p = ProbeBox {
        number: a.int("Number")?.unwrap_or(0).max(0) as u32,
        probe_type: ProbeType::from_id(a.int("Type")?.unwrap_or(0)),
        norm_dir: match a.int("NormDir")? {
            Some(d) if d >= 0 => Some(axis_attr(d, "NormDir")?),
            _ => None,
        },
        ..Default::default()
    };
    if let Some(v) = a.scalar("Weight") {
        p.weight = v;
    }
    if let Some(v) = a.scalar("StartTime") {
        p.start_time = v;
    }
    if let Some(v) = a.scalar("StopTime") {
        p.stop_time = v;
    }
    if let Some(fd) = el.child("FD_Samples") {
        let text = fd.text();
        for t in text
            .split([',', ';', ' ', '\t', '\n'])
            .filter(|t| !t.trim().is_empty())
        {
            p.frequencies.push(number(t, "FD_Samples")?);
        }
    }
    Ok(p)
}

/// A coordinate that must be numeric in the model (polyhedron vertices): plain numbers, or expressions evaluated once against the
/// file's parameters.
fn coordinate(s: &str, ps: &ParameterSet, what: &str) -> Result<f64> {
    match parse_number_lenient(s) {
        Some(v) => Ok(v),
        None => ps
            .eval(term(s))
            .map_err(|e| Error::Format(format!("{what}: {s:?}: {e}"))),
    }
}

/// Polygon vertices; expressions are kept and evaluated by the structure
/// update that follows reading.
fn vertices_2d(el: &Element) -> PolyVertices {
    el.children_named("Vertex")
        .map(|v| {
            let c = |n: &str| v.attr(n).map_or(Scalar::new(0.0), scalar);
            [c("X1"), c("X2")]
        })
        .collect()
}

fn axis(a: &Attrs<'_>, names: &[&str], default: usize) -> Result<usize> {
    match a.get_any(names) {
        Some(v) => {
            let d = int(v, names[0])?;
            usize::try_from(d)
                .ok()
                .filter(|&d| d < 3)
                .ok_or_else(|| Error::Format(format!("{} {d} is not an axis", names[0])))
        }
        None => Ok(default),
    }
}

fn polyhedron(el: &Element, ps: &ParameterSet) -> Result<Polyhedron> {
    let mut verts = Vec::new();
    for v in el.children_named("Vertex") {
        let p = if v.attr("X").is_some() {
            let c = |n: &str| v.attr(n).map_or(Ok(0.0), |s| coordinate(s, ps, n));
            Vec3::new(c("X")?, c("Y")?, c("Z")?)
        } else {
            let t = v.text();
            let c: Vec<f64> = split_list(&t)
                .into_iter()
                .map(|s| coordinate(s, ps, "Vertex"))
                .collect::<Result<_>>()?;
            if c.len() != 3 {
                return Err(Error::Format(format!("polyhedron vertex {t:?}")));
            }
            Vec3::new(c[0], c[1], c[2])
        };
        verts.push(p);
    }
    let mut faces = Vec::new();
    for f in el.children_named("Face") {
        let t = f.text();
        let idx = split_list(&t)
            .into_iter()
            .map(|s| {
                let v = number(s, "Face")?;
                if v >= 0.0 && v <= u32::MAX as f64 && v.fract() == 0.0 {
                    Ok(v as u32)
                } else {
                    Err(Error::Format(format!("polyhedron face index {s:?}")))
                }
            })
            .collect::<Result<Vec<u32>>>()?;
        faces.push(idx);
    }
    Polyhedron::new(verts, faces)
}

fn primitive(el: &Element, ps: &ParameterSet) -> Result<Primitive> {
    let a = Attrs::new(el);
    // Child element a `<Point>` took its coordinates from, if any.
    let mut point_child: Option<&str> = None;
    let need = |name: &str| -> Result<Scalar> {
        a.scalar(name)
            .ok_or_else(|| Error::Format(format!("<{}> without {name}", el.name)))
    };
    let shape = match el.name.as_str() {
        "Box" => Shape::Box(BoxShape {
            start: child_point(el, "P1")?,
            stop: child_point(el, "P2")?,
        }),
        "MultiBox" => {
            a.get("QtyBox");
            let starts: Vec<_> = el.children_named("StartP").map(point).collect();
            let ends: Vec<_> = el.children_named("EndP").map(point).collect();
            if starts.len() != ends.len() {
                return Err(Error::Format("MultiBox with unpaired StartP/EndP".into()));
            }
            Shape::MultiBox(MultiBox {
                boxes: starts.into_iter().zip(ends).collect(),
            })
        }
        "Sphere" => Shape::Sphere(Sphere {
            center: child_point(el, "Center")?,
            radius: need("Radius")?,
        }),
        "SphericalShell" => Shape::SphericalShell(SphericalShell {
            center: child_point(el, "Center")?,
            radius: need("Radius")?,
            shell_width: need("ShellWidth")?,
        }),
        "Cylinder" => Shape::Cylinder(Cylinder {
            start: child_point(el, "P1")?,
            stop: child_point(el, "P2")?,
            radius: need("Radius")?,
        }),
        "CylindricalShell" => Shape::CylindricalShell(CylindricalShell {
            start: child_point(el, "P1")?,
            stop: child_point(el, "P2")?,
            radius: need("Radius")?,
            shell_width: need("ShellWidth")?,
        }),
        "Polygon" => {
            a.get("QtyVertices");
            Shape::Polygon(Polygon {
                norm_dir: axis(&a, &["NormDir"], 2)?,
                elevation: a.scalar("Elevation").unwrap_or_default(),
                vertices: vertices_2d(el),
            })
        }
        "LinPoly" => {
            a.get("QtyVertices");
            Shape::LinPoly(LinPoly {
                norm_dir: axis(&a, &["NormDir"], 2)?,
                elevation: a.scalar("Elevation").unwrap_or_default(),
                length: a.scalar("Length").unwrap_or_default(),
                vertices: vertices_2d(el),
            })
        }
        "RotPoly" => {
            a.get("QtyVertices");
            let mut r = RotPoly {
                norm_dir: axis(&a, &["NormDir"], 0)?,
                rot_axis_dir: axis(&a, &["RotAxisDir", "RotAxis", "RotationAxis"], 2)?,
                elevation: a.scalar("Elevation").unwrap_or_default(),
                vertices: vertices_2d(el),
                ..Default::default()
            };
            if let Some(ang) = el.child("Angles") {
                if let Some(v) = ang.attr("Start") {
                    r.start_angle = scalar(v);
                }
                if let Some(v) = ang.attr("Stop") {
                    r.stop_angle = scalar(v);
                }
            }
            if let Some(v) = a.scalar("StartAngle") {
                r.start_angle = v;
            }
            if let Some(v) = a.scalar("StopAngle") {
                r.stop_angle = v;
            }
            Shape::RotPoly(r)
        }
        "Curve" => Shape::Curve(Curve {
            points: el.children_named("Vertex").map(point).collect(),
        }),
        "Wire" => Shape::Wire(Wire {
            points: el.children_named("Vertex").map(point).collect(),
            radius: need("WireRadius")?,
        }),
        "Point" => {
            let p = if a.get("X").is_some() || a.get("Y").is_some() || a.get("Z").is_some() {
                point(el)
            } else if let Some(c) = el.child("Coordinates").or_else(|| el.child("P")) {
                point_child = Some(c.name.as_str());
                point(c)
            } else {
                ScalarPoint::default()
            };
            Shape::Point(PointShape { point: p })
        }
        "Polyhedron" => Shape::Polyhedron(polyhedron(el, ps)?),
        "PolyhedronReader" => {
            let file_name = a.get("FileName").unwrap_or("").to_owned();
            let file_type = match a.get("FileType") {
                Some(t) => t.parse::<MeshFileType>()?,
                None => MeshFileType::from_path(&file_name).unwrap_or_default(),
            };
            Shape::PolyhedronReader(PolyhedronReader::new(file_name).with_file_type(file_type))
        }
        "UserDefined" => {
            let f = el
                .child("Equation")
                .or_else(|| el.child("Function"))
                .map(Element::text)
                .or_else(|| a.get_any(&["Equation", "Function"]).map(str::to_owned))
                .unwrap_or_default();
            let mut u = UserDefined::new(f.trim());
            if let Some(b) = el.child("CoordShift").or_else(|| el.child("P0")) {
                u.base = point(b);
            }
            Shape::UserDefined(u)
        }
        other => return Err(Error::Format(format!("unknown primitive <{other}>"))),
    };
    let mut prim = Primitive::new(shape);
    prim.priority = a.int("Priority")?.unwrap_or(0);
    if let Some(cs) = el.attr("CoordSystem")
        && let Some(cs) = coord_system(cs)
    {
        a.get("CoordSystem");
        prim.coord_system = Some(cs);
    }
    if let Some(t) = el.child("Transformation") {
        prim.transform = transform(t)?;
    }
    let specific: &[&str] = match el.name.as_str() {
        "Box" | "Cylinder" | "CylindricalShell" => &["P1", "P2"],
        "MultiBox" => &["StartP", "EndP"],
        "Sphere" | "SphericalShell" => &["Center"],
        "Polygon" | "LinPoly" | "Curve" | "Wire" => &["Vertex"],
        "RotPoly" => &["Vertex", "Angles"],
        "Polyhedron" => &["Vertex", "Face"],
        "UserDefined" => &["Equation", "Function", "CoordShift", "P0"],
        _ => &[],
    };
    let mut known = vec!["Transformation"];
    known.extend_from_slice(specific);
    known.extend(point_child);
    prim.extra_elements = unknown_children(el, &known);
    prim.extra_attributes = a.rest();
    Ok(prim)
}

fn transform(el: &Element) -> Result<Transform> {
    let mut t = Transform::new();
    for op in el.elements() {
        let arg = op.attr("Argument").unwrap_or("");
        let args: Vec<Scalar> = split_list(arg).into_iter().map(scalar).collect();
        t.push(TransformOp::from_name_args(&op.name, args)?);
    }
    Ok(t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lists() {
        assert_eq!(split_list("1, pow(x,2) ,3"), vec!["1", "pow(x,2)", "3"]);
        assert_eq!(split_list(""), Vec::<&str>::new());
        assert_eq!(split_list("a"), vec!["a"]);
        assert_eq!(order_of("EpsilonDelta_3", "EpsilonDelta"), Some(3));
        assert_eq!(order_of("EpsilonDelta", "EpsilonDelta"), Some(1));
        assert_eq!(order_of("EpsilonDeltaX", "EpsilonDelta"), None);
        assert_eq!(order_of("EpsilonDelta_0", "EpsilonDelta"), None);
        assert!(boolean("1") && !boolean("0") && boolean("true") && !boolean("0.000000e+00"));
        assert_eq!(int("0.000000e+00", "x").unwrap(), 0);
        assert_eq!(coord_system("nan"), None);
        assert_eq!(coord_system("1.5"), None);
        assert_eq!(coord_system("1"), Some(CoordSystem::Cylindrical));
    }
}
