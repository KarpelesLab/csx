//! Data model → CSXCAD XML.

use super::dom::Element;
use crate::param::{ParameterKind, ParameterSet, Scalar, ScalarPoint, format_f64};
use crate::primitives::{MeshFileType, PolyVertices, Primitive, Shape};
use crate::properties::{
    Aniso, Material, MaterialQuantity, ProbeBox, Property, PropertyKind, Rgba, WeightFn,
};
use crate::structure::ContinuousStructure;
use crate::transform::{Transform, TransformOp};

fn join<I: IntoIterator<Item = String>>(it: I) -> String {
    it.into_iter().collect::<Vec<_>>().join(",")
}

fn s3(v: &[Scalar; 3]) -> String {
    join(v.iter().map(Scalar::to_text))
}

fn w3(v: &[WeightFn; 3]) -> String {
    join(v.iter().map(|w| w.source().to_owned()))
}

fn point(name: &str, p: &ScalarPoint) -> Element {
    Element::new(name)
        .with_attr("X", p.0[0].to_text())
        .with_attr("Y", p.0[1].to_text())
        .with_attr("Z", p.0[2].to_text())
}

fn color(name: &str, c: Rgba) -> Element {
    Element::new(name)
        .with_attr("R", c.r.to_string())
        .with_attr("G", c.g.to_string())
        .with_attr("B", c.b.to_string())
        .with_attr("a", c.a.to_string())
}

fn b(v: bool) -> &'static str {
    if v { "1" } else { "0" }
}

fn extend(el: &mut Element, extra: &[(String, String)]) {
    for (k, v) in extra {
        if el.attr(k).is_none() {
            el.attrs.push((k.clone(), v.clone()));
        }
    }
}

/// Builds the `<ContinuousStructure>` element.
pub(super) fn structure(s: &ContinuousStructure) -> Element {
    let mut el = Element::new("ContinuousStructure")
        .with_attr("CoordSystem", s.coord_system.id().to_string());
    extend(&mut el, &s.extra_attributes);

    let g = &s.grid;
    let mut ge = Element::new("RectilinearGrid")
        .with_attr("DeltaUnit", format_f64(g.delta_unit))
        .with_attr("CoordSystem", g.coord_system.id().to_string());
    for (dir, name) in ["XLines", "YLines", "ZLines"].into_iter().enumerate() {
        let lines = g.lines(dir);
        let mut l = Element::new(name).with_attr("Qty", lines.len().to_string());
        if !lines.is_empty() {
            l.push_text(join(lines.iter().map(|&v| format_f64(v))));
        }
        ge.push(l);
    }
    el.push(ge);

    let bg = &s.background;
    el.push(
        Element::new("BackgroundMaterial")
            .with_attr("Epsilon", bg.epsilon.to_text())
            .with_attr("Mue", bg.mue.to_text())
            .with_attr("Kappa", bg.kappa.to_text())
            .with_attr("Sigma", bg.sigma.to_text()),
    );
    el.push(parameters(&s.parameters));

    let mut props = Element::new("Properties");
    for (i, p) in s.properties().iter().enumerate() {
        props.push(property(p, i));
    }
    el.push(props);
    el
}

fn parameters(ps: &ParameterSet) -> Element {
    let mut el = Element::new("ParameterSet");
    for p in ps.iter() {
        let mut e = Element::new("Parameter");
        match &p.kind {
            ParameterKind::Const => e.set_attr("Type", "Const"),
            ParameterKind::Linear { .. } => e.set_attr("Type", "Linear"),
        }
        e.set_attr("name", p.name.clone());
        e.set_attr("value", format_f64(p.value));
        e.set_attr("Sweep", b(p.sweep));
        if let ParameterKind::Linear { start, stop, step } = p.kind {
            e.set_attr("Start", format_f64(start));
            e.set_attr("Stop", format_f64(stop));
            e.set_attr("Step", format_f64(step));
        }
        el.push(e);
    }
    el
}

fn property(p: &Property, index: usize) -> Element {
    let mut el = Element::new(match &p.kind {
        PropertyKind::Unknown {
            element: Some(e), ..
        } if super::dom::is_valid_name(e) => e.as_str(),
        PropertyKind::Unknown { .. } => "Unknown",
        k => k.type_name(),
    });
    el.set_attr("ID", index.to_string());
    el.set_attr("Name", p.name.clone());
    if !p.visible {
        el.set_attr("Visible", "0");
    }
    let mut children = vec![
        color("FillColor", p.fill_color),
        color("EdgeColor", p.edge_color),
    ];

    match &p.kind {
        PropertyKind::Unknown { property, element } => {
            // A foreign element keeps its own name; `Property` is only
            // needed when it carries information beyond that name.
            if element.as_deref() != Some(property.as_str()) {
                el.set_attr("Property", property.clone());
            }
        }
        PropertyKind::Material(m) => material(&mut el, &mut children, m, &[]),
        PropertyKind::LorentzMaterial(m) => {
            let names = [
                "EpsilonPlasmaFrequency",
                "EpsilonLorPoleFrequency",
                "EpsilonRelaxTime",
                "MuePlasmaFrequency",
                "MueLorPoleFrequency",
                "MueRelaxTime",
            ];
            let mut disp = Vec::new();
            for (n, pole) in m.poles.iter().enumerate() {
                let q = [
                    &pole.eps_plasma_frequency,
                    &pole.eps_lorentz_pole_frequency,
                    &pole.eps_relax_time,
                    &pole.mue_plasma_frequency,
                    &pole.mue_lorentz_pole_frequency,
                    &pole.mue_relax_time,
                ];
                for (name, a) in names.iter().zip(q) {
                    disp.push((format!("{name}_{}", n + 1), a));
                }
            }
            material(&mut el, &mut children, &m.material, &disp);
        }
        PropertyKind::DebyeMaterial(m) => {
            let mut disp = Vec::new();
            for (n, pole) in m.poles.iter().enumerate() {
                disp.push((format!("EpsilonDelta_{}", n + 1), &pole.eps_delta));
                disp.push((format!("EpsilonRelaxTime_{}", n + 1), &pole.eps_relax_time));
            }
            material(&mut el, &mut children, &m.material, &disp);
        }
        PropertyKind::DiscMaterial(d) => {
            material(&mut el, &mut children, &d.material, &[]);
            el.set_attr("File", d.file.clone());
            el.set_attr("Scale", d.scale.to_text());
            extend(&mut el, &d.attributes);
        }
        PropertyKind::Metal => {}
        PropertyKind::ConductingSheet(c) => {
            el.set_attr("Conductivity", c.conductivity.to_text());
            el.set_attr("Thickness", c.thickness.to_text());
        }
        PropertyKind::Excitation(e) => {
            el.set_attr("Number", e.number.to_string());
            el.set_attr("Enabled", b(e.enabled));
            el.set_attr("Frequency", e.frequency.to_text());
            el.set_attr("Delay", e.delay.to_text());
            el.set_attr("Type", e.excitation_type.to_string());
            el.set_attr("Excite", s3(&e.excite));
            el.set_attr("PropDir", s3(&e.prop_dir));
            children.push(
                Element::new("Weight")
                    .with_attr("X", e.weight[0].source())
                    .with_attr("Y", e.weight[1].source())
                    .with_attr("Z", e.weight[2].source()),
            );
        }
        PropertyKind::ProbeBox(pb) => probe(&mut el, &mut children, pb),
        PropertyKind::DumpBox(d) => {
            probe(&mut el, &mut children, &d.probe);
            el.set_attr("DumpType", d.dump_type.to_string());
            el.set_attr("DumpMode", d.dump_mode.to_string());
            el.set_attr("FileType", d.file_type.to_string());
            el.set_attr("MultiGridLevel", d.multi_grid_level.to_string());
            if let Some(s) = d.sub_sampling {
                el.set_attr("SubSampling", join(s.iter().map(u32::to_string)));
            }
            if let Some(r) = &d.opt_resolution {
                el.set_attr("OptResolution", s3(r));
            }
        }
        PropertyKind::LumpedElement(l) => {
            el.set_attr("Direction", l.direction.to_string());
            el.set_attr("Caps", b(l.caps));
            el.set_attr("R", l.r.to_text());
            el.set_attr("C", l.c.to_text());
            el.set_attr("L", l.l.to_text());
            if let Some(t) = l.le_type {
                el.set_attr("LEtype", t.to_string());
            }
        }
        PropertyKind::AbsorbingBC(a) => {
            el.set_attr("NormalSignPositive", b(a.normal_sign_positive));
            el.set_attr("PhaseVelocity", a.phase_velocity.to_text());
            el.set_attr("AbsorbingBoundaryType", a.abc_type.to_string());
        }
    }
    extend(&mut el, &p.extra_attributes);

    if !p.attributes.is_empty() {
        let mut a = Element::new("Attributes");
        a.attrs = p.attributes.clone();
        children.push(a);
    }
    let mut prims = Element::new("Primitives");
    for q in &p.primitives {
        prims.push(primitive(q));
    }
    children.push(prims);
    for c in children {
        el.push(c);
    }
    el
}

fn material(
    el: &mut Element,
    children: &mut Vec<Element>,
    m: &Material,
    disp: &[(String, &Aniso)],
) {
    el.set_attr("Isotropy", b(m.isotropic));
    let mut prop = Element::new("Property");
    let mut weight = Element::new("Weight");
    for q in MaterialQuantity::ALL {
        let a = m.quantity(q);
        prop.set_attr(q.name(), s3(&a.value));
        weight.set_attr(q.name(), w3(&a.weight));
    }
    prop.set_attr("Density", m.density.to_text());
    weight.set_attr("Density", m.density_weight.source());
    for (name, a) in disp {
        prop.set_attr(name.clone(), s3(&a.value));
        weight.set_attr(name.clone(), w3(&a.weight));
    }
    children.push(prop);
    children.push(weight);
}

fn probe(el: &mut Element, children: &mut Vec<Element>, p: &ProbeBox) {
    el.set_attr("Number", p.number.to_string());
    el.set_attr("Type", p.probe_type.to_string());
    el.set_attr("Weight", p.weight.to_text());
    el.set_attr("NormDir", p.norm_dir.unwrap_or(-1).to_string());
    el.set_attr("StartTime", p.start_time.to_text());
    el.set_attr("StopTime", p.stop_time.to_text());
    if !p.frequencies.is_empty() {
        let mut fd = Element::new("FD_Samples");
        fd.push_text(join(p.frequencies.iter().map(|&f| format_f64(f))));
        children.push(fd);
    }
}

fn vertices_2d(el: &mut Element, v: &PolyVertices) {
    el.set_attr("QtyVertices", v.len().to_string());
    for [a, c] in v.points() {
        el.push(
            Element::new("Vertex")
                .with_attr("X1", a.to_text())
                .with_attr("X2", c.to_text()),
        );
    }
}

fn primitive(p: &Primitive) -> Element {
    let mut el = Element::new(p.shape.type_name());
    el.set_attr("Priority", p.priority.to_string());
    if let Some(cs) = p.coord_system {
        el.set_attr("CoordSystem", cs.id().to_string());
    }
    if !p.transform.is_empty() {
        el.push(transform(&p.transform));
    }
    match &p.shape {
        Shape::Box(s) => {
            el.push(point("P1", &s.start));
            el.push(point("P2", &s.stop));
        }
        Shape::MultiBox(s) => {
            for (a, c) in &s.boxes {
                el.push(point("StartP", a));
                el.push(point("EndP", c));
            }
        }
        Shape::Sphere(s) => {
            el.set_attr("Radius", s.radius.to_text());
            el.push(point("Center", &s.center));
        }
        Shape::SphericalShell(s) => {
            el.set_attr("Radius", s.radius.to_text());
            el.set_attr("ShellWidth", s.shell_width.to_text());
            el.push(point("Center", &s.center));
        }
        Shape::Cylinder(s) => {
            el.set_attr("Radius", s.radius.to_text());
            el.push(point("P1", &s.start));
            el.push(point("P2", &s.stop));
        }
        Shape::CylindricalShell(s) => {
            el.set_attr("Radius", s.radius.to_text());
            el.set_attr("ShellWidth", s.shell_width.to_text());
            el.push(point("P1", &s.start));
            el.push(point("P2", &s.stop));
        }
        Shape::Polygon(s) => {
            el.set_attr("NormDir", s.norm_dir.to_string());
            el.set_attr("Elevation", s.elevation.to_text());
            vertices_2d(&mut el, &s.vertices);
        }
        Shape::LinPoly(s) => {
            el.set_attr("NormDir", s.norm_dir.to_string());
            el.set_attr("Elevation", s.elevation.to_text());
            el.set_attr("Length", s.length.to_text());
            vertices_2d(&mut el, &s.vertices);
        }
        Shape::RotPoly(s) => {
            el.set_attr("NormDir", s.norm_dir.to_string());
            el.set_attr("Elevation", s.elevation.to_text());
            el.set_attr("RotAxisDir", s.rot_axis_dir.to_string());
            vertices_2d(&mut el, &s.vertices);
            el.push(
                Element::new("Angles")
                    .with_attr("Start", s.start_angle.to_text())
                    .with_attr("Stop", s.stop_angle.to_text()),
            );
        }
        Shape::Curve(s) => {
            for v in &s.points {
                el.push(point("Vertex", v));
            }
        }
        Shape::Wire(s) => {
            el.set_attr("WireRadius", s.radius.to_text());
            for v in &s.points {
                el.push(point("Vertex", v));
            }
        }
        Shape::Point(s) => {
            el.set_attr("X", s.point.0[0].to_text());
            el.set_attr("Y", s.point.0[1].to_text());
            el.set_attr("Z", s.point.0[2].to_text());
        }
        Shape::Polyhedron(s) => {
            for v in s.vertices() {
                let mut e = Element::new("Vertex");
                e.push_text(join(v.0.iter().map(|&c| format_f64(c))));
                el.push(e);
            }
            for f in s.faces() {
                let mut e = Element::new("Face");
                e.push_text(join(f.iter().map(u32::to_string)));
                el.push(e);
            }
        }
        Shape::PolyhedronReader(s) => {
            el.set_attr("FileName", s.file_name.clone());
            el.set_attr(
                "FileType",
                match s.file_type {
                    MeshFileType::Stl => "STL",
                    MeshFileType::Ply => "PLY",
                },
            );
        }
        Shape::UserDefined(s) => {
            el.push(point("CoordShift", &s.base));
            let mut e = Element::new("Equation");
            e.push_text(s.function.clone());
            el.push(e);
        }
    }
    extend(&mut el, &p.extra_attributes);
    el
}

/// Converts a degree angle to its radian text form.
fn rad(a: &Scalar, deg: bool) -> String {
    if !deg {
        return a.to_text();
    }
    match a.expression() {
        Some(e) => format!("({e})*pi/180"),
        None => format_f64(a.value().to_radians()),
    }
}

fn transform(t: &Transform) -> Element {
    let deg = t.angles_in_degrees();
    let mut el = Element::new("Transformation");
    for op in t.ops() {
        let arg = match op {
            TransformOp::RotateX(a) | TransformOp::RotateY(a) | TransformOp::RotateZ(a) => {
                rad(a, deg)
            }
            TransformOp::RotateOrigin(v, a) => format!("{},{}", s3(v), rad(a, deg)),
            _ => join(op.args().into_iter().map(Scalar::to_text)),
        };
        el.push(Element::new(op.name()).with_attr("Argument", arg));
    }
    el
}
