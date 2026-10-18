//! `<FDTD>` element ↔ [`FdtdSettings`].

use super::dom::{Element, Node};
use super::read::Attrs;
use crate::error::{Error, Result};
use crate::expr::parse_number_lenient;
use crate::fdtd::{
    Boundary, BoundaryConditions, ExcitationSignal, FACE_NAMES, FdtdSettings, TimeStepMethod,
};
use crate::param::format_f64;

fn number(s: &str, what: &str) -> Result<f64> {
    parse_number_lenient(s)
        .ok_or_else(|| Error::Format(format!("FDTD {what}: {s:?} is not a number")))
}

fn opt_number(a: &Attrs<'_>, names: &[&str]) -> Result<Option<f64>> {
    names
        .iter()
        .find_map(|n| a.get(n))
        .map(|v| number(v, names[0]))
        .transpose()
}

fn opt_int(a: &Attrs<'_>, names: &[&str]) -> Result<Option<i64>> {
    match opt_number(a, names)? {
        Some(v) if v.is_finite() && v.abs() < 9.0e18 => Ok(Some(v.round() as i64)),
        Some(v) => Err(Error::Format(format!(
            "FDTD {}: {v} is not an integer",
            names[0]
        ))),
        None => Ok(None),
    }
}

fn to_i32(v: Option<i64>, what: &str) -> Result<Option<i32>> {
    v.map(|x| i32::try_from(x).map_err(|_| Error::Format(format!("FDTD {what}: {x} out of range"))))
        .transpose()
}

fn to_u32(v: Option<i64>, what: &str) -> Result<Option<u32>> {
    v.map(|x| u32::try_from(x).map_err(|_| Error::Format(format!("FDTD {what}: {x} out of range"))))
        .transpose()
}

impl FdtdSettings {
    /// Parses an `<FDTD>` element.
    pub fn from_xml_element(el: &Element) -> Result<Self> {
        let a = Attrs::new(el);
        let mut s = FdtdSettings {
            number_of_timesteps: opt_int(&a, &["NumberOfTimesteps"])?.map(|v| v.max(0) as u64),
            end_criteria: opt_number(&a, &["endCriteria", "EndCriteria"])?,
            f_max: opt_number(&a, &["f_max"])?,
            over_sampling: to_u32(opt_int(&a, &["OverSampling"])?, "OverSampling")?,
            time_step_method: to_i32(opt_int(&a, &["TimeStepMethod"])?, "TimeStepMethod")?
                .map(TimeStepMethod::from_id),
            max_time: opt_number(&a, &["MaxTime"])?,
            time_step: opt_number(&a, &["TimeStep"])?,
            time_step_factor: opt_number(&a, &["TimeStepFactor"])?,
            cell_constant_material: opt_number(&a, &["CellConstantMaterial"])?.map(|v| v != 0.0),
            coord_system: to_i32(opt_int(&a, &["CoordSystem"])?, "CoordSystem")?,
            multi_grid: a.get("MultiGrid").map(str::to_owned),
            ..Default::default()
        };
        s.extra_attributes = a.rest();
        if let Some(e) = el.child("Excitation") {
            let (sig, extra) = excitation(e, s.f_max)?;
            s.excitation = Some(sig);
            s.excitation_extra_attributes = extra;
        }
        if let Some(b) = el.child("BoundaryCond") {
            s.boundary = Some(boundary(b)?);
        }
        Ok(s)
    }

    /// Builds an `<FDTD>` element (with `<Excitation>` and `<BoundaryCond>`
    /// children when set).
    pub fn to_xml_element(&self) -> Element {
        let mut e = Element::new("FDTD");
        if let Some(n) = self.number_of_timesteps {
            e.set_attr("NumberOfTimesteps", n.to_string());
        }
        let mut num = |name: &str, v: Option<f64>| {
            if let Some(v) = v {
                e.set_attr(name, format_f64(v));
            }
        };
        num("endCriteria", self.end_criteria);
        num("f_max", self.f_max);
        num("MaxTime", self.max_time);
        num("TimeStep", self.time_step);
        num("TimeStepFactor", self.time_step_factor);
        if let Some(v) = self.over_sampling {
            e.set_attr("OverSampling", v.to_string());
        }
        if let Some(v) = self.time_step_method {
            e.set_attr("TimeStepMethod", v.id().to_string());
        }
        if let Some(v) = self.cell_constant_material {
            e.set_attr("CellConstantMaterial", if v { "1" } else { "0" });
        }
        if let Some(v) = self.coord_system {
            e.set_attr("CoordSystem", v.to_string());
        }
        if let Some(v) = &self.multi_grid {
            e.set_attr("MultiGrid", v.clone());
        }
        for (k, v) in &self.extra_attributes {
            e.set_attr(k.clone(), v.clone());
        }
        if let Some(sig) = &self.excitation {
            let mut x = Element::new("Excitation").with_attr("Type", sig.type_id().to_string());
            match sig {
                ExcitationSignal::Gaussian { f0, fc } => {
                    x.set_attr("f0", format_f64(*f0));
                    x.set_attr("fc", format_f64(*fc));
                }
                ExcitationSignal::Sinus { f0 } => x.set_attr("f0", format_f64(*f0)),
                ExcitationSignal::Dirac { .. } | ExcitationSignal::Step { .. } => {}
                ExcitationSignal::Custom { function, f0, .. } => {
                    x.set_attr("f0", format_f64(*f0));
                    x.set_attr("Function", function.clone());
                }
            }
            for (k, v) in &self.excitation_extra_attributes {
                x.set_attr(k.clone(), v.clone());
            }
            e.push(x);
        }
        if let Some(b) = &self.boundary {
            let mut x = Element::new("BoundaryCond");
            for (name, f) in FACE_NAMES.iter().zip(&b.faces) {
                x.set_attr(*name, f.to_text());
            }
            if let Some(v) = b.mur_phase_velocity {
                x.set_attr("MUR_PhaseVelocity", format_f64(v));
            }
            if let Some(g) = &b.pml_grading {
                x.set_attr("PML_Grading", g.clone());
            }
            for (k, v) in &b.extra_attributes {
                x.set_attr(k.clone(), v.clone());
            }
            e.push(x);
        }
        e
    }

    /// Like [`FdtdSettings::to_xml_element`] but keeps the unrelated child
    /// elements (anything but `<Excitation>` / `<BoundaryCond>`) of `base`.
    pub fn merge_into_element(&self, base: &Element) -> Element {
        let mut e = self.to_xml_element();
        for c in &base.children {
            if let Node::Element(c) = c
                && c.name != "Excitation"
                && c.name != "BoundaryCond"
            {
                e.push(c.clone());
            }
        }
        e
    }
}

fn excitation(
    el: &Element,
    f_max: Option<f64>,
) -> Result<(ExcitationSignal, Vec<(String, String)>)> {
    let a = Attrs::new(el);
    let ty = opt_int(&a, &["Type"])?.unwrap_or(0);
    let f0 = opt_number(&a, &["f0"])?.unwrap_or(0.0);
    let fc = opt_number(&a, &["fc"])?;
    let fmax = f_max.unwrap_or(0.0);
    let sig = match ty {
        0 => ExcitationSignal::Gaussian {
            f0,
            fc: fc.unwrap_or(0.0),
        },
        1 => ExcitationSignal::Sinus { f0 },
        2 => ExcitationSignal::Dirac { f_max: fmax },
        3 => ExcitationSignal::Step { f_max: fmax },
        10 => ExcitationSignal::Custom {
            function: a.get("Function").unwrap_or("0").to_owned(),
            f0,
            f_max: fmax,
        },
        other => {
            return Err(Error::Format(format!(
                "unknown excitation signal type {other}"
            )));
        }
    };
    // f0 / fc are part of the typed signal even when the type does not use
    // them; only keep them verbatim when they are not.
    let mut extra = a.rest();
    if let ExcitationSignal::Dirac { .. } | ExcitationSignal::Step { .. } = sig {
        for (k, v) in el.attrs.iter().filter(|(k, _)| k == "f0" || k == "fc") {
            extra.push((k.clone(), v.clone()));
        }
    } else if !matches!(sig, ExcitationSignal::Gaussian { .. })
        && let Some(v) = el.attr("fc")
    {
        extra.push(("fc".into(), v.into()));
    }
    Ok((sig, extra))
}

fn boundary(el: &Element) -> Result<BoundaryConditions> {
    let a = Attrs::new(el);
    let mut b = BoundaryConditions::default();
    for (i, name) in FACE_NAMES.iter().enumerate() {
        if let Some(v) = a.get(name) {
            b.faces[i] = Boundary::parse(v)
                .ok_or_else(|| Error::Format(format!("unknown boundary condition {name}={v:?}")))?;
        }
    }
    b.mur_phase_velocity = opt_number(&a, &["MUR_PhaseVelocity"])?;
    b.pml_grading = a.get("PML_Grading").map(str::to_owned);
    b.extra_attributes = a.rest();
    Ok(b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::xml::OpenEmsDocument;
    use crate::io::xml::dom::parse;

    const DOC: &str = r#"<openEMS>
  <FDTD NumberOfTimesteps="30000" endCriteria="1e-05" f_max="1.5e+09" OverSampling="4" Custom="x">
    <Excitation Type="0" f0="1e+09" fc="5e+08"/>
    <BoundaryCond xmin="PML_8" xmax="PML 12" ymin="1" ymax="MUR" zmin="PEC" zmax="3" MUR_PhaseVelocity="2e8"/>
    <Other A="1"/>
  </FDTD>
  <ContinuousStructure CoordSystem="0"><Properties/></ContinuousStructure>
</openEMS>"#;

    #[test]
    fn parse_settings() {
        let doc = OpenEmsDocument::from_xml_str(DOC).unwrap();
        let s = doc.fdtd_settings().unwrap().unwrap();
        assert_eq!(s.number_of_timesteps, Some(30000));
        assert_eq!(s.end_criteria, Some(1e-5));
        assert_eq!(s.over_sampling, Some(4));
        assert_eq!(
            s.extra_attributes,
            vec![("Custom".to_string(), "x".to_string())]
        );
        assert_eq!(
            s.excitation,
            Some(ExcitationSignal::Gaussian { f0: 1e9, fc: 5e8 })
        );
        let b = s.boundary.as_ref().unwrap();
        assert_eq!(
            b.faces,
            [
                Boundary::Pml(8),
                Boundary::Pml(12),
                Boundary::Pmc,
                Boundary::Mur,
                Boundary::Pec,
                Boundary::Pml(8)
            ]
        );
        assert_eq!(b.mur_phase_velocity, Some(2e8));
    }

    #[test]
    fn round_trip() {
        let mut doc = OpenEmsDocument::from_xml_str(DOC).unwrap();
        let s = doc.fdtd_settings().unwrap().unwrap();
        doc.set_fdtd_settings(&s);
        // unrelated children survive
        assert!(doc.fdtd.as_ref().unwrap().child("Other").is_some());
        let again = OpenEmsDocument::from_xml_str(&doc.to_xml_string()).unwrap();
        assert_eq!(again.fdtd_settings().unwrap().unwrap(), s);
        assert_eq!(again.to_xml_string(), doc.to_xml_string());
    }

    #[test]
    fn signal_variants() {
        for (xml, sig) in [
            (
                r#"<FDTD f_max="2e9"><Excitation Type="3"/></FDTD>"#,
                ExcitationSignal::Step { f_max: 2e9 },
            ),
            (
                r#"<FDTD f_max="2e9"><Excitation Type="2"/></FDTD>"#,
                ExcitationSignal::Dirac { f_max: 2e9 },
            ),
            (
                r#"<FDTD><Excitation Type="1" f0="5e8"/></FDTD>"#,
                ExcitationSignal::Sinus { f0: 5e8 },
            ),
            (
                r#"<FDTD f_max="3e9"><Excitation Type="10" f0="1e9" Function="sin(t)"/></FDTD>"#,
                ExcitationSignal::Custom {
                    function: "sin(t)".into(),
                    f0: 1e9,
                    f_max: 3e9,
                },
            ),
        ] {
            let s = FdtdSettings::from_xml_element(&parse(xml).unwrap()).unwrap();
            assert_eq!(s.excitation.as_ref(), Some(&sig), "{xml}");
            let back = FdtdSettings::from_xml_element(&s.to_xml_element()).unwrap();
            assert_eq!(back, s);
        }
        assert!(
            FdtdSettings::from_xml_element(
                &parse(r#"<FDTD><Excitation Type="7"/></FDTD>"#).unwrap()
            )
            .is_err()
        );
        assert!(
            FdtdSettings::from_xml_element(
                &parse(r#"<FDTD><BoundaryCond xmin="wall"/></FDTD>"#).unwrap()
            )
            .is_err()
        );
        // a fresh document without raw FDTD element
        let mut doc = OpenEmsDocument::default();
        assert!(doc.fdtd_settings().is_none());
        doc.set_fdtd_settings(&FdtdSettings::with_defaults(1e9, 5e8));
        let x = doc.to_xml_string();
        assert!(x.contains(r#"xmin="PML_8""#), "{x}");
        assert!(x.contains(r#"NumberOfTimesteps="1000000000""#), "{x}");
    }
}
