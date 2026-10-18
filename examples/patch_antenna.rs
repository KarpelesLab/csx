//! Builds a parametric patch antenna, meshes it, writes an openEMS file,
//! reads it back and prints a summary.
//!
//! ```text
//! cargo run --example patch_antenna -- [output.xml] [output.vtk]
//! ```

use csx::fdtd::FdtdSettings;
use csx::io::export::{VtkOptions, vtk_ascii};
use csx::mesh::{MeshBuilder, MeshOptions};
use csx::primitives::BoxShape;
use csx::properties::{DumpBox, Excitation, LumpedElement, Material, ProbeBox, Rgba};
use csx::{
    ContinuousStructure, CoordSystem, OpenEmsDocument, Primitive, Property, PropertyKind, Scalar,
    ScalarPoint, TypeMask, Vec3,
};

/// A box whose corners are given as expressions.
fn pbox(a: [&str; 3], b: [&str; 3]) -> BoxShape {
    let p = |v: [&str; 3]| ScalarPoint(v.map(Scalar::expr));
    BoxShape {
        start: p(a),
        stop: p(b),
    }
}

fn colored(mut p: Property, seed: u32) -> Property {
    p.fill_color = Rgba::from_seed(seed);
    p.edge_color = p.fill_color;
    p
}

fn build() -> csx::Result<ContinuousStructure> {
    let mut s = ContinuousStructure::new();
    // Dimensions in mm.
    for (name, v) in [
        ("patch_w", 32.86),
        ("patch_l", 41.37),
        ("sub_w", 60.0),
        ("sub_l", 60.0),
        ("sub_h", 1.524),
        ("feed_pos", -6.0),
        ("eps_r", 3.38),
    ] {
        s.parameters.set(name, v);
    }

    let mut sub = Material::default();
    sub.epsilon.set("eps_r");
    // Loss tangent 1e-3 at 2.45 GHz.
    sub.kappa
        .set(1e-3 * 2.0 * core::f64::consts::PI * 2.45e9 * 8.854e-12 * 3.38);
    let id = s.add_property(colored(Property::material("substrate", sub), 1));
    s.add_primitive(
        id,
        Primitive::new(pbox(
            ["-sub_w/2", "-sub_l/2", "0"],
            ["sub_w/2", "sub_l/2", "sub_h"],
        )),
    )?;

    let gnd = s.add_property(colored(Property::metal("ground"), 2));
    s.add_primitive(
        gnd,
        Primitive::new(pbox(
            ["-sub_w/2", "-sub_l/2", "0"],
            ["sub_w/2", "sub_l/2", "0"],
        ))
        .with_priority(10),
    )?;

    let patch = s.add_property(colored(Property::metal("patch"), 3));
    s.add_primitive(
        patch,
        Primitive::new(pbox(
            ["-patch_w/2", "-patch_l/2", "sub_h"],
            ["patch_w/2", "patch_l/2", "sub_h"],
        ))
        .with_priority(10),
    )?;

    // Lumped 50 Ω port between ground and patch.
    let port = ["feed_pos", "-0.5", "0"];
    let port_end = ["feed_pos+1", "0.5", "sub_h"];
    let le = LumpedElement {
        direction: csx::Axis::Z,
        r: Scalar::new(50.0),
        ..Default::default()
    };
    let r = s.add_property(Property::new(
        "port_resist",
        PropertyKind::LumpedElement(le),
    ));
    s.add_primitive(r, Primitive::new(pbox(port, port_end)).with_priority(5))?;

    let ex = Excitation {
        excite: [0.0.into(), 0.0.into(), (-1.0).into()],
        ..Default::default()
    };
    let e = s.add_property(Property::new("port_excite", PropertyKind::Excitation(ex)));
    s.add_primitive(e, Primitive::new(pbox(port, port_end)).with_priority(5))?;

    let ut = ProbeBox {
        probe_type: csx::properties::ProbeType::Voltage,
        weight: Scalar::new(-1.0),
        ..Default::default()
    };
    let u = s.add_property(Property::new("port_ut", PropertyKind::ProbeBox(ut)));
    s.add_primitive(
        u,
        Primitive::new(pbox(
            ["feed_pos+0.5", "0", "0"],
            ["feed_pos+0.5", "0", "sub_h"],
        )),
    )?;
    let it = ProbeBox {
        probe_type: csx::properties::ProbeType::Current,
        norm_dir: Some(csx::Axis::Z),
        ..Default::default()
    };
    let i = s.add_property(Property::new("port_it", PropertyKind::ProbeBox(it)));
    s.add_primitive(
        i,
        Primitive::new(pbox(
            ["feed_pos", "-0.5", "sub_h/2"],
            ["feed_pos+1", "0.5", "sub_h/2"],
        )),
    )?;

    let d = s.add_property(Property::new(
        "Et",
        PropertyKind::DumpBox(DumpBox::default()),
    ));
    s.add_primitive(
        d,
        Primitive::new(pbox(
            ["-sub_w/2", "-sub_l/2", "sub_h/2"],
            ["sub_w/2", "sub_l/2", "sub_h/2"],
        )),
    )?;

    s.update()?;
    s.grid.delta_unit = 1e-3;
    MeshBuilder::new(
        MeshOptions::new(2.0)
            .with_max_res([2.0, 2.0, 1.0])
            .with_ratio(1.4)
            .with_margin(30.0)
            .with_mask(TypeMask::PHYSICAL),
    )
    .with_lines(2, [0.508, 1.016])
    .apply(&mut s)?;
    Ok(s)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let out = args
        .next()
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("csx_patch_antenna.xml"));

    let structure = build()?;
    let mut doc = OpenEmsDocument {
        structure,
        ..Default::default()
    };
    doc.set_fdtd_settings(&FdtdSettings::gaussian(2.0e9, 1.0e9));
    doc.write_file(&out)?;
    println!("wrote {}", out.display());

    let back = OpenEmsDocument::read_file(&out)?;
    let s = &back.structure;
    println!(
        "grid: {:?} lines, {} cells, delta unit {}",
        s.grid.num_lines(),
        s.grid.num_cells(),
        s.grid.delta_unit
    );
    if let Some(f) = back.fdtd_settings() {
        println!("excitation: {:?}", f?.excitation);
    }
    for p in s.properties() {
        println!("{:>14}  {}", p.kind.type_name(), p.name);
    }
    let probe = Vec3::new(0.0, 0.0, 0.5);
    if let Some(m) = s.medium_at(probe, CoordSystem::Cartesian) {
        println!("medium at {:?}: eps = {:?}", probe.0, m.epsilon);
    }

    if let Some(vtk) = args.next() {
        std::fs::write(&vtk, vtk_ascii(s, &VtkOptions::default()))?;
        println!("wrote {vtk}");
    }
    Ok(())
}
