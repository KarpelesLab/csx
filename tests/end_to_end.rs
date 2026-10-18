//! End-to-end tests: build → update → mesh → discretise → XML round-trip →
//! tessellate / export.

use csx::discretize::{DiscretizeOptions, discretize};
use csx::fdtd::{Boundary, ExcitationSignal, FdtdSettings};
use csx::io::export::{VtkOptions, export_vtk, ply_ascii, stl_binary};
use csx::mesh::{MeshBuilder, MeshOptions, detect_edges};
use csx::primitives::{BoxShape, Cylinder, Sphere, Wire};
use csx::properties::{Excitation, LumpedElement, Material, MaterialQuantity, ProbeBox};
use csx::tessellate::{TessellationOptions, tessellate_structure};
use csx::{
    ContinuousStructure, CoordSystem, OpenEmsDocument, Primitive, Property, PropertyKind,
    PropertyType, Scalar, ScalarPoint, Shape, TypeMask, Vec3,
};

fn sp(v: [&str; 3]) -> ScalarPoint {
    ScalarPoint(v.map(Scalar::expr))
}

/// A parametric microstrip line: substrate, ground, strip, port, probe.
fn microstrip() -> ContinuousStructure {
    let mut s = ContinuousStructure::new();
    s.parameters.set("w", 3.0);
    s.parameters.set("h", 1.0);
    s.parameters.set("eps", 4.0);
    s.grid.delta_unit = 1e-3;

    let mut m = Material::default();
    m.epsilon.set("eps");
    let sub = s.add_property(Property::material("sub", m));
    s.add_primitive(
        sub,
        BoxShape {
            start: sp(["-10", "-10", "0"]),
            stop: sp(["10", "10", "h"]),
        },
    )
    .unwrap();
    let gnd = s.add_property(Property::metal("gnd"));
    s.add_primitive(
        gnd,
        Primitive::new(BoxShape::new([-10.0, -10.0, 0.0], [10.0, 10.0, 0.0])).with_priority(10),
    )
    .unwrap();
    let strip = s.add_property(Property::metal("strip"));
    s.add_primitive(
        strip,
        Primitive::new(BoxShape {
            start: sp(["-10", "-w/2", "h"]),
            stop: sp(["10", "w/2", "h"]),
        })
        .with_priority(10),
    )
    .unwrap();
    let ex = Excitation {
        excite: [0.0.into(), 0.0.into(), 1.0.into()],
        ..Default::default()
    };
    let e = s.add_property(Property::new("port", PropertyKind::Excitation(ex)));
    s.add_primitive(
        e,
        BoxShape {
            start: sp(["-10", "-w/2", "0"]),
            stop: sp(["-10", "w/2", "h"]),
        },
    )
    .unwrap();
    let le = LumpedElement {
        direction: csx::Axis::Z,
        r: 50.0.into(),
        ..Default::default()
    };
    let r = s.add_property(Property::new("load", PropertyKind::LumpedElement(le)));
    s.add_primitive(
        r,
        BoxShape {
            start: sp(["10", "-w/2", "0"]),
            stop: sp(["10", "w/2", "h"]),
        },
    )
    .unwrap();
    let p = s.add_property(Property::new(
        "ut",
        PropertyKind::ProbeBox(ProbeBox::default()),
    ));
    s.add_primitive(p, BoxShape::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]))
        .unwrap();
    s.update().unwrap();
    s
}

fn mesh(s: &mut ContinuousStructure) {
    MeshBuilder::new(
        MeshOptions::new(1.0)
            .with_max_res([1.0, 0.5, 0.25])
            .with_mask(TypeMask::PHYSICAL),
    )
    .apply(s)
    .unwrap();
}

/// Compares two structures property by property, tolerating NaN scalars
/// (unset lumped element values) by comparing the XML text instead.
fn assert_same(a: &ContinuousStructure, b: &ContinuousStructure) {
    assert_eq!(a.properties().len(), b.properties().len());
    assert_eq!(a.grid, b.grid);
    assert_eq!(a.parameters, b.parameters);
    assert_eq!(a.coord_system, b.coord_system);
    for (p, q) in a.properties().iter().zip(b.properties()) {
        if matches!(p.kind, PropertyKind::LumpedElement(_)) {
            assert_eq!(p.name, q.name);
            assert_eq!(p.primitives, q.primitives);
        } else {
            assert_eq!(p, q, "{}", p.name);
        }
    }
    assert_eq!(a.to_xml_string(), b.to_xml_string());
}

#[test]
fn full_pipeline() {
    let mut s = microstrip();
    let edges = detect_edges(&s, TypeMask::PHYSICAL, CoordSystem::Cartesian);
    assert!(edges[1].contains(&-1.5) && edges[1].contains(&1.5));
    assert!(edges[2].contains(&0.0) && edges[2].contains(&1.0));
    mesh(&mut s);
    let g = &s.grid;
    assert!(g.snap(2, 1.0).is_some() && g.snap(1, 1.5).is_some());
    for d in 0..3 {
        let (_, max) = g.cell_size_range(d).unwrap();
        assert!(max <= [1.0, 0.5, 0.25][d] * (1.0 + 1e-9), "axis {d}: {max}");
    }

    // Discretisation.
    let mg = discretize(&s, &DiscretizeOptions::default()).unwrap();
    assert!(mg.pec_count() > 0);
    let eps_max = mg.epsilon[0].iter().copied().fold(0.0, f64::max);
    assert!((eps_max - 4.0).abs() < 1e-12, "{eps_max}");
    assert_eq!(mg.excitations.len(), 1);
    assert!(!mg.excitations[0].amplitudes.is_empty());
    assert!(!mg.lumped.is_empty());
    assert!(mg.lumped.iter().all(|l| l.edge.dir == 2));

    // XML round-trip with FDTD settings.
    let mut doc = OpenEmsDocument {
        structure: s.clone(),
        ..Default::default()
    };
    let fdtd = FdtdSettings::with_defaults(5e9, 4e9);
    doc.set_fdtd_settings(&fdtd);
    let xml = doc.to_xml_string();
    let back = OpenEmsDocument::from_xml_str(&xml).unwrap();
    assert_eq!(back.fdtd_settings().unwrap().unwrap(), fdtd);
    assert_same(&s, &back.structure);
    assert_eq!(back.to_xml_string(), xml);

    // Re-discretising the re-read structure yields the same result.
    let mg2 = discretize(&back.structure, &DiscretizeOptions::default()).unwrap();
    assert_eq!(mg.epsilon, mg2.epsilon);
    assert_eq!(mg.pec, mg2.pec);

    // Tessellation and exports.
    let parts = tessellate_structure(&s, TypeMask::ANY, &TessellationOptions::default());
    assert_eq!(parts.len(), 6);
    let sub = &parts[0].mesh;
    assert!((sub.signed_volume() - 400.0).abs() < 1e-9);
    let stl = stl_binary(sub);
    let reread = csx::io::stl::parse(&stl).unwrap();
    assert!(reread.is_closed());
    assert!(ply_ascii(sub).starts_with("ply"));
    let vtk = export_vtk(&s, &VtkOptions::default());
    assert!(vtk.contains("POLYDATA") && vtk.contains("property_index"));
}

#[test]
fn file_roundtrip() {
    let s = microstrip();
    let dir = std::env::temp_dir().join(format!("csx-e2e-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("model.xml");
    s.write_file(&path).unwrap();
    let r = ContinuousStructure::read_file(&path).unwrap();
    assert_eq!(r.base_dir.as_deref(), Some(dir.as_path()));
    assert_same(&s, &r);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn parametric_update_changes_results() {
    let mut s = microstrip();
    let p = Vec3::new(0.0, 2.0, 1.0);
    let hit = |s: &ContinuousStructure| {
        s.property_at(p, CoordSystem::Cartesian, TypeMask::of(PropertyType::Metal))
            .map(|h| h.property.name.clone())
    };
    assert_eq!(hit(&s), None);
    s.parameters.set("w", 5.0);
    s.update().unwrap();
    assert_eq!(hit(&s).as_deref(), Some("strip"));

    let eps_at = |s: &ContinuousStructure| {
        s.medium_at(Vec3::new(5.0, 5.0, 0.5), CoordSystem::Cartesian)
            .unwrap()
            .epsilon[0]
    };
    assert_eq!(eps_at(&s), 4.0);
    s.parameters.set("eps", 2.2);
    s.update().unwrap();
    assert_eq!(eps_at(&s), 2.2);

    // The parametric expressions survive the XML round-trip.
    let r = ContinuousStructure::from_xml_str(&s.to_xml_string()).unwrap();
    let m = r.property_by_name("sub").unwrap().kind.material().unwrap();
    assert_eq!(m.epsilon.value[0].expression(), Some("eps"));
    assert_eq!(m.base_value(MaterialQuantity::Epsilon, 0), 2.2);

    // A thicker substrate moves the mesh.
    s.parameters.set("h", 2.0);
    s.update().unwrap();
    mesh(&mut s);
    assert!(s.grid.snap(2, 2.0).is_some());
}

#[test]
fn cylindrical_model() {
    let mut s = ContinuousStructure::new();
    s.coord_system = CoordSystem::Cylindrical;
    s.grid.coord_system = CoordSystem::Cylindrical;
    let pi = core::f64::consts::PI;
    // Coaxial line: inner conductor, dielectric annulus, outer shell sector.
    let inner = s.add_property(Property::metal("inner"));
    s.add_primitive(
        inner,
        Primitive::new(BoxShape::new([0.0, -pi, 0.0], [1.0, pi, 10.0])).with_priority(10),
    )
    .unwrap();
    let diel = s.add_property(Property::material("ptfe", Material::new(2.1, 0.0)));
    s.add_primitive(diel, BoxShape::new([1.0, -pi, 0.0], [3.0, pi, 10.0]))
        .unwrap();
    s.update().unwrap();

    // Queries work in both coordinate systems.
    let name = |p: Vec3, cs| {
        s.property_at(p, cs, TypeMask::ANY)
            .map(|h| h.property.name.as_str())
    };
    assert_eq!(
        name(Vec3::new(2.0, 1.0, 5.0), CoordSystem::Cylindrical),
        Some("ptfe")
    );
    assert_eq!(
        name(Vec3::new(0.0, -2.0, 5.0), CoordSystem::Cartesian),
        Some("ptfe")
    );
    assert_eq!(
        name(Vec3::new(0.5, 0.0, 5.0), CoordSystem::Cartesian),
        Some("inner")
    );
    assert_eq!(
        name(Vec3::new(4.0, 0.0, 5.0), CoordSystem::Cylindrical),
        None
    );

    // Full revolutions contribute no angular lines; give the angular range
    // explicitly.
    MeshBuilder::new(MeshOptions::new(0.5).with_max_res([0.5, 2.0 * pi / 16.0, 1.0]))
        .add_lines(1, [-pi, pi])
        .apply(&mut s)
        .unwrap();
    assert_eq!(s.grid.coord_system, CoordSystem::Cylindrical);
    assert!(s.grid.num_lines()[1] >= 17);
    assert!(s.grid.snap(0, 1.0).is_some() && s.grid.snap(0, 3.0).is_some());

    let mg = discretize(&s, &DiscretizeOptions::default()).unwrap();
    assert_eq!(mg.coord_system, CoordSystem::Cylindrical);
    assert!(mg.pec_count() > 0);
    let eps_max = mg.epsilon[2].iter().copied().fold(0.0, f64::max);
    assert!((eps_max - 2.1).abs() < 1e-12);

    let r = ContinuousStructure::from_xml_str(&s.to_xml_string()).unwrap();
    assert_eq!(r.coord_system, CoordSystem::Cylindrical);
    assert_same(&s, &r);

    // Tessellated sectors are closed solids.
    let parts = tessellate_structure(&s, TypeMask::ANY, &TessellationOptions::default());
    for p in &parts {
        assert!(p.mesh.to_polyhedron().unwrap().is_closed());
    }
}

#[test]
fn round_shapes_and_wires() {
    let mut s = ContinuousStructure::new();
    let m = s.add_property(Property::metal("m"));
    s.add_primitive(m, Sphere::new([0.0, 0.0, 0.0], 2.0))
        .unwrap();
    s.add_primitive(m, Cylinder::new([5.0, 0.0, 0.0], [5.0, 0.0, 4.0], 1.0))
        .unwrap();
    s.add_primitive(
        m,
        Wire::new([[10.0, 0.0, 0.0], [10.0, 0.0, 3.0], [12.0, 0.0, 3.0]], 0.2),
    )
    .unwrap();
    s.update().unwrap();
    let r = ContinuousStructure::from_xml_str(&s.to_xml_string()).unwrap();
    assert_same(&s, &r);
    let parts = tessellate_structure(&r, TypeMask::ANY, &TessellationOptions::default());
    assert_eq!(parts.len(), 3);
    assert!(parts.iter().all(|p| !p.mesh.is_empty()));
    assert_eq!(parts[2].polyline.len(), 3);
}

/// A hand-written openEMS file with every property and primitive type.
const ALL_TYPES: &str = r#"<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>
<openEMS>
  <FDTD NumberOfTimesteps="10000" endCriteria="1e-4" f_max="3e9">
    <Excitation Type="0" f0="1.5e9" fc="1.5e9"/>
    <BoundaryCond xmin="PML_8" xmax="PML 8" ymin="MUR" ymax="PEC" zmin="1" zmax="3"/>
  </FDTD>
  <ContinuousStructure CoordSystem="0">
    <Properties>
      <Material ID="0" Name="mat" Isotropy="0">
        <FillColor R="10" G="20" B="30" a="255"/>
        <EdgeColor R="1" G="2" B="3" a="128"/>
        <Property Epsilon="2,3,4" Mue="1" Kappa="0.1" Sigma="0" Density="1000"/>
        <Weight Epsilon="1,max(1,x),1" Mue="1,1,1" Kappa="1" Sigma="1" Density="1"/>
        <Primitives>
          <Box Priority="1"><P1 X="0" Y="0" Z="0"/><P2 X="1" Y="1" Z="1"/></Box>
          <MultiBox Priority="1">
            <StartP X="2" Y="0" Z="0"/><EndP X="3" Y="1" Z="1"/>
            <StartP X="4" Y="0" Z="0"/><EndP X="5" Y="1" Z="1"/>
          </MultiBox>
          <Sphere Priority="2" Radius="1"><Center X="10" Y="0" Z="0"/></Sphere>
          <SphericalShell Priority="2" Radius="2" ShellWidth="0.2"><Center X="20" Y="0" Z="0"/></SphericalShell>
          <Cylinder Priority="2" Radius="0.5"><P1 X="30" Y="0" Z="0"/><P2 X="30" Y="0" Z="2"/></Cylinder>
          <CylindricalShell Priority="2" Radius="1" ShellWidth="0.1"><P1 X="40" Y="0" Z="0"/><P2 X="40" Y="0" Z="2"/></CylindricalShell>
        </Primitives>
      </Material>
      <LorentzMaterial ID="1" Name="drude" Isotropy="1">
        <Property Epsilon="1" EpsilonPlasmaFrequency="5e9" EpsilonLorPoleFrequency="0" EpsilonRelaxTime="1e-9"/>
        <Primitives>
          <Polygon Priority="3" NormDir="2" Elevation="0.5" QtyVertices="3">
            <Vertex X1="50" X2="0"/><Vertex X1="52" X2="0"/><Vertex X1="50" X2="2"/>
          </Polygon>
        </Primitives>
      </LorentzMaterial>
      <DebyeMaterial ID="2" Name="debye" Isotropy="1">
        <Property Epsilon="2" EpsilonDelta_1="10" EpsilonRelaxTime_1="1e-10" EpsilonDelta_2="5" EpsilonRelaxTime_2="1e-9"/>
        <Primitives>
          <LinPoly Priority="3" NormDir="1" Elevation="0" Length="2">
            <Vertex X1="0" X2="60"/><Vertex X1="1" X2="60"/><Vertex X1="1" X2="61"/><Vertex X1="0" X2="61"/>
          </LinPoly>
        </Primitives>
      </DebyeMaterial>
      <Metal ID="3" Name="pec">
        <Primitives>
          <RotPoly Priority="5" NormDir="1" RotAxisDir="2">
            <Vertex X1="0" X2="1"/><Vertex X1="1" X2="1"/><Vertex X1="1" X2="2"/><Vertex X1="0" X2="2"/>
            <Angles Start="0" Stop="6.283185307179586"/>
            <Transformation><Translate Argument="70,0,0"/></Transformation>
          </RotPoly>
          <Curve Priority="5"><Vertex X="80" Y="0" Z="0"/><Vertex X="81" Y="0" Z="0"/></Curve>
          <Wire Priority="5" WireRadius="0.1"><Vertex X="90" Y="0" Z="0"/><Vertex X="91" Y="0" Z="1"/></Wire>
          <Point Priority="5" X="100" Y="0" Z="0"/>
          <Polyhedron Priority="5">
            <Vertex>110,0,0</Vertex><Vertex>111,0,0</Vertex><Vertex>110,1,0</Vertex><Vertex>110,0,1</Vertex>
            <Face>0,2,1</Face><Face>0,1,3</Face><Face>0,3,2</Face><Face>1,2,3</Face>
          </Polyhedron>
          <UserDefined Priority="5"><Equation>(x-120)^2+y^2+z^2&lt;1</Equation></UserDefined>
          <Box Priority="5">
            <Transformation><Rotate_Z Argument="0.7853981633974483"/><Scale Argument="2"/></Transformation>
            <P1 X="0" Y="0" Z="200"/><P2 X="1" Y="1" Z="201"/>
          </Box>
        </Primitives>
      </Metal>
      <ConductingSheet ID="4" Name="copper" Conductivity="5.8e7" Thickness="35e-6">
        <Primitives><Box Priority="4"><P1 X="0" Y="0" Z="300"/><P2 X="1" Y="1" Z="300"/></Box></Primitives>
      </ConductingSheet>
      <Excitation ID="5" Name="exc" Number="1" Enabled="1" Type="0" Excite="0,0,1" PropDir="0,0,0" Delay="0" Frequency="0">
        <Weight X="1" Y="1" Z="cos(x)"/>
        <Primitives><Box Priority="0"><P1 X="0" Y="0" Z="400"/><P2 X="0" Y="1" Z="401"/></Box></Primitives>
      </Excitation>
      <ProbeBox ID="6" Name="probe" Visible="0" Number="0" Type="1" Weight="-1" NormDir="2" StartTime="0" StopTime="0">
        <FD_Samples>1e9,2e9</FD_Samples>
        <Attributes ModeFunctionX="sin(x)" ModeFunctionY="0" ModeFunctionZ="0"/>
        <Primitives><Box Priority="0"><P1 X="0" Y="0" Z="500"/><P2 X="1" Y="1" Z="500"/></Box></Primitives>
      </ProbeBox>
      <DumpBox ID="7" Name="dump" DumpType="10" DumpMode="2" FileType="1" OptResolution="1,1,1" SubSampling="2,2,2">
        <FD_Samples>1e9</FD_Samples>
        <Primitives><Box Priority="0"><P1 X="0" Y="0" Z="600"/><P2 X="10" Y="10" Z="600"/></Box></Primitives>
      </DumpBox>
      <LumpedElement ID="8" Name="rlc" Direction="2" Caps="1" R="50" C="1e-12" L="nan" LEtype="1">
        <Primitives><Box Priority="0"><P1 X="0" Y="0" Z="700"/><P2 X="1" Y="1" Z="701"/></Box></Primitives>
      </LumpedElement>
      <AbsorbingBC ID="9" Name="abc" NormalSignPositive="1" PhaseVelocity="3e8" AbsorbingBoundaryType="0">
        <Primitives><Box Priority="0"><P1 X="0" Y="0" Z="800"/><P2 X="10" Y="10" Z="800"/></Box></Primitives>
      </AbsorbingBC>
      <DiscMaterial ID="10" Name="voxels" File="model.h5" Scale="0.001">
        <Primitives><Box Priority="0"><P1 X="0" Y="0" Z="900"/><P2 X="1" Y="1" Z="901"/></Box></Primitives>
      </DiscMaterial>
      <Unknown ID="11" Name="misc" Property="whatever">
        <Primitives/>
      </Unknown>
    </Properties>
    <RectilinearGrid DeltaUnit="0.001" CoordSystem="0">
      <XLines>0,1,2</XLines>
      <YLines>0,1</YLines>
      <ZLines>0,1</ZLines>
    </RectilinearGrid>
    <BackgroundMaterial Epsilon="1" Mue="1" Kappa="0" Sigma="0"/>
    <ParameterSet/>
  </ContinuousStructure>
</openEMS>
"#;

#[test]
fn every_type_from_fixture() {
    let doc = OpenEmsDocument::from_xml_str(ALL_TYPES).unwrap();
    let s = &doc.structure;
    let kinds: Vec<&str> = s.properties().iter().map(|p| p.kind.type_name()).collect();
    assert_eq!(
        kinds,
        [
            "Material",
            "LorentzMaterial",
            "DebyeMaterial",
            "Metal",
            "ConductingSheet",
            "Excitation",
            "ProbeBox",
            "DumpBox",
            "LumpedElement",
            "AbsorbingBC",
            "DiscMaterial",
            "Unknown"
        ]
    );
    let shapes: Vec<&str> = s.primitives().map(|(_, p)| p.shape.type_name()).collect();
    for t in [
        "Box",
        "MultiBox",
        "Sphere",
        "SphericalShell",
        "Cylinder",
        "CylindricalShell",
        "Polygon",
        "LinPoly",
        "RotPoly",
        "Curve",
        "Wire",
        "Point",
        "Polyhedron",
        "UserDefined",
    ] {
        assert!(shapes.contains(&t), "missing {t}");
    }

    // FDTD settings.
    let f = doc.fdtd_settings().unwrap().unwrap();
    assert_eq!(f.number_of_timesteps, Some(10000));
    assert_eq!(
        f.excitation,
        Some(ExcitationSignal::Gaussian {
            f0: 1.5e9,
            fc: 1.5e9
        })
    );
    let b = f.boundary.unwrap().faces;
    assert_eq!(b[0], Boundary::Pml(8));
    assert_eq!(b[1], Boundary::Pml(8));
    assert_eq!(b[2], Boundary::Mur);
    assert_eq!(b[4], Boundary::Pmc);

    // Material details.
    let mat = s.property_by_name("mat").unwrap();
    assert_eq!(mat.fill_color.g, 20);
    assert_eq!(mat.edge_color.a, 128);
    let m = mat.kind.material().unwrap();
    assert!(!m.isotropic);
    assert_eq!(
        m.values_at(MaterialQuantity::Epsilon, Vec3::new(5.0, 0.0, 0.0)),
        [2.0, 15.0, 4.0]
    );
    assert_eq!(m.density.value(), 1000.0);
    match &s.property_by_name("debye").unwrap().kind {
        PropertyKind::DebyeMaterial(d) => assert_eq!(d.poles.len(), 2),
        k => panic!("{k:?}"),
    }
    match &s.property_by_name("rlc").unwrap().kind {
        PropertyKind::LumpedElement(l) => {
            assert_eq!(l.c.value(), 1e-12);
            assert!(l.l.value().is_nan());
            assert_eq!(l.le_type, Some(csx::properties::LumpedElementType::Series));
        }
        k => panic!("{k:?}"),
    }
    let probe = s.property_by_name("probe").unwrap();
    assert!(!probe.visible);
    assert_eq!(probe.attribute("ModeFunctionX"), Some("sin(x)"));

    // Geometry queries on several primitive types.
    let at = |p: Vec3| {
        s.property_at(p, CoordSystem::Cartesian, TypeMask::ANY)
            .map(|h| (h.property.name.as_str(), h.primitive.shape.type_name()))
    };
    assert_eq!(at(Vec3::new(4.5, 0.5, 0.5)), Some(("mat", "MultiBox")));
    assert_eq!(at(Vec3::new(10.0, 0.5, 0.0)), Some(("mat", "Sphere")));
    assert_eq!(
        at(Vec3::new(20.0, 2.05, 0.0)),
        Some(("mat", "SphericalShell"))
    );
    assert_eq!(at(Vec3::new(50.5, 0.5, 0.5)), Some(("drude", "Polygon")));
    assert_eq!(at(Vec3::new(60.5, 1.0, 0.5)), Some(("debye", "LinPoly")));
    // RotPoly: washer of radius 1..2 around z, translated to x = 70.
    assert_eq!(at(Vec3::new(70.0, 1.5, 0.5)), Some(("pec", "RotPoly")));
    assert_eq!(at(Vec3::new(70.0, 0.5, 0.5)), None);
    assert_eq!(at(Vec3::new(110.1, 0.1, 0.1)), Some(("pec", "Polyhedron")));
    assert_eq!(at(Vec3::new(120.2, 0.2, 0.2)), Some(("pec", "UserDefined")));
    // Rotated by 45° and scaled by 2: local (1,1) → (0, 2√2).
    assert_eq!(at(Vec3::new(0.0, 2.0, 401.0)), Some(("pec", "Box")));

    // Round-trip of the whole document.
    let xml = doc.to_xml_string();
    let again = OpenEmsDocument::from_xml_str(&xml).unwrap();
    assert_eq!(again.to_xml_string(), xml);
    assert_same(s, &again.structure);
    assert_eq!(
        again.fdtd_settings().unwrap().unwrap(),
        doc.fdtd_settings().unwrap().unwrap()
    );

    // Every closed solid tessellates to a closed mesh; nothing panics.
    let parts = tessellate_structure(s, TypeMask::ANY, &TessellationOptions::default());
    assert_eq!(parts.len(), s.primitives().count());
    let vtk = export_vtk(s, &VtkOptions::default());
    assert!(vtk.contains("LINES") && vtk.contains("POLYGONS"));
    for (p, (prop, prim)) in parts.iter().zip(s.primitives()) {
        if prim.dimension() == 3 && !matches!(prim.shape, Shape::UserDefined(_) | Shape::Wire(_)) {
            assert!(
                p.mesh.to_polyhedron().unwrap().is_closed(),
                "{} {}",
                prop.name,
                prim.shape.type_name()
            );
        }
    }

    // Discretisation over the fixture's small grid works.
    let mg = discretize(s, &DiscretizeOptions::default()).unwrap();
    assert_eq!(mg.dims, [3, 2, 2]);
}

/// Mirrors the example in README.md.
#[test]
fn readme_example() {
    use csx::primitives::BoxShape;
    use csx::properties::Material;
    use csx::{ContinuousStructure, CoordSystem, Property, TypeMask, Vec3};

    let mut model = ContinuousStructure::new();
    let fr4 = model.add_property(Property::material("FR4", Material::new(4.3, 0.0)));
    model
        .add_primitive(fr4, BoxShape::new([-10.0, -10.0, 0.0], [10.0, 10.0, 1.6]))
        .unwrap();

    // Write and read back a CSXCAD XML file (feature `xml`).
    let xml = model.to_xml_string();
    let model = ContinuousStructure::from_xml_str(&xml).unwrap();

    let hit = model.property_at(
        Vec3::new(0.0, 0.0, 0.5),
        CoordSystem::Cartesian,
        TypeMask::ANY,
    );
    assert_eq!(hit.unwrap().property.name, "FR4");
}
