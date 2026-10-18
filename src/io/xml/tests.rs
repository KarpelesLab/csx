use super::*;
use crate::coords::CoordSystem;
use crate::math::Vec3;
use crate::param::Scalar;
use crate::primitives::Shape;
use crate::properties::{DumpMode, DumpType, LumpedElementType};
use crate::properties::{MaterialQuantity, PropertyKind};
use crate::{PropertyType, TypeMask};

const SAMPLE: &str = r#"<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>
<openEMS>
  <FDTD NumberOfTimesteps="30000" endCriteria="1e-05" f_max="7e9">
    <Excitation Type="0" f0="0" fc="7e9"/>
    <BoundaryCond xmin="PML_8" xmax="PML_8" ymin="MUR" ymax="PEC" zmin="PEC" zmax="PMC"/>
  </FDTD>
  <ContinuousStructure CoordSystem="0">
    <Properties>
      <Material Name="sub" Isotropy="1" Custom="kept">
        <Property Epsilon="3.38" Kappa="0.0004"/>
        <Weight Epsilon="1+z, pow(x,2) ,1"/>
        <Primitives>
          <Box Priority="1">
            <P1 X="-10" Y="-10" Z="0"/>
            <P2 X="10" Y="10" Z="h"/>
          </Box>
        </Primitives>
      </Material>
      <Metal Name="patch">
        <FillColor R="200" G="10" B="1" a="128"/>
        <Primitives>
          <Polygon Priority="10" Elevation="1.5" NormDir="2" CoordSystem="0">
            <Vertex X1="0" X2="0"/><Vertex X1="4" X2="0"/><Vertex X1="4" X2="3"/>
          </Polygon>
          <Sphere Priority="11" Radius="1">
            <Center X="0" Y="0" Z="0"/>
            <Transformation>
              <Scale Argument="2,1,1"/>
              <Translate Argument="20, 0, 0"/>
              <Rotate_Z Argument="0.5"/>
            </Transformation>
          </Sphere>
        </Primitives>
      </Metal>
      <LumpedElement Name="R1" Direction="2" Caps="1" R="50" C="nan" L="-nan(ind)" LEtype="0.000000e+00">
        <Primitives><Box Priority="5"><P1 X="0" Y="0" Z="0"/><P2 X="0" Y="0" Z="1.5"/></Box></Primitives>
      </LumpedElement>
      <Excitation Name="ex" Type="0" Excite="0,0,-1" Number="0" Enabled="1" Delay="0">
        <Weight X="1" Y="1" Z="sin(x)"/>
        <Primitives><Box Priority="5"><P1 X="0" Y="0" Z="0"/><P2 X="0" Y="0" Z="1.5"/></Box></Primitives>
      </Excitation>
      <ProbeBox Name="ut" Type="0" Weight="-1" NormDir="-1" Visible="0">
        <FD_Samples>1e9,2e9</FD_Samples>
        <Attributes ModeFunctionX="sin(x)" ModeFunctionY="0" ModeFunctionZ="0"/>
        <Primitives><Box Priority="0"><P1 X="0" Y="0" Z="0"/><P2 X="0" Y="0" Z="1.5"/></Box></Primitives>
      </ProbeBox>
      <DumpBox Name="Et" DumpType="10" DumpMode="2" SubSampling="2,2,4" OptResolution="1">
        <Primitives><Box Priority="0"><P1 X="-10" Y="0" Z="-10"/><P2 X="10" Y="0" Z="30"/></Box></Primitives>
      </DumpBox>
      <DebyeMaterial Name="skin" Isotropy="1">
        <Property Epsilon="21.6,21.6,21.6" EpsilonDelta_1="19.3,19.3,19.3" EpsilonRelaxTime_1="7e-12" EpsilonDelta_2="2" EpsilonRelaxTime_2="1e-9"/>
        <Primitives/>
      </DebyeMaterial>
      <LorentzMaterial Name="drude" Isotropy="1">
        <Property Epsilon="4.2" EpsilonPlasmaFrequency="5e9" EpsilonLorPoleFrequency="12.57e9" EpsilonRelaxTime="2.294e-11"/>
        <Primitives/>
      </LorentzMaterial>
      <ConductingSheet Name="cu" Conductivity="5.8e7" Thickness="35e-6">
        <Primitives>
          <MultiBox Priority="1"><StartP X="0" Y="0" Z="0"/><EndP X="1" Y="1" Z="0"/><StartP X="2" Y="2" Z="0"/><EndP X="3" Y="3" Z="0"/></MultiBox>
        </Primitives>
      </ConductingSheet>
      <Metal Name="misc">
        <Primitives>
          <Cylinder Priority="1" Radius="0.5"><P1 X="0" Y="0" Z="0"/><P2 X="0" Y="0" Z="2"/></Cylinder>
          <CylindricalShell Priority="1" Radius="1" ShellWidth="0.1"><P1 X="0" Y="0" Z="0"/><P2 X="0" Y="0" Z="2"/></CylindricalShell>
          <SphericalShell Priority="1" Radius="1" ShellWidth="0.1"><Center X="5" Y="5" Z="5"/></SphericalShell>
          <LinPoly Priority="1" Elevation="1" Length="3" NormDir="2"><Vertex X1="10" X2="10"/><Vertex X1="10" X2="20"/><Vertex X1="30" X2="20"/></LinPoly>
          <RotPoly Priority="1" NormDir="1" RotAxisDir="2"><Vertex X1="0" X2="1"/><Vertex X1="1" X2="1"/><Vertex X1="1" X2="2"/><Angles Start="0" Stop="3.14"/></RotPoly>
          <Curve Priority="0"><Vertex X="0" Y="0" Z="0"></Vertex><Vertex X="1" Y="0" Z="0"/></Curve>
          <Wire Priority="0" WireRadius="0.1"><Vertex X="0" Y="0" Z="0"/><Vertex X="1" Y="0" Z="0"/></Wire>
          <Point Priority="1" X="1" Y="1" Z="1"/>
          <Polyhedron Priority="2">
            <Vertex>0,0,0</Vertex><Vertex X="1" Y="0" Z="0"/><Vertex>0,1,0</Vertex><Vertex>0,0,1</Vertex>
            <Face>0,2,1</Face><Face>0,1,3</Face><Face>0,3,2</Face><Face>1,2,3</Face>
          </Polyhedron>
          <PolyhedronReader Priority="10" FileName="model.stl" FileType="STL"/>
          <UserDefined Priority="3"><CoordShift X="1" Y="0" Z="0"/><Equation>x^2+y^2 &lt; 1</Equation></UserDefined>
          <Box Priority="1" CoordSystem="1"><P1 X="1" Y="0" Z="0"/><P2 X="2" Y="1.5" Z="1"/></Box>
        </Primitives>
      </Metal>
      <ResBox Name="weird" Foo="bar"><Primitives/></ResBox>
    </Properties>
    <RectilinearGrid DeltaUnit="0.001" CoordSystem="0">
      <XLines>-10,0,10,h*2</XLines>
      <YLines Qty="2">-10, 10</YLines>
      <ZLines/>
    </RectilinearGrid>
    <BackgroundMaterial Epsilon="1" Mue="1" Kappa="0" Sigma="0"/>
    <ParameterSet>
      <Parameter Type="Const" name="h" value="1.524" Sweep="1"/>
    </ParameterSet>
  </ContinuousStructure>
</openEMS>
"#;

#[test]
fn read_sample() {
    let doc = OpenEmsDocument::from_xml_str(SAMPLE).unwrap();
    assert_eq!(doc.fdtd.as_ref().unwrap().attr("f_max"), Some("7e9"));
    let s = &doc.structure;
    assert_eq!(s.grid.lines(0), &[-10.0, 0.0, 3.048, 10.0]);
    assert_eq!(s.grid.lines(1), &[-10.0, 10.0]);
    assert!(s.grid.lines(2).is_empty());
    assert_eq!(s.grid.delta_unit, 0.001);

    let sub = s.property_by_name("sub").unwrap();
    assert_eq!(
        sub.extra_attributes,
        vec![("Custom".to_owned(), "kept".to_owned())]
    );
    let m = sub.kind.material().unwrap();
    assert_eq!(m.base_value(MaterialQuantity::Epsilon, 2), 3.38);
    assert_eq!(
        m.quantity(MaterialQuantity::Epsilon).weight[1].source(),
        "pow(x,2)"
    );
    assert_eq!(
        m.value_at(MaterialQuantity::Epsilon, 0, Vec3::new(0.0, 0.0, 1.0)),
        6.76
    );
    // box top at parameter h
    let hit = s.property_at(
        Vec3::new(0.0, 0.0, 1.5),
        CoordSystem::Cartesian,
        TypeMask::MATERIALS,
    );
    assert_eq!(hit.unwrap().property.name, "sub");

    let patch = s.property_by_name("patch").unwrap();
    assert_eq!(patch.fill_color.r, 200);
    assert_eq!(patch.fill_color.a, 128);
    // transformed sphere: scaled ×2 in x, moved to x=20, rotated by 0.5 rad
    let c = Vec3::new(20.0 * 0.5f64.cos(), 20.0 * 0.5f64.sin(), 0.0);
    assert!(
        s.property_at(c, CoordSystem::Cartesian, TypeMask::ANY)
            .is_some()
    );

    let PropertyKind::LumpedElement(l) = &s.property_by_name("R1").unwrap().kind else {
        panic!()
    };
    assert_eq!(l.r.value(), 50.0);
    assert!(l.c.value().is_nan() && l.l.value().is_nan());
    assert_eq!(l.le_type, Some(LumpedElementType::Parallel));

    let PropertyKind::Excitation(e) = &s.property_by_name("ex").unwrap().kind else {
        panic!()
    };
    assert_eq!(e.excite[2].value(), -1.0);
    assert_eq!(e.weight[2].source(), "sin(x)");

    let ut = s.property_by_name("ut").unwrap();
    assert!(!ut.visible);
    assert_eq!(ut.attribute("ModeFunctionX"), Some("sin(x)"));
    let PropertyKind::ProbeBox(p) = &ut.kind else {
        panic!()
    };
    assert_eq!(p.frequencies, vec![1e9, 2e9]);
    assert_eq!(p.weight.value(), -1.0);
    assert_eq!(p.norm_dir, None);

    let PropertyKind::DumpBox(d) = &s.property_by_name("Et").unwrap().kind else {
        panic!()
    };
    assert_eq!(
        (d.dump_type, d.dump_mode, d.sub_sampling),
        (DumpType::EFieldFd, DumpMode::Cell, Some([2, 2, 4]))
    );
    assert_eq!(d.opt_resolution.as_ref().unwrap()[2].value(), 1.0);

    let PropertyKind::DebyeMaterial(db) = &s.property_by_name("skin").unwrap().kind else {
        panic!()
    };
    assert_eq!(db.poles.len(), 2);
    assert_eq!(db.poles[1].eps_relax_time.base(0), 1e-9);
    let PropertyKind::LorentzMaterial(lz) = &s.property_by_name("drude").unwrap().kind else {
        panic!()
    };
    assert_eq!(lz.poles.len(), 1);
    assert_eq!(lz.poles[0].eps_plasma_frequency.base(1), 5e9);

    let misc = s.property_by_name("misc").unwrap();
    let names: Vec<_> = misc
        .primitives
        .iter()
        .map(|p| p.shape.type_name())
        .collect();
    assert_eq!(
        names,
        [
            "Cylinder",
            "CylindricalShell",
            "SphericalShell",
            "LinPoly",
            "RotPoly",
            "Curve",
            "Wire",
            "Point",
            "Polyhedron",
            "PolyhedronReader",
            "UserDefined",
            "Box"
        ]
    );
    let Shape::Polyhedron(ph) = &misc.primitives[8].shape else {
        panic!()
    };
    assert!(ph.is_closed());
    assert_eq!(
        misc.primitives[11].coord_system,
        Some(CoordSystem::Cylindrical)
    );
    assert!(misc.primitives[10].contains(Vec3::new(1.5, 0.0, 0.0), CoordSystem::Cartesian));

    let weird = s.property_by_name("weird").unwrap();
    assert_eq!(weird.property_type(), PropertyType::Unknown);
    // Regression: foreign property elements were renamed to <Unknown>.
    let out = s.to_xml_string();
    assert!(
        out.contains("<ResBox ") && out.contains("Foo=\"bar\""),
        "{out}"
    );
    assert!(!out.contains("Property=\"ResBox\""));
    assert_eq!(s.parameters.get("h"), Some(1.524));
}

#[test]
fn roundtrip_is_stable() {
    let doc = OpenEmsDocument::from_xml_str(SAMPLE).unwrap();
    let first = doc.to_xml_string();
    let again = OpenEmsDocument::from_xml_str(&first).unwrap();
    let second = again.to_xml_string();
    assert_eq!(first, second);
    assert_eq!(again.fdtd, doc.fdtd);
    assert_eq!(
        again.structure.properties().len(),
        doc.structure.properties().len()
    );
    // Structures compare equal except for NaN-valued scalars.
    for (a, b) in doc
        .structure
        .properties()
        .iter()
        .zip(again.structure.properties())
    {
        if !matches!(a.kind, PropertyKind::LumpedElement(_)) {
            assert_eq!(a, b, "{}", a.name);
        }
    }
    assert_eq!(doc.structure.grid, again.structure.grid);
    assert_eq!(doc.structure.parameters, again.structure.parameters);
}

#[test]
fn structure_root_and_files() {
    let s = ContinuousStructure::from_xml_str(SAMPLE).unwrap();
    let xml = s.to_xml_string();
    assert!(xml.starts_with("<?xml"));
    assert!(xml.contains("<ContinuousStructure CoordSystem=\"0\">"));
    let t = ContinuousStructure::from_xml_str(&xml).unwrap();
    assert_eq!(t.to_xml_string(), xml);

    let dir = std::env::temp_dir().join(format!("csx-xml-test-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("model.xml");
    t.write_file(&path).unwrap();
    let r = ContinuousStructure::read_file(&path).unwrap();
    assert_eq!(r.base_dir.as_deref(), Some(dir.as_path()));
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn degree_transforms_written_in_radians() {
    let mut s = ContinuousStructure::new();
    let id = s.add_property(crate::Property::metal("m"));
    let mut t = crate::Transform::new();
    t.set_angles_in_degrees(true);
    t.rotate_axis(2, 90.0);
    let prim = crate::Primitive::new(crate::primitives::BoxShape::new(
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0],
    ))
    .with_transform(t);
    s.add_primitive(id, prim).unwrap();
    let r = ContinuousStructure::from_xml_str(&s.to_xml_string()).unwrap();
    let p = &r.properties()[0].primitives[0];
    assert!(!p.transform.angles_in_degrees());
    let q = p.transform.apply(Vec3::unit(0));
    assert!((q - Vec3::unit(1)).norm() < 1e-12);
}

#[test]
fn vertex_expressions() {
    // Regression: expression-valued vertices made the whole file unreadable.
    let xml = r#"<ContinuousStructure CoordSystem="0">
      <Properties>
        <Metal Name="m"><Primitives>
          <Polygon Priority="1" NormDir="2" Elevation="0">
            <Vertex X1="0" X2="0"/><Vertex X1="w/2" X2="0"/><Vertex X1="w/2" X2="w"/>
          </Polygon>
          <Polyhedron><Vertex>0,0,0</Vertex><Vertex>w,0,0</Vertex><Vertex>0,w,0</Vertex>
            <Face>0,1,2</Face></Polyhedron>
        </Primitives></Metal>
      </Properties>
      <ParameterSet><Parameter Type="Const" name="w" value="4"/></ParameterSet>
    </ContinuousStructure>"#;
    let s = ContinuousStructure::from_xml_str(xml).unwrap();
    let p = &s.property_by_name("m").unwrap().primitives;
    let Shape::Polygon(poly) = &p[0].shape else {
        panic!()
    };
    assert_eq!(poly.vertices.values(), [(0.0, 0.0), (2.0, 0.0), (2.0, 4.0)]);
    assert_eq!(poly.vertices.points()[2][1].expression(), Some("w"));
    let Shape::Polyhedron(ph) = &p[1].shape else {
        panic!()
    };
    assert_eq!(ph.vertices()[1], Vec3::new(4.0, 0.0, 0.0));
    assert!(ContinuousStructure::from_xml_str(&xml.replace("w/2", "q/2")).is_err());
}

#[test]
fn errors() {
    assert!(ContinuousStructure::from_xml_str("<foo/>").is_err());
    assert!(ContinuousStructure::from_xml_str("<ContinuousStructure CoordSystem=\"7\"/>").is_err());
    let bad_prim = "<ContinuousStructure><Properties><Metal Name=\"m\"><Primitives><Blob/></Primitives></Metal></Properties></ContinuousStructure>";
    assert!(ContinuousStructure::from_xml_str(bad_prim).is_err());
    let bad_box = "<ContinuousStructure><Properties><Metal Name=\"m\"><Primitives><Box><P1 X=\"0\"/></Box></Primitives></Metal></Properties></ContinuousStructure>";
    assert!(ContinuousStructure::from_xml_str(bad_box).is_err());
    let bad_face = "<ContinuousStructure><Properties><Metal Name=\"m\"><Primitives><Polyhedron><Vertex>0,0,0</Vertex><Face>0,1,2</Face></Polyhedron></Primitives></Metal></Properties></ContinuousStructure>";
    assert!(ContinuousStructure::from_xml_str(bad_face).is_err());
    let empty = ContinuousStructure::from_xml_str("<ContinuousStructure/>").unwrap();
    assert!(empty.properties().is_empty());
}

#[test]
fn typed_ids_round_trip() {
    use crate::properties::{AbcType, Axis, ExcitationType, ProbeType};
    let xml = r#"<ContinuousStructure><Properties>
        <Excitation Name="ex" Type="7"/>
        <ProbeBox Name="p" Type="1.000000e+00" NormDir="2"/>
        <DumpBox Name="d" DumpType="99" DumpMode="0" FileType="5"/>
        <LumpedElement Name="l" Direction="1.000000e+00" LEtype="1.000000e+00" R="50"/>
        <AbsorbingBC Name="a" AbsorbingBoundaryType="MUR_1ST_SA"/>
    </Properties></ContinuousStructure>"#;
    let s = ContinuousStructure::from_xml_str(xml).unwrap();
    let kind = |n: &str| s.property_by_name(n).unwrap().kind.clone();
    let PropertyKind::Excitation(e) = kind("ex") else {
        panic!()
    };
    assert_eq!(e.excitation_type, ExcitationType::Other(7));
    let PropertyKind::ProbeBox(p) = kind("p") else {
        panic!()
    };
    assert_eq!(
        (p.probe_type, p.norm_dir),
        (ProbeType::Current, Some(Axis::Z))
    );
    let PropertyKind::DumpBox(d) = kind("d") else {
        panic!()
    };
    assert_eq!(
        (d.dump_type, d.dump_mode, d.file_type),
        (
            DumpType::Other(99),
            DumpMode::NoInterpolation,
            crate::properties::DumpFileType::Other(5)
        )
    );
    let PropertyKind::LumpedElement(l) = kind("l") else {
        panic!()
    };
    assert_eq!(
        (l.direction, l.le_type),
        (Axis::Y, Some(LumpedElementType::Series))
    );
    let PropertyKind::AbsorbingBC(a) = kind("a") else {
        panic!()
    };
    assert_eq!(a.abc_type, AbcType::Mur1stSa);

    let out = s.to_xml_string();
    assert!(out.contains("DumpType=\"99\"") && out.contains("FileType=\"5\""));
    let again = ContinuousStructure::from_xml_str(&out).unwrap();
    assert_eq!(again.to_xml_string(), out);

    let bad = r#"<ContinuousStructure><Properties><LumpedElement Name="l" Direction="3"/></Properties></ContinuousStructure>"#;
    assert!(ContinuousStructure::from_xml_str(bad).is_err());
    let bad = r#"<ContinuousStructure><Properties><ProbeBox Name="p" NormDir="5"/></Properties></ContinuousStructure>"#;
    assert!(ContinuousStructure::from_xml_str(bad).is_err());
}

/// Parses every sample file found in `CSX_XML_SAMPLES` (a directory of
/// third-party openEMS/CSXCAD files, not distributed with the crate) and
/// checks that writing and re-reading is stable.
#[test]
fn external_samples() {
    let Some(dir) = std::env::var_os("CSX_XML_SAMPLES") else {
        return;
    };
    let mut n = 0;
    for entry in std::fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
        if ext != "xml" && ext != "csx" {
            continue;
        }
        let doc =
            OpenEmsDocument::read_file(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let first = doc.to_xml_string();
        let again = OpenEmsDocument::from_xml_str(&first)
            .unwrap_or_else(|e| panic!("{} (rewritten): {e}", path.display()));
        assert_eq!(first, again.to_xml_string(), "{}", path.display());
        let settings = doc
            .fdtd_settings()
            .unwrap_or_else(|e| panic!("{} FDTD: {e}", path.display()));
        if let Some(settings) = settings {
            let mut typed = doc.clone();
            typed.set_fdtd_settings(&settings);
            let reread = OpenEmsDocument::from_xml_str(&typed.to_xml_string()).unwrap();
            assert_eq!(
                reread.fdtd_settings().unwrap().unwrap(),
                settings,
                "{}",
                path.display()
            );
        }
        n += 1;
    }
    assert!(n > 0);
}

#[test]
fn parametric_polygon_round_trip() {
    use crate::primitives::{LinPoly, PolyVertices};
    let mut s = ContinuousStructure::new();
    s.parameters.set("w", 3.0);
    let id = s.add_property(crate::Property::metal("m"));
    let verts: PolyVertices = [
        [Scalar::new(0.0), Scalar::new(0.0)],
        [Scalar::expr("w"), Scalar::new(0.0)],
        [Scalar::expr("w"), Scalar::expr("w*2")],
    ]
    .into_iter()
    .collect();
    s.add_primitive(id, LinPoly::new(2, 0.0, 1.0, verts))
        .unwrap();
    s.update().unwrap();
    let xml = s.to_xml_string();
    assert!(xml.contains(r#"X2="w*2""#), "{xml}");
    let r = ContinuousStructure::from_xml_str(&xml).unwrap();
    let Shape::LinPoly(lp) = &r.property_by_name("m").unwrap().primitives[0].shape else {
        panic!()
    };
    assert_eq!(lp.vertices.values()[2], (3.0, 6.0));
    assert_eq!(lp.vertices.points()[2][1].expression(), Some("w*2"));
    assert_eq!(r.to_xml_string(), xml);
    // changing the parameter changes the geometry after update()
    let mut r = r;
    r.parameters.set("w", 10.0);
    r.update().unwrap();
    let hit = r.property_at(
        Vec3::new(8.0, 1.0, 0.5),
        crate::CoordSystem::Cartesian,
        crate::TypeMask::ANY,
    );
    assert!(hit.is_some());
}

#[test]
fn unknown_children_survive_roundtrip() {
    let xml = r#"<ContinuousStructure CoordSystem="0">
  <Properties>
    <Material Name="m">
      <Property Epsilon="2"/>
      <CustomData a="1"><Nested b="2">text &amp; more</Nested></CustomData>
      <Primitives>
        <Box Priority="1">
          <P1 X="0" Y="0" Z="0"/><P2 X="1" Y="1" Z="1"/>
          <Annotation note="primitive level"/>
        </Box>
        <Point X="1" Y="2" Z="3"><P X="9" Y="9" Z="9"/></Point>
      </Primitives>
    </Material>
    <ResBox Name="foreign"><Inner/></ResBox>
    <ProbeBox Name="p"><Weight X="1"/><Primitives/></ProbeBox>
  </Properties>
  <Extension version="2"><Deep><Deeper/></Deep></Extension>
</ContinuousStructure>"#;
    let s = ContinuousStructure::from_xml_str(xml).unwrap();
    assert_eq!(s.extra_elements.len(), 1);
    assert_eq!(s.extra_elements[0].name, "Extension");
    let m = s.property_by_name("m").unwrap();
    assert_eq!(m.extra_elements.len(), 1);
    assert_eq!(m.extra_elements[0].name, "CustomData");
    assert_eq!(m.primitives[0].extra_elements[0].name, "Annotation");
    // `<P>` is unknown when the point uses attributes.
    assert_eq!(m.primitives[1].extra_elements[0].name, "P");
    assert_eq!(
        s.property_by_name("foreign").unwrap().extra_elements[0].name,
        "Inner"
    );
    // Excitation-only children are foreign on other kinds.
    assert_eq!(
        s.property_by_name("p").unwrap().extra_elements[0].name,
        "Weight"
    );

    let out = s.to_xml_string();
    let r = ContinuousStructure::from_xml_str(&out).unwrap();
    assert_eq!(r, s);
    assert_eq!(r.to_xml_string(), out);
    let nested = &r.property_by_name("m").unwrap().extra_elements[0];
    assert_eq!(nested.attr("a"), Some("1"));
    assert_eq!(nested.child("Nested").unwrap().text(), "text & more");
}

#[test]
fn grid_and_background_extra_attributes_survive() {
    let xml = r#"<ContinuousStructure CoordSystem="0">
  <RectilinearGrid DeltaUnit="0.001" CoordSystem="0" Vendor="x">
    <XLines>0,1</XLines><YLines>0,1</YLines><ZLines>0,1</ZLines>
  </RectilinearGrid>
  <BackgroundMaterial Epsilon="2" Mue="1" Kappa="0" Sigma="0" Custom="7"/>
  <Properties/>
</ContinuousStructure>"#;
    let s = ContinuousStructure::from_xml_str(xml).unwrap();
    assert_eq!(
        s.grid.extra_attributes,
        vec![("Vendor".to_owned(), "x".to_owned())]
    );
    assert_eq!(
        s.background.extra_attributes,
        vec![("Custom".to_owned(), "7".to_owned())]
    );
    let out = s.to_xml_string();
    let r = ContinuousStructure::from_xml_str(&out).unwrap();
    assert_eq!(r, s);
}

/// Hand-written fixture in the layout CSXCAD 0.7 itself writes (verified
/// black-box with its Python API).
const CSXCAD_LAYOUT: &str = r#"<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>
<ContinuousStructure CoordSystem="0">
    <RectilinearGrid DeltaUnit="0.001" CoordSystem="0">
        <XLines Qty="2">0,1</XLines>
        <YLines Qty="2">0,1</YLines>
        <ZLines Qty="2">0,1</ZLines>
    </RectilinearGrid>
    <BackgroundMaterial Epsilon="1" Mue="1" Kappa="0" Sigma="0" />
    <ParameterSet>
        <Parameter Type="Const" name="w" Sweep="1" value="2" />
        <Parameter Type="Linear" name="h" Sweep="1" value="1" min="0" max="2" step="0.5" />
    </ParameterSet>
    <Properties>
        <Metal ID="0" Name="m">
            <FillColor R="1" G="2" B="3" a="255" />
            <EdgeColor R="1" G="2" B="3" a="255" />
            <Primitives>
                <Box Priority="0">
                    <Transformation>
                        <Scale3 Argument="1,2,3" />
                        <Translate Argument="5,6,7" />
                        <Scale Argument="2" />
                    </Transformation>
                    <P1 X="0.000000e+00" Y="0.000000e+00" Z="0.000000e+00" />
                    <P2 X="1.000000e+00" Y="1.000000e+00" Z="1.000000e+00" />
                </Box>
                <Box Priority="0">
                    <Transformation>
                        <Matrix Argument="0,-1,0,10,1,0,0,20,0,0,1,30,0,0,0,1" />
                    </Transformation>
                    <P1 X="1" Y="0" Z="0" />
                    <P2 X="2" Y="0.5" Z="0.5" />
                </Box>
                <MultiBox Priority="0" QtyBox="1">
                    <StartP X="0" Y="0" Z="0" />
                    <EndP X="1" Y="1" Z="1" />
                </MultiBox>
            </Primitives>
        </Metal>
        <LorentzMaterial ID="1" Name="lor" Isotropy="1">
            <FillColor R="1" G="2" B="3" a="123" />
            <EdgeColor R="1" G="2" B="3" a="123" />
            <Primitives />
            <Property Epsilon="2,1,1" Mue="1,1,1" Kappa="0,0,0" Sigma="0,0,0" Density="0" EpsilonPlasmaFrequency_1="1e9,2e9,3e9" MuePlasmaFrequency_1="0,0,0" EpsilonLorPoleFrequency_1="4e9,0,0" MueLorPoleFrequency_1="0,0,0" EpsilonRelaxTime_1="1e-9,0,0" MueRelaxTime_1="0,0,0" EpsilonPlasmaFrequency_2="5e9,0,0" MuePlasmaFrequency_2="0,0,0" EpsilonLorPoleFrequency_2="0,0,0" MueLorPoleFrequency_2="0,0,0" EpsilonRelaxTime_2="0,0,0" MueRelaxTime_2="0,0,0" />
            <Weight Epsilon="1,1,1" Mue="1,1,1" Kappa="1,1,1" Sigma="1,1,1" Density="1" EpsilonPlasmaFrequency_1="x,1,1" MuePlasmaFrequency_1="1,1,1" EpsilonLorPoleFrequency_1="1,1,1" MueLorPoleFrequency_1="1,1,1" EpsilonRelaxTime_1="1,1,1" MueRelaxTime_1="1,1,1" EpsilonPlasmaFrequency_2="1,1,1" MuePlasmaFrequency_2="1,1,1" EpsilonLorPoleFrequency_2="1,1,1" MueLorPoleFrequency_2="1,1,1" EpsilonRelaxTime_2="y,1,1" MueRelaxTime_2="1,1,1" />
        </LorentzMaterial>
        <DiscMaterial ID="2" Name="disc" Isotropy="1">
            <FillColor R="1" G="2" B="3" a="123" />
            <EdgeColor R="1" G="2" B="3" a="123" />
            <Primitives />
            <Property Epsilon="1,1,1" Mue="1,1,1" Kappa="0,0,0" Sigma="0,0,0" Density="0" />
            <Weight Epsilon="1,1,1" Mue="1,1,1" Kappa="1,1,1" Sigma="1,1,1" Density="1" />
            <DiscFile Type="0" File="mat.h5" UseDBBackground="0" Scale="0.001" />
        </DiscMaterial>
        <LumpedElement ID="3" Name="le" Direction="2" Caps="1" R="5.000000e+01" C="1.000000e-12" LEtype="0.000000e+00">
            <FillColor R="1" G="2" B="3" a="255" />
            <EdgeColor R="1" G="2" B="3" a="255" />
            <Primitives />
        </LumpedElement>
    </Properties>
</ContinuousStructure>
"#;

/// Index of `needle` in `hay` at or after `from`.
fn find_from(hay: &str, from: usize, needle: &str) -> usize {
    from + hay[from..]
        .find(needle)
        .unwrap_or_else(|| panic!("{needle} missing"))
}

#[test]
fn csxcad_layout() {
    use crate::param::ParameterKind;
    let s = ContinuousStructure::from_xml_str(CSXCAD_LAYOUT).unwrap();
    let h = s.parameters.parameter("h").unwrap();
    assert_eq!(
        h.kind,
        ParameterKind::Linear {
            start: 0.0,
            stop: 2.0,
            step: 0.5
        }
    );
    assert_eq!(s.parameters.get("w"), Some(2.0));

    let m = &s.properties()[0];
    // Scale3(1,2,3), Translate(5,6,7), Scale(2) in list order.
    let t = &m.primitives[0].transform;
    assert_eq!(t.ops()[0].name(), "Scale3");
    assert_eq!(
        t.apply(Vec3::new(1.0, 1.0, 1.0)),
        Vec3::new(12.0, 16.0, 20.0)
    );
    // Matrix: row-major, acting on column vectors.
    let q = m.primitives[1].transform.apply(Vec3::new(1.0, 0.0, 0.0));
    assert!((q - Vec3::new(10.0, 21.0, 30.0)).norm() < 1e-12);
    assert!(
        s.property_at(
            Vec3::new(9.75, 21.5, 30.25),
            CoordSystem::Cartesian,
            TypeMask::ANY
        )
        .is_some()
    );
    assert!(
        m.primitives[2].extra_attributes.is_empty(),
        "QtyBox consumed"
    );

    let PropertyKind::LorentzMaterial(l) = &s.properties()[1].kind else {
        panic!("lorentz expected")
    };
    assert_eq!(l.poles.len(), 2);
    assert_eq!(l.poles[0].eps_plasma_frequency.base(2), 3e9);
    assert_eq!(l.poles[0].eps_lorentz_pole_frequency.base(0), 4e9);
    assert_eq!(l.poles[0].eps_plasma_frequency.weight[0].source(), "x");
    assert_eq!(l.poles[1].eps_relax_time.weight[0].source(), "y");
    assert_eq!(l.poles[1].eps_plasma_frequency.base(0), 5e9);

    let PropertyKind::DiscMaterial(d) = &s.properties()[2].kind else {
        panic!("disc expected")
    };
    assert_eq!(d.file, "mat.h5");
    assert_eq!(d.scale.value(), 0.001);
    assert!(s.properties()[2].extra_attributes.is_empty());

    let PropertyKind::LumpedElement(le) = &s.properties()[3].kind else {
        panic!("lumped expected")
    };
    assert!(le.l.value().is_nan());

    // The writer reproduces the CSXCAD layout.
    let out = s.to_xml_string();
    for needle in [
        r#"<Parameter Type="Linear" name="h" Sweep="1" value="1" min="0" max="2" step="0.5""#,
        r#"<Scale3 Argument="1,2,3""#,
        r#"<Matrix Argument="0,-1,0,10,1,0,0,20,0,0,1,30,0,0,0,1""#,
        r#"<MultiBox Priority="0" QtyBox="1">"#,
        r#"<DiscFile Type="0" File="mat.h5" Scale="1e-3" UseDBBackground="0""#,
        "EpsilonRelaxTime_2=\"y,1,1\"",
    ] {
        assert!(out.contains(needle), "{needle} missing in\n{out}");
    }
    assert!(!out.contains(" L=\""), "unset inductance is omitted");
    for prop in ["<LorentzMaterial", "<DiscMaterial"] {
        let start = find_from(&out, 0, prop);
        let prims = find_from(&out, start, "<Primitives");
        let attrs = find_from(&out, start, "<Property ");
        let weight = find_from(&out, start, "<Weight ");
        assert!(prims < attrs && attrs < weight, "{prop}: child order");
    }
    let disc = find_from(&out, 0, "<DiscMaterial");
    assert!(find_from(&out, disc, "<Weight ") < find_from(&out, disc, "<DiscFile"));
    let again = ContinuousStructure::from_xml_str(&out).unwrap();
    assert_eq!(again.to_xml_string(), out);
}
