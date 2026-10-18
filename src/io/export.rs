//! Mesh export: STL, PLY and legacy VTK.

use core::fmt::Write as _;

use crate::math::Vec3;
use crate::properties::TypeMask;
use crate::structure::ContinuousStructure;
use crate::tessellate::{TessellationOptions, TriMesh, tessellate_structure};

/// Binary STL of a mesh.
pub fn stl_binary(mesh: &TriMesh) -> Vec<u8> {
    super::stl::write_binary(&mesh.triangle_soup())
}

/// ASCII STL of a mesh; `name` is written after `solid`.
pub fn stl_ascii(mesh: &TriMesh, name: &str) -> String {
    let name: String = name.chars().filter(|c| !c.is_control()).collect();
    let mut s = format!("solid {name}\n");
    for [a, b, c] in mesh.triangle_soup() {
        let n = (b - a).cross(c - a).normalized();
        let _ = writeln!(s, "  facet normal {} {} {}", n[0], n[1], n[2]);
        s.push_str("    outer loop\n");
        for v in [a, b, c] {
            let _ = writeln!(s, "      vertex {} {} {}", v[0], v[1], v[2]);
        }
        s.push_str("    endloop\n  endfacet\n");
    }
    let _ = writeln!(s, "endsolid {name}");
    s
}

/// ASCII PLY of a mesh.
pub fn ply_ascii(mesh: &TriMesh) -> String {
    let mut s = String::new();
    let _ = write!(
        s,
        "ply\nformat ascii 1.0\ncomment written by csx\nelement vertex {}\nproperty double x\nproperty double y\nproperty double z\nelement face {}\nproperty list uchar uint vertex_indices\nend_header\n",
        mesh.vertices.len(),
        mesh.triangles.len()
    );
    for v in &mesh.vertices {
        let _ = writeln!(s, "{} {} {}", v[0], v[1], v[2]);
    }
    for t in &mesh.triangles {
        let _ = writeln!(s, "3 {} {} {}", t[0], t[1], t[2]);
    }
    s
}

/// Options for [`export_vtk`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VtkOptions {
    /// Tessellation resolution.
    pub tessellation: TessellationOptions,
    /// Which properties to export.
    pub mask: TypeMask,
    /// Scale applied to coordinates (e.g. the grid's delta unit to export in
    /// metres).
    pub scale: f64,
}

impl Default for VtkOptions {
    fn default() -> Self {
        VtkOptions {
            tessellation: TessellationOptions::default(),
            mask: TypeMask::ANY,
            scale: 1.0,
        }
    }
}

/// Exports the tessellated geometry of a structure as a legacy ASCII VTK
/// `POLYDATA` file. Surfaces become `POLYGONS`, curves and wires `LINES`.
/// Cell data: `property_index` (index into
/// [`ContinuousStructure::properties`]) and `priority`; field data lists the
/// property names (whitespace replaced by `_`) in index order.
pub fn export_vtk(s: &ContinuousStructure, opts: &VtkOptions) -> String {
    let parts = tessellate_structure(s, opts.mask, &opts.tessellation);
    let mut points: Vec<Vec3> = Vec::new();
    let mut lines: Vec<(Vec<u32>, usize, i32)> = Vec::new();
    let mut polys: Vec<([u32; 3], usize, i32)> = Vec::new();
    for p in &parts {
        let off = points.len() as u32;
        points.extend(p.mesh.vertices.iter().map(|v| *v * opts.scale));
        polys.extend(
            p.mesh
                .triangles
                .iter()
                .map(|t| (t.map(|i| i + off), p.property_index, p.priority)),
        );
        if p.polyline.len() >= 2 {
            let off = points.len() as u32;
            points.extend(p.polyline.iter().map(|v| *v * opts.scale));
            lines.push((
                (0..p.polyline.len() as u32).map(|i| i + off).collect(),
                p.property_index,
                p.priority,
            ));
        }
    }
    let mut out =
        String::from("# vtk DataFile Version 3.0\ncsx structure\nASCII\nDATASET POLYDATA\n");
    let _ = writeln!(out, "POINTS {} double", points.len());
    for v in &points {
        let _ = writeln!(out, "{} {} {}", v[0], v[1], v[2]);
    }
    if !lines.is_empty() {
        let size: usize = lines.iter().map(|l| l.0.len() + 1).sum();
        let _ = writeln!(out, "LINES {} {}", lines.len(), size);
        for (l, _, _) in &lines {
            let _ = write!(out, "{}", l.len());
            for i in l {
                let _ = write!(out, " {i}");
            }
            out.push('\n');
        }
    }
    if !polys.is_empty() {
        let _ = writeln!(out, "POLYGONS {} {}", polys.len(), polys.len() * 4);
        for (t, _, _) in &polys {
            let _ = writeln!(out, "3 {} {} {}", t[0], t[1], t[2]);
        }
    }
    // Cell order in legacy POLYDATA: vertices, lines, polygons, strips.
    let cells: Vec<(usize, i32)> = lines
        .iter()
        .map(|l| (l.1, l.2))
        .chain(polys.iter().map(|p| (p.1, p.2)))
        .collect();
    let _ = writeln!(out, "CELL_DATA {}", cells.len());
    if !cells.is_empty() {
        out.push_str("SCALARS property_index int 1\nLOOKUP_TABLE default\n");
        for (i, _) in &cells {
            let _ = writeln!(out, "{i}");
        }
        out.push_str("SCALARS priority int 1\nLOOKUP_TABLE default\n");
        for (_, p) in &cells {
            let _ = writeln!(out, "{p}");
        }
    }
    let names: Vec<String> = s
        .properties()
        .iter()
        .map(|p| {
            let n: String = p
                .name
                .chars()
                .map(|c| {
                    if c.is_whitespace() || c.is_control() {
                        '_'
                    } else {
                        c
                    }
                })
                .collect();
            if n.is_empty() { "_".to_owned() } else { n }
        })
        .collect();
    if !names.is_empty() {
        let _ = writeln!(
            out,
            "FIELD property_names 1\nnames {} 1 string",
            names.len()
        );
        for n in &names {
            let _ = writeln!(out, "{n}");
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{BoxShape, Curve, Sphere};
    use crate::{CoordSystem, Primitive, Property};

    fn cube() -> TriMesh {
        crate::tessellate::tessellate(
            &Primitive::new(BoxShape::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])),
            CoordSystem::Cartesian,
            &TessellationOptions::default(),
        )
    }

    #[test]
    fn stl_roundtrip() {
        let m = cube();
        let a = super::super::stl::parse(stl_ascii(&m, "cube").as_bytes()).unwrap();
        assert!(a.is_closed());
        assert_eq!(a.vertices().len(), 8);
        let b = super::super::stl::parse(&stl_binary(&m)).unwrap();
        assert!(b.contains(Vec3::new(0.5, 0.5, 0.5)));
    }

    #[test]
    fn ply_roundtrip() {
        let m = cube();
        let p = super::super::ply::parse(ply_ascii(&m).as_bytes()).unwrap();
        assert_eq!(p.vertices(), &m.vertices[..]);
        assert!(p.is_closed());
    }

    #[test]
    fn vtk() {
        let mut s = ContinuousStructure::new();
        let a = s.add_property(Property::metal("pec plate"));
        s.add_primitive(a, BoxShape::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]))
            .unwrap();
        let b = s.add_property(Property::metal("wire"));
        s.add_primitive(
            b,
            Curve::new([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]),
        )
        .unwrap();
        s.add_primitive(
            b,
            Primitive::new(Sphere::new([3.0, 0.0, 0.0], 1.0)).with_priority(4),
        )
        .unwrap();
        let v = export_vtk(
            &s,
            &VtkOptions {
                scale: 1e-3,
                ..Default::default()
            },
        );
        assert!(v.starts_with("# vtk DataFile Version 3.0\n"));
        assert!(v.contains("LINES 1 4\n"));
        assert!(v.contains("POLYGONS "));
        let lines: Vec<&str> = v.lines().collect();
        let cd = lines
            .iter()
            .position(|l| l.starts_with("CELL_DATA"))
            .unwrap();
        let n: usize = lines[cd]["CELL_DATA ".len()..].parse().unwrap();
        // first cell is the line (property 1), then the cube (property 0)
        assert_eq!(lines[cd + 3], "1");
        assert_eq!(lines[cd + 4], "0");
        assert!(n > 13);
        assert!(v.contains("pec_plate"));
        // empty structure still valid
        assert!(
            export_vtk(&ContinuousStructure::new(), &VtkOptions::default())
                .contains("POINTS 0 double")
        );
    }
}
