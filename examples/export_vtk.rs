//! Exports the geometry of a CSXCAD / openEMS XML file as a VTK surface file
//! and one binary STL file per property.
//!
//! ```text
//! cargo run --example export_vtk -- model.xml out_dir
//! ```

use csx::io::export::{VtkOptions, export_vtk, stl_binary};
use csx::tessellate::{TessellationOptions, TriMesh, tessellate_structure};
use csx::{ContinuousStructure, TypeMask};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let input = args
        .next()
        .ok_or("usage: export_vtk <file.xml> [out_dir]")?;
    let out = std::path::PathBuf::from(args.next().unwrap_or_else(|| ".".into()));
    std::fs::create_dir_all(&out)?;

    let s = ContinuousStructure::read_file(&input)?;
    let vtk = out.join("geometry.vtk");
    std::fs::write(&vtk, export_vtk(&s, &VtkOptions::default()))?;
    println!("wrote {}", vtk.display());

    let parts = tessellate_structure(&s, TypeMask::ANY, &TessellationOptions::default());
    let mut per_prop: Vec<TriMesh> = vec![TriMesh::new(); s.properties().len()];
    for p in &parts {
        per_prop[p.property_index].append(&p.mesh);
    }
    for (k, (prop, mesh)) in s.properties().iter().zip(&per_prop).enumerate() {
        if mesh.is_empty() {
            continue;
        }
        let name: String = prop
            .name
            .chars()
            .map(|c| {
                if c.is_ascii_alphanumeric() || c == '-' {
                    c
                } else {
                    '_'
                }
            })
            .collect();
        let path = out.join(format!("{k:03}_{name}.stl"));
        std::fs::write(&path, stl_binary(mesh))?;
        println!(
            "wrote {} ({} triangles)",
            path.display(),
            mesh.triangles.len()
        );
    }
    Ok(())
}
