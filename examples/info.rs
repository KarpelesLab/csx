//! Prints a summary of a CSXCAD / openEMS XML file.
//!
//! ```text
//! cargo run --example info -- model.xml
//! ```

use csx::{ContinuousStructure, CoordSystem};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args().nth(1).ok_or("usage: info <file.xml>")?;
    let s = ContinuousStructure::read_file(&path)?;
    let cs = s.coord_system;
    println!("coordinate system: {cs:?}");
    let g = &s.grid;
    println!(
        "grid: {:?} lines ({} cells), delta unit {} m",
        g.num_lines(),
        g.num_cells(),
        g.delta_unit
    );
    for p in s.properties() {
        println!(
            "{:>16} {:<24} {} primitive(s)",
            p.kind.type_name(),
            p.name,
            p.primitives.len()
        );
    }
    let b = s.bounding_box(CoordSystem::Cartesian);
    println!("bounding box: {:?} .. {:?}", b.min.0, b.max.0);
    for w in s.validate() {
        println!("warning: {w}");
    }
    Ok(())
}
