//! Discretises a CSXCAD / openEMS XML file onto its mesh and prints
//! statistics.
//!
//! ```text
//! cargo run --release --example discretize -- model.xml
//! ```

use csx::ContinuousStructure;
use csx::discretize::{DiscretizeOptions, discretize};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("usage: discretize <file.xml>")?;
    let mut s = ContinuousStructure::read_file(&path)?;
    if let Err(e) = s.update() {
        eprintln!("warning: update failed ({e}); using values as read");
    }
    let t = std::time::Instant::now();
    let g = discretize(&s, &DiscretizeOptions::default())?;
    println!(
        "grid {:?} ({} nodes) in {:.2?}",
        g.dims,
        g.len(),
        t.elapsed()
    );
    for d in 0..3 {
        let pec = g.pec[d].iter().filter(|&&b| b).count();
        let sheet = g.sheet[d].iter().filter(|&&n| n != 0).count();
        let (lo, hi) = g.epsilon[d]
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| {
                (a.min(v), b.max(v))
            });
        println!("axis {d}: {pec} PEC edges, {sheet} sheet edges, eps in [{lo}, {hi}]");
    }
    for e in &g.excitations {
        println!("excitation {:?}: {} amplitudes", e.name, e.amplitudes.len());
    }
    println!("{} lumped element edges", g.lumped.len());
    for w in &g.warnings {
        println!("warning: {w}");
    }
    Ok(())
}
