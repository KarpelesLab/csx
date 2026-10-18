# csx

[![CI](https://github.com/KarpelesLab/csx/actions/workflows/ci.yml/badge.svg)](https://github.com/KarpelesLab/csx/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/csx.svg)](https://crates.io/crates/csx)
[![docs.rs](https://img.shields.io/docsrs/csx)](https://docs.rs/csx)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A pure Rust library for continuous-structure geometry and material modeling,
compatible with [CSXCAD](https://openems.de) and openEMS.

* Clean-room implementation (no CSXCAD/openEMS code was consulted), MIT licensed.
* Zero runtime dependencies; XML support is an optional (default) feature.
* Edition 2024, MSRV 1.95.

Status: early development — see [PLAN.md](PLAN.md).

## Overview

* Geometric primitives: box, multi-box, sphere, spherical shell, cylinder,
  cylindrical shell, polygon, extruded polygon, revolved polygon, curve, wire,
  point, polyhedron (inline or STL/PLY file) and expression-defined regions,
  with transformations and Cartesian/cylindrical coordinates.
* Properties: materials (anisotropic, position-weighted), Drude/Lorentz and
  Debye dispersive materials, metals, conducting sheets, excitations,
  probes, field dumps, lumped elements.
* Parametric models: every value may be an expression over named parameters.
* Priority-based point queries (`property_at`, `medium_at`) and spatial index.
* Rectilinear grids with graded smoothing and automatic mesh generation from
  the geometry (`mesh`).
* Read/write CSXCAD and openEMS XML files, including typed FDTD settings
  (feature `xml`).
* Tessellation of primitives and STL/PLY/VTK export; STL/PLY import.
* Discretisation of materials, metals, excitations and lumped elements onto
  the FDTD grid (`discretize`).

## Example

```rust
use csx::primitives::BoxShape;
use csx::properties::Material;
use csx::{ContinuousStructure, CoordSystem, Property, TypeMask, Vec3};

let mut model = ContinuousStructure::new();
let fr4 = model.add_property(Property::material("FR4", Material::new(4.3, 0.0)));
model.add_primitive(fr4, BoxShape::new([-10.0, -10.0, 0.0], [10.0, 10.0, 1.6])).unwrap();

// Write and read back a CSXCAD XML file (feature `xml`).
let xml = model.to_xml_string();
let model = ContinuousStructure::from_xml_str(&xml).unwrap();

let hit = model.property_at(Vec3::new(0.0, 0.0, 0.5), CoordSystem::Cartesian, TypeMask::ANY);
assert_eq!(hit.unwrap().property.name, "FR4");
```

More examples:

* `cargo run --example info -- file.xml` — summary of a CSXCAD/openEMS file
* `cargo run --example patch_antenna -- out.xml [out.vtk]` — build, mesh and
  write a parametric patch antenna
* `cargo run --example export_vtk -- file.xml out_dir` — VTK + STL export
* `cargo run --release --example discretize -- file.xml` — discretise onto
  the mesh and print statistics

## Features

| feature | default | description |
|---|---|---|
| `xml` | yes | read/write CSXCAD & openEMS XML files |

## License

MIT — see [LICENSE](LICENSE).
