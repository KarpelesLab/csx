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

```rust,ignore
let model = csx::ContinuousStructure::read_file("patch.xml")?;
let hit = model.property_at(csx::Vec3::new(0.0, 0.0, 0.5), csx::CoordSystem::Cartesian, csx::TypeMask::ANY);
```

See `cargo run --example info -- file.xml`.

## Features

| feature | default | description |
|---|---|---|
| `xml` | yes | read/write CSXCAD & openEMS XML files |

## License

MIT — see [LICENSE](LICENSE).
