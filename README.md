# csx

[![CI](https://github.com/KarpelesLab/csx/actions/workflows/ci.yml/badge.svg)](https://github.com/KarpelesLab/csx/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/csx.svg)](https://crates.io/crates/csx)
[![docs.rs](https://img.shields.io/docsrs/csx)](https://docs.rs/csx)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A pure Rust library for continuous-structure geometry and material modeling,
compatible with [CSXCAD](https://openems.de) and openEMS.

* Clean-room implementation (no CSXCAD/openEMS code was consulted), MIT licensed.
  Behaviour is derived from public documentation, the observable file
  format and black-box runs of CSXCAD.
* Zero runtime dependencies; XML support is an optional (default) feature.
* Edition 2024, MSRV 1.95.

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

## Compatibility notes

Behaviour confirmed from public documentation, published XML output and
black-box runs of CSXCAD 0.7 (its Python API, never its source):

* Coordinate systems: a `Box` in cylindrical coordinates is an annular
  sector; spheres, cylinders, their shells, curves, wires and points convert
  their defining points from `(ρ, α, z)`; multi-boxes, polygons, extruded and
  revolved polygons and polyhedra are always Cartesian. Multi-boxes ignore
  their own `CoordSystem` and compare raw coordinates in the structure's
  system (no angle wrapping). Primitives without `CoordSystem` inherit the
  structure's; structure queries take coordinates in the structure's system.
  A cylindrical box compares cylindrical queries as given (negative `ρ` is
  outside, `α` modulo 2π); every other shape folds `(-ρ, α)` through the
  Cartesian conversion. Verified with 140k random cylindrical-structure
  queries. Known divergence: CSXCAD clips cylindrical wires to an
  inconsistently computed bounding box; csx keeps the true geometry.
* Curves and points never contain a point.
* Equal priorities: the earlier property wins.
* Transformations are applied in list order; angles are radians in files;
  `Matrix` is row-major acting on column vectors (translation in the last
  column); vector scaling is written `Scale3`.
* Polygon vertices `(X1, X2)` map to axes `(n+1)%3, (n+2)%3` for normal `n`;
  revolved polygons rotate right-handed about `RotAxisDir` starting at the
  remaining (radial) axis.
* Material values are written as three comma separated values; with
  `Isotropy="1"` only the first is used. Dispersive terms use 1-based `_n`
  suffixes in `Property` and `Weight`.
* Weighting function variables: `x y z rho a r t`, in drawing units.
* Expressions: any non-zero value is true.
* `UserDefined` primitives are not supported by CSXCAD 0.7; csx reads and
  writes them in its own format (`<Equation>` text, `<CoordShift>` origin).

## Features

| feature | default | description |
|---|---|---|
| `xml` | yes | read/write CSXCAD & openEMS XML files |

## License

MIT — see [LICENSE](LICENSE).
