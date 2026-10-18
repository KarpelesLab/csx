# csx — implementation plan

`csx` is a pure-Rust, clean-room implementation of a *continuous structure*
geometry and material model, compatible with the CSXCAD data model and the
XML files consumed/produced by CSXCAD and openEMS.

## Ground rules

* **Clean room.** No CSXCAD / openEMS source code is read or copied. Behaviour
  is derived from public documentation, the published API reference and the
  observable XML file format (interoperability). Where the documentation is
  silent, behaviour is chosen deliberately and documented in rustdoc.
* **License:** MIT.
* **Toolchain:** edition 2024, MSRV 1.95 (`rust-version = "1.95"`).
* **Minimum dependencies:** zero runtime dependencies. The XML reader/writer
  is a small in-tree implementation, gated behind the optional `xml` cargo
  feature (enabled by default).
* `#![forbid(unsafe_code)]`, no panics on malformed input (errors instead).
* Versioning starts at `0.0.0`.

## Architecture (single crate, modules)

| Module | Responsibility |
|---|---|
| `math` | `Vec3`, 4×4 affine matrices, `BoundingBox`, numeric helpers |
| `coords` | `CoordSystem` (Cartesian / Cylindrical), conversions, angle wrapping |
| `expr` | fparser-compatible expression parser + evaluator (constant folding, variables) |
| `param` | `ParameterSet`, `Scalar` (constant or expression with cached value) |
| `transform` | Ordered transformation list (translate, rotate origin/x/y/z, scale, matrix) |
| `primitives` | Box, MultiBox, Sphere, SphericalShell, Cylinder, CylindricalShell, Polygon, LinPoly, RotPoly, Curve, Wire, Point, Polyhedron, PolyhedronReader, UserDefined |
| `properties` | Material, Lorentz/Debye dispersive materials, DiscMaterial (metadata), Metal, ConductingSheet, Excitation, ProbeBox, DumpBox, LumpedElement, AbsorbingBC, Unknown |
| `grid` | `RectilinearGrid`: lines, delta unit, lookup/snap, helpers for line generation/smoothing |
| `structure` | `ContinuousStructure`: owns properties, grid, parameters, background material; priority-based spatial queries; `update()` re-evaluation |
| `io::xml` (feature `xml`) | minimal XML DOM parser + writer; CSXCAD/openEMS reader & writer with round-trip of unknown attributes |
| `io::stl`, `io::ply` | STL (ascii/binary) and PLY (ascii/binary LE/BE) mesh readers for polyhedra |

### Key design decisions

* **Data model.** Properties own their primitives. Spatial queries return
  `(&Property, &Primitive)` pairs. Properties and primitives carry stable ids.
* **Values.** Every user-facing number that CSXCAD allows to be an expression
  is a `Scalar` = `{ expr: Option<String>, value: f64 }`. `structure.update()`
  re-evaluates all expressions against the `ParameterSet` (like CSXCAD's
  update step); hot paths only read the cached `f64`.
* **Coordinates.** All queries take a point *and* the coordinate system it is
  expressed in. Primitives may be defined in either system; a cylindrical
  Box is an annular sector (ρ, α, z ranges).
* **Inside test.** `point → (convert CS) → inverse transform → local test`.
  Polyhedra use a generalised winding number (robust for closed meshes) with
  bounding-box rejection; a BVH can be added later.
* **Priority.** `structure.property_at(point, cs, filter)` returns the
  property whose containing primitive has highest priority (ties: last
  inserted wins, documented).
* **Material evaluation.** `Material::epsilon_at(point)` etc. apply the
  weighting expressions (variables x, y, z, rho, r, a, t) per direction.
* **Round-trip.** Unknown XML attributes/elements are preserved in
  `extra` maps so read→write does not lose information.

## Milestones

1. **Scaffold** — Cargo.toml (0.0.0, MSRV, features), LICENSE, README, CI
   (fmt, clippy, test on stable + 1.95, `--no-default-features`).
2. **Core math & coordinates** — Vec3, Mat4, BoundingBox, cylindrical
   conversions, tests.
3. **Expressions & parameters** — tokenizer, Pratt parser, AST, evaluator,
   fparser function set, `ParameterSet`, `Scalar`; property tests.
4. **Transformations** — op list, composition, inverse, (de)serialisable
   argument form.
5. **Primitives** — all types: `is_inside`, `bounding_box`, `dimension`,
   update; extensive tests per primitive (incl. cylindrical CS, transforms).
6. **Properties** — all types, colours, visibility, attributes, material
   weighting, dispersive parameters.
7. **ContinuousStructure** — container API, ids, lookup by name/type,
   priority queries, background material, bounding box of the whole model.
8. **Rectilinear grid** — lines management, lookup, snapping, simple
   smooth-mesh generation helpers.
9. **XML I/O** (feature `xml`) — minimal XML parser/writer, CSXCAD reader
   and writer, openEMS document wrapper (FDTD section kept typed-or-raw),
   round-trip tests with hand-written fixtures.
10. **Mesh file readers** — STL/PLY for `PolyhedronReader`.
11. **Polish** — docs, examples, benches, fuzz-style robustness tests,
    performance (BVH for polyhedra, per-property bbox culling).

## Status

| Milestone | State |
|---|---|
| 1. Scaffold | done |
| 2. Core math & coordinates | done |
| 3. Expressions & parameters | done (incl. Bessel extensions, depth limits) |
| 4. Transformations | done |
| 5. Primitives | done (BVH-accelerated polyhedra) |
| 6. Properties | done |
| 7. ContinuousStructure | done (BVH spatial index) |
| 8. Rectilinear grid | done (graded smoothing, automatic mesh generation in `mesh`) |
| 9. XML I/O | done (round-trips 39 real-world files; typed `<FDTD>` settings in `fdtd`) |
| 10. Mesh file readers (STL/PLY) | done |
| 11. Polish | ongoing |

Beyond the original plan:

* `tessellate` + `io::export`: triangle meshes of every primitive, STL/PLY/VTK export.
* `discretize`: sampling materials, PEC/sheets, excitations and lumped elements
  onto the Yee grid.
* `fdtd`: typed openEMS simulation settings, excitation signals, boundaries.

### Compatibility notes

Facts confirmed from public documentation, published XML output and
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

## Open items before 1.0

All six items identified after the first implementation pass are resolved:

1. ✅ Polygon vertices are parametric (`PolyVertices`, cached evaluation).
2. ✅ Stale-prone caches are private behind invalidating setters
   (`UserDefined`, `PolyhedronReader`); the `update()` contract is documented
   on `ContinuousStructure`.
3. ✅ Type ids are enums with `Other(i32)` fallbacks (`ExcitationType`,
   `ProbeType`, `DumpType`, `DumpMode`, `DumpFileType`, `LumpedElementType`,
   `AbcType`, `TimeStepMethod`, `Axis`).
4. ✅ Unknown child elements and attributes survive XML round-trips
   (structure, properties, primitives, grid, background material).
5. ✅ Expressions compile to iterative bytecode; no nesting limit.
6. ✅ File-format guesses verified black-box against CSXCAD 0.7 (see
   compatibility notes). `UserDefined` is not supported by CSXCAD 0.7 and
   keeps this crate's own format.

Also done: cylindrical-coordinate primitives are honoured in `mesh`,
`discretize` and `tessellate`, and cylindrical query semantics match CSXCAD
(differentially tested, 0 mismatches over 140k random points).

Known divergence: CSXCAD clips wires defined in cylindrical coordinates to an
inconsistently computed bounding box; this crate uses the true geometry.
