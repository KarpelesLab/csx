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

Facts confirmed from public documentation / published XML output that shaped
the implementation:

* Only `Box` (and `MultiBox`) interpret their coordinates in the primitive's
  coordinate system (cylindrical boxes are annular sectors); all other
  primitives are always Cartesian.
* Transformations are applied in list order; rotation angles are radians in
  files.
* Polygon vertices `(X1, X2)` map to axes `(n+1)%3, (n+2)%3` for normal `n`.
* Material values are written as one or three comma separated values; with
  `Isotropy="1"` only the first is used. Debye terms use 1-based `_n`
  suffixes.
* Weighting function variables: `x y z rho a r t`, in drawing units.
* Expressions: any non-zero value is true.

## Open items before 1.0

1. Polygon vertices are plain `f64`; parametric vertices are evaluated on read
   and lost on write. Consider `Vec<[Scalar; 2]>`.
2. Public fields that can go stale (`UserDefined::function`, cached
   `PolyhedronReader::mesh`, ...): make private with setters, or keep
   `update()` as the documented contract everywhere.
3. Raw `i32` type ids (excitation/probe/dump types, lumped direction) →
   enums with an `Other(i32)` fallback.
4. Unknown child elements of properties/primitives are dropped on XML
   round-trip (attributes are preserved).
5. `MAX_DEPTH` still bounds nested parentheses / unary chains; iterative
   evaluation would remove the limit.
6. Unconfirmed file-format guesses (Matrix convention, RotPoly, Polyhedron,
   ParameterSet, UserDefined attribute names) need verification against
   genuine CSXCAD output.
