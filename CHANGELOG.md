# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.0.2](https://github.com/KarpelesLab/csx/compare/v0.0.1...v0.0.2) - 2026-10-18

### Other

- Document Axis-to-index conversion in the axis convention
- Document properties, structure and grid; opaque PropertyId; structured validation
- Document primitives: module guide, doctests on every shape, Errors/Panics sections
- Document mesh, discretize, tessellate, fdtd and io modules
- Document core modules with examples, Errors sections and conventions
- Add missing constructors for consistency across shapes and properties
- Add prelude and root re-exports; operator impls; idiomatic fdtd_settings
- Add #[must_use] to builders and value-producing functions
- Mark growable public types non_exhaustive; option builders; misc impls
- Make public API names consistent
- Remove panic paths from public functions; hide Bessel helpers
- Remove completed PLAN.md; move compatibility notes to README
- Use core EULER_GAMMA constant (newer clippy approx_constant lint)
- Record cylindrical follow-up in plan
- Honour cylindrical primitive coordinates everywhere; match CSXCAD query semantics
- Update open items in plan
- Match the XML layout CSXCAD writes and reads (verified black-box)
- Align geometry semantics with CSXCAD (verified black-box)
- Compile expressions to iterative postfix bytecode; remove MAX_DEPTH
- Encapsulate stale-prone caches; keep grid/background extra attributes
- Replace raw i32 type ids with typed enums
- Preserve unknown XML child elements on round-trip
- Make polygon vertices parametric
- Record open items before 1.0
- Add convenience constructors for properties, excitations, probes and points
- Preserve foreign property elements and accept expression vertices in XML
- Return exact material values when all averaged samples agree
- Mesh the full circle when a cylindrical model defines no angles
- Exact graded smoothing and resource limits on grid-sized work
- Treat NaN consistently: false in conditions, equal in Scalar comparisons
- Fix singular/unbounded transforms, Matrix convention and cylindrical alpha bounds
- Add patch antenna, VTK export and discretize examples; end-to-end tests
- Update plan status and README feature list
- Add discretisation of structures onto the rectilinear grid
- Add typed openEMS FDTD settings, excitation signals and boundary conditions
- Add primitive tessellation and STL/PLY/VTK export
- Store operator chains flat so long sums are not limited by MAX_DEPTH
- Fix panics on invalid polygon axes and huge sector angles; add primitive property tests
- Harden expression parser and PLY reader; add Bessel functions
- Add info example that summarises a CSXCAD/openEMS file
- Add XML support: minimal DOM parser/writer and CSXCAD/openEMS reader/writer
- Accelerate polyhedron and spatial index queries with a BVH
- Add automatic mesh-line generation

## [0.0.1](https://github.com/KarpelesLab/csx/compare/v0.0.0...v0.0.1) - 2026-10-18

### Other

- Document crate overview with a doctest; record status and compatibility notes
- Add properties, rectilinear grid with smoothing, and ContinuousStructure
- Add primitives (box, multibox, sphere, shells, cylinder, polygons, rotpoly, curve, wire, point, polyhedron, user-defined) and STL/PLY readers
- Add parameter set, parametric scalars and transformations
