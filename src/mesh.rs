//! Automatic mesh-line generation.
//!
//! Mesh generation works in two steps:
//!
//! 1. [`detect_edges`] collects *fixed* lines from the geometry: box faces,
//!    polygon vertices and elevations, curve vertices, the extents of round
//!    and polyhedral shapes, ... Primitives with a rotating transformation
//!    contribute the extents of their bounding box only.
//! 2. [`MeshBuilder`] merges lines that are closer than a minimum
//!    resolution, optionally pads the model with an air margin and fills the
//!    gaps with [`smooth_lines`] so that no cell exceeds the maximum
//!    resolution and neighbouring cells grow by at most a given ratio.
//!
//! Optionally the *thirds rule* is applied to metal edges: instead of a line
//! exactly on the edge of a metal sheet (where the field is singular), one
//! line is placed `d/3` inside the metal and one `2d/3` outside, where `d` is
//! the local resolution.

use crate::coords::CoordSystem;
use crate::error::{Error, Result};
use crate::grid::{RectilinearGrid, smooth_lines};
use crate::math::{BoundingBox, Mat4, Vec3};
use crate::primitives::{Primitive, Shape, plane_axes};
use crate::properties::{PropertyType, TypeMask};
use crate::structure::ContinuousStructure;

/// Options of the automatic mesher.
#[derive(Debug, Clone, PartialEq)]
pub struct MeshOptions {
    /// Largest allowed cell size per axis (drawing units; radians for the
    /// angular axis of a cylindrical mesh).
    pub max_res: [f64; 3],
    /// Largest allowed size ratio of neighbouring cells (≥ 1).
    pub ratio: f64,
    /// Lines closer than this (per axis) are merged before smoothing.
    pub min_res: Option<[f64; 3]>,
    /// Enables the thirds rule for metal edges with the given resolution
    /// per axis (the distance `d` between the two lines around an edge is
    /// `d/3 + 2d/3 = d`).
    pub thirds: Option<[f64; 3]>,
    /// Padding `[below, above]` added around the detected lines per axis
    /// (air box). The radial axis of a cylindrical mesh is clamped at 0.
    pub margin: [[f64; 2]; 3],
    /// Properties whose primitives contribute lines.
    pub mask: TypeMask,
}

impl MeshOptions {
    /// Options with the same maximum resolution on every axis, a growth
    /// ratio of 1.5, no thirds rule and no margin.
    pub fn new(max_res: f64) -> Self {
        MeshOptions {
            max_res: [max_res; 3],
            ratio: 1.5,
            min_res: None,
            thirds: None,
            margin: [[0.0; 2]; 3],
            mask: TypeMask::ANY,
        }
    }

    /// Builder: per-axis maximum resolution.
    pub fn with_max_res(mut self, max_res: [f64; 3]) -> Self {
        self.max_res = max_res;
        self
    }

    /// Builder: growth ratio.
    pub fn with_ratio(mut self, ratio: f64) -> Self {
        self.ratio = ratio;
        self
    }

    /// Builder: minimum line distance (merging threshold).
    pub fn with_min_res(mut self, min_res: [f64; 3]) -> Self {
        self.min_res = Some(min_res);
        self
    }

    /// Builder: enables the thirds rule with resolution `d` per axis.
    pub fn with_thirds(mut self, d: [f64; 3]) -> Self {
        self.thirds = Some(d);
        self
    }

    /// Builder: the same margin on both sides of every axis.
    pub fn with_margin(mut self, m: f64) -> Self {
        self.margin = [[m; 2]; 3];
        self
    }

    /// Builder: property filter.
    pub fn with_mask(mut self, mask: impl Into<TypeMask>) -> Self {
        self.mask = mask.into();
        self
    }
}

/// Collects candidate mesh lines from all primitives of properties matching
/// `mask`, expressed in the coordinate system `cs` (`ρ, α, z` lines for
/// cylindrical meshes). The result is sorted and de-duplicated.
pub fn detect_edges(
    s: &ContinuousStructure,
    mask: impl Into<TypeMask>,
    cs: CoordSystem,
) -> [Vec<f64>; 3] {
    let mask = mask.into();
    let mut lines: [Vec<f64>; 3] = Default::default();
    for (prop, prim) in s.primitives() {
        if mask.contains(prop.property_type()) {
            primitive_lines(prim, s.coord_system, cs, &mut lines);
        }
    }
    lines.iter_mut().for_each(normalize);
    lines
}

fn normalize(v: &mut Vec<f64>) {
    v.retain(|x| x.is_finite());
    v.sort_by(f64::total_cmp);
    let scale = v.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    let eps = scale * 1e-13;
    v.dedup_by(|b, a| (*b - *a).abs() <= eps);
}

/// True if the linear part of `m` keeps the coordinate axes aligned
/// (translation and per-axis scaling only).
fn axis_aligned(m: &Mat4) -> bool {
    (0..3).all(|i| (0..3).all(|j| i == j || m.0[i][j] == 0.0))
}

/// Adds the coordinates of `pts` (local, primitive frame) to `lines`.
fn add_points(prim: &Primitive, cs: CoordSystem, pts: &[Vec3], lines: &mut [Vec<f64>; 3]) {
    for &p in pts {
        let w = prim.transform.apply(p);
        let q = cs.from_cartesian(w);
        for d in 0..3 {
            if cs == CoordSystem::Cylindrical && d == 1 && q[0] == 0.0 {
                continue; // the angle is meaningless on the axis
            }
            lines[d].push(q[d]);
        }
    }
}

fn add_bbox(b: &BoundingBox, cs: CoordSystem, lines: &mut [Vec<f64>; 3]) {
    if b.is_empty() {
        return;
    }
    for (d, l) in lines.iter_mut().enumerate() {
        if cs == CoordSystem::Cylindrical
            && d == 1
            && b.max[1] - b.min[1] >= core::f64::consts::TAU - 1e-12
        {
            continue; // full revolution: no angular edges
        }
        l.push(b.min[d]);
        l.push(b.max[d]);
    }
}

fn primitive_lines(
    prim: &Primitive,
    default_cs: CoordSystem,
    cs: CoordSystem,
    lines: &mut [Vec<f64>; 3],
) {
    let pcs = prim.effective_cs(default_cs);
    let aligned = prim.transform.is_empty() || axis_aligned(prim.transform.matrix());
    // Exact bounds in the mesh coordinate system always contribute.
    add_bbox(&prim.bounding_box_in(cs, default_cs), cs, lines);
    // Boxes and multi-boxes whose values are not Cartesian (sectors, ρ/α/z
    // ranges) are fully described by their bounds; everything else may
    // contribute interior lines when the mapping keeps axes aligned.
    // Point-defined shapes are converted to Cartesian first.
    if !aligned
        || (prim.shape_cs(default_cs) != CoordSystem::Cartesian
            && matches!(prim.shape, Shape::Box(_) | Shape::MultiBox(_)))
    {
        return;
    }
    let shape = prim.shape.to_cartesian(pcs);
    let mut pts: Vec<Vec3> = Vec::new();
    match &*shape {
        Shape::Box(b) => pts.extend([b.start.value(), b.stop.value()]),
        Shape::MultiBox(m) => m
            .boxes
            .iter()
            .for_each(|(a, b)| pts.extend([a.value(), b.value()])),
        Shape::Polygon(p) => {
            let n = p.norm_dir;
            let (u, v) = plane_axes(n);
            for &(a, c) in &p.vertices {
                let mut q = Vec3::ZERO;
                q[u] = a;
                q[v] = c;
                q[n] = p.elevation.value();
                pts.push(q);
            }
        }
        Shape::LinPoly(p) => {
            let n = p.norm_dir;
            let (u, v) = plane_axes(n);
            for &(a, c) in &p.vertices {
                for h in [p.elevation.value(), p.elevation.value() + p.length.value()] {
                    let mut q = Vec3::ZERO;
                    q[u] = a;
                    q[v] = c;
                    q[n] = h;
                    pts.push(q);
                }
            }
        }
        Shape::Curve(c) => pts.extend(c.points.iter().map(|p| p.value())),
        Shape::Wire(w) => pts.extend(w.points.iter().map(|p| p.value())),
        Shape::Point(p) => pts.push(p.point.value()),
        Shape::Sphere(s) => pts.push(s.center.value()),
        Shape::SphericalShell(s) => pts.push(s.center.value()),
        Shape::Cylinder(c) => pts.extend([c.start.value(), c.stop.value()]),
        Shape::CylindricalShell(c) => pts.extend([c.start.value(), c.stop.value()]),
        Shape::RotPoly(_)
        | Shape::Polyhedron(_)
        | Shape::PolyhedronReader(_)
        | Shape::UserDefined(_) => {}
    }
    add_points(prim, cs, &pts, lines);
}

/// Lines produced by the thirds rule for metal primitives, together with
/// the edge coordinates they replace.
fn thirds_lines(
    s: &ContinuousStructure,
    mask: TypeMask,
    cs: CoordSystem,
    d: [f64; 3],
) -> ([Vec<f64>; 3], [Vec<f64>; 3]) {
    let mut add: [Vec<f64>; 3] = Default::default();
    let mut remove: [Vec<f64>; 3] = Default::default();
    for (prop, prim) in s.primitives() {
        let t = prop.property_type();
        let metal = TypeMask::of(PropertyType::Metal) | TypeMask::of(PropertyType::ConductingSheet);
        if !mask.contains(t) || !metal.contains(t) || !matches!(prim.shape, Shape::Box(_)) {
            continue;
        }
        let b = prim.bounding_box_in(cs, s.coord_system);
        if b.is_empty() {
            continue;
        }
        for dir in 0..3 {
            let (lo, hi) = (b.min[dir], b.max[dir]);
            let dd = d[dir];
            // Only edges of extents wide enough to hold both inner lines.
            if dd.is_nan() || dd <= 0.0 || hi - lo <= 2.0 * dd / 3.0 {
                continue;
            }
            add[dir].extend([
                lo + dd / 3.0,
                lo - 2.0 * dd / 3.0,
                hi - dd / 3.0,
                hi + 2.0 * dd / 3.0,
            ]);
            remove[dir].extend([lo, hi]);
        }
    }
    (add, remove)
}

/// Merges lines closer than `min` (keeping the first of each cluster, but
/// preferring lines in `keep`).
fn merge_close(v: &[f64], min: f64, keep: &[f64]) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::with_capacity(v.len());
    for &x in v {
        match out.last_mut() {
            Some(last) if x - *last < min => {
                if keep.contains(&x) && !keep.contains(last) {
                    *last = x;
                }
            }
            _ => out.push(x),
        }
    }
    out
}

/// Builds a [`RectilinearGrid`] for a structure.
#[derive(Debug, Clone, PartialEq)]
pub struct MeshBuilder {
    /// Options.
    pub options: MeshOptions,
    extra: [Vec<f64>; 3],
}

impl MeshBuilder {
    /// Builder with the given options.
    pub fn new(options: MeshOptions) -> Self {
        MeshBuilder {
            options,
            extra: Default::default(),
        }
    }

    /// Adds user-defined fixed lines along `dir`.
    pub fn with_lines(mut self, dir: usize, lines: impl IntoIterator<Item = f64>) -> Self {
        self.extra[dir].extend(lines);
        self
    }

    /// Fixed lines (before smoothing) for the structure, in the structure's
    /// coordinate system.
    pub fn fixed_lines(&self, s: &ContinuousStructure) -> Result<[Vec<f64>; 3]> {
        let o = &self.options;
        let cs = s.coord_system;
        let mut lines = detect_edges(s, o.mask, cs);
        if let Some(d) = o.thirds {
            let (add, remove) = thirds_lines(s, o.mask, cs, d);
            // An edge is removed unless a material boundary needs it too.
            let others = detect_edges(s, TypeMask(o.mask.0 & TypeMask::MATERIALS.0), cs);
            for dir in 0..3 {
                lines[dir].retain(|x| !remove[dir].contains(x) || others[dir].contains(x));
                lines[dir].extend(&add[dir]);
            }
        }
        // A cylindrical model whose primitives all span full revolutions
        // defines no angular lines; mesh the full circle then.
        if cs == CoordSystem::Cylindrical && lines[1].len() < 2 && self.extra[1].len() < 2 {
            lines[1].extend([-core::f64::consts::PI, core::f64::consts::PI]);
        }
        for dir in 0..3 {
            lines[dir].extend(self.extra[dir].iter().copied());
            normalize(&mut lines[dir]);
            if let (Some(&lo), Some(&hi)) = (lines[dir].first(), lines[dir].last()) {
                let [below, above] = o.margin[dir];
                let angular = cs == CoordSystem::Cylindrical && dir == 1;
                if below > 0.0 {
                    let mut l = lo - below;
                    if cs == CoordSystem::Cylindrical && dir == 0 {
                        l = l.max(0.0);
                    }
                    if angular {
                        // never exceed a full turn
                        l = l.max(hi - core::f64::consts::TAU);
                    }
                    if l < lo {
                        lines[dir].insert(0, l);
                    }
                }
                if above > 0.0 {
                    let mut h = hi + above;
                    if angular {
                        h = h.min(lines[dir][0] + core::f64::consts::TAU);
                    }
                    if h > hi {
                        lines[dir].push(h);
                    }
                }
            }
            if let Some(min) = o.min_res {
                let keep = self.extra[dir].clone();
                lines[dir] = merge_close(&lines[dir], min[dir], &keep);
            }
        }
        Ok(lines)
    }

    /// Generates the grid: fixed lines, merged and padded, then smoothed.
    /// The grid takes the structure's coordinate system and delta unit.
    pub fn build(&self, s: &ContinuousStructure) -> Result<RectilinearGrid> {
        let o = &self.options;
        if o.max_res.iter().any(|r| r.is_nan() || *r <= 0.0) {
            return Err(Error::Invalid(format!(
                "mesh: invalid max_res {:?}",
                o.max_res
            )));
        }
        let fixed = self.fixed_lines(s)?;
        let mut g = RectilinearGrid::new();
        g.delta_unit = s.grid.delta_unit;
        g.coord_system = s.coord_system;
        for (dir, f) in fixed.iter().enumerate() {
            g.set_lines(dir, smooth_lines(f, o.max_res[dir], o.ratio)?)?;
        }
        Ok(g)
    }

    /// Builds the grid and stores it in the structure.
    pub fn apply(&self, s: &mut ContinuousStructure) -> Result<()> {
        s.grid = self.build(s)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{BoxShape, Cylinder, Polygon, Sphere, Wire};
    use crate::properties::{Excitation, Material, Property, PropertyKind};
    use core::f64::consts::PI;

    fn microstrip() -> ContinuousStructure {
        let mut s = ContinuousStructure::new();
        let sub = s.add_property(Property::material("sub", Material::new(3.66, 0.0)));
        s.add_primitive(sub, BoxShape::new([-20.0, -10.0, 0.0], [20.0, 10.0, 0.5]))
            .unwrap();
        let m = s.add_property(Property::metal("strip"));
        s.add_primitive(
            m,
            Primitive::new(BoxShape::new([-20.0, -0.5, 0.5], [20.0, 0.5, 0.5])).with_priority(10),
        )
        .unwrap();
        let p = s.add_property(Property::new(
            "port",
            PropertyKind::Excitation(Excitation::default()),
        ));
        s.add_primitive(p, BoxShape::new([-20.0, -0.5, 0.0], [-20.0, 0.5, 0.5]))
            .unwrap();
        s
    }

    fn has(v: &[f64], x: f64) -> bool {
        v.iter().any(|&y| (y - x).abs() < 1e-9)
    }

    #[test]
    fn edges() {
        let s = microstrip();
        let e = detect_edges(&s, TypeMask::ANY, CoordSystem::Cartesian);
        assert_eq!(e[0], vec![-20.0, 20.0]);
        assert_eq!(e[1], vec![-10.0, -0.5, 0.5, 10.0]);
        assert_eq!(e[2], vec![0.0, 0.5]);
        // property masks
        let m = detect_edges(&s, PropertyType::Metal, CoordSystem::Cartesian);
        assert_eq!(m[1], vec![-0.5, 0.5]);
        assert_eq!(m[2], vec![0.5]);
        let x = detect_edges(&s, TypeMask::MATERIALS, CoordSystem::Cartesian);
        assert_eq!(x[1], vec![-10.0, 10.0]);
    }

    #[test]
    fn build_grid() {
        let s = microstrip();
        let g = MeshBuilder::new(MeshOptions::new(1.0).with_ratio(1.4).with_margin(5.0))
            .build(&s)
            .unwrap();
        for dir in 0..3 {
            let l = g.lines(dir);
            let (_, max) = g.cell_size_range(dir).unwrap();
            assert!(max <= 1.0 + 1e-9, "dir {dir}: {max}");
            for w in l.windows(3) {
                let (a, b) = (w[1] - w[0], w[2] - w[1]);
                assert!(
                    a.max(b) / a.min(b) <= 1.4 * 1.05,
                    "dir {dir} ratio {}",
                    a.max(b) / a.min(b)
                );
            }
        }
        assert!(has(g.lines(0), -25.0) && has(g.lines(0), 25.0));
        assert!(has(g.lines(1), -0.5) && has(g.lines(1), 0.5));
        assert!(has(g.lines(2), 0.5) && has(g.lines(2), -5.0));
    }

    #[test]
    fn thirds_rule() {
        let s = microstrip();
        let b = MeshBuilder::new(MeshOptions::new(1.0).with_thirds([0.3; 3]));
        let f = b.fixed_lines(&s).unwrap();
        // strip edges in y replaced by thirds lines
        assert!(!has(&f[1], 0.5) && !has(&f[1], -0.5));
        assert!(has(&f[1], 0.4) && has(&f[1], 0.7));
        assert!(has(&f[1], -0.4) && has(&f[1], -0.7));
        // x edges are shared with the substrate and stay
        assert!(has(&f[0], 20.0) && has(&f[0], -20.0));
        // the zero-thickness z extent is not subject to the rule
        assert!(has(&f[2], 0.5));
    }

    #[test]
    fn merging_and_user_lines() {
        let s = microstrip();
        let b =
            MeshBuilder::new(MeshOptions::new(1.0).with_min_res([0.05; 3])).with_lines(1, [0.52]);
        let f = b.fixed_lines(&s).unwrap();
        assert!(has(&f[1], 0.52) && !has(&f[1], 0.5));
    }

    #[test]
    fn cylindrical() {
        let mut s = ContinuousStructure::new();
        s.coord_system = CoordSystem::Cylindrical;
        let m = s.add_property(Property::metal("inner"));
        s.add_primitive(m, Cylinder::new([0.0, 0.0, 0.0], [0.0, 0.0, 10.0], 1.0))
            .unwrap();
        let d = s.add_property(Property::material("diel", Material::new(2.0, 0.0)));
        s.add_primitive(
            d,
            BoxShape::new([1.0, -PI / 4.0, 0.0], [3.0, PI / 4.0, 10.0]),
        )
        .unwrap();
        let e = detect_edges(&s, TypeMask::ANY, CoordSystem::Cylindrical);
        assert!(has(&e[0], 0.0) && has(&e[0], 1.0) && has(&e[0], 3.0));
        assert!(has(&e[1], -PI / 4.0) && has(&e[1], PI / 4.0));
        assert!(has(&e[2], 0.0) && has(&e[2], 10.0));
        let g = MeshBuilder::new(
            MeshOptions::new(0.5)
                .with_max_res([0.5, 0.1, 1.0])
                .with_margin(2.0),
        )
        .build(&s)
        .unwrap();
        assert_eq!(g.coord_system, CoordSystem::Cylindrical);
        assert_eq!(g.lines(0)[0], 0.0); // clamped radial margin
        assert!(g.cell_size_range(1).unwrap().1 <= 0.1 + 1e-9);
    }

    #[test]
    fn cylindrical_full_revolution() {
        // Only full-revolution primitives: no angular lines are detected,
        // the mesher falls back to the full circle (regression: the alpha
        // axis was left empty and discretize failed).
        let mut s = ContinuousStructure::new();
        s.coord_system = CoordSystem::Cylindrical;
        let m = s.add_property(Property::metal("inner"));
        s.add_primitive(m, Cylinder::new([0.0, 0.0, 0.0], [0.0, 0.0, 10.0], 1.0))
            .unwrap();
        let g = MeshBuilder::new(
            MeshOptions::new(0.5)
                .with_max_res([0.5, 0.2, 1.0])
                .with_margin(1.0),
        )
        .build(&s)
        .unwrap();
        let a = g.lines(1);
        assert!(a.len() > 2);
        assert!((a[0] + PI).abs() < 1e-12 && (a[a.len() - 1] - PI).abs() < 1e-12);
        s.grid = g;
        assert!(crate::discretize::discretize(&s, &Default::default()).is_ok());
    }

    #[test]
    fn polygon_and_rotation() {
        let mut s = ContinuousStructure::new();
        let m = s.add_property(Property::metal("tri"));
        s.add_primitive(
            m,
            Polygon::new(2, 1.0, vec![(0.0, 0.0), (4.0, 0.0), (2.0, 3.0)]),
        )
        .unwrap();
        let e = detect_edges(&s, TypeMask::ANY, CoordSystem::Cartesian);
        assert_eq!(e[0], vec![0.0, 2.0, 4.0]);
        assert_eq!(e[1], vec![0.0, 3.0]);
        let mut t = crate::Transform::new();
        t.rotate_axis(2, PI / 4.0);
        let r = s.add_property(Property::metal("rot"));
        s.add_primitive(
            r,
            Primitive::new(BoxShape::new([0.0, 0.0, 5.0], [1.0, 1.0, 6.0])).with_transform(t),
        )
        .unwrap();
        let e = detect_edges(&s, PropertyType::Metal, CoordSystem::Cartesian);
        // rotated box contributes its bounding box only
        assert!(has(&e[0], -(0.5f64).sqrt()) && has(&e[1], 2f64.sqrt()));
        assert!(MeshBuilder::new(MeshOptions::new(0.0)).build(&s).is_err());
    }

    #[test]
    fn cylindrical_point_defined_shapes() {
        // Cartesian mesh, primitives defined in cylindrical coordinates.
        let cyl = CoordSystem::Cylindrical;
        let mut s = ContinuousStructure::new();
        let m = s.add_property(Property::metal("m"));
        s.add_primitive(
            m,
            Primitive::new(Sphere::new([2.0, PI / 2.0, 1.0], 0.5)).with_coord_system(cyl),
        )
        .unwrap();
        s.add_primitive(
            m,
            Primitive::new(Cylinder::new([3.0, PI, 0.0], [3.0, PI, 4.0], 0.25))
                .with_coord_system(cyl),
        )
        .unwrap();
        s.add_primitive(
            m,
            Primitive::new(Wire::new([[1.0, 0.0, 0.0], [1.0, PI / 2.0, 0.0]], 0.0))
                .with_coord_system(cyl),
        )
        .unwrap();
        let l = detect_edges(&s, TypeMask::ANY, CoordSystem::Cartesian);
        let has = |d: usize, v: f64| l[d].iter().any(|x| (x - v).abs() < 1e-9);
        // sphere centre (0, 2, 1)
        assert!(has(0, 0.0) && has(1, 2.0) && has(2, 1.0));
        // cylinder axis at (-3, 0)
        assert!(has(0, -3.0) && has(1, 0.0) && has(2, 4.0));
        // wire vertices (1, 0, 0) and (0, 1, 0)
        assert!(has(0, 1.0) && has(1, 1.0));
        // no raw (rho, alpha) values leaked in
        assert!(!has(1, PI / 2.0) && !has(1, PI));
    }
}
