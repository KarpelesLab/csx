//! Surface tessellation of primitives.
//!
//! [`tessellate`] turns a [`Primitive`] into a triangle mesh in world
//! Cartesian coordinates (transformations applied). Solids produce closed,
//! outward oriented meshes; sheets (polygons, flat boxes) produce open
//! surfaces; curves, points and expression-defined primitives produce empty
//! meshes (see [`polyline`] for curves).

use std::collections::HashMap;

use core::f64::consts::{PI, TAU};

use crate::coords::{CoordSystem, cylindrical_to_cartesian};
use crate::error::Result;
use crate::math::{BoundingBox, Mat4, Vec3};
use crate::primitives::{Polyhedron, Primitive, RotPoly, Shape, plane_axes};
use crate::properties::TypeMask;
use crate::structure::ContinuousStructure;

/// Tessellation resolution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TessellationOptions {
    /// Number of segments used for a full circle (partial arcs use a
    /// proportional number, at least one).
    pub segments: usize,
    /// Number of latitude bands of a sphere (pole to pole).
    pub sphere_rings: usize,
}

impl Default for TessellationOptions {
    fn default() -> Self {
        TessellationOptions {
            segments: 32,
            sphere_rings: 16,
        }
    }
}

impl TessellationOptions {
    fn arc_segments(&self, span: f64) -> usize {
        let full = self.segments.max(3) as f64;
        // Arcs beyond a full turn overlap themselves; cap them (also keeps
        // huge or NaN spans from producing absurd segment counts).
        let span = if span.is_nan() {
            0.0
        } else {
            span.abs().min(TAU)
        };
        ((full * span / TAU).ceil() as usize).max(1)
    }
}

/// An indexed triangle mesh.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TriMesh {
    /// Vertex positions.
    pub vertices: Vec<Vec3>,
    /// Triangles as vertex indices (counter-clockwise seen from outside).
    pub triangles: Vec<[u32; 3]>,
}

impl TriMesh {
    /// Empty mesh.
    pub fn new() -> Self {
        Self::default()
    }

    /// Builds a mesh from a triangle soup, merging bitwise identical
    /// vertices and dropping degenerate triangles (repeated vertices).
    pub fn from_soup(soup: &[[Vec3; 3]]) -> Self {
        let mut index: HashMap<[u64; 3], u32> = HashMap::new();
        let mut m = TriMesh::new();
        for t in soup {
            let ids = t.map(|v| {
                let key = v.0.map(|c| (c + 0.0).to_bits());
                *index.entry(key).or_insert_with(|| {
                    m.vertices.push(Vec3(v.0.map(|c| c + 0.0)));
                    (m.vertices.len() - 1) as u32
                })
            });
            if ids[0] != ids[1] && ids[1] != ids[2] && ids[0] != ids[2] {
                m.triangles.push(ids);
            }
        }
        m
    }

    /// True if there are no triangles.
    pub fn is_empty(&self) -> bool {
        self.triangles.is_empty()
    }

    /// Appends another mesh (no vertex merging).
    pub fn append(&mut self, o: &TriMesh) {
        let off = self.vertices.len() as u32;
        self.vertices.extend_from_slice(&o.vertices);
        self.triangles
            .extend(o.triangles.iter().map(|t| t.map(|i| i + off)));
    }

    /// Reverses the orientation of every triangle.
    pub fn flip(&mut self) {
        for t in &mut self.triangles {
            t.swap(1, 2);
        }
    }

    /// Applies an affine transformation; orientation is preserved for
    /// mirroring matrices.
    pub fn transform(&mut self, m: &Mat4) {
        for v in &mut self.vertices {
            *v = m.transform_point(*v);
        }
        let a = &m.0;
        let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
        if det < 0.0 {
            self.flip();
        }
    }

    /// The triangles as vertex position triples.
    pub fn triangle_soup(&self) -> Vec<[Vec3; 3]> {
        self.triangles
            .iter()
            .map(|t| t.map(|i| self.vertices[i as usize]))
            .collect()
    }

    /// Signed enclosed volume (positive for outward oriented closed meshes).
    pub fn signed_volume(&self) -> f64 {
        self.triangle_soup()
            .iter()
            .map(|[a, b, c]| a.dot(b.cross(*c)))
            .sum::<f64>()
            / 6.0
    }

    /// Total surface area.
    pub fn area(&self) -> f64 {
        self.triangle_soup()
            .iter()
            .map(|[a, b, c]| (*b - *a).cross(*c - *a).norm() / 2.0)
            .sum()
    }

    /// Bounding box of the vertices.
    pub fn bounding_box(&self) -> BoundingBox {
        BoundingBox::from_points(self.vertices.iter().copied())
    }

    /// Converts into a [`Polyhedron`] (for inside tests).
    pub fn to_polyhedron(&self) -> Result<Polyhedron> {
        Polyhedron::new(
            self.vertices.clone(),
            self.triangles.iter().map(|t| t.to_vec()).collect(),
        )
    }
}

/// Collects triangles before welding.
#[derive(Default)]
struct Soup(Vec<[Vec3; 3]>);

impl Soup {
    fn quad(&mut self, a: Vec3, b: Vec3, c: Vec3, d: Vec3) {
        self.0.push([a, b, c]);
        self.0.push([a, c, d]);
    }

    /// Grid surface `f(i, j)` for `i ∈ 0..=nu`, `j ∈ 0..=nv`; quads are
    /// emitted so that the normal is `∂u × ∂v` (reversed if `flip`).
    fn grid(&mut self, nu: usize, nv: usize, flip: bool, f: impl Fn(usize, usize) -> Vec3) {
        for i in 0..nu {
            for j in 0..nv {
                let (a, b, c, d) = (f(i, j), f(i + 1, j), f(i + 1, j + 1), f(i, j + 1));
                if flip {
                    self.quad(a, d, c, b)
                } else {
                    self.quad(a, b, c, d)
                }
            }
        }
    }

    /// Surface of the parameter box `[lo, hi]` mapped through `map`, with
    /// `divs` subdivisions per axis. A periodic axis wraps (its end faces are
    /// omitted). Faces of zero-thickness axes are emitted once.
    fn param_box(
        &mut self,
        lo: Vec3,
        hi: Vec3,
        divs: [usize; 3],
        periodic: [bool; 3],
        map: impl Fn(Vec3) -> Vec3,
    ) {
        let coord = |axis: usize, i: usize| {
            let n = divs[axis];
            if periodic[axis] && i == n {
                lo[axis]
            } else if i == n {
                hi[axis]
            } else {
                lo[axis] + (hi[axis] - lo[axis]) * i as f64 / n as f64
            }
        };
        for k in 0..3 {
            if periodic[k] {
                continue;
            }
            let (u, v) = ((k + 1) % 3, (k + 2) % 3);
            for (side, val) in [(false, lo[k]), (true, hi[k])] {
                if side && hi[k] == lo[k] {
                    continue;
                }
                self.grid(divs[u], divs[v], !side, |i, j| {
                    let mut p = Vec3::ZERO;
                    p[k] = val;
                    p[u] = coord(u, i);
                    p[v] = coord(v, j);
                    map(p)
                });
            }
        }
    }

    fn finish(self) -> TriMesh {
        TriMesh::from_soup(&self.0)
    }
}

/// Orthonormal frame `(e1, e2, d)` (right handed) around direction `d`.
fn frame(d: Vec3) -> (Vec3, Vec3) {
    let a = if d[0].abs() < 0.9 {
        Vec3::unit(0)
    } else {
        Vec3::unit(1)
    };
    let e1 = d.cross(a).normalized();
    (e1, d.cross(e1))
}

fn sphere(s: &mut Soup, c: Vec3, r: f64, opts: &TessellationOptions, flip: bool) {
    if r <= 0.0 {
        return;
    }
    let nt = opts.sphere_rings.max(2);
    let np = opts.segments.max(3);
    s.grid(nt, np, flip, |i, j| {
        let (st, ct) = match i {
            0 => (0.0, 1.0),
            _ if i == nt => (0.0, -1.0),
            _ => (PI * i as f64 / nt as f64).sin_cos(),
        };
        let (sp, cp) = (TAU * (j % np) as f64 / np as f64).sin_cos();
        c + Vec3::new(st * cp, st * sp, ct) * r
    });
}

/// Cylinder (shell) around segment `a`–`b` with radii `[r0, r1]`.
fn cylinder(s: &mut Soup, a: Vec3, b: Vec3, r0: f64, r1: f64, opts: &TessellationOptions) {
    let axis = b - a;
    let len = axis.norm();
    if len == 0.0 || r1 <= 0.0 {
        return;
    }
    let d = axis / len;
    let (e1, e2) = frame(d);
    let np = opts.segments.max(3);
    s.param_box(
        Vec3::new(r0.max(0.0), 0.0, 0.0),
        Vec3::new(r1, TAU, len),
        [1, np, 1],
        [false, true, false],
        |p| {
            let (sn, cs) = p[1].sin_cos();
            a + d * p[2] + (e1 * cs + e2 * sn) * p[0]
        },
    );
}

/// Removes consecutive duplicates (and a closing duplicate).
fn clean_polygon(v: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut out: Vec<(f64, f64)> = Vec::with_capacity(v.len());
    for &p in v {
        if out.last() != Some(&p) {
            out.push(p);
        }
    }
    while out.len() > 1 && out.first() == out.last() {
        out.pop();
    }
    out
}

fn signed_area(v: &[(f64, f64)]) -> f64 {
    let n = v.len();
    (0..n)
        .map(|i| {
            let (a, b) = (v[i], v[(i + 1) % n]);
            a.0 * b.1 - b.0 * a.1
        })
        .sum::<f64>()
        / 2.0
}

/// Ear-clipping triangulation of a simple polygon given counter-clockwise.
/// Returns triangles (CCW) as indices into `v`.
pub fn triangulate_polygon(v: &[(f64, f64)]) -> Vec<[usize; 3]> {
    let n = v.len();
    let mut out = Vec::with_capacity(n.saturating_sub(2));
    if n < 3 {
        return out;
    }
    let scale = v
        .iter()
        .fold(0.0f64, |m, p| m.max(p.0.abs()).max(p.1.abs()))
        .max(1e-300);
    let eps = 1e-12 * scale * scale;
    let cross = |a: (f64, f64), b: (f64, f64), c: (f64, f64)| {
        (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
    };
    let mut idx: Vec<usize> = (0..n).collect();
    let mut guard = 0;
    while idx.len() > 3 && guard < n * n {
        guard += 1;
        let m = idx.len();
        let mut clipped = false;
        for k in 0..m {
            let (ip, ic, inx) = (idx[(k + m - 1) % m], idx[k], idx[(k + 1) % m]);
            let (a, b, c) = (v[ip], v[ic], v[inx]);
            if cross(a, b, c) <= eps {
                continue;
            }
            let blocked = idx.iter().any(|&o| {
                if o == ip || o == ic || o == inx {
                    return false;
                }
                let p = v[o];
                if p == a || p == b || p == c {
                    return false;
                }
                cross(a, b, p) >= -eps && cross(b, c, p) >= -eps && cross(c, a, p) >= -eps
            });
            if !blocked {
                out.push([ip, ic, inx]);
                idx.remove(k);
                clipped = true;
                break;
            }
        }
        if !clipped {
            // Degenerate or self-intersecting input: drop a flat vertex if
            // any, otherwise fan the remainder.
            if let Some(k) = (0..m).find(|&k| {
                cross(v[idx[(k + m - 1) % m]], v[idx[k]], v[idx[(k + 1) % m]]).abs() <= eps
            }) {
                idx.remove(k);
            } else {
                for k in 1..m - 1 {
                    out.push([idx[0], idx[k], idx[k + 1]]);
                }
                return out;
            }
        }
    }
    if idx.len() == 3 && cross(v[idx[0]], v[idx[1]], v[idx[2]]) > eps {
        out.push([idx[0], idx[1], idx[2]]);
    }
    out
}

/// Cleans `v` and returns it counter-clockwise with its triangulation.
fn ccw_polygon(v: &[(f64, f64)]) -> (Vec<(f64, f64)>, Vec<[usize; 3]>) {
    let mut p = clean_polygon(v);
    if signed_area(&p) < 0.0 {
        p.reverse();
    }
    let t = triangulate_polygon(&p);
    (p, t)
}

fn plane_point(n: usize, (a, b): (f64, f64), h: f64) -> Vec3 {
    let (u, v) = plane_axes(n);
    let mut p = Vec3::ZERO;
    p[u] = a;
    p[v] = b;
    p[n] = h;
    p
}

fn extruded_polygon(s: &mut Soup, n: usize, verts: &[(f64, f64)], lo: f64, hi: f64) {
    let (p, tris) = ccw_polygon(verts);
    if p.len() < 3 {
        return;
    }
    let at = |i: usize, h: f64| plane_point(n, p[i], h);
    for t in &tris {
        s.0.push([at(t[0], hi), at(t[1], hi), at(t[2], hi)]);
        if hi != lo {
            s.0.push([at(t[0], lo), at(t[2], lo), at(t[1], lo)]);
        }
    }
    if hi != lo {
        for i in 0..p.len() {
            let j = (i + 1) % p.len();
            s.quad(at(i, lo), at(j, lo), at(j, hi), at(i, hi));
        }
    }
}

fn rotpoly(s: &mut Soup, r: &RotPoly, opts: &TessellationOptions) {
    let a = r.rot_axis_dir;
    let q = 3 - a - r.norm_dir;
    let w = Vec3::unit(a).cross(Vec3::unit(q));
    let (u, _) = plane_axes(r.norm_dir);
    // Polygon in (radial, axial) coordinates.
    let rh: Vec<(f64, f64)> = r
        .vertices
        .iter()
        .map(|&(c0, c1)| if u == q { (c0, c1) } else { (c1, c0) })
        .collect();
    let (p, tris) = ccw_polygon(&rh);
    if p.len() < 3 {
        return;
    }
    let (start, stop) = (r.start_angle.value(), r.stop_angle.value());
    let span = stop - start;
    let full = span.abs() >= TAU - 1e-12;
    let n = opts.arc_segments(if full { TAU } else { span });
    let angle = |k: usize| {
        if full {
            start + TAU * (k % n) as f64 / n as f64
        } else {
            start + span * k as f64 / n as f64
        }
    };
    let at = |vi: usize, k: usize| {
        let (rad, h) = p[vi];
        let (sn, cs) = angle(k).sin_cos();
        Vec3::unit(a) * h + (Vec3::unit(q) * cs + w * sn) * rad
    };
    for i in 0..p.len() {
        let j = (i + 1) % p.len();
        for k in 0..n {
            s.quad(at(i, k), at(i, k + 1), at(j, k + 1), at(j, k));
        }
    }
    if !full {
        for t in &tris {
            s.0.push([at(t[0], 0), at(t[1], 0), at(t[2], 0)]);
            s.0.push([at(t[0], n), at(t[2], n), at(t[1], n)]);
        }
    }
}

fn box_shape(s: &mut Soup, a: Vec3, b: Vec3, cs: CoordSystem, opts: &TessellationOptions) {
    let (lo, hi) = (a.min(b), a.max(b));
    match cs {
        CoordSystem::Cartesian => s.param_box(lo, hi, [1, 1, 1], [false; 3], |p| p),
        CoordSystem::Cylindrical => {
            let span = hi[1] - lo[1];
            let full = span >= TAU - 1e-12;
            let (lo, hi) = if full {
                (lo, Vec3::new(hi[0], lo[1] + TAU, hi[2]))
            } else {
                (lo, hi)
            };
            let n = opts.arc_segments(hi[1] - lo[1]);
            s.param_box(
                lo,
                hi,
                [1, n, 1],
                [false, full, false],
                cylindrical_to_cartesian,
            );
        }
    }
}

/// Tessellates a primitive in world Cartesian coordinates. `default_cs` is
/// the structure's coordinate system (used by boxes without their own).
pub fn tessellate(
    prim: &Primitive,
    default_cs: CoordSystem,
    opts: &TessellationOptions,
) -> TriMesh {
    let mut s = Soup::default();
    let cs = prim.effective_cs(default_cs);
    // Point-defined shapes are converted to Cartesian; boxes keep `cs`.
    let shape = prim.shape.to_cartesian(cs);
    let mut mesh = match &*shape {
        Shape::Box(b) => {
            box_shape(&mut s, b.start.value(), b.stop.value(), cs, opts);
            s.finish()
        }
        Shape::MultiBox(m) => {
            let mut out = TriMesh::new();
            for (a, b) in &m.boxes {
                let mut t = Soup::default();
                // Multi-boxes use the structure's coordinate system.
                box_shape(
                    &mut t,
                    a.value(),
                    b.value(),
                    prim.shape_cs(default_cs),
                    opts,
                );
                out.append(&t.finish());
            }
            out
        }
        Shape::Sphere(sp) => {
            sphere(&mut s, sp.center.value(), sp.radius.value(), opts, false);
            s.finish()
        }
        Shape::SphericalShell(sh) => {
            let (r, w) = (sh.radius.value(), sh.shell_width.value() / 2.0);
            let c = sh.center.value();
            sphere(&mut s, c, r + w, opts, false);
            if w > 0.0 {
                sphere(&mut s, c, r - w, opts, true);
            }
            s.finish()
        }
        Shape::Cylinder(c) => {
            cylinder(
                &mut s,
                c.start.value(),
                c.stop.value(),
                0.0,
                c.radius.value(),
                opts,
            );
            s.finish()
        }
        Shape::CylindricalShell(c) => {
            let (r, w) = (c.radius.value(), c.shell_width.value() / 2.0);
            cylinder(&mut s, c.start.value(), c.stop.value(), r - w, r + w, opts);
            s.finish()
        }
        Shape::Polygon(p) => {
            let e = p.elevation.value();
            extruded_polygon(&mut s, p.norm_dir.min(2), &p.vertices, e, e);
            s.finish()
        }
        Shape::LinPoly(p) => {
            let (e, l) = (p.elevation.value(), p.length.value());
            extruded_polygon(
                &mut s,
                p.norm_dir.min(2),
                &p.vertices,
                e.min(e + l),
                e.max(e + l),
            );
            s.finish()
        }
        Shape::RotPoly(r) => {
            if r.norm_dir < 3 && r.rot_axis_dir < 3 && r.norm_dir != r.rot_axis_dir {
                rotpoly(&mut s, r, opts);
            }
            let mut m = s.finish();
            if m.signed_volume() < 0.0 {
                m.flip();
            }
            m
        }
        Shape::Wire(w) => {
            let r = w.radius.value();
            let pts: Vec<Vec3> = w.points.iter().map(|p| p.value()).collect();
            let mut out = TriMesh::new();
            if r > 0.0 {
                for seg in pts.windows(2) {
                    let mut t = Soup::default();
                    cylinder(&mut t, seg[0], seg[1], 0.0, r, opts);
                    out.append(&t.finish());
                }
                for p in pts.iter().skip(1).take(pts.len().saturating_sub(2)) {
                    let mut t = Soup::default();
                    sphere(&mut t, *p, r, opts, false);
                    out.append(&t.finish());
                }
            }
            out
        }
        Shape::Polyhedron(p) => TriMesh {
            vertices: p.vertices().to_vec(),
            triangles: p.triangles().to_vec(),
        },
        Shape::PolyhedronReader(r) => r
            .mesh()
            .map(|p| TriMesh {
                vertices: p.vertices().to_vec(),
                triangles: p.triangles().to_vec(),
            })
            .unwrap_or_default(),
        Shape::Curve(_) | Shape::Point(_) | Shape::UserDefined(_) => TriMesh::new(),
    };
    if !prim.transform.is_empty() {
        mesh.transform(prim.transform.matrix());
    }
    mesh
}

/// World-space polyline of a curve or wire primitive (empty for others).
/// `default_cs` is the structure's coordinate system; vertices of
/// primitives in a cylindrical system are converted from `(ρ, α, z)`.
pub fn polyline(prim: &Primitive, default_cs: CoordSystem) -> Vec<Vec3> {
    let shape = prim.cartesian_shape(default_cs);
    let pts: Vec<Vec3> = match &*shape {
        Shape::Curve(c) => c.points.iter().map(|p| p.value()).collect(),
        Shape::Wire(w) => w.points.iter().map(|p| p.value()).collect(),
        _ => return Vec::new(),
    };
    if prim.transform.is_empty() {
        pts
    } else {
        pts.into_iter().map(|p| prim.transform.apply(p)).collect()
    }
}

/// One tessellated primitive of a structure.
#[derive(Debug, Clone, PartialEq)]
pub struct TessellatedPrimitive {
    /// Index of the owning property in [`ContinuousStructure::properties`].
    pub property_index: usize,
    /// Index of the primitive within the property.
    pub primitive_index: usize,
    /// Priority of the primitive.
    pub priority: i32,
    /// Surface mesh.
    pub mesh: TriMesh,
    /// Polyline (curves and wires).
    pub polyline: Vec<Vec3>,
}

/// Tessellates every primitive of properties matching `mask`.
pub fn tessellate_structure(
    s: &ContinuousStructure,
    mask: impl Into<TypeMask>,
    opts: &TessellationOptions,
) -> Vec<TessellatedPrimitive> {
    let mask = mask.into();
    let mut out = Vec::new();
    for (pi, prop) in s.properties().iter().enumerate() {
        if !mask.contains(prop.property_type()) {
            continue;
        }
        for (k, prim) in prop.primitives.iter().enumerate() {
            out.push(TessellatedPrimitive {
                property_index: pi,
                primitive_index: k,
                priority: prim.priority,
                mesh: tessellate(prim, s.coord_system, opts),
                polyline: polyline(prim, s.coord_system),
            });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{
        BoxShape, Cylinder, CylindricalShell, LinPoly, MultiBox, Polygon, Sphere, SphericalShell,
        Wire,
    };
    use crate::transform::Transform;
    use core::f64::consts::FRAC_PI_2;

    /// Deterministic LCG in [0, 1).
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    fn opts() -> TessellationOptions {
        TessellationOptions {
            segments: 48,
            sphere_rings: 24,
        }
    }

    /// Every directed edge appears once and its reverse once.
    fn closed_and_consistent(m: &TriMesh) -> bool {
        let mut edges: HashMap<(u32, u32), i32> = HashMap::new();
        for t in &m.triangles {
            for k in 0..3 {
                *edges.entry((t[k], t[(k + 1) % 3])).or_default() += 1;
            }
        }
        !m.is_empty()
            && edges
                .iter()
                .all(|(&(a, b), &n)| n == 1 && edges.get(&(b, a)) == Some(&1))
    }

    /// Compares inside tests of the primitive and its tessellation at random
    /// points, ignoring points within `margin` of the surface (judged by the
    /// primitive itself on a slightly shrunk/expanded probe).
    fn agree(prim: &Primitive, margin: f64) {
        let m = tessellate(prim, CoordSystem::Cartesian, &opts());
        assert!(
            closed_and_consistent(&m),
            "{} not closed",
            prim.shape.type_name()
        );
        assert!(m.signed_volume() > 0.0);
        let poly = m.to_polyhedron().unwrap();
        assert!(poly.is_closed());
        let bb = prim
            .bounding_box(CoordSystem::Cartesian)
            .expanded(margin * 2.0);
        let mut rng = Lcg(7);
        let mut checked = 0;
        for _ in 0..3000 {
            let p = Vec3::new(
                bb.min[0] + rng.next() * (bb.max[0] - bb.min[0]),
                bb.min[1] + rng.next() * (bb.max[1] - bb.min[1]),
                bb.min[2] + rng.next() * (bb.max[2] - bb.min[2]),
            );
            let inside = prim.contains(p, CoordSystem::Cartesian);
            // Skip points near the boundary.
            let near = [Vec3::unit(0), Vec3::unit(1), Vec3::unit(2)]
                .iter()
                .any(|d| {
                    prim.contains(p + *d * margin, CoordSystem::Cartesian) != inside
                        || prim.contains(p - *d * margin, CoordSystem::Cartesian) != inside
                });
            if near {
                continue;
            }
            checked += 1;
            assert_eq!(
                poly.contains(p),
                inside,
                "{} at {p:?}",
                prim.shape.type_name()
            );
        }
        assert!(checked > 1000);
    }

    #[test]
    fn boxes() {
        let b = Primitive::new(BoxShape::new([1.0, 2.0, 3.0], [-1.0, 0.0, 0.0]));
        let m = tessellate(&b, CoordSystem::Cartesian, &opts());
        assert_eq!(m.vertices.len(), 8);
        assert_eq!(m.triangles.len(), 12);
        assert!((m.signed_volume() - 12.0).abs() < 1e-12);
        agree(&b, 1e-6);
        // flat box → single rectangle
        let f = Primitive::new(BoxShape::new([0.0, 0.0, 1.0], [2.0, 3.0, 1.0]));
        let m = tessellate(&f, CoordSystem::Cartesian, &opts());
        assert_eq!(m.triangles.len(), 2);
        assert!((m.area() - 6.0).abs() < 1e-12);
        // line → empty
        let l = Primitive::new(BoxShape::new([0.0, 0.0, 1.0], [2.0, 0.0, 1.0]));
        assert!(tessellate(&l, CoordSystem::Cartesian, &opts()).is_empty());
        let mb = Primitive::new(MultiBox {
            boxes: vec![
                ([0.0, 0.0, 0.0].into(), [1.0, 1.0, 1.0].into()),
                ([2.0, 0.0, 0.0].into(), [3.0, 1.0, 1.0].into()),
            ],
        });
        let m = tessellate(&mb, CoordSystem::Cartesian, &opts());
        assert!(closed_and_consistent(&m));
        assert!((m.signed_volume() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn cylindrical_boxes() {
        let sector = Primitive::new(BoxShape::new([1.0, -0.5, 0.0], [2.0, 1.0, 1.0]))
            .with_coord_system(CoordSystem::Cylindrical);
        agree(&sector, 0.01);
        let exact = 0.5 * 1.5 * (4.0 - 1.0);
        let v = tessellate(&sector, CoordSystem::Cartesian, &opts()).signed_volume();
        assert!((v - exact).abs() / exact < 0.01, "{v} vs {exact}");
        // full annulus and full disc (inner radius 0)
        for r0 in [0.5, 0.0] {
            let ring = Primitive::new(BoxShape::new([r0, 0.0, 0.0], [1.0, TAU, 1.0]))
                .with_coord_system(CoordSystem::Cylindrical);
            let m = tessellate(&ring, CoordSystem::Cartesian, &opts());
            assert!(closed_and_consistent(&m), "r0 = {r0}");
            let exact = PI * (1.0 - r0 * r0);
            assert!((m.signed_volume() - exact).abs() / exact < 0.01);
        }
        // the structure default coordinate system is honoured
        let b = Primitive::new(BoxShape::new([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]));
        let m = tessellate(&b, CoordSystem::Cylindrical, &opts());
        assert!(m.vertices.len() > 8);
    }

    #[test]
    fn round_solids() {
        let s = Primitive::new(Sphere::new([1.0, -2.0, 0.5], 2.0));
        agree(&s, 0.05);
        let v = tessellate(&s, CoordSystem::Cartesian, &opts()).signed_volume();
        let exact = 4.0 / 3.0 * PI * 8.0;
        assert!((v - exact).abs() / exact < 0.02);
        agree(
            &Primitive::new(SphericalShell {
                center: [0.0, 0.0, 0.0].into(),
                radius: 2.0.into(),
                shell_width: 1.0.into(),
            }),
            0.05,
        );
        agree(
            &Primitive::new(Cylinder::new([0.0, 0.0, 0.0], [1.0, 2.0, 3.0], 0.7)),
            0.03,
        );
        agree(
            &Primitive::new(CylindricalShell {
                start: [0.0, 0.0, -1.0].into(),
                stop: [0.0, 0.0, 1.0].into(),
                radius: 1.0.into(),
                shell_width: 0.4.into(),
            }),
            0.03,
        );
        // zero-width shell → open mantle surface
        let m = tessellate(
            &Primitive::new(CylindricalShell {
                start: [0.0, 0.0, 0.0].into(),
                stop: [0.0, 0.0, 1.0].into(),
                radius: 1.0.into(),
                shell_width: 0.0.into(),
            }),
            CoordSystem::Cartesian,
            &opts(),
        );
        assert!(!m.is_empty());
        assert!((m.area() - TAU).abs() / TAU < 0.01);
    }

    #[test]
    fn polygons() {
        // concave L shape given clockwise
        let l = vec![
            (0.0, 0.0),
            (0.0, 2.0),
            (1.0, 2.0),
            (1.0, 1.0),
            (2.0, 1.0),
            (2.0, 0.0),
            (0.0, 0.0),
        ];
        let p = Primitive::new(Polygon::new(2, 0.5, l.clone()));
        let m = tessellate(&p, CoordSystem::Cartesian, &opts());
        assert_eq!(m.triangles.len(), 4);
        assert!((m.area() - 3.0).abs() < 1e-12);
        // normal is +z
        for [a, b, c] in m.triangle_soup() {
            assert!((b - a).cross(c - a)[2] > 0.0);
        }
        for n in 0..3 {
            agree(&Primitive::new(LinPoly::new(n, -0.5, 1.5, l.clone())), 1e-6);
            agree(&Primitive::new(LinPoly::new(n, 0.5, -1.5, l.clone())), 1e-6);
        }
        // collinear vertices keep the extrusion closed
        let sq = vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)];
        let m = tessellate(
            &Primitive::new(LinPoly::new(2, 0.0, 1.0, sq)),
            CoordSystem::Cartesian,
            &opts(),
        );
        assert!(closed_and_consistent(&m));
        assert!((m.signed_volume() - 4.0).abs() < 1e-12);
    }

    #[test]
    fn triangulation_random_star() {
        let mut rng = Lcg(3);
        for _ in 0..50 {
            let n = 5 + (rng.next() * 20.0) as usize;
            let v: Vec<(f64, f64)> = (0..n)
                .map(|i| {
                    let a = TAU * i as f64 / n as f64;
                    let r = 0.3 + rng.next();
                    (r * a.cos(), r * a.sin())
                })
                .collect();
            let t = triangulate_polygon(&v);
            assert_eq!(t.len(), n - 2);
            let area: f64 = t
                .iter()
                .map(|t| signed_area(&[v[t[0]], v[t[1]], v[t[2]]]))
                .sum();
            assert!((area - signed_area(&v)).abs() < 1e-9);
        }
    }

    #[test]
    fn revolved() {
        // washer: square [1,2]x[0,1] in x/z plane revolved around z
        let full = RotPoly::new(1, 2, vec![(0.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0)]);
        agree(&Primitive::new(full.clone()), 0.03);
        let v = tessellate(
            &Primitive::new(full.clone()),
            CoordSystem::Cartesian,
            &opts(),
        )
        .signed_volume();
        let exact = PI * 3.0;
        assert!((v - exact).abs() / exact < 0.01);
        let mut part = full.clone();
        part.start_angle = 0.3.into();
        part.stop_angle = 2.0.into();
        agree(&Primitive::new(part), 0.03);
        // polygon touching the axis (a cone-like solid)
        let cone = RotPoly::new(0, 2, vec![(0.0, 0.0), (1.0, 0.0), (0.0, 2.0)]);
        // norm x: u = y (radial), v = z (axial)
        let m = tessellate(&Primitive::new(cone), CoordSystem::Cartesian, &opts());
        assert!(closed_and_consistent(&m));
        let exact = PI * 1.0 * 2.0 / 3.0;
        assert!((m.signed_volume() - exact).abs() / exact < 0.01);
    }

    #[test]
    fn wires_curves_transforms() {
        let w = Primitive::new(Wire::new(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
            0.1,
        ));
        let m = tessellate(&w, CoordSystem::Cartesian, &opts());
        assert!(!m.is_empty());
        assert_eq!(polyline(&w, CoordSystem::Cartesian).len(), 3);
        let mut t = Transform::new();
        t.scale(2.0).translate(Vec3::new(5.0, 0.0, 0.0));
        let b = Primitive::new(BoxShape::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])).with_transform(t);
        let m = tessellate(&b, CoordSystem::Cartesian, &opts());
        assert!((m.signed_volume() - 8.0).abs() < 1e-12);
        assert_eq!(m.bounding_box().min, Vec3::new(5.0, 0.0, 0.0));
        // mirroring keeps outward orientation
        let mut t = Transform::new();
        t.matrix_op(&Mat4::scaling(Vec3::new(-1.0, 1.0, 1.0)));
        let b = Primitive::new(BoxShape::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])).with_transform(t);
        assert!(tessellate(&b, CoordSystem::Cartesian, &opts()).signed_volume() > 0.0);
        let mut rot = Transform::new();
        rot.rotate(Vec3::new(1.0, 1.0, 0.0), 0.7)
            .translate(Vec3::new(0.0, 0.0, 3.0));
        agree(
            &Primitive::new(Cylinder::new([0.0, 0.0, 0.0], [0.0, 0.0, 2.0], 0.5))
                .with_transform(rot),
            0.03,
        );
    }

    #[test]
    fn structure() {
        let mut s = ContinuousStructure::new();
        let id = s.add_property(crate::Property::metal("m"));
        s.add_primitive(id, BoxShape::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]))
            .unwrap();
        s.add_primitive(
            id,
            crate::primitives::Curve::new([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]),
        )
        .unwrap();
        let t = tessellate_structure(&s, TypeMask::ANY, &opts());
        assert_eq!(t.len(), 2);
        assert_eq!(t[0].mesh.triangles.len(), 12);
        assert_eq!(t[1].polyline.len(), 2);
    }

    #[test]
    fn huge_angles_are_capped() {
        let o = TessellationOptions::default();
        assert_eq!(o.arc_segments(1e300), o.arc_segments(TAU));
        assert_eq!(o.arc_segments(f64::NAN), 1);
    }

    #[test]
    fn cylindrical_point_defined_shapes() {
        let cyl = CoordSystem::Cylindrical;
        let o = TessellationOptions::default();
        // Cylinder along z at (rho=1, alpha=pi/2) → axis at (0, 1).
        let c = Primitive::new(Cylinder::new(
            [1.0, FRAC_PI_2, 0.0],
            [1.0, FRAC_PI_2, 2.0],
            0.5,
        ))
        .with_coord_system(cyl);
        let b = tessellate(&c, CoordSystem::Cartesian, &o).bounding_box();
        assert!((b.min - Vec3::new(-0.5, 0.5, 0.0)).norm() < 1e-9, "{b:?}");
        assert!((b.max - Vec3::new(0.5, 1.5, 2.0)).norm() < 1e-9, "{b:?}");
        // Sphere with the inherited (cylindrical) system.
        let s = Primitive::new(Sphere::new([2.0, core::f64::consts::PI, 0.0], 1.0));
        let b = tessellate(&s, cyl, &o).bounding_box();
        assert!(
            (b.center() - Vec3::new(-2.0, 0.0, 0.0)).norm() < 1e-9,
            "{b:?}"
        );
        assert!(tessellate(&s, cyl, &o).to_polyhedron().unwrap().is_closed());
        // Spherical shell and cylindrical shell move the same way.
        let sh = Primitive::new(SphericalShell {
            center: [2.0, FRAC_PI_2, 0.0].into(),
            radius: 1.0.into(),
            shell_width: 0.2.into(),
        })
        .with_coord_system(cyl);
        let b = tessellate(&sh, CoordSystem::Cartesian, &o).bounding_box();
        assert!((b.center() - Vec3::new(0.0, 2.0, 0.0)).norm() < 1e-9);
        let csh = Primitive::new(CylindricalShell {
            start: [2.0, 0.0, 0.0].into(),
            stop: [2.0, 0.0, 1.0].into(),
            radius: 0.5.into(),
            shell_width: 0.1.into(),
        })
        .with_coord_system(cyl);
        let b = tessellate(&csh, CoordSystem::Cartesian, &o).bounding_box();
        assert!((b.center() - Vec3::new(2.0, 0.0, 0.5)).norm() < 1e-9);
        // Wire tube and polyline.
        let w = Primitive::new(Wire::new([[1.0, 0.0, 0.0], [1.0, FRAC_PI_2, 0.0]], 0.1))
            .with_coord_system(cyl);
        let pl = polyline(&w, CoordSystem::Cartesian);
        assert!((pl[1] - Vec3::new(0.0, 1.0, 0.0)).norm() < 1e-12);
        let b = tessellate(&w, CoordSystem::Cartesian, &o).bounding_box();
        let e = 1.0 + 0.1 * core::f64::consts::FRAC_1_SQRT_2;
        assert!((b.max - Vec3::new(e, e, 0.1)).norm() < 1e-9, "{b:?}");
    }
}
