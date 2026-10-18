//! Polyhedral surfaces and solids.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::math::{BoundingBox, Vec3};

use super::bvh::Bvh;
use super::geom::REL_TOL;

/// Fixed, "irrational" ray directions used for the parity test. Choosing
/// directions that are not aligned with any axis or diagonal makes grazing
/// hits on typical CAD meshes unlikely.
const RAY_DIRS: [[f64; 3]; 3] = [
    [
        0.577_215_664_901_532_9,
        0.618_033_988_749_894_8,
        0.533_402_669_220_018_2,
    ],
    [
        -0.734_513_220_164_721_4,
        0.271_828_182_845_904_5,
        0.621_498_663_213_558_7,
    ],
    [
        0.198_062_264_195_162_1,
        -0.881_373_587_019_543,
        0.429_449_408_670_437_8,
    ],
];

/// A polyhedron defined by vertices and polygonal faces (vertex index
/// lists). Faces with more than three vertices are fan-triangulated.
///
/// A *closed* polyhedron (every edge shared by exactly two faces) is a solid;
/// a point is inside when a ray from it crosses the surface an odd number
/// of times (several rays vote; degenerate cases fall back to the
/// generalised winding number). An open polyhedron is a surface and only
/// points on it are "inside". A bounding volume hierarchy over the
/// triangles keeps queries roughly logarithmic in the face count.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Polyhedron {
    vertices: Vec<Vec3>,
    faces: Vec<Vec<u32>>,
    tris: Vec<[u32; 3]>,
    bbox: BoundingBox,
    closed: bool,
    bvh: Bvh,
}

impl Polyhedron {
    /// Builds a polyhedron, validating face indices.
    pub fn new(vertices: Vec<Vec3>, faces: Vec<Vec<u32>>) -> Result<Self> {
        let nv = vertices.len();
        let mut tris = Vec::new();
        for (i, f) in faces.iter().enumerate() {
            if f.len() < 3 {
                return Err(Error::Invalid(format!(
                    "face {i} has fewer than 3 vertices"
                )));
            }
            if let Some(&bad) = f.iter().find(|&&k| k as usize >= nv) {
                return Err(Error::Invalid(format!(
                    "face {i} references vertex {bad} of {nv}"
                )));
            }
            for k in 1..f.len() - 1 {
                tris.push([f[0], f[k], f[k + 1]]);
            }
        }
        if let Some(v) = vertices.iter().find(|v| !v.is_finite()) {
            return Err(Error::Invalid(format!("non-finite vertex {v:?}")));
        }
        let bbox = BoundingBox::from_points(vertices.iter().copied());
        let closed = is_closed(&faces);
        let boxes: Vec<BoundingBox> = tris
            .iter()
            .map(|t| BoundingBox::from_points(t.map(|i| vertices[i as usize])))
            .collect();
        let bvh = Bvh::build(&boxes);
        Ok(Polyhedron {
            vertices,
            faces,
            tris,
            bbox,
            closed,
            bvh,
        })
    }

    /// Builds a polyhedron from a triangle soup, merging identical vertices.
    pub fn from_triangles(triangles: &[[Vec3; 3]]) -> Result<Self> {
        let mut index: HashMap<[u64; 3], u32> = HashMap::new();
        let mut vertices = Vec::new();
        let mut faces = Vec::with_capacity(triangles.len());
        for t in triangles {
            let mut f = Vec::with_capacity(3);
            for v in t {
                let key = v.0.map(|c| (c + 0.0).to_bits());
                let id = *index.entry(key).or_insert_with(|| {
                    vertices.push(*v);
                    (vertices.len() - 1) as u32
                });
                f.push(id);
            }
            if f[0] != f[1] && f[1] != f[2] && f[0] != f[2] {
                faces.push(f);
            }
        }
        Polyhedron::new(vertices, faces)
    }

    /// Vertex coordinates.
    pub fn vertices(&self) -> &[Vec3] {
        &self.vertices
    }

    /// Faces as vertex index lists.
    pub fn faces(&self) -> &[Vec<u32>] {
        &self.faces
    }

    /// Triangulated faces.
    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.tris
    }

    /// True if the surface is closed (watertight).
    pub fn is_closed(&self) -> bool {
        self.closed
    }

    /// Axis aligned bounding box of the vertices.
    pub fn bounding_box(&self) -> BoundingBox {
        self.bbox
    }

    fn tri(&self, t: &[u32; 3]) -> [Vec3; 3] {
        t.map(|i| self.vertices[i as usize])
    }

    /// Generalised winding number of `p` (≈1 inside a closed, outward
    /// oriented mesh, ≈0 outside; sign flips with orientation).
    pub fn winding_number(&self, p: Vec3) -> f64 {
        let mut total = 0.0;
        for t in &self.tris {
            let [a, b, c] = self.tri(t).map(|v| v - p);
            let (la, lb, lc) = (a.norm(), b.norm(), c.norm());
            let det = a.dot(b.cross(c));
            let den = la * lb * lc + a.dot(b) * lc + a.dot(c) * lb + b.dot(c) * la;
            total += 2.0 * det.atan2(den);
        }
        total / (4.0 * core::f64::consts::PI)
    }

    /// Absolute tolerance used for the on-surface test.
    fn surface_eps(&self) -> f64 {
        let scale = self
            .bbox
            .size()
            .norm()
            .max(self.bbox.max.norm())
            .max(self.bbox.min.norm());
        REL_TOL * scale.max(1e-300)
    }

    /// True if `p` lies on the surface (within tolerance).
    pub fn on_surface(&self, p: Vec3) -> bool {
        let eps = self.surface_eps();
        let mut found = false;
        self.bvh.visit(
            |b| b.expanded(eps).contains(p),
            |i| {
                let [a, b, c] = self.tri(&self.tris[i as usize]);
                found = (closest_on_triangle(p, a, b, c) - p).norm() <= eps;
                !found
            },
        );
        found
    }

    /// Counts the crossings of the ray `p + t·d` (`t > 0`) with the surface.
    /// Returns `None` when a hit is too close to an edge or vertex for the
    /// parity to be trusted.
    fn ray_crossings(&self, p: Vec3, d: Vec3) -> Option<usize> {
        let inv = Vec3::new(1.0 / d[0], 1.0 / d[1], 1.0 / d[2]);
        let mut count = 0usize;
        let mut reliable = true;
        self.bvh.visit(
            |b| ray_hits_box(p, inv, b),
            |i| {
                let [a, b, c] = self.tri(&self.tris[i as usize]);
                match ray_triangle(p, d, a, b, c) {
                    RayHit::Miss => {}
                    RayHit::Hit => count += 1,
                    RayHit::Degenerate => reliable = false,
                }
                reliable
            },
        );
        reliable.then_some(count)
    }

    /// Parity based inside test for closed meshes (majority vote of up to
    /// three rays, winding number fallback).
    fn inside_closed(&self, p: Vec3) -> bool {
        let mut votes = [0usize; 2];
        for d in RAY_DIRS {
            if let Some(n) = self.ray_crossings(p, Vec3(d)) {
                votes[n % 2] += 1;
                if votes[0] == 2 || votes[1] == 2 {
                    break;
                }
            }
        }
        match votes {
            [0, 0] => self.winding_number(p).abs() > 0.5,
            [o, i] => i > o || (i == o && self.winding_number(p).abs() > 0.5),
        }
    }

    /// Inside test: points on the surface count as inside.
    pub fn contains(&self, p: Vec3) -> bool {
        if !self
            .bbox
            .expanded(REL_TOL * self.bbox.size().norm())
            .contains(p)
        {
            return false;
        }
        if self.on_surface(p) {
            return true;
        }
        self.closed && self.inside_closed(p)
    }
}

/// Slab test of a ray (origin `o`, inverse direction `inv`, `t ≥ 0`)
/// against a box.
#[inline]
fn ray_hits_box(o: Vec3, inv: Vec3, b: &BoundingBox) -> bool {
    let mut t0 = 0.0f64;
    let mut t1 = f64::INFINITY;
    for k in 0..3 {
        let a = (b.min[k] - o[k]) * inv[k];
        let c = (b.max[k] - o[k]) * inv[k];
        let (lo, hi) = if a <= c { (a, c) } else { (c, a) };
        t0 = t0.max(lo);
        t1 = t1.min(hi);
    }
    t0 <= t1
}

enum RayHit {
    Miss,
    Hit,
    /// The ray passes (nearly) through an edge or vertex.
    Degenerate,
}

/// Möller–Trumbore ray/triangle intersection for `t > 0`.
fn ray_triangle(o: Vec3, d: Vec3, a: Vec3, b: Vec3, c: Vec3) -> RayHit {
    const EPS: f64 = 1e-9;
    let e1 = b - a;
    let e2 = c - a;
    let h = d.cross(e2);
    let det = e1.dot(h);
    let scale = e1.norm() * e2.norm();
    if scale == 0.0 || det.abs() <= 1e-12 * scale {
        // Degenerate triangle, or ray parallel to its plane: a ray lying in
        // the plane is caught by the edge hits of the neighbouring faces.
        return RayHit::Miss;
    }
    let f = 1.0 / det;
    let s = o - a;
    let u = f * s.dot(h);
    if !(-EPS..=1.0 + EPS).contains(&u) {
        return RayHit::Miss;
    }
    let q = s.cross(e1);
    let v = f * d.dot(q);
    if v < -EPS || u + v > 1.0 + EPS {
        return RayHit::Miss;
    }
    let t = f * e2.dot(q);
    if t <= 0.0 {
        return RayHit::Miss;
    }
    if u <= EPS || v <= EPS || u + v >= 1.0 - EPS {
        return RayHit::Degenerate;
    }
    RayHit::Hit
}

fn is_closed(faces: &[Vec<u32>]) -> bool {
    if faces.is_empty() {
        return false;
    }
    let mut edges: HashMap<(u32, u32), i32> = HashMap::new();
    for f in faces {
        for k in 0..f.len() {
            let (a, b) = (f[k], f[(k + 1) % f.len()]);
            *edges.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    edges.values().all(|&n| n == 2)
}

/// Closest point to `p` on triangle `abc` (Ericson, Real-Time Collision
/// Detection, §5.1.5).
fn closest_on_triangle(p: Vec3, a: Vec3, b: Vec3, c: Vec3) -> Vec3 {
    let ab = b - a;
    let ac = c - a;
    let ap = p - a;
    let d1 = ab.dot(ap);
    let d2 = ac.dot(ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = p - b;
    let d3 = ab.dot(bp);
    let d4 = ac.dot(bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return a + ab * (d1 / (d1 - d3));
    }
    let cp = p - c;
    let d5 = ab.dot(cp);
    let d6 = ac.dot(cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return a + ac * (d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    let denom = va + vb + vc;
    if denom == 0.0 {
        return a;
    }
    a + ab * (vb / denom) + ac * (vc / denom)
}

#[cfg(test)]
pub(crate) fn unit_cube() -> Polyhedron {
    let v = (0..8)
        .map(|i| Vec3::new((i & 1) as f64, ((i >> 1) & 1) as f64, ((i >> 2) & 1) as f64))
        .collect();
    let faces = vec![
        vec![0, 2, 3, 1],
        vec![4, 5, 7, 6],
        vec![0, 1, 5, 4],
        vec![2, 6, 7, 3],
        vec![0, 4, 6, 2],
        vec![1, 3, 7, 5],
    ];
    Polyhedron::new(v, faces).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cube() {
        let c = unit_cube();
        assert!(c.is_closed());
        assert_eq!(c.triangles().len(), 12);
        assert!((c.winding_number(Vec3::new(0.5, 0.5, 0.5)).abs() - 1.0).abs() < 1e-9);
        assert!(c.contains(Vec3::new(0.2, 0.7, 0.9)));
        assert!(c.contains(Vec3::new(1.0, 0.5, 0.5)));
        assert!(c.contains(Vec3::new(1.0, 1.0, 1.0)));
        assert!(!c.contains(Vec3::new(1.2, 0.5, 0.5)));
        assert!(!c.contains(Vec3::new(-0.1, 0.5, 0.5)));
    }

    #[test]
    fn open_surface() {
        let p = Polyhedron::new(
            vec![Vec3::ZERO, Vec3::unit(0), Vec3::unit(1)],
            vec![vec![0, 1, 2]],
        )
        .unwrap();
        assert!(!p.is_closed());
        assert!(p.contains(Vec3::new(0.2, 0.2, 0.0)));
        assert!(!p.contains(Vec3::new(0.2, 0.2, 0.1)));
        assert!(Polyhedron::new(vec![Vec3::ZERO], vec![vec![0, 1, 2]]).is_err());
    }

    #[test]
    fn soup() {
        let c = unit_cube();
        let tris: Vec<[Vec3; 3]> = c.triangles().iter().map(|t| c.tri(t)).collect();
        let s = Polyhedron::from_triangles(&tris).unwrap();
        assert_eq!(s.vertices().len(), 8);
        assert!(s.is_closed());
    }

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

    /// UV sphere (outward oriented).
    fn sphere(n: u32, r: f64) -> Polyhedron {
        let mut v = vec![Vec3::new(0.0, 0.0, r)];
        for i in 1..n {
            let th = core::f64::consts::PI * i as f64 / n as f64;
            for j in 0..n {
                let ph = core::f64::consts::TAU * j as f64 / n as f64;
                v.push(Vec3::new(th.sin() * ph.cos(), th.sin() * ph.sin(), th.cos()) * r);
            }
        }
        v.push(Vec3::new(0.0, 0.0, -r));
        let south = (v.len() - 1) as u32;
        let idx = |i: u32, j: u32| 1 + (i - 1) * n + (j % n);
        let mut f = Vec::new();
        for j in 0..n {
            f.push(vec![0, idx(1, j), idx(1, j + 1)]);
            f.push(vec![south, idx(n - 1, j + 1), idx(n - 1, j)]);
        }
        for i in 1..n - 1 {
            for j in 0..n {
                f.push(vec![idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)]);
                f.push(vec![idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)]);
            }
        }
        Polyhedron::new(v, f).unwrap()
    }

    /// Non-convex closed mesh: an L-shaped prism (extruded along z).
    fn l_prism() -> Polyhedron {
        let outline = [
            (0.0, 0.0),
            (3.0, 0.0),
            (3.0, 1.0),
            (1.0, 1.0),
            (1.0, 3.0),
            (0.0, 3.0),
        ];
        let n = outline.len() as u32;
        let mut v: Vec<Vec3> = outline.iter().map(|&(x, y)| Vec3::new(x, y, 0.0)).collect();
        v.extend(outline.iter().map(|&(x, y)| Vec3::new(x, y, 2.0)));
        let mut f = vec![(0..n).rev().collect::<Vec<u32>>(), (n..2 * n).collect()];
        for k in 0..n {
            let k1 = (k + 1) % n;
            f.push(vec![k, k1, n + k1, n + k]);
        }
        // The fan triangulation of the (non-convex) caps is wrong for an L
        // outline, so split the caps explicitly into convex parts.
        f.drain(..2);
        f.push(vec![0, 3, 2, 1]);
        f.push(vec![0, 5, 4, 3]);
        f.push(vec![n, n + 1, n + 2, n + 3]);
        f.push(vec![n, n + 3, n + 4, n + 5]);
        Polyhedron::new(v, f).unwrap()
    }

    fn brute_on_surface(p: &Polyhedron, q: Vec3) -> bool {
        let eps = p.surface_eps();
        p.tris.iter().any(|t| {
            let [a, b, c] = p.tri(t);
            (closest_on_triangle(q, a, b, c) - q).norm() <= eps
        })
    }

    fn cross_check(p: &Polyhedron, lo: f64, hi: f64, n: usize, seed: u64) {
        assert!(p.is_closed());
        let mut rng = Lcg(seed);
        for _ in 0..n {
            let q = Vec3::new(
                lo + (hi - lo) * rng.next(),
                lo + (hi - lo) * rng.next(),
                lo + (hi - lo) * rng.next(),
            );
            let w = p.winding_number(q).abs() > 0.5 || brute_on_surface(p, q);
            assert_eq!(p.contains(q), w, "point {q:?}");
            assert_eq!(p.on_surface(q), brute_on_surface(p, q), "point {q:?}");
        }
    }

    #[test]
    fn bvh_matches_winding_number() {
        cross_check(&sphere(24, 1.0), -1.2, 1.2, 3000, 1);
        cross_check(&unit_cube(), -0.2, 1.2, 3000, 2);
        cross_check(&l_prism(), -0.5, 3.5, 3000, 3);
    }

    #[test]
    fn grid_aligned_points() {
        // Points on a regular lattice hit edges/vertices of axis-aligned
        // meshes in the axis directions; the skewed rays must not care.
        let c = l_prism();
        for i in 0..=14 {
            for j in 0..=14 {
                for k in 0..=10 {
                    let q = Vec3::new(
                        i as f64 * 0.25 - 0.25,
                        j as f64 * 0.25 - 0.25,
                        k as f64 * 0.25 - 0.25,
                    );
                    let w = c.winding_number(q).abs() > 0.5 || brute_on_surface(&c, q);
                    assert_eq!(c.contains(q), w, "point {q:?}");
                }
            }
        }
    }

    #[test]
    fn surface_points() {
        let s = sphere(16, 2.0);
        let mut rng = Lcg(9);
        for t in s.triangles().iter().take(200) {
            let [a, b, c] = s.tri(t);
            let (u, v) = (rng.next() * 0.5, rng.next() * 0.5);
            let q = a + (b - a) * u + (c - a) * v;
            assert!(s.on_surface(q));
            assert!(s.contains(q));
        }
    }
}
