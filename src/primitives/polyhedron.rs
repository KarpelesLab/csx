//! Polyhedral surfaces and solids.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::math::{BoundingBox, Vec3};

use super::geom::REL_TOL;

/// A polyhedron defined by vertices and polygonal faces (vertex index
/// lists). Faces with more than three vertices are fan-triangulated.
///
/// A *closed* polyhedron (every edge shared by exactly two faces) is a solid
/// and its interior is determined with a generalised winding number; an
/// open polyhedron is a surface and only points on it are "inside".
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Polyhedron {
    vertices: Vec<Vec3>,
    faces: Vec<Vec<u32>>,
    tris: Vec<[u32; 3]>,
    bbox: BoundingBox,
    closed: bool,
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
        Ok(Polyhedron {
            vertices,
            faces,
            tris,
            bbox,
            closed,
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

    /// True if `p` lies on the surface (within tolerance).
    pub fn on_surface(&self, p: Vec3) -> bool {
        let scale = self
            .bbox
            .size()
            .norm()
            .max(self.bbox.max.norm())
            .max(self.bbox.min.norm());
        let eps = REL_TOL * scale.max(1e-300);
        self.tris.iter().any(|t| {
            let [a, b, c] = self.tri(t);
            (closest_on_triangle(p, a, b, c) - p).norm() <= eps
        })
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
        self.closed && self.winding_number(p).abs() > 0.5
    }
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
}
