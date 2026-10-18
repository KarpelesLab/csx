//! Planar polygons and the solids derived from them.
//!
//! A polygon lies in a plane normal to axis `norm_dir` (0 = x, 1 = y,
//! 2 = z). Its 2-D vertices `[u, v]` are coordinates along the two remaining
//! axes in cyclic order: `u` along `(norm_dir + 1) % 3` and `v` along
//! `(norm_dir + 2) % 3` (for `norm_dir = 2`: `u = x`, `v = y`).
//!
//! Polygon based primitives are always defined in Cartesian coordinates.

use core::f64::consts::{PI, TAU};

use crate::error::{Error, Result};
use crate::math::{BoundingBox, Vec3};
use crate::param::{ParameterSet, Scalar};

use super::geom::{angle_ok, in_range, point_in_polygon, sector_bbox};

/// The two in-plane axes for a plane normal to `n`.
#[inline]
pub fn plane_axes(n: usize) -> (usize, usize) {
    ((n + 1) % 3, (n + 2) % 3)
}

fn check_dir(d: usize, what: &str) -> Result<()> {
    if d < 3 {
        Ok(())
    } else {
        Err(Error::Invalid(format!("{what} {d} is not an axis (0..=2)")))
    }
}

fn verts_bbox(n: usize, verts: &[(f64, f64)], lo: f64, hi: f64) -> BoundingBox {
    let (u, v) = plane_axes(n);
    let mut b = BoundingBox::EMPTY;
    for &(a, c) in verts {
        for h in [lo, hi] {
            let mut p = Vec3::ZERO;
            p[u] = a;
            p[v] = c;
            p[n] = h;
            b.include(p);
        }
    }
    b
}

/// A planar polygon (2-D sheet).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Polygon {
    /// Normal axis (0, 1, 2).
    pub norm_dir: usize,
    /// Position of the plane along `norm_dir`.
    pub elevation: Scalar,
    /// Vertices `(u, v)`; see the module docs for the axis mapping.
    pub vertices: Vec<(f64, f64)>,
}

impl Polygon {
    /// Polygon in the plane `axis[norm_dir] = elevation`.
    pub fn new(norm_dir: usize, elevation: impl Into<Scalar>, vertices: Vec<(f64, f64)>) -> Self {
        Polygon {
            norm_dir,
            elevation: elevation.into(),
            vertices,
        }
    }

    pub(crate) fn inside(&self, p: Vec3) -> bool {
        let n = self.norm_dir;
        let (u, v) = plane_axes(n);
        let e = self.elevation.value();
        in_range(p[n], e, e) && point_in_polygon((p[u], p[v]), &self.vertices)
    }

    fn native(&self) -> BoundingBox {
        let e = self.elevation.value();
        verts_bbox(self.norm_dir, &self.vertices, e, e)
    }

    pub(crate) fn bbox(&self) -> BoundingBox {
        self.native()
    }

    pub(crate) fn dimension(&self) -> usize {
        if self.vertices.len() >= 3 {
            2
        } else {
            self.native().dimension()
        }
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        check_dir(self.norm_dir, "normal direction")?;
        self.elevation.update(p)
    }
}

/// A polygon linearly extruded along its normal axis from `elevation` to
/// `elevation + length`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LinPoly {
    /// Normal (extrusion) axis.
    pub norm_dir: usize,
    /// Start of the extrusion along `norm_dir`.
    pub elevation: Scalar,
    /// Extrusion length (may be negative).
    pub length: Scalar,
    /// Vertices `(u, v)`.
    pub vertices: Vec<(f64, f64)>,
}

impl LinPoly {
    /// Extruded polygon.
    pub fn new(
        norm_dir: usize,
        elevation: impl Into<Scalar>,
        length: impl Into<Scalar>,
        vertices: Vec<(f64, f64)>,
    ) -> Self {
        LinPoly {
            norm_dir,
            elevation: elevation.into(),
            length: length.into(),
            vertices,
        }
    }

    pub(crate) fn inside(&self, p: Vec3) -> bool {
        let n = self.norm_dir;
        let (u, v) = plane_axes(n);
        let e = self.elevation.value();
        let l = self.length.value();
        in_range(p[n], e, e + l) && point_in_polygon((p[u], p[v]), &self.vertices)
    }

    fn native(&self) -> BoundingBox {
        let e = self.elevation.value();
        verts_bbox(self.norm_dir, &self.vertices, e, e + self.length.value())
    }

    pub(crate) fn bbox(&self) -> BoundingBox {
        self.native()
    }

    pub(crate) fn dimension(&self) -> usize {
        if self.length.value() != 0.0 { 3 } else { 2 }
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        check_dir(self.norm_dir, "normal direction")?;
        self.elevation.update(p)?;
        self.length.update(p)
    }
}

/// A polygon revolved about a coordinate axis lying in its plane.
///
/// The polygon lies in the plane normal to `norm_dir` that contains the
/// rotation axis `rot_axis_dir`. Let `q` be the third axis (the *radial*
/// axis). The polygon is swept by a right-handed rotation about
/// `rot_axis_dir` from `start_angle` to `stop_angle` (radians), where angle
/// `0` is the polygon's own position (radial direction `+q`).
#[derive(Debug, Clone, PartialEq)]
pub struct RotPoly {
    /// Normal axis of the polygon plane.
    pub norm_dir: usize,
    /// Rotation axis (must differ from `norm_dir`).
    pub rot_axis_dir: usize,
    /// Plane offset (stored for round-tripping; the rotation axis passes
    /// through the origin).
    pub elevation: Scalar,
    /// Start angle (radians).
    pub start_angle: Scalar,
    /// Stop angle (radians).
    pub stop_angle: Scalar,
    /// Vertices `(u, v)`.
    pub vertices: Vec<(f64, f64)>,
}

impl Default for RotPoly {
    fn default() -> Self {
        RotPoly {
            norm_dir: 0,
            rot_axis_dir: 2,
            elevation: Scalar::new(0.0),
            start_angle: Scalar::new(0.0),
            stop_angle: Scalar::new(TAU),
            vertices: Vec::new(),
        }
    }
}

impl RotPoly {
    /// Full (2π) revolution of a polygon.
    pub fn new(norm_dir: usize, rot_axis_dir: usize, vertices: Vec<(f64, f64)>) -> Self {
        RotPoly {
            norm_dir,
            rot_axis_dir,
            vertices,
            ..Default::default()
        }
    }

    /// Radial axis and the unit vector of positive rotation from it.
    fn frame(&self) -> (usize, Vec3) {
        let a = self.rot_axis_dir;
        let q = 3 - a - self.norm_dir;
        (q, Vec3::unit(a).cross(Vec3::unit(q)))
    }

    fn to_poly(&self, radial: f64, axial: f64, q: usize) -> (f64, f64) {
        let (u, _) = plane_axes(self.norm_dir);
        if u == q {
            (radial, axial)
        } else {
            (axial, radial)
        }
    }

    pub(crate) fn inside(&self, p: Vec3) -> bool {
        let a = self.rot_axis_dir;
        let (q, w) = self.frame();
        let h = p[a];
        let x = p[q];
        let y = p.dot(w);
        let r = x.hypot(y);
        let phi = y.atan2(x);
        let (s, e) = (self.start_angle.value(), self.stop_angle.value());
        if r == 0.0 {
            return point_in_polygon(self.to_poly(0.0, h, q), &self.vertices);
        }
        (angle_ok(phi, s, e) && point_in_polygon(self.to_poly(r, h, q), &self.vertices))
            || (angle_ok(phi + PI, s, e)
                && point_in_polygon(self.to_poly(-r, h, q), &self.vertices))
    }

    pub(crate) fn bbox(&self) -> BoundingBox {
        let a = self.rot_axis_dir;
        let (q, w) = self.frame();
        let (u, _) = plane_axes(self.norm_dir);
        let (mut rpos, mut rneg) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
        let (mut hmin, mut hmax) = (f64::INFINITY, f64::NEG_INFINITY);
        for &(c0, c1) in &self.vertices {
            let (r, h) = if u == q { (c0, c1) } else { (c1, c0) };
            rpos = rpos.max(r);
            rneg = rneg.max(-r);
            hmin = hmin.min(h);
            hmax = hmax.max(h);
        }
        if self.vertices.is_empty() {
            return BoundingBox::EMPTY;
        }
        let (s, e) = (self.start_angle.value(), self.stop_angle.value());
        let mut local = BoundingBox::EMPTY;
        if rpos > 0.0 {
            local = local.union(&sector_bbox((0.0, rpos), (s, e), (hmin, hmax)));
        }
        if rneg > 0.0 {
            local = local.union(&sector_bbox((0.0, rneg), (s + PI, e + PI), (hmin, hmax)));
        }
        if local.is_empty() {
            // Polygon degenerate on the axis.
            local = BoundingBox::from_corners(Vec3::new(0.0, 0.0, hmin), Vec3::new(0.0, 0.0, hmax));
        }
        // local (x', y', z') = (q, w, a) → global.
        BoundingBox::from_points(
            local
                .corners()
                .map(|c| Vec3::unit(q) * c[0] + w * c[1] + Vec3::unit(a) * c[2]),
        )
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        check_dir(self.norm_dir, "normal direction")?;
        check_dir(self.rot_axis_dir, "rotation axis")?;
        if self.norm_dir == self.rot_axis_dir {
            return Err(Error::Invalid(
                "rotation axis must lie in the polygon plane".into(),
            ));
        }
        self.elevation.update(p)?;
        self.start_angle.update(p)?;
        self.stop_angle.update(p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square() -> Vec<(f64, f64)> {
        vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    }

    #[test]
    fn polygon_axes() {
        // norm z: u = x, v = y
        let p = Polygon::new(2, 0.5, vec![(0.0, 0.0), (2.0, 0.0), (2.0, 1.0)]);
        assert!(p.inside(Vec3::new(1.5, 0.5, 0.5)));
        assert!(!p.inside(Vec3::new(0.5, 1.5, 0.5)));
        assert!(!p.inside(Vec3::new(1.5, 0.5, 0.6)));
        // norm x: u = y, v = z
        let p = Polygon::new(0, 0.0, vec![(0.0, 0.0), (2.0, 0.0), (2.0, 1.0)]);
        assert!(p.inside(Vec3::new(0.0, 1.5, 0.5)));
        assert_eq!(p.dimension(), 2);
    }

    #[test]
    fn linpoly() {
        let l = LinPoly::new(1, 1.0, -2.0, square());
        // norm y: u = z, v = x
        assert!(l.inside(Vec3::new(0.25, 0.0, 0.75)));
        assert!(!l.inside(Vec3::new(0.25, 1.5, 0.75)));
        let b = l.bbox();
        assert_eq!(b.min, Vec3::new(0.0, -1.0, 0.0));
        assert_eq!(b.max, Vec3::new(1.0, 1.0, 1.0));
    }

    #[test]
    fn rotpoly_torus_like() {
        // Square [1,2]x[0,1] in the x/z plane (normal y), revolved around z:
        // a washer of inner radius 1, outer radius 2, height 1.
        // norm y: u = z, v = x  → vertices given as (z, x).
        let r = RotPoly::new(1, 2, vec![(0.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0)]);
        for ang in [0.0f64, 1.0, 2.5, -2.0] {
            let p = Vec3::new(1.5 * ang.cos(), 1.5 * ang.sin(), 0.5);
            assert!(r.inside(p), "{ang}");
        }
        assert!(!r.inside(Vec3::new(0.5, 0.0, 0.5)));
        let b = r.bbox();
        assert!((b.min - Vec3::new(-2.0, -2.0, 0.0)).norm() < 1e-12);
        assert!((b.max - Vec3::new(2.0, 2.0, 1.0)).norm() < 1e-12);

        let mut half = r.clone();
        half.stop_angle = Scalar::new(PI / 2.0);
        assert!(half.inside(Vec3::new(1.0, 1.0, 0.5)));
        assert!(!half.inside(Vec3::new(-1.0, 1.0, 0.5)));
    }
}
