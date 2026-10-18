//! Planar polygons and the solids derived from them.
//!
//! A polygon lies in a plane normal to axis `norm_dir` (0 = x, 1 = y,
//! 2 = z). Its 2-D vertices `[u, v]` are coordinates along the two remaining
//! axes in cyclic order: `u` along `(norm_dir + 1) % 3` and `v` along
//! `(norm_dir + 2) % 3` (for `norm_dir = 2`: `u = x`, `v = y`).
//!
//! | `norm_dir` | `u` (`X1`) | `v` (`X2`) |
//! |---|---|---|
//! | 0 (x) | y | z |
//! | 1 (y) | z | x |
//! | 2 (z) | x | y |
//!
//! Polygon based primitives are always defined in Cartesian coordinates.
//! Polygons are implicitly closed and may be concave; points on an edge
//! count as inside.

use core::f64::consts::{PI, TAU};

use crate::error::{Error, Result};
use crate::math::{BoundingBox, Vec3};
use crate::param::{ParameterSet, Scalar};

use super::geom::{angle_ok, in_range, point_in_polygon, sector_bbox};

/// Polygon vertices: parametric `[u, v]` coordinates plus a cache of their
/// evaluated values.
///
/// The cache is what geometric queries use; it is rebuilt whenever the
/// vertices are modified through this type's methods and by
/// [`PolyVertices::update`] (called by the owning structure's `update()`).
/// `PolyVertices` dereferences to the evaluated `(u, v)` slice. Vertices given
/// as expressions evaluate to `NaN` until updated, unless they are
/// parameter-free.
///
/// ```
/// use csx::primitives::PolyVertices;
/// use csx::ParameterSet;
///
/// // Plain numbers convert directly.
/// let square: PolyVertices = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)].into();
/// assert_eq!(square.len(), 4);
/// assert_eq!(square[2], (1.0, 1.0));
///
/// // Parametric vertices keep their expressions.
/// let mut tri = PolyVertices::new();
/// tri.push(0.0, 0.0);
/// tri.push("w", 0.0);
/// tri.push("w", "w/2");
/// assert!(tri[1].0.is_nan()); // not evaluated yet
/// let mut ps = ParameterSet::new();
/// ps.set("w", 4.0);
/// tri.update(&ps)?;
/// assert_eq!(tri[2], (4.0, 2.0));
/// assert_eq!(tri.points()[2][1].expression(), Some("w/2"));
/// # Ok::<(), csx::Error>(())
/// ```
#[derive(Debug, Clone, Default)]
pub struct PolyVertices {
    points: Vec<[Scalar; 2]>,
    values: Vec<(f64, f64)>,
}

impl PartialEq for PolyVertices {
    fn eq(&self, o: &Self) -> bool {
        self.points == o.points
    }
}

impl PolyVertices {
    /// No vertices.
    pub fn new() -> Self {
        Self::default()
    }

    /// Vertices from parametric coordinates.
    pub fn from_scalars(points: impl IntoIterator<Item = [Scalar; 2]>) -> Self {
        let points: Vec<[Scalar; 2]> = points.into_iter().collect();
        let values = points.iter().map(|[a, b]| (a.value(), b.value())).collect();
        PolyVertices { points, values }
    }

    /// The parametric coordinates.
    pub fn points(&self) -> &[[Scalar; 2]] {
        &self.points
    }

    /// The evaluated coordinates (same as dereferencing).
    pub fn values(&self) -> &[(f64, f64)] {
        &self.values
    }

    /// Appends a vertex.
    pub fn push(&mut self, u: impl Into<Scalar>, v: impl Into<Scalar>) {
        let p = [u.into(), v.into()];
        self.values.push((p[0].value(), p[1].value()));
        self.points.push(p);
    }

    /// Replaces vertex `i`.
    ///
    /// # Panics
    ///
    /// If `i` is out of bounds (like slice indexing).
    pub fn set(&mut self, i: usize, u: impl Into<Scalar>, v: impl Into<Scalar>) {
        let p = [u.into(), v.into()];
        self.values[i] = (p[0].value(), p[1].value());
        self.points[i] = p;
    }

    /// Removes all vertices.
    pub fn clear(&mut self) {
        self.points.clear();
        self.values.clear();
    }

    /// Re-evaluates the parametric coordinates against `params`.
    ///
    /// # Errors
    ///
    /// [`Error::UnknownVariable`](crate::Error::UnknownVariable) or
    /// [`Error::Parse`](crate::Error::Parse) when a vertex expression cannot
    /// be evaluated, [`Error::Invalid`](crate::Error::Invalid) when it
    /// evaluates to NaN. Vertices before the failing one are already updated.
    pub fn update(&mut self, params: &ParameterSet) -> Result<()> {
        for (p, v) in self.points.iter_mut().zip(&mut self.values) {
            p[0].update(params)?;
            p[1].update(params)?;
            *v = (p[0].value(), p[1].value());
        }
        Ok(())
    }
}

impl core::ops::Deref for PolyVertices {
    type Target = [(f64, f64)];
    fn deref(&self) -> &[(f64, f64)] {
        &self.values
    }
}

impl<'a> IntoIterator for &'a PolyVertices {
    type Item = &'a (f64, f64);
    type IntoIter = core::slice::Iter<'a, (f64, f64)>;
    fn into_iter(self) -> Self::IntoIter {
        self.values.iter()
    }
}

impl From<Vec<(f64, f64)>> for PolyVertices {
    fn from(v: Vec<(f64, f64)>) -> Self {
        let points = v
            .iter()
            .map(|&(a, b)| [Scalar::new(a), Scalar::new(b)])
            .collect();
        PolyVertices { points, values: v }
    }
}

impl FromIterator<[Scalar; 2]> for PolyVertices {
    fn from_iter<I: IntoIterator<Item = [Scalar; 2]>>(it: I) -> Self {
        PolyVertices::from_scalars(it)
    }
}

impl FromIterator<(f64, f64)> for PolyVertices {
    fn from_iter<I: IntoIterator<Item = (f64, f64)>>(it: I) -> Self {
        PolyVertices::from(it.into_iter().collect::<Vec<_>>())
    }
}

/// The two in-plane axes `(u, v)` for a plane normal to axis `n` (0..=2):
/// `((n + 1) % 3, (n + 2) % 3)`.
///
/// ```
/// assert_eq!(csx::primitives::plane_axes(1), (2, 0)); // normal y: u = z, v = x
/// ```
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
    if n > 2 {
        return BoundingBox::EMPTY;
    }
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

/// A planar polygon (2-D sheet) in the plane `axis[norm_dir] = elevation`.
///
/// ```
/// use csx::primitives::{Polygon, Primitive};
/// use csx::{CoordSystem, Vec3};
///
/// // An L-shaped sheet at z = 0.5 (normal z: vertices are (x, y)).
/// let l = vec![(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0)];
/// let p = Primitive::new(Polygon::new(2, 0.5, l));
/// assert_eq!(p.dimension(), 2);
/// assert!(p.contains(Vec3::new(0.5, 1.5, 0.5), CoordSystem::Cartesian));
/// assert!(!p.contains(Vec3::new(1.5, 1.5, 0.5), CoordSystem::Cartesian)); // the notch
/// assert!(!p.contains(Vec3::new(0.5, 1.5, 0.6), CoordSystem::Cartesian)); // off-plane
/// ```
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Polygon {
    /// Normal axis (0, 1, 2).
    pub norm_dir: usize,
    /// Position of the plane along `norm_dir`.
    pub elevation: Scalar,
    /// Vertices `(u, v)`; see the module docs for the axis mapping.
    pub vertices: PolyVertices,
}

impl Polygon {
    /// Polygon in the plane `axis[norm_dir] = elevation`.
    pub fn new(
        norm_dir: usize,
        elevation: impl Into<Scalar>,
        vertices: impl Into<PolyVertices>,
    ) -> Self {
        Polygon {
            norm_dir,
            elevation: elevation.into(),
            vertices: vertices.into(),
        }
    }

    pub(crate) fn inside(&self, p: Vec3) -> bool {
        let n = self.norm_dir;
        if n > 2 {
            return false;
        }
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
        self.vertices.update(p)?;
        self.elevation.update(p)
    }
}

/// A polygon linearly extruded along its normal axis from `elevation` to
/// `elevation + length` (`LinPoly` in CSXCAD files). A zero length makes it
/// a sheet.
///
/// ```
/// use csx::primitives::{LinPoly, Primitive};
/// use csx::{CoordSystem, Vec3};
///
/// // A triangle in the x/y plane extruded from z = 0 to z = 1.6.
/// let p = Primitive::new(LinPoly::new(2, 0.0, 1.6, vec![(0.0, 0.0), (4.0, 0.0), (0.0, 3.0)]));
/// assert_eq!(p.dimension(), 3);
/// assert!(p.contains(Vec3::new(1.0, 1.0, 0.8), CoordSystem::Cartesian));
/// assert!(!p.contains(Vec3::new(3.0, 3.0, 0.8), CoordSystem::Cartesian));
/// assert_eq!(p.bounding_box(CoordSystem::Cartesian).max, Vec3::new(4.0, 3.0, 1.6));
/// ```
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LinPoly {
    /// Normal (extrusion) axis.
    pub norm_dir: usize,
    /// Start of the extrusion along `norm_dir`.
    pub elevation: Scalar,
    /// Extrusion length (may be negative).
    pub length: Scalar,
    /// Vertices `(u, v)`.
    pub vertices: PolyVertices,
}

impl LinPoly {
    /// Extruded polygon.
    pub fn new(
        norm_dir: usize,
        elevation: impl Into<Scalar>,
        length: impl Into<Scalar>,
        vertices: impl Into<PolyVertices>,
    ) -> Self {
        LinPoly {
            norm_dir,
            elevation: elevation.into(),
            length: length.into(),
            vertices: vertices.into(),
        }
    }

    pub(crate) fn inside(&self, p: Vec3) -> bool {
        let n = self.norm_dir;
        if n > 2 {
            return false;
        }
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
        self.vertices.update(p)?;
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
/// `0` is the polygon's own position (radial direction `+q`). Both
/// orientations of the angle range describe the same sweep; a range of 2π
/// or more is a full revolution. The XML form is `<RotPoly NormDir
/// RotAxisDir>` with `<Angles Start Stop/>` (verified against CSXCAD).
///
/// ```
/// use csx::primitives::{Primitive, RotPoly};
/// use csx::{CoordSystem, Vec3};
///
/// // A washer: the square ρ ∈ [1, 2], z ∈ [0, 1] in the x/z plane (normal y,
/// // vertices are (z, x)) revolved about the z axis.
/// let washer = RotPoly::new(1, 2, vec![(0.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0)]);
/// let p = Primitive::new(washer);
/// assert!(p.contains(Vec3::new(0.0, -1.5, 0.5), CoordSystem::Cartesian));
/// assert!(!p.contains(Vec3::new(0.5, 0.0, 0.5), CoordSystem::Cartesian)); // the hole
/// assert_eq!(p.bounding_box(CoordSystem::Cartesian).max, Vec3::new(2.0, 2.0, 1.0));
/// ```
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
    pub vertices: PolyVertices,
}

impl Default for RotPoly {
    fn default() -> Self {
        RotPoly {
            norm_dir: 0,
            rot_axis_dir: 2,
            elevation: Scalar::new(0.0),
            start_angle: Scalar::new(0.0),
            stop_angle: Scalar::new(TAU),
            vertices: PolyVertices::new(),
        }
    }
}

impl RotPoly {
    /// Full (2π) revolution of a polygon.
    pub fn new(norm_dir: usize, rot_axis_dir: usize, vertices: impl Into<PolyVertices>) -> Self {
        RotPoly {
            norm_dir,
            rot_axis_dir,
            vertices: vertices.into(),
            ..Default::default()
        }
    }

    /// True if the axes are valid (both in 0..=2 and distinct).
    fn axes_valid(&self) -> bool {
        self.norm_dir < 3 && self.rot_axis_dir < 3 && self.norm_dir != self.rot_axis_dir
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
        if !self.axes_valid() {
            return false;
        }
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
        if !self.axes_valid() {
            return BoundingBox::EMPTY;
        }
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
        self.vertices.update(p)?;
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

    #[test]
    fn parametric_vertices() {
        let mut v = PolyVertices::new();
        v.push(0.0, 0.0);
        v.push("w", 0.0);
        v.push("w", "w/2");
        v.push(0.0, "2*3");
        // parameter-free expressions evaluate immediately
        assert_eq!(v[3], (0.0, 6.0));
        assert!(v[1].0.is_nan());
        let mut p = Polygon::new(2, 0.0, v);
        let mut ps = ParameterSet::new();
        ps.set("w", 2.0);
        p.update(&ps).unwrap();
        assert_eq!(p.vertices.values()[2], (2.0, 1.0));
        assert!(p.inside(Vec3::new(1.5, 0.5, 0.0)));
        assert!(!p.inside(Vec3::new(3.0, 0.5, 0.0)));
        ps.set("w", 4.0);
        p.update(&ps).unwrap();
        assert!(p.inside(Vec3::new(3.0, 0.5, 0.0)));
        assert_eq!(p.vertices.points()[1][0].expression(), Some("w"));
        // equality ignores the cache, compares definitions
        let q = Polygon::new(2, 0.0, p.vertices.clone());
        assert_eq!(p, q);
        p.vertices.set(0, -1.0, -1.0);
        assert_eq!(p.vertices[0], (-1.0, -1.0));
        assert!(p.update(&ParameterSet::new()).is_err());
    }
}
