//! Boxes, spheres, cylinders and their shells.

use crate::coords::CoordSystem;
use crate::error::{Error, Result};
use crate::math::{BoundingBox, Vec3};
use crate::param::{ParameterSet, Scalar, ScalarPoint};

use super::geom::{angle_ok, in_range, sector_bbox, tol};

/// Axis-aligned box between two corners (`P1`, `P2`). In a cylindrical
/// coordinate system the corners are `(ρ, α, z)` and the box is an annular
/// sector.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct BoxShape {
    /// First corner.
    pub start: ScalarPoint,
    /// Opposite corner.
    pub stop: ScalarPoint,
}

impl BoxShape {
    /// Box between two corners.
    pub fn new(start: impl Into<ScalarPoint>, stop: impl Into<ScalarPoint>) -> Self {
        BoxShape {
            start: start.into(),
            stop: stop.into(),
        }
    }

    pub(crate) fn inside(&self, p: Vec3, cs: CoordSystem) -> bool {
        inside_box(self.start.value(), self.stop.value(), p, cs)
    }

    pub(crate) fn bbox(&self, cs: CoordSystem) -> BoundingBox {
        box_bbox(self.start.value(), self.stop.value(), cs)
    }

    pub(crate) fn native(&self) -> BoundingBox {
        BoundingBox::from_corners(self.start.value(), self.stop.value())
    }

    pub(crate) fn dimension(&self) -> usize {
        self.native().dimension()
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.start.update(p)?;
        self.stop.update(p)
    }
}

/// Inside test for a box whose corners `a`, `b` are expressed in `cs`, for
/// the Cartesian point `p`.
pub(crate) fn inside_box(a: Vec3, b: Vec3, p: Vec3, cs: CoordSystem) -> bool {
    match cs {
        CoordSystem::Cartesian => (0..3).all(|i| in_range(p[i], a[i], b[i])),
        CoordSystem::Cylindrical => {
            let c = cs.from_cartesian(p);
            if !in_range(c[2], a[2], b[2]) || !in_range(c[0], a[0], b[0]) {
                return false;
            }
            // On the axis every angle is valid.
            c[0] <= tol(a[0], b[0]) || angle_ok(c[1], a[1], b[1])
        }
    }
}

pub(crate) fn box_bbox(a: Vec3, b: Vec3, cs: CoordSystem) -> BoundingBox {
    match cs {
        CoordSystem::Cartesian => BoundingBox::from_corners(a, b),
        CoordSystem::Cylindrical => sector_bbox((a[0], b[0]), (a[1], b[1]), (a[2], b[2])),
    }
}

/// A union of axis-aligned boxes.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MultiBox {
    /// Pairs of opposite corners.
    pub boxes: Vec<(ScalarPoint, ScalarPoint)>,
}

impl MultiBox {
    pub(crate) fn inside(&self, p: Vec3, cs: CoordSystem) -> bool {
        self.boxes
            .iter()
            .any(|(a, b)| inside_box(a.value(), b.value(), p, cs))
    }

    pub(crate) fn bbox(&self, cs: CoordSystem) -> BoundingBox {
        self.boxes.iter().fold(BoundingBox::EMPTY, |acc, (a, b)| {
            acc.union(&box_bbox(a.value(), b.value(), cs))
        })
    }

    pub(crate) fn native(&self) -> BoundingBox {
        self.boxes.iter().fold(BoundingBox::EMPTY, |acc, (a, b)| {
            acc.union(&BoundingBox::from_corners(a.value(), b.value()))
        })
    }

    pub(crate) fn dimension(&self) -> usize {
        self.boxes
            .iter()
            .map(|(a, b)| BoundingBox::from_corners(a.value(), b.value()).dimension())
            .max()
            .unwrap_or(0)
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        for (a, b) in &mut self.boxes {
            a.update(p)?;
            b.update(p)?;
        }
        Ok(())
    }
}

fn check_radius(r: &Scalar, what: &str) -> Result<()> {
    if r.value() < 0.0 {
        Err(Error::Invalid(format!("negative {what} {}", r.value())))
    } else {
        Ok(())
    }
}

/// Solid sphere. Like all primitives except boxes, its values are always
/// Cartesian.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Sphere {
    /// Center (in the primitive's coordinate system).
    pub center: ScalarPoint,
    /// Radius.
    pub radius: Scalar,
}

impl Sphere {
    /// Sphere from center and radius.
    pub fn new(center: impl Into<ScalarPoint>, radius: impl Into<Scalar>) -> Self {
        Sphere {
            center: center.into(),
            radius: radius.into(),
        }
    }

    pub(crate) fn inside(&self, p: Vec3) -> bool {
        let c = self.center.value();
        let r = self.radius.value();
        (p - c).norm() <= r + tol(r, 0.0)
    }

    pub(crate) fn bbox(&self) -> BoundingBox {
        let c = self.center.value();
        BoundingBox::from_corners(c, c).expanded(self.radius.value())
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.center.update(p)?;
        self.radius.update(p)?;
        check_radius(&self.radius, "radius")
    }
}

/// Spherical shell: points whose distance from the center is within
/// `shell_width / 2` of `radius`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SphericalShell {
    /// Center.
    pub center: ScalarPoint,
    /// Mid-surface radius.
    pub radius: Scalar,
    /// Total shell thickness.
    pub shell_width: Scalar,
}

impl SphericalShell {
    pub(crate) fn inside(&self, p: Vec3) -> bool {
        let c = self.center.value();
        let (r, w) = (self.radius.value(), self.shell_width.value() / 2.0);
        ((p - c).norm() - r).abs() <= w + tol(r, w)
    }

    pub(crate) fn bbox(&self) -> BoundingBox {
        let c = self.center.value();
        BoundingBox::from_corners(c, c)
            .expanded(self.radius.value() + self.shell_width.value() / 2.0)
    }

    pub(crate) fn dimension(&self) -> usize {
        if self.shell_width.value() > 0.0 { 3 } else { 2 }
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.center.update(p)?;
        self.radius.update(p)?;
        self.shell_width.update(p)?;
        check_radius(&self.radius, "radius")?;
        check_radius(&self.shell_width, "shell width")
    }
}

/// Axial coordinate `t ∈ [0, 1]` (unclamped) and radial distance of `p`
/// relative to the segment `a`–`b`.
fn axial(a: Vec3, b: Vec3, p: Vec3) -> Option<(f64, f64, f64)> {
    let d = b - a;
    let l2 = d.norm_sq();
    if l2 == 0.0 {
        return None;
    }
    let t = (p - a).dot(d) / l2;
    let r = (p - (a + d * t)).norm();
    Some((t, r, l2.sqrt()))
}

fn cylinder_bbox(a: Vec3, b: Vec3, r: f64) -> BoundingBox {
    let d = (b - a).normalized();
    let e = Vec3::new(
        r * (1.0 - d[0] * d[0]).max(0.0).sqrt(),
        r * (1.0 - d[1] * d[1]).max(0.0).sqrt(),
        r * (1.0 - d[2] * d[2]).max(0.0).sqrt(),
    );
    let bb = BoundingBox::from_corners(a, b);
    BoundingBox {
        min: bb.min - e,
        max: bb.max + e,
    }
}

/// Solid circular cylinder around the axis segment `start`–`stop`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Cylinder {
    /// Axis start point.
    pub start: ScalarPoint,
    /// Axis end point.
    pub stop: ScalarPoint,
    /// Radius.
    pub radius: Scalar,
}

impl Cylinder {
    /// Cylinder from axis end points and radius.
    pub fn new(
        start: impl Into<ScalarPoint>,
        stop: impl Into<ScalarPoint>,
        radius: impl Into<Scalar>,
    ) -> Self {
        Cylinder {
            start: start.into(),
            stop: stop.into(),
            radius: radius.into(),
        }
    }

    pub(crate) fn inside(&self, p: Vec3) -> bool {
        let (a, b) = (self.start.value(), self.stop.value());
        let r = self.radius.value();
        match axial(a, b, p) {
            Some((t, d, len)) => {
                let tt = tol(len, 0.0) / len;
                t >= -tt && t <= 1.0 + tt && d <= r + tol(r, 0.0)
            }
            None => false,
        }
    }

    pub(crate) fn bbox(&self) -> BoundingBox {
        let (a, b) = (self.start.value(), self.stop.value());
        cylinder_bbox(a, b, self.radius.value())
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.start.update(p)?;
        self.stop.update(p)?;
        self.radius.update(p)?;
        check_radius(&self.radius, "radius")
    }
}

/// Cylindrical shell: points within `shell_width / 2` of the cylinder
/// mantle of radius `radius`, between the two end planes.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CylindricalShell {
    /// Axis start point.
    pub start: ScalarPoint,
    /// Axis end point.
    pub stop: ScalarPoint,
    /// Mid-surface radius.
    pub radius: Scalar,
    /// Total shell thickness.
    pub shell_width: Scalar,
}

impl CylindricalShell {
    pub(crate) fn inside(&self, p: Vec3) -> bool {
        let (a, b) = (self.start.value(), self.stop.value());
        let (r, w) = (self.radius.value(), self.shell_width.value() / 2.0);
        match axial(a, b, p) {
            Some((t, d, len)) => {
                let tt = tol(len, 0.0) / len;
                t >= -tt && t <= 1.0 + tt && (d - r).abs() <= w + tol(r, w)
            }
            None => false,
        }
    }

    pub(crate) fn bbox(&self) -> BoundingBox {
        let (a, b) = (self.start.value(), self.stop.value());
        cylinder_bbox(a, b, self.radius.value() + self.shell_width.value() / 2.0)
    }

    pub(crate) fn dimension(&self) -> usize {
        if self.shell_width.value() > 0.0 { 3 } else { 2 }
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.start.update(p)?;
        self.stop.update(p)?;
        self.radius.update(p)?;
        self.shell_width.update(p)?;
        check_radius(&self.radius, "radius")?;
        check_radius(&self.shell_width, "shell width")
    }
}
