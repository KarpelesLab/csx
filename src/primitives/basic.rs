//! Boxes, spheres, cylinders and their shells.

use crate::coords::CoordSystem;
use crate::error::{Error, Result};
use crate::math::{BoundingBox, Vec3};
use crate::param::{ParameterSet, Scalar, ScalarPoint};

use super::geom::{angle_ok, in_range, sector_bbox, tol};

/// Axis-aligned box between two corners (`<Box>` with `P1`, `P2` in CSXCAD
/// files). The corners may be given in any order; faces count as inside.
///
/// In a cylindrical coordinate system the corners are `(ρ, α, z)` and the
/// box is an annular sector (see the
/// [module docs](super#coordinate-systems)). A box with a zero extent along
/// one or two axes is a sheet or a line ([`Shape::dimension`](super::Shape::dimension)).
///
/// ```
/// use csx::primitives::{BoxShape, Primitive};
/// use csx::{CoordSystem, ScalarPoint, Vec3};
///
/// let b = Primitive::new(BoxShape::new([2.0, 0.0, 0.0], [0.0, 1.0, 1.0]));
/// assert!(b.contains(Vec3::new(1.0, 0.5, 0.5), CoordSystem::Cartesian));
/// assert!(b.contains(Vec3::new(2.0, 1.0, 1.0), CoordSystem::Cartesian)); // corner
/// assert_eq!(b.bounding_box(CoordSystem::Cartesian).max, Vec3::new(2.0, 1.0, 1.0));
///
/// // Corners may be expressions over parameters (evaluated on update).
/// let p = BoxShape::new([0.0, 0.0, 0.0], ScalarPoint::from(["w", "w/2", "h"]));
/// assert_eq!(p.stop.0[0].expression(), Some("w"));
/// ```
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

    /// Inside test for a point already expressed in `cs`.
    pub(crate) fn inside_native(&self, c: Vec3, cs: CoordSystem) -> bool {
        inside_box_native(self.start.value(), self.stop.value(), c, cs)
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
    inside_box_native(a, b, cs.from_cartesian(p), cs)
}

/// Inside test for a box whose corners `a`, `b` are expressed in `cs`, for a
/// point `c` given in the same coordinate system. In cylindrical
/// coordinates the radius is compared as given (a negative `ρ` is outside)
/// and the angle modulo 2π.
pub(crate) fn inside_box_native(a: Vec3, b: Vec3, c: Vec3, cs: CoordSystem) -> bool {
    match cs {
        CoordSystem::Cartesian => (0..3).all(|i| in_range(c[i], a[i], b[i])),
        CoordSystem::Cylindrical => {
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
///
/// Unlike [`BoxShape`], a multi-box ignores the primitive's own coordinate
/// system (as in CSXCAD): its corners are in the *structure's* coordinate
/// system and are compared component-wise with the query coordinates,
/// without any angle wrapping (in a cylindrical structure a multi-box is a
/// set of plain `(ρ, α, z)` ranges).
///
/// ```
/// use csx::primitives::{MultiBox, Primitive};
/// use csx::{CoordSystem, Vec3};
///
/// let mut m = MultiBox::new([([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])]);
/// m.push([3.0, 0.0, 0.0], [4.0, 1.0, 1.0]);
/// let p = Primitive::new(m);
/// assert!(p.contains(Vec3::new(3.5, 0.5, 0.5), CoordSystem::Cartesian));
/// assert!(!p.contains(Vec3::new(2.0, 0.5, 0.5), CoordSystem::Cartesian)); // the gap
/// assert_eq!(p.bounding_box(CoordSystem::Cartesian).max, Vec3::new(4.0, 1.0, 1.0));
/// ```
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MultiBox {
    /// Pairs of opposite corners.
    pub boxes: Vec<(ScalarPoint, ScalarPoint)>,
}

impl MultiBox {
    /// Union of the boxes given by pairs of opposite corners.
    pub fn new<A: Into<ScalarPoint>, B: Into<ScalarPoint>>(
        boxes: impl IntoIterator<Item = (A, B)>,
    ) -> Self {
        MultiBox {
            boxes: boxes
                .into_iter()
                .map(|(a, b)| (a.into(), b.into()))
                .collect(),
        }
    }

    /// Adds a box.
    pub fn push(&mut self, start: impl Into<ScalarPoint>, stop: impl Into<ScalarPoint>) {
        self.boxes.push((start.into(), stop.into()));
    }

    /// Inside test for the Cartesian point `p`, corners in `cs`.
    pub(crate) fn inside(&self, p: Vec3, cs: CoordSystem) -> bool {
        self.inside_native(cs.from_cartesian(p))
    }

    /// Inside test for a point expressed in the corners' coordinate system.
    pub(crate) fn inside_native(&self, c: Vec3) -> bool {
        self.boxes.iter().any(|(a, b)| {
            let (a, b) = (a.value(), b.value());
            (0..3).all(|i| in_range(c[i], a[i], b[i]))
        })
    }

    pub(crate) fn bbox(&self, cs: CoordSystem) -> BoundingBox {
        self.boxes.iter().fold(BoundingBox::EMPTY, |acc, (a, b)| {
            acc.union(&box_bbox(a.value(), b.value(), cs))
        })
    }

    /// Bounds in the corners' own coordinate system.
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

/// Solid sphere. In a cylindrical coordinate system the center is given as
/// `(ρ, α, z)` and converted; the shape itself stays a Euclidean sphere.
///
/// ```
/// use csx::primitives::{Primitive, Sphere};
/// use csx::{CoordSystem, Vec3};
///
/// let s = Primitive::new(Sphere::new([1.0, 2.0, 3.0], 0.5));
/// assert!(s.contains(Vec3::new(1.0, 2.0, 3.4), CoordSystem::Cartesian));
/// assert!(!s.contains(Vec3::new(1.0, 2.0, 3.6), CoordSystem::Cartesian));
/// let b = s.bounding_box(CoordSystem::Cartesian);
/// assert_eq!((b.min, b.max), (Vec3::new(0.5, 1.5, 2.5), Vec3::new(1.5, 2.5, 3.5)));
/// ```
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

    pub(crate) fn inside(&self, p: Vec3, cs: CoordSystem) -> bool {
        let c = cs.to_cartesian(self.center.value());
        let r = self.radius.value();
        (p - c).norm() <= r + tol(r, 0.0)
    }

    pub(crate) fn bbox(&self, cs: CoordSystem) -> BoundingBox {
        let c = cs.to_cartesian(self.center.value());
        BoundingBox::from_corners(c, c).expanded(self.radius.value())
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.center.update(p)?;
        self.radius.update(p)?;
        check_radius(&self.radius, "radius")
    }
}

/// Spherical shell: points whose distance from the center is within
/// `|shell_width| / 2` of `radius` (the sign of the width is ignored). The
/// inner radius is `radius − width/2`, the outer `radius + width/2`. A zero
/// width makes it a surface (dimension 2).
///
/// ```
/// use csx::primitives::{Primitive, SphericalShell};
/// use csx::{CoordSystem, Vec3};
///
/// let s = Primitive::new(SphericalShell::new([0.0, 0.0, 0.0], 1.0, 0.2));
/// assert!(s.contains(Vec3::new(0.0, 0.0, 1.05), CoordSystem::Cartesian));
/// assert!(!s.contains(Vec3::new(0.0, 0.0, 0.5), CoordSystem::Cartesian)); // hollow
/// assert_eq!(s.bounding_box(CoordSystem::Cartesian).max, Vec3::new(1.1, 1.1, 1.1));
/// ```
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
    /// Shell around `center` with mid-surface radius `radius` and total
    /// thickness `shell_width`.
    pub fn new(
        center: impl Into<ScalarPoint>,
        radius: impl Into<Scalar>,
        shell_width: impl Into<Scalar>,
    ) -> Self {
        SphericalShell {
            center: center.into(),
            radius: radius.into(),
            shell_width: shell_width.into(),
        }
    }

    pub(crate) fn inside(&self, p: Vec3, cs: CoordSystem) -> bool {
        let c = cs.to_cartesian(self.center.value());
        let (r, w) = (self.radius.value(), self.shell_width.value().abs() / 2.0);
        ((p - c).norm() - r).abs() <= w + tol(r, w)
    }

    pub(crate) fn bbox(&self, cs: CoordSystem) -> BoundingBox {
        let c = cs.to_cartesian(self.center.value());
        BoundingBox::from_corners(c, c)
            .expanded(self.radius.value() + self.shell_width.value().abs() / 2.0)
    }

    pub(crate) fn dimension(&self) -> usize {
        if self.shell_width.value() != 0.0 {
            3
        } else {
            2
        }
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.center.update(p)?;
        self.radius.update(p)?;
        self.shell_width.update(p)?;
        check_radius(&self.radius, "radius")
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

/// Solid circular cylinder around the axis segment `start`–`stop` (the
/// centres of the two end faces). The axis may point in any direction. In a
/// cylindrical coordinate system the end points are given as `(ρ, α, z)`.
///
/// ```
/// use csx::primitives::{Cylinder, Primitive};
/// use csx::{CoordSystem, Vec3};
///
/// // A diagonal rod of radius 0.1.
/// let c = Primitive::new(Cylinder::new([0.0, 0.0, 0.0], [1.0, 1.0, 0.0], 0.1));
/// assert!(c.contains(Vec3::new(0.5, 0.5, 0.05), CoordSystem::Cartesian));
/// assert!(!c.contains(Vec3::new(0.5, 0.5, 0.2), CoordSystem::Cartesian));
/// assert!(!c.contains(Vec3::new(1.2, 1.2, 0.0), CoordSystem::Cartesian)); // past the end
/// ```
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

    pub(crate) fn inside(&self, p: Vec3, cs: CoordSystem) -> bool {
        let (a, b) = (
            cs.to_cartesian(self.start.value()),
            cs.to_cartesian(self.stop.value()),
        );
        let r = self.radius.value();
        match axial(a, b, p) {
            Some((t, d, len)) => {
                let tt = tol(len, 0.0) / len;
                t >= -tt && t <= 1.0 + tt && d <= r + tol(r, 0.0)
            }
            None => false,
        }
    }

    pub(crate) fn bbox(&self, cs: CoordSystem) -> BoundingBox {
        let (a, b) = (
            cs.to_cartesian(self.start.value()),
            cs.to_cartesian(self.stop.value()),
        );
        cylinder_bbox(a, b, self.radius.value())
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.start.update(p)?;
        self.stop.update(p)?;
        self.radius.update(p)?;
        check_radius(&self.radius, "radius")
    }
}

/// Cylindrical shell: points within `|shell_width| / 2` of the cylinder
/// mantle of radius `radius`, between the two end planes (a tube; the ends
/// are open). A zero width makes it a surface (dimension 2).
///
/// ```
/// use csx::primitives::{CylindricalShell, Primitive};
/// use csx::{CoordSystem, Vec3};
///
/// // Outer conductor of a coax: radius 2.25 ± 0.25 along z.
/// let tube = Primitive::new(CylindricalShell::new([0.0, 0.0, 0.0], [0.0, 0.0, 10.0], 2.25, 0.5));
/// assert!(tube.contains(Vec3::new(2.4, 0.0, 5.0), CoordSystem::Cartesian));
/// assert!(!tube.contains(Vec3::new(1.0, 0.0, 5.0), CoordSystem::Cartesian));
/// ```
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
    /// Shell around the axis `start`–`stop` with mid-surface radius
    /// `radius` and total thickness `shell_width`.
    pub fn new(
        start: impl Into<ScalarPoint>,
        stop: impl Into<ScalarPoint>,
        radius: impl Into<Scalar>,
        shell_width: impl Into<Scalar>,
    ) -> Self {
        CylindricalShell {
            start: start.into(),
            stop: stop.into(),
            radius: radius.into(),
            shell_width: shell_width.into(),
        }
    }

    pub(crate) fn inside(&self, p: Vec3, cs: CoordSystem) -> bool {
        let (a, b) = (
            cs.to_cartesian(self.start.value()),
            cs.to_cartesian(self.stop.value()),
        );
        let (r, w) = (self.radius.value(), self.shell_width.value().abs() / 2.0);
        match axial(a, b, p) {
            Some((t, d, len)) => {
                let tt = tol(len, 0.0) / len;
                t >= -tt && t <= 1.0 + tt && (d - r).abs() <= w + tol(r, w)
            }
            None => false,
        }
    }

    pub(crate) fn bbox(&self, cs: CoordSystem) -> BoundingBox {
        let (a, b) = (
            cs.to_cartesian(self.start.value()),
            cs.to_cartesian(self.stop.value()),
        );
        cylinder_bbox(
            a,
            b,
            self.radius.value() + self.shell_width.value().abs() / 2.0,
        )
    }

    pub(crate) fn dimension(&self) -> usize {
        if self.shell_width.value() != 0.0 {
            3
        } else {
            2
        }
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.start.update(p)?;
        self.stop.update(p)?;
        self.radius.update(p)?;
        self.shell_width.update(p)?;
        check_radius(&self.radius, "radius")
    }
}
