//! Coordinate systems and conversions.
//!
//! CSXCAD knows two coordinate systems: Cartesian `(x, y, z)` and cylindrical
//! `(ρ, α, z)` where `α` is the azimuth in radians measured from the +x axis
//! towards +y, and `z` is shared with the Cartesian system. Both are stored
//! in a [`Vec3`]; the meaning of its components depends on the system.
//!
//! # Crate-wide conventions
//!
//! * **Axis indices** are `usize`: `0` = x (or ρ), `1` = y (or α),
//!   `2` = z. They index [`Vec3`] components, grid line sets, material axes,
//!   polygon normals and so on. [`crate::Axis`] is the typed equivalent used
//!   by a few property fields.
//! * **Lengths** are in *drawing units*; the mesh's
//!   [`crate::RectilinearGrid::delta_unit`] gives the size of one unit in
//!   metres. Material parameters are in SI units.
//! * **Angles** are radians everywhere (cylindrical `α`, rotation
//!   arguments, revolved polygon angles), except for transforms that opt into
//!   degrees with [`crate::Transform::set_angles_in_degrees`].
//! * **Queries** always state the coordinate system of the query point
//!   explicitly; nothing is converted implicitly.
//!
//! ```
//! use core::f64::consts::FRAC_PI_2;
//! use csx::{CoordSystem, Vec3};
//!
//! // ρ = 2, α = 90°, z = 1 is the Cartesian point (0, 2, 1).
//! let cyl = Vec3::new(2.0, FRAC_PI_2, 1.0);
//! let cart = CoordSystem::Cylindrical.to_cartesian(cyl);
//! assert!((cart - Vec3::new(0.0, 2.0, 1.0)).norm() < 1e-12);
//! let back = CoordSystem::Cartesian.convert(cart, CoordSystem::Cylindrical);
//! assert!((back - cyl).norm() < 1e-12);
//! ```

use core::f64::consts::{PI, TAU};

use crate::math::Vec3;

/// A coordinate system.
///
/// Marked `#[non_exhaustive]`: CSXCAD reserves a (currently unimplemented)
/// spherical system, which may be added later.
///
/// ```
/// use csx::CoordSystem;
///
/// assert_eq!(CoordSystem::default(), CoordSystem::Cartesian);
/// assert_eq!(CoordSystem::from_id(1), Some(CoordSystem::Cylindrical));
/// assert_eq!(CoordSystem::Cylindrical.id(), 1);
/// assert_eq!(CoordSystem::from_id(7), None);
/// ```
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum CoordSystem {
    /// `(x, y, z)`; numeric id `0`.
    #[default]
    Cartesian,
    /// `(ρ, α, z)`; numeric id `1`.
    Cylindrical,
}

impl CoordSystem {
    /// Numeric id as used in CSXCAD files.
    pub const fn id(self) -> i32 {
        match self {
            CoordSystem::Cartesian => 0,
            CoordSystem::Cylindrical => 1,
        }
    }

    /// Parses a numeric id as used in CSXCAD files.
    pub const fn from_id(id: i32) -> Option<Self> {
        match id {
            0 => Some(CoordSystem::Cartesian),
            1 => Some(CoordSystem::Cylindrical),
            _ => None,
        }
    }

    /// Converts `p`, expressed in `self`, to Cartesian coordinates.
    #[inline]
    pub fn to_cartesian(self, p: Vec3) -> Vec3 {
        match self {
            CoordSystem::Cartesian => p,
            CoordSystem::Cylindrical => cylindrical_to_cartesian(p),
        }
    }

    /// Converts the Cartesian point `p` into `self`.
    #[inline]
    pub fn from_cartesian(self, p: Vec3) -> Vec3 {
        match self {
            CoordSystem::Cartesian => p,
            CoordSystem::Cylindrical => cartesian_to_cylindrical(p),
        }
    }

    /// Converts `p`, expressed in `self`, into `target` (via Cartesian
    /// coordinates; a no-op when both systems are equal).
    #[inline]
    pub fn convert(self, p: Vec3, target: CoordSystem) -> Vec3 {
        if self == target {
            p
        } else {
            target.from_cartesian(self.to_cartesian(p))
        }
    }
}

/// `(ρ, α, z)` → `(x, y, z)`. Negative `ρ` is allowed and mirrors the
/// point through the z axis.
///
/// ```
/// use csx::coords::cylindrical_to_cartesian;
/// use csx::Vec3;
///
/// let p = cylindrical_to_cartesian(Vec3::new(1.0, core::f64::consts::PI, 3.0));
/// assert!((p - Vec3::new(-1.0, 0.0, 3.0)).norm() < 1e-12);
/// ```
#[inline]
pub fn cylindrical_to_cartesian(p: Vec3) -> Vec3 {
    let (s, c) = p[1].sin_cos();
    Vec3::new(p[0] * c, p[0] * s, p[2])
}

/// `(x, y, z)` → `(ρ, α, z)` with `ρ ≥ 0` and `α ∈ (-π, π]` (`α = 0` on
/// the z axis).
///
/// ```
/// use csx::coords::cartesian_to_cylindrical;
/// use csx::Vec3;
///
/// let c = cartesian_to_cylindrical(Vec3::new(0.0, -2.0, 5.0));
/// assert_eq!(c.x(), 2.0);
/// assert!((c.y() + core::f64::consts::FRAC_PI_2).abs() < 1e-12);
/// ```
#[inline]
pub fn cartesian_to_cylindrical(p: Vec3) -> Vec3 {
    Vec3::new(p[0].hypot(p[1]), p[1].atan2(p[0]), p[2])
}

/// Spherical radius `r`, polar angle `θ` (from +z, in `[0, π]`) and
/// azimuth `α` (in `(-π, π]`) of a Cartesian point. These are the `r`, `t`
/// and `a` variables of weighting functions.
///
/// ```
/// use csx::coords::cartesian_to_spherical;
/// use csx::Vec3;
///
/// let (r, theta, alpha) = cartesian_to_spherical(Vec3::new(0.0, 0.0, -3.0));
/// assert_eq!((r, alpha), (3.0, 0.0));
/// assert!((theta - core::f64::consts::PI).abs() < 1e-12);
/// ```
pub fn cartesian_to_spherical(p: Vec3) -> (f64, f64, f64) {
    let r = p.norm();
    let theta = if r > 0.0 {
        (p[2] / r).clamp(-1.0, 1.0).acos()
    } else {
        0.0
    };
    (r, theta, p[1].atan2(p[0]))
}

/// Wraps an angle into `[-π, π)`.
///
/// ```
/// use core::f64::consts::PI;
/// use csx::coords::wrap_angle;
///
/// assert!((wrap_angle(1.5 * PI) + 0.5 * PI).abs() < 1e-12);
/// assert_eq!(wrap_angle(PI), -PI);
/// ```
#[inline]
pub fn wrap_angle(a: f64) -> f64 {
    let r = (a + PI).rem_euclid(TAU) - PI;
    if r >= PI { r - TAU } else { r }
}

/// Tests whether angle `a` lies within the angular interval going
/// counter-clockwise from `start` to `stop` (inclusive, tolerance `eps`).
/// Intervals spanning `2π` or more contain every angle; `stop < start`
/// describes the same interval as `[stop, start]`.
///
/// ```
/// use core::f64::consts::PI;
/// use csx::coords::angle_in_range;
///
/// // An interval crossing the ±π branch cut.
/// assert!(angle_in_range(-PI + 0.05, 3.0, 3.3, 0.0));
/// assert!(!angle_in_range(0.0, 3.0, 3.3, 0.0));
/// ```
pub fn angle_in_range(a: f64, start: f64, stop: f64, eps: f64) -> bool {
    let (lo, hi) = if stop < start {
        (stop, start)
    } else {
        (start, stop)
    };
    let span = hi - lo;
    if span >= TAU - eps {
        return true;
    }
    let d = (a - lo).rem_euclid(TAU);
    d <= span + eps || d >= TAU - eps
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() {
        let p = Vec3::new(-1.0, 2.0, 3.0);
        let c = cartesian_to_cylindrical(p);
        assert!((cylindrical_to_cartesian(c) - p).norm() < 1e-12);
        assert_eq!(
            CoordSystem::Cylindrical.convert(c, CoordSystem::Cartesian),
            cylindrical_to_cartesian(c)
        );
    }

    #[test]
    fn angles() {
        assert!((wrap_angle(3.0 * PI) + PI).abs() < 1e-12);
        assert!(angle_in_range(0.1, -0.5, 0.5, 0.0));
        assert!(angle_in_range(PI, 3.0, 3.3, 0.0));
        assert!(angle_in_range(-PI + 0.1, 3.0, 3.3, 0.0));
        assert!(!angle_in_range(0.0, 3.0, 3.3, 0.0));
        assert!(angle_in_range(1.0, 0.0, 7.0, 0.0));
    }
}
