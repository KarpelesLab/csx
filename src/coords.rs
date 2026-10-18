//! Coordinate systems and conversions.
//!
//! CSXCAD knows two coordinate systems: Cartesian `(x, y, z)` and cylindrical
//! `(ρ, α, z)` where `α` is the azimuth in radians measured from the +x axis
//! towards +y, and `z` is shared with the Cartesian system.

use core::f64::consts::{PI, TAU};

use crate::math::Vec3;

/// A coordinate system.
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

    /// Converts `p`, expressed in `self`, into `target`.
    #[inline]
    pub fn convert(self, p: Vec3, target: CoordSystem) -> Vec3 {
        if self == target {
            p
        } else {
            target.from_cartesian(self.to_cartesian(p))
        }
    }
}

/// `(ρ, α, z)` → `(x, y, z)`.
#[inline]
pub fn cylindrical_to_cartesian(p: Vec3) -> Vec3 {
    let (s, c) = p[1].sin_cos();
    Vec3::new(p[0] * c, p[0] * s, p[2])
}

/// `(x, y, z)` → `(ρ, α, z)` with `α ∈ (-π, π]`.
#[inline]
pub fn cartesian_to_cylindrical(p: Vec3) -> Vec3 {
    Vec3::new(p[0].hypot(p[1]), p[1].atan2(p[0]), p[2])
}

/// Spherical radius `r`, polar angle `θ` (from +z) and azimuth `α` of a
/// Cartesian point.
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
#[inline]
pub fn wrap_angle(a: f64) -> f64 {
    let r = (a + PI).rem_euclid(TAU) - PI;
    if r >= PI { r - TAU } else { r }
}

/// Tests whether angle `a` lies within the angular interval going
/// counter-clockwise from `start` to `stop` (inclusive, tolerance `eps`).
/// Intervals spanning `2π` or more contain every angle; `stop < start`
/// describes the same interval as `[stop, start]`.
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
