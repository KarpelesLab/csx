//! Geometric helper routines shared by several primitives.

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use crate::coords::angle_in_range;
use crate::math::{BoundingBox, Vec3};

/// Relative tolerance used for closed-interval comparisons.
pub(crate) const REL_TOL: f64 = 1e-12;

/// Absolute tolerance scaled from the magnitudes involved.
#[inline]
pub(crate) fn tol(a: f64, b: f64) -> f64 {
    REL_TOL * a.abs().max(b.abs()).max(1e-300)
}

/// `lo <= v <= hi` (bounds in any order) with a tiny relative tolerance.
#[inline]
pub(crate) fn in_range(v: f64, a: f64, b: f64) -> bool {
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    let t = tol(lo, hi).max(tol(v, 0.0));
    v >= lo - t && v <= hi + t
}

/// Bounding box (Cartesian) of the cylindrical sector
/// `ρ ∈ [r1, r2], α ∈ [a1, a2], z ∈ [z1, z2]`.
pub(crate) fn sector_bbox(r: (f64, f64), a: (f64, f64), z: (f64, f64)) -> BoundingBox {
    let (r1, r2) = (r.0.min(r.1), r.0.max(r.1));
    let (a1, a2) = (a.0.min(a.1), a.0.max(a.1));
    let mut b = BoundingBox::EMPTY;
    let full = a2 - a1 >= TAU;
    let mut add = |ang: f64| {
        let (s, c) = ang.sin_cos();
        for rr in [r1, r2] {
            for zz in [z.0, z.1] {
                b.include(Vec3::new(rr * c, rr * s, zz));
            }
        }
    };
    if full {
        for k in 0..4 {
            add(k as f64 * FRAC_PI_2);
        }
    } else {
        add(a1);
        add(a2);
        // Cardinal directions inside the angular interval.
        let first = (a1 / FRAC_PI_2).ceil() as i64;
        let last = (a2 / FRAC_PI_2).floor() as i64;
        for k in first..=last.min(first.saturating_add(4)) {
            add(k as f64 * FRAC_PI_2);
        }
    }
    b
}

/// Conservative cylindrical bounds `(ρ, α, z)` of a Cartesian box.
pub(crate) fn cartesian_bbox_to_cylindrical(b: &BoundingBox) -> BoundingBox {
    if b.is_empty() {
        return *b;
    }
    let (x0, x1, y0, y1) = (b.min[0], b.max[0], b.min[1], b.max[1]);
    let cx = 0.0f64.clamp(x0, x1);
    let cy = 0.0f64.clamp(y0, y1);
    let rmin = cx.hypot(cy);
    let corners = [(x0, y0), (x0, y1), (x1, y0), (x1, y1)];
    let rmax = corners.iter().map(|&(x, y)| x.hypot(y)).fold(0.0, f64::max);
    let (amin, amax) = if x0 <= 0.0 && x1 >= 0.0 && y0 <= 0.0 && y1 >= 0.0 {
        (-PI, PI)
    } else if x1 < 0.0 && y0 < 0.0 && y1 >= 0.0 {
        // Crosses the ±π branch cut: measure angles in [0, 2π).
        let angs = corners.map(|(x, y)| y.atan2(x).rem_euclid(TAU));
        let lo = angs.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = angs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        (lo, hi)
    } else {
        let angs = corners.map(|(x, y)| y.atan2(x));
        let lo = angs.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = angs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        (lo, hi)
    };
    BoundingBox {
        min: Vec3::new(rmin, amin, b.min[2]),
        max: Vec3::new(rmax, amax, b.max[2]),
    }
}

/// Point-in-polygon test (even–odd rule) with boundary points counted as
/// inside.
pub(crate) fn point_in_polygon(p: (f64, f64), poly: &[(f64, f64)]) -> bool {
    let n = poly.len();
    if n == 0 {
        return false;
    }
    let scale = poly
        .iter()
        .fold(0.0f64, |m, &(x, y)| m.max(x.abs()).max(y.abs()));
    let eps = REL_TOL * scale.max(1e-300);
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = poly[i];
        let (xj, yj) = poly[j];
        if dist_point_segment_2d(p, (xi, yi), (xj, yj)) <= eps {
            return true;
        }
        if (yi > p.1) != (yj > p.1) {
            let x = xi + (p.1 - yi) * (xj - xi) / (yj - yi);
            if p.0 < x {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

fn dist_point_segment_2d(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (dx, dy) = (b.0 - a.0, b.1 - a.1);
    let l2 = dx * dx + dy * dy;
    let t = if l2 > 0.0 {
        (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (p.0 - a.0 - t * dx).hypot(p.1 - a.1 - t * dy)
}

/// Distance from `p` to the segment `a`–`b`.
pub(crate) fn dist_point_segment(p: Vec3, a: Vec3, b: Vec3) -> f64 {
    let d = b - a;
    let l2 = d.norm_sq();
    let t = if l2 > 0.0 {
        ((p - a).dot(d) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (p - (a + d * t)).norm()
}

/// Angle test helper re-exported for primitives.
#[inline]
pub(crate) fn angle_ok(a: f64, start: f64, stop: f64) -> bool {
    angle_in_range(a, start, stop, 1e-12)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polygon() {
        let sq = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        assert!(point_in_polygon((0.5, 0.5), &sq));
        assert!(point_in_polygon((1.0, 0.5), &sq));
        assert!(point_in_polygon((0.0, 0.0), &sq));
        assert!(!point_in_polygon((1.5, 0.5), &sq));
        let l = [
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (1.0, 1.0),
            (1.0, 2.0),
            (0.0, 2.0),
        ];
        assert!(!point_in_polygon((1.5, 1.5), &l));
        assert!(point_in_polygon((0.5, 1.5), &l));
    }

    #[test]
    fn sectors() {
        let b = sector_bbox((1.0, 2.0), (0.0, FRAC_PI_2), (0.0, 1.0));
        assert!((b.min - Vec3::new(0.0, 0.0, 0.0)).norm() < 1e-12);
        assert!((b.max - Vec3::new(2.0, 2.0, 1.0)).norm() < 1e-12);
        let c = cartesian_bbox_to_cylindrical(&BoundingBox::from_corners(
            Vec3::new(-2.0, -1.0, 0.0),
            Vec3::new(-1.0, 1.0, 1.0),
        ));
        assert!(c.min[1] > FRAC_PI_2 && c.max[1] < 3.0 * FRAC_PI_2);
        assert!((c.min[0] - 1.0).abs() < 1e-12);
    }
}
