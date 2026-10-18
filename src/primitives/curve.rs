//! One- and zero-dimensional primitives: curves, wires and points.

use crate::error::{Error, Result};
use crate::math::{BoundingBox, Vec3};
use crate::param::{ParameterSet, Scalar, ScalarPoint};

use super::geom::{dist_point_segment, tol};

fn polyline_dist(points: &[ScalarPoint], p: Vec3) -> f64 {
    let pts: Vec<Vec3> = points.iter().map(|q| q.value()).collect();
    match pts.len() {
        0 => f64::INFINITY,
        1 => (p - pts[0]).norm(),
        _ => pts
            .windows(2)
            .map(|w| dist_point_segment(p, w[0], w[1]))
            .fold(f64::INFINITY, f64::min),
    }
}

fn polyline_bbox(points: &[ScalarPoint]) -> BoundingBox {
    BoundingBox::from_points(points.iter().map(|q| q.value()))
}

fn polyline_scale(points: &[ScalarPoint]) -> f64 {
    points
        .iter()
        .flat_map(|q| q.value().0)
        .fold(0.0, |m, c| m.max(c.abs()))
}

/// A polyline through a list of points (1-D).
///
/// A point is "inside" a curve when it lies on it (within a tiny relative
/// tolerance).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Curve {
    /// Polyline vertices.
    pub points: Vec<ScalarPoint>,
}

impl Curve {
    /// Curve through `points`.
    pub fn new<P: Into<ScalarPoint>>(points: impl IntoIterator<Item = P>) -> Self {
        Curve {
            points: points.into_iter().map(Into::into).collect(),
        }
    }

    pub(crate) fn inside(&self, p: Vec3) -> bool {
        polyline_dist(&self.points, p) <= tol(polyline_scale(&self.points), 0.0)
    }

    pub(crate) fn bbox(&self) -> BoundingBox {
        polyline_bbox(&self.points)
    }

    pub(crate) fn dimension(&self) -> usize {
        if self.points.len() >= 2 { 1 } else { 0 }
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.points.iter_mut().try_for_each(|q| q.update(p))
    }
}

/// A curve with a circular cross section of radius `radius`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Wire {
    /// Polyline vertices.
    pub points: Vec<ScalarPoint>,
    /// Wire radius.
    pub radius: Scalar,
}

impl Wire {
    /// Wire through `points` with radius `radius`.
    pub fn new<P: Into<ScalarPoint>>(
        points: impl IntoIterator<Item = P>,
        radius: impl Into<Scalar>,
    ) -> Self {
        Wire {
            points: points.into_iter().map(Into::into).collect(),
            radius: radius.into(),
        }
    }

    pub(crate) fn inside(&self, p: Vec3) -> bool {
        let r = self.radius.value();
        polyline_dist(&self.points, p) <= r + tol(polyline_scale(&self.points), r)
    }

    pub(crate) fn bbox(&self) -> BoundingBox {
        let b = polyline_bbox(&self.points);
        if b.is_empty() {
            b
        } else {
            b.expanded(self.radius.value())
        }
    }

    pub(crate) fn dimension(&self) -> usize {
        if self.radius.value() > 0.0 {
            3
        } else if self.points.len() >= 2 {
            1
        } else {
            0
        }
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.points.iter_mut().try_for_each(|q| q.update(p))?;
        self.radius.update(p)?;
        if self.radius.value() < 0.0 {
            return Err(Error::Invalid(format!(
                "negative wire radius {}",
                self.radius.value()
            )));
        }
        Ok(())
    }
}

/// A single point (0-D).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PointShape {
    /// Location.
    pub point: ScalarPoint,
}

impl PointShape {
    /// Point primitive at `p`.
    pub fn new(p: impl Into<ScalarPoint>) -> Self {
        PointShape { point: p.into() }
    }

    pub(crate) fn inside(&self, p: Vec3) -> bool {
        let q = self.point.value();
        (p - q).norm() <= tol(q.norm(), 0.0)
    }

    pub(crate) fn bbox(&self) -> BoundingBox {
        let q = self.point.value();
        BoundingBox::from_corners(q, q)
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.point.update(p)
    }
}
