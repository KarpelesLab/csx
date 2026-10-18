//! Affine transformations of primitives.
//!
//! A [`Transform`] is an ordered list of elementary operations. Operations
//! are applied to the primitive's geometry in list order (the first entry is
//! applied first). Angles are in radians unless the operation was created
//! with [`Transform::set_angles_in_degrees`] semantics (see below).

use crate::error::{Error, Result};
use crate::math::{Mat4, Vec3};
use crate::param::{ParameterSet, Scalar};

/// One elementary transformation.
#[derive(Debug, Clone, PartialEq)]
pub enum TransformOp {
    /// Translation by a vector.
    Translate([Scalar; 3]),
    /// Rotation about an axis through the origin: axis vector then angle.
    RotateOrigin([Scalar; 3], Scalar),
    /// Rotation about the x axis.
    RotateX(Scalar),
    /// Rotation about the y axis.
    RotateY(Scalar),
    /// Rotation about the z axis.
    RotateZ(Scalar),
    /// Uniform scaling.
    Scale(Scalar),
    /// Non-uniform scaling.
    Scale3([Scalar; 3]),
    /// Arbitrary 4×4 matrix (row major).
    Matrix(Box<[Scalar; 16]>),
}

impl TransformOp {
    /// XML element name of the operation.
    pub fn name(&self) -> &'static str {
        match self {
            TransformOp::Translate(_) => "Translate",
            TransformOp::RotateOrigin(..) => "Rotate_Origin",
            TransformOp::RotateX(_) => "Rotate_X",
            TransformOp::RotateY(_) => "Rotate_Y",
            TransformOp::RotateZ(_) => "Rotate_Z",
            TransformOp::Scale(_) | TransformOp::Scale3(_) => "Scale",
            TransformOp::Matrix(_) => "Matrix",
        }
    }

    /// Arguments in their textual order.
    pub fn args(&self) -> Vec<&Scalar> {
        match self {
            TransformOp::Translate(v) | TransformOp::Scale3(v) => v.iter().collect(),
            TransformOp::RotateOrigin(v, a) => v.iter().chain(core::iter::once(a)).collect(),
            TransformOp::RotateX(a)
            | TransformOp::RotateY(a)
            | TransformOp::RotateZ(a)
            | TransformOp::Scale(a) => vec![a],
            TransformOp::Matrix(m) => m.iter().collect(),
        }
    }

    /// Builds an operation from its element name and arguments.
    pub fn from_name_args(name: &str, args: Vec<Scalar>) -> Result<TransformOp> {
        let n = args.len();
        let bad = || Error::Format(format!("transformation {name} with {n} argument(s)"));
        let mut it = args.into_iter();
        let next3 = |it: &mut std::vec::IntoIter<Scalar>| -> [Scalar; 3] {
            [it.next().unwrap(), it.next().unwrap(), it.next().unwrap()]
        };
        Ok(match (name, n) {
            ("Translate", 3) => TransformOp::Translate(next3(&mut it)),
            ("Rotate_Origin", 4) => {
                let v = next3(&mut it);
                TransformOp::RotateOrigin(v, it.next().unwrap())
            }
            ("Rotate_X", 1) => TransformOp::RotateX(it.next().unwrap()),
            ("Rotate_Y", 1) => TransformOp::RotateY(it.next().unwrap()),
            ("Rotate_Z", 1) => TransformOp::RotateZ(it.next().unwrap()),
            ("Scale", 1) => TransformOp::Scale(it.next().unwrap()),
            ("Scale", 3) => TransformOp::Scale3(next3(&mut it)),
            ("Matrix", 16) => {
                let v: Vec<Scalar> = it.collect();
                TransformOp::Matrix(Box::new(v.try_into().map_err(|_| bad())?))
            }
            _ => return Err(bad()),
        })
    }

    fn update(&mut self, p: &ParameterSet) -> Result<()> {
        match self {
            TransformOp::Translate(v) | TransformOp::Scale3(v) => {
                v.iter_mut().try_for_each(|s| s.update(p))
            }
            TransformOp::RotateOrigin(v, a) => {
                v.iter_mut().try_for_each(|s| s.update(p))?;
                a.update(p)
            }
            TransformOp::RotateX(a)
            | TransformOp::RotateY(a)
            | TransformOp::RotateZ(a)
            | TransformOp::Scale(a) => a.update(p),
            TransformOp::Matrix(m) => m.iter_mut().try_for_each(|s| s.update(p)),
        }
    }

    fn matrix(&self, deg: bool) -> Mat4 {
        let ang = |a: &Scalar| {
            if deg {
                a.value().to_radians()
            } else {
                a.value()
            }
        };
        let v3 = |v: &[Scalar; 3]| Vec3::new(v[0].value(), v[1].value(), v[2].value());
        match self {
            TransformOp::Translate(v) => Mat4::translation(v3(v)),
            TransformOp::RotateOrigin(v, a) => Mat4::rotation(v3(v), ang(a)),
            TransformOp::RotateX(a) => Mat4::rotation(Vec3::unit(0), ang(a)),
            TransformOp::RotateY(a) => Mat4::rotation(Vec3::unit(1), ang(a)),
            TransformOp::RotateZ(a) => Mat4::rotation(Vec3::unit(2), ang(a)),
            TransformOp::Scale(s) => {
                let s = s.value();
                Mat4::scaling(Vec3::new(s, s, s))
            }
            TransformOp::Scale3(v) => Mat4::scaling(v3(v)),
            TransformOp::Matrix(m) => Mat4::from_row_major(&core::array::from_fn(|i| m[i].value())),
        }
    }
}

/// An ordered list of transformations with cached forward and inverse
/// matrices.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Transform {
    ops: Vec<TransformOp>,
    degrees: bool,
    matrix: Mat4,
    inverse: Mat4,
}

impl Transform {
    /// Identity transform.
    pub fn new() -> Self {
        Self::default()
    }

    /// Whether rotation angles are interpreted in degrees (default: radians).
    pub fn angles_in_degrees(&self) -> bool {
        self.degrees
    }

    /// Selects degree (`true`) or radian (`false`) rotation angles.
    pub fn set_angles_in_degrees(&mut self, deg: bool) {
        self.degrees = deg;
        self.rebuild();
    }

    /// The operations in application order.
    pub fn ops(&self) -> &[TransformOp] {
        &self.ops
    }

    /// True if no operation is present.
    pub fn is_empty(&self) -> bool {
        self.ops.is_empty()
    }

    /// Appends an operation (applied after the existing ones).
    pub fn push(&mut self, op: TransformOp) -> &mut Self {
        self.ops.push(op);
        self.rebuild();
        self
    }

    /// Appends a translation.
    pub fn translate(&mut self, v: Vec3) -> &mut Self {
        self.push(TransformOp::Translate(v.0.map(Scalar::new)))
    }

    /// Appends a rotation about an arbitrary axis through the origin.
    pub fn rotate(&mut self, axis: Vec3, angle: f64) -> &mut Self {
        self.push(TransformOp::RotateOrigin(
            axis.0.map(Scalar::new),
            angle.into(),
        ))
    }

    /// Appends a rotation about the coordinate axis `dir` (0, 1, 2).
    pub fn rotate_axis(&mut self, dir: usize, angle: f64) -> &mut Self {
        let a = Scalar::new(angle);
        self.push(match dir {
            0 => TransformOp::RotateX(a),
            1 => TransformOp::RotateY(a),
            _ => TransformOp::RotateZ(a),
        })
    }

    /// Appends a uniform scaling.
    pub fn scale(&mut self, s: f64) -> &mut Self {
        self.push(TransformOp::Scale(s.into()))
    }

    /// Appends an arbitrary matrix.
    pub fn matrix_op(&mut self, m: &Mat4) -> &mut Self {
        self.push(TransformOp::Matrix(Box::new(
            m.to_row_major().map(Scalar::new),
        )))
    }

    /// Removes all operations.
    pub fn clear(&mut self) {
        self.ops.clear();
        self.rebuild();
    }

    /// Re-evaluates parametric arguments and rebuilds the matrices.
    pub fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.ops.iter_mut().try_for_each(|o| o.update(p))?;
        self.rebuild();
        if self.matrix.inverse().is_none() {
            return Err(Error::Invalid("singular transformation".into()));
        }
        Ok(())
    }

    fn rebuild(&mut self) {
        self.matrix = self
            .ops
            .iter()
            .fold(Mat4::IDENTITY, |acc, op| op.matrix(self.degrees).mul(&acc));
        self.inverse = self.matrix.inverse().unwrap_or(Mat4::IDENTITY);
    }

    /// Forward matrix (local → world).
    pub fn matrix(&self) -> &Mat4 {
        &self.matrix
    }

    /// Inverse matrix (world → local).
    pub fn inverse(&self) -> &Mat4 {
        &self.inverse
    }

    /// Maps a local point to world coordinates.
    #[inline]
    pub fn apply(&self, p: Vec3) -> Vec3 {
        self.matrix.transform_point(p)
    }

    /// Maps a world point back to local coordinates.
    #[inline]
    pub fn apply_inverse(&self, p: Vec3) -> Vec3 {
        self.inverse.transform_point(p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::FRAC_PI_2;

    #[test]
    fn order() {
        let mut t = Transform::new();
        t.translate(Vec3::new(1.0, 0.0, 0.0))
            .rotate_axis(2, FRAC_PI_2);
        // translate first, then rotate: (0,0,0) -> (1,0,0) -> (0,1,0)
        let p = t.apply(Vec3::ZERO);
        assert!((p - Vec3::new(0.0, 1.0, 0.0)).norm() < 1e-12);
        assert!((t.apply_inverse(p)).norm() < 1e-12);
    }

    #[test]
    fn degrees_and_params() {
        let mut t = Transform::new();
        t.set_angles_in_degrees(true);
        t.push(TransformOp::RotateZ(Scalar::expr("ang")));
        let mut ps = ParameterSet::new();
        ps.set("ang", 90.0);
        t.update(&ps).unwrap();
        assert!((t.apply(Vec3::unit(0)) - Vec3::unit(1)).norm() < 1e-12);
        assert!(TransformOp::from_name_args("Translate", vec![1.0.into()]).is_err());
        let op = TransformOp::from_name_args(
            "Rotate_Origin",
            vec![0.0.into(), 0.0.into(), 1.0.into(), 2.0.into()],
        )
        .unwrap();
        assert_eq!(op.args().len(), 4);
    }
}
