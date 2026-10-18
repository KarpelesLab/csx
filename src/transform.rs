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
    /// Arbitrary 4×4 matrix: 16 values in row-major order acting on column
    /// vectors, `p' = M · [x y z 1]ᵀ` (the translation is in the last
    /// column). This matches the behaviour of CSXCAD.
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
            TransformOp::Scale(_) => "Scale",
            TransformOp::Scale3(_) => "Scale3",
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
            ("Scale" | "Scale3", 3) => TransformOp::Scale3(next3(&mut it)),
            ("Scale3", 1) => TransformOp::Scale(it.next().unwrap()),
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

    /// The operation as a matrix acting on column vectors (`p' = M·p`),
    /// with rotation angles in radians (`degrees = false`) or degrees.
    pub fn to_mat4(&self, degrees: bool) -> Mat4 {
        self.matrix(degrees)
    }
}

/// An ordered list of transformations with cached forward and inverse
/// matrices.
#[derive(Debug, Clone, PartialEq)]
pub struct Transform {
    ops: Vec<TransformOp>,
    degrees: bool,
    matrix: Mat4,
    inverse: Mat4,
    invertible: bool,
}

impl Default for Transform {
    fn default() -> Self {
        Transform {
            ops: Vec::new(),
            degrees: false,
            matrix: Mat4::IDENTITY,
            inverse: Mat4::IDENTITY,
            invertible: true,
        }
    }
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

    /// Appends an arbitrary matrix acting on column vectors (`p' = m·p`).
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
        if !self.invertible {
            return Err(Error::Invalid("singular transformation".into()));
        }
        Ok(())
    }

    fn rebuild(&mut self) {
        self.matrix = self
            .ops
            .iter()
            .fold(Mat4::IDENTITY, |acc, op| op.matrix(self.degrees).mul(&acc));
        match self.matrix.inverse() {
            Some(inv) => {
                self.inverse = inv;
                self.invertible = true;
            }
            None => {
                self.inverse = Mat4::IDENTITY;
                self.invertible = false;
            }
        }
    }

    /// False if the combined matrix is singular (e.g. a zero scale); such a
    /// transform maps every shape to something without interior, and
    /// inside tests of primitives using it return `false`.
    pub fn is_invertible(&self) -> bool {
        self.invertible
    }

    /// Forward matrix (local → world).
    pub fn matrix(&self) -> &Mat4 {
        &self.matrix
    }

    /// Inverse matrix (world → local); the identity when the transform is
    /// not invertible (see [`Transform::is_invertible`]).
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

    #[test]
    fn matrix_convention() {
        // Verified against CSXCAD: translation in the last column.
        let col = [
            1.0, 0.0, 0.0, 5.0, 0.0, 1.0, 0.0, 6.0, 0.0, 0.0, 1.0, 7.0, 0.0, 0.0, 0.0, 1.0,
        ];
        let op = TransformOp::from_name_args("Matrix", col.map(Scalar::new).to_vec()).unwrap();
        assert_eq!(
            op.to_mat4(false).transform_point(Vec3::ZERO),
            Vec3::new(5.0, 6.0, 7.0)
        );
        // [[0,-1],[1,0]] rotates +x onto +y.
        let rot = [
            0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        ];
        let op = TransformOp::from_name_args("Matrix", rot.map(Scalar::new).to_vec()).unwrap();
        assert!((op.to_mat4(false).transform_point(Vec3::unit(0)) - Vec3::unit(1)).norm() < 1e-12);
        let r = Mat4::rotation(Vec3::unit(2), core::f64::consts::FRAC_PI_2);
        let mut t = Transform::new();
        t.matrix_op(&r);
        assert!((t.apply(Vec3::unit(0)) - Vec3::unit(1)).norm() < 1e-12);
        let mut t2 = Transform::new();
        t2.matrix_op(&Mat4::translation(Vec3::new(1.0, 2.0, 3.0)));
        assert_eq!(t2.apply(Vec3::ZERO), Vec3::new(1.0, 2.0, 3.0));
        // Vector scaling is written as `Scale3`; both names are accepted.
        let s3 = TransformOp::from_name_args("Scale3", vec![1.0.into(), 2.0.into(), 3.0.into()])
            .unwrap();
        assert_eq!(s3.name(), "Scale3");
        assert!(
            TransformOp::from_name_args("Scale", vec![1.0.into(), 2.0.into(), 3.0.into()]).is_ok()
        );
    }
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
        let mut z = Transform::new();
        z.scale(0.0);
        assert!(!z.is_invertible());
        assert!(Transform::new().is_invertible());
        let op = TransformOp::from_name_args(
            "Rotate_Origin",
            vec![0.0.into(), 0.0.into(), 1.0.into(), 2.0.into()],
        )
        .unwrap();
        assert_eq!(op.args().len(), 4);
    }
}
