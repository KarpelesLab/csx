//! Basic linear algebra: 3-vectors, affine 4×4 matrices and bounding boxes.

use core::ops::{Add, AddAssign, Div, Index, IndexMut, Mul, Neg, Sub, SubAssign};

/// A 3-component vector of `f64`, used for points and directions.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3(pub [f64; 3]);

impl Vec3 {
    /// The zero vector.
    pub const ZERO: Vec3 = Vec3([0.0; 3]);

    /// Creates a vector from its components.
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3([x, y, z])
    }

    /// Unit vector along axis `dir` (0 = x, 1 = y, 2 = z).
    #[inline]
    pub fn unit(dir: usize) -> Self {
        let mut v = Vec3::ZERO;
        v.0[dir] = 1.0;
        v
    }

    /// First component.
    #[inline]
    pub const fn x(&self) -> f64 {
        self.0[0]
    }
    /// Second component.
    #[inline]
    pub const fn y(&self) -> f64 {
        self.0[1]
    }
    /// Third component.
    #[inline]
    pub const fn z(&self) -> f64 {
        self.0[2]
    }

    /// Dot product.
    #[inline]
    pub fn dot(self, o: Vec3) -> f64 {
        self.0[0] * o.0[0] + self.0[1] * o.0[1] + self.0[2] * o.0[2]
    }

    /// Cross product.
    #[inline]
    pub fn cross(self, o: Vec3) -> Vec3 {
        Vec3::new(
            self.0[1] * o.0[2] - self.0[2] * o.0[1],
            self.0[2] * o.0[0] - self.0[0] * o.0[2],
            self.0[0] * o.0[1] - self.0[1] * o.0[0],
        )
    }

    /// Squared Euclidean norm.
    #[inline]
    pub fn norm_sq(self) -> f64 {
        self.dot(self)
    }

    /// Euclidean norm.
    #[inline]
    pub fn norm(self) -> f64 {
        self.norm_sq().sqrt()
    }

    /// Returns the vector scaled to unit length (or zero if it has no length).
    pub fn normalized(self) -> Vec3 {
        let n = self.norm();
        if n > 0.0 { self / n } else { Vec3::ZERO }
    }

    /// Component-wise minimum.
    #[inline]
    pub fn min(self, o: Vec3) -> Vec3 {
        Vec3::new(
            self.0[0].min(o.0[0]),
            self.0[1].min(o.0[1]),
            self.0[2].min(o.0[2]),
        )
    }

    /// Component-wise maximum.
    #[inline]
    pub fn max(self, o: Vec3) -> Vec3 {
        Vec3::new(
            self.0[0].max(o.0[0]),
            self.0[1].max(o.0[1]),
            self.0[2].max(o.0[2]),
        )
    }

    /// True if all components are finite.
    pub fn is_finite(self) -> bool {
        self.0.iter().all(|c| c.is_finite())
    }
}

impl From<[f64; 3]> for Vec3 {
    fn from(a: [f64; 3]) -> Self {
        Vec3(a)
    }
}

impl From<Vec3> for [f64; 3] {
    fn from(v: Vec3) -> Self {
        v.0
    }
}

impl Index<usize> for Vec3 {
    type Output = f64;
    #[inline]
    fn index(&self, i: usize) -> &f64 {
        &self.0[i]
    }
}

impl IndexMut<usize> for Vec3 {
    #[inline]
    fn index_mut(&mut self, i: usize) -> &mut f64 {
        &mut self.0[i]
    }
}

impl Add for Vec3 {
    type Output = Vec3;
    #[inline]
    fn add(self, o: Vec3) -> Vec3 {
        Vec3::new(self.0[0] + o.0[0], self.0[1] + o.0[1], self.0[2] + o.0[2])
    }
}

impl AddAssign for Vec3 {
    #[inline]
    fn add_assign(&mut self, o: Vec3) {
        *self = *self + o;
    }
}

impl Sub for Vec3 {
    type Output = Vec3;
    #[inline]
    fn sub(self, o: Vec3) -> Vec3 {
        Vec3::new(self.0[0] - o.0[0], self.0[1] - o.0[1], self.0[2] - o.0[2])
    }
}

impl SubAssign for Vec3 {
    #[inline]
    fn sub_assign(&mut self, o: Vec3) {
        *self = *self - o;
    }
}

impl Mul<f64> for Vec3 {
    type Output = Vec3;
    #[inline]
    fn mul(self, s: f64) -> Vec3 {
        Vec3::new(self.0[0] * s, self.0[1] * s, self.0[2] * s)
    }
}

impl Div<f64> for Vec3 {
    type Output = Vec3;
    #[inline]
    fn div(self, s: f64) -> Vec3 {
        Vec3::new(self.0[0] / s, self.0[1] / s, self.0[2] / s)
    }
}

impl Neg for Vec3 {
    type Output = Vec3;
    #[inline]
    fn neg(self) -> Vec3 {
        Vec3::new(-self.0[0], -self.0[1], -self.0[2])
    }
}

/// A 4×4 matrix in row-major order, used for affine transformations of
/// homogeneous coordinates (`p' = M · [p, 1]`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mat4(pub [[f64; 4]; 4]);

impl Default for Mat4 {
    fn default() -> Self {
        Mat4::IDENTITY
    }
}

impl Mat4 {
    /// The identity matrix.
    pub const IDENTITY: Mat4 = Mat4([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]);

    /// Builds a matrix from 16 values in row-major order.
    pub fn from_row_major(v: &[f64; 16]) -> Mat4 {
        let mut m = [[0.0; 4]; 4];
        for (i, row) in m.iter_mut().enumerate() {
            row.copy_from_slice(&v[i * 4..i * 4 + 4]);
        }
        Mat4(m)
    }

    /// Returns the 16 values in row-major order.
    pub fn to_row_major(&self) -> [f64; 16] {
        let mut out = [0.0; 16];
        for i in 0..4 {
            out[i * 4..i * 4 + 4].copy_from_slice(&self.0[i]);
        }
        out
    }

    /// Translation by `t`.
    pub fn translation(t: Vec3) -> Mat4 {
        let mut m = Mat4::IDENTITY;
        m.0[0][3] = t[0];
        m.0[1][3] = t[1];
        m.0[2][3] = t[2];
        m
    }

    /// Non-uniform scaling.
    pub fn scaling(s: Vec3) -> Mat4 {
        let mut m = Mat4::IDENTITY;
        m.0[0][0] = s[0];
        m.0[1][1] = s[1];
        m.0[2][2] = s[2];
        m
    }

    /// Right-handed rotation by `angle` radians about `axis` through the origin.
    pub fn rotation(axis: Vec3, angle: f64) -> Mat4 {
        let a = axis.normalized();
        let (s, c) = angle.sin_cos();
        let t = 1.0 - c;
        let (x, y, z) = (a[0], a[1], a[2]);
        Mat4([
            [t * x * x + c, t * x * y - s * z, t * x * z + s * y, 0.0],
            [t * x * y + s * z, t * y * y + c, t * y * z - s * x, 0.0],
            [t * x * z - s * y, t * y * z + s * x, t * z * z + c, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }

    /// Matrix product `self · o` (applies `o` first, then `self`).
    pub fn mul(&self, o: &Mat4) -> Mat4 {
        let mut r = [[0.0; 4]; 4];
        for (i, row) in r.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = (0..4).map(|k| self.0[i][k] * o.0[k][j]).sum();
            }
        }
        Mat4(r)
    }

    /// Transposed matrix.
    pub fn transpose(&self) -> Mat4 {
        Mat4(core::array::from_fn(|i| {
            core::array::from_fn(|j| self.0[j][i])
        }))
    }

    /// Transforms a point (w = 1).
    #[inline]
    pub fn transform_point(&self, p: Vec3) -> Vec3 {
        let m = &self.0;
        let x = m[0][0] * p[0] + m[0][1] * p[1] + m[0][2] * p[2] + m[0][3];
        let y = m[1][0] * p[0] + m[1][1] * p[1] + m[1][2] * p[2] + m[1][3];
        let z = m[2][0] * p[0] + m[2][1] * p[1] + m[2][2] * p[2] + m[2][3];
        let w = m[3][0] * p[0] + m[3][1] * p[1] + m[3][2] * p[2] + m[3][3];
        if w != 1.0 && w != 0.0 {
            Vec3::new(x / w, y / w, z / w)
        } else {
            Vec3::new(x, y, z)
        }
    }

    /// Transforms a direction (w = 0, translation ignored).
    #[inline]
    pub fn transform_vector(&self, v: Vec3) -> Vec3 {
        let m = &self.0;
        Vec3::new(
            m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
        )
    }

    /// General 4×4 inverse via Gauss–Jordan elimination with partial
    /// pivoting. Returns `None` for (numerically) singular matrices.
    pub fn inverse(&self) -> Option<Mat4> {
        let mut a = self.0;
        let mut inv = Mat4::IDENTITY.0;
        for col in 0..4 {
            let pivot = (col..4).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
            if a[pivot][col].abs() < 1e-300 {
                return None;
            }
            a.swap(col, pivot);
            inv.swap(col, pivot);
            let d = a[col][col];
            for k in 0..4 {
                a[col][k] /= d;
                inv[col][k] /= d;
            }
            for row in 0..4 {
                if row != col {
                    let f = a[row][col];
                    if f != 0.0 {
                        for k in 0..4 {
                            a[row][k] -= f * a[col][k];
                            inv[row][k] -= f * inv[col][k];
                        }
                    }
                }
            }
        }
        Some(Mat4(inv))
    }

    /// True if this is (exactly) the identity matrix.
    pub fn is_identity(&self) -> bool {
        *self == Mat4::IDENTITY
    }
}

/// An axis-aligned bounding box. An *empty* box has `min > max`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundingBox {
    /// Lower corner.
    pub min: Vec3,
    /// Upper corner.
    pub max: Vec3,
}

impl Default for BoundingBox {
    fn default() -> Self {
        BoundingBox::EMPTY
    }
}

impl BoundingBox {
    /// The empty box (neutral element for [`BoundingBox::union`]).
    pub const EMPTY: BoundingBox = BoundingBox {
        min: Vec3([f64::INFINITY; 3]),
        max: Vec3([f64::NEG_INFINITY; 3]),
    };

    /// Box spanned by two arbitrary corners.
    pub fn from_corners(a: Vec3, b: Vec3) -> Self {
        BoundingBox {
            min: a.min(b),
            max: a.max(b),
        }
    }

    /// Smallest box containing all `points`.
    pub fn from_points<I: IntoIterator<Item = Vec3>>(points: I) -> Self {
        let mut b = BoundingBox::EMPTY;
        for p in points {
            b.include(p);
        }
        b
    }

    /// True if the box contains no point.
    pub fn is_empty(&self) -> bool {
        (0..3).any(|i| {
            self.min[i]
                .partial_cmp(&self.max[i])
                .is_none_or(|o| o.is_gt())
        })
    }

    /// Grows the box to contain `p`.
    pub fn include(&mut self, p: Vec3) {
        self.min = self.min.min(p);
        self.max = self.max.max(p);
    }

    /// Union of two boxes.
    pub fn union(&self, o: &BoundingBox) -> BoundingBox {
        if self.is_empty() {
            return *o;
        }
        if o.is_empty() {
            return *self;
        }
        BoundingBox {
            min: self.min.min(o.min),
            max: self.max.max(o.max),
        }
    }

    /// Intersection of two boxes (possibly empty).
    pub fn intersection(&self, o: &BoundingBox) -> BoundingBox {
        BoundingBox {
            min: self.min.max(o.min),
            max: self.max.min(o.max),
        }
    }

    /// Box enlarged by `d` on every side.
    pub fn expanded(&self, d: f64) -> BoundingBox {
        let e = Vec3::new(d, d, d);
        BoundingBox {
            min: self.min - e,
            max: self.max + e,
        }
    }

    /// Point containment (closed box).
    #[inline]
    pub fn contains(&self, p: Vec3) -> bool {
        (0..3).all(|i| p[i] >= self.min[i] && p[i] <= self.max[i])
    }

    /// True if the two (closed) boxes overlap.
    pub fn intersects(&self, o: &BoundingBox) -> bool {
        (0..3).all(|i| self.min[i] <= o.max[i] && o.min[i] <= self.max[i])
    }

    /// Size along each axis.
    pub fn size(&self) -> Vec3 {
        self.max - self.min
    }

    /// Center point.
    pub fn center(&self) -> Vec3 {
        (self.min + self.max) * 0.5
    }

    /// The eight corners of the box.
    pub fn corners(&self) -> [Vec3; 8] {
        let (a, b) = (self.min, self.max);
        core::array::from_fn(|i| {
            Vec3::new(
                if i & 1 == 0 { a[0] } else { b[0] },
                if i & 2 == 0 { a[1] } else { b[1] },
                if i & 4 == 0 { a[2] } else { b[2] },
            )
        })
    }

    /// Bounding box of this box after an affine transformation.
    pub fn transformed(&self, m: &Mat4) -> BoundingBox {
        if self.is_empty() {
            return *self;
        }
        BoundingBox::from_points(self.corners().map(|c| m.transform_point(c)))
    }

    /// Number of axes along which the box has non-zero extent.
    pub fn dimension(&self) -> usize {
        (0..3).filter(|&i| self.max[i] > self.min[i]).count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: Vec3, b: Vec3) -> bool {
        (a - b).norm() < 1e-12
    }

    #[test]
    fn vector_ops() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        assert_eq!(a + b, Vec3::new(5.0, 7.0, 9.0));
        assert_eq!(a.dot(b), 32.0);
        assert_eq!(Vec3::unit(0).cross(Vec3::unit(1)), Vec3::unit(2));
    }

    #[test]
    fn rotation_and_inverse() {
        let r = Mat4::rotation(Vec3::unit(2), core::f64::consts::FRAC_PI_2);
        assert!(close(r.transform_point(Vec3::unit(0)), Vec3::unit(1)));
        let m = Mat4::translation(Vec3::new(1.0, 2.0, 3.0))
            .mul(&r)
            .mul(&Mat4::scaling(Vec3::new(2.0, 3.0, 4.0)));
        let inv = m.inverse().unwrap();
        let p = Vec3::new(0.3, -1.2, 7.0);
        assert!(close(inv.transform_point(m.transform_point(p)), p));
        assert!(Mat4::scaling(Vec3::new(1.0, 0.0, 1.0)).inverse().is_none());
    }

    #[test]
    fn bbox() {
        let b = BoundingBox::from_corners(Vec3::new(1.0, 1.0, 1.0), Vec3::new(0.0, 2.0, 1.0));
        assert_eq!(b.dimension(), 2);
        assert!(b.contains(Vec3::new(0.5, 1.5, 1.0)));
        assert!(BoundingBox::EMPTY.is_empty());
        assert_eq!(BoundingBox::EMPTY.union(&b), b);
    }
}
