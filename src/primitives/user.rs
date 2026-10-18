//! Primitives defined by a user expression.

use crate::coords::cartesian_to_spherical;
use crate::error::Result;
use crate::expr::{CompiledExpr, Expr};
use crate::math::{BoundingBox, Vec3};
use crate::param::{ParameterSet, ScalarPoint};

/// Variable names available in [`UserDefined`] functions, in slot order.
pub const USER_VARS: [&str; 7] = ["x", "y", "z", "rho", "a", "r", "t"];

/// A primitive whose interior is the set of points where `function`
/// evaluates to a non-zero (and non-NaN) value.
///
/// Variables: `x, y, z` (Cartesian, relative to `base`), `rho` and `a`
/// (cylindrical radius and azimuth), `r` and `t` (spherical radius and polar
/// angle). The extent is unbounded.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct UserDefined {
    /// The defining expression.
    pub function: String,
    /// Origin of the local coordinates.
    pub base: ScalarPoint,
    compiled: Option<CompiledExpr>,
}

impl UserDefined {
    /// Creates the primitive; parameter-free functions are compiled
    /// immediately, others when the owning structure is updated.
    pub fn new(function: impl Into<String>) -> Self {
        let function = function.into();
        let compiled = Expr::parse(&function)
            .and_then(|e| e.compile(&USER_VARS, |_| None))
            .ok();
        UserDefined {
            function,
            base: ScalarPoint::default(),
            compiled,
        }
    }

    pub(crate) fn inside(&self, p: Vec3) -> bool {
        let Some(f) = &self.compiled else {
            return false;
        };
        let q = p - self.base.value();
        let (r, t, a) = cartesian_to_spherical(q);
        let rho = q[0].hypot(q[1]);
        let v = f.eval(&[q[0], q[1], q[2], rho, a, r, t]);
        v != 0.0 && !v.is_nan()
    }

    pub(crate) fn bbox(&self) -> BoundingBox {
        BoundingBox {
            min: Vec3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
            max: Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
        }
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.base.update(p)?;
        self.compiled = Some(Expr::parse(&self.function)?.compile(&USER_VARS, |n| p.get(n))?);
        Ok(())
    }
}
