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
///
/// The function is private so its compiled form can never go stale; change
/// it with [`UserDefined::set_function`]. Functions that reference
/// parameters are compiled when the owning structure is updated.
///
/// CSXCAD 0.7 itself does not support this primitive (it drops it when
/// reading); this crate reads and writes it as `<UserDefined>` with an
/// `<Equation>` text and a `<CoordShift X Y Z/>` origin.
///
/// ```
/// use csx::primitives::{Primitive, UserDefined};
/// use csx::{CoordSystem, Vec3};
///
/// // A cone with its apex at (0, 0, 5), opening downwards.
/// let cone = UserDefined::new("rho < -z & z < 0 & z > -5").with_base([0.0, 0.0, 5.0]);
/// assert!(cone.is_compiled());
/// let p = Primitive::new(cone);
/// assert!(p.contains(Vec3::new(0.5, 0.0, 2.0), CoordSystem::Cartesian));
/// assert!(!p.contains(Vec3::new(4.0, 0.0, 2.0), CoordSystem::Cartesian));
/// // The extent is unbounded.
/// assert!(p.bounding_box(CoordSystem::Cartesian).max.x().is_infinite());
/// ```
#[derive(Debug, Clone, Default)]
pub struct UserDefined {
    function: String,
    /// Origin of the local coordinates.
    pub base: ScalarPoint,
    compiled: Option<CompiledExpr>,
}

impl PartialEq for UserDefined {
    /// Compares the definitions (function text and base point).
    fn eq(&self, o: &Self) -> bool {
        self.function == o.function && self.base == o.base
    }
}

impl UserDefined {
    /// Creates the primitive; parameter-free functions are compiled
    /// immediately, others when the owning structure is updated.
    pub fn new(function: impl Into<String>) -> Self {
        let mut u = UserDefined::default();
        u.set_function(function);
        u
    }

    /// Builder: sets the origin of the local coordinates.
    #[must_use]
    pub fn with_base(mut self, base: impl Into<ScalarPoint>) -> Self {
        self.base = base.into();
        self
    }

    /// The defining expression.
    pub fn function(&self) -> &str {
        &self.function
    }

    /// Replaces the defining expression. Parameter-free functions are
    /// compiled immediately; others are compiled by the next update, and
    /// until then the primitive contains no point.
    pub fn set_function(&mut self, function: impl Into<String>) {
        self.function = function.into();
        self.compiled = Expr::parse(&self.function)
            .and_then(|e| e.compile(&USER_VARS, |_| None))
            .ok();
    }

    /// True once the function has been compiled successfully.
    pub fn is_compiled(&self) -> bool {
        self.compiled.is_some()
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
