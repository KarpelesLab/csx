//! Position-dependent weighting functions (re-exported from
//! [`crate::properties`]).

use crate::coords::cartesian_to_spherical;
use crate::error::Result;
use crate::expr::{CompiledExpr, Expr, parse_number};
use crate::math::Vec3;
use crate::param::{ParameterSet, Scalar};

/// Variables available in weighting functions, in slot order: Cartesian
/// `x, y, z`, cylindrical `rho, a` (azimuth) and spherical `r, t` (polar
/// angle).
pub const WEIGHT_VARS: [&str; 7] = ["x", "y", "z", "rho", "a", "r", "t"];

/// Computes the weighting-function variable slots ([`WEIGHT_VARS`] order)
/// for a Cartesian point.
///
/// ```
/// use csx::properties::weight_vars;
/// use csx::Vec3;
///
/// let [x, y, z, rho, a, r, t] = weight_vars(Vec3::new(0.0, 3.0, 4.0));
/// assert_eq!((x, y, z, rho, r), (0.0, 3.0, 4.0, 3.0, 5.0));
/// assert!((a - std::f64::consts::FRAC_PI_2).abs() < 1e-15);
/// assert!((t - (3.0f64).atan2(4.0)).abs() < 1e-15);
/// ```
#[inline]
pub fn weight_vars(p: Vec3) -> [f64; 7] {
    let (r, t, a) = cartesian_to_spherical(p);
    [p[0], p[1], p[2], p[0].hypot(p[1]), a, r, t]
}

/// A weighting function: an expression over [`WEIGHT_VARS`] that multiplies
/// a base value. The default function is the constant `1`.
///
/// The expression may also reference parameters; those are resolved by
/// [`WeightFn::update`] (normally through
/// [`ContinuousStructure::update`](crate::ContinuousStructure::update)).
/// Until then, and if the expression is malformed, the function evaluates
/// to `1`.
///
/// ```
/// use csx::properties::WeightFn;
/// use csx::{ParameterSet, Vec3};
///
/// let w = WeightFn::new("exp(-(x^2 + y^2))");
/// assert_eq!(w.eval(Vec3::ZERO), 1.0);
/// assert!(w.eval(Vec3::new(1.0, 0.0, 0.0)) < 0.37);
///
/// let mut taper = WeightFn::new("1 - rho/R");
/// let mut params = ParameterSet::new();
/// params.set("R", 4.0);
/// taper.update(&params)?;
/// assert_eq!(taper.eval(Vec3::new(0.0, 2.0, 0.0)), 0.5);
/// # Ok::<(), csx::Error>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct WeightFn {
    src: String,
    constant: Option<f64>,
    compiled: Option<CompiledExpr>,
}

impl Default for WeightFn {
    fn default() -> Self {
        WeightFn {
            src: "1".into(),
            constant: Some(1.0),
            compiled: None,
        }
    }
}

impl WeightFn {
    /// Creates a weighting function; parameter-free expressions are compiled
    /// immediately, others on [`WeightFn::update`].
    pub fn new(src: impl Into<String>) -> Self {
        let src = src.into();
        if let Some(v) = parse_number(&src) {
            return WeightFn {
                src,
                constant: Some(v),
                compiled: None,
            };
        }
        let compiled = Expr::parse(&src)
            .and_then(|e| e.compile(&WEIGHT_VARS, |_| None))
            .ok();
        let constant = compiled.as_ref().and_then(|c| c.as_constant());
        WeightFn {
            src,
            constant,
            compiled,
        }
    }

    /// Expression source.
    pub fn source(&self) -> &str {
        &self.src
    }

    /// True for the default (constant 1) function.
    pub fn is_unity(&self) -> bool {
        self.constant == Some(1.0)
    }

    /// Evaluates at precomputed variables (see [`weight_vars`]). Functions
    /// that failed to compile evaluate to `1`.
    #[inline]
    pub fn eval_vars(&self, vars: &[f64; 7]) -> f64 {
        if let Some(c) = self.constant {
            return c;
        }
        self.compiled.as_ref().map_or(1.0, |c| c.eval(vars))
    }

    /// Evaluates at a Cartesian point.
    pub fn eval(&self, p: Vec3) -> f64 {
        if let Some(c) = self.constant {
            return c;
        }
        self.eval_vars(&weight_vars(p))
    }

    /// Recompiles the expression against `params`.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`](crate::Error::Parse) for a malformed expression,
    /// [`Error::UnknownVariable`](crate::Error::UnknownVariable) for a name
    /// that is neither a weighting variable nor a parameter.
    pub fn update(&mut self, params: &ParameterSet) -> Result<()> {
        if parse_number(&self.src).is_some() {
            return Ok(());
        }
        let c = Expr::parse(&self.src)?.compile(&WEIGHT_VARS, |n| params.get(n))?;
        self.constant = c.as_constant();
        self.compiled = Some(c);
        Ok(())
    }
}

impl From<&str> for WeightFn {
    fn from(s: &str) -> Self {
        WeightFn::new(s)
    }
}

/// A possibly anisotropic quantity: one value and one weighting function per
/// axis.
///
/// ```
/// use csx::properties::{Aniso, weight_vars};
/// use csx::{Scalar, Vec3};
///
/// let mut a = Aniso::iso(2.0);
/// assert!(a.is_isotropic() && !a.is_weighted());
/// a.value[2] = Scalar::new(3.0);
/// a.set_weight("1 + z");
/// assert!(!a.is_isotropic() && a.is_weighted());
/// let vars = weight_vars(Vec3::new(0.0, 0.0, 1.0));
/// assert_eq!(a.at_vars(0, &vars), 4.0);
/// assert_eq!(a.at_vars(2, &vars), 6.0);
/// assert_eq!(a.base(2), 3.0);
/// ```
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Aniso {
    /// Base values per axis.
    pub value: [Scalar; 3],
    /// Weighting functions per axis.
    pub weight: [WeightFn; 3],
}

impl Aniso {
    /// Isotropic quantity with value `v` on every axis.
    pub fn iso(v: f64) -> Self {
        Aniso {
            value: [Scalar::new(v), Scalar::new(v), Scalar::new(v)],
            weight: Default::default(),
        }
    }

    /// Sets the same value on every axis.
    pub fn set(&mut self, v: impl Into<Scalar>) {
        let v = v.into();
        self.value = [v.clone(), v.clone(), v];
    }

    /// Sets the same weighting function on every axis.
    pub fn set_weight(&mut self, w: &str) {
        self.weight = [WeightFn::new(w), WeightFn::new(w), WeightFn::new(w)];
    }

    /// Unweighted value along `dir` (0..=2).
    #[inline]
    pub fn base(&self, dir: usize) -> f64 {
        self.value[dir].value()
    }

    /// Weighted value along `dir` at the given weighting variables.
    #[inline]
    pub fn at_vars(&self, dir: usize, vars: &[f64; 7]) -> f64 {
        self.value[dir].value() * self.weight[dir].eval_vars(vars)
    }

    /// True if all three axes share the same value and weight.
    pub fn is_isotropic(&self) -> bool {
        self.value[0] == self.value[1]
            && self.value[1] == self.value[2]
            && self.weight[0] == self.weight[1]
            && self.weight[1] == self.weight[2]
    }

    /// True if any axis carries a non-trivial weighting function.
    pub fn is_weighted(&self) -> bool {
        self.weight.iter().any(|w| !w.is_unity())
    }

    pub(crate) fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.value.iter_mut().try_for_each(|v| v.update(p))?;
        self.weight.iter_mut().try_for_each(|w| w.update(p))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weights() {
        let w = WeightFn::new("2*x + rho");
        assert_eq!(w.eval(Vec3::new(3.0, 4.0, 0.0)), 11.0);
        let mut w = WeightFn::new("k*z");
        assert_eq!(w.eval(Vec3::new(0.0, 0.0, 2.0)), 1.0);
        let mut p = ParameterSet::new();
        p.set("k", 3.0);
        w.update(&p).unwrap();
        assert_eq!(w.eval(Vec3::new(0.0, 0.0, 2.0)), 6.0);
        assert!(WeightFn::default().is_unity());
        let mut a = Aniso::iso(4.0);
        a.weight[2] = WeightFn::new("0.5");
        assert!(!a.is_isotropic());
        assert_eq!(a.at_vars(2, &weight_vars(Vec3::ZERO)), 2.0);
    }
}
