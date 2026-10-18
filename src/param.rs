//! Parameters and parametric scalar values.
//!
//! Most numbers in a CSXCAD model may be given either as a constant or as an
//! expression referring to named parameters (e.g. `"substrate_h/2"`). Such a
//! value is a [`Scalar`]: it remembers the expression and caches its last
//! evaluated value. [`crate::ContinuousStructure::update`] re-evaluates every
//! scalar against the structure's [`ParameterSet`].

use core::fmt;

use crate::error::{Error, Result};
use crate::expr::{Expr, parse_number, parse_number_lenient};

/// How a parameter is swept (informational; used by sweep front-ends).
#[derive(Debug, Clone, PartialEq, Default)]
pub enum ParameterKind {
    /// A constant value.
    #[default]
    Const,
    /// A linearly swept value: `start, start+step, ..., stop`.
    Linear {
        /// First value.
        start: f64,
        /// Last value.
        stop: f64,
        /// Increment.
        step: f64,
    },
}

/// A named parameter.
#[derive(Debug, Clone, PartialEq)]
pub struct Parameter {
    /// Name, usable in expressions.
    pub name: String,
    /// Current value.
    pub value: f64,
    /// Whether the parameter takes part in sweeps.
    pub sweep: bool,
    /// Sweep definition.
    pub kind: ParameterKind,
}

impl Parameter {
    /// A constant parameter.
    pub fn new(name: impl Into<String>, value: f64) -> Self {
        Parameter {
            name: name.into(),
            value,
            sweep: true,
            kind: ParameterKind::Const,
        }
    }
}

/// An ordered collection of uniquely named parameters.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ParameterSet {
    params: Vec<Parameter>,
}

impl ParameterSet {
    /// Empty set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Inserts or replaces a parameter (by name).
    pub fn insert(&mut self, p: Parameter) {
        match self.params.iter_mut().find(|q| q.name == p.name) {
            Some(q) => *q = p,
            None => self.params.push(p),
        }
    }

    /// Sets the value of a parameter, creating a constant one if needed.
    pub fn set(&mut self, name: &str, value: f64) {
        match self.params.iter_mut().find(|q| q.name == name) {
            Some(q) => q.value = value,
            None => self.params.push(Parameter::new(name, value)),
        }
    }

    /// Value of a parameter.
    pub fn get(&self, name: &str) -> Option<f64> {
        self.params.iter().find(|p| p.name == name).map(|p| p.value)
    }

    /// Parameter by name.
    pub fn parameter(&self, name: &str) -> Option<&Parameter> {
        self.params.iter().find(|p| p.name == name)
    }

    /// Removes a parameter, returning it.
    pub fn remove(&mut self, name: &str) -> Option<Parameter> {
        let i = self.params.iter().position(|p| p.name == name)?;
        Some(self.params.remove(i))
    }

    /// Iterates over all parameters in insertion order.
    pub fn iter(&self) -> impl Iterator<Item = &Parameter> {
        self.params.iter()
    }

    /// Number of parameters.
    pub fn len(&self) -> usize {
        self.params.len()
    }

    /// True if there are no parameters.
    pub fn is_empty(&self) -> bool {
        self.params.is_empty()
    }

    /// Evaluates an expression string against this set.
    pub fn eval(&self, expr: &str) -> Result<f64> {
        if let Some(v) = parse_number(expr) {
            return Ok(v);
        }
        Expr::parse(expr)?.eval_with(|n| self.get(n))
    }
}

/// A number that is either a constant or an expression over parameters.
///
/// The cached [`Scalar::value`] is what geometric queries use; for
/// expressions it is refreshed by [`Scalar::update`].
#[derive(Debug, Clone, PartialEq)]
pub struct Scalar {
    expr: Option<Box<str>>,
    value: f64,
}

impl Default for Scalar {
    fn default() -> Self {
        Scalar::new(0.0)
    }
}

impl Scalar {
    /// A constant.
    pub const fn new(value: f64) -> Self {
        Scalar { expr: None, value }
    }

    /// An expression. Plain numbers are stored as constants. The cached value
    /// is computed immediately when the expression does not depend on any
    /// parameter, otherwise it is `NaN` until [`Scalar::update`] is called.
    pub fn expr(src: &str) -> Self {
        if let Some(v) = parse_number_lenient(src) {
            return Scalar::new(v);
        }
        let value = Expr::parse(src)
            .and_then(|e| e.eval_const())
            .unwrap_or(f64::NAN);
        Scalar {
            expr: Some(src.trim().into()),
            value,
        }
    }

    /// Like [`Scalar::expr`] but validates the syntax.
    pub fn try_expr(src: &str) -> Result<Self> {
        if parse_number_lenient(src).is_none() {
            Expr::parse(src)?;
        }
        Ok(Scalar::expr(src))
    }

    /// Cached numeric value.
    #[inline]
    pub fn value(&self) -> f64 {
        self.value
    }

    /// Expression source, if this scalar is parametric.
    pub fn expression(&self) -> Option<&str> {
        self.expr.as_deref()
    }

    /// True for parametric scalars.
    pub fn is_expr(&self) -> bool {
        self.expr.is_some()
    }

    /// Replaces the value with a constant.
    pub fn set(&mut self, v: f64) {
        *self = Scalar::new(v);
    }

    /// Re-evaluates the expression (if any) against `params`.
    pub fn update(&mut self, params: &ParameterSet) -> Result<()> {
        if let Some(e) = &self.expr {
            let v = params.eval(e)?;
            if v.is_nan() {
                return Err(Error::Invalid(format!("expression {e:?} evaluates to NaN")));
            }
            self.value = v;
        }
        Ok(())
    }

    /// The textual form used in files: the expression or the number.
    pub fn to_text(&self) -> String {
        match &self.expr {
            Some(e) => e.to_string(),
            None => format_f64(self.value),
        }
    }
}

impl From<f64> for Scalar {
    fn from(v: f64) -> Self {
        Scalar::new(v)
    }
}

impl From<&str> for Scalar {
    fn from(s: &str) -> Self {
        Scalar::expr(s)
    }
}

impl fmt::Display for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_text())
    }
}

/// Formats a float with the shortest representation that round-trips.
pub fn format_f64(v: f64) -> String {
    if v.is_nan() {
        return "nan".into();
    }
    if v.is_infinite() {
        return if v > 0.0 { "inf".into() } else { "-inf".into() };
    }
    if v == 0.0 {
        return "0".into();
    }
    let s = format!("{v}");
    // Very large/small magnitudes are more readable in exponent notation.
    let e = format!("{v:e}");
    if e.len() < s.len() { e } else { s }
}

/// A point whose three coordinates are [`Scalar`]s.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ScalarPoint(pub [Scalar; 3]);

impl ScalarPoint {
    /// Constant point.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        ScalarPoint([x.into(), y.into(), z.into()])
    }

    /// Cached numeric value.
    #[inline]
    pub fn value(&self) -> crate::Vec3 {
        crate::Vec3::new(self.0[0].value(), self.0[1].value(), self.0[2].value())
    }

    /// Re-evaluates all coordinates.
    pub fn update(&mut self, params: &ParameterSet) -> Result<()> {
        self.0.iter_mut().try_for_each(|s| s.update(params))
    }
}

impl From<crate::Vec3> for ScalarPoint {
    fn from(v: crate::Vec3) -> Self {
        ScalarPoint::new(v[0], v[1], v[2])
    }
}

impl From<[f64; 3]> for ScalarPoint {
    fn from(v: [f64; 3]) -> Self {
        ScalarPoint::new(v[0], v[1], v[2])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scalars() {
        let mut ps = ParameterSet::new();
        ps.set("h", 1.5);
        let mut s = Scalar::expr("h*2");
        assert!(s.value().is_nan());
        s.update(&ps).unwrap();
        assert_eq!(s.value(), 3.0);
        assert_eq!(s.to_text(), "h*2");
        let c = Scalar::expr(" 2.5 ");
        assert!(!c.is_expr());
        assert_eq!(Scalar::expr("2*3").value(), 6.0);
        assert!(Scalar::expr("q").update(&ps).is_err());
        assert!(Scalar::try_expr("1+").is_err());
    }

    #[test]
    fn formatting() {
        assert_eq!(format_f64(0.5), "0.5");
        assert_eq!(format_f64(1e-12), "1e-12");
        assert_eq!(format_f64(100.0), "100");
        assert_eq!(format_f64(299792458.0), "299792458");
        for v in [0.1, 1.0 / 3.0, 6.02e23, -1e-300] {
            assert_eq!(format_f64(v).parse::<f64>().unwrap(), v);
        }
    }
}
