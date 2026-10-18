//! Mathematical expressions.
//!
//! The accepted syntax follows the conventions of the classic "function
//! parser" dialect used by CSXCAD files:
//!
//! * numbers: `1`, `2.5`, `.5`, `1e-3`
//! * operators, lowest to highest precedence:
//!   `|`, `&`, comparisons (`=`, `==`, `!=`, `<`, `<=`, `>`, `>=`), `+ -`,
//!   `* / %`, unary `-`/`!`, `^` (right associative, binds tighter than unary
//!   minus so `-2^2 = -4`)
//! * functions: `abs acos acosh asin asinh atan atan2 atanh cbrt ceil cos
//!   cosh cot csc exp exp2 floor hypot if int log log2 log10 max min not pow
//!   sec sin sinh sqrt tan tanh trunc`
//! * constants `pi` and `e` (overridable by variables / parameters)
//!
//! Booleans are `1` / `0`; any non-zero value is *true*.

use std::collections::BTreeSet;

use crate::error::{Error, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BinOp {
    Or,
    And,
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
    Add,
    Sub,
    Mul,
    Div,
    Mod,
    Pow,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum UnOp {
    Neg,
    Not,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Func {
    Abs,
    Acos,
    Acosh,
    Asin,
    Asinh,
    Atan,
    Atan2,
    Atanh,
    Cbrt,
    Ceil,
    Cos,
    Cosh,
    Cot,
    Csc,
    Exp,
    Exp2,
    Floor,
    Hypot,
    If,
    Int,
    Log,
    Log2,
    Log10,
    Max,
    Min,
    Not,
    Pow,
    Sec,
    Sin,
    Sinh,
    Sqrt,
    Tan,
    Tanh,
    Trunc,
}

impl Func {
    fn lookup(name: &str) -> Option<(Func, usize)> {
        use Func::*;
        Some(match name {
            "abs" => (Abs, 1),
            "acos" => (Acos, 1),
            "acosh" => (Acosh, 1),
            "asin" => (Asin, 1),
            "asinh" => (Asinh, 1),
            "atan" => (Atan, 1),
            "atan2" => (Atan2, 2),
            "atanh" => (Atanh, 1),
            "cbrt" => (Cbrt, 1),
            "ceil" => (Ceil, 1),
            "cos" => (Cos, 1),
            "cosh" => (Cosh, 1),
            "cot" => (Cot, 1),
            "csc" => (Csc, 1),
            "exp" => (Exp, 1),
            "exp2" => (Exp2, 1),
            "floor" => (Floor, 1),
            "hypot" => (Hypot, 2),
            "if" => (If, 3),
            "int" => (Int, 1),
            "log" => (Log, 1),
            "log2" => (Log2, 1),
            "log10" => (Log10, 1),
            "max" => (Max, 2),
            "min" => (Min, 2),
            "not" => (Not, 1),
            "pow" => (Pow, 2),
            "sec" => (Sec, 1),
            "sin" => (Sin, 1),
            "sinh" => (Sinh, 1),
            "sqrt" => (Sqrt, 1),
            "tan" => (Tan, 1),
            "tanh" => (Tanh, 1),
            "trunc" => (Trunc, 1),
            _ => return None,
        })
    }

    fn apply(self, a: &[f64]) -> f64 {
        use Func::*;
        match self {
            Abs => a[0].abs(),
            Acos => a[0].acos(),
            Acosh => a[0].acosh(),
            Asin => a[0].asin(),
            Asinh => a[0].asinh(),
            Atan => a[0].atan(),
            Atan2 => a[0].atan2(a[1]),
            Atanh => a[0].atanh(),
            Cbrt => a[0].cbrt(),
            Ceil => a[0].ceil(),
            Cos => a[0].cos(),
            Cosh => a[0].cosh(),
            Cot => 1.0 / a[0].tan(),
            Csc => 1.0 / a[0].sin(),
            Exp => a[0].exp(),
            Exp2 => a[0].exp2(),
            Floor => a[0].floor(),
            Hypot => a[0].hypot(a[1]),
            If => unreachable!("if() is evaluated lazily"),
            Int => a[0].round(),
            Log => a[0].ln(),
            Log2 => a[0].log2(),
            Log10 => a[0].log10(),
            Max => a[0].max(a[1]),
            Min => a[0].min(a[1]),
            Not => bool_val(!truthy(a[0])),
            Pow => a[0].powf(a[1]),
            Sec => 1.0 / a[0].cos(),
            Sin => a[0].sin(),
            Sinh => a[0].sinh(),
            Sqrt => a[0].sqrt(),
            Tan => a[0].tan(),
            Tanh => a[0].tanh(),
            Trunc => a[0].trunc(),
        }
    }
}

#[inline]
fn truthy(v: f64) -> bool {
    v != 0.0
}

#[inline]
fn bool_val(b: bool) -> f64 {
    if b { 1.0 } else { 0.0 }
}

#[derive(Debug, Clone, PartialEq)]
enum Node {
    Num(f64),
    Name(String),
    Slot(usize),
    Unary(UnOp, Box<Node>),
    Binary(BinOp, Box<Node>, Box<Node>),
    Call(Func, Vec<Node>),
}

impl Node {
    fn eval(&self, slots: &[f64], lookup: &mut dyn FnMut(&str) -> Option<f64>) -> Result<f64> {
        Ok(match self {
            Node::Num(v) => *v,
            Node::Slot(i) => slots[*i],
            Node::Name(n) => match lookup(n) {
                Some(v) => v,
                None => builtin_constant(n).ok_or_else(|| Error::UnknownVariable(n.clone()))?,
            },
            Node::Unary(op, a) => {
                let a = a.eval(slots, lookup)?;
                match op {
                    UnOp::Neg => -a,
                    UnOp::Not => bool_val(!truthy(a)),
                }
            }
            Node::Binary(op, a, b) => {
                let a = a.eval(slots, lookup)?;
                // Short-circuit logic operators.
                match op {
                    BinOp::And if !truthy(a) => return Ok(0.0),
                    BinOp::Or if truthy(a) => return Ok(1.0),
                    _ => {}
                }
                let b = b.eval(slots, lookup)?;
                binop(*op, a, b)
            }
            Node::Call(Func::If, args) => {
                if truthy(args[0].eval(slots, lookup)?) {
                    args[1].eval(slots, lookup)?
                } else {
                    args[2].eval(slots, lookup)?
                }
            }
            Node::Call(f, args) => {
                let mut vals = [0.0; 3];
                for (v, a) in vals.iter_mut().zip(args) {
                    *v = a.eval(slots, lookup)?;
                }
                f.apply(&vals[..args.len()])
            }
        })
    }

    /// Evaluates a node that contains only numbers and slots.
    fn eval_fast(&self, slots: &[f64]) -> f64 {
        match self {
            Node::Num(v) => *v,
            Node::Slot(i) => slots[*i],
            Node::Name(_) => f64::NAN,
            Node::Unary(UnOp::Neg, a) => -a.eval_fast(slots),
            Node::Unary(UnOp::Not, a) => bool_val(!truthy(a.eval_fast(slots))),
            Node::Binary(op, a, b) => {
                let a = a.eval_fast(slots);
                match op {
                    BinOp::And if !truthy(a) => 0.0,
                    BinOp::Or if truthy(a) => 1.0,
                    _ => binop(*op, a, b.eval_fast(slots)),
                }
            }
            Node::Call(Func::If, args) => {
                if truthy(args[0].eval_fast(slots)) {
                    args[1].eval_fast(slots)
                } else {
                    args[2].eval_fast(slots)
                }
            }
            Node::Call(f, args) => {
                let mut vals = [0.0; 3];
                for (v, a) in vals.iter_mut().zip(args) {
                    *v = a.eval_fast(slots);
                }
                f.apply(&vals[..args.len()])
            }
        }
    }

    fn collect_names<'a>(&'a self, out: &mut BTreeSet<&'a str>) {
        match self {
            Node::Name(n) => {
                out.insert(n);
            }
            Node::Unary(_, a) => a.collect_names(out),
            Node::Binary(_, a, b) => {
                a.collect_names(out);
                b.collect_names(out);
            }
            Node::Call(_, args) => args.iter().for_each(|a| a.collect_names(out)),
            Node::Num(_) | Node::Slot(_) => {}
        }
    }

    /// Replaces names by slots or constants, folding constant sub-trees.
    fn resolve(&self, vars: &[&str], lookup: &mut dyn FnMut(&str) -> Option<f64>) -> Result<Node> {
        Ok(match self {
            Node::Num(v) => Node::Num(*v),
            Node::Slot(i) => Node::Slot(*i),
            Node::Name(n) => {
                if let Some(i) = vars.iter().position(|v| v == n) {
                    Node::Slot(i)
                } else if let Some(v) = lookup(n).or_else(|| builtin_constant(n)) {
                    Node::Num(v)
                } else {
                    return Err(Error::UnknownVariable(n.clone()));
                }
            }
            Node::Unary(op, a) => fold(Node::Unary(*op, Box::new(a.resolve(vars, lookup)?))),
            Node::Binary(op, a, b) => fold(Node::Binary(
                *op,
                Box::new(a.resolve(vars, lookup)?),
                Box::new(b.resolve(vars, lookup)?),
            )),
            Node::Call(f, args) => fold(Node::Call(
                *f,
                args.iter()
                    .map(|a| a.resolve(vars, lookup))
                    .collect::<Result<_>>()?,
            )),
        })
    }
}

fn fold(n: Node) -> Node {
    let constant = match &n {
        Node::Unary(_, a) => matches!(**a, Node::Num(_)),
        Node::Binary(_, a, b) => matches!(**a, Node::Num(_)) && matches!(**b, Node::Num(_)),
        Node::Call(_, args) => args.iter().all(|a| matches!(a, Node::Num(_))),
        _ => false,
    };
    if constant {
        Node::Num(n.eval_fast(&[]))
    } else {
        n
    }
}

fn binop(op: BinOp, a: f64, b: f64) -> f64 {
    match op {
        BinOp::Or => bool_val(truthy(a) || truthy(b)),
        BinOp::And => bool_val(truthy(a) && truthy(b)),
        BinOp::Eq => bool_val(a == b),
        BinOp::Ne => bool_val(a != b),
        BinOp::Lt => bool_val(a < b),
        BinOp::Le => bool_val(a <= b),
        BinOp::Gt => bool_val(a > b),
        BinOp::Ge => bool_val(a >= b),
        BinOp::Add => a + b,
        BinOp::Sub => a - b,
        BinOp::Mul => a * b,
        BinOp::Div => a / b,
        BinOp::Mod => a % b,
        BinOp::Pow => a.powf(b),
    }
}

fn builtin_constant(name: &str) -> Option<f64> {
    match name {
        "pi" => Some(core::f64::consts::PI),
        "e" => Some(core::f64::consts::E),
        _ => None,
    }
}

/// A parsed expression.
#[derive(Debug, Clone, PartialEq)]
pub struct Expr {
    src: String,
    root: Node,
}

impl Expr {
    /// Parses `src`.
    pub fn parse(src: &str) -> Result<Expr> {
        let mut p = Parser { src, pos: 0 };
        let root = p.parse_expr(0)?;
        p.skip_ws();
        if p.pos != src.len() {
            return Err(p.err("unexpected trailing input"));
        }
        Ok(Expr {
            src: src.to_owned(),
            root,
        })
    }

    /// Original source text.
    pub fn source(&self) -> &str {
        &self.src
    }

    /// Names referenced by the expression (excluding function names).
    pub fn names(&self) -> BTreeSet<&str> {
        let mut s = BTreeSet::new();
        self.root.collect_names(&mut s);
        s
    }

    /// Evaluates the expression; names are resolved by `lookup`, falling
    /// back to the built-in constants.
    pub fn eval_with(&self, mut lookup: impl FnMut(&str) -> Option<f64>) -> Result<f64> {
        self.root.eval(&[], &mut lookup)
    }

    /// Evaluates an expression without any variables.
    pub fn eval_const(&self) -> Result<f64> {
        self.eval_with(|_| None)
    }

    /// Compiles the expression for fast repeated evaluation. Names listed in
    /// `vars` become positional arguments of [`CompiledExpr::eval`]; all other
    /// names are resolved once through `lookup` (then built-in constants).
    pub fn compile(
        &self,
        vars: &[&str],
        mut lookup: impl FnMut(&str) -> Option<f64>,
    ) -> Result<CompiledExpr> {
        Ok(CompiledExpr {
            root: self.root.resolve(vars, &mut lookup)?,
            arity: vars.len(),
        })
    }
}

impl core::str::FromStr for Expr {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self> {
        Expr::parse(s)
    }
}

impl core::fmt::Display for Expr {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(&self.src)
    }
}

/// An expression whose free names were bound to positional slots or folded
/// to constants; see [`Expr::compile`].
#[derive(Debug, Clone, PartialEq)]
pub struct CompiledExpr {
    root: Node,
    arity: usize,
}

impl CompiledExpr {
    /// Evaluates with the given variable values (in the order passed to
    /// [`Expr::compile`]). Missing values are treated as `0`.
    pub fn eval(&self, vars: &[f64]) -> f64 {
        if vars.len() >= self.arity {
            self.root.eval_fast(vars)
        } else {
            let mut v = vars.to_vec();
            v.resize(self.arity, 0.0);
            self.root.eval_fast(&v)
        }
    }

    /// The constant value, if the expression does not depend on any variable.
    pub fn as_constant(&self) -> Option<f64> {
        match self.root {
            Node::Num(v) => Some(v),
            _ => None,
        }
    }
}

struct Parser<'a> {
    src: &'a str,
    pos: usize,
}

impl Parser<'_> {
    fn err(&self, msg: &str) -> Error {
        Error::Parse {
            expr: self.src.to_owned(),
            pos: self.pos,
            msg: msg.to_owned(),
        }
    }

    fn skip_ws(&mut self) {
        while let Some(c) = self.peek() {
            if c.is_ascii_whitespace() {
                self.pos += 1;
            } else {
                break;
            }
        }
    }

    fn peek(&self) -> Option<u8> {
        self.src.as_bytes().get(self.pos).copied()
    }

    fn peek2(&self) -> Option<u8> {
        self.src.as_bytes().get(self.pos + 1).copied()
    }

    /// Reads a binary operator (without consuming it) and its binding power.
    fn peek_binop(&mut self) -> Option<(BinOp, u8, usize)> {
        self.skip_ws();
        let c = self.peek()?;
        let c2 = self.peek2();
        Some(match (c, c2) {
            (b'|', _) => (BinOp::Or, 1, 1),
            (b'&', _) => (BinOp::And, 2, 1),
            (b'=', Some(b'=')) => (BinOp::Eq, 3, 2),
            (b'=', _) => (BinOp::Eq, 3, 1),
            (b'!', Some(b'=')) => (BinOp::Ne, 3, 2),
            (b'<', Some(b'=')) => (BinOp::Le, 3, 2),
            (b'<', _) => (BinOp::Lt, 3, 1),
            (b'>', Some(b'=')) => (BinOp::Ge, 3, 2),
            (b'>', _) => (BinOp::Gt, 3, 1),
            (b'+', _) => (BinOp::Add, 4, 1),
            (b'-', _) => (BinOp::Sub, 4, 1),
            (b'*', _) => (BinOp::Mul, 5, 1),
            (b'/', _) => (BinOp::Div, 5, 1),
            (b'%', _) => (BinOp::Mod, 5, 1),
            (b'^', _) => (BinOp::Pow, 7, 1),
            _ => return None,
        })
    }

    fn parse_expr(&mut self, min_bp: u8) -> Result<Node> {
        let mut lhs = self.parse_unary()?;
        while let Some((op, bp, len)) = self.peek_binop() {
            if bp < min_bp {
                break;
            }
            self.pos += len;
            // `^` is right associative; everything else left associative.
            let rhs = if op == BinOp::Pow {
                self.parse_pow_rhs()?
            } else {
                self.parse_expr(bp + 1)?
            };
            lhs = Node::Binary(op, Box::new(lhs), Box::new(rhs));
        }
        Ok(lhs)
    }

    /// Right-hand side of `^`: allows a unary sign (`2^-1`) and is right
    /// associative.
    fn parse_pow_rhs(&mut self) -> Result<Node> {
        self.skip_ws();
        match self.peek() {
            Some(b'-') => {
                self.pos += 1;
                Ok(Node::Unary(UnOp::Neg, Box::new(self.parse_pow_rhs()?)))
            }
            Some(b'+') => {
                self.pos += 1;
                self.parse_pow_rhs()
            }
            _ => {
                let base = self.parse_atom()?;
                self.maybe_pow(base)
            }
        }
    }

    fn maybe_pow(&mut self, base: Node) -> Result<Node> {
        self.skip_ws();
        if self.peek() == Some(b'^') {
            self.pos += 1;
            let e = self.parse_pow_rhs()?;
            Ok(Node::Binary(BinOp::Pow, Box::new(base), Box::new(e)))
        } else {
            Ok(base)
        }
    }

    fn parse_unary(&mut self) -> Result<Node> {
        self.skip_ws();
        match self.peek() {
            Some(b'-') => {
                self.pos += 1;
                Ok(Node::Unary(UnOp::Neg, Box::new(self.parse_unary()?)))
            }
            Some(b'+') => {
                self.pos += 1;
                self.parse_unary()
            }
            Some(b'!') if self.peek2() != Some(b'=') => {
                self.pos += 1;
                Ok(Node::Unary(UnOp::Not, Box::new(self.parse_unary()?)))
            }
            _ => {
                let a = self.parse_atom()?;
                self.maybe_pow(a)
            }
        }
    }

    fn parse_atom(&mut self) -> Result<Node> {
        self.skip_ws();
        let start = self.pos;
        match self.peek() {
            None => Err(self.err("unexpected end of expression")),
            Some(b'(') => {
                self.pos += 1;
                let e = self.parse_expr(0)?;
                self.expect(b')')?;
                Ok(e)
            }
            Some(c) if c.is_ascii_digit() || c == b'.' => self.parse_number(start),
            Some(c) if c.is_ascii_alphabetic() || c == b'_' => {
                while matches!(self.peek(), Some(c) if c.is_ascii_alphanumeric() || c == b'_') {
                    self.pos += 1;
                }
                let name = &self.src[start..self.pos];
                self.skip_ws();
                if self.peek() == Some(b'(') {
                    let Some((f, arity)) = Func::lookup(name) else {
                        self.pos = start;
                        return Err(self.err(&format!("unknown function {name:?}")));
                    };
                    self.pos += 1;
                    let mut args = Vec::with_capacity(arity);
                    self.skip_ws();
                    if self.peek() != Some(b')') {
                        loop {
                            args.push(self.parse_expr(0)?);
                            self.skip_ws();
                            if self.peek() == Some(b',') {
                                self.pos += 1;
                            } else {
                                break;
                            }
                        }
                    }
                    self.expect(b')')?;
                    if args.len() != arity {
                        self.pos = start;
                        return Err(self.err(&format!(
                            "function {name:?} expects {arity} argument(s), got {}",
                            args.len()
                        )));
                    }
                    Ok(Node::Call(f, args))
                } else {
                    Ok(Node::Name(name.to_owned()))
                }
            }
            Some(_) => Err(self.err("unexpected character")),
        }
    }

    fn parse_number(&mut self, start: usize) -> Result<Node> {
        let b = self.src.as_bytes();
        let mut i = self.pos;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
        if i < b.len() && b[i] == b'.' {
            i += 1;
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
            }
        }
        if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
            let mut j = i + 1;
            if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
                j += 1;
            }
            if j < b.len() && b[j].is_ascii_digit() {
                while j < b.len() && b[j].is_ascii_digit() {
                    j += 1;
                }
                i = j;
            }
        }
        self.pos = i;
        self.src[start..i]
            .parse::<f64>()
            .map(Node::Num)
            .map_err(|_| {
                self.pos = start;
                self.err("invalid number")
            })
    }

    fn expect(&mut self, c: u8) -> Result<()> {
        self.skip_ws();
        if self.peek() == Some(c) {
            self.pos += 1;
            Ok(())
        } else {
            Err(self.err(&format!("expected '{}'", c as char)))
        }
    }
}

/// Parses `s` as a plain floating point number, accepting surrounding
/// whitespace. Returns `None` when `s` is not a bare number.
pub fn parse_number(s: &str) -> Option<f64> {
    let t = s.trim();
    let first = *t.as_bytes().first()?;
    if !(first.is_ascii_digit() || matches!(first, b'.' | b'-' | b'+')) {
        return None;
    }
    t.parse::<f64>().ok().filter(|v| v.is_finite())
}

/// Like [`parse_number`] but also accepts the textual NaN / infinity forms
/// written by C and C++ runtimes (`nan`, `-nan(ind)`, `inf`, `1.#QNAN`, ...).
pub fn parse_number_lenient(s: &str) -> Option<f64> {
    if let Some(v) = parse_number(s) {
        return Some(v);
    }
    let t = s.trim().to_ascii_lowercase();
    let (neg, body) = match t.strip_prefix('-') {
        Some(b) => (true, b),
        None => (false, t.strip_prefix('+').unwrap_or(&t)),
    };
    if body.starts_with("nan") || body.contains("#qnan") || body.contains("#ind") {
        return Some(f64::NAN);
    }
    if body == "inf" || body == "infinity" || body.contains("#inf") {
        return Some(if neg {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        });
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ev(s: &str) -> f64 {
        Expr::parse(s).unwrap().eval_const().unwrap()
    }

    #[test]
    fn arithmetic() {
        assert_eq!(ev("1+2*3"), 7.0);
        assert_eq!(ev("(1+2)*3"), 9.0);
        assert_eq!(ev("10-4-3"), 3.0);
        assert_eq!(ev("8/4/2"), 1.0);
        assert_eq!(ev("2^3^2"), 512.0);
        assert_eq!(ev("-2^2"), -4.0);
        assert_eq!(ev("2^-1"), 0.5);
        assert_eq!(ev("7%4"), 3.0);
        assert_eq!(ev(" 1.5e3 + .5 "), 1500.5);
        assert_eq!(ev("--3"), 3.0);
    }

    #[test]
    fn logic_and_functions() {
        assert_eq!(ev("1<2 & 3>=3"), 1.0);
        assert_eq!(ev("1=2 | 0"), 0.0);
        assert_eq!(ev("!0"), 1.0);
        assert_eq!(ev("1!=2"), 1.0);
        assert_eq!(ev("if(1>0, 5, 6)"), 5.0);
        assert_eq!(ev("max(3, min(10, 4))"), 4.0);
        assert!((ev("sin(pi/2)") - 1.0).abs() < 1e-15);
        assert_eq!(ev("atan2(0, 1)"), 0.0);
        assert_eq!(ev("sqrt(16)+abs(-1)"), 5.0);
    }

    #[test]
    fn variables() {
        let e = Expr::parse("a*x + b").unwrap();
        assert_eq!(
            e.names().into_iter().collect::<Vec<_>>(),
            vec!["a", "b", "x"]
        );
        let c = e
            .compile(&["x"], |n| match n {
                "a" => Some(2.0),
                "b" => Some(1.0),
                _ => None,
            })
            .unwrap();
        assert_eq!(c.eval(&[3.0]), 7.0);
        assert!(e.compile(&["x"], |_| None).is_err());
        assert_eq!(
            Expr::parse("2*pi")
                .unwrap()
                .compile(&[], |_| None)
                .unwrap()
                .as_constant(),
            Some(core::f64::consts::TAU)
        );
        // parameters shadow built-in constants
        assert_eq!(
            Expr::parse("e").unwrap().eval_with(|_| Some(3.0)).unwrap(),
            3.0
        );
    }

    #[test]
    fn errors() {
        for s in ["", "1+", "(1", "foo(1)", "sin(1,2)", "1 2", "#"] {
            assert!(Expr::parse(s).is_err(), "{s}");
        }
    }

    #[test]
    fn numbers() {
        assert_eq!(parse_number(" -1e-3 "), Some(-1e-3));
        assert_eq!(parse_number("a"), None);
        assert_eq!(parse_number("inf"), None);
        assert!(parse_number_lenient("-nan(ind)").unwrap().is_nan());
        assert!(parse_number_lenient("nan").unwrap().is_nan());
        assert_eq!(parse_number_lenient("-inf"), Some(f64::NEG_INFINITY));
        assert_eq!(parse_number_lenient("x"), None);
    }
}
