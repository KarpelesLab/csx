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
//!   sec sin sinh sqrt tan tanh trunc`, plus the Bessel functions `j0 j1
//!   jn y0 y1 yn` (`jn(n, x)`, `yn(n, x)` with integer order `n`)
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
    J0,
    J1,
    Jn,
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
    Y0,
    Y1,
    Yn,
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
            "j0" => (J0, 1),
            "j1" => (J1, 1),
            "jn" => (Jn, 2),
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
            "y0" => (Y0, 1),
            "y1" => (Y1, 1),
            "yn" => (Yn, 2),
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
            J0 => bessel::j0(a[0]),
            J1 => bessel::j1(a[0]),
            Jn => bessel::jn(a[0], a[1]),
            Y0 => bessel::y0(a[0]),
            Y1 => bessel::y1(a[0]),
            Yn => bessel::yn(a[0], a[1]),
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
        let mut p = Parser {
            src,
            pos: 0,
            rec: 0,
        };
        let (root, _) = p.parse_expr(0)?;
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

/// Maximum nesting depth of expressions (parser recursion and tree depth).
///
/// Parsing, evaluation and dropping of expression trees are recursive, so
/// untrusted input is bounded to keep stack use small (well below 1 MiB
/// even in debug builds). Each parenthesis level counts twice; each term of
/// a left-associative chain such as `a+b+c` counts once.
pub const MAX_DEPTH: usize = 256;

struct Parser<'a> {
    src: &'a str,
    pos: usize,
    rec: usize,
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

    /// Guards parser recursion (parentheses, unary chains, ...).
    fn enter(&mut self) -> Result<()> {
        self.rec += 1;
        if self.rec > MAX_DEPTH {
            return Err(self.err("expression nested too deeply"));
        }
        Ok(())
    }

    /// Builds a node of depth `depth`, rejecting trees deeper than
    /// [`MAX_DEPTH`] (evaluation and drop are recursive).
    fn node(&self, node: Node, depth: usize) -> Result<(Node, usize)> {
        if depth > MAX_DEPTH {
            return Err(self.err("expression nested too deeply"));
        }
        Ok((node, depth))
    }

    fn parse_expr(&mut self, min_bp: u8) -> Result<(Node, usize)> {
        self.enter()?;
        let (mut lhs, mut d) = self.parse_unary()?;
        while let Some((op, bp, len)) = self.peek_binop() {
            if bp < min_bp {
                break;
            }
            self.pos += len;
            // `^` is right associative; everything else left associative.
            let (rhs, rd) = if op == BinOp::Pow {
                self.parse_pow_rhs()?
            } else {
                self.parse_expr(bp + 1)?
            };
            (lhs, d) = self.node(
                Node::Binary(op, Box::new(lhs), Box::new(rhs)),
                d.max(rd) + 1,
            )?;
        }
        self.rec -= 1;
        Ok((lhs, d))
    }

    /// Right-hand side of `^`: allows a unary sign (`2^-1`) and is right
    /// associative.
    fn parse_pow_rhs(&mut self) -> Result<(Node, usize)> {
        self.enter()?;
        self.skip_ws();
        let r = match self.peek() {
            Some(b'-') => {
                self.pos += 1;
                let (n, d) = self.parse_pow_rhs()?;
                self.node(Node::Unary(UnOp::Neg, Box::new(n)), d + 1)
            }
            Some(b'+') => {
                self.pos += 1;
                self.parse_pow_rhs()
            }
            _ => {
                let base = self.parse_atom()?;
                self.maybe_pow(base)
            }
        };
        self.rec -= 1;
        r
    }

    fn maybe_pow(&mut self, base: (Node, usize)) -> Result<(Node, usize)> {
        self.skip_ws();
        if self.peek() == Some(b'^') {
            self.pos += 1;
            let (e, d) = self.parse_pow_rhs()?;
            self.node(
                Node::Binary(BinOp::Pow, Box::new(base.0), Box::new(e)),
                base.1.max(d) + 1,
            )
        } else {
            Ok(base)
        }
    }

    fn parse_unary(&mut self) -> Result<(Node, usize)> {
        self.enter()?;
        self.skip_ws();
        let r = match self.peek() {
            Some(b'-') => {
                self.pos += 1;
                let (n, d) = self.parse_unary()?;
                self.node(Node::Unary(UnOp::Neg, Box::new(n)), d + 1)
            }
            Some(b'+') => {
                self.pos += 1;
                self.parse_unary()
            }
            Some(b'!') if self.peek2() != Some(b'=') => {
                self.pos += 1;
                let (n, d) = self.parse_unary()?;
                self.node(Node::Unary(UnOp::Not, Box::new(n)), d + 1)
            }
            _ => {
                let a = self.parse_atom()?;
                self.maybe_pow(a)
            }
        };
        self.rec -= 1;
        r
    }

    fn parse_atom(&mut self) -> Result<(Node, usize)> {
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
            Some(c) if c.is_ascii_digit() || c == b'.' => Ok((self.parse_number(start)?, 1)),
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
                    let mut depth = 0;
                    self.skip_ws();
                    if self.peek() != Some(b')') {
                        loop {
                            let (a, d) = self.parse_expr(0)?;
                            args.push(a);
                            depth = depth.max(d);
                            self.skip_ws();
                            if self.peek() == Some(b',') && args.len() < arity {
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
                    self.node(Node::Call(f, args), depth + 1)
                } else {
                    Ok((Node::Name(name.to_owned()), 1))
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

/// Bessel functions of the first and second kind for integer orders.
///
/// Small arguments use the ascending power series; large arguments use the
/// Hankel asymptotic expansion. Higher orders use recurrences (backward,
/// normalised recurrence for `J_n` when `n > |x|`, forward otherwise).
mod bessel {
    use core::f64::consts::{FRAC_PI_4, PI};

    /// Euler–Mascheroni constant.
    const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;
    /// Switch-over between power series and asymptotic expansion.
    const SERIES_LIMIT: f64 = 12.0;

    /// Hankel asymptotic expansion for order `nu` (0 or 1), `x > 0`:
    /// returns `(J_nu(x), Y_nu(x))`.
    fn asymptotic(nu: f64, x: f64) -> (f64, f64) {
        let mu = 4.0 * nu * nu;
        let (mut p, mut q) = (1.0, 0.0);
        let mut term = 1.0f64;
        let mut prev = f64::INFINITY;
        for k in 1..60 {
            let odd = (2 * k - 1) as f64;
            term *= (mu - odd * odd) / (k as f64 * 8.0 * x);
            if term.abs() >= prev || term.abs() < 1e-17 {
                break;
            }
            prev = term.abs();
            // k odd → Q with sign (-1)^((k-1)/2); k even → P with (-1)^(k/2)
            let sign = if (k / 2) % 2 == 0 { 1.0 } else { -1.0 };
            if k % 2 == 1 {
                q += sign * term;
            } else {
                p += sign * term;
            }
        }
        let chi = x - (nu / 2.0) * PI - FRAC_PI_4;
        let (s, c) = chi.sin_cos();
        let amp = (2.0 / (PI * x)).sqrt();
        (amp * (p * c - q * s), amp * (p * s + q * c))
    }

    /// Power series: `(J0, J1, Y0, Y1)` for `0 < x < SERIES_LIMIT`.
    fn series(x: f64) -> (f64, f64, f64, f64) {
        let h = x / 2.0;
        let h2 = -h * h;
        // J0 = Σ (-h²)^k/(k!)²; Y0 helper Σ (-1)^(k+1) H_k (h²)^k/(k!)²
        let (mut j0, mut y0s) = (1.0, 0.0);
        // J1 = h Σ (-h²)^k/(k!(k+1)!); Y1 helper Σ (-1)^k (ψ(k+1)+ψ(k+2)) h^(2k+1)/(k!(k+1)!)
        let (mut j1, mut y1s) = (h, h * (1.0 - 2.0 * EULER_GAMMA));
        let mut t0 = 1.0; // (-h²)^k/(k!)²
        let mut t1 = h; // h (-h²)^k/(k!(k+1)!)
        let mut hk = 0.0; // harmonic number H_k
        for k in 1..200 {
            let kf = k as f64;
            t0 *= h2 / (kf * kf);
            t1 *= h2 / (kf * (kf + 1.0));
            hk += 1.0 / kf;
            j0 += t0;
            y0s -= hk * t0;
            j1 += t1;
            // ψ(k+1) + ψ(k+2) = 2 H_k + 1/(k+1) - 2γ
            y1s += t1 * (2.0 * hk + 1.0 / (kf + 1.0) - 2.0 * EULER_GAMMA);
            if t0.abs() < 1e-17 * j0.abs().max(1e-300) && t1.abs() < 1e-17 {
                break;
            }
        }
        let ln = h.ln();
        let y0 = 2.0 / PI * ((ln + EULER_GAMMA) * j0 + y0s);
        let y1 = 2.0 / PI * ln * j1 - 2.0 / (PI * x) - y1s / PI;
        (j0, j1, y0, y1)
    }

    fn j01(x: f64) -> (f64, f64) {
        let ax = x.abs();
        let (j0, j1) = if ax == 0.0 {
            (1.0, 0.0)
        } else if ax < SERIES_LIMIT {
            let (a, b, _, _) = series(ax);
            (a, b)
        } else if ax.is_finite() {
            (asymptotic(0.0, ax).0, asymptotic(1.0, ax).0)
        } else if ax.is_infinite() {
            (0.0, 0.0)
        } else {
            (f64::NAN, f64::NAN)
        };
        (j0, if x < 0.0 { -j1 } else { j1 })
    }

    fn y01(x: f64) -> (f64, f64) {
        if x.is_nan() || x < 0.0 {
            (f64::NAN, f64::NAN)
        } else if x == 0.0 {
            (f64::NEG_INFINITY, f64::NEG_INFINITY)
        } else if x < SERIES_LIMIT {
            let (_, _, a, b) = series(x);
            (a, b)
        } else if x.is_infinite() {
            (0.0, 0.0)
        } else {
            (asymptotic(0.0, x).1, asymptotic(1.0, x).1)
        }
    }

    /// `J0(x)`.
    pub fn j0(x: f64) -> f64 {
        j01(x).0
    }

    /// `J1(x)`.
    pub fn j1(x: f64) -> f64 {
        j01(x).1
    }

    /// `Y0(x)` (NaN for negative arguments).
    pub fn y0(x: f64) -> f64 {
        y01(x).0
    }

    /// `Y1(x)` (NaN for negative arguments).
    pub fn y1(x: f64) -> f64 {
        y01(x).1
    }

    /// Converts a float order to an integer (rounded), rejecting absurd
    /// magnitudes.
    fn order(n: f64) -> Option<i64> {
        (n.is_finite() && n.abs() <= 1e6).then(|| n.round() as i64)
    }

    /// `J_n(x)` for integer `n`.
    pub fn jn(n: f64, x: f64) -> f64 {
        let Some(n) = order(n) else { return f64::NAN };
        if x.is_nan() {
            return f64::NAN;
        }
        let sign_n = if n < 0 && n % 2 != 0 { -1.0 } else { 1.0 };
        let m = n.unsigned_abs();
        let sign_x = if x < 0.0 && m % 2 == 1 { -1.0 } else { 1.0 };
        let ax = x.abs();
        let v = match m {
            0 => j0(ax),
            1 => j1(ax),
            _ if ax == 0.0 || ax.is_infinite() => 0.0,
            _ if ax > m as f64 => {
                // Forward recurrence is stable while k < x.
                let (mut a, mut b) = (j0(ax), j1(ax));
                for k in 1..m {
                    let c = 2.0 * k as f64 / ax * b - a;
                    a = b;
                    b = c;
                }
                b
            }
            _ => {
                // Miller's backward recurrence, normalised with
                // J0 + 2 Σ J_2k = 1.
                let start = 2 * ((m + 16 + (40.0 * m as f64).sqrt() as u64 + ax as u64) / 2);
                let (mut next, mut cur) = (0.0f64, 1e-300f64);
                let (mut sum, mut result) = (0.0f64, 0.0f64);
                for k in (1..=start).rev() {
                    let prev = 2.0 * k as f64 / ax * cur - next;
                    next = cur;
                    cur = prev; // cur = J_{k-1} (unnormalised)
                    if k - 1 == m {
                        result = cur;
                    }
                    if (k - 1) % 2 == 0 && k - 1 > 0 {
                        sum += cur;
                    }
                    if cur.abs() > 1e250 {
                        cur *= 1e-250;
                        next *= 1e-250;
                        sum *= 1e-250;
                        result *= 1e-250;
                    }
                }
                let norm = cur + 2.0 * sum; // cur = J0
                result / norm
            }
        };
        sign_n * sign_x * v
    }

    /// `Y_n(x)` for integer `n` (NaN for negative `x`).
    pub fn yn(n: f64, x: f64) -> f64 {
        let Some(n) = order(n) else { return f64::NAN };
        let sign_n = if n < 0 && n % 2 != 0 { -1.0 } else { 1.0 };
        let m = n.unsigned_abs();
        let v = match m {
            0 => y0(x),
            1 => y1(x),
            _ if x.is_nan() || x < 0.0 => f64::NAN,
            _ if x == 0.0 => f64::NEG_INFINITY,
            _ if x.is_infinite() => 0.0,
            _ => {
                let (mut a, mut b) = (y0(x), y1(x));
                for k in 1..m {
                    let c = 2.0 * k as f64 / x * b - a;
                    a = b;
                    b = c;
                    if b.is_infinite() {
                        break;
                    }
                }
                b
            }
        };
        sign_n * v
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        fn close(a: f64, b: f64, tol: f64) -> bool {
            (a - b).abs() <= tol * b.abs().max(1.0)
        }

        #[test]
        fn reference_values() {
            // Reference values (Abramowitz & Stegun tables / standard libm).
            let cases: [(f64, f64, f64, f64, f64); 8] = [
                // x, J0, J1, Y0, Y1
                (
                    0.5,
                    0.938_469_807_240_813,
                    0.242_268_457_674_874,
                    -0.444_518_733_506_707,
                    -1.471_472_392_670_243,
                ),
                (
                    1.0,
                    0.765_197_686_557_967,
                    0.440_050_585_744_934,
                    0.088_256_964_215_677,
                    -0.781_212_821_300_289,
                ),
                (
                    2.5,
                    -0.048_383_776_468_198,
                    0.497_094_102_464_274,
                    0.498_070_359_615_232,
                    0.145_918_137_966_786,
                ),
                (
                    5.0,
                    -0.177_596_771_314_338,
                    -0.327_579_137_591_465,
                    -0.308_517_625_249_034,
                    0.147_863_143_391_227,
                ),
                (
                    10.0,
                    -0.245_935_764_451_348,
                    0.043_472_746_168_862,
                    0.055_671_167_283_599,
                    0.249_015_424_206_954,
                ),
                (
                    11.9,
                    0.025_049_441_699_590,
                    -0.228_983_249_661_924,
                    -0.229_833_213_943_375,
                    -0.034_711_498_334_030,
                ),
                (
                    12.1,
                    0.069_666_773_606_808,
                    -0.215_748_973_376_925,
                    -0.218_438_380_550_925,
                    -0.078_736_931_451_396,
                ),
                (
                    50.0,
                    0.055_812_327_669_252,
                    -0.097_511_828_125_175,
                    -0.098_064_995_470_077,
                    -0.056_795_668_562_015,
                ),
            ];
            for (x, a, b, c, d) in cases {
                assert!(close(j0(x), a, 1e-9), "j0({x}) = {} want {a}", j0(x));
                assert!(close(j1(x), b, 1e-9), "j1({x}) = {} want {b}", j1(x));
                assert!(close(y0(x), c, 1e-9), "y0({x}) = {} want {c}", y0(x));
                assert!(close(y1(x), d, 1e-9), "y1({x}) = {} want {d}", y1(x));
            }
            assert_eq!(j0(0.0), 1.0);
            assert_eq!(j1(0.0), 0.0);
            assert!(close(j1(-1.0), -0.440_050_585_744_934, 1e-12));
            assert!(y0(-1.0).is_nan());
            assert_eq!(y0(0.0), f64::NEG_INFINITY);
        }

        #[test]
        fn integer_orders() {
            // J_2(1), J_5(1), J_3(10), J_20(10), Y_2(1), Y_3(10)
            assert!(close(jn(2.0, 1.0), 0.114_903_484_931_900, 1e-9));
            assert!(close(jn(5.0, 1.0), 2.497_577_302_112_344e-4, 1e-9));
            assert!(close(jn(3.0, 10.0), 0.058_379_379_305_187, 1e-9));
            assert!(close(jn(20.0, 10.0), 1.151_336_924_781_34e-5, 1e-8));
            assert!(close(yn(2.0, 1.0), -1.650_682_606_816_254, 1e-9));
            assert!(close(yn(3.0, 10.0), -0.251_362_657_183_837, 1e-9));
            // symmetries
            assert!(close(jn(-3.0, 2.0), -jn(3.0, 2.0), 1e-14));
            assert!(close(jn(3.0, -2.0), -jn(3.0, 2.0), 1e-14));
            assert!(close(yn(-3.0, 2.0), -yn(3.0, 2.0), 1e-14));
            assert_eq!(jn(4.0, 0.0), 0.0);
            assert!(jn(f64::NAN, 1.0).is_nan());
            assert!(jn(1e300, 1.0).is_nan());
            // Recurrence identity J_{n-1} + J_{n+1} = 2n/x J_n
            for &x in &[0.3, 3.7, 17.0, 40.0] {
                for n in 1..15 {
                    let n = n as f64;
                    let l = jn(n - 1.0, x) + jn(n + 1.0, x);
                    let r = 2.0 * n / x * jn(n, x);
                    assert!((l - r).abs() < 1e-9, "n={n} x={x}: {l} vs {r}");
                }
            }
        }
    }
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
