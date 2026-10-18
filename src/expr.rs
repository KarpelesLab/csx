//! Mathematical expressions.
//!
//! Values in CSXCAD models (coordinates, material parameters, weighting
//! functions, ...) may be written as expressions over named parameters and,
//! for weighting functions, over coordinate variables. The accepted syntax
//! follows the "function parser" dialect used by CSXCAD files.
//!
//! [`Expr::parse`] produces an [`Expr`] that can be evaluated with any name
//! lookup; [`Expr::compile`] binds names once and yields a [`CompiledExpr`]
//! for fast repeated evaluation (e.g. per mesh cell).
//!
//! ```
//! use csx::expr::Expr;
//!
//! let e = Expr::parse("a * sin(x) + 1")?;
//! assert_eq!(e.names().into_iter().collect::<Vec<_>>(), ["a", "x"]);
//!
//! // Ad-hoc evaluation with a closure resolving names.
//! let v = e.eval_with(|n| match n { "a" => Some(2.0), "x" => Some(0.0), _ => None })?;
//! assert_eq!(v, 1.0);
//!
//! // Compile once: `x` becomes argument 0, `a` is folded to a constant.
//! let f = e.compile(&["x"], |n| (n == "a").then_some(2.0))?;
//! assert!((f.eval(&[core::f64::consts::FRAC_PI_2]) - 3.0).abs() < 1e-12);
//! # Ok::<(), csx::Error>(())
//! ```
//!
//! # Syntax
//!
//! **Numbers**: `1`, `2.5`, `.5`, `1e-3`, `6.02E23`.
//!
//! **Operators**, from lowest to highest precedence:
//!
//! | precedence | operators | notes |
//! |---|---|---|
//! | 1 | `a \| b` | logical or (short-circuit) |
//! | 2 | `a & b` | logical and (short-circuit) |
//! | 3 | `=` / `==`, `!=`, `<`, `<=`, `>`, `>=` | comparisons, result `1` or `0` |
//! | 4 | `+`, `-` | left associative |
//! | 5 | `*`, `/`, `%` | left associative; `%` is the floating remainder |
//! | 6 | unary `-`, `+`, `!` | `!x` is logical not |
//! | 7 | `^` | power, right associative, binds tighter than unary minus: `-2^2 = -4`, `2^-1 = 0.5` |
//!
//! Parentheses group as usual; there is no nesting limit.
//!
//! **Functions** (`name(args)`, fixed arity):
//!
//! | functions | meaning |
//! |---|---|
//! | `abs ceil floor trunc` | rounding / magnitude |
//! | `int(x)` | round to nearest |
//! | `sqrt cbrt exp exp2 log log2 log10` | `log` is the natural logarithm |
//! | `pow(x, y)`, `hypot(x, y)` | `xʸ`, `√(x²+y²)` |
//! | `min(a, b)`, `max(a, b)` | |
//! | `sin cos tan cot sec csc` | trigonometry, radians |
//! | `asin acos atan atan2(y, x)` | inverse trigonometry |
//! | `sinh cosh tanh asinh acosh atanh` | hyperbolic |
//! | `j0 j1 jn(n, x)` | Bessel functions of the first kind |
//! | `y0 y1 yn(n, x)` | Bessel functions of the second kind (NaN for `x < 0`) |
//! | `if(c, a, b)` | `a` if `c` is true else `b`; only the taken branch is evaluated |
//! | `not(x)` | logical not |
//!
//! **Constants**: `pi` and `e`. They are only used when the name is not
//! resolved otherwise, so a parameter or variable named `e` shadows Euler's
//! number.
//!
//! **Truth values**: logical operators return `1` or `0`. Any non-zero
//! value is *true*; `0` and NaN are *false*.
//!
//! **NaN**: arithmetic follows IEEE 754 (e.g. `0/0` and `sqrt(-1)` are NaN).
//! Contexts that need a definite number reject NaN results (see
//! [`crate::Scalar::update`]).
//!
//! # Implementation
//!
//! Expressions compile to flat postfix programs evaluated with an explicit
//! value stack, so there is no limit on nesting depth (memory use is linear
//! in the input length) and evaluating a [`CompiledExpr`] does not allocate
//! for typical expressions.

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
    // NaN (e.g. 0/0) is false, not "non-zero".
    v != 0.0 && !v.is_nan()
}

#[inline]
fn bool_val(b: bool) -> f64 {
    if b { 1.0 } else { 0.0 }
}

/// One instruction of the postfix program an expression compiles to.
///
/// Expressions are stored as flat programs evaluated with a value stack, so
/// parsing, evaluation and dropping never recurse and nesting depth is only
/// limited by memory.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Op {
    /// Push a constant.
    Num(f64),
    /// Push the value of a named variable (index into the name table).
    Name(u32),
    /// Push a positional variable of a compiled expression.
    Slot(u32),
    /// Unary operator on the top of the stack.
    Unary(UnOp),
    /// Binary operator on the two topmost values.
    Bin(BinOp),
    /// Function call on the topmost `arity` values (never `if`).
    Call(Func),
    /// `&`: if the top is false, replace it by `0` and jump; else pop it.
    AndJump(u32),
    /// `|`: if the top is true, replace it by `1` and jump; else pop it.
    OrJump(u32),
    /// Replace the top by its truth value (`1` / `0`).
    ToBool,
    /// Pop the top and jump if it is false (`if` condition).
    JumpIfFalse(u32),
    /// Unconditional jump.
    Jump(u32),
}

impl Func {
    fn arity(self) -> usize {
        match self {
            Func::Atan2 | Func::Hypot | Func::Jn | Func::Max | Func::Min | Func::Pow | Func::Yn => {
                2
            }
            Func::If => 3,
            _ => 1,
        }
    }
}

/// A postfix program plus the stack depth it needs.
#[derive(Debug, Clone, PartialEq, Default)]
struct Program {
    code: Vec<Op>,
    max_stack: usize,
}

/// Values kept on the machine stack before falling back to the heap.
const INLINE_STACK: usize = 32;

impl Program {
    /// Computes the maximum stack depth. Code is structured (jumps only go
    /// forward and every branch leaves the same depth), so a single pass
    /// suffices.
    fn finish(code: Vec<Op>) -> Program {
        let mut at: Vec<Option<usize>> = vec![None; code.len() + 1];
        let (mut d, mut max) = (0usize, 0usize);
        for (i, op) in code.iter().enumerate() {
            if let Some(v) = at[i] {
                d = v;
            }
            match *op {
                Op::Num(_) | Op::Name(_) | Op::Slot(_) => d += 1,
                Op::Unary(_) | Op::ToBool => {}
                Op::Bin(_) => d -= 1,
                Op::Call(f) => d = d + 1 - f.arity(),
                Op::AndJump(t) | Op::OrJump(t) => {
                    at[t as usize] = Some(d);
                    d -= 1;
                }
                Op::JumpIfFalse(t) => {
                    d -= 1;
                    at[t as usize] = Some(d);
                }
                Op::Jump(t) => at[t as usize] = Some(d),
            }
            max = max.max(d);
        }
        Program {
            code,
            max_stack: max,
        }
    }

    /// Runs the program; `name` resolves [`Op::Name`] operands and aborts
    /// the evaluation by returning `None`.
    #[inline]
    fn run<F>(&self, slots: &[f64], name: F) -> Option<f64>
    where
        F: FnMut(u32) -> Option<f64>,
    {
        // Small programs get a small buffer: zeroing it is part of every
        // evaluation.
        if self.max_stack <= 8 {
            let mut buf = [0.0f64; 8];
            exec(&self.code, slots, &mut buf, name)
        } else if self.max_stack <= INLINE_STACK {
            let mut buf = [0.0f64; INLINE_STACK];
            exec(&self.code, slots, &mut buf, name)
        } else {
            let mut buf = vec![0.0f64; self.max_stack];
            exec(&self.code, slots, &mut buf, name)
        }
    }
}

#[inline]
fn exec<F>(code: &[Op], slots: &[f64], st: &mut [f64], mut name: F) -> Option<f64>
where
    F: FnMut(u32) -> Option<f64>,
{
    let mut sp = 0usize;
    let mut pc = 0usize;
    while let Some(op) = code.get(pc) {
        match *op {
            Op::Num(v) => {
                st[sp] = v;
                sp += 1;
            }
            Op::Slot(i) => {
                st[sp] = slots.get(i as usize).copied().unwrap_or(0.0);
                sp += 1;
            }
            Op::Name(i) => {
                st[sp] = name(i)?;
                sp += 1;
            }
            Op::Unary(UnOp::Neg) => st[sp - 1] = -st[sp - 1],
            Op::Unary(UnOp::Not) => st[sp - 1] = bool_val(!truthy(st[sp - 1])),
            Op::Bin(b) => {
                sp -= 1;
                st[sp - 1] = binop(b, st[sp - 1], st[sp]);
            }
            Op::Call(f) => {
                let n = f.arity();
                let r = f.apply(&st[sp - n..sp]);
                sp -= n - 1;
                st[sp - 1] = r;
            }
            Op::AndJump(t) => {
                if !truthy(st[sp - 1]) {
                    st[sp - 1] = 0.0;
                    pc = t as usize;
                    continue;
                }
                sp -= 1;
            }
            Op::OrJump(t) => {
                if truthy(st[sp - 1]) {
                    st[sp - 1] = 1.0;
                    pc = t as usize;
                    continue;
                }
                sp -= 1;
            }
            Op::ToBool => st[sp - 1] = bool_val(truthy(st[sp - 1])),
            Op::JumpIfFalse(t) => {
                sp -= 1;
                if !truthy(st[sp]) {
                    pc = t as usize;
                    continue;
                }
            }
            Op::Jump(t) => {
                pc = t as usize;
                continue;
            }
        }
        pc += 1;
    }
    Some(if sp > 0 { st[sp - 1] } else { f64::NAN })
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

/// A parsed expression; see the [module documentation](self) for the
/// syntax.
///
/// `Expr` implements [`FromStr`](core::str::FromStr) and displays as its
/// source text.
///
/// ```
/// use csx::expr::Expr;
///
/// let e: Expr = "if(x > 0, x, -x)".parse()?;
/// assert_eq!(e.eval_with(|_| Some(-3.0))?, 3.0);
/// assert_eq!(e.to_string(), "if(x > 0, x, -x)");
/// # Ok::<(), csx::Error>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Expr {
    src: String,
    prog: Program,
    names: Vec<String>,
}

impl Expr {
    /// Parses `src`. There is no limit on nesting depth: parsing and
    /// evaluation use explicit stacks, so memory use is linear in the input
    /// length. Names are not resolved at this point.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] with the byte offset of the first syntax error
    /// (unexpected character or end of input, unknown function, wrong
    /// argument count, unbalanced parentheses).
    ///
    /// ```
    /// use csx::{expr::Expr, Error};
    ///
    /// assert!(matches!(Expr::parse("sin(1, 2)"), Err(Error::Parse { .. })));
    /// assert!(matches!(Expr::parse("(1 + 2"), Err(Error::Parse { pos: 6, .. })));
    /// ```
    pub fn parse(src: &str) -> Result<Expr> {
        let mut names: Vec<String> = Vec::new();
        let code = Parser::new(src).parse(&mut |n: &str| {
            let i = match names.iter().position(|m| m == n) {
                Some(i) => i,
                None => {
                    names.push(n.to_owned());
                    names.len() - 1
                }
            };
            Ok(Op::Name(i as u32))
        })?;
        Ok(Expr {
            src: src.to_owned(),
            prog: Program::finish(code),
            names,
        })
    }

    /// Original source text.
    pub fn source(&self) -> &str {
        &self.src
    }

    /// Names referenced by the expression (excluding function names), in
    /// sorted order.
    pub fn names(&self) -> BTreeSet<&str> {
        self.names.iter().map(String::as_str).collect()
    }

    /// Evaluates the expression; names are resolved by `lookup`, falling
    /// back to the built-in constants. Names are only resolved when they are
    /// actually evaluated (e.g. in the taken branch of `if`).
    ///
    /// # Errors
    ///
    /// [`Error::UnknownVariable`] naming an evaluated name that neither
    /// `lookup` nor the built-in constants resolve.
    ///
    /// ```
    /// use csx::expr::Expr;
    ///
    /// let e = Expr::parse("if(1, 2, unknown)")?;
    /// assert_eq!(e.eval_with(|_| None)?, 2.0); // `unknown` is never evaluated
    /// assert!(Expr::parse("unknown")?.eval_with(|_| None).is_err());
    /// # Ok::<(), csx::Error>(())
    /// ```
    pub fn eval_with(&self, mut lookup: impl FnMut(&str) -> Option<f64>) -> Result<f64> {
        let mut missing = 0;
        self.prog
            .run(&[], |i| {
                let n = &self.names[i as usize];
                let v = lookup(n).or_else(|| builtin_constant(n));
                if v.is_none() {
                    missing = i;
                }
                v
            })
            .ok_or_else(|| Error::UnknownVariable(self.names[missing as usize].clone()))
    }

    /// Evaluates an expression without any variables (only the built-in
    /// constants are available).
    ///
    /// ```
    /// let v = csx::expr::Expr::parse("2^10 - max(3, 4)")?.eval_const()?;
    /// assert_eq!(v, 1020.0);
    /// # Ok::<(), csx::Error>(())
    /// ```
    ///
    /// # Errors
    ///
    /// [`Error::UnknownVariable`] if the expression evaluates any name other
    /// than `pi` or `e`.
    pub fn eval_const(&self) -> Result<f64> {
        self.eval_with(|_| None)
    }

    /// Compiles the expression for fast repeated evaluation. Names listed in
    /// `vars` become positional arguments of [`CompiledExpr::eval`]; all other
    /// names are resolved once through `lookup` (then built-in constants) and
    /// constant sub-expressions are folded. Every name must resolve, even in
    /// branches that are never taken.
    ///
    /// ```
    /// use csx::expr::Expr;
    ///
    /// let e = Expr::parse("k * (x + y)")?;
    /// let f = e.compile(&["x", "y"], |n| (n == "k").then_some(10.0))?;
    /// assert_eq!(f.eval(&[1.0, 2.0]), 30.0);
    /// assert_eq!(f.as_constant(), None);
    ///
    /// let c = Expr::parse("k * 2")?.compile(&[], |_| Some(4.0))?;
    /// assert_eq!(c.as_constant(), Some(8.0));
    /// # Ok::<(), csx::Error>(())
    /// ```
    ///
    /// # Errors
    ///
    /// [`Error::UnknownVariable`] for a name that is neither in `vars`, nor
    /// resolved by `lookup`, nor a built-in constant.
    pub fn compile(
        &self,
        vars: &[&str],
        mut lookup: impl FnMut(&str) -> Option<f64>,
    ) -> Result<CompiledExpr> {
        let code = Parser::new(&self.src).parse(&mut |n: &str| {
            if let Some(i) = vars.iter().position(|v| *v == n) {
                Ok(Op::Slot(i as u32))
            } else if let Some(v) = lookup(n).or_else(|| builtin_constant(n)) {
                Ok(Op::Num(v))
            } else {
                Err(Error::UnknownVariable(n.to_owned()))
            }
        })?;
        Ok(CompiledExpr {
            prog: Program::finish(code),
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
/// to constants; see [`Expr::compile`]. Evaluation is infallible: the result
/// is a plain `f64` (possibly NaN).
#[derive(Debug, Clone, PartialEq)]
pub struct CompiledExpr {
    prog: Program,
}

impl CompiledExpr {
    /// Evaluates with the given variable values (in the order passed to
    /// [`Expr::compile`]). Missing values are treated as `0`. Does not
    /// allocate unless the expression needs more than 32 stack slots.
    #[inline]
    pub fn eval(&self, vars: &[f64]) -> f64 {
        if let [Op::Num(v)] = self.prog.code[..] {
            return v;
        }
        self.prog.run(vars, |_| Some(f64::NAN)).unwrap_or(f64::NAN)
    }

    /// The constant value, if the expression does not depend on any variable.
    pub fn as_constant(&self) -> Option<f64> {
        match self.prog.code[..] {
            [Op::Num(v)] => Some(v),
            _ => None,
        }
    }
}

/// Pending operator on the parser's operator stack.
enum Pending {
    /// Binary operator; for `&` / `|` the index of the emitted jump and
    /// whether the left operand was a foldable constant.
    Bin {
        op: BinOp,
        prec: u8,
        jump: usize,
        fold_lhs: bool,
    },
    Unary(UnOp),
    Paren,
    /// Function call in progress.
    Call {
        f: Func,
        start: usize,
        argc: usize,
        /// `if`: index of the conditional jump / end jump, and (when the
        /// condition is a foldable constant) the folding barrier before it.
        jumps: [usize; 2],
        fold_cond: Option<usize>,
    },
}

/// Precedence of the prefix operators `-` and `!` (between `*` and `^`).
const UNARY_PREC: u8 = 6;

/// Code emitter with constant folding.
struct Emitter {
    code: Vec<Op>,
    /// No folding may consume instructions before this index (a jump lands
    /// there).
    barrier: usize,
}

impl Emitter {
    fn num_at(&self, i: usize) -> Option<f64> {
        if i < self.barrier {
            return None;
        }
        match self.code.get(i) {
            Some(Op::Num(v)) => Some(*v),
            _ => None,
        }
    }

    fn unary(&mut self, op: UnOp) {
        let n = self.code.len();
        if let Some(a) = n.checked_sub(1).and_then(|i| self.num_at(i)) {
            self.code[n - 1] = Op::Num(match op {
                UnOp::Neg => -a,
                UnOp::Not => bool_val(!truthy(a)),
            });
        } else {
            self.code.push(Op::Unary(op));
        }
    }

    fn bin(&mut self, op: BinOp) {
        let n = self.code.len();
        if n >= 2
            && let (Some(a), Some(b)) = (self.num_at(n - 2), self.num_at(n - 1))
        {
            self.code.truncate(n - 2);
            self.code.push(Op::Num(binop(op, a, b)));
        } else {
            self.code.push(Op::Bin(op));
        }
    }

    fn call(&mut self, f: Func) {
        let n = self.code.len();
        let k = f.arity();
        if n >= k {
            let mut args = [0.0; 3];
            if (0..k).all(|j| match self.num_at(n - k + j) {
                Some(v) => {
                    args[j] = v;
                    true
                }
                None => false,
            }) {
                self.code.truncate(n - k);
                self.code.push(Op::Num(f.apply(&args[..k])));
                return;
            }
        }
        self.code.push(Op::Call(f));
    }

    /// Points the jump at `at` to the current end.
    fn patch(&mut self, at: usize) {
        let here = self.code.len() as u32;
        match &mut self.code[at] {
            Op::AndJump(t) | Op::OrJump(t) | Op::JumpIfFalse(t) | Op::Jump(t) => *t = here,
            _ => unreachable!("not a jump"),
        }
        self.barrier = self.barrier.max(self.code.len());
    }

    /// Completes `lhs & rhs` / `lhs | rhs` whose jump is at `jump`.
    fn finish_logic(&mut self, op: BinOp, jump: usize, fold_lhs: bool) {
        let n = self.code.len();
        if fold_lhs
            && n == jump + 2
            && let (Some(Op::Num(a)), Some(b)) =
                (self.code.get(jump - 1).copied(), self.num_at(n - 1))
        {
            self.code.truncate(jump - 1);
            self.code.push(Op::Num(binop(op, a, b)));
            return;
        }
        self.code.push(Op::ToBool);
        self.patch(jump);
    }
}

struct Parser<'a> {
    src: &'a str,
    pos: usize,
}

impl<'a> Parser<'a> {
    fn new(src: &'a str) -> Self {
        Parser { src, pos: 0 }
    }

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

    /// Reads a binary operator (without consuming it) and its precedence.
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

    /// Emits the code for a popped operator.
    fn reduce(em: &mut Emitter, p: Pending) {
        match p {
            Pending::Unary(op) => em.unary(op),
            Pending::Bin {
                op: op @ (BinOp::And | BinOp::Or),
                jump,
                fold_lhs,
                ..
            } => em.finish_logic(op, jump, fold_lhs),
            Pending::Bin { op, .. } => em.bin(op),
            Pending::Paren | Pending::Call { .. } => unreachable!("not an operator"),
        }
    }

    /// Pops operators that bind at least as tightly as an incoming binary
    /// operator of precedence `prec` (`^` is right associative).
    fn reduce_for(em: &mut Emitter, ops: &mut Vec<Pending>, prec: u8, right_assoc: bool) {
        while let Some(top) = ops.last() {
            let top_prec = match top {
                Pending::Bin { prec, .. } => *prec,
                Pending::Unary(_) => UNARY_PREC,
                Pending::Paren | Pending::Call { .. } => break,
            };
            if top_prec > prec || (top_prec == prec && !right_assoc) {
                let p = ops.pop().expect("non-empty");
                Self::reduce(em, p);
            } else {
                break;
            }
        }
    }

    /// Pops all operators down to the nearest parenthesis / call marker.
    fn reduce_group(em: &mut Emitter, ops: &mut Vec<Pending>) {
        while matches!(ops.last(), Some(Pending::Bin { .. } | Pending::Unary(_))) {
            let p = ops.pop().expect("non-empty");
            Self::reduce(em, p);
        }
    }

    /// Parses the whole input into postfix code; `name` maps identifiers to
    /// an operand instruction.
    fn parse(mut self, name: &mut dyn FnMut(&str) -> Result<Op>) -> Result<Vec<Op>> {
        let mut em = Emitter {
            code: Vec::new(),
            barrier: 0,
        };
        let mut ops: Vec<Pending> = Vec::new();
        let mut expect_operand = true;
        loop {
            self.skip_ws();
            if expect_operand {
                let start = self.pos;
                match self.peek() {
                    None => return Err(self.err("unexpected end of expression")),
                    Some(b'(') => {
                        self.pos += 1;
                        ops.push(Pending::Paren);
                    }
                    Some(b'-') => {
                        self.pos += 1;
                        ops.push(Pending::Unary(UnOp::Neg));
                    }
                    Some(b'+') => self.pos += 1,
                    Some(b'!') if self.peek2() != Some(b'=') => {
                        self.pos += 1;
                        ops.push(Pending::Unary(UnOp::Not));
                    }
                    Some(c) if c.is_ascii_digit() || c == b'.' => {
                        let v = self.parse_number(start)?;
                        em.code.push(Op::Num(v));
                        expect_operand = false;
                    }
                    Some(c) if c.is_ascii_alphabetic() || c == b'_' => {
                        while matches!(self.peek(), Some(c) if c.is_ascii_alphanumeric() || c == b'_')
                        {
                            self.pos += 1;
                        }
                        let id = &self.src[start..self.pos];
                        self.skip_ws();
                        if self.peek() == Some(b'(') {
                            let Some((f, _)) = Func::lookup(id) else {
                                self.pos = start;
                                return Err(self.err(&format!("unknown function {id:?}")));
                            };
                            self.pos += 1;
                            ops.push(Pending::Call {
                                f,
                                start,
                                argc: 0,
                                jumps: [0; 2],
                                fold_cond: None,
                            });
                            self.skip_ws();
                            if self.peek() == Some(b')') {
                                self.pos = start;
                                let arity = f.arity();
                                return Err(self.err(&format!(
                                    "function {id:?} expects {arity} argument(s), got 0"
                                )));
                            }
                        } else {
                            em.code.push(name(id)?);
                            expect_operand = false;
                        }
                    }
                    Some(_) => return Err(self.err("unexpected character")),
                }
                continue;
            }
            // Operator position.
            if let Some((op, prec, len)) = self.peek_binop() {
                self.pos += len;
                Self::reduce_for(&mut em, &mut ops, prec, op == BinOp::Pow);
                let mut pending = Pending::Bin {
                    op,
                    prec,
                    jump: 0,
                    fold_lhs: false,
                };
                if let Pending::Bin { jump, fold_lhs, .. } = &mut pending
                    && matches!(op, BinOp::And | BinOp::Or)
                {
                    let n = em.code.len();
                    *fold_lhs = n >= 1 && em.num_at(n - 1).is_some();
                    *jump = n;
                    em.code.push(if op == BinOp::And {
                        Op::AndJump(0)
                    } else {
                        Op::OrJump(0)
                    });
                }
                ops.push(pending);
                expect_operand = true;
                continue;
            }
            match self.peek() {
                None => break,
                Some(b',') => {
                    Self::reduce_group(&mut em, &mut ops);
                    let Some(Pending::Call {
                        f,
                        argc,
                        jumps,
                        fold_cond,
                        ..
                    }) = ops.last_mut()
                    else {
                        return Err(self.err("unexpected trailing input"));
                    };
                    *argc += 1;
                    if *argc >= f.arity() {
                        return Err(self.err("expected ')'"));
                    }
                    if *f == Func::If {
                        let n = em.code.len();
                        if *argc == 1 {
                            *fold_cond =
                                (n >= 1 && em.num_at(n - 1).is_some()).then_some(em.barrier);
                            jumps[0] = n;
                            em.code.push(Op::JumpIfFalse(0));
                        } else {
                            jumps[1] = n;
                            em.code.push(Op::Jump(0));
                            let j = jumps[0];
                            em.patch(j);
                        }
                    }
                    self.pos += 1;
                    expect_operand = true;
                }
                Some(b')') => {
                    Self::reduce_group(&mut em, &mut ops);
                    match ops.pop() {
                        Some(Pending::Paren) => {}
                        Some(Pending::Call {
                            f,
                            start,
                            argc,
                            jumps,
                            fold_cond,
                        }) => {
                            let argc = argc + 1;
                            if argc != f.arity() {
                                self.pos = start;
                                return Err(self.err(&format!(
                                    "function expects {} argument(s), got {argc}",
                                    f.arity()
                                )));
                            }
                            if f == Func::If {
                                Self::finish_if(&mut em, jumps, fold_cond);
                            } else {
                                em.call(f);
                            }
                        }
                        _ => return Err(self.err("unexpected ')'")),
                    }
                    self.pos += 1;
                }
                Some(_) => return Err(self.err("unexpected trailing input")),
            }
        }
        Self::reduce_group(&mut em, &mut ops);
        if !ops.is_empty() {
            return Err(self.err("expected ')'"));
        }
        Ok(em.code)
    }

    /// Completes `if(c, a, b)`; folds it when all three parts are constant.
    fn finish_if(em: &mut Emitter, jumps: [usize; 2], fold_cond: Option<usize>) {
        let [jz, jmp] = jumps;
        let n = em.code.len();
        if let Some(barrier) = fold_cond
            && jmp == jz + 2
            && n == jmp + 2
            && let (Some(Op::Num(c)), Some(Op::Num(a)), Some(Op::Num(b))) = (
                em.code.get(jz - 1).copied(),
                em.code.get(jz + 1).copied(),
                em.code.get(jmp + 1).copied(),
            )
        {
            em.code.truncate(jz - 1);
            em.code.push(Op::Num(if truthy(c) { a } else { b }));
            // The patched conditional jump raised the barrier past the
            // folded code; restore the barrier from before the `if`.
            em.barrier = barrier;
            return;
        }
        em.patch(jmp);
    }

    fn parse_number(&mut self, start: usize) -> Result<f64> {
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
        self.src[start..i].parse::<f64>().map_err(|_| {
            self.pos = start;
            self.err("invalid number")
        })
    }
}

/// Parses `s` as a plain finite floating point number, accepting
/// surrounding whitespace. Returns `None` when `s` is not a bare number.
///
/// ```
/// use csx::expr::parse_number;
///
/// assert_eq!(parse_number(" -1.5e3 "), Some(-1500.0));
/// assert_eq!(parse_number("2*3"), None);
/// assert_eq!(parse_number("inf"), None);
/// ```
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
///
/// ```
/// use csx::expr::parse_number_lenient;
///
/// assert!(parse_number_lenient("-nan(ind)").unwrap().is_nan());
/// assert_eq!(parse_number_lenient("-inf"), Some(f64::NEG_INFINITY));
/// assert_eq!(parse_number_lenient("7"), Some(7.0));
/// ```
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
    use core::f64::consts::EULER_GAMMA;
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
    pub(crate) fn j0(x: f64) -> f64 {
        j01(x).0
    }

    /// `J1(x)`.
    pub(crate) fn j1(x: f64) -> f64 {
        j01(x).1
    }

    /// `Y0(x)` (NaN for negative arguments).
    pub(crate) fn y0(x: f64) -> f64 {
        y01(x).0
    }

    /// `Y1(x)` (NaN for negative arguments).
    pub(crate) fn y1(x: f64) -> f64 {
        y01(x).1
    }

    /// Converts a float order to an integer (rounded), rejecting absurd
    /// magnitudes.
    fn order(n: f64) -> Option<i64> {
        (n.is_finite() && n.abs() <= 1e6).then(|| n.round() as i64)
    }

    /// `J_n(x)` for integer `n`.
    pub(crate) fn jn(n: f64, x: f64) -> f64 {
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
    pub(crate) fn yn(n: f64, x: f64) -> f64 {
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
        assert_eq!(ev("1*2+3*4-5/5"), 13.0);
        assert_eq!(ev("2-3+4"), 3.0);
        let long = vec!["1"; 5000].join("+");
        assert_eq!(ev(&long), 5000.0);
        assert_eq!(ev("0 & (1/0) | 2"), 1.0);
        // NaN is not true
        assert_eq!(ev("if(0/0, 1, 2)"), 2.0);
        assert_eq!(ev("!(0/0)"), 1.0);
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
    fn precedence_corner_cases() {
        assert_eq!(ev("2^-1*3"), 1.5);
        assert_eq!(ev("-2^-2"), -0.25);
        assert_eq!(ev("2^-1^2"), 0.5);
        assert_eq!(ev("!1=0"), 1.0);
        assert_eq!(ev("1<2<3"), 1.0);
        assert_eq!(ev("3>2>1"), 0.0);
        assert_eq!(ev("-3%2"), -1.0);
        assert_eq!(ev("1 - -1"), 2.0);
        assert_eq!(ev("+-+1"), -1.0);
        assert_eq!(ev("1|0&0"), 1.0);
        assert_eq!(ev("(1|0)&0"), 0.0);
        assert_eq!(ev("2*if(1,3,4)"), 6.0);
        assert_eq!(ev("if(0,1,2)+3"), 5.0);
        assert_eq!(ev("if(1,if(0,5,6),7)"), 6.0);
        assert_eq!(ev("max(if(0,1,2), 1 & 3)"), 2.0);
    }

    /// Branches and short-circuit operators mixed with variables: folding
    /// must never merge code across a jump target.
    #[test]
    fn folding_respects_jumps() {
        let cases: &[(&str, f64, f64)] = &[
            ("if(x,1,2)+3", 0.0, 5.0),
            ("if(x,1,2)+3", 1.0, 4.0),
            ("3+if(x,1,2)", 0.0, 5.0),
            ("-if(x,1,2)", 1.0, -1.0),
            ("x&3", 2.0, 1.0),
            ("(x&3)+1", 0.0, 1.0),
            ("(x|0)*5", 0.0, 0.0),
            ("(x|0)*5", 4.0, 5.0),
            ("max(if(x,1,2),3)", 1.0, 3.0),
            ("2*if(x,3,4)", 0.0, 8.0),
            ("if(x,1+1,2*3)", 0.0, 6.0),
            ("if(1,x,2)", 7.0, 7.0),
            ("if(0,x,2)*x", 7.0, 14.0),
            ("1 & x | 0", 0.0, 0.0),
            ("0 & x | 2", 9.0, 1.0),
        ];
        for &(src, x, want) in cases {
            let e = Expr::parse(src).unwrap();
            assert_eq!(e.eval_with(|_| Some(x)).unwrap(), want, "{src} eval_with");
            let c = e.compile(&["x"], |_| None).unwrap();
            assert_eq!(c.eval(&[x]), want, "{src} compiled");
            let k = e.compile(&[], |_| Some(x)).unwrap();
            assert_eq!(k.as_constant(), Some(want), "{src} folded");
        }
    }

    /// Lazily evaluated names: only the taken branch is resolved.
    #[test]
    fn lazy_names() {
        let e = Expr::parse("if(1, 2, missing)").unwrap();
        assert_eq!(e.eval_const().unwrap(), 2.0);
        assert!(e.compile(&[], |_| None).is_err());
        assert_eq!(
            Expr::parse("0 & missing").unwrap().eval_const().unwrap(),
            0.0
        );
        assert!(Expr::parse("1 & missing").unwrap().eval_const().is_err());
    }

    /// Random expressions: `eval_with`, compiled and fully folded evaluation
    /// agree.
    #[test]
    fn differential_random() {
        let mut seed = 0x1234_5678_9abc_def0u64;
        let mut rnd = |n: u64| {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (seed >> 33) % n
        };
        let atoms = ["x", "y", "0", "1", "2.5", "pi"];
        let bins = ["+", "-", "*", "/", "^", "<", ">=", "=", "!=", "&", "|", "%"];
        for _ in 0..3000 {
            // Build a random expression with a small explicit stack.
            let mut parts: Vec<String> = Vec::new();
            for _ in 0..(2 + rnd(6)) {
                let a = atoms[rnd(atoms.len() as u64) as usize].to_string();
                let t = match rnd(6) {
                    0 => format!("-{a}"),
                    1 => format!("!{a}"),
                    2 => format!("sin({a})"),
                    _ => a,
                };
                parts.push(t);
                while parts.len() >= 2 && rnd(3) == 0 {
                    let b = parts.pop().unwrap();
                    let a = parts.pop().unwrap();
                    let op = bins[rnd(bins.len() as u64) as usize];
                    parts.push(match rnd(4) {
                        0 => format!("if({a},{b},{a}{op}{b})"),
                        1 => format!("max({a},{b})"),
                        _ => format!("({a}{op}{b})"),
                    });
                }
            }
            let src = parts.join("+");
            let e = Expr::parse(&src).unwrap();
            let (x, y) = (rnd(5) as f64 - 2.0, rnd(3) as f64 * 0.5);
            let look = |n: &str| match n {
                "x" => Some(x),
                "y" => Some(y),
                _ => None,
            };
            let a = e.eval_with(look).unwrap();
            let b = e.compile(&["x", "y"], |_| None).unwrap().eval(&[x, y]);
            let c = e.compile(&[], look).unwrap().as_constant().unwrap();
            let same = |p: f64, q: f64| p == q || (p.is_nan() && q.is_nan());
            assert!(same(a, b) && same(a, c), "{src}: {a} {b} {c}");
        }
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
