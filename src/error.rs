//! Error type.
//!
//! Every fallible operation of the crate returns [`Result`], whose error is
//! the single [`Error`] enum. Errors carry human readable messages; match on
//! the variant to distinguish causes.
//!
//! ```
//! use csx::{Error, ParameterSet};
//!
//! let params = ParameterSet::new();
//! match params.eval("width * 2") {
//!     Err(Error::UnknownVariable(name)) => assert_eq!(name, "width"),
//!     other => panic!("unexpected {other:?}"),
//! }
//! assert!(matches!(params.eval("1 +"), Err(Error::Parse { .. })));
//! ```

use core::fmt;

/// Errors produced by this crate.
///
/// `#[non_exhaustive]`: new variants may be added; match with a wildcard
/// arm.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// An expression could not be parsed.
    Parse {
        /// The expression source.
        expr: String,
        /// Byte offset of the problem.
        pos: usize,
        /// Description.
        msg: String,
    },
    /// An expression references a name that is neither a variable nor a
    /// parameter.
    UnknownVariable(String),
    /// A value is out of its valid domain (negative radius, ...).
    Invalid(String),
    /// Malformed or unsupported file content.
    Format(String),
    /// I/O failure (message of the underlying [`std::io::Error`], which is
    /// converted by `?`).
    Io(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Parse { expr, pos, msg } => {
                write!(f, "cannot parse expression {expr:?} at offset {pos}: {msg}")
            }
            Error::UnknownVariable(v) => write!(f, "unknown variable or parameter {v:?}"),
            Error::Invalid(m) => write!(f, "invalid value: {m}"),
            Error::Format(m) => write!(f, "format error: {m}"),
            Error::Io(m) => write!(f, "I/O error: {m}"),
        }
    }
}

impl std::error::Error for Error {}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e.to_string())
    }
}

/// `Result` with [`Error`] as the default error type.
pub type Result<T, E = Error> = core::result::Result<T, E>;
