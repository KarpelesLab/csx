//! File formats.
//!
//! | Module | Purpose |
//! |---|---|
//! | [`xml`] (feature `xml`) | read / write CSXCAD and openEMS XML files |
//! | [`dom`] | the minimal XML document model used by [`xml`]; also holds unknown elements preserved on round-trips |
//! | [`stl`], [`ply`] | mesh readers for polyhedra (`PolyhedronReader`) |
//! | [`export`] | STL / PLY / VTK writers for tessellated geometry |
//!
//! All parsers work on in-memory data, never panic on malformed input and
//! bound their allocations by the input size.

pub mod dom;
pub mod export;
pub mod ply;
pub mod stl;
#[cfg(feature = "xml")]
pub mod xml;

#[cfg(feature = "xml")]
pub use xml::OpenEmsDocument;
