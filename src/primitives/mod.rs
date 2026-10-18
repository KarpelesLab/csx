//! Geometric primitives.
//!
//! A [`Primitive`] combines a [`Shape`] with a priority, an optional
//! coordinate system and an optional [`Transform`]. Primitives are owned by
//! properties (see [`crate::properties::Property`]); when several primitives
//! overlap, the one with the highest priority determines the property at a
//! point.
//!
//! Inside tests take a point together with the coordinate system it is
//! expressed in. The point is converted to Cartesian coordinates, mapped
//! through the inverse transformation into the primitive's local frame and
//! tested against the shape.
//!
//! # Coordinate systems
//!
//! The semantics below match CSXCAD 0.7 (verified black-box with random
//! models and point queries):
//!
//! * A [`BoxShape`] in cylindrical coordinates is an annular sector. Queried
//!   with cylindrical coordinates it compares them as given: a negative `ρ`
//!   is outside and `α` is taken modulo 2π.
//! * A [`MultiBox`] ignores its own coordinate system: its corners are in the
//!   structure's system and are compared component-wise, without angle
//!   wrapping.
//! * Point-defined shapes (spheres, cylinders, their shells, curves, wires,
//!   points) convert their defining points from `(ρ, α, z)` when their
//!   coordinate system is cylindrical; the geometry itself is Euclidean
//!   (straight wire segments, round spheres). See [`Shape::to_cartesian`].
//! * Polygon based shapes and polyhedra are always Cartesian.
//! * Shapes other than boxes see query points through the Cartesian
//!   conversion, so `(-ρ, α)` is the same point as `(ρ, α + π)`.
//!
//! Known divergence: CSXCAD 0.7 clips cylindrical-coordinate wires to a
//! bounding box that it computes inconsistently (it can even be inverted),
//! so parts of such wires are invisible to CSXCAD; this crate uses the true
//! geometry.

mod basic;
pub(crate) mod bvh;
mod curve;
mod geom;
mod poly;
mod polyhedron;
mod user;

pub use basic::{BoxShape, Cylinder, CylindricalShell, MultiBox, Sphere, SphericalShell};
pub use curve::{Curve, PointShape, Wire};
pub use poly::{LinPoly, PolyVertices, Polygon, RotPoly, plane_axes};
pub use polyhedron::Polyhedron;
pub use user::{USER_VARS, UserDefined};

pub(crate) use geom::cartesian_bbox_to_cylindrical;

use std::borrow::Cow;

use crate::coords::CoordSystem;
use crate::error::Result;
use crate::math::{BoundingBox, Vec3};
use crate::param::{ParameterSet, ScalarPoint};
use crate::transform::Transform;

/// Mesh file formats understood by [`PolyhedronReader`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MeshFileType {
    /// STL (ASCII or binary); file type id `1`.
    #[default]
    Stl,
    /// PLY (ASCII or binary); file type id `2`.
    Ply,
}

impl MeshFileType {
    /// Numeric id as used in files.
    pub fn id(self) -> i32 {
        match self {
            MeshFileType::Stl => 1,
            MeshFileType::Ply => 2,
        }
    }

    /// Guesses the type from a file name extension.
    pub fn from_path(path: &str) -> Option<Self> {
        path.rsplit('.').next()?.parse().ok()
    }
}

impl core::str::FromStr for MeshFileType {
    type Err = crate::Error;

    /// Parses a numeric id (`1`, `2`) or a type name (`STL`, `PLY`, case
    /// insensitive).
    fn from_str(s: &str) -> Result<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "1" | "stl" => Ok(MeshFileType::Stl),
            "2" | "ply" => Ok(MeshFileType::Ply),
            _ => Err(crate::Error::Format(format!(
                "unknown mesh file type {s:?}"
            ))),
        }
    }
}

impl core::fmt::Display for MeshFileType {
    /// The type name used in files (`STL`, `PLY`).
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            MeshFileType::Stl => "STL",
            MeshFileType::Ply => "PLY",
        })
    }
}

/// A polyhedron loaded from an external mesh file.
///
/// The loaded mesh is a cache of the file contents: changing the file name
/// or type through the setters discards it, and
/// [`crate::ContinuousStructure::update`] (or [`PolyhedronReader::load`])
/// loads it again. Equality compares the definition only.
#[derive(Debug, Clone, Default)]
pub struct PolyhedronReader {
    file_name: String,
    file_type: MeshFileType,
    mesh: Option<Polyhedron>,
}

impl PartialEq for PolyhedronReader {
    fn eq(&self, o: &Self) -> bool {
        self.file_name == o.file_name && self.file_type == o.file_type
    }
}

impl PolyhedronReader {
    /// Reader for `file_name`; the type is guessed from the extension.
    pub fn new(file_name: impl Into<String>) -> Self {
        let file_name = file_name.into();
        let file_type = MeshFileType::from_path(&file_name).unwrap_or_default();
        PolyhedronReader {
            file_name,
            file_type,
            mesh: None,
        }
    }

    /// Builder: overrides the file type.
    #[must_use]
    pub fn with_file_type(mut self, t: MeshFileType) -> Self {
        self.set_file_type(t);
        self
    }

    /// File name as written in the model.
    pub fn file_name(&self) -> &str {
        &self.file_name
    }

    /// File format.
    pub fn file_type(&self) -> MeshFileType {
        self.file_type
    }

    /// Changes the file name, discarding any loaded mesh.
    pub fn set_file_name(&mut self, name: impl Into<String>) {
        self.file_name = name.into();
        self.mesh = None;
    }

    /// Changes the file type, discarding any loaded mesh.
    pub fn set_file_type(&mut self, t: MeshFileType) {
        if t != self.file_type {
            self.file_type = t;
            self.mesh = None;
        }
    }

    /// The loaded mesh, if any.
    pub fn mesh(&self) -> Option<&Polyhedron> {
        self.mesh.as_ref()
    }

    /// True once the mesh has been loaded.
    pub fn is_loaded(&self) -> bool {
        self.mesh.is_some()
    }

    /// Loads the mesh from in-memory file contents (format per
    /// [`PolyhedronReader::file_type`]).
    pub fn load_bytes(&mut self, data: &[u8]) -> Result<()> {
        self.mesh = Some(match self.file_type {
            MeshFileType::Stl => crate::io::stl::parse(data)?,
            MeshFileType::Ply => crate::io::ply::parse(data)?,
        });
        Ok(())
    }

    /// Loads the mesh from disk. Relative paths are resolved against
    /// `base_dir` when given.
    pub fn load(&mut self, base_dir: Option<&std::path::Path>) -> Result<()> {
        let mut path = std::path::PathBuf::from(&self.file_name);
        if path.is_relative()
            && let Some(b) = base_dir
        {
            path = b.join(path);
        }
        let data = std::fs::read(&path)?;
        self.load_bytes(&data)
    }
}

/// The geometric definition of a primitive.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum Shape {
    /// Axis aligned box.
    Box(BoxShape),
    /// Union of boxes.
    MultiBox(MultiBox),
    /// Sphere.
    Sphere(Sphere),
    /// Spherical shell.
    SphericalShell(SphericalShell),
    /// Cylinder.
    Cylinder(Cylinder),
    /// Cylindrical shell.
    CylindricalShell(CylindricalShell),
    /// Planar polygon.
    Polygon(Polygon),
    /// Extruded polygon.
    LinPoly(LinPoly),
    /// Revolved polygon.
    RotPoly(RotPoly),
    /// Polyline.
    Curve(Curve),
    /// Polyline with radius.
    Wire(Wire),
    /// Single point.
    Point(PointShape),
    /// Polyhedron given inline.
    Polyhedron(Polyhedron),
    /// Polyhedron read from a file.
    PolyhedronReader(PolyhedronReader),
    /// Expression-defined region.
    UserDefined(UserDefined),
}

impl Shape {
    /// True for shapes defined by points (centres, axis end points, polyline
    /// vertices) that are converted from `(ρ, α, z)` when the primitive's
    /// coordinate system is cylindrical: spheres, cylinders, their shells,
    /// curves, wires and points.
    pub fn is_point_defined(&self) -> bool {
        matches!(
            self,
            Shape::Sphere(_)
                | Shape::SphericalShell(_)
                | Shape::Cylinder(_)
                | Shape::CylindricalShell(_)
                | Shape::Curve(_)
                | Shape::Wire(_)
                | Shape::Point(_)
        )
    }

    /// The shape with its defining points expressed in Cartesian
    /// coordinates, given the coordinate system `cs` they are written in.
    ///
    /// Point-defined shapes (see [`Shape::is_point_defined`]) in a
    /// cylindrical system get their points converted (expressions are
    /// replaced by their current values); every other shape is returned
    /// unchanged — boxes keep interpreting their corners in `cs`, polygon
    /// based shapes and polyhedra are always Cartesian. After this, every
    /// shape except boxes can be read as Cartesian geometry.
    pub fn to_cartesian(&self, cs: CoordSystem) -> Cow<'_, Shape> {
        if cs == CoordSystem::Cartesian || !self.is_point_defined() {
            return Cow::Borrowed(self);
        }
        let conv = |p: &ScalarPoint| ScalarPoint::from(cs.to_cartesian(p.value()));
        let convs = |v: &[ScalarPoint]| v.iter().map(conv).collect::<Vec<_>>();
        Cow::Owned(match self {
            Shape::Sphere(s) => Shape::Sphere(Sphere {
                center: conv(&s.center),
                ..s.clone()
            }),
            Shape::SphericalShell(s) => Shape::SphericalShell(SphericalShell {
                center: conv(&s.center),
                ..s.clone()
            }),
            Shape::Cylinder(c) => Shape::Cylinder(Cylinder {
                start: conv(&c.start),
                stop: conv(&c.stop),
                ..c.clone()
            }),
            Shape::CylindricalShell(c) => Shape::CylindricalShell(CylindricalShell {
                start: conv(&c.start),
                stop: conv(&c.stop),
                ..c.clone()
            }),
            Shape::Curve(c) => Shape::Curve(Curve {
                points: convs(&c.points),
            }),
            Shape::Wire(w) => Shape::Wire(Wire {
                points: convs(&w.points),
                ..w.clone()
            }),
            Shape::Point(p) => Shape::Point(PointShape {
                point: conv(&p.point),
            }),
            other => other.clone(),
        })
    }

    /// Element name used in CSXCAD files.
    pub fn type_name(&self) -> &'static str {
        match self {
            Shape::Box(_) => "Box",
            Shape::MultiBox(_) => "MultiBox",
            Shape::Sphere(_) => "Sphere",
            Shape::SphericalShell(_) => "SphericalShell",
            Shape::Cylinder(_) => "Cylinder",
            Shape::CylindricalShell(_) => "CylindricalShell",
            Shape::Polygon(_) => "Polygon",
            Shape::LinPoly(_) => "LinPoly",
            Shape::RotPoly(_) => "RotPoly",
            Shape::Curve(_) => "Curve",
            Shape::Wire(_) => "Wire",
            Shape::Point(_) => "Point",
            Shape::Polyhedron(_) => "Polyhedron",
            Shape::PolyhedronReader(_) => "PolyhedronReader",
            Shape::UserDefined(_) => "UserDefined",
        }
    }

    /// Inside test for a local Cartesian point; `cs` is the coordinate system
    /// the shape's values are expressed in.
    pub fn contains_local(&self, p: Vec3, cs: CoordSystem) -> bool {
        match self {
            Shape::Box(s) => s.inside(p, cs),
            Shape::MultiBox(s) => s.inside(p, cs),
            Shape::Sphere(s) => s.inside(p, cs),
            Shape::SphericalShell(s) => s.inside(p, cs),
            Shape::Cylinder(s) => s.inside(p, cs),
            Shape::CylindricalShell(s) => s.inside(p, cs),
            Shape::Polygon(s) => s.inside(p),
            Shape::LinPoly(s) => s.inside(p),
            Shape::RotPoly(s) => s.inside(p),
            Shape::Curve(s) => s.inside(p),
            Shape::Wire(s) => s.inside(p, cs),
            Shape::Point(s) => s.inside(p),
            Shape::Polyhedron(s) => s.contains(p),
            Shape::PolyhedronReader(s) => s.mesh.as_ref().is_some_and(|m| m.contains(p)),
            Shape::UserDefined(s) => s.inside(p),
        }
    }

    /// Local Cartesian bounding box.
    pub fn local_bounding_box(&self, cs: CoordSystem) -> BoundingBox {
        match self {
            Shape::Box(s) => s.bbox(cs),
            Shape::MultiBox(s) => s.bbox(cs),
            Shape::Sphere(s) => s.bbox(cs),
            Shape::SphericalShell(s) => s.bbox(cs),
            Shape::Cylinder(s) => s.bbox(cs),
            Shape::CylindricalShell(s) => s.bbox(cs),
            Shape::Polygon(s) => s.bbox(),
            Shape::LinPoly(s) => s.bbox(),
            Shape::RotPoly(s) => s.bbox(),
            Shape::Curve(s) => s.bbox(cs),
            Shape::Wire(s) => s.bbox(cs),
            Shape::Point(s) => s.bbox(cs),
            Shape::Polyhedron(s) => s.bounding_box(),
            Shape::PolyhedronReader(s) => s
                .mesh
                .as_ref()
                .map_or(BoundingBox::EMPTY, |m| m.bounding_box()),
            Shape::UserDefined(s) => s.bbox(),
        }
    }

    /// Bounds expressed in the shape's own coordinate system, for the shapes
    /// whose interpretation depends on it (boxes and multi-boxes).
    fn native_bbox(&self) -> Option<BoundingBox> {
        match self {
            Shape::Box(s) => Some(s.native()),
            Shape::MultiBox(s) => Some(s.native()),
            _ => None,
        }
    }

    /// Topological dimension: 0 (point), 1 (curve), 2 (sheet) or 3 (solid).
    pub fn dimension(&self) -> usize {
        match self {
            Shape::Box(s) => s.dimension(),
            Shape::MultiBox(s) => s.dimension(),
            Shape::Sphere(s) => usize::from(s.radius.value() > 0.0) * 3,
            Shape::SphericalShell(s) => s.dimension(),
            Shape::Cylinder(s) => usize::from(s.radius.value() > 0.0) * 3,
            Shape::CylindricalShell(s) => s.dimension(),
            Shape::Polygon(s) => s.dimension(),
            Shape::LinPoly(s) => s.dimension(),
            Shape::RotPoly(_) => 3,
            Shape::Curve(s) => s.dimension(),
            Shape::Wire(s) => s.dimension(),
            Shape::Point(_) => 0,
            Shape::Polyhedron(s) => 2 + usize::from(s.is_closed()),
            Shape::PolyhedronReader(s) => s
                .mesh
                .as_ref()
                .map_or(0, |m| 2 + usize::from(m.is_closed())),
            Shape::UserDefined(_) => 3,
        }
    }

    /// Re-evaluates parametric values.
    pub fn update(&mut self, p: &ParameterSet) -> Result<()> {
        match self {
            Shape::Box(s) => s.update(p),
            Shape::MultiBox(s) => s.update(p),
            Shape::Sphere(s) => s.update(p),
            Shape::SphericalShell(s) => s.update(p),
            Shape::Cylinder(s) => s.update(p),
            Shape::CylindricalShell(s) => s.update(p),
            Shape::Polygon(s) => s.update(p),
            Shape::LinPoly(s) => s.update(p),
            Shape::RotPoly(s) => s.update(p),
            Shape::Curve(s) => s.update(p),
            Shape::Wire(s) => s.update(p),
            Shape::Point(s) => s.update(p),
            Shape::Polyhedron(_) | Shape::PolyhedronReader(_) => Ok(()),
            Shape::UserDefined(s) => s.update(p),
        }
    }
}

macro_rules! impl_from_shape {
    ($($t:ty => $v:ident),* $(,)?) => {$(
        impl From<$t> for Shape {
            fn from(s: $t) -> Shape { Shape::$v(s) }
        }
        impl From<$t> for Primitive {
            fn from(s: $t) -> Primitive { Primitive::new(Shape::$v(s)) }
        }
    )*};
}

impl_from_shape! {
    BoxShape => Box, MultiBox => MultiBox, Sphere => Sphere,
    SphericalShell => SphericalShell, Cylinder => Cylinder,
    CylindricalShell => CylindricalShell, Polygon => Polygon, LinPoly => LinPoly,
    RotPoly => RotPoly, Curve => Curve, Wire => Wire, PointShape => Point,
    Polyhedron => Polyhedron, PolyhedronReader => PolyhedronReader,
    UserDefined => UserDefined,
}

/// A shape placed in a model.
#[derive(Debug, Clone, PartialEq)]
pub struct Primitive {
    /// Overlap resolution priority (higher wins).
    pub priority: i32,
    /// Coordinate system of the shape's values; `None` inherits the
    /// structure's coordinate system.
    pub coord_system: Option<CoordSystem>,
    /// Placement transformation (applied in Cartesian space).
    pub transform: Transform,
    /// The geometry.
    pub shape: Shape,
    /// Unrecognised attributes preserved from files.
    pub extra_attributes: Vec<(String, String)>,
    /// Unrecognised child elements preserved from files.
    pub extra_elements: Vec<crate::io::dom::Element>,
}

impl Primitive {
    /// Primitive with priority 0, inherited coordinate system and no
    /// transformation.
    pub fn new(shape: impl Into<Shape>) -> Self {
        Primitive {
            priority: 0,
            coord_system: None,
            transform: Transform::new(),
            shape: shape.into(),
            extra_attributes: Vec::new(),
            extra_elements: Vec::new(),
        }
    }

    /// Builder: sets the priority.
    #[must_use]
    pub fn with_priority(mut self, priority: i32) -> Self {
        self.priority = priority;
        self
    }

    /// Builder: sets the coordinate system.
    #[must_use]
    pub fn with_coord_system(mut self, cs: CoordSystem) -> Self {
        self.coord_system = Some(cs);
        self
    }

    /// Builder: sets the transformation.
    #[must_use]
    pub fn with_transform(mut self, t: Transform) -> Self {
        self.transform = t;
        self
    }

    /// Effective coordinate system given the structure default.
    #[inline]
    pub fn effective_cs(&self, default: CoordSystem) -> CoordSystem {
        self.coord_system.unwrap_or(default)
    }

    /// The shape in local Cartesian form (see [`Shape::to_cartesian`]),
    /// given the structure's default coordinate system.
    pub fn cartesian_shape(&self, default_cs: CoordSystem) -> Cow<'_, Shape> {
        self.shape.to_cartesian(self.effective_cs(default_cs))
    }

    /// Inside test for `p` expressed in `p_cs`, with `default_cs` used when
    /// the primitive does not define its own coordinate system.
    ///
    /// As in CSXCAD, an untransformed box (or multi-box) whose values are in
    /// the query's (cylindrical) coordinate system compares the query
    /// coordinates as given: a negative radius `ρ < 0` is outside a
    /// cylindrical box. Every other shape receives the point converted to
    /// Cartesian coordinates, where `(-ρ, α)` is the same point as
    /// `(ρ, α + π)`.
    pub fn contains_in(&self, p: Vec3, p_cs: CoordSystem, default_cs: CoordSystem) -> bool {
        let scs = self.shape_cs(default_cs);
        if p_cs == scs && scs != CoordSystem::Cartesian && self.transform.is_empty() {
            match &self.shape {
                Shape::Box(b) => return b.inside_native(p, scs),
                Shape::MultiBox(m) => return m.inside_native(p),
                _ => {}
            }
        }
        let mut c = p_cs.to_cartesian(p);
        if !self.transform.is_empty() {
            if !self.transform.is_invertible() {
                return false;
            }
            c = self.transform.apply_inverse(c);
        }
        self.shape.contains_local(c, scs)
    }

    /// The coordinate system in which the shape's values are interpreted:
    /// the effective coordinate system, except for multi-boxes, which always
    /// use the structure's (`default_cs`).
    pub fn shape_cs(&self, default_cs: CoordSystem) -> CoordSystem {
        match self.shape {
            Shape::MultiBox(_) => default_cs,
            _ => self.effective_cs(default_cs),
        }
    }

    /// Inside test for `p` expressed in `p_cs` (default coordinate system:
    /// Cartesian).
    pub fn contains(&self, p: Vec3, p_cs: CoordSystem) -> bool {
        self.contains_in(p, p_cs, CoordSystem::Cartesian)
    }

    /// World-space Cartesian bounding box.
    pub fn bounding_box(&self, default_cs: CoordSystem) -> BoundingBox {
        let b = self.shape.local_bounding_box(self.shape_cs(default_cs));
        if self.transform.is_empty() || b.is_empty() || !(b.min.is_finite() && b.max.is_finite()) {
            // Unbounded shapes stay unbounded (transforming infinite corners
            // would produce NaN).
            b
        } else {
            b.transformed(self.transform.matrix())
        }
    }

    /// Bounding box expressed in `target` coordinates. For cylindrical
    /// targets the result is `(ρ, α, z)` ranges; it is exact for boxes and
    /// multi-boxes defined in cylindrical coordinates and conservative
    /// otherwise.
    pub fn bounding_box_in(&self, target: CoordSystem, default_cs: CoordSystem) -> BoundingBox {
        let cs = self.shape_cs(default_cs);
        if cs == target
            && self.transform.is_empty()
            && let Some(b) = self.shape.native_bbox()
        {
            // Multi-boxes compare angles without wrapping, so their ranges
            // are kept as written.
            return match target {
                CoordSystem::Cylindrical if matches!(self.shape, Shape::Box(_)) => {
                    normalize_alpha(b)
                }
                _ => b,
            };
        }
        let b = self.bounding_box(default_cs);
        match target {
            CoordSystem::Cartesian => b,
            CoordSystem::Cylindrical => cartesian_bbox_to_cylindrical(&b),
        }
    }

    /// Topological dimension of the shape (0..=3).
    pub fn dimension(&self) -> usize {
        self.shape.dimension()
    }

    /// Re-evaluates parametric values of the shape and the transformation.
    pub fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.shape.update(p)?;
        self.transform.update(p)
    }
}

/// Shifts the angular range of a cylindrical box so that it starts in
/// `[-π, π)`; ranges of a full turn or more become `[-π, π]`.
fn normalize_alpha(mut b: BoundingBox) -> BoundingBox {
    use core::f64::consts::{PI, TAU};
    if b.is_empty() || !b.min[1].is_finite() || !b.max[1].is_finite() {
        return b;
    }
    if b.max[1] - b.min[1] >= TAU {
        b.min[1] = -PI;
        b.max[1] = PI;
        return b;
    }
    let shift = ((b.min[1] + PI) / TAU).floor() * TAU;
    b.min[1] -= shift;
    b.max[1] -= shift;
    b
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::{FRAC_PI_2, PI, TAU};

    #[test]
    fn cartesian_box() {
        let b = Primitive::new(BoxShape::new([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]));
        assert!(b.contains(Vec3::new(0.5, 1.0, 1.0), CoordSystem::Cartesian));
        assert!(b.contains(Vec3::new(1.0, 2.0, 3.0), CoordSystem::Cartesian));
        assert!(!b.contains(Vec3::new(1.1, 1.0, 1.0), CoordSystem::Cartesian));
        assert_eq!(b.dimension(), 3);
        // query in cylindrical coordinates: rho=1, alpha=pi/4 → (0.707,0.707)
        assert!(b.contains(Vec3::new(1.0, PI / 4.0, 1.0), CoordSystem::Cylindrical));
    }

    #[test]
    fn cylindrical_box() {
        let b = Primitive::new(BoxShape::new([1.0, -0.5, 0.0], [2.0, 0.5, 1.0]))
            .with_coord_system(CoordSystem::Cylindrical);
        assert!(b.contains(Vec3::new(1.5, 0.0, 0.5), CoordSystem::Cartesian));
        assert!(!b.contains(Vec3::new(0.0, 1.5, 0.5), CoordSystem::Cartesian));
        assert!(b.contains(Vec3::new(1.5, 0.4, 0.5), CoordSystem::Cylindrical));
        let n = b.bounding_box_in(CoordSystem::Cylindrical, CoordSystem::Cartesian);
        assert_eq!(n.min, Vec3::new(1.0, -0.5, 0.0));
        let c = b.bounding_box(CoordSystem::Cartesian);
        assert!((c.max[0] - 2.0).abs() < 1e-12);
        assert!((c.min[0] - 0.5f64.cos()).abs() < 1e-12);
        // wrapping sector around ±π
        let w = Primitive::new(BoxShape::new([1.0, 3.0, 0.0], [2.0, 3.3, 1.0]))
            .with_coord_system(CoordSystem::Cylindrical);
        assert!(w.contains(Vec3::new(-1.5, 0.0, 0.5), CoordSystem::Cartesian));
        assert!(w.contains(Vec3::new(-1.5, -0.1, 0.5), CoordSystem::Cartesian));
    }

    #[test]
    fn cylindrical_bbox_alpha_normalised() {
        let b = Primitive::new(BoxShape::new([1.0, 6.5, 0.0], [2.0, 7.0, 1.0]))
            .with_coord_system(CoordSystem::Cylindrical);
        let n = b.bounding_box_in(CoordSystem::Cylindrical, CoordSystem::Cartesian);
        assert!((n.min[1] - (6.5 - TAU)).abs() < 1e-12);
        assert!((n.max[1] - (7.0 - TAU)).abs() < 1e-12);
        let full = Primitive::new(BoxShape::new([1.0, 0.0, 0.0], [2.0, 10.0, 1.0]))
            .with_coord_system(CoordSystem::Cylindrical);
        let n = full.bounding_box_in(CoordSystem::Cylindrical, CoordSystem::Cartesian);
        assert_eq!((n.min[1], n.max[1]), (-PI, PI));
    }

    #[test]
    fn singular_and_unbounded_transforms() {
        let mut t = Transform::new();
        t.scale(0.0);
        let b = Primitive::new(BoxShape::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])).with_transform(t);
        // Previously the identity was used as the inverse of a singular
        // transform, so the untransformed box was reported.
        assert!(!b.contains(Vec3::new(0.5, 0.5, 0.5), CoordSystem::Cartesian));
        let mut t = Transform::new();
        t.rotate_axis(2, 0.3);
        let u = Primitive::new(UserDefined::new("1")).with_transform(t);
        let bb = u.bounding_box(CoordSystem::Cartesian);
        assert!(bb.min[0] == f64::NEG_INFINITY && bb.max[2] == f64::INFINITY);
    }

    #[test]
    fn transformed() {
        let mut t = Transform::new();
        t.rotate_axis(2, FRAC_PI_2)
            .translate(Vec3::new(10.0, 0.0, 0.0));
        let b = Primitive::new(BoxShape::new([0.0, 0.0, 0.0], [2.0, 1.0, 1.0])).with_transform(t);
        // local (1.5, 0.5) → rotate → (-0.5, 1.5) → translate → (9.5, 1.5)
        assert!(b.contains(Vec3::new(9.5, 1.5, 0.5), CoordSystem::Cartesian));
        assert!(!b.contains(Vec3::new(1.5, 0.5, 0.5), CoordSystem::Cartesian));
        let bb = b.bounding_box(CoordSystem::Cartesian);
        assert!((bb.min - Vec3::new(9.0, 0.0, 0.0)).norm() < 1e-12);
        assert!((bb.max - Vec3::new(10.0, 2.0, 1.0)).norm() < 1e-12);
    }

    #[test]
    fn round_shapes() {
        let s = Primitive::new(Sphere::new([1.0, 0.0, 0.0], 1.0));
        assert!(s.contains(Vec3::new(1.5, 0.5, 0.5), CoordSystem::Cartesian));
        assert!(!s.contains(Vec3::new(1.5, 0.8, 0.8), CoordSystem::Cartesian));
        let c = Primitive::new(Cylinder::new([0.0, 0.0, 0.0], [0.0, 0.0, 2.0], 0.5));
        assert!(c.contains(Vec3::new(0.3, 0.3, 2.0), CoordSystem::Cartesian));
        assert!(!c.contains(Vec3::new(0.3, 0.3, 2.1), CoordSystem::Cartesian));
        assert!(!c.contains(Vec3::new(0.4, 0.4, 1.0), CoordSystem::Cartesian));
        let bb = c.bounding_box(CoordSystem::Cartesian);
        assert_eq!(bb.min, Vec3::new(-0.5, -0.5, 0.0));
        let sh = Primitive::new(CylindricalShell {
            start: [0.0, 0.0, 0.0].into(),
            stop: [1.0, 0.0, 0.0].into(),
            radius: 1.0.into(),
            shell_width: 0.2.into(),
        });
        assert!(sh.contains(Vec3::new(0.5, 0.0, 1.05), CoordSystem::Cartesian));
        assert!(!sh.contains(Vec3::new(0.5, 0.0, 0.5), CoordSystem::Cartesian));
        let ss = Primitive::new(SphericalShell {
            center: [0.0, 0.0, 0.0].into(),
            radius: 1.0.into(),
            shell_width: 0.2.into(),
        });
        assert!(ss.contains(Vec3::new(0.0, 0.95, 0.0), CoordSystem::Cartesian));
        assert!(!ss.contains(Vec3::new(0.0, 0.85, 0.0), CoordSystem::Cartesian));
    }

    #[test]
    fn wires_and_points() {
        let w = Primitive::new(Wire::new(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
            0.1,
        ));
        assert!(w.contains(Vec3::new(1.05, 0.5, 0.0), CoordSystem::Cartesian));
        assert!(!w.contains(Vec3::new(0.5, 0.5, 0.0), CoordSystem::Cartesian));
        // Curves and points never contain anything (verified in CSXCAD).
        let c = Primitive::new(Curve::new([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]));
        assert!(!c.contains(Vec3::new(0.5, 0.0, 0.0), CoordSystem::Cartesian));
        assert_eq!(c.dimension(), 1);
        let p = Primitive::new(PointShape::new([1.0, 2.0, 3.0]));
        assert!(!p.contains(Vec3::new(1.0, 2.0, 3.0), CoordSystem::Cartesian));
        assert_eq!(p.dimension(), 0);
    }

    /// Coordinate system semantics observed in CSXCAD 0.7.
    #[test]
    fn cylindrical_definitions() {
        let cyl = CoordSystem::Cylindrical;
        let half = core::f64::consts::FRAC_PI_2;
        let cart = |x: f64, y: f64, z: f64| Vec3::new(x, y, z);
        // Defining points of round shapes are converted: center (1, π/2, 0)
        // is the Cartesian point (0, 1, 0).
        let s = Primitive::new(Sphere::new([1.0, half, 0.0], 0.1)).with_coord_system(cyl);
        assert!(s.contains(cart(0.0, 1.0, 0.0), CoordSystem::Cartesian));
        assert!(!s.contains(cart(1.0, half, 0.0), CoordSystem::Cartesian));
        let c = Primitive::new(Cylinder::new([1.0, half, 0.0], [1.0, half, 1.0], 0.1))
            .with_coord_system(cyl);
        assert!(c.contains(cart(0.0, 1.0, 0.5), CoordSystem::Cartesian));
        let w = Primitive::new(Wire::new([[1.0, half, 0.0], [1.0, half, 1.0]], 0.1))
            .with_coord_system(cyl);
        assert!(w.contains(cart(0.0, 1.0, 0.5), CoordSystem::Cartesian));
        let p = Primitive::new(PointShape::new([2.0, half, 1.0])).with_coord_system(cyl);
        let b = p.bounding_box(CoordSystem::Cartesian);
        assert!((b.min - cart(0.0, 2.0, 1.0)).norm() < 1e-12);
        // Multi-boxes stay Cartesian: (1.9, 0.6) is outside the box
        // [1,2]×[-0.5,0.5] even though it lies in the sector.
        let m = Primitive::new(MultiBox {
            boxes: vec![([1.0, -0.5, 0.0].into(), [2.0, 0.5, 1.0].into())],
        })
        .with_coord_system(cyl);
        assert!(!m.contains(cart(1.9, 0.6, 0.5), CoordSystem::Cartesian));
        assert!(m.contains(cart(1.5, 0.0, 0.5), CoordSystem::Cartesian));
        // A box is a sector.
        let bx =
            Primitive::new(BoxShape::new([1.0, -0.5, 0.0], [2.0, 0.5, 1.0])).with_coord_system(cyl);
        assert!(bx.contains(cart(1.9, 0.6, 0.5), CoordSystem::Cartesian));
        // Primitives without an explicit system inherit the structure's.
        let implicit = Primitive::new(BoxShape::new([1.0, -0.5, 0.0], [2.0, 0.5, 1.0]));
        assert!(implicit.contains_in(Vec3::new(1.9, 0.3, 0.5), cyl, cyl));
        assert!(!implicit.contains_in(Vec3::new(1.9, 0.6, 0.5), cyl, cyl));
    }

    #[test]
    fn user_defined() {
        let mut u = UserDefined::new("x^2 + y^2 < R^2 & z > 0");
        assert!(!u.inside(Vec3::new(0.1, 0.1, 1.0)));
        let mut ps = ParameterSet::new();
        ps.set("R", 1.0);
        u.update(&ps).unwrap();
        let p = Primitive::new(u);
        assert!(p.contains(Vec3::new(0.1, 0.1, 1.0), CoordSystem::Cartesian));
        assert!(!p.contains(Vec3::new(0.1, 0.1, -1.0), CoordSystem::Cartesian));
        assert!(p.contains(Vec3::new(0.5, 0.3, 1.0), CoordSystem::Cylindrical));
    }

    #[test]
    fn polyhedron_primitive() {
        let mut t = Transform::new();
        t.scale(2.0);
        let p = Primitive::new(polyhedron::unit_cube()).with_transform(t);
        assert!(p.contains(Vec3::new(1.5, 1.5, 1.5), CoordSystem::Cartesian));
        assert_eq!(p.dimension(), 3);
        assert_eq!(
            p.bounding_box(CoordSystem::Cartesian).max,
            Vec3::new(2.0, 2.0, 2.0)
        );
    }

    #[test]
    fn cartesian_shape_conversion() {
        let cyl = CoordSystem::Cylindrical;
        // Sphere centred at rho=2, alpha=pi/2 → (0, 2, 1).
        let s = Primitive::new(Sphere::new([2.0, FRAC_PI_2, 1.0], 0.5)).with_coord_system(cyl);
        let c = s.cartesian_shape(CoordSystem::Cartesian);
        let Shape::Sphere(sp) = &*c else { panic!() };
        assert!((sp.center.value() - Vec3::new(0.0, 2.0, 1.0)).norm() < 1e-12);
        assert!(s.contains(Vec3::new(0.0, 2.0, 1.2), CoordSystem::Cartesian));
        // Inherited coordinate system.
        let w = Primitive::new(Wire::new([[1.0, 0.0, 0.0], [1.0, PI, 0.0]], 0.1));
        let Shape::Wire(cw) = &*w.cartesian_shape(cyl) else {
            panic!()
        };
        assert!((cw.points[1].value() - Vec3::new(-1.0, 0.0, 0.0)).norm() < 1e-12);
        assert_eq!(cw.radius.value(), 0.1);
        // Boxes, polygons and Cartesian shapes are borrowed unchanged.
        let b =
            Primitive::new(BoxShape::new([1.0, 0.0, 0.0], [2.0, 1.0, 1.0])).with_coord_system(cyl);
        assert!(matches!(b.cartesian_shape(cyl), Cow::Borrowed(_)));
        assert!(matches!(
            s.cartesian_shape(CoordSystem::Cartesian),
            Cow::Owned(_)
        ));
        let cart = Primitive::new(Sphere::new([2.0, 0.0, 0.0], 1.0));
        assert!(matches!(
            cart.cartesian_shape(CoordSystem::Cartesian),
            Cow::Borrowed(_)
        ));
    }

    /// Raw versus converted cylindrical queries (verified black-box against
    /// CSXCAD 0.7).
    #[test]
    fn cylindrical_query_semantics() {
        let cyl = CoordSystem::Cylindrical;
        // A cylindrical box compares the query's (ρ, α, z) as given: a
        // negative radius is outside even if (-ρ, α) folds into the sector.
        let sector =
            Primitive::new(BoxShape::new([0.5, -0.5, 0.0], [2.0, 0.5, 1.0])).with_coord_system(cyl);
        assert!(sector.contains_in(Vec3::new(1.0, 0.0, 0.5), cyl, cyl));
        assert!(!sector.contains_in(Vec3::new(-1.0, PI, 0.5), cyl, cyl));
        assert!(!sector.contains_in(Vec3::new(-1.0, 0.0, 0.5), cyl, cyl));
        // ... with the angle taken modulo 2π.
        assert!(sector.contains_in(Vec3::new(1.0, 2.0 * PI, 0.5), cyl, cyl));
        // Cartesian queries go through the conversion.
        let cart = CoordSystem::Cartesian;
        assert!(sector.contains_in(Vec3::new(1.0, 0.0, 0.5), cart, cyl));
        // Every other shape folds a negative radius through the Cartesian
        // conversion: (-1, 0) is the point (-1, 0, z).
        let b = Primitive::new(BoxShape::new([-2.0, -1.0, 0.0], [0.0, 1.0, 1.0]))
            .with_coord_system(cart);
        assert!(b.contains_in(Vec3::new(-1.0, 0.0, 0.5), cyl, cyl));
        assert!(b.contains_in(Vec3::new(1.0, PI, 0.5), cyl, cyl));
        let s = Primitive::new(Sphere::new([-1.0, 0.0, 0.0], 0.5)).with_coord_system(cart);
        assert!(s.contains_in(Vec3::new(-1.0, 0.0, 0.0), cyl, cyl));
        let all = Primitive::new(UserDefined::new("1"));
        assert!(all.contains_in(Vec3::new(-1.0, 0.3, 0.0), cyl, cyl));
        // A multi-box uses the structure's coordinate system whatever its
        // own, and compares without angle wrapping.
        let m = Primitive::new(MultiBox {
            boxes: vec![([1.0, -1.0, 0.0].into(), [2.0, 1.0, 1.0].into())],
        })
        .with_coord_system(cart);
        assert!(m.contains_in(Vec3::new(1.5, 0.0, 0.5), cyl, cyl));
        assert!(!m.contains_in(Vec3::new(1.5, 2.0 * PI, 0.5), cyl, cyl));
        assert!(!m.contains_in(Vec3::new(-1.5, PI, 0.5), cyl, cyl));
        // In a Cartesian structure the same multi-box is Cartesian, even
        // when flagged cylindrical.
        let m = m.with_coord_system(cyl);
        assert!(m.contains_in(Vec3::new(1.5, 0.5, 0.5), cart, cart));
        assert!(!m.contains_in(Vec3::new(1.5, 1.5, 0.5), cart, cart));
        assert_eq!(m.shape_cs(cart), cart);
        assert_eq!(m.shape_cs(cyl), cyl);
        // Its cylindrical bounds are the raw ranges.
        let nb = m.bounding_box_in(cyl, cyl);
        assert_eq!(nb.min, Vec3::new(1.0, -1.0, 0.0));
        assert_eq!(nb.max, Vec3::new(2.0, 1.0, 1.0));
    }

    #[test]
    fn constructors() {
        let m = MultiBox::new([
            ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            ([2.0, 0.0, 0.0], [3.0, 1.0, 1.0]),
        ]);
        let mut n = MultiBox::default();
        n.push([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        n.push([2.0, 0.0, 0.0], [3.0, 1.0, 1.0]);
        assert_eq!(m, n);
        assert!(Primitive::new(m).contains(Vec3::new(2.5, 0.5, 0.5), CoordSystem::Cartesian));
        let s = SphericalShell::new([0.0, 0.0, 0.0], 1.0, 0.2);
        assert_eq!(s.shell_width.value(), 0.2);
        let c = CylindricalShell::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0, 0.2);
        assert!(Primitive::new(c).contains(Vec3::new(1.05, 0.0, 0.5), CoordSystem::Cartesian));
        let u = UserDefined::new("x < 0").with_base([1.0, 0.0, 0.0]);
        assert!(Primitive::new(u).contains(Vec3::new(0.5, 0.0, 0.0), CoordSystem::Cartesian));
    }
}
