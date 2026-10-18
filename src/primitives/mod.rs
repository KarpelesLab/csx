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
//! tested against the shape. Only boxes (and multi-boxes) interpret their
//! values in the primitive's coordinate system — in cylindrical coordinates
//! they describe annular sectors; every other shape is always Cartesian.

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

use crate::coords::CoordSystem;
use crate::error::Result;
use crate::math::{BoundingBox, Vec3};
use crate::param::ParameterSet;
use crate::transform::Transform;

/// Mesh file formats understood by [`PolyhedronReader`].
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

    /// Parses a numeric id or a file type name.
    pub fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "1" | "stl" => Some(MeshFileType::Stl),
            "2" | "ply" => Some(MeshFileType::Ply),
            _ => None,
        }
    }

    /// Guesses the type from a file name extension.
    pub fn from_path(path: &str) -> Option<Self> {
        let ext = path.rsplit('.').next()?;
        MeshFileType::parse(ext)
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
            Shape::MultiBox(s) => s.inside(p),
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
    pub fn local_bbox(&self, cs: CoordSystem) -> BoundingBox {
        match self {
            Shape::Box(s) => s.bbox(cs),
            Shape::MultiBox(s) => s.bbox(),
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
    /// whose interpretation depends on it (boxes).
    fn native_bbox(&self) -> Option<BoundingBox> {
        match self {
            Shape::Box(s) => Some(s.native()),
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
    pub fn with_priority(mut self, priority: i32) -> Self {
        self.priority = priority;
        self
    }

    /// Builder: sets the coordinate system.
    pub fn with_coord_system(mut self, cs: CoordSystem) -> Self {
        self.coord_system = Some(cs);
        self
    }

    /// Builder: sets the transformation.
    pub fn with_transform(mut self, t: Transform) -> Self {
        self.transform = t;
        self
    }

    /// Effective coordinate system given the structure default.
    #[inline]
    pub fn effective_cs(&self, default: CoordSystem) -> CoordSystem {
        self.coord_system.unwrap_or(default)
    }

    /// Inside test for `p` expressed in `p_cs`, with `default_cs` used when
    /// the primitive does not define its own coordinate system.
    pub fn contains_in(&self, p: Vec3, p_cs: CoordSystem, default_cs: CoordSystem) -> bool {
        let mut c = p_cs.to_cartesian(p);
        if !self.transform.is_empty() {
            if !self.transform.is_invertible() {
                return false;
            }
            c = self.transform.apply_inverse(c);
        }
        self.shape.contains_local(c, self.effective_cs(default_cs))
    }

    /// Inside test for `p` expressed in `p_cs` (default coordinate system:
    /// Cartesian).
    pub fn contains(&self, p: Vec3, p_cs: CoordSystem) -> bool {
        self.contains_in(p, p_cs, CoordSystem::Cartesian)
    }

    /// World-space Cartesian bounding box.
    pub fn bounding_box(&self, default_cs: CoordSystem) -> BoundingBox {
        let b = self.shape.local_bbox(self.effective_cs(default_cs));
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
        let cs = self.effective_cs(default_cs);
        if cs == target
            && self.transform.is_empty()
            && let Some(b) = self.shape.native_bbox()
        {
            return match target {
                CoordSystem::Cartesian => b,
                CoordSystem::Cylindrical => normalize_alpha(b),
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
}
