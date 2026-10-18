//! The continuous structure: the complete geometry and material model.

use std::path::PathBuf;

use crate::coords::CoordSystem;
use crate::error::{Error, Result};
use crate::grid::RectilinearGrid;
use crate::math::{BoundingBox, Vec3};
use crate::param::{ParameterSet, Scalar};
use crate::primitives::{Primitive, Shape};
use crate::properties::{MaterialQuantity, Property, PropertyKind, TypeMask};

/// Material filling all space not covered by any material primitive.
#[derive(Debug, Clone, PartialEq)]
pub struct BackgroundMaterial {
    /// Relative permittivity.
    pub epsilon: Scalar,
    /// Relative permeability.
    pub mue: Scalar,
    /// Electric conductivity.
    pub kappa: Scalar,
    /// Magnetic conductivity.
    pub sigma: Scalar,
}

impl Default for BackgroundMaterial {
    fn default() -> Self {
        BackgroundMaterial {
            epsilon: Scalar::new(1.0),
            mue: Scalar::new(1.0),
            kappa: Scalar::new(0.0),
            sigma: Scalar::new(0.0),
        }
    }
}

impl BackgroundMaterial {
    fn update(&mut self, p: &ParameterSet) -> Result<()> {
        self.epsilon.update(p)?;
        self.mue.update(p)?;
        self.kappa.update(p)?;
        self.sigma.update(p)
    }

    /// Value of a quantity.
    pub fn value(&self, q: MaterialQuantity) -> f64 {
        match q {
            MaterialQuantity::Epsilon => self.epsilon.value(),
            MaterialQuantity::Mue => self.mue.value(),
            MaterialQuantity::Kappa => self.kappa.value(),
            MaterialQuantity::Sigma => self.sigma.value(),
        }
    }
}

/// Opaque handle of a property inside a [`ContinuousStructure`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct PropertyId(pub u32);

/// The result of a priority query: a property and the winning primitive.
#[derive(Debug, Clone, Copy)]
pub struct Hit<'a> {
    /// The property.
    pub property: &'a Property,
    /// The primitive that contains the point.
    pub primitive: &'a Primitive,
    /// Index of the primitive within the property.
    pub primitive_index: usize,
}

/// Electromagnetic medium at a point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Medium {
    /// Relative permittivity per axis.
    pub epsilon: [f64; 3],
    /// Relative permeability per axis.
    pub mue: [f64; 3],
    /// Electric conductivity per axis.
    pub kappa: [f64; 3],
    /// Magnetic conductivity per axis.
    pub sigma: [f64; 3],
}

/// The complete model: properties with their primitives, the mesh, the
/// parameters and the background material.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ContinuousStructure {
    properties: Vec<Property>,
    next_id: u32,
    /// Mesh.
    pub grid: RectilinearGrid,
    /// Named parameters used by expressions.
    pub parameters: ParameterSet,
    /// Background material.
    pub background: BackgroundMaterial,
    /// Default coordinate system of primitives.
    pub coord_system: CoordSystem,
    /// Directory against which relative file names (mesh files) are
    /// resolved.
    pub base_dir: Option<PathBuf>,
    /// Unrecognised attributes preserved from files.
    pub extra_attributes: Vec<(String, String)>,
}

impl ContinuousStructure {
    /// Empty structure.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a property and returns its handle. Parametric values that do not
    /// depend on parameters are usable immediately; call
    /// [`ContinuousStructure::update`] after adding parametric content.
    pub fn add_property(&mut self, mut p: Property) -> PropertyId {
        self.next_id += 1;
        p.id = self.next_id;
        self.properties.push(p);
        PropertyId(self.next_id)
    }

    /// Adds a primitive to the property `id`.
    pub fn add_primitive(
        &mut self,
        id: PropertyId,
        prim: impl Into<Primitive>,
    ) -> Result<&mut Primitive> {
        let p = self
            .property_mut(id)
            .ok_or_else(|| Error::Invalid(format!("no property {id:?}")))?;
        Ok(p.add_primitive(prim))
    }

    /// Removes a property, returning it.
    pub fn remove_property(&mut self, id: PropertyId) -> Option<Property> {
        let i = self.properties.iter().position(|p| p.id == id.0)?;
        let mut p = self.properties.remove(i);
        p.id = 0;
        Some(p)
    }

    /// All properties in insertion order.
    pub fn properties(&self) -> &[Property] {
        &self.properties
    }

    /// Mutable iterator over all properties.
    pub fn properties_mut(&mut self) -> impl Iterator<Item = &mut Property> {
        self.properties.iter_mut()
    }

    /// Properties matching `mask`.
    pub fn properties_of(&self, mask: impl Into<TypeMask>) -> impl Iterator<Item = &Property> {
        let mask = mask.into();
        self.properties
            .iter()
            .filter(move |p| mask.contains(p.property_type()))
    }

    /// Property by handle.
    pub fn property(&self, id: PropertyId) -> Option<&Property> {
        self.properties.iter().find(|p| p.id == id.0)
    }

    /// Mutable property by handle.
    pub fn property_mut(&mut self, id: PropertyId) -> Option<&mut Property> {
        self.properties.iter_mut().find(|p| p.id == id.0)
    }

    /// First property with the given name.
    pub fn property_by_name(&self, name: &str) -> Option<&Property> {
        self.properties.iter().find(|p| p.name == name)
    }

    /// Mutable first property with the given name.
    pub fn property_by_name_mut(&mut self, name: &str) -> Option<&mut Property> {
        self.properties.iter_mut().find(|p| p.name == name)
    }

    /// Handle of the first property with the given name.
    pub fn property_id(&self, name: &str) -> Option<PropertyId> {
        self.property_by_name(name).map(|p| PropertyId(p.id))
    }

    /// Iterates over all `(property, primitive)` pairs.
    pub fn primitives(&self) -> impl Iterator<Item = (&Property, &Primitive)> {
        self.properties
            .iter()
            .flat_map(|p| p.primitives.iter().map(move |q| (p, q)))
    }

    /// Re-evaluates every parametric value against [`Self::parameters`] and
    /// loads external mesh files of `PolyhedronReader` primitives.
    pub fn update(&mut self) -> Result<()> {
        let params = &self.parameters;
        self.background.update(params)?;
        for p in &mut self.properties {
            p.update(params)?;
            for prim in &mut p.primitives {
                if let Shape::PolyhedronReader(r) = &mut prim.shape
                    && r.mesh.is_none()
                {
                    r.load(self.base_dir.as_deref())?;
                }
            }
        }
        Ok(())
    }

    /// Finds the highest-priority primitive containing `p` (expressed in
    /// `cs`) among properties matching `mask`. On equal priority the
    /// primitive added last wins.
    pub fn property_at(
        &self,
        p: Vec3,
        cs: CoordSystem,
        mask: impl Into<TypeMask>,
    ) -> Option<Hit<'_>> {
        let mask = mask.into();
        let mut best: Option<Hit<'_>> = None;
        for prop in self
            .properties
            .iter()
            .filter(|q| mask.contains(q.property_type()))
        {
            for (i, prim) in prop.primitives.iter().enumerate() {
                if best.is_some_and(|b| b.primitive.priority > prim.priority) {
                    continue;
                }
                if prim.contains_in(p, cs, self.coord_system) {
                    best = Some(Hit {
                        property: prop,
                        primitive: prim,
                        primitive_index: i,
                    });
                }
            }
        }
        best
    }

    /// All properties (matching `mask`) having at least one primitive that
    /// contains `p`, highest priority first.
    pub fn properties_at(
        &self,
        p: Vec3,
        cs: CoordSystem,
        mask: impl Into<TypeMask>,
    ) -> Vec<Hit<'_>> {
        let mask = mask.into();
        let mut hits: Vec<Hit<'_>> = Vec::new();
        for prop in self
            .properties
            .iter()
            .filter(|q| mask.contains(q.property_type()))
        {
            for (i, prim) in prop.primitives.iter().enumerate() {
                if prim.contains_in(p, cs, self.coord_system) {
                    hits.push(Hit {
                        property: prop,
                        primitive: prim,
                        primitive_index: i,
                    });
                }
            }
        }
        // stable sort keeps insertion order; reverse so later ones win ties
        hits.reverse();
        hits.sort_by_key(|h| core::cmp::Reverse(h.primitive.priority));
        hits
    }

    /// The physical medium at `p`: the highest-priority material/metal
    /// primitive, or the background material. Metals are reported as
    /// `None`.
    pub fn medium_at(&self, p: Vec3, cs: CoordSystem) -> Option<Medium> {
        let hit = self.property_at(
            p,
            cs,
            TypeMask::MATERIALS | TypeMask::of(crate::PropertyType::Metal),
        );
        let cart = cs.to_cartesian(p);
        match hit.map(|h| &h.property.kind) {
            Some(PropertyKind::Metal) => None,
            Some(k) => {
                let m = k.material().expect("mask only matches materials");
                Some(Medium {
                    epsilon: m.values_at(MaterialQuantity::Epsilon, cart),
                    mue: m.values_at(MaterialQuantity::Mue, cart),
                    kappa: m.values_at(MaterialQuantity::Kappa, cart),
                    sigma: m.values_at(MaterialQuantity::Sigma, cart),
                })
            }
            None => {
                let b = &self.background;
                Some(Medium {
                    epsilon: [b.epsilon.value(); 3],
                    mue: [b.mue.value(); 3],
                    kappa: [b.kappa.value(); 3],
                    sigma: [b.sigma.value(); 3],
                })
            }
        }
    }

    /// Union of the bounding boxes of all primitives (in `cs`).
    pub fn bounding_box(&self, cs: CoordSystem) -> BoundingBox {
        self.primitives()
            .map(|(_, q)| q.bounding_box_in(cs, self.coord_system))
            .filter(|b| b.min.is_finite() && b.max.is_finite())
            .fold(BoundingBox::EMPTY, |a, b| a.union(&b))
    }

    /// Primitives (of properties matching `mask`) whose bounding box in `cs`
    /// intersects `region`.
    pub fn primitives_in<'a>(
        &'a self,
        region: &'a BoundingBox,
        cs: CoordSystem,
        mask: impl Into<TypeMask>,
    ) -> impl Iterator<Item = (&'a Property, &'a Primitive)> + 'a {
        let mask = mask.into();
        let dcs = self.coord_system;
        self.primitives()
            .filter(move |(p, _)| mask.contains(p.property_type()))
            .filter(move |(_, q)| q.bounding_box_in(cs, dcs).intersects(region))
    }

    /// Builds an index for fast repeated point queries.
    pub fn index(&self, cs: CoordSystem, mask: impl Into<TypeMask>) -> SpatialIndex<'_> {
        SpatialIndex::new(self, cs, mask.into())
    }

    /// Consistency checks: duplicate names, properties without primitives,
    /// primitives with empty or non-finite bounds. Returns human readable
    /// warnings.
    pub fn validate(&self) -> Vec<String> {
        let mut w = Vec::new();
        for (i, p) in self.properties.iter().enumerate() {
            if self.properties[..i].iter().any(|q| q.name == p.name) {
                w.push(format!("duplicate property name {:?}", p.name));
            }
            if p.primitives.is_empty() {
                w.push(format!("property {:?} has no primitives", p.name));
            }
            for (k, q) in p.primitives.iter().enumerate() {
                let b = q.bounding_box(self.coord_system);
                if b.is_empty() {
                    w.push(format!(
                        "property {:?}: primitive #{k} ({}) is empty",
                        p.name,
                        q.shape.type_name()
                    ));
                }
            }
        }
        w
    }
}

struct Entry<'a> {
    bbox: BoundingBox,
    hit: Hit<'a>,
}

/// Point-query accelerator over a borrowed structure: primitives are
/// pre-sorted by priority and culled with their bounding boxes.
pub struct SpatialIndex<'a> {
    entries: Vec<Entry<'a>>,
    cs: CoordSystem,
    default_cs: CoordSystem,
}

impl<'a> SpatialIndex<'a> {
    fn new(s: &'a ContinuousStructure, cs: CoordSystem, mask: TypeMask) -> Self {
        let mut entries: Vec<Entry<'a>> = Vec::new();
        for prop in s
            .properties
            .iter()
            .filter(|q| mask.contains(q.property_type()))
        {
            for (i, prim) in prop.primitives.iter().enumerate() {
                let bbox = prim
                    .bounding_box(s.coord_system)
                    .expanded(1e-9 * prim.bounding_box(s.coord_system).size().norm());
                entries.push(Entry {
                    bbox,
                    hit: Hit {
                        property: prop,
                        primitive: prim,
                        primitive_index: i,
                    },
                });
            }
        }
        entries.reverse();
        entries.sort_by_key(|e| core::cmp::Reverse(e.hit.primitive.priority));
        SpatialIndex {
            entries,
            cs,
            default_cs: s.coord_system,
        }
    }

    /// Same semantics as [`ContinuousStructure::property_at`].
    pub fn query(&self, p: Vec3) -> Option<Hit<'a>> {
        let c = self.cs.to_cartesian(p);
        self.entries
            .iter()
            .filter(|e| e.bbox.contains(c))
            .find(|e| {
                e.hit
                    .primitive
                    .contains_in(c, CoordSystem::Cartesian, self.default_cs)
            })
            .map(|e| e.hit)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{BoxShape, Sphere};
    use crate::properties::Material;

    fn model() -> ContinuousStructure {
        let mut s = ContinuousStructure::new();
        let sub = s.add_property(Property::material("sub", Material::new(4.0, 0.0)));
        s.add_primitive(
            sub,
            Primitive::new(BoxShape::new([0.0, 0.0, 0.0], [10.0, 10.0, 1.0])).with_priority(1),
        )
        .unwrap();
        let m = s.add_property(Property::metal("pec"));
        s.add_primitive(
            m,
            Primitive::new(Sphere::new([5.0, 5.0, 1.0], 0.5)).with_priority(5),
        )
        .unwrap();
        s
    }

    #[test]
    fn priority() {
        let s = model();
        let h = s
            .property_at(
                Vec3::new(5.0, 5.0, 0.8),
                CoordSystem::Cartesian,
                TypeMask::ANY,
            )
            .unwrap();
        assert_eq!(h.property.name, "pec");
        let h = s
            .property_at(
                Vec3::new(1.0, 1.0, 0.5),
                CoordSystem::Cartesian,
                TypeMask::ANY,
            )
            .unwrap();
        assert_eq!(h.property.name, "sub");
        assert!(
            s.property_at(
                Vec3::new(1.0, 1.0, 2.0),
                CoordSystem::Cartesian,
                TypeMask::ANY
            )
            .is_none()
        );
        assert_eq!(
            s.properties_at(
                Vec3::new(5.0, 5.0, 0.8),
                CoordSystem::Cartesian,
                TypeMask::ANY
            )
            .len(),
            2
        );
        let idx = s.index(CoordSystem::Cartesian, TypeMask::ANY);
        for p in [
            Vec3::new(5.0, 5.0, 0.8),
            Vec3::new(1.0, 1.0, 0.5),
            Vec3::new(1.0, 1.0, 2.0),
        ] {
            assert_eq!(
                idx.query(p).map(|h| &h.property.name),
                s.property_at(p, CoordSystem::Cartesian, TypeMask::ANY)
                    .map(|h| &h.property.name)
            );
        }
    }

    #[test]
    fn media() {
        let s = model();
        assert_eq!(
            s.medium_at(Vec3::new(1.0, 1.0, 0.5), CoordSystem::Cartesian)
                .unwrap()
                .epsilon,
            [4.0; 3]
        );
        assert_eq!(
            s.medium_at(Vec3::new(1.0, 1.0, 5.0), CoordSystem::Cartesian)
                .unwrap()
                .epsilon,
            [1.0; 3]
        );
        assert!(
            s.medium_at(Vec3::new(5.0, 5.0, 1.0), CoordSystem::Cartesian)
                .is_none()
        );
        let b = s.bounding_box(CoordSystem::Cartesian);
        assert_eq!(b.max, Vec3::new(10.0, 10.0, 1.5));
    }

    #[test]
    fn parametric() {
        let mut s = ContinuousStructure::new();
        s.parameters.set("w", 2.0);
        let id = s.add_property(Property::metal("m"));
        s.add_primitive(
            id,
            BoxShape {
                start: [0.0, 0.0, 0.0].into(),
                stop: crate::ScalarPoint(["w".into(), "w*2".into(), 1.0.into()]),
            },
        )
        .unwrap();
        s.update().unwrap();
        assert!(
            s.property_at(
                Vec3::new(1.9, 3.9, 0.5),
                CoordSystem::Cartesian,
                TypeMask::ANY
            )
            .is_some()
        );
        s.parameters.set("w", 1.0);
        s.update().unwrap();
        assert!(
            s.property_at(
                Vec3::new(1.9, 3.9, 0.5),
                CoordSystem::Cartesian,
                TypeMask::ANY
            )
            .is_none()
        );
        assert!(s.validate().is_empty());
        assert!(s.remove_property(id).is_some());
    }
}
