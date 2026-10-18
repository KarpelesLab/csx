//! Property-style (deterministic random) tests of the primitives.

use core::f64::consts::{PI, TAU};

use csx::coords::{CoordSystem, cartesian_to_cylindrical};
use csx::primitives::{
    BoxShape, Curve, Cylinder, CylindricalShell, LinPoly, MultiBox, PointShape, Polygon,
    Polyhedron, Primitive, RotPoly, Sphere, SphericalShell, UserDefined, Wire,
};
use csx::{BoundingBox, Scalar, Transform, Vec3};

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 11
    }
    /// Uniform in [0, 1).
    fn unit(&mut self) -> f64 {
        self.next() as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.unit()
    }
    fn vec(&mut self, a: f64, b: f64) -> Vec3 {
        Vec3::new(self.range(a, b), self.range(a, b), self.range(a, b))
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// A random simple (star-shaped) polygon around `c`.
fn star_polygon(rng: &mut Lcg, c: (f64, f64), n: usize) -> Vec<(f64, f64)> {
    let mut angles: Vec<f64> = (0..n).map(|_| rng.range(0.0, TAU)).collect();
    angles.sort_by(f64::total_cmp);
    angles
        .iter()
        .map(|&a| {
            (
                c.0 + rng.range(0.3, 2.0) * a.cos(),
                c.1 + rng.range(0.3, 2.0) * a.sin(),
            )
        })
        .collect()
}

fn random_cube(rng: &mut Lcg) -> Polyhedron {
    let o = rng.vec(-2.0, 2.0);
    let s = rng.range(0.5, 2.0);
    let v = (0..8)
        .map(|i| o + Vec3::new((i & 1) as f64, ((i >> 1) & 1) as f64, ((i >> 2) & 1) as f64) * s)
        .collect();
    let faces = vec![
        vec![0, 2, 3, 1],
        vec![4, 5, 7, 6],
        vec![0, 1, 5, 4],
        vec![2, 6, 7, 3],
        vec![0, 4, 6, 2],
        vec![1, 3, 7, 5],
    ];
    Polyhedron::new(v, faces).unwrap()
}

fn random_shape(rng: &mut Lcg, kind: usize) -> Primitive {
    match kind {
        0 => Primitive::new(BoxShape::new(rng.vec(-3.0, 3.0), rng.vec(-3.0, 3.0))),
        1 => {
            // Cylindrical box: rho in [r1, r2], arbitrary (possibly wrapping)
            // angle range, z range.
            let r1 = rng.range(0.0, 2.0);
            let a1 = rng.range(-2.0 * PI, 2.0 * PI);
            let a2 = a1 + rng.range(0.1, 7.0);
            Primitive::new(BoxShape::new(
                [r1, a1, rng.range(-2.0, 0.0)],
                [r1 + rng.range(0.1, 2.0), a2, rng.range(0.0, 2.0)],
            ))
            .with_coord_system(CoordSystem::Cylindrical)
        }
        2 => Primitive::new(MultiBox {
            boxes: (0..3)
                .map(|_| (rng.vec(-3.0, 3.0).into(), rng.vec(-3.0, 3.0).into()))
                .collect(),
        }),
        3 => Primitive::new(Sphere::new(rng.vec(-2.0, 2.0), rng.range(0.1, 2.0))),
        4 => Primitive::new(SphericalShell {
            center: rng.vec(-2.0, 2.0).into(),
            radius: rng.range(0.5, 2.0).into(),
            shell_width: rng.range(0.1, 0.8).into(),
        }),
        5 => Primitive::new(Cylinder::new(
            rng.vec(-2.0, 2.0),
            rng.vec(-2.0, 2.0),
            rng.range(0.1, 1.5),
        )),
        6 => Primitive::new(CylindricalShell {
            start: rng.vec(-2.0, 2.0).into(),
            stop: rng.vec(-2.0, 2.0).into(),
            radius: rng.range(0.5, 1.5).into(),
            shell_width: rng.range(0.1, 0.5).into(),
        }),
        7 => {
            let n = 3 + rng.below(6);
            let c = (rng.range(-1.0, 1.0), rng.range(-1.0, 1.0));
            let poly = star_polygon(rng, c, n);
            Primitive::new(LinPoly::new(
                rng.below(3),
                rng.range(-1.0, 1.0),
                rng.range(-2.0, 2.0),
                poly,
            ))
        }
        8 => {
            let norm = rng.below(3);
            let mut axis = rng.below(3);
            if axis == norm {
                axis = (axis + 1) % 3;
            }
            let n = 3 + rng.below(6);
            let c = (rng.range(-1.0, 1.0), rng.range(-1.0, 1.0));
            let poly = star_polygon(rng, c, n);
            let mut r = RotPoly::new(norm, axis, poly);
            let a = rng.range(-PI, PI);
            r.start_angle = Scalar::new(a);
            r.stop_angle = Scalar::new(a + rng.range(0.2, TAU));
            Primitive::new(r)
        }
        9 => {
            let pts: Vec<Vec3> = (0..2 + rng.below(4)).map(|_| rng.vec(-2.0, 2.0)).collect();
            Primitive::new(Wire::new(pts, rng.range(0.05, 0.8)))
        }
        10 => Primitive::new(random_cube(rng)),
        _ => unreachable!(),
    }
}

const SOLID_KINDS: usize = 11;

fn random_transform(rng: &mut Lcg) -> Transform {
    let mut t = Transform::new();
    for _ in 0..1 + rng.below(3) {
        match rng.below(4) {
            0 => {
                t.translate(rng.vec(-3.0, 3.0));
            }
            1 => {
                t.rotate(
                    rng.vec(-1.0, 1.0) + Vec3::new(0.0, 0.0, 1e-3),
                    rng.range(-PI, PI),
                );
            }
            2 => {
                t.rotate_axis(rng.below(3), rng.range(-PI, PI));
            }
            _ => {
                t.scale(rng.range(0.5, 2.0));
            }
        }
    }
    t
}

fn contains_with_tol(b: &BoundingBox, p: Vec3) -> bool {
    b.expanded(1e-9 * (1.0 + b.size().norm())).contains(p)
}

/// Cylindrical bounds of a cylindrical box are reported in the box's own
/// (unwrapped) angle range, so compare modulo 2π.
fn angle_in_cyl_box(b: &BoundingBox, a: f64) -> bool {
    (-3..=3).any(|k| {
        let x = a + k as f64 * TAU;
        x >= b.min[1] - 1e-9 && x <= b.max[1] + 1e-9
    })
}

#[test]
fn inside_points_lie_in_bounding_boxes() {
    let mut rng = Lcg(11);
    let mut hits = 0usize;
    for iter in 0..1_200 {
        let kind = iter % SOLID_KINDS;
        let mut prim = random_shape(&mut rng, kind);
        if rng.below(2) == 0 {
            prim.transform = random_transform(&mut rng);
        }
        let bb = prim.bounding_box(CoordSystem::Cartesian);
        assert!(!bb.is_empty(), "kind {kind}");
        let cyl = prim.bounding_box_in(CoordSystem::Cylindrical, CoordSystem::Cartesian);
        let probe = bb.expanded(0.5);
        for _ in 0..300 {
            let p = Vec3::new(
                rng.range(probe.min[0], probe.max[0]),
                rng.range(probe.min[1], probe.max[1]),
                rng.range(probe.min[2], probe.max[2]),
            );
            if prim.contains(p, CoordSystem::Cartesian) {
                hits += 1;
                assert!(
                    contains_with_tol(&bb, p),
                    "kind {kind}: {p:?} outside {bb:?}"
                );
                let c = cartesian_to_cylindrical(p);
                assert!(
                    c[0] >= cyl.min[0] - 1e-9 && c[0] <= cyl.max[0] + 1e-9,
                    "kind {kind}: rho {} outside {cyl:?}",
                    c[0]
                );
                assert!(
                    angle_in_cyl_box(&cyl, c[1]),
                    "kind {kind}: alpha {} outside {cyl:?}",
                    c[1]
                );
                assert!(
                    c[2] >= cyl.min[2] - 1e-9 && c[2] <= cyl.max[2] + 1e-9,
                    "kind {kind}: z"
                );
            }
        }
    }
    assert!(
        hits > 10_000,
        "too few inside samples ({hits}); test is not exercising anything"
    );
}

#[test]
fn cartesian_and_cylindrical_queries_agree() {
    let mut rng = Lcg(12);
    for iter in 0..600 {
        let kind = iter % SOLID_KINDS;
        let mut prim = random_shape(&mut rng, kind);
        if rng.below(3) == 0 {
            prim.transform = random_transform(&mut rng);
        }
        let probe = prim.bounding_box(CoordSystem::Cartesian).expanded(0.5);
        for _ in 0..200 {
            let p = Vec3::new(
                rng.range(probe.min[0], probe.max[0]),
                rng.range(probe.min[1], probe.max[1]),
                rng.range(probe.min[2], probe.max[2]),
            );
            let c = cartesian_to_cylindrical(p);
            let a = prim.contains(p, CoordSystem::Cartesian);
            let b = prim.contains(c, CoordSystem::Cylindrical);
            // Same alpha shifted by a full turn is the same point.
            let d = prim.contains(Vec3::new(c[0], c[1] + TAU, c[2]), CoordSystem::Cylindrical);
            assert_eq!(a, b, "kind {kind} at {p:?}");
            assert_eq!(a, d, "kind {kind} at {p:?} (+2π)");
        }
    }
}

/// Re-expresses the defining points of a round shape in cylindrical
/// coordinates (identity for shapes that are always Cartesian).
fn to_cylindrical_points(prim: &Primitive) -> Primitive {
    use csx::Shape;
    let c =
        |p: &csx::ScalarPoint| -> csx::ScalarPoint { cartesian_to_cylindrical(p.value()).into() };
    let mut out = prim.clone();
    match &mut out.shape {
        Shape::Sphere(s) => s.center = c(&s.center),
        Shape::SphericalShell(s) => s.center = c(&s.center),
        Shape::Cylinder(s) => {
            s.start = c(&s.start);
            s.stop = c(&s.stop);
        }
        Shape::CylindricalShell(s) => {
            s.start = c(&s.start);
            s.stop = c(&s.stop);
        }
        Shape::Wire(w) => w.points = w.points.iter().map(c).collect(),
        _ => {}
    }
    out.with_coord_system(CoordSystem::Cylindrical)
}

#[test]
fn primitive_coord_system_semantics() {
    // As observed in CSXCAD: round shapes (spheres, cylinders, their shells,
    // wires) take their defining points in the primitive's coordinate
    // system; polygons and polyhedra are always Cartesian and multi-boxes
    // use the structure's coordinate system (Cartesian here).
    let mut rng = Lcg(13);
    for iter in 0..400 {
        let kind = 2 + iter % (SOLID_KINDS - 2);
        let prim = random_shape(&mut rng, kind);
        let cyl = to_cylindrical_points(&prim);
        let probe = prim.bounding_box(CoordSystem::Cartesian).expanded(0.3);
        let pb = cyl.bounding_box(CoordSystem::Cartesian).expanded(0.3);
        assert!(
            (probe.min - pb.min).norm() < 1e-9 && (probe.max - pb.max).norm() < 1e-9,
            "kind {kind}: {probe:?} vs {pb:?}"
        );
        for _ in 0..100 {
            let p = Vec3::new(
                rng.range(probe.min[0], probe.max[0]),
                rng.range(probe.min[1], probe.max[1]),
                rng.range(probe.min[2], probe.max[2]),
            );
            let a = prim.contains(p, CoordSystem::Cartesian);
            if a != cyl.contains(p, CoordSystem::Cartesian) {
                // Only points within rounding distance of the surface may
                // differ (the converted points are not bit-identical).
                let near = [1e-9, -1e-9]
                    .iter()
                    .any(|&e| prim.contains(p * (1.0 + e), CoordSystem::Cartesian) != a);
                assert!(near, "kind {kind} at {p:?}");
            }
        }
    }
}

#[test]
fn transformed_matches_inverse_mapped_query() {
    let mut rng = Lcg(14);
    for iter in 0..400 {
        let kind = iter % SOLID_KINDS;
        let base = random_shape(&mut rng, kind);
        let t = random_transform(&mut rng);
        let moved = base.clone().with_transform(t.clone());
        let probe = base.bounding_box(CoordSystem::Cartesian).expanded(0.3);
        for _ in 0..100 {
            let local = Vec3::new(
                rng.range(probe.min[0], probe.max[0]),
                rng.range(probe.min[1], probe.max[1]),
                rng.range(probe.min[2], probe.max[2]),
            );
            let world = t.apply(local);
            let back = t.apply_inverse(world);
            // Only compare when the round trip is exact enough not to cross
            // a boundary.
            if (back - local).norm() < 1e-12 {
                assert_eq!(
                    base.contains(local, CoordSystem::Cartesian),
                    moved.contains(world, CoordSystem::Cartesian),
                    "kind {kind}"
                );
            }
        }
    }
}

#[test]
fn rotpoly_partial_angles() {
    let mut rng = Lcg(15);
    for _ in 0..300 {
        // Square in the (radial, axial) half plane revolved about z.
        // norm_dir = 1 (y), rot axis = 2 (z): u = z, v = x → vertices (z, x).
        let (r0, r1) = (rng.range(0.2, 1.0), rng.range(1.2, 2.0));
        let (h0, h1) = (rng.range(-1.0, 0.0), rng.range(0.1, 1.0));
        let mut rp = RotPoly::new(1, 2, vec![(h0, r0), (h1, r0), (h1, r1), (h0, r1)]);
        let start = rng.range(-PI, PI);
        let stop = start + rng.range(0.1, TAU - 0.1);
        rp.start_angle = Scalar::new(start);
        rp.stop_angle = Scalar::new(stop);
        let prim = Primitive::new(rp);
        for _ in 0..100 {
            let phi = rng.range(-2.0 * PI, 2.0 * PI);
            let r = rng.range(r0 + 1e-6, r1 - 1e-6);
            let h = rng.range(h0 + 1e-6, h1 - 1e-6);
            let p = Vec3::new(r * phi.cos(), r * phi.sin(), h);
            let d = (phi - start).rem_euclid(TAU);
            let expected = d <= stop - start;
            // Avoid points within rounding distance of the end planes.
            let margin = 1e-9;
            if (d - (stop - start)).abs() > margin && d > margin && TAU - d > margin {
                assert_eq!(
                    prim.contains(p, CoordSystem::Cartesian),
                    expected,
                    "phi {phi} in [{start}, {stop}]"
                );
            }
        }
    }
}

#[test]
fn degenerate_inputs_do_not_panic() {
    let pts = [
        Vec3::ZERO,
        Vec3::new(1.0, 1.0, 1.0),
        Vec3::new(-0.0, 0.0, 0.0),
        Vec3::new(1e300, -1e300, 0.0),
    ];
    let mut prims: Vec<Primitive> = vec![
        Primitive::new(BoxShape::new([0.0; 3], [0.0; 3])),
        Primitive::new(BoxShape::new([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]))
            .with_coord_system(CoordSystem::Cylindrical),
        Primitive::new(BoxShape::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]))
            .with_coord_system(CoordSystem::Cylindrical),
        Primitive::new(BoxShape::new([1.0, 1e300, 0.0], [2.0, 1e300, 1.0]))
            .with_coord_system(CoordSystem::Cylindrical),
        Primitive::new(BoxShape::new([1.0, f64::NAN, 0.0], [2.0, 1.0, 1.0]))
            .with_coord_system(CoordSystem::Cylindrical),
        Primitive::new(MultiBox::default()),
        Primitive::new(Sphere::new([0.0; 3], 0.0)),
        Primitive::new(SphericalShell::default()),
        Primitive::new(Cylinder::new([0.0; 3], [0.0; 3], 1.0)),
        Primitive::new(Cylinder::new([0.0; 3], [0.0, 0.0, 1.0], 0.0)),
        Primitive::new(CylindricalShell::default()),
        Primitive::new(Polygon::new(2, 0.0, vec![])),
        Primitive::new(Polygon::new(2, 0.0, vec![(0.0, 0.0)])),
        Primitive::new(Polygon::new(
            2,
            0.0,
            vec![(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)],
        )),
        Primitive::new(LinPoly::new(0, 0.0, 0.0, vec![])),
        Primitive::new(LinPoly::new(1, 0.0, 1.0, vec![(1.0, 1.0), (1.0, 1.0)])),
        Primitive::new(RotPoly::new(0, 2, vec![])),
        Primitive::new(RotPoly::new(0, 2, vec![(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)])),
        Primitive::new(Curve::new(Vec::<Vec3>::new())),
        Primitive::new(Curve::new([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])),
        Primitive::new(Wire::new(Vec::<Vec3>::new(), 1.0)),
        Primitive::new(Wire::new([[0.0; 3], [0.0; 3]], 0.0)),
        Primitive::new(PointShape::new([0.0; 3])),
        Primitive::new(UserDefined::new("")),
        Primitive::new(UserDefined::new("1/0")),
        Primitive::new(UserDefined::new("undefined_param > 0")),
        Primitive::new(Polyhedron::default()),
    ];
    // Invalid directions on polygons must not panic either (update() reports
    // them as errors).
    let mut bad = Polygon::new(2, 0.0, vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]);
    bad.norm_dir = 7;
    prims.push(Primitive::new(bad));
    let mut bad = RotPoly::new(1, 1, vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]);
    bad.start_angle = Scalar::new(f64::NAN);
    prims.push(Primitive::new(bad));
    let params = csx::ParameterSet::new();
    for prim in &mut prims {
        for &p in &pts {
            let _ = prim.contains(p, CoordSystem::Cartesian);
            let _ = prim.contains(p, CoordSystem::Cylindrical);
            let _ = prim.contains(Vec3::new(f64::NAN, 0.0, 0.0), CoordSystem::Cartesian);
        }
        let _ = prim.bounding_box(CoordSystem::Cartesian);
        let _ = prim.bounding_box_in(CoordSystem::Cylindrical, CoordSystem::Cartesian);
        let _ = prim.dimension();
        let _ = prim.update(&params);
    }
    // Singular transformations are reported, not panicking.
    let mut t = Transform::new();
    t.scale(0.0);
    let mut p = Primitive::new(Sphere::new([0.0; 3], 1.0)).with_transform(t);
    assert!(p.update(&params).is_err());
    let _ = p.contains(Vec3::ZERO, CoordSystem::Cartesian);
}
