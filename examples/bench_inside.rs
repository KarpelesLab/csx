//! Dependency-free micro benchmark of point-inside queries.
//!
//! Run with `cargo run --release --example bench_inside`.

use std::time::Instant;

use csx::primitives::{BoxShape, Polyhedron};
use csx::{ContinuousStructure, CoordSystem, Material, Primitive, Property, TypeMask, Vec3};

/// Small deterministic LCG so the benchmark needs no dependencies.
struct Lcg(u64);

impl Lcg {
    fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// UV sphere of radius 1 with `n` latitude bands and `n` longitude segments.
fn uv_sphere(n: u32) -> Polyhedron {
    let mut v = vec![Vec3::new(0.0, 0.0, 1.0)];
    for i in 1..n {
        let th = std::f64::consts::PI * i as f64 / n as f64;
        for j in 0..n {
            let ph = std::f64::consts::TAU * j as f64 / n as f64;
            v.push(Vec3::new(
                th.sin() * ph.cos(),
                th.sin() * ph.sin(),
                th.cos(),
            ));
        }
    }
    v.push(Vec3::new(0.0, 0.0, -1.0));
    let south = (v.len() - 1) as u32;
    let idx = |i: u32, j: u32| 1 + (i - 1) * n + (j % n);
    let mut f = Vec::new();
    for j in 0..n {
        f.push(vec![0, idx(1, j), idx(1, j + 1)]);
        f.push(vec![south, idx(n - 1, j + 1), idx(n - 1, j)]);
    }
    for i in 1..n - 1 {
        for j in 0..n {
            f.push(vec![idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)]);
            f.push(vec![idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)]);
        }
    }
    Polyhedron::new(v, f).expect("valid sphere")
}

fn main() {
    let mut rng = Lcg(42);

    let sphere = uv_sphere(100);
    println!(
        "polyhedron: {} triangles, closed = {}",
        sphere.triangles().len(),
        sphere.is_closed()
    );
    let pts: Vec<Vec3> = (0..2000)
        .map(|_| {
            Vec3::new(
                rng.next_f64() * 2.4 - 1.2,
                rng.next_f64() * 2.4 - 1.2,
                rng.next_f64() * 2.4 - 1.2,
            )
        })
        .collect();
    let t = Instant::now();
    let inside = pts.iter().filter(|&&p| sphere.contains(p)).count();
    let dt = t.elapsed();
    println!(
        "polyhedron contains: {} queries, {} inside, {:.2} µs/query",
        pts.len(),
        inside,
        dt.as_secs_f64() * 1e6 / pts.len() as f64
    );

    let mut s = ContinuousStructure::new();
    for k in 0..1000 {
        let id = s.add_property(Property::material(
            format!("m{k}"),
            Material::new(1.0 + k as f64, 0.0),
        ));
        let c = Vec3::new(
            rng.next_f64() * 100.0,
            rng.next_f64() * 100.0,
            rng.next_f64() * 100.0,
        );
        let d = Vec3::new(
            1.0 + rng.next_f64() * 4.0,
            1.0 + rng.next_f64() * 4.0,
            1.0 + rng.next_f64() * 4.0,
        );
        s.add_primitive(
            id,
            Primitive::new(BoxShape::new(c - d, c + d)).with_priority(k % 7),
        )
        .expect("property exists");
    }
    let qs: Vec<Vec3> = (0..100_000)
        .map(|_| {
            Vec3::new(
                rng.next_f64() * 100.0,
                rng.next_f64() * 100.0,
                rng.next_f64() * 100.0,
            )
        })
        .collect();

    let t = Instant::now();
    let idx = s.index(CoordSystem::Cartesian, TypeMask::ANY);
    let build = t.elapsed();
    let t = Instant::now();
    let hits = qs.iter().filter(|&&p| idx.query(p).is_some()).count();
    let dt = t.elapsed();
    println!(
        "spatial index (1000 boxes): build {:.2} ms, {} queries, {} hits, {:.3} µs/query",
        build.as_secs_f64() * 1e3,
        qs.len(),
        hits,
        dt.as_secs_f64() * 1e6 / qs.len() as f64
    );

    let t = Instant::now();
    let hits2 = qs[..10_000]
        .iter()
        .filter(|&&p| {
            s.property_at(p, CoordSystem::Cartesian, TypeMask::ANY)
                .is_some()
        })
        .count();
    let dt = t.elapsed();
    println!(
        "linear property_at: 10000 queries, {} hits, {:.3} µs/query",
        hits2,
        dt.as_secs_f64() * 1e6 / 10_000.0
    );
}
