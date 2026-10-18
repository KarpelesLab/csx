//! Deterministic fuzz-style robustness tests: parsers must never panic on
//! arbitrary input and must not allocate based on untrusted counts.

use csx::expr::{Expr, MAX_DEPTH};
use csx::io::{ply, stl};

/// Small deterministic PRNG (64-bit LCG, high bits).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn bytes(&mut self, len: usize) -> Vec<u8> {
        (0..len).map(|_| self.next() as u8).collect()
    }
}

const EXPR_ALPHABET: &[&str] = &[
    "1", "2.5", ".5", "1e3", "1e", "x", "y", "pi", "e", "(", ")", ",", "+", "-", "*", "/", "%",
    "^", "!", "=", "==", "!=", "<", "<=", ">", ">=", "&", "|", " ", "sin(", "atan2(", "if(", "jn(",
    "yn(", "max(", "foo(", "#", "é",
];

#[test]
fn expr_random_tokens_never_panic() {
    let mut rng = Lcg(1);
    for _ in 0..20_000 {
        let n = rng.below(24);
        let s: String = (0..n)
            .map(|_| EXPR_ALPHABET[rng.below(EXPR_ALPHABET.len())])
            .collect();
        if let Ok(e) = Expr::parse(&s) {
            // Evaluation never panics either (values may be NaN/inf).
            let _ = e.eval_with(|_| Some(0.7));
            if let Ok(c) = e.compile(&["x", "y"], |_| Some(1.0)) {
                let _ = c.eval(&[0.3]);
                let _ = c.eval(&[0.3, -2.0]);
            }
            // Display / source round-trip.
            let again = Expr::parse(&e.to_string()).expect("re-parse of displayed expression");
            assert_eq!(again, e, "{s:?}");
            assert_eq!(e.source(), s);
        }
    }
}

#[test]
fn expr_random_bytes_never_panic() {
    let mut rng = Lcg(2);
    for _ in 0..20_000 {
        let len = rng.below(40);
        let b = rng.bytes(len);
        let s = String::from_utf8_lossy(&b);
        let _ = Expr::parse(&s);
        // ASCII-only variant exercises more of the grammar.
        let a: String = b.iter().map(|&c| (c % 95 + 32) as char).collect();
        if let Ok(e) = Expr::parse(&a) {
            let _ = e.eval_with(|_| Some(1.0));
        }
    }
}

#[test]
fn expr_deep_nesting_is_rejected_not_overflowed() {
    let n = 100_000;
    for s in [
        format!("{}1{}", "(".repeat(n), ")".repeat(n)),
        format!("{}1", "-".repeat(n)),
        format!("{}1", "!".repeat(n)),
        "1+".repeat(n) + "1",
        "2^".repeat(n) + "2",
        "2^-".repeat(n) + "2",
        "sin(".repeat(n) + "1" + &")".repeat(n),
    ] {
        assert!(Expr::parse(&s).is_err(), "{}", &s[..20]);
    }
    // Expressions close to the limit still parse and evaluate on a default
    // (2 MiB) test thread stack. Each parenthesis level costs two levels of
    // parser recursion; a left-associative sum adds one tree level per term.
    let d = MAX_DEPTH / 2 - 2;
    let s = format!("{}x{}", "(".repeat(d), ")".repeat(d));
    let e = Expr::parse(&s).unwrap();
    assert_eq!(e.eval_with(|_| Some(1.0)).unwrap(), 1.0);
    assert_eq!(e.compile(&["x"], |_| None).unwrap().eval(&[2.0]), 2.0);
    let d = MAX_DEPTH - 2;
    let s = "x+".repeat(d) + "1";
    let e = Expr::parse(&s).unwrap();
    assert_eq!(e.eval_with(|_| Some(1.0)).unwrap(), (d + 1) as f64);
    assert_eq!(
        e.compile(&["x"], |_| None).unwrap().eval(&[1.0]),
        (d + 1) as f64
    );
    let e = Expr::parse(&("1+".repeat(d) + "1")).unwrap();
    assert_eq!(
        e.compile(&[], |_| None).unwrap().as_constant(),
        Some((d + 1) as f64)
    );
}

#[test]
fn bessel_functions_in_expressions() {
    let v = |s: &str| Expr::parse(s).unwrap().eval_const().unwrap();
    assert!((v("j0(1)") - 0.765_197_686_557_966_5).abs() < 1e-12);
    assert!((v("jn(2, 1)") - 0.114_903_484_931_900_5).abs() < 1e-12);
    assert!((v("y1(2.5)") - 0.145_918_137_966_785_77).abs() < 1e-12);
    assert!(Expr::parse("jn(1)").is_err());
}

fn ply_header(format: &str, body_elements: &str) -> String {
    format!("ply\nformat {format} 1.0\n{body_elements}end_header\n")
}

#[test]
fn ply_untrusted_counts_are_bounded() {
    // Huge element counts with tiny bodies must fail fast.
    for fmt in ["ascii", "binary_little_endian", "binary_big_endian"] {
        let h = ply_header(
            fmt,
            "element vertex 18446744073709551615\nproperty float x\nproperty float y\nproperty float z\n",
        );
        assert!(ply::parse(h.as_bytes()).is_err());
        let h = ply_header(fmt, "element vertex 1000000000000\n");
        assert!(ply::parse(h.as_bytes()).is_err());
        let mut d = ply_header(
            fmt,
            "element face 1\nproperty list uint uint vertex_indices\n",
        )
        .into_bytes();
        d.extend_from_slice(if fmt == "ascii" {
            b"9999999 1 2 3"
        } else {
            &[0x7f, 0x96, 0x98, 0, 1, 0, 0, 0]
        });
        assert!(ply::parse(&d).is_err());
    }
    // Elements without properties that we ignore are fine.
    let ok = ply_header(
        "ascii",
        "element vertex 3\nproperty float x\nproperty float y\nproperty float z\nelement junk 1000000000\nelement face 1\nproperty list uchar int vertex_indices\n",
    ) + "0 0 0\n1 0 0\n0 1 0\n3 0 1 2\n";
    assert_eq!(ply::parse(ok.as_bytes()).unwrap().faces().len(), 1);
}

#[test]
fn ply_garbage_and_truncation_never_panic() {
    let good = ply_header(
        "binary_little_endian",
        "element vertex 4\nproperty double x\nproperty double y\nproperty double z\nelement face 4\nproperty list uchar uint vertex_indices\n",
    );
    let mut d = good.into_bytes();
    for v in [
        [0.0f64, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ] {
        for c in v {
            d.extend_from_slice(&c.to_le_bytes());
        }
    }
    for f in [[0u32, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]] {
        d.push(3);
        for i in f {
            d.extend_from_slice(&i.to_le_bytes());
        }
    }
    assert!(ply::parse(&d).is_ok());
    for cut in 0..d.len() {
        let _ = ply::parse(&d[..cut]);
    }
    let mut rng = Lcg(3);
    for _ in 0..3_000 {
        let mut m = d.clone();
        for _ in 0..1 + rng.below(4) {
            let i = rng.below(m.len());
            m[i] = rng.next() as u8;
        }
        let _ = ply::parse(&m);
        let len = rng.below(300);
        let _ = ply::parse(&rng.bytes(len));
        let mut h = b"ply\nformat ascii 1.0\n".to_vec();
        h.extend(rng.bytes(len));
        let _ = ply::parse(&h);
    }
}

#[test]
fn stl_garbage_and_truncation_never_panic() {
    let ascii = "solid t\nfacet normal 0 0 1\nouter loop\nvertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\nendloop\nendfacet\nendsolid t\n";
    let tris = stl::parse_triangles(ascii.as_bytes()).unwrap();
    let bin = stl::write_binary(&tris);
    for cut in 0..bin.len() {
        let _ = stl::parse(&bin[..cut]);
    }
    for cut in 0..ascii.len() {
        let _ = stl::parse(&ascii.as_bytes()[..cut]);
    }
    // A header claiming billions of triangles must not allocate them.
    let mut huge = vec![0u8; 84];
    huge[80..84].copy_from_slice(&u32::MAX.to_le_bytes());
    assert!(stl::parse(&huge).is_err());
    let mut rng = Lcg(4);
    for _ in 0..3_000 {
        let mut m = bin.clone();
        let i = rng.below(m.len());
        m[i] = rng.next() as u8;
        let _ = stl::parse(&m);
        let len = rng.below(400);
        let _ = stl::parse(&rng.bytes(len));
        let mut a = b"solid x\n".to_vec();
        a.extend((0..len).map(|_| b" vertex0123456789.e-nan\nendloop"[rng.below(30)]));
        let _ = stl::parse(&a);
    }
    // Non-finite coordinates are rejected rather than producing a bogus mesh.
    let nan = ascii.replace("vertex 1 0 0", "vertex nan 0 0");
    assert!(stl::parse(nan.as_bytes()).is_err());
}
