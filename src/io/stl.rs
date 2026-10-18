//! STL mesh reader (ASCII and binary).

use crate::error::{Error, Result};
use crate::math::Vec3;
use crate::primitives::Polyhedron;

/// Parses STL data (binary or ASCII) into a polyhedron, merging coincident
/// vertices.
pub fn parse(data: &[u8]) -> Result<Polyhedron> {
    Polyhedron::from_triangles(&parse_triangles(data)?)
}

/// Parses STL data into a list of triangles.
pub fn parse_triangles(data: &[u8]) -> Result<Vec<[Vec3; 3]>> {
    if is_binary(data) {
        parse_binary(data)
    } else {
        parse_ascii(data)
    }
}

fn is_binary(data: &[u8]) -> bool {
    if data.len() < 84 {
        return false;
    }
    let n = u32::from_le_bytes([data[80], data[81], data[82], data[83]]) as usize;
    let expected = n.checked_mul(50).and_then(|x| x.checked_add(84));
    if expected == Some(data.len()) {
        return true;
    }
    // Binary files whose header happens to start with "solid" still match
    // the size check above; everything else must look like ASCII.
    !data.trim_ascii_start().starts_with(b"solid")
}

fn parse_binary(data: &[u8]) -> Result<Vec<[Vec3; 3]>> {
    if data.len() < 84 {
        return Err(Error::Format("STL: truncated binary header".into()));
    }
    let n = u32::from_le_bytes([data[80], data[81], data[82], data[83]]) as usize;
    let body = &data[84..];
    if body.len() / 50 < n {
        return Err(Error::Format(format!(
            "STL: expected {n} triangles, file is truncated"
        )));
    }
    let f = |b: &[u8], o: usize| f32::from_le_bytes([b[o], b[o + 1], b[o + 2], b[o + 3]]) as f64;
    Ok(body
        .as_chunks::<50>()
        .0
        .iter()
        .take(n)
        .map(|rec| {
            core::array::from_fn(|k| {
                let o = 12 + 12 * k;
                Vec3::new(f(rec, o), f(rec, o + 4), f(rec, o + 8))
            })
        })
        .collect())
}

fn parse_ascii(data: &[u8]) -> Result<Vec<[Vec3; 3]>> {
    let text =
        core::str::from_utf8(data).map_err(|_| Error::Format("STL: invalid ASCII data".into()))?;
    let mut tris = Vec::new();
    let mut cur: Vec<Vec3> = Vec::with_capacity(3);
    let mut tok = text.split_ascii_whitespace();
    while let Some(t) = tok.next() {
        match t {
            "vertex" => {
                let mut c = [0.0; 3];
                for v in &mut c {
                    *v = tok
                        .next()
                        .and_then(|s| s.parse().ok())
                        .ok_or_else(|| Error::Format("STL: bad vertex".into()))?;
                }
                cur.push(Vec3(c));
            }
            "endloop" => {
                if cur.len() < 3 {
                    return Err(Error::Format(
                        "STL: facet with fewer than 3 vertices".into(),
                    ));
                }
                for k in 1..cur.len() - 1 {
                    tris.push([cur[0], cur[k], cur[k + 1]]);
                }
                cur.clear();
            }
            _ => {}
        }
    }
    Ok(tris)
}

/// Serialises triangles as binary STL.
pub fn write_binary(tris: &[[Vec3; 3]]) -> Vec<u8> {
    let mut out = vec![0u8; 80];
    out.extend_from_slice(&(tris.len() as u32).to_le_bytes());
    for t in tris {
        let n = (t[1] - t[0]).cross(t[2] - t[0]).normalized();
        for v in core::iter::once(n).chain(t.iter().copied()) {
            for c in v.0 {
                out.extend_from_slice(&(c as f32).to_le_bytes());
            }
        }
        out.extend_from_slice(&[0, 0]);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const TETRA: &str = "solid t
facet normal 0 0 -1
 outer loop
  vertex 0 0 0
  vertex 0 1 0
  vertex 1 0 0
 endloop
endfacet
facet normal 0 -1 0
 outer loop
  vertex 0 0 0
  vertex 1 0 0
  vertex 0 0 1
 endloop
endfacet
facet normal -1 0 0
 outer loop
  vertex 0 0 0
  vertex 0 0 1
  vertex 0 1 0
 endloop
endfacet
facet normal 1 1 1
 outer loop
  vertex 1 0 0
  vertex 0 1 0
  vertex 0 0 1
 endloop
endfacet
endsolid t
";

    #[test]
    fn ascii_and_binary() {
        let p = parse(TETRA.as_bytes()).unwrap();
        assert_eq!(p.vertices().len(), 4);
        assert!(p.is_closed());
        assert!(p.contains(Vec3::new(0.1, 0.1, 0.1)));
        assert!(!p.contains(Vec3::new(0.5, 0.5, 0.5)));
        let bin = write_binary(&parse_triangles(TETRA.as_bytes()).unwrap());
        let q = parse(&bin).unwrap();
        assert_eq!(q.faces().len(), 4);
        assert!(q.contains(Vec3::new(0.1, 0.1, 0.1)));
        assert!(parse(&bin[..100]).is_err());
    }
}
