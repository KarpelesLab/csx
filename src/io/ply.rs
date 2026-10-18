//! PLY mesh reader (ASCII, binary little and big endian).

use crate::error::{Error, Result};
use crate::math::Vec3;
use crate::primitives::Polyhedron;

#[derive(Debug, Clone, Copy, PartialEq)]
enum Format {
    Ascii,
    BinLe,
    BinBe,
}

#[derive(Debug, Clone, Copy)]
enum Ty {
    I8,
    U8,
    I16,
    U16,
    I32,
    U32,
    F32,
    F64,
}

impl Ty {
    fn parse(s: &str) -> Result<Ty> {
        Ok(match s {
            "char" | "int8" => Ty::I8,
            "uchar" | "uint8" => Ty::U8,
            "short" | "int16" => Ty::I16,
            "ushort" | "uint16" => Ty::U16,
            "int" | "int32" => Ty::I32,
            "uint" | "uint32" => Ty::U32,
            "float" | "float32" => Ty::F32,
            "double" | "float64" => Ty::F64,
            _ => return Err(Error::Format(format!("PLY: unknown type {s:?}"))),
        })
    }

    fn size(self) -> usize {
        match self {
            Ty::I8 | Ty::U8 => 1,
            Ty::I16 | Ty::U16 => 2,
            Ty::I32 | Ty::U32 | Ty::F32 => 4,
            Ty::F64 => 8,
        }
    }
}

#[derive(Debug)]
struct Prop {
    name: String,
    ty: Ty,
    /// Count type for list properties.
    list: Option<Ty>,
}

#[derive(Debug)]
struct Element {
    name: String,
    count: usize,
    props: Vec<Prop>,
}

struct Reader<'a> {
    fmt: Format,
    data: &'a [u8],
    pos: usize,
    tokens: core::str::SplitAsciiWhitespace<'a>,
}

impl Reader<'_> {
    fn read(&mut self, ty: Ty) -> Result<f64> {
        let eof = || Error::Format("PLY: unexpected end of data".into());
        if self.fmt == Format::Ascii {
            let t = self.tokens.next().ok_or_else(eof)?;
            return t
                .parse::<f64>()
                .map_err(|_| Error::Format(format!("PLY: bad number {t:?}")));
        }
        let n = ty.size();
        let b = self.data.get(self.pos..self.pos + n).ok_or_else(eof)?;
        self.pos += n;
        let mut a = [0u8; 8];
        a[..n].copy_from_slice(b);
        if self.fmt == Format::BinBe {
            a[..n].reverse();
        }
        Ok(match ty {
            Ty::I8 => a[0] as i8 as f64,
            Ty::U8 => a[0] as f64,
            Ty::I16 => i16::from_le_bytes([a[0], a[1]]) as f64,
            Ty::U16 => u16::from_le_bytes([a[0], a[1]]) as f64,
            Ty::I32 => i32::from_le_bytes([a[0], a[1], a[2], a[3]]) as f64,
            Ty::U32 => u32::from_le_bytes([a[0], a[1], a[2], a[3]]) as f64,
            Ty::F32 => f32::from_le_bytes([a[0], a[1], a[2], a[3]]) as f64,
            Ty::F64 => f64::from_le_bytes(a),
        })
    }
}

/// Parses PLY data into a polyhedron. Uses the `x`, `y`, `z` properties of
/// the `vertex` element and the `vertex_indices` (or `vertex_index`) list of
/// the `face` element; everything else is skipped.
///
/// ```
/// use csx::io::ply;
///
/// let data = "ply\nformat ascii 1.0\nelement vertex 3\nproperty float x\n\
///             property float y\nproperty float z\nelement face 1\n\
///             property list uchar int vertex_indices\nend_header\n\
///             0 0 0\n1 0 0\n0 1 0\n3 0 1 2\n";
/// let p = ply::parse(data.as_bytes())?;
/// assert_eq!(p.faces().len(), 1);
/// assert!(!p.is_closed()); // a single triangle is an open surface
/// # Ok::<(), csx::Error>(())
/// ```
///
/// # Errors
///
/// [`Error::Format`] for a missing or malformed header, unknown types,
/// counts that cannot fit in the data, or truncated / malformed bodies;
/// [`Error::Invalid`] if faces reference missing vertices.
pub fn parse(data: &[u8]) -> Result<Polyhedron> {
    let end =
        find_header_end(data).ok_or_else(|| Error::Format("PLY: missing end_header".into()))?;
    let header = core::str::from_utf8(&data[..end.0])
        .map_err(|_| Error::Format("PLY: bad header".into()))?;
    let mut lines = header.lines();
    if lines.next().map(str::trim) != Some("ply") {
        return Err(Error::Format("PLY: missing magic".into()));
    }
    let mut fmt = None;
    let mut elements: Vec<Element> = Vec::new();
    for line in lines {
        let w: Vec<&str> = line.split_ascii_whitespace().collect();
        match w.as_slice() {
            ["format", f, _] => {
                fmt = Some(match *f {
                    "ascii" => Format::Ascii,
                    "binary_little_endian" => Format::BinLe,
                    "binary_big_endian" => Format::BinBe,
                    _ => return Err(Error::Format(format!("PLY: unknown format {f:?}"))),
                })
            }
            ["element", name, count] => elements.push(Element {
                name: (*name).to_owned(),
                count: count
                    .parse()
                    .map_err(|_| Error::Format("PLY: bad element count".into()))?,
                props: Vec::new(),
            }),
            ["property", "list", ct, ty, name] => elements
                .last_mut()
                .ok_or_else(|| Error::Format("PLY: property before element".into()))?
                .props
                .push(Prop {
                    name: (*name).to_owned(),
                    ty: Ty::parse(ty)?,
                    list: Some(Ty::parse(ct)?),
                }),
            ["property", ty, name] => elements
                .last_mut()
                .ok_or_else(|| Error::Format("PLY: property before element".into()))?
                .props
                .push(Prop {
                    name: (*name).to_owned(),
                    ty: Ty::parse(ty)?,
                    list: None,
                }),
            _ => {}
        }
    }
    let fmt = fmt.ok_or_else(|| Error::Format("PLY: missing format line".into()))?;
    let body = &data[end.1..];
    let text = if fmt == Format::Ascii {
        core::str::from_utf8(body).map_err(|_| Error::Format("PLY: bad ASCII body".into()))?
    } else {
        ""
    };
    let mut r = Reader {
        fmt,
        data: body,
        pos: 0,
        tokens: text.split_ascii_whitespace(),
    };
    let mut verts = Vec::new();
    let mut faces = Vec::new();
    for el in &elements {
        if el.props.is_empty() {
            // Nothing to read; only meaningful for elements we ignore.
            if el.name == "vertex" && el.count > 0 {
                return Err(Error::Format(
                    "PLY: vertex element without properties".into(),
                ));
            }
            continue;
        }
        // Every instance consumes at least one byte, which bounds the
        // allocations driven by the (untrusted) header counts.
        if el.count > body.len() {
            return Err(Error::Format(format!(
                "PLY: element {:?} count {} exceeds the data size",
                el.name, el.count
            )));
        }
        if el.name == "vertex" {
            verts.reserve(el.count);
        }
        for _ in 0..el.count {
            let mut xyz = [0.0; 3];
            for p in &el.props {
                if let Some(ct) = p.list {
                    let n = r.read(ct)?;
                    if !(0.0..=1e7).contains(&n) {
                        return Err(Error::Format("PLY: bad list length".into()));
                    }
                    let mut idx = Vec::with_capacity((n as usize).min(body.len()));
                    for _ in 0..n as usize {
                        idx.push(r.read(p.ty)?);
                    }
                    if el.name == "face" && (p.name == "vertex_indices" || p.name == "vertex_index")
                    {
                        let f = idx
                            .into_iter()
                            .map(|v| {
                                if v >= 0.0 && v <= u32::MAX as f64 {
                                    Ok(v as u32)
                                } else {
                                    Err(Error::Format("PLY: bad index".into()))
                                }
                            })
                            .collect::<Result<Vec<u32>>>()?;
                        faces.push(f);
                    }
                } else {
                    let v = r.read(p.ty)?;
                    if el.name == "vertex" {
                        match p.name.as_str() {
                            "x" => xyz[0] = v,
                            "y" => xyz[1] = v,
                            "z" => xyz[2] = v,
                            _ => {}
                        }
                    }
                }
            }
            if el.name == "vertex" {
                verts.push(Vec3(xyz));
            }
        }
    }
    Polyhedron::new(verts, faces)
}

/// Returns (header length, body start).
fn find_header_end(data: &[u8]) -> Option<(usize, usize)> {
    let pat = b"end_header";
    let i = data.windows(pat.len()).position(|w| w == pat)?;
    let mut j = i + pat.len();
    if data.get(j) == Some(&b'\r') {
        j += 1;
    }
    if data.get(j) == Some(&b'\n') {
        j += 1;
    }
    Some((i, j))
}

#[cfg(test)]
mod tests {
    use super::*;

    const ASCII: &str = "ply
format ascii 1.0
comment tetrahedron
element vertex 4
property float x
property float y
property float z
element face 4
property list uchar int vertex_indices
end_header
0 0 0
1 0 0
0 1 0
0 0 1
3 0 2 1
3 0 1 3
3 0 3 2
3 1 2 3
";

    #[test]
    fn ascii() {
        let p = parse(ASCII.as_bytes()).unwrap();
        assert_eq!(p.vertices().len(), 4);
        assert!(p.is_closed());
        assert!(p.contains(Vec3::new(0.1, 0.1, 0.1)));
    }

    #[test]
    fn binary() {
        for be in [false, true] {
            let mut d = format!(
                "ply\nformat {} 1.0\nelement vertex 4\nproperty double x\nproperty double y\nproperty double z\nproperty uchar red\nelement face 4\nproperty list uchar uint vertex_indices\nend_header\n",
                if be { "binary_big_endian" } else { "binary_little_endian" }
            )
            .into_bytes();
            let verts = [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ];
            for v in verts {
                for c in v {
                    d.extend_from_slice(&if be {
                        f64::to_be_bytes(c)
                    } else {
                        f64::to_le_bytes(c)
                    });
                }
                d.push(255);
            }
            for f in [[0u32, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]] {
                d.push(3);
                for i in f {
                    d.extend_from_slice(&if be { i.to_be_bytes() } else { i.to_le_bytes() });
                }
            }
            let p = parse(&d).unwrap();
            assert!(p.is_closed());
            assert!(p.contains(Vec3::new(0.2, 0.2, 0.2)));
            assert!(parse(&d[..d.len() - 3]).is_err());
        }
    }
}
