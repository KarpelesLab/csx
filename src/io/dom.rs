//! A minimal XML document model with a parser and a pretty printer.
//!
//! Supports elements, attributes (single or double quoted), text, CDATA,
//! comments, processing instructions (skipped), `<!DOCTYPE>` (skipped) and
//! the predefined and numeric character entities. Whitespace-only text
//! between elements is dropped. Namespaces and DTD entities are not
//! interpreted.
//!
//! This module is always available (it has no dependencies): the data model
//! uses [`Element`] to keep unknown XML content of properties, primitives
//! and structures so it survives a read/write round-trip. The `xml` feature
//! adds the CSXCAD/openEMS reader and writer on top of it.

use core::fmt::Write as _;

use crate::error::{Error, Result};

/// Maximum element nesting depth accepted by the parser.
pub const MAX_DEPTH: usize = 256;

/// A node of the document tree.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Node {
    /// A child element.
    Element(Element),
    /// Character data (entities decoded).
    Text(String),
}

/// An XML element.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Element {
    /// Tag name.
    pub name: String,
    /// Attributes in document order.
    pub attrs: Vec<(String, String)>,
    /// Child nodes in document order.
    pub children: Vec<Node>,
}

impl Element {
    /// Empty element named `name`.
    pub fn new(name: impl Into<String>) -> Self {
        Element {
            name: name.into(),
            attrs: Vec::new(),
            children: Vec::new(),
        }
    }

    /// Value of the first attribute named `name`.
    pub fn attr(&self, name: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|(k, _)| k == name)
            .map(|(_, v)| v.as_str())
    }

    /// Sets (or replaces) an attribute.
    pub fn set_attr(&mut self, name: impl Into<String>, value: impl Into<String>) {
        let (name, value) = (name.into(), value.into());
        match self.attrs.iter_mut().find(|(k, _)| *k == name) {
            Some(a) => a.1 = value,
            None => self.attrs.push((name, value)),
        }
    }

    /// Builder variant of [`Element::set_attr`].
    pub fn with_attr(mut self, name: impl Into<String>, value: impl Into<String>) -> Self {
        self.set_attr(name, value);
        self
    }

    /// Appends a child element.
    pub fn push(&mut self, child: Element) {
        self.children.push(Node::Element(child));
    }

    /// Appends a text node.
    pub fn push_text(&mut self, text: impl Into<String>) {
        self.children.push(Node::Text(text.into()));
    }

    /// Iterates over child elements.
    pub fn elements(&self) -> impl Iterator<Item = &Element> {
        self.children.iter().filter_map(|n| match n {
            Node::Element(e) => Some(e),
            Node::Text(_) => None,
        })
    }

    /// First child element named `name`.
    pub fn child(&self, name: &str) -> Option<&Element> {
        self.elements().find(|e| e.name == name)
    }

    /// Child elements named `name`.
    pub fn children_named<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a Element> + 'a {
        self.elements().filter(move |e| e.name == name)
    }

    /// Concatenated text content of the direct text children.
    pub fn text(&self) -> String {
        let mut s = String::new();
        for n in &self.children {
            if let Node::Text(t) = n {
                s.push_str(t);
            }
        }
        s
    }

    /// Serialises the element (pretty printed, 4-space indentation).
    pub fn to_xml(&self) -> String {
        let mut out = String::new();
        self.write_into(&mut out, 0);
        out
    }

    /// Serialises the element preceded by an XML declaration.
    pub fn to_document(&self) -> String {
        let mut out =
            String::from("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>\n");
        self.write_into(&mut out, 0);
        out
    }

    fn write_into(&self, out: &mut String, depth: usize) {
        indent(out, depth);
        out.push('<');
        out.push_str(&self.name);
        for (k, v) in &self.attrs {
            let _ = write!(out, " {k}=\"");
            escape_into(out, v, true);
            out.push('"');
        }
        if self.children.is_empty() {
            out.push_str(" />\n");
            return;
        }
        out.push('>');
        let text_only = self.children.iter().all(|n| matches!(n, Node::Text(_)));
        if text_only {
            escape_into(out, &self.text(), false);
        } else {
            out.push('\n');
            for n in &self.children {
                match n {
                    Node::Element(e) => e.write_into(out, depth + 1),
                    Node::Text(t) => {
                        indent(out, depth + 1);
                        escape_into(out, t.trim(), false);
                        out.push('\n');
                    }
                }
            }
            indent(out, depth);
        }
        let _ = writeln!(out, "</{}>", self.name);
    }
}

fn indent(out: &mut String, depth: usize) {
    for _ in 0..depth {
        out.push_str("    ");
    }
}

fn escape_into(out: &mut String, s: &str, attr: bool) {
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' if attr => out.push_str("&quot;"),
            '\n' if attr => out.push_str("&#10;"),
            '\r' if attr => out.push_str("&#13;"),
            '\t' if attr => out.push_str("&#9;"),
            c => out.push(c),
        }
    }
}

/// True if `name` can be written as an element or attribute name (as
/// accepted by [`parse`]).
pub fn is_valid_name(name: &str) -> bool {
    !name.is_empty()
        && !name.bytes().any(|c| {
            c.is_ascii_whitespace()
                || matches!(
                    c,
                    b'/' | b'>' | b'=' | b'<' | b'"' | b'\'' | b'&' | b'!' | b'?'
                )
        })
}

/// Parses a document and returns its root element.
pub fn parse(src: &str) -> Result<Element> {
    let mut p = Parser {
        s: src.as_bytes(),
        src,
        pos: 0,
    };
    p.parse_document()
}

struct Parser<'a> {
    s: &'a [u8],
    src: &'a str,
    pos: usize,
}

impl Parser<'_> {
    fn err(&self, msg: impl core::fmt::Display) -> Error {
        let pos = self.pos.min(self.s.len());
        let before = &self.s[..pos];
        let line = before.iter().filter(|&&b| b == b'\n').count() + 1;
        let col = pos
            - before
                .iter()
                .rposition(|&b| b == b'\n')
                .map_or(0, |i| i + 1)
            + 1;
        Error::Format(format!("XML line {line}, column {col}: {msg}"))
    }

    fn starts_with(&self, pat: &str) -> bool {
        self.s[self.pos..].starts_with(pat.as_bytes())
    }

    fn skip_ws(&mut self) {
        while self.pos < self.s.len() && self.s[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
    }

    /// Skips past the next occurrence of `pat`.
    fn skip_past(&mut self, pat: &str, what: &str) -> Result<()> {
        match find(&self.s[self.pos..], pat.as_bytes()) {
            Some(i) => {
                self.pos += i + pat.len();
                Ok(())
            }
            None => Err(self.err(format!("unterminated {what}"))),
        }
    }

    /// Skips comments, processing instructions, DOCTYPE and whitespace
    /// outside the root element.
    fn skip_misc(&mut self) -> Result<()> {
        loop {
            self.skip_ws();
            if self.starts_with("<?") {
                self.skip_past("?>", "processing instruction")?;
            } else if self.starts_with("<!--") {
                self.skip_past("-->", "comment")?;
            } else if self.starts_with("<!DOCTYPE") || self.starts_with("<!doctype") {
                self.skip_doctype()?;
            } else {
                return Ok(());
            }
        }
    }

    fn skip_doctype(&mut self) -> Result<()> {
        let mut depth = 0i32;
        while self.pos < self.s.len() {
            match self.s[self.pos] {
                b'[' => depth += 1,
                b']' => depth -= 1,
                b'>' if depth <= 0 => {
                    self.pos += 1;
                    return Ok(());
                }
                _ => {}
            }
            self.pos += 1;
        }
        Err(self.err("unterminated DOCTYPE"))
    }

    fn parse_document(&mut self) -> Result<Element> {
        // Optional UTF-8 byte order mark.
        if self.starts_with("\u{feff}") {
            self.pos += 3;
        }
        self.skip_misc()?;
        if self.pos >= self.s.len() || self.s[self.pos] != b'<' {
            return Err(self.err("expected root element"));
        }
        let root = self.parse_element()?;
        self.skip_misc()?;
        if self.pos < self.s.len() {
            return Err(self.err("content after root element"));
        }
        Ok(root)
    }

    fn parse_name(&mut self) -> Result<String> {
        let start = self.pos;
        while self.pos < self.s.len() {
            let c = self.s[self.pos];
            if c.is_ascii_whitespace() || matches!(c, b'/' | b'>' | b'=' | b'<' | b'"' | b'\'') {
                break;
            }
            self.pos += 1;
        }
        if self.pos == start {
            return Err(self.err("expected a name"));
        }
        Ok(self.src[start..self.pos].to_owned())
    }

    /// Parses an element starting at `<`, using an explicit stack instead of
    /// recursion.
    fn parse_element(&mut self) -> Result<Element> {
        let mut stack: Vec<Element> = Vec::new();
        loop {
            // At a start tag '<'.
            self.pos += 1;
            let mut el = Element::new(self.parse_name()?);
            let self_closing = self.parse_attrs(&mut el)?;
            if self_closing {
                match stack.last_mut() {
                    Some(parent) => parent.children.push(Node::Element(el)),
                    None => return Ok(el),
                }
            } else {
                if stack.len() >= MAX_DEPTH {
                    return Err(self.err(format!("nesting deeper than {MAX_DEPTH}")));
                }
                stack.push(el);
            }
            // Content until the next child start tag.
            loop {
                let top = stack.last_mut().expect("stack is non-empty here");
                let text_start = self.pos;
                while self.pos < self.s.len() && self.s[self.pos] != b'<' {
                    self.pos += 1;
                }
                if self.pos > text_start {
                    let raw = &self.src[text_start..self.pos];
                    if !raw.trim().is_empty() {
                        push_text(top, decode_entities(raw));
                    }
                }
                if self.pos >= self.s.len() {
                    let name = top.name.clone();
                    return Err(self.err(format!("unclosed element <{name}>")));
                }
                if self.starts_with("</") {
                    self.pos += 2;
                    let name = self.parse_name()?;
                    self.skip_ws();
                    if self.s.get(self.pos) != Some(&b'>') {
                        return Err(self.err("expected '>'"));
                    }
                    self.pos += 1;
                    let done = stack.pop().expect("stack is non-empty here");
                    if done.name != name {
                        return Err(
                            self.err(format!("mismatched </{name}>, expected </{}>", done.name))
                        );
                    }
                    match stack.last_mut() {
                        Some(parent) => parent.children.push(Node::Element(done)),
                        None => return Ok(done),
                    }
                } else if self.starts_with("<!--") {
                    self.skip_past("-->", "comment")?;
                } else if self.starts_with("<![CDATA[") {
                    self.pos += 9;
                    let start = self.pos;
                    self.skip_past("]]>", "CDATA section")?;
                    let t = self.src[start..self.pos - 3].to_owned();
                    let top = stack.last_mut().expect("stack is non-empty here");
                    push_text(top, t);
                } else if self.starts_with("<?") {
                    self.skip_past("?>", "processing instruction")?;
                } else if self.starts_with("<!") {
                    return Err(self.err("unexpected markup declaration"));
                } else {
                    break;
                }
            }
        }
    }

    /// Parses attributes up to and including `>` or `/>`. Returns true for
    /// self-closing tags.
    fn parse_attrs(&mut self, el: &mut Element) -> Result<bool> {
        loop {
            self.skip_ws();
            match self.s.get(self.pos) {
                None => return Err(self.err(format!("unterminated tag <{}>", el.name))),
                Some(b'>') => {
                    self.pos += 1;
                    return Ok(false);
                }
                Some(b'/') => {
                    if self.s.get(self.pos + 1) == Some(&b'>') {
                        self.pos += 2;
                        return Ok(true);
                    }
                    return Err(self.err("expected '/>'"));
                }
                Some(_) => {
                    let name = self.parse_name()?;
                    self.skip_ws();
                    if self.s.get(self.pos) != Some(&b'=') {
                        return Err(self.err(format!("expected '=' after attribute {name}")));
                    }
                    self.pos += 1;
                    self.skip_ws();
                    let q = match self.s.get(self.pos) {
                        Some(&q @ (b'"' | b'\'')) => q,
                        _ => return Err(self.err("expected quoted attribute value")),
                    };
                    self.pos += 1;
                    let start = self.pos;
                    while self.pos < self.s.len() && self.s[self.pos] != q {
                        self.pos += 1;
                    }
                    if self.pos >= self.s.len() {
                        return Err(self.err("unterminated attribute value"));
                    }
                    let v = decode_entities(&self.src[start..self.pos]);
                    self.pos += 1;
                    el.attrs.push((name, v));
                }
            }
        }
    }
}

fn push_text(el: &mut Element, t: String) {
    if let Some(Node::Text(prev)) = el.children.last_mut() {
        prev.push_str(&t);
    } else {
        el.children.push(Node::Text(t));
    }
}

fn find(hay: &[u8], pat: &[u8]) -> Option<usize> {
    hay.windows(pat.len()).position(|w| w == pat)
}

/// Decodes predefined and numeric entities and normalises CR LF to LF.
/// Unknown or malformed entities are kept verbatim.
fn decode_entities(s: &str) -> String {
    if !s.contains('&') && !s.contains('\r') {
        return s.to_owned();
    }
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(i) = rest.find(['&', '\r']) {
        out.push_str(&rest[..i]);
        rest = &rest[i..];
        if let Some(r) = rest.strip_prefix('\r') {
            out.push('\n');
            rest = r.strip_prefix('\n').unwrap_or(r);
            continue;
        }
        let window = rest
            .char_indices()
            .take(16)
            .last()
            .map_or(0, |(i, c)| i + c.len_utf8());
        let Some(end) = rest[..window].find(';') else {
            out.push('&');
            rest = &rest[1..];
            continue;
        };
        let ent = &rest[1..end];
        let c = match ent {
            "lt" => Some('<'),
            "gt" => Some('>'),
            "amp" => Some('&'),
            "quot" => Some('"'),
            "apos" => Some('\''),
            _ => {
                let code =
                    if let Some(h) = ent.strip_prefix("#x").or_else(|| ent.strip_prefix("#X")) {
                        u32::from_str_radix(h, 16).ok()
                    } else if let Some(d) = ent.strip_prefix('#') {
                        d.parse::<u32>().ok()
                    } else {
                        None
                    };
                code.and_then(char::from_u32)
            }
        };
        match c {
            Some(c) => {
                out.push(c);
                rest = &rest[end + 1..];
            }
            None => {
                out.push('&');
                rest = &rest[1..];
            }
        }
    }
    out.push_str(rest);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_basic() {
        let doc = "\u{feff}<?xml version=\"1.0\"?>\r\n<!DOCTYPE x [<!ENTITY a \"b\">]>\n<!-- c -->\n<root a='1' b=\"x &amp; &lt;y&gt; &#65;&#x42;\">\n  <child/>\n  <t>hello <![CDATA[<raw>]]> &unknown; & world</t>\n  <!-- inner --><e k = \"v\" ></e>\n</root>\n";
        let r = parse(doc).unwrap();
        assert_eq!(r.name, "root");
        assert_eq!(r.attr("a"), Some("1"));
        assert_eq!(r.attr("b"), Some("x & <y> AB"));
        assert_eq!(r.elements().count(), 3);
        assert_eq!(
            r.child("t").unwrap().text(),
            "hello <raw> &unknown; & world"
        );
        assert_eq!(r.child("e").unwrap().attr("k"), Some("v"));
    }

    #[test]
    fn errors() {
        for bad in [
            "",
            "<a>",
            "<a></b>",
            "<a b></a>",
            "<a b=1></a>",
            "<a b=\"1></a>",
            "<a/><b/>",
            "text",
            "<a><!-- </a>",
            "<a><![CDATA[x</a>",
            "<",
            "<a",
        ] {
            assert!(parse(bad).is_err(), "{bad:?}");
        }
        let deep = "<a>".repeat(MAX_DEPTH + 1) + &"</a>".repeat(MAX_DEPTH + 1);
        assert!(parse(&deep).is_err());
        let ok = "<a>".repeat(MAX_DEPTH) + &"</a>".repeat(MAX_DEPTH);
        assert!(parse(&ok).is_ok());
        let e = parse("<a>\n<b c=x/></a>").unwrap_err().to_string();
        assert!(e.contains("line 2"), "{e}");
    }

    #[test]
    fn roundtrip() {
        let mut r = Element::new("r").with_attr("q", "a\"b<\n");
        let mut c = Element::new("c");
        c.push_text("1,2 & 3");
        r.push(c);
        r.push(Element::new("d"));
        let s = r.to_document();
        assert_eq!(parse(&s).unwrap(), r);
    }

    #[test]
    fn no_panic_on_truncation() {
        let seeds = [
            "<a b='1'>x&amp;</a>",
            "<?xml?><r><s/>é</r>",
            "<a>&#xFFFFFFFF;&#;</a>",
        ];
        for s in seeds {
            for cut in 0..=s.len() {
                if s.is_char_boundary(cut) {
                    let _ = parse(&s[..cut]);
                }
            }
        }
        assert_eq!(parse("<a>&#xFFFFFFFF;</a>").unwrap().text(), "&#xFFFFFFFF;");
        assert_eq!(
            parse("<a>&éééééééééééééééééé;</a>").unwrap().text(),
            "&éééééééééééééééééé;"
        );
    }
}
