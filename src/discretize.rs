//! Sampling of a structure onto its rectilinear grid (FDTD discretisation).
//!
//! [`discretize`] turns the continuous model into per-edge / per-face
//! arrays suitable for a Yee-scheme solver.
//!
//! # Conventions
//!
//! * **Nodes** are the intersections of the mesh lines; node `(i, j, k)` has
//!   the flat index `i + nx·(j + ny·k)` (see [`MaterialGrid::index`]). All
//!   per-component arrays have one entry per node.
//! * **E edges.** Component `d` at node `n` is the edge from `n` to the next
//!   node along `d`, located at the edge centre. It exists when
//!   `n[d] < n_lines[d] − 1`; entries of non-existing edges hold background
//!   values.
//! * **H faces.** Component `d` at node `n` is the face normal to `d` whose
//!   lowest corner is `n`, located at the face centre. It exists when both
//!   transverse indices are below their last line.
//! * **ε and κ** of an edge are the volume-weighted arithmetic mean over the
//!   (up to four) quarter cells around the edge: along the edge its full
//!   length, transversally half of each adjacent cell. Each quarter is
//!   sampled at the centres of an `N×N×N` sub-grid
//!   ([`DiscretizeOptions::subsampling`]).
//! * **μ** of a face is the volume-weighted *harmonic* mean over the two half
//!   cells on either side of the face (the normal flux density is continuous
//!   through the layers, which therefore act in series); **σ** (magnetic
//!   conductivity) uses the arithmetic mean.
//! * Samples falling into metal are left out of the averages (the edge is
//!   handled by the PEC mask instead); if every sample is metal, background
//!   values are used.
//! * In a **cylindrical grid** the three axes are `(ρ, α, z)`, material
//!   components apply to those axes, and volume weights use the Jacobian
//!   (`ρ·dα` is the length of an angular step).
//!
//! # Metal and sheets
//!
//! An edge is PEC when
//!
//! * its centre lies in a solid metal primitive (the highest-priority
//!   material or metal at that point is a metal),
//! * it lies in the plane of a planar metal sheet (zero-thickness box,
//!   polygon, ...) snapped to a mesh plane and its centre is on the sheet,
//! * it is part of the staircase approximation of a metal curve or wire
//!   (curves are walked through the nearest mesh nodes).
//!
//! Conducting sheets are detected like metal sheets but reported through
//! [`MaterialGrid::sheet`] / [`MaterialGrid::sheets`] instead of the PEC
//! mask. Sheets, lumped elements and line-like primitives that do not
//! coincide with mesh lines produce warnings in [`MaterialGrid::warnings`].

use crate::coords::CoordSystem;
use crate::error::{Error, Result};
use crate::grid::RectilinearGrid;
use crate::math::{BoundingBox, Vec3};
use crate::primitives::{Primitive, Shape};
use crate::properties::{
    ExcitationType, LumpedElementType, MaterialQuantity, PropertyKind, PropertyType, TypeMask,
};
use crate::structure::{ContinuousStructure, SpatialIndex};

/// Options of [`discretize`].
#[derive(Debug, Clone, PartialEq)]
pub struct DiscretizeOptions {
    /// Number of sub-samples per axis inside each quarter / half cell
    /// (`N`, total `N³` samples). At least 1.
    pub subsampling: usize,
    /// Worker threads (0 = available parallelism).
    pub threads: usize,
    /// Relative snapping tolerance: a sheet or line coordinate is on a mesh
    /// line when it is within `snap_tolerance × local cell size`.
    pub snap_tolerance: f64,
}

impl Default for DiscretizeOptions {
    fn default() -> Self {
        DiscretizeOptions {
            subsampling: 1,
            threads: 0,
            snap_tolerance: 1e-6,
        }
    }
}

/// A conducting sheet referenced from [`MaterialGrid::sheet`].
#[derive(Debug, Clone, PartialEq)]
pub struct SheetInfo {
    /// Property name.
    pub name: String,
    /// Conductivity (S/m).
    pub conductivity: f64,
    /// Thickness (m).
    pub thickness: f64,
}

/// A grid edge: component `dir` at flat node index `index`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Edge {
    /// Axis of the edge.
    pub dir: usize,
    /// Flat index of the start node.
    pub index: usize,
}

/// Excitation amplitudes of one excitation property.
#[derive(Debug, Clone, PartialEq)]
pub struct ExcitationField {
    /// Property name.
    pub name: String,
    /// Excitation number.
    pub number: u32,
    /// Excitation type (E or H field source).
    pub excitation_type: ExcitationType,
    /// Delay (s).
    pub delay: f64,
    /// Non-zero amplitudes: E edges for E types, H faces for H types.
    pub amplitudes: Vec<(Edge, f64)>,
}

/// One edge of a lumped element. The element's values are distributed over
/// its edges: with `s` edges in series (along the element direction) and `p`
/// in parallel, each edge gets `R·p/s`, `L·p/s` and `C·s/p`.
#[derive(Debug, Clone, PartialEq)]
pub struct LumpedEdge {
    /// Property name.
    pub name: String,
    /// The edge.
    pub edge: Edge,
    /// Resistance of this edge (NaN if unset).
    pub r: f64,
    /// Capacitance of this edge (NaN if unset).
    pub c: f64,
    /// Inductance of this edge (NaN if unset).
    pub l: f64,
    /// How R, L and C are connected, if given.
    pub le_type: Option<LumpedElementType>,
    /// Whether PEC caps are requested.
    pub caps: bool,
}

/// Result of [`discretize`]; see the module docs for the conventions.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialGrid {
    /// Number of mesh lines per axis.
    pub dims: [usize; 3],
    /// Coordinate system of the grid.
    pub coord_system: CoordSystem,
    /// Relative permittivity per E edge.
    pub epsilon: [Vec<f64>; 3],
    /// Electric conductivity per E edge.
    pub kappa: [Vec<f64>; 3],
    /// Relative permeability per H face.
    pub mue: [Vec<f64>; 3],
    /// Magnetic conductivity per H face.
    pub sigma: [Vec<f64>; 3],
    /// PEC flag per E edge.
    pub pec: [Vec<bool>; 3],
    /// Conducting sheet per E edge: 0 = none, `n` = `sheets[n - 1]`.
    pub sheet: [Vec<u32>; 3],
    /// Conducting sheets referenced by [`MaterialGrid::sheet`].
    pub sheets: Vec<SheetInfo>,
    /// Excitations.
    pub excitations: Vec<ExcitationField>,
    /// Lumped element edges.
    pub lumped: Vec<LumpedEdge>,
    /// Discretisation warnings (unsnapped sheets, empty primitives, ...).
    pub warnings: Vec<String>,
}

impl MaterialGrid {
    /// Total number of nodes.
    pub fn len(&self) -> usize {
        self.dims.iter().product()
    }

    /// True if the grid has no node.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Flat index of node `(i, j, k)`.
    #[inline]
    pub fn index(&self, n: [usize; 3]) -> usize {
        n[0] + self.dims[0] * (n[1] + self.dims[1] * n[2])
    }

    /// Node of a flat index.
    #[inline]
    pub fn node(&self, idx: usize) -> [usize; 3] {
        let i = idx % self.dims[0];
        let r = idx / self.dims[0];
        [i, r % self.dims[1], r / self.dims[1]]
    }

    /// True if E edge `dir` at `n` exists.
    pub fn edge_exists(&self, dir: usize, n: [usize; 3]) -> bool {
        (0..3).all(|a| n[a] < self.dims[a]) && n[dir] + 1 < self.dims[dir]
    }

    /// True if H face `dir` at `n` exists.
    pub fn face_exists(&self, dir: usize, n: [usize; 3]) -> bool {
        (0..3).all(|a| n[a] < self.dims[a])
            && (0..3)
                .filter(|&a| a != dir)
                .all(|a| n[a] + 1 < self.dims[a])
    }

    /// Number of PEC edges.
    pub fn pec_count(&self) -> usize {
        self.pec
            .iter()
            .map(|v| v.iter().filter(|&&b| b).count())
            .sum()
    }
}

/// Largest per-axis sub-sampling factor accepted by [`discretize`].
pub const MAX_SUBSAMPLING: usize = 64;

/// Largest number of grid nodes [`discretize`] accepts.
pub const MAX_NODES: usize = 1 << 32;

/// Discretises the materials, metals, sheets, excitations and lumped
/// elements of `s` onto `s.grid`. The structure should be up to date (see
/// [`ContinuousStructure::update`]).
pub fn discretize(s: &ContinuousStructure, opts: &DiscretizeOptions) -> Result<MaterialGrid> {
    let dims = s.grid.num_lines();
    if dims.contains(&0) {
        return Err(Error::Invalid(
            "discretize: the grid has an axis without lines".into(),
        ));
    }
    if opts.subsampling == 0 || opts.subsampling > MAX_SUBSAMPLING {
        return Err(Error::Invalid(format!(
            "discretize: subsampling must be in 1..={MAX_SUBSAMPLING}"
        )));
    }
    let cs = s.grid.coord_system;
    let ctx = Ctx {
        s,
        g: &s.grid,
        cs,
        dims,
        index: s.index(cs, TypeMask::MATERIALS | TypeMask::of(PropertyType::Metal)),
        bg: [
            s.background.epsilon.value(),
            s.background.mue.value(),
            s.background.kappa.value(),
            s.background.sigma.value(),
        ],
        n_sub: opts.subsampling,
        snap_tol: opts.snap_tolerance,
    };
    let threads = if opts.threads == 0 {
        std::thread::available_parallelism().map_or(1, |n| n.get())
    } else {
        opts.threads
    };
    let total = dims
        .iter()
        .try_fold(1usize, |acc, &n| acc.checked_mul(n))
        .filter(|&t| t <= MAX_NODES)
        .ok_or_else(|| {
            Error::Invalid(format!(
                "discretize: grid {dims:?} exceeds {MAX_NODES} nodes"
            ))
        })?;

    let mut out = MaterialGrid {
        dims,
        coord_system: cs,
        epsilon: Default::default(),
        kappa: Default::default(),
        mue: Default::default(),
        sigma: Default::default(),
        pec: core::array::from_fn(|_| vec![false; total]),
        sheet: core::array::from_fn(|_| vec![0; total]),
        sheets: Vec::new(),
        excitations: Vec::new(),
        lumped: Vec::new(),
        warnings: Vec::new(),
    };

    for d in 0..3 {
        let e = ctx.par_map(threads, |idx| ctx.edge_material(d, idx));
        out.epsilon[d] = e.iter().map(|v| v.0).collect();
        out.kappa[d] = e.iter().map(|v| v.1).collect();
        let h = ctx.par_map(threads, |idx| ctx.face_material(d, idx));
        out.mue[d] = h.iter().map(|v| v.0).collect();
        out.sigma[d] = h.iter().map(|v| v.1).collect();
        let pec = ctx.par_map(threads, |idx| ctx.edge_in_solid_metal(d, idx));
        out.pec[d] = pec;
    }

    let mut warnings = Vec::new();
    for prop in s.properties() {
        for (pi, prim) in prop.primitives.iter().enumerate() {
            let label = || {
                format!(
                    "property {:?} primitive #{pi} ({})",
                    prop.name,
                    prim.shape.type_name()
                )
            };
            match &prop.kind {
                PropertyKind::Metal => {
                    // Solid metal is handled by the point query above; add
                    // sheets and lines (and thin wires) explicitly.
                    let dim = prim.dimension();
                    let is_wire = matches!(prim.shape, Shape::Wire(_));
                    if dim < 3 || is_wire {
                        for e in ctx.select_edges(prim, &label, &mut warnings, is_wire) {
                            out.pec[e.dir][e.index] = true;
                        }
                    }
                }
                PropertyKind::ConductingSheet(cs_) => {
                    let edges = ctx.select_edges(prim, &label, &mut warnings, false);
                    if !edges.is_empty() {
                        out.sheets.push(SheetInfo {
                            name: prop.name.clone(),
                            conductivity: cs_.conductivity.value(),
                            thickness: cs_.thickness.value(),
                        });
                        let id = out.sheets.len() as u32;
                        for e in edges {
                            out.sheet[e.dir][e.index] = id;
                        }
                    }
                }
                PropertyKind::Excitation(ex) => {
                    let h_type = ex.excitation_type.is_magnetic();
                    if !h_type && !ex.excitation_type.is_electric() {
                        warnings.push(format!(
                            "{}: excitation type {} is not discretised",
                            label(),
                            ex.excitation_type
                        ));
                        continue;
                    }
                    let edges = if h_type {
                        ctx.select_faces(prim)
                    } else {
                        ctx.select_edges(prim, &label, &mut warnings, false)
                    };
                    let mut amplitudes = Vec::new();
                    for e in edges {
                        let p = if h_type {
                            ctx.face_center(e.dir, e.index)
                        } else {
                            ctx.edge_mid(e.dir, e.index)
                        };
                        let a = ex.amplitude_at(e.dir, cs.to_cartesian(p));
                        if a != 0.0 {
                            amplitudes.push((e, a));
                        }
                    }
                    if amplitudes.is_empty() {
                        warnings.push(format!("{}: excitation covers no grid edge", label()));
                    }
                    let entry = out.excitations.iter_mut().find(|f| f.name == prop.name);
                    match entry {
                        Some(f) => f.amplitudes.extend(amplitudes),
                        None => out.excitations.push(ExcitationField {
                            name: prop.name.clone(),
                            number: ex.number,
                            excitation_type: ex.excitation_type,
                            delay: ex.delay.value(),
                            amplitudes,
                        }),
                    }
                }
                PropertyKind::LumpedElement(le) => {
                    let dir = le.direction.index();
                    let edges: Vec<Edge> = ctx
                        .select_edges(prim, &label, &mut warnings, false)
                        .into_iter()
                        .filter(|e| e.dir == dir)
                        .collect();
                    if edges.is_empty() {
                        warnings.push(format!(
                            "{}: lumped element covers no edge along axis {dir}",
                            label()
                        ));
                        continue;
                    }
                    // Series edges: distinct positions along `dir`; parallel
                    // edges: distinct transverse positions.
                    let mut along: Vec<usize> =
                        edges.iter().map(|e| ctx.node(e.index)[dir]).collect();
                    along.sort_unstable();
                    along.dedup();
                    let mut across: Vec<[usize; 2]> = edges
                        .iter()
                        .map(|e| {
                            let n = ctx.node(e.index);
                            [n[(dir + 1) % 3], n[(dir + 2) % 3]]
                        })
                        .collect();
                    across.sort_unstable();
                    across.dedup();
                    let (ns, np) = (along.len() as f64, across.len() as f64);
                    for e in edges {
                        out.lumped.push(LumpedEdge {
                            name: prop.name.clone(),
                            edge: e,
                            r: le.r.value() * np / ns,
                            c: le.c.value() * ns / np,
                            l: le.l.value() * np / ns,
                            le_type: le.le_type,
                            caps: le.caps,
                        });
                    }
                }
                _ => {}
            }
        }
    }
    out.warnings = warnings;
    Ok(out)
}

struct Ctx<'a> {
    s: &'a ContinuousStructure,
    g: &'a RectilinearGrid,
    cs: CoordSystem,
    dims: [usize; 3],
    index: SpatialIndex<'a>,
    /// Background epsilon, mue, kappa, sigma.
    bg: [f64; 4],
    n_sub: usize,
    snap_tol: f64,
}

/// Tracks whether all averaged sample pairs are identical.
#[derive(Default)]
struct Uniform {
    first: Option<(f64, f64)>,
    mixed: bool,
}

impl Uniform {
    fn see(&mut self, a: f64, b: f64) {
        match self.first {
            None => self.first = Some((a, b)),
            Some(f) => self.mixed |= f != (a, b),
        }
    }

    fn value(&self) -> Option<(f64, f64)> {
        if self.mixed { None } else { self.first }
    }
}

/// Sampled material: `[epsilon, mue, kappa, sigma]` per axis.
type Sample = [[f64; 3]; 4];

impl Ctx<'_> {
    #[inline]
    fn node(&self, idx: usize) -> [usize; 3] {
        let i = idx % self.dims[0];
        let r = idx / self.dims[0];
        [i, r % self.dims[1], r / self.dims[1]]
    }

    #[inline]
    fn flat(&self, n: [usize; 3]) -> usize {
        n[0] + self.dims[0] * (n[1] + self.dims[1] * n[2])
    }

    #[inline]
    fn line(&self, d: usize, i: usize) -> f64 {
        self.g.lines(d)[i]
    }

    /// Runs `f` for every node index, splitting the work over threads.
    fn par_map<T: Send>(&self, threads: usize, f: impl Fn(usize) -> T + Sync) -> Vec<T> {
        let total: usize = self.dims.iter().product();
        let threads = threads.clamp(1, total.max(1));
        if threads == 1 || total < 4096 {
            return (0..total).map(f).collect();
        }
        let chunk = total.div_ceil(threads);
        let f = &f;
        std::thread::scope(|sc| {
            let handles: Vec<_> = (0..threads)
                .map(|t| {
                    let lo = (t * chunk).min(total);
                    let hi = ((t + 1) * chunk).min(total);
                    sc.spawn(move || (lo..hi).map(f).collect::<Vec<T>>())
                })
                .collect();
            handles
                .into_iter()
                .flat_map(|h| h.join().expect("discretize worker panicked"))
                .collect()
        })
    }

    /// Material at a grid-coordinate point; `None` inside metal.
    fn sample(&self, p: Vec3) -> Option<Sample> {
        match self.index.query(p) {
            None => Some(self.bg.map(|v| [v; 3])),
            Some(h) => match &h.property.kind {
                PropertyKind::Metal => None,
                k => {
                    let m = k.material()?;
                    let c = self.cs.to_cartesian(p);
                    Some([
                        m.values_at(MaterialQuantity::Epsilon, c),
                        m.values_at(MaterialQuantity::Mue, c),
                        m.values_at(MaterialQuantity::Kappa, c),
                        m.values_at(MaterialQuantity::Sigma, c),
                    ])
                }
            },
        }
    }

    /// Physical length of a coordinate step `h` along `dir` at point `p`.
    #[inline]
    fn plen(&self, dir: usize, h: f64, p: Vec3) -> f64 {
        if self.cs == CoordSystem::Cylindrical && dir == 1 {
            h * p[0].abs()
        } else {
            h
        }
    }

    /// Integrates the samples of the box `[lo, hi]` (grid coordinates):
    /// returns `Σw`, and `Σw·f(sample)` for the given accumulator.
    fn integrate(&self, lo: Vec3, hi: Vec3, acc: &mut impl FnMut(f64, &Sample)) -> f64 {
        let n = self.n_sub;
        let step = (hi - lo) / n as f64;
        let mut wsum = 0.0;
        for a in 0..n {
            for b in 0..n {
                for c in 0..n {
                    let f = Vec3::new(a as f64 + 0.5, b as f64 + 0.5, c as f64 + 0.5);
                    let p = Vec3::new(
                        lo[0] + step[0] * f[0],
                        lo[1] + step[1] * f[1],
                        lo[2] + step[2] * f[2],
                    );
                    let w = (0..3)
                        .map(|d| self.plen(d, step[d].abs(), p))
                        .product::<f64>();
                    if let Some(s) = self.sample(p) {
                        wsum += w;
                        acc(w, &s);
                    }
                }
            }
        }
        wsum
    }

    fn edge_mid(&self, d: usize, idx: usize) -> Vec3 {
        let n = self.node(idx);
        let mut p = Vec3::new(self.line(0, n[0]), self.line(1, n[1]), self.line(2, n[2]));
        if n[d] + 1 < self.dims[d] {
            p[d] = 0.5 * (self.line(d, n[d]) + self.line(d, n[d] + 1));
        }
        p
    }

    fn face_center(&self, d: usize, idx: usize) -> Vec3 {
        let n = self.node(idx);
        let mut p = Vec3::new(self.line(0, n[0]), self.line(1, n[1]), self.line(2, n[2]));
        for a in (0..3).filter(|&a| a != d) {
            if n[a] + 1 < self.dims[a] {
                p[a] = 0.5 * (self.line(a, n[a]) + self.line(a, n[a] + 1));
            }
        }
        p
    }

    /// (ε, κ) of E edge `d` at node `idx`.
    fn edge_material(&self, d: usize, idx: usize) -> (f64, f64) {
        let n = self.node(idx);
        let bg = (self.bg[0], self.bg[2]);
        if n[d] + 1 >= self.dims[d] {
            return bg;
        }
        let (t1, t2) = ((d + 1) % 3, (d + 2) % 3);
        let mut lo = Vec3::ZERO;
        let mut hi = Vec3::ZERO;
        lo[d] = self.line(d, n[d]);
        hi[d] = self.line(d, n[d] + 1);
        let (mut eps, mut kap, mut wsum) = (0.0, 0.0, 0.0);
        let mut uniform = Uniform::default();
        for s1 in [-1i64, 1] {
            let Some(h1) = self.half(t1, n[t1], s1) else {
                continue;
            };
            for s2 in [-1i64, 1] {
                let Some(h2) = self.half(t2, n[t2], s2) else {
                    continue;
                };
                lo[t1] = self.line(t1, n[t1]);
                hi[t1] = lo[t1] + h1;
                lo[t2] = self.line(t2, n[t2]);
                hi[t2] = lo[t2] + h2;
                wsum += self.integrate(lo, hi, &mut |w, s| {
                    eps += w * s[0][d];
                    kap += w * s[2][d];
                    uniform.see(s[0][d], s[2][d]);
                });
            }
        }
        if wsum > 0.0 {
            // Identical samples average to themselves exactly.
            return uniform.value().unwrap_or((eps / wsum, kap / wsum));
        }
        match self.sample(self.edge_mid(d, idx)) {
            Some(s) => (s[0][d], s[2][d]),
            None => bg,
        }
    }

    /// Signed half cell size next to line `i` on side `s` along `a`.
    fn half(&self, a: usize, i: usize, s: i64) -> Option<f64> {
        if s < 0 {
            (i >= 1).then(|| -(self.line(a, i) - self.line(a, i - 1)) / 2.0)
        } else {
            (i + 1 < self.dims[a]).then(|| (self.line(a, i + 1) - self.line(a, i)) / 2.0)
        }
    }

    /// (μ, σ) of H face `d` at node `idx`.
    fn face_material(&self, d: usize, idx: usize) -> (f64, f64) {
        let n = self.node(idx);
        let bg = (self.bg[1], self.bg[3]);
        let (t1, t2) = ((d + 1) % 3, (d + 2) % 3);
        if n[t1] + 1 >= self.dims[t1] || n[t2] + 1 >= self.dims[t2] {
            return bg;
        }
        let mut lo = Vec3::ZERO;
        let mut hi = Vec3::ZERO;
        for t in [t1, t2] {
            lo[t] = self.line(t, n[t]);
            hi[t] = self.line(t, n[t] + 1);
        }
        let (mut inv_mu, mut sig, mut wsum) = (0.0, 0.0, 0.0);
        let mut uniform = Uniform::default();
        for s in [-1i64, 1] {
            let Some(h) = self.half(d, n[d], s) else {
                continue;
            };
            lo[d] = self.line(d, n[d]);
            hi[d] = lo[d] + h;
            wsum += self.integrate(lo, hi, &mut |w, smp| {
                inv_mu += w / smp[1][d];
                sig += w * smp[3][d];
                uniform.see(smp[1][d], smp[3][d]);
            });
        }
        if wsum > 0.0 {
            return uniform.value().unwrap_or((wsum / inv_mu, sig / wsum));
        }
        match self.sample(self.face_center(d, idx)) {
            Some(s) => (s[1][d], s[3][d]),
            None => bg,
        }
    }

    fn edge_in_solid_metal(&self, d: usize, idx: usize) -> bool {
        let n = self.node(idx);
        if n[d] + 1 >= self.dims[d] {
            return false;
        }
        matches!(self.index.query(self.edge_mid(d, idx)), Some(h) if matches!(h.property.kind, PropertyKind::Metal))
    }

    /// Local cell size around line `i` along `a` (the smaller neighbour).
    fn local_cell(&self, a: usize, i: usize) -> f64 {
        let l = self.g.lines(a);
        let mut h = f64::INFINITY;
        if i >= 1 {
            h = h.min(l[i] - l[i - 1]);
        }
        if i + 1 < l.len() {
            h = h.min(l[i + 1] - l[i]);
        }
        if h.is_finite() { h } else { 1.0 }
    }

    /// Mesh line index that `v` lies on (within the snapping tolerance).
    fn snap(&self, a: usize, v: f64) -> Option<usize> {
        let i = self.g.nearest_line(a, v)?;
        let tol = self.snap_tol * self.local_cell(a, i) + 1e-12 * v.abs();
        ((self.line(a, i) - v).abs() <= tol).then_some(i)
    }

    /// Index range `[lo, hi]` of lines within `[min, max]` (with tolerance).
    fn line_range(&self, a: usize, min: f64, max: f64) -> Option<(usize, usize)> {
        let l = self.g.lines(a);
        let lo = match self.snap(a, min) {
            Some(i) => i,
            None => l.partition_point(|&x| x < min),
        };
        let hi = match self.snap(a, max) {
            Some(i) => i,
            None => l.partition_point(|&x| x <= max).checked_sub(1)?,
        };
        (lo < l.len() && lo <= hi).then_some((lo, hi))
    }

    /// Range of edge start indices along `a` whose edges overlap
    /// `[min, max]`.
    fn edge_range(&self, a: usize, min: f64, max: f64) -> Option<(usize, usize)> {
        let l = self.g.lines(a);
        if l.len() < 2 || max < l[0] || min > l[l.len() - 1] {
            return None;
        }
        let lo = l.partition_point(|&x| x <= min).saturating_sub(1);
        let hi = l
            .partition_point(|&x| x < max)
            .saturating_sub(1)
            .min(l.len() - 2);
        (lo <= hi).then_some((lo, hi))
    }

    fn bbox(&self, prim: &Primitive) -> BoundingBox {
        prim.bounding_box_in(self.cs, self.s.coord_system)
    }

    #[inline]
    fn contains(&self, prim: &Primitive, p: Vec3) -> bool {
        prim.contains_in(p, self.cs, self.s.coord_system)
    }

    /// Edges covered by a primitive, depending on its dimension (see module
    /// docs). `wire` additionally adds the staircase of a wire's centre line.
    fn select_edges(
        &self,
        prim: &Primitive,
        label: &dyn Fn() -> String,
        warnings: &mut Vec<String>,
        wire: bool,
    ) -> Vec<Edge> {
        let b = self.bbox(prim);
        if b.is_empty() {
            warnings.push(format!("{}: empty primitive", label()));
            return Vec::new();
        }
        let flat: Vec<usize> = (0..3)
            .filter(|&a| b.max[a] - b.min[a] <= self.flat_tol(a, b.min[a]))
            .collect();
        let dim = prim.dimension();
        let mut out = Vec::new();
        match (dim, flat.len()) {
            (0, _) => warnings.push(format!("{}: point primitives cover no edge", label())),
            (_, 3) => warnings.push(format!("{}: primitive is a single point", label())),
            (_, 2) => {
                // Axis aligned line along `d`.
                let d = (0..3).find(|a| !flat.contains(a)).expect("one free axis");
                let (a1, a2) = (flat[0], flat[1]);
                let (Some(i1), Some(i2)) = (self.snap(a1, b.min[a1]), self.snap(a2, b.min[a2]))
                else {
                    warnings.push(format!("{}: line not on mesh lines", label()));
                    return out;
                };
                if let Some((lo, hi)) = self.edge_range(d, b.min[d], b.max[d]) {
                    for i in lo..=hi {
                        let mut n = [0; 3];
                        n[d] = i;
                        n[a1] = i1;
                        n[a2] = i2;
                        let mut p = self.edge_mid(d, self.flat(n));
                        p[a1] = b.min[a1];
                        p[a2] = b.min[a2];
                        if self.contains(prim, p) {
                            out.push(Edge {
                                dir: d,
                                index: self.flat(n),
                            });
                        }
                    }
                }
            }
            _ if dim == 1 || wire => {
                out.extend(self.staircase(prim, label, warnings));
                if wire && dim == 3 {
                    out.extend(self.volume_edges(prim, &b));
                }
            }
            (2, 1) => {
                let nrm = flat[0];
                let plane = b.min[nrm];
                let Some(k) = self.snap(nrm, plane) else {
                    warnings.push(format!(
                        "{}: sheet at {plane} is not on a mesh plane (axis {nrm})",
                        label()
                    ));
                    return out;
                };
                for d in (0..3).filter(|&d| d != nrm) {
                    let t = 3 - d - nrm;
                    let (Some((lo, hi)), Some((tlo, thi))) = (
                        self.edge_range(d, b.min[d], b.max[d]),
                        self.line_range(t, b.min[t], b.max[t]),
                    ) else {
                        continue;
                    };
                    for i in lo..=hi {
                        for j in tlo..=thi {
                            let mut n = [0; 3];
                            n[d] = i;
                            n[t] = j;
                            n[nrm] = k;
                            let idx = self.flat(n);
                            let mut p = self.edge_mid(d, idx);
                            p[nrm] = plane;
                            if self.contains(prim, p) {
                                out.push(Edge { dir: d, index: idx });
                            }
                        }
                    }
                }
                if out.is_empty() {
                    warnings.push(format!("{}: sheet covers no grid edge", label()));
                }
            }
            (2, _) => {
                warnings.push(format!(
                    "{}: non axis-aligned sheet cannot be discretised",
                    label()
                ));
            }
            _ => {
                out = self.volume_edges(prim, &b);
                if out.is_empty() {
                    warnings.push(format!("{}: primitive covers no grid edge", label()));
                }
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    fn flat_tol(&self, a: usize, v: f64) -> f64 {
        let i = self.g.nearest_line(a, v).unwrap_or(0);
        self.snap_tol * self.local_cell(a, i) + 1e-12 * v.abs()
    }

    /// Edges whose centre lies inside a volumetric primitive.
    fn volume_edges(&self, prim: &Primitive, b: &BoundingBox) -> Vec<Edge> {
        let mut out = Vec::new();
        for d in 0..3 {
            let Some((lo, hi)) = self.edge_range(d, b.min[d], b.max[d]) else {
                continue;
            };
            let (t1, t2) = ((d + 1) % 3, (d + 2) % 3);
            let (Some((l1, h1)), Some((l2, h2))) = (
                self.line_range(t1, b.min[t1], b.max[t1]),
                self.line_range(t2, b.min[t2], b.max[t2]),
            ) else {
                continue;
            };
            for i in lo..=hi {
                for j in l1..=h1 {
                    for k in l2..=h2 {
                        let mut n = [0; 3];
                        n[d] = i;
                        n[t1] = j;
                        n[t2] = k;
                        let idx = self.flat(n);
                        if self.contains(prim, self.edge_mid(d, idx)) {
                            out.push(Edge { dir: d, index: idx });
                        }
                    }
                }
            }
        }
        out
    }

    /// Faces (H components) whose centre lies in the primitive.
    fn select_faces(&self, prim: &Primitive) -> Vec<Edge> {
        let b = self.bbox(prim);
        let mut out = Vec::new();
        if b.is_empty() {
            return out;
        }
        for d in 0..3 {
            let (t1, t2) = ((d + 1) % 3, (d + 2) % 3);
            let (Some((lo, hi)), Some((l1, h1)), Some((l2, h2))) = (
                self.line_range(d, b.min[d], b.max[d]),
                self.edge_range(t1, b.min[t1], b.max[t1]),
                self.edge_range(t2, b.min[t2], b.max[t2]),
            ) else {
                continue;
            };
            for i in lo..=hi {
                for j in l1..=h1 {
                    for k in l2..=h2 {
                        let mut n = [0; 3];
                        n[d] = i;
                        n[t1] = j;
                        n[t2] = k;
                        let idx = self.flat(n);
                        if self.contains(prim, self.face_center(d, idx)) {
                            out.push(Edge { dir: d, index: idx });
                        }
                    }
                }
            }
        }
        out
    }

    /// Staircase approximation of a curve or wire centre line through the
    /// nearest mesh nodes.
    fn staircase(
        &self,
        prim: &Primitive,
        label: &dyn Fn() -> String,
        warnings: &mut Vec<String>,
    ) -> Vec<Edge> {
        let pts: Vec<Vec3> = match &prim.shape {
            Shape::Curve(c) => c.points.iter().map(|p| p.value()).collect(),
            Shape::Wire(w) => w.points.iter().map(|p| p.value()).collect(),
            Shape::Box(bx) => vec![bx.start.value(), bx.stop.value()],
            _ => {
                warnings.push(format!(
                    "{}: one-dimensional shape cannot be discretised",
                    label()
                ));
                return Vec::new();
            }
        };
        // Points are Cartesian except for boxes, which use the primitive's
        // coordinate system.
        let pcs = if matches!(prim.shape, Shape::Box(_)) {
            prim.effective_cs(self.s.coord_system)
        } else {
            CoordSystem::Cartesian
        };
        let world: Vec<Vec3> = pts
            .iter()
            .map(|&p| {
                let c = pcs.to_cartesian(p);
                if prim.transform.is_empty() {
                    c
                } else {
                    prim.transform.apply(c)
                }
            })
            .collect();
        let mut out = Vec::new();
        let mut off_mesh = false;
        let to_node = |c: Vec3| -> Option<[usize; 3]> {
            let g = self.cs.from_cartesian(c);
            Some([
                self.g.nearest_line(0, g[0])?,
                self.g.nearest_line(1, g[1])?,
                self.g.nearest_line(2, g[2])?,
            ])
        };
        for w in world.windows(2) {
            let (Some(na), Some(nb)) = (to_node(w[0]), to_node(w[1])) else {
                continue;
            };
            let span: usize = (0..3).map(|a| na[a].abs_diff(nb[a])).sum();
            let steps = 4 * span.max(1) + 2;
            let mut cur = na;
            for s in 1..=steps {
                let t = s as f64 / steps as f64;
                let Some(next) = to_node(w[0] + (w[1] - w[0]) * t) else {
                    continue;
                };
                while cur != next {
                    let a = (0..3).find(|&a| cur[a] != next[a]).expect("differs");
                    let mut start = cur;
                    if next[a] > cur[a] {
                        cur[a] += 1;
                    } else {
                        cur[a] -= 1;
                        start = cur;
                    }
                    out.push(Edge {
                        dir: a,
                        index: self.flat(start),
                    });
                }
            }
            for p in [w[0], w[1]] {
                let g = self.cs.from_cartesian(p);
                if (0..3).any(|a| self.snap(a, g[a]).is_none()) {
                    off_mesh = true;
                }
            }
        }
        if off_mesh {
            warnings.push(format!(
                "{}: curve vertices are not on mesh nodes (staircased)",
                label()
            ));
        }
        if out.is_empty() {
            warnings.push(format!("{}: curve covers no grid edge", label()));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{BoxShape, Curve, Polygon};
    use crate::properties::{ConductingSheet, Excitation, LumpedElement, Material, Property};
    use crate::{Primitive, Scalar};

    fn grid(s: &mut ContinuousStructure, lines: [&[f64]; 3]) {
        for (d, l) in lines.iter().enumerate() {
            s.grid.set_lines(d, l.iter().copied()).unwrap();
        }
    }

    fn lines(n: usize) -> Vec<f64> {
        (0..=n).map(|i| i as f64).collect()
    }

    #[test]
    fn uniform_block() {
        let mut s = ContinuousStructure::new();
        let l = lines(4);
        grid(&mut s, [&l, &l, &l]);
        let id = s.add_property(Property::material("d", Material::new(4.0, 0.5)));
        s.add_primitive(id, BoxShape::new([-1.0, -1.0, -1.0], [5.0, 5.0, 5.0]))
            .unwrap();
        let g = discretize(&s, &DiscretizeOptions::default()).unwrap();
        for d in 0..3 {
            for idx in 0..g.len() {
                let n = g.node(idx);
                // non-existing edges hold background values
                let (e, k) = if g.edge_exists(d, n) {
                    (4.0, 0.5)
                } else {
                    (1.0, 0.0)
                };
                assert_eq!(g.epsilon[d][idx], e);
                assert_eq!(g.kappa[d][idx], k);
                assert_eq!(g.mue[d][idx], 1.0);
            }
        }
        assert_eq!(g.pec_count(), 0);
        assert!(g.warnings.is_empty(), "{:?}", g.warnings);
    }

    #[test]
    fn uniform_values_are_exact() {
        // Regression: weighted averaging of identical samples returned
        // e.g. 3.380000000000001.
        let mut s = ContinuousStructure::new();
        let l = [0.0, 0.1, 0.35, 1.0, 1.7];
        grid(&mut s, [&l, &l, &l]);
        let id = s.add_property(Property::material("d", Material::new(3.38, 0.0021)));
        s.add_primitive(id, BoxShape::new([-1.0, -1.0, -1.0], [5.0, 5.0, 5.0]))
            .unwrap();
        let opts = DiscretizeOptions {
            subsampling: 3,
            ..Default::default()
        };
        let g = discretize(&s, &opts).unwrap();
        for d in 0..3 {
            for idx in 0..g.len() {
                if g.edge_exists(d, g.node(idx)) {
                    assert_eq!(g.epsilon[d][idx], 3.38);
                    assert_eq!(g.kappa[d][idx], 0.0021);
                }
            }
        }
    }

    #[test]
    fn half_filled_averaging() {
        // Dielectric fills z < 1; the z = 1 plane is an interface.
        let mut s = ContinuousStructure::new();
        let l = lines(2);
        grid(&mut s, [&l, &l, &l]);
        let id = s.add_property(Property::material("d", Material::new(3.0, 0.0)));
        let mut m = Material::new(1.0, 0.0);
        m.mue.set(4.0);
        let id2 = s.add_property(Property::material("mu", m));
        s.add_primitive(id, BoxShape::new([-1.0, -1.0, -1.0], [3.0, 3.0, 1.0]))
            .unwrap();
        // permeable region x < 1 (priority above the dielectric)
        s.add_primitive(
            id2,
            Primitive::new(BoxShape::new([-1.0, -1.0, 1.5], [1.0, 3.0, 3.0])).with_priority(1),
        )
        .unwrap();
        let opts = DiscretizeOptions {
            subsampling: 2,
            threads: 2,
            ..Default::default()
        };
        let g = discretize(&s, &opts).unwrap();
        // E_x edge on the interface at (0,1,1): two quarters below (eps 3),
        // two above (eps 1).
        let e = g.epsilon[0][g.index([0, 1, 1])];
        assert!((e - 2.0).abs() < 1e-12, "{e}");
        // edge fully inside the dielectric
        assert_eq!(g.epsilon[0][g.index([0, 1, 0])], 3.0);
        // Harmonic mean along z: half cell z∈[1.5,2] has mu 4 only for
        // x<1, H_z face at node (0,0,2): face x∈[0,1], y∈[0,1]; halves
        // z∈[1.5,2] (mu=4) and z∈[2,2.5] (mu=4) → 4.
        assert!((g.mue[2][g.index([0, 0, 2])] - 4.0).abs() < 1e-12);
        // H_z at node (1,0,2): x∈[1,2] → mu 1.
        assert!((g.mue[2][g.index([1, 0, 2])] - 1.0).abs() < 1e-12);
        // H_x at x = 1, node (1,0,1): face y∈[0,1], z∈[1,2]; 2×2×2 samples
        // per half cell. Only the samples at z = 1.75 with x < 1
        // have mu = 4 (4 of 16) → harmonic mean 16 / (12 + 4/4).
        let h = g.mue[0][g.index([1, 0, 1])];
        assert!((h - 16.0 / 13.0).abs() < 1e-12, "{h}");
    }

    #[test]
    fn pec_sheet_on_mesh_plane() {
        let mut s = ContinuousStructure::new();
        let l = lines(4);
        grid(&mut s, [&l, &l, &l]);
        let id = s.add_property(Property::metal("m"));
        s.add_primitive(id, BoxShape::new([1.0, 1.0, 2.0], [3.0, 2.0, 2.0]))
            .unwrap();
        let g = discretize(&s, &DiscretizeOptions::default()).unwrap();
        // x edges: x∈{1,2} starts, y∈{1,2}: 4; y edges: x∈{1,2,3}, y start 1: 3
        assert_eq!(g.pec[0].iter().filter(|&&b| b).count(), 4);
        assert_eq!(g.pec[1].iter().filter(|&&b| b).count(), 3);
        assert_eq!(g.pec[2].iter().filter(|&&b| b).count(), 0);
        assert!(g.pec[0][g.index([1, 1, 2])]);
        assert!(g.warnings.is_empty(), "{:?}", g.warnings);

        // A sheet between mesh planes is reported.
        let mut s2 = s.clone();
        let id = s2.add_property(Property::metal("off"));
        s2.add_primitive(
            id,
            Polygon::new(2, 2.5, vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]),
        )
        .unwrap();
        let g = discretize(&s2, &DiscretizeOptions::default()).unwrap();
        assert!(
            g.warnings.iter().any(|w| w.contains("not on a mesh plane")),
            "{:?}",
            g.warnings
        );
    }

    #[test]
    fn solid_metal_and_conducting_sheet() {
        let mut s = ContinuousStructure::new();
        let l = lines(4);
        grid(&mut s, [&l, &l, &l]);
        let id = s.add_property(Property::metal("m"));
        s.add_primitive(id, BoxShape::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]))
            .unwrap();
        let cs = s.add_property(Property::new(
            "cu",
            PropertyKind::ConductingSheet(ConductingSheet {
                conductivity: 5.8e7.into(),
                thickness: 35e-6.into(),
            }),
        ));
        s.add_primitive(cs, BoxShape::new([2.0, 2.0, 3.0], [3.0, 3.0, 3.0]))
            .unwrap();
        let g = discretize(&s, &DiscretizeOptions::default()).unwrap();
        // the 12 edges of the unit cube (midpoints on its surface)
        assert_eq!(g.pec_count(), 12);
        assert_eq!(g.sheets.len(), 1);
        assert_eq!(g.sheet[0].iter().filter(|&&v| v == 1).count(), 2);
        assert_eq!(g.sheets[0].conductivity, 5.8e7);
    }

    #[test]
    fn wire_edges() {
        let mut s = ContinuousStructure::new();
        let l = lines(4);
        grid(&mut s, [&l, &l, &l]);
        let id = s.add_property(Property::metal("w"));
        // L shaped curve on mesh nodes plus a diagonal segment.
        s.add_primitive(
            id,
            Curve::new([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 2.0, 0.0]]),
        )
        .unwrap();
        let g = discretize(&s, &DiscretizeOptions::default()).unwrap();
        assert_eq!(g.pec[0].iter().filter(|&&b| b).count(), 2);
        assert_eq!(g.pec[1].iter().filter(|&&b| b).count(), 2);
        assert!(g.pec[1][g.index([2, 1, 0])]);
        assert!(g.warnings.is_empty(), "{:?}", g.warnings);

        let mut s = ContinuousStructure::new();
        grid(&mut s, [&l, &l, &l]);
        let id = s.add_property(Property::metal("d"));
        s.add_primitive(id, Curve::new([[0.0, 0.0, 0.0], [3.0, 3.0, 0.0]]))
            .unwrap();
        let g = discretize(&s, &DiscretizeOptions::default()).unwrap();
        // a connected staircase: 3 x steps and 3 y steps
        assert_eq!(g.pec[0].iter().filter(|&&b| b).count(), 3);
        assert_eq!(g.pec[1].iter().filter(|&&b| b).count(), 3);
    }

    #[test]
    fn cylindrical_grid() {
        let mut s = ContinuousStructure::new();
        s.coord_system = CoordSystem::Cylindrical;
        s.grid.coord_system = CoordSystem::Cylindrical;
        let r = [0.5, 1.0, 1.5, 2.0, 2.5];
        let a: Vec<f64> = (0..9).map(|i| -0.4 + 0.1 * i as f64).collect();
        let z = [0.0, 1.0, 2.0];
        grid(&mut s, [&r, &a, &z]);
        // Dielectric annulus 1 ≤ ρ ≤ 2 (cylindrical box covering all α).
        let id = s.add_property(Property::material("d", Material::new(2.0, 0.0)));
        s.add_primitive(
            id,
            Primitive::new(BoxShape::new([1.0, -3.2, -1.0], [2.0, 3.2, 3.0]))
                .with_coord_system(CoordSystem::Cylindrical),
        )
        .unwrap();
        // Metal sheet at ρ = 2.5? no: a cylindrical sheet at ρ = 1.5.
        let m = s.add_property(Property::metal("m"));
        s.add_primitive(
            m,
            Primitive::new(BoxShape::new([1.5, -0.2, 0.0], [1.5, 0.2, 2.0]))
                .with_coord_system(CoordSystem::Cylindrical),
        )
        .unwrap();
        let g = discretize(&s, &DiscretizeOptions::default()).unwrap();
        // E_alpha edge inside the annulus
        assert_eq!(g.epsilon[1][g.index([2, 3, 0])], 2.0);
        // E_rho edge from ρ=0.5 to 1.0 at the inner boundary? fully outside
        assert_eq!(g.epsilon[0][g.index([0, 3, 1])], 1.0);
        // E_z edge on ρ = 1 (half inside, weights ∝ ρ): between 1 and 2,
        // closer to 2 because the outer quarter has larger radius.
        let e = g.epsilon[2][g.index([1, 3, 0])];
        assert!(e > 1.5 && e < 2.0, "{e}");
        // sheet edges: α edges from -0.2 to 0.2 (4) × 3 z lines + z edges 5 α lines × 2
        assert_eq!(g.pec[1].iter().filter(|&&b| b).count(), 12);
        assert_eq!(g.pec[2].iter().filter(|&&b| b).count(), 10);
        assert!(g.warnings.is_empty(), "{:?}", g.warnings);
    }

    #[test]
    fn excitation_and_lumped() {
        let mut s = ContinuousStructure::new();
        let l = lines(4);
        grid(&mut s, [&l, &l, &l]);
        let ex = Excitation {
            excite: [Scalar::new(0.0), Scalar::new(0.0), Scalar::new(2.0)],
            weight: [Default::default(), Default::default(), "x".into()],
            ..Default::default()
        };
        let e = s.add_property(Property::new("port", PropertyKind::Excitation(ex)));
        s.add_primitive(e, BoxShape::new([1.0, 1.0, 1.0], [2.0, 1.0, 3.0]))
            .unwrap();
        let le = LumpedElement {
            direction: crate::Axis::Z,
            r: Scalar::new(50.0),
            ..Default::default()
        };
        let r = s.add_property(Property::new("R", PropertyKind::LumpedElement(le)));
        s.add_primitive(r, BoxShape::new([1.0, 1.0, 1.0], [2.0, 1.0, 3.0]))
            .unwrap();
        let g = discretize(&s, &DiscretizeOptions::default()).unwrap();
        let f = &g.excitations[0];
        // z edges at x∈{1,2}, y=1, z starts {1,2}: 4 edges, amplitude 2·x
        assert_eq!(f.amplitudes.len(), 4);
        for (edge, a) in &f.amplitudes {
            assert_eq!(edge.dir, 2);
            let n = g.node(edge.index);
            assert_eq!(*a, 2.0 * n[0] as f64);
        }
        // lumped: 2 in series × 2 in parallel → each edge R = 50
        assert_eq!(g.lumped.len(), 4);
        assert!(
            g.lumped
                .iter()
                .all(|e| (e.r - 50.0).abs() < 1e-12 && e.c.is_nan())
        );
        assert!(g.warnings.is_empty(), "{:?}", g.warnings);
    }

    #[test]
    fn threads_agree() {
        let mut s = ContinuousStructure::new();
        let l = lines(19);
        grid(&mut s, [&l, &l, &l]);
        let id = s.add_property(Property::material("d", Material::new(3.0, 0.1)));
        s.add_primitive(id, crate::primitives::Sphere::new([9.5, 9.5, 9.5], 6.3))
            .unwrap();
        let m = s.add_property(Property::metal("m"));
        s.add_primitive(m, BoxShape::new([2.0, 2.0, 2.0], [5.0, 4.0, 3.0]))
            .unwrap();
        let one = discretize(
            &s,
            &DiscretizeOptions {
                threads: 1,
                ..Default::default()
            },
        )
        .unwrap();
        let many = discretize(
            &s,
            &DiscretizeOptions {
                threads: 4,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(one, many);
        assert!(one.pec_count() > 0);
    }

    #[test]
    fn index_helpers_and_errors() {
        let mut s = ContinuousStructure::new();
        assert!(discretize(&s, &DiscretizeOptions::default()).is_err());
        let l = lines(2);
        grid(&mut s, [&l, &lines(3), &lines(4)]);
        let opts = DiscretizeOptions {
            subsampling: 0,
            ..Default::default()
        };
        assert!(discretize(&s, &opts).is_err());
        let opts = DiscretizeOptions {
            subsampling: 1000,
            ..Default::default()
        };
        assert!(discretize(&s, &opts).is_err());
        let g = discretize(&s, &DiscretizeOptions::default()).unwrap();
        assert_eq!(g.dims, [3, 4, 5]);
        // Grids too large to allocate are rejected instead of aborting.
        let mut big = s.clone();
        for d in 0..3 {
            big.grid.set_lines(d, (0..2048).map(|i| i as f64)).unwrap();
        }
        assert!(discretize(&big, &DiscretizeOptions::default()).is_err());
        assert_eq!(g.node(g.index([2, 3, 4])), [2, 3, 4]);
        assert!(g.edge_exists(0, [1, 3, 4]));
        assert!(!g.edge_exists(0, [2, 0, 0]));
        assert!(g.face_exists(2, [1, 2, 4]));
        assert!(!g.face_exists(2, [2, 2, 4]));
    }
}
