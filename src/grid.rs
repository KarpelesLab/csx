//! Rectilinear grids (mesh lines).
//!
//! A [`RectilinearGrid`] holds three sorted lists of mesh line positions.
//! Positions are in drawing units; multiply by [`RectilinearGrid::delta_unit`]
//! to obtain metres. In a cylindrical grid the three axes are `ρ`, `α`
//! (radians) and `z`.

use crate::coords::CoordSystem;
use crate::error::{Error, Result};
use crate::math::{BoundingBox, Vec3};

/// Three sets of mesh lines.
#[derive(Debug, Clone, PartialEq)]
pub struct RectilinearGrid {
    lines: [Vec<f64>; 3],
    /// Size of one drawing unit in metres (e.g. `1e-3` for millimetres).
    pub delta_unit: f64,
    /// Coordinate system of the grid.
    pub coord_system: CoordSystem,
}

impl Default for RectilinearGrid {
    fn default() -> Self {
        RectilinearGrid {
            lines: Default::default(),
            delta_unit: 1.0,
            coord_system: CoordSystem::Cartesian,
        }
    }
}

impl RectilinearGrid {
    /// Empty Cartesian grid with `delta_unit = 1`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Mesh lines along `dir`, sorted ascending without duplicates.
    pub fn lines(&self, dir: usize) -> &[f64] {
        &self.lines[dir]
    }

    /// Number of lines along each axis.
    pub fn num_lines(&self) -> [usize; 3] {
        [
            self.lines[0].len(),
            self.lines[1].len(),
            self.lines[2].len(),
        ]
    }

    /// Total number of cells (product of line count minus one per axis).
    pub fn num_cells(&self) -> usize {
        self.num_lines()
            .iter()
            .map(|&n| n.saturating_sub(1))
            .product()
    }

    /// Replaces the lines along `dir` (sorted and de-duplicated).
    pub fn set_lines(&mut self, dir: usize, lines: impl IntoIterator<Item = f64>) -> Result<()> {
        let v: Vec<f64> = lines.into_iter().collect();
        if let Some(bad) = v.iter().find(|x| !x.is_finite()) {
            return Err(Error::Invalid(format!("non-finite mesh line {bad}")));
        }
        self.lines[dir] = v;
        self.normalize(dir);
        Ok(())
    }

    /// Adds lines along `dir`.
    pub fn add_lines(&mut self, dir: usize, lines: impl IntoIterator<Item = f64>) -> Result<()> {
        for l in lines {
            if !l.is_finite() {
                return Err(Error::Invalid(format!("non-finite mesh line {l}")));
            }
            self.lines[dir].push(l);
        }
        self.normalize(dir);
        Ok(())
    }

    /// Adds a single line.
    pub fn add_line(&mut self, dir: usize, l: f64) -> Result<()> {
        self.add_lines(dir, [l])
    }

    /// Removes all lines along `dir`.
    pub fn clear_lines(&mut self, dir: usize) {
        self.lines[dir].clear();
    }

    /// Removes all lines.
    pub fn clear(&mut self) {
        self.lines.iter_mut().for_each(Vec::clear);
    }

    fn normalize(&mut self, dir: usize) {
        let v = &mut self.lines[dir];
        v.sort_by(f64::total_cmp);
        let scale = v.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        let eps = scale * 1e-13;
        v.dedup_by(|b, a| (*b - *a).abs() <= eps);
    }

    /// Bounds of the mesh (first and last line per axis).
    pub fn bounds(&self) -> BoundingBox {
        if self.lines.iter().any(Vec::is_empty) {
            return BoundingBox::EMPTY;
        }
        BoundingBox {
            min: Vec3::new(self.lines[0][0], self.lines[1][0], self.lines[2][0]),
            max: Vec3::new(
                *self.lines[0].last().unwrap(),
                *self.lines[1].last().unwrap(),
                *self.lines[2].last().unwrap(),
            ),
        }
    }

    /// Index of the line nearest to `v` along `dir`.
    pub fn nearest_line(&self, dir: usize, v: f64) -> Option<usize> {
        let l = &self.lines[dir];
        if l.is_empty() {
            return None;
        }
        let i = l.partition_point(|&x| x < v);
        if i == 0 {
            return Some(0);
        }
        if i == l.len() {
            return Some(l.len() - 1);
        }
        Some(if v - l[i - 1] <= l[i] - v { i - 1 } else { i })
    }

    /// Index of the line exactly at `v` (within a relative tolerance).
    pub fn snap(&self, dir: usize, v: f64) -> Option<usize> {
        let i = self.nearest_line(dir, v)?;
        let x = self.lines[dir][i];
        let tol = 1e-12
            * x.abs()
                .max(v.abs())
                .max(self.delta_cell(dir, i).unwrap_or(0.0));
        ((x - v).abs() <= tol).then_some(i)
    }

    /// Index `i` of the cell `[line[i], line[i+1])` containing `v`; the last
    /// line belongs to the last cell.
    pub fn cell_index(&self, dir: usize, v: f64) -> Option<usize> {
        let l = &self.lines[dir];
        if l.len() < 2 || v < l[0] || v > l[l.len() - 1] {
            return None;
        }
        let i = l.partition_point(|&x| x <= v);
        Some(i.saturating_sub(1).min(l.len() - 2))
    }

    /// Size of cell `i` along `dir`.
    pub fn delta_cell(&self, dir: usize, i: usize) -> Option<f64> {
        let l = &self.lines[dir];
        Some(l.get(i + 1)? - l.get(i)?)
    }

    /// Mesh node position `(lines[0][i], lines[1][j], lines[2][k])`.
    pub fn node(&self, i: usize, j: usize, k: usize) -> Option<Vec3> {
        Some(Vec3::new(
            *self.lines[0].get(i)?,
            *self.lines[1].get(j)?,
            *self.lines[2].get(k)?,
        ))
    }

    /// Smallest and largest cell size along `dir`.
    pub fn cell_size_range(&self, dir: usize) -> Option<(f64, f64)> {
        let l = &self.lines[dir];
        if l.len() < 2 {
            return None;
        }
        Some(
            l.windows(2)
                .map(|w| w[1] - w[0])
                .fold((f64::INFINITY, 0.0f64), |(a, b), d| (a.min(d), b.max(d))),
        )
    }

    /// Smooths the lines along `dir` in place; see [`smooth_lines`].
    pub fn smooth(&mut self, dir: usize, max_res: f64, ratio: f64) -> Result<()> {
        let s = smooth_lines(&self.lines[dir], max_res, ratio)?;
        self.lines[dir] = s;
        Ok(())
    }
}

/// Inserts lines between the given fixed lines so that no cell is larger
/// than `max_res` and neighbouring cells differ in size by roughly at most
/// `ratio` (≥ 1). Fixed lines are always kept.
///
/// The algorithm assigns every gap between fixed lines a target size
/// function that grows linearly (i.e. geometrically per cell) away from the
/// neighbouring cells, capped at `max_res`, and distributes points so each
/// cell spans one unit of `∫ dx / h(x)`. Neighbour sizes are refined
/// iteratively until stable.
pub fn smooth_lines(fixed: &[f64], max_res: f64, ratio: f64) -> Result<Vec<f64>> {
    if max_res.is_nan() || !ratio.is_finite() || max_res <= 0.0 || ratio < 1.0 {
        return Err(Error::Invalid(format!(
            "smooth_lines: max_res {max_res}, ratio {ratio}"
        )));
    }
    let mut f: Vec<f64> = fixed.iter().copied().filter(|x| x.is_finite()).collect();
    f.sort_by(f64::total_cmp);
    f.dedup();
    if f.len() < 2 {
        return Ok(f);
    }
    let n = f.len() - 1;
    // Edge cell sizes of every gap: (first cell, last cell).
    let mut gaps: Vec<Vec<f64>> = vec![Vec::new(); n];
    let mut edge: Vec<(f64, f64)> = vec![(max_res, max_res); n];
    for _ in 0..50 {
        let mut changed = false;
        for g in 0..n {
            let hl = if g > 0 { edge[g - 1].1 } else { max_res };
            let hr = if g + 1 < n { edge[g + 1].0 } else { max_res };
            let pts = fill_gap(f[g], f[g + 1], hl, hr, max_res, ratio)?;
            let first = pts.first().copied().unwrap_or(f[g + 1]) - f[g];
            let last = f[g + 1] - pts.last().copied().unwrap_or(f[g]);
            let new = (first, last);
            if (new.0 - edge[g].0).abs() > 1e-12 * new.0
                || (new.1 - edge[g].1).abs() > 1e-12 * new.1
            {
                changed = true;
            }
            edge[g] = new;
            gaps[g] = pts;
        }
        if !changed {
            break;
        }
    }
    let count = f.len() + gaps.iter().map(Vec::len).sum::<usize>();
    if count > MAX_SMOOTH_LINES {
        return Err(Error::Invalid(format!(
            "smooth_lines: more than {MAX_SMOOTH_LINES} lines required"
        )));
    }
    let mut out = Vec::with_capacity(count);
    for g in 0..n {
        out.push(f[g]);
        out.extend_from_slice(&gaps[g]);
    }
    out.push(f[n]);
    Ok(out)
}

/// Maximum number of lines [`smooth_lines`] may produce.
pub const MAX_SMOOTH_LINES: usize = 10_000_000;

/// Interior points for the gap `[a, b]` given neighbour cell sizes.
///
/// The target size function is `h(x) = min(hmax, A + k(x−a), B + k(b−x))`,
/// piecewise linear; `∫dx/h` is integrated and inverted exactly on each
/// linear piece so that cells follow the grading even when `h` varies by
/// orders of magnitude across the gap.
fn fill_gap(a: f64, b: f64, hl: f64, hr: f64, hmax: f64, r: f64) -> Result<Vec<f64>> {
    let len = b - a;
    if len <= hmax && len <= hl * r && len <= hr * r {
        return Ok(Vec::new());
    }
    // With h(x) = h0 + k·x and one cell per unit of ∫dx/h, consecutive cells
    // grow by exactly e^k; choose k = ln r and h0 so the first cell is h·r.
    let k = r.ln();
    let c = if k > 1e-12 { k / (r - 1.0) } else { 1.0 };
    let (ha, hb) = (c * hl * r, c * hr * r);
    let h = |x: f64| hmax.min(ha + k * (x - a)).min(hb + k * (b - x));
    // Breakpoints of the piecewise linear h.
    let mut xs = vec![a, b];
    if k > 0.0 {
        for x in [
            a + (hmax - ha) / k,
            b - (hmax - hb) / k,
            0.5 * (a + b + (hb - ha) / k),
        ] {
            if x > a && x < b {
                xs.push(x);
            }
        }
    }
    xs.sort_by(f64::total_cmp);
    // Pieces (x0, x1, h(x0), slope, ∫ over the piece).
    let mut pieces = Vec::with_capacity(xs.len());
    let mut total = 0.0;
    for w in xs.windows(2) {
        let (x0, x1) = (w[0], w[1]);
        if x1 <= x0 {
            continue;
        }
        let (h0, h1) = (h(x0), h(x1));
        let slope = (h1 - h0) / (x1 - x0);
        let integral = if slope.abs() <= 1e-12 * h0.max(h1) / (x1 - x0) {
            (x1 - x0) / h0
        } else {
            (h1 / h0).ln() / slope
        };
        total += integral;
        pieces.push((x0, x1, h0, slope, integral));
    }
    if !total.is_finite() || total > MAX_SMOOTH_LINES as f64 {
        return Err(Error::Invalid(format!(
            "smooth_lines: gap [{a}, {b}] would need more than {MAX_SMOOTH_LINES} lines"
        )));
    }
    let cells = (total.ceil() as usize).max(1);
    let mut pts = Vec::with_capacity(cells - 1);
    let mut piece = 0;
    let mut before = 0.0;
    for n in 1..cells {
        let target = total * n as f64 / cells as f64;
        while piece + 1 < pieces.len() && before + pieces[piece].4 < target {
            before += pieces[piece].4;
            piece += 1;
        }
        let (x0, x1, h0, slope, _) = pieces[piece];
        let f = target - before;
        let x = if slope.abs() <= 1e-12 * h0 / (x1 - x0) {
            x0 + f * h0
        } else {
            x0 + h0 * ((slope * f).exp() - 1.0) / slope
        };
        pts.push(x.clamp(x0, x1));
    }
    Ok(pts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lines() {
        let mut g = RectilinearGrid::new();
        g.set_lines(0, [3.0, 1.0, 2.0, 1.0]).unwrap();
        assert_eq!(g.lines(0), &[1.0, 2.0, 3.0]);
        g.add_line(0, 2.5).unwrap();
        assert_eq!(g.nearest_line(0, 2.4), Some(2));
        assert_eq!(g.snap(0, 2.5), Some(2));
        assert_eq!(g.snap(0, 2.4), None);
        assert_eq!(g.cell_index(0, 1.5), Some(0));
        assert_eq!(g.cell_index(0, 3.0), Some(2));
        assert_eq!(g.cell_index(0, 0.5), None);
        assert!(g.set_lines(1, [f64::NAN]).is_err());
    }

    #[test]
    fn smoothing() {
        let fixed = [0.0, 0.1, 10.0];
        let s = smooth_lines(&fixed, 1.0, 1.5).unwrap();
        for f in fixed {
            assert!(s.contains(&f));
        }
        let d: Vec<f64> = s.windows(2).map(|w| w[1] - w[0]).collect();
        assert!(d.iter().all(|&x| x <= 1.0 + 1e-9), "{d:?}");
        for w in d.windows(2) {
            let q = w[1].max(w[0]) / w[1].min(w[0]);
            assert!(q <= 1.5 * 1.1, "ratio {q} in {d:?}");
        }
        assert!(smooth_lines(&fixed, 0.0, 1.5).is_err());
    }

    fn check_graded(s: &[f64], max_res: f64, ratio: f64) {
        let d: Vec<f64> = s.windows(2).map(|w| w[1] - w[0]).collect();
        assert!(
            d.iter().all(|&x| x > 0.0 && x <= max_res * (1.0 + 1e-9)),
            "{d:?}"
        );
        for w in d.windows(2) {
            let q = w[1].max(w[0]) / w[1].min(w[0]);
            assert!(q <= ratio * 1.05, "ratio {q} in {d:?}");
        }
    }

    #[test]
    fn small_gap_next_to_large_one() {
        // Regression: the numeric integration was too coarse when the size
        // function varied by orders of magnitude, giving ratios above 3.
        let s = smooth_lines(&[0.0, 1e-3, 10.0], 1.0, 1.5).unwrap();
        // the tiny fixed cell is kept; grading starts after it
        assert!((s[1] - 1e-3).abs() < 1e-15);
        check_graded(&s[1..], 1.0, 1.5);
        let s = smooth_lines(&[-50.0, 0.0, 3.5e-5, 50.0], 2.0, 1.3).unwrap();
        let i = s.iter().position(|&x| x == 3.5e-5).unwrap();
        check_graded(&s[i..], 2.0, 1.3);
        let z = s.iter().position(|&x| x == 0.0).unwrap();
        check_graded(&s[..=z], 2.0, 1.3);
    }

    #[test]
    fn smoothing_limits() {
        // ratio 1 is uniform at the smallest neighbour size
        let s = smooth_lines(&[0.0, 0.5, 3.0], 1.0, 1.0).unwrap();
        check_graded(&s, 1.0, 1.0);
        // absurd counts are rejected instead of exhausting memory
        assert!(smooth_lines(&[0.0, 1e12], 1e-3, 1.5).is_err());
        assert!(smooth_lines(&[0.0, 1e-9, 1e9], 1.0, 1.0).is_err());
        assert!(smooth_lines(&[0.0, 1.0], 0.1, f64::INFINITY).is_err());
        assert!(smooth_lines(&[0.0, 1.0], f64::NAN, 1.5).is_err());
        // infinite max_res: only grading constraints apply
        let s = smooth_lines(&[0.0, 1.0, 100.0], f64::INFINITY, 2.0).unwrap();
        assert!(s.len() > 3);
    }
}
