//! A small bounding volume hierarchy over axis-aligned boxes.

use crate::math::BoundingBox;

/// Maximum number of items stored in a leaf.
const LEAF_SIZE: usize = 4;

#[derive(Debug, Clone, PartialEq)]
struct Node {
    bbox: BoundingBox,
    /// Leaf: first index into `items`. Internal: index of the right child
    /// (the left child is always the next node).
    start: u32,
    /// Number of items for leaves, 0 for internal nodes.
    count: u32,
}

/// Bounding volume hierarchy; items are referred to by their index in the
/// slice of boxes passed to [`Bvh::build`].
#[derive(Debug, Clone, PartialEq, Default)]
pub(crate) struct Bvh {
    nodes: Vec<Node>,
    items: Vec<u32>,
}

impl Bvh {
    /// Builds the hierarchy (median split along the longest centroid axis).
    /// Empty boxes are left out.
    pub(crate) fn build(boxes: &[BoundingBox]) -> Bvh {
        let mut items: Vec<u32> = (0..boxes.len() as u32)
            .filter(|&i| !boxes[i as usize].is_empty())
            .collect();
        let mut nodes = Vec::with_capacity(2 * items.len() / LEAF_SIZE + 1);
        if !items.is_empty() {
            let centers: Vec<[f64; 3]> = boxes.iter().map(|b| b.center().0).collect();
            build_node(boxes, &centers, &mut items, 0, &mut nodes);
        }
        Bvh { nodes, items }
    }

    /// Calls `f` for every item in a leaf whose node boxes all satisfy
    /// `overlaps`. Traversal stops early when `f` returns `false`.
    pub(crate) fn visit(
        &self,
        overlaps: impl Fn(&BoundingBox) -> bool,
        mut f: impl FnMut(u32) -> bool,
    ) {
        if self.nodes.is_empty() {
            return;
        }
        let mut stack: Vec<u32> = Vec::with_capacity(64);
        stack.push(0);
        while let Some(n) = stack.pop() {
            let node = &self.nodes[n as usize];
            if !overlaps(&node.bbox) {
                continue;
            }
            if node.count > 0 {
                let s = node.start as usize;
                for &it in &self.items[s..s + node.count as usize] {
                    if !f(it) {
                        return;
                    }
                }
            } else {
                stack.push(node.start);
                stack.push(n + 1);
            }
        }
    }
}

fn build_node(
    boxes: &[BoundingBox],
    centers: &[[f64; 3]],
    items: &mut [u32],
    offset: usize,
    nodes: &mut Vec<Node>,
) -> u32 {
    let bbox = items
        .iter()
        .fold(BoundingBox::EMPTY, |a, &i| a.union(&boxes[i as usize]));
    let me = nodes.len() as u32;
    nodes.push(Node {
        bbox,
        start: offset as u32,
        count: items.len() as u32,
    });
    if items.len() <= LEAF_SIZE {
        return me;
    }
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for &i in items.iter() {
        let c = centers[i as usize];
        for k in 0..3 {
            // Unbounded boxes have non-finite centers; keep them out of the
            // split heuristic.
            if c[k].is_finite() {
                lo[k] = lo[k].min(c[k]);
                hi[k] = hi[k].max(c[k]);
            }
        }
    }
    let axis = (0..3)
        .max_by(|&a, &b| (hi[a] - lo[a]).total_cmp(&(hi[b] - lo[b])))
        .unwrap_or(0);
    let mid = items.len() / 2;
    items.select_nth_unstable_by(mid, |&a, &b| {
        centers[a as usize][axis].total_cmp(&centers[b as usize][axis])
    });
    let (left, right) = items.split_at_mut(mid);
    build_node(boxes, centers, left, offset, nodes);
    let r = build_node(boxes, centers, right, offset + mid, nodes);
    let node = &mut nodes[me as usize];
    node.start = r;
    node.count = 0;
    me
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Vec3;

    #[test]
    fn point_queries_match_brute_force() {
        let mut seed = 7u64;
        let mut rnd = || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let boxes: Vec<BoundingBox> = (0..500)
            .map(|_| {
                let c = Vec3::new(rnd() * 10.0, rnd() * 10.0, rnd() * 10.0);
                let d = Vec3::new(rnd(), rnd(), rnd());
                BoundingBox::from_corners(c - d, c + d)
            })
            .collect();
        let bvh = Bvh::build(&boxes);
        for _ in 0..2000 {
            let p = Vec3::new(rnd() * 10.0, rnd() * 10.0, rnd() * 10.0);
            let mut got = Vec::new();
            bvh.visit(
                |b| b.contains(p),
                |i| {
                    if boxes[i as usize].contains(p) {
                        got.push(i);
                    }
                    true
                },
            );
            got.sort_unstable();
            let want: Vec<u32> = (0..boxes.len() as u32)
                .filter(|&i| boxes[i as usize].contains(p))
                .collect();
            assert_eq!(got, want);
        }
        assert!(Bvh::build(&[]).nodes.is_empty());
    }
}
