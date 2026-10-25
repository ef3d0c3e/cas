//! Multifactor Hensel lifting over ℤ (CAS-181) \[vzGG13 §15.4–15.5\]: from
//! `f ≡ lc(f) g₁ ⋯ gᵣ (mod p)` to the same modulo `pᵏ`, through a tree of
//! the factors whose inner nodes carry Bézout coefficients, by quadratic
//! Hensel steps.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of the Hensel step: f = g h, s g + t h = 1, e, q, r, c, d"
)]

use cas_num::Integer;
use cas_num::algebra::{EuclideanDomain, Integers, Ring, WordPrimeField};
use cas_num::modular::{ModRing, Residue};

use crate::cancel::{Cancel, Cancelled};
use crate::dense::{Dense, DenseRing};

/// A node of the factor tree: the product of the factors at its leaves,
/// monic, modulo the current modulus, its coefficients the representatives
/// in `[0, m)`.
struct Node {
    product: Dense<Integer>,
    split: Option<Box<Split>>,
}

/// An inner node's children, and `s`, `t` with `s · left + t · right ≡ 1`,
/// `deg s < deg right` and `deg t < deg left`.
struct Split {
    left: Node,
    right: Node,
    s: Dense<Integer>,
    t: Dense<Integer>,
}

/// The polynomial over ℤ of the representatives in `[0, p)`.
fn from_words(p: &Dense<u64>) -> Dense<Integer> {
    DenseRing::new(Integers).poly(p.coeffs().iter().map(|&c| Integer::from(c)).collect())
}

/// `p` modulo the modulus of `ring`.
fn residues(ring: &DenseRing<ModRing>, p: &Dense<Integer>) -> Dense<Residue> {
    let m = ring.base();
    ring.poly(p.coeffs().iter().map(|c| m.elem_integer(c)).collect())
}

/// The polynomial over ℤ of the representatives in `[0, m)`.
fn representatives(ring: &DenseRing<ModRing>, p: &Dense<Residue>) -> Dense<Integer> {
    let m = ring.base();
    DenseRing::new(Integers).poly(
        p.coeffs()
            .iter()
            .map(|c| Integer::from(m.to_natural(c)))
            .collect(),
    )
}

/// The tree over `factors`, split in halves, with the Bézout coefficients
/// modulo `p`, and the product of the factors.
fn build(fp: &DenseRing<WordPrimeField>, factors: &[Dense<u64>]) -> (Node, Dense<u64>) {
    if let [g] = factors {
        let leaf = Node {
            product: from_words(g),
            split: None,
        };
        return (leaf, g.clone());
    }
    let (left, a) = build(fp, &factors[..factors.len() / 2]);
    let (right, b) = build(fp, &factors[factors.len() / 2..]);
    let (g, s, t) = EuclideanDomain::gcd_ext(fp, &a, &b);
    debug_assert!(fp.is_one(&g), "coprime factors");
    let product = fp.mul(&a, &b);
    let node = Node {
        product: from_words(&product),
        split: Some(Box::new(Split {
            left,
            right,
            s: from_words(&s),
            t: from_words(&t),
        })),
    };
    (node, product)
}

/// One Hensel step [vzGG13 Algorithm 15.10] below `node`, whose product
/// becomes `v`, monic: from `v ≡ g h`, `s g + t h ≡ 1` modulo `m` to the same
/// modulo `next`, the modulus of `ring`, with `m | next | m²`. Then the
/// children lift to `g` and `h`.
fn lift_node(
    node: &mut Node,
    v: Dense<Integer>,
    ring: &DenseRing<ModRing>,
    cancel: &dyn Cancel,
) -> Result<(), Cancelled> {
    if let Some(split) = &mut node.split {
        let r = |p: &Dense<Integer>| residues(ring, p);
        let (f, g, h) = (r(&v), r(&split.left.product), r(&split.right.product));
        let (s, t) = (r(&split.s), r(&split.t));
        // e = f - g h; s e = q h + r; g* = g + t e + q g; h* = h + r.
        let e = ring.sub(&f, &ring.mul_with(&g, &h, cancel)?);
        let (q, rem) = ring.pseudo_div_rem(&ring.mul_with(&s, &e, cancel)?, &h);
        let te = ring.mul_with(&t, &e, cancel)?;
        let g1 = ring.add(&g, &ring.add(&te, &ring.mul_with(&q, &g, cancel)?));
        let h1 = ring.add(&h, &rem);
        // b = s g* + t h* - 1; s b = c h* + d; s* = s - d; t* = t - t b - c g*.
        let sg = ring.mul_with(&s, &g1, cancel)?;
        let b = ring.sub(
            &ring.add(&sg, &ring.mul_with(&t, &h1, cancel)?),
            &ring.one(),
        );
        let (c, d) = ring.pseudo_div_rem(&ring.mul_with(&s, &b, cancel)?, &h1);
        let s1 = ring.sub(&s, &d);
        let tb = ring.mul_with(&t, &b, cancel)?;
        let t1 = ring.sub(&t, &ring.add(&tb, &ring.mul_with(&c, &g1, cancel)?));
        split.s = representatives(ring, &s1);
        split.t = representatives(ring, &t1);
        lift_node(&mut split.left, representatives(ring, &g1), ring, cancel)?;
        lift_node(&mut split.right, representatives(ring, &h1), ring, cancel)?;
    }
    node.product = v;
    Ok(())
}

/// The leaves' products, from left to right.
fn leaves(node: Node, out: &mut Vec<Dense<Integer>>) {
    match node.split {
        None => out.push(node.product),
        Some(split) => {
            leaves(split.left, out);
            leaves(split.right, out);
        }
    }
}

/// The exponents of the quadratic steps to `k`: `⌈k/2ⁱ⌉` from the smallest
/// above 1 up to `k` itself, each at most twice the one before.
fn exponents(k: u64) -> Vec<u64> {
    let mut out = Vec::new();
    let mut e = k;
    while e > 1 {
        out.push(e);
        e = e.div_ceil(2);
    }
    out.reverse();
    out
}

/// The lifts of the local factors of `f` to `pᵏ`, `p` the characteristic of
/// `fp`. For `f ≡ lc(f) g₁ ⋯ gᵣ (mod p)`, `p` not dividing `lc(f)` and the
/// `gᵢ` monic and pairwise coprime modulo `p`, they are the monic
/// `Gᵢ ≡ gᵢ (mod p)` with `f ≡ lc(f) G₁ ⋯ Gᵣ (mod pᵏ)`, in the same order,
/// their coefficients in `[0, pᵏ)`. The lifts are unique, since the `gᵢ` are
/// coprime.
///
/// # Panics
///
/// If `factors` is empty, or `p` divides `lc(f)`.
pub(crate) fn lift(
    f: &Dense<Integer>,
    factors: &[Dense<u64>],
    fp: &DenseRing<WordPrimeField>,
    k: u64,
    cancel: &dyn Cancel,
) -> Result<Vec<Dense<Integer>>, Cancelled> {
    assert!(!factors.is_empty(), "a factor to lift");
    let p = fp.characteristic();
    let lead = f.lead().expect("a non-zero polynomial");
    let (mut tree, _) = build(fp, factors);
    for e in exponents(k) {
        let ring = DenseRing::new(ModRing::new(p.pow(e)));
        let m = ring.base();
        let inverse = m
            .inv(&m.elem_integer(lead))
            .expect("p does not divide the leading coefficient");
        let target = ring.scale(&residues(&ring, f), &inverse);
        lift_node(&mut tree, representatives(&ring, &target), &ring, cancel)?;
    }
    let mut out = Vec::with_capacity(factors.len());
    leaves(tree, &mut out);
    Ok(out)
}
