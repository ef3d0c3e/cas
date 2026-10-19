//! Tests of the handle, the interner and node release.
//!
//! The interner is global and tests run in parallel, so tests that look at
//! interner entries use symbols of their own and count entries by hash.

#![allow(
    clippy::many_single_char_names,
    reason = "f, g, x, y mirror the notation"
)]

use std::sync::{Arc, Barrier, Weak};

use proptest::prelude::*;

use crate::expr::Repr;
use crate::intern::{self, entries_with_hash};
use crate::node::Key;
use crate::{Attributes, Complex, Expr, Integer, Rational, Symbol, View};

fn sym(name: &str) -> Expr {
    Symbol::new(name).into()
}

fn node(e: &Expr) -> &Arc<crate::node::Node> {
    match &e.0 {
        Repr::Node(n) => n,
        _ => panic!("{e:?} is not a node"),
    }
}

/// The hashes of every node in `e`, outermost first.
fn node_hashes(e: &Expr, out: &mut Vec<u64>) {
    if let Repr::Node(n) = &e.0 {
        out.push(n.hash);
    }
    if let Some(h) = e.head() {
        node_hashes(h, out);
    }
    for a in e.args() {
        node_hashes(a, out);
    }
}

#[test]
fn handle_is_16_bytes() {
    fn send_sync<T: Send + Sync>() {}
    assert_eq!(size_of::<Expr>(), 16);
    assert_eq!(size_of::<Option<Expr>>(), 16);
    send_sync::<Expr>();
    send_sync::<Symbol>();
}

#[test]
fn symbols_intern_by_name() {
    let a = Symbol::new("symbols_intern_by_name_a");
    assert_eq!(a, Symbol::new("symbols_intern_by_name_a"));
    assert_ne!(a, Symbol::new("symbols_intern_by_name_b"));
    assert_eq!(a.name(), "symbols_intern_by_name_a");
    assert_eq!(format!("{a:?}"), "symbols_intern_by_name_a");
}

#[test]
fn concurrent_symbol_creation_agrees() {
    const THREADS: usize = 8;
    let barrier = Barrier::new(THREADS);
    let results: Vec<Vec<Symbol>> = std::thread::scope(|s| {
        let handles: Vec<_> = (0..THREADS)
            .map(|_| {
                s.spawn(|| {
                    barrier.wait();
                    (0..2000)
                        .map(|i| Symbol::new(&format!("concurrent_symbol_{i}")))
                        .collect()
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for r in &results[1..] {
        assert_eq!(&results[0], r);
    }
    for (i, s) in results[0].iter().enumerate() {
        assert_eq!(s.name(), format!("concurrent_symbol_{i}"));
    }
}

#[test]
fn same_contents_share_a_node() {
    let (f, g, x) = (sym("f"), sym("g"), sym("x"));
    let build = || {
        let inner = Expr::apply(&g, &[x.clone(), Expr::from(-7)]);
        Expr::apply(&f, &[x.clone(), Expr::from(1), Expr::string("a"), inner])
    };
    let (a, b) = (build(), build());
    assert_eq!(a, b);
    assert!(Arc::ptr_eq(node(&a), node(&b)));
    assert!(Arc::ptr_eq(node(&a.args()[3]), node(&b.args()[3])));
    assert_eq!(a.structural_hash(), b.structural_hash());
    assert_eq!(format!("{a:?}"), r#"f[x, 1, "a", g[x, -7]]"#);
}

#[test]
fn different_contents_differ() {
    let (f, g, x, y) = (sym("f"), sym("g"), sym("x"), sym("y"));
    let fx = Expr::apply(&f, std::slice::from_ref(&x));
    let distinct = [
        Expr::apply(&f, &[x.clone(), y.clone()]),
        Expr::apply(&f, &[y.clone(), x.clone()]),
        fx.clone(),
        Expr::apply(&g, std::slice::from_ref(&x)),
        Expr::apply(&f, &[]),
        f.clone(),
        x.clone(),
        Expr::string("x"),
        Expr::string(""),
        Expr::from(0),
        Expr::from(1),
        Expr::apply(&fx, std::slice::from_ref(&y)),
        Expr::apply(&f, &[Expr::from(0)]),
        Expr::apply(&f, &[Expr::string("0")]),
        Expr::apply(&f, std::slice::from_ref(&fx)),
    ];
    for (i, a) in distinct.iter().enumerate() {
        for (j, b) in distinct.iter().enumerate() {
            assert_eq!(i == j, a == b, "{a:?} vs {b:?}");
        }
    }
}

#[test]
fn views_and_accessors() {
    let (f, x) = (sym("f"), sym("x"));
    let e = Expr::apply(&f, &[x.clone(), Expr::from(i64::MIN), Expr::string("s")]);
    let View::Apply { head, args } = e.view() else {
        panic!("not an application");
    };
    assert_eq!(head, &f);
    assert_eq!(args.len(), 3);
    assert_eq!(e.head(), Some(&f));
    assert!(!e.is_atom());
    assert_eq!(args[0].as_symbol(), Some(Symbol::new("x")));
    assert_eq!(args[1].as_i64(), Some(i64::MIN));
    assert_eq!(args[2].as_str(), Some("s"));
    for a in args {
        assert!(a.is_atom());
        assert!(a.head().is_none());
        assert!(a.args().is_empty());
    }
    assert_eq!(x.as_i64(), None);
    assert_eq!(x.as_str(), None);
    assert_eq!(Expr::from(3).as_symbol(), None);
    assert_eq!(Expr::apply(&f, &[]).args(), &[] as &[Expr]);
}

#[test]
fn metadata() {
    let (f, g, x, y) = (sym("meta_f"), sym("meta_g"), sym("meta_x"), sym("meta_y"));
    let gx = Expr::apply(&g, &[x.clone(), Expr::from(2)]);
    let e = Expr::apply(&f, &[gx.clone(), gx.clone(), Expr::string("s")]);
    assert_eq!((x.depth(), x.size()), (1, 1));
    assert_eq!((gx.depth(), gx.size()), (2, 4));
    // Shared subterms count once per occurrence.
    assert_eq!((e.depth(), e.size()), (3, 1 + 4 + 4 + 1 + 1));
    for s in [&f, &g, &x] {
        assert!(e.may_contain(s.as_symbol().unwrap()));
    }
    // A symbol whose Bloom bit gx does not have is reported absent.
    let absent = (0..10_000)
        .map(|i| Symbol::new(&format!("meta_absent_{i}")))
        .find(|s| Expr::from(*s).meta().symbols & gx.meta().symbols == 0)
        .unwrap();
    assert!(!gx.may_contain(absent));
    assert!(!Expr::from(1).may_contain(y.as_symbol().unwrap()));
    assert!(Expr::from(1).is_numeric());
    assert!(!x.is_numeric() && !gx.is_numeric() && !Expr::string("1").is_numeric());
}

#[test]
fn builtins_have_fixed_ids() {
    let builtins: Vec<Symbol> = Symbol::builtins().collect();
    assert!(builtins.len() > 20);
    for (i, &s) in builtins.iter().enumerate() {
        assert_eq!(s.id() as usize, i);
        assert!(s.is_builtin());
        assert_eq!(Symbol::new(s.name()), s, "{s:?} interns to its constant");
    }
    let names: std::collections::HashSet<&str> = builtins.iter().map(|s| s.name()).collect();
    assert_eq!(names.len(), builtins.len(), "names are unique");
    // The language escapes user names that start with an ASCII capital
    // (CAS-DOC-7 §5), which is safe only while every builtin's name does.
    for name in names {
        assert!(
            name.starts_with(|c: char| c.is_ascii_uppercase())
                && name.chars().all(|c| c.is_ascii_alphanumeric()),
            "{name}: a builtin's name is an ASCII capital, then ASCII letters and digits"
        );
    }
    assert_eq!(Symbol::PLUS.name(), "Plus");
    assert_eq!(Symbol::BLANK_NULL_SEQUENCE.name(), "BlankNullSequence");
    assert!(!Symbol::new("builtins_user").is_builtin());
}

#[test]
fn builtin_attributes() {
    use Attributes as A;
    let attrs = |s: Symbol| s.builtin_attributes().unwrap();
    for s in Symbol::builtins() {
        assert!(attrs(s).contains(A::PROTECTED), "{s:?} is protected");
    }
    let arith =
        A::FLAT | A::ORDERLESS | A::ONE_IDENTITY | A::LISTABLE | A::NUMERIC_FUNCTION | A::PROTECTED;
    assert_eq!(attrs(Symbol::PLUS), arith);
    assert_eq!(attrs(Symbol::TIMES), arith);
    assert_eq!(
        attrs(Symbol::POWER),
        A::ONE_IDENTITY | A::LISTABLE | A::NUMERIC_FUNCTION | A::PROTECTED
    );
    assert_eq!(attrs(Symbol::PI), A::CONSTANT | A::PROTECTED);
    assert_eq!(attrs(Symbol::LIST), A::PROTECTED);
    assert_eq!(attrs(Symbol::HOLD), A::HOLD_ALL | A::PROTECTED);
    assert_eq!(attrs(Symbol::PATTERN), A::HOLD_FIRST | A::PROTECTED);
    assert_eq!(attrs(Symbol::PATTERN_TEST), A::HOLD_REST | A::PROTECTED);
    assert_eq!(A::HOLD_ALL, A::HOLD_FIRST | A::HOLD_REST);
    assert_eq!(
        Symbol::new("builtin_attributes_f").builtin_attributes(),
        None
    );
}

#[test]
fn attributes_set_operations_and_debug() {
    use Attributes as A;
    let a = A::FLAT | A::HOLD_FIRST;
    assert!(a.contains(A::FLAT) && !a.contains(A::HOLD_ALL));
    assert!(a.contains(A::NONE) && A::NONE.is_empty() && A::default().is_empty());
    assert_eq!(a & A::HOLD_ALL, A::HOLD_FIRST);
    let mut b = a;
    b |= A::HOLD_REST;
    assert!(b.contains(A::HOLD_ALL));
    assert_eq!(format!("{:?}", A::NONE), "None");
    assert_eq!(format!("{a:?}"), "Flat | HoldFirst");
    assert_eq!(format!("{b:?}"), "Flat | HoldAll");
    assert_eq!(
        format!("{:?}", Symbol::PLUS.builtin_attributes().unwrap()),
        "Flat | Orderless | OneIdentity | Listable | NumericFunction | Protected"
    );
}

#[test]
fn numeric_propagates_through_builtin_functions() {
    let (plus, sin, pi) = (
        Symbol::PLUS.into(),
        Symbol::SIN.into(),
        Expr::from(Symbol::PI),
    );
    let (f, x) = (sym("numeric_f"), sym("numeric_x"));
    let one = Expr::from(1);
    assert!(pi.is_numeric() && Expr::from(Symbol::E).is_numeric());
    assert!(!Expr::from(Symbol::PLUS).is_numeric());
    let sum = Expr::apply(&plus, &[one.clone(), pi.clone()]);
    assert!(sum.is_numeric());
    assert!(Expr::apply(&sin, std::slice::from_ref(&sum)).is_numeric());
    assert!(!Expr::apply(&plus, &[one.clone(), x.clone()]).is_numeric());
    assert!(!Expr::apply(&plus, &[one.clone(), Expr::string("1")]).is_numeric());
    assert!(
        !Expr::apply(&f, std::slice::from_ref(&one)).is_numeric(),
        "user heads are unknown"
    );
    assert!(!Expr::apply(&plus, &[]).is_numeric(), "no arguments");
    let list = Expr::apply(&Symbol::LIST.into(), std::slice::from_ref(&one));
    assert!(!list.is_numeric());
    // A compound head is not a numeric function, even if numeric itself.
    assert!(!Expr::apply(&sum, &[one]).is_numeric());
}

#[test]
fn has_pattern_propagates() {
    let blank = Expr::apply(&Symbol::BLANK.into(), &[]);
    let (f, x) = (sym("pattern_f"), sym("pattern_x"));
    assert!(blank.has_pattern());
    assert!(
        !Expr::from(Symbol::BLANK).has_pattern(),
        "the bare symbol is literal"
    );
    assert!(!Expr::apply(&f, &[Expr::from(Symbol::BLANK)]).has_pattern());
    let px = Expr::apply(&Symbol::PATTERN.into(), &[x.clone(), blank.clone()]);
    let deep = Expr::apply(&f, &[x.clone(), Expr::apply(&f, std::slice::from_ref(&px))]);
    assert!(px.has_pattern() && deep.has_pattern());
    assert!(
        Expr::apply(&blank, std::slice::from_ref(&x)).has_pattern(),
        "through the head"
    );
    let alt = Expr::apply(&Symbol::ALTERNATIVES.into(), &[x.clone(), f.clone()]);
    assert!(alt.has_pattern());
    assert!(!Expr::apply(&f, &[x, Expr::from(1), Expr::string("_")]).has_pattern());
}

#[test]
fn bloom_bits_spread() {
    let bits: std::collections::HashSet<u64> = (0..512)
        .map(|i| sym(&format!("bloom_{i}")).meta().symbols)
        .collect();
    assert_eq!(bits.len(), 64, "every bit is used");
}

#[test]
fn dropping_removes_entries() {
    let (f, x) = (sym("dropping_f"), sym("dropping_x"));
    let inner = Expr::apply(&f, &[x.clone(), Expr::string("dropping")]);
    let outer = Expr::apply(&f, &[inner.clone(), inner]);
    let mut hashes = Vec::new();
    node_hashes(&outer, &mut hashes);
    assert_eq!(hashes.len(), 5); // outer, inner twice, the string twice
    for &h in &hashes {
        assert_eq!(entries_with_hash(h), 1);
    }
    drop(outer);
    for &h in &hashes {
        assert_eq!(entries_with_hash(h), 0);
    }
}

#[test]
fn shared_node_outlives_its_parent() {
    let (f, x) = (sym("outlives_f"), sym("outlives_x"));
    let inner = Expr::apply(&f, std::slice::from_ref(&x));
    let outer = Expr::apply(&f, std::slice::from_ref(&inner));
    let (hi, ho) = (inner.structural_hash(), outer.structural_hash());
    drop(outer);
    assert_eq!((entries_with_hash(hi), entries_with_hash(ho)), (1, 0));
    assert_eq!(Expr::apply(&f, std::slice::from_ref(&x)), inner);
    drop(inner);
    assert_eq!(entries_with_hash(hi), 0);
}

#[test]
fn deep_drop_uses_constant_stack() {
    const DEPTH: u32 = 1_000_000;
    // Far too small for a recursive drop of DEPTH levels.
    std::thread::Builder::new()
        .stack_size(64 * 1024)
        .spawn(|| {
            let (f, x) = (sym("deep_f"), sym("deep_x"));
            let mut e = x;
            for _ in 0..DEPTH {
                e = Expr::apply(&f, &[e]);
            }
            assert_eq!(e.depth(), u16::MAX);
            assert_eq!(e.size(), 2 * DEPTH + 1);
            let h = e.structural_hash();
            drop(e);
            assert_eq!(entries_with_hash(h), 0);
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn wide_drop() {
    let f = sym("wide_f");
    let args: Vec<Expr> = (0..100_000)
        .map(|i| Expr::apply(&f, &[Expr::from(i)]))
        .collect();
    let probe = args[777].structural_hash();
    let e = Expr::apply(&f, &args);
    drop(args);
    assert_eq!(entries_with_hash(probe), 1);
    drop(e);
    assert_eq!(entries_with_hash(probe), 0);
}

#[test]
fn concurrent_construction_agrees() {
    const THREADS: usize = 8;
    let (f, g) = (sym("concurrent_f"), sym("concurrent_g"));
    let build = |f: &Expr, g: &Expr| -> Vec<Expr> {
        (0..2000)
            .map(|i| {
                let inner = Expr::apply(
                    g,
                    &[Expr::from(i % 37), Expr::string(&(i % 11).to_string())],
                );
                Expr::apply(f, &[inner, Expr::from(i)])
            })
            .collect()
    };
    let barrier = Barrier::new(THREADS);
    let results: Vec<Vec<Expr>> = std::thread::scope(|s| {
        let handles: Vec<_> = (0..THREADS)
            .map(|_| {
                s.spawn(|| {
                    barrier.wait();
                    build(&f, &g)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for r in &results[1..] {
        for (a, b) in results[0].iter().zip(r) {
            assert!(Arc::ptr_eq(node(a), node(b)));
        }
    }
}

/// Threads build and drop the same few expressions as fast as they can, so
/// lookups keep meeting nodes whose last reference is being dropped.
#[test]
fn concurrent_churn() {
    const THREADS: usize = 8;
    const KEYS: usize = 16;
    let (f, x) = (sym("churn_f"), sym("churn_x"));
    let key = |i: usize| Expr::apply(&f, &[x.clone(), Expr::from(i64::try_from(i).unwrap())]);
    // Half of the keys stay alive throughout; the others come and go.
    let held: Vec<Expr> = (0..KEYS / 2).map(key).collect();
    let hashes: Vec<u64> = (0..KEYS).map(|i| key(i).structural_hash()).collect();
    std::thread::scope(|s| {
        for t in 0..THREADS {
            let (key, held) = (&key, &held);
            s.spawn(move || {
                let mut alive = Vec::new();
                for n in 0..20_000 {
                    let i = (n * 7 + t) % KEYS;
                    let e = key(i);
                    if i < KEYS / 2 {
                        assert_eq!(e, held[i]);
                    }
                    // Two live handles to the same key are always one node.
                    if let Some(other) = alive.iter().find(|a: &&Expr| a.args()[1] == e.args()[1]) {
                        assert_eq!(*other, e);
                    }
                    alive.push(e);
                    if alive.len() > 3 {
                        alive.remove(0);
                    }
                }
            });
        }
    });
    for (i, &h) in hashes.iter().enumerate() {
        let expected = usize::from(i < held.len());
        assert_eq!(entries_with_hash(h), expected, "key {i}");
    }
    drop(held);
    for &h in &hashes {
        assert_eq!(entries_with_hash(h), 0);
    }
}

/// A full 64-bit hash collision: the lookup must not return the other node,
/// and must keep the upgraded candidate alive until after the lock.
#[test]
fn colliding_hashes_stay_apart() {
    let (f, x) = (sym("collide_f"), sym("collide_x"));
    let a = Expr::apply(&f, std::slice::from_ref(&x));
    let b = Expr::apply(&f, &[x.clone(), x.clone()]);
    let (ka, kb) = (
        Key::Apply {
            head: &f,
            args: a.args(),
        },
        Key::Apply {
            head: &f,
            args: b.args(),
        },
    );
    let mut shard = intern::Shard::new_for_test();
    let mut mismatched = Vec::new();
    let got = shard.get_or_insert(0, ka, &mut mismatched, || node(&a).clone());
    assert!(Arc::ptr_eq(&got, node(&a)) && mismatched.is_empty());
    let got = shard.get_or_insert(0, kb, &mut mismatched, || node(&b).clone());
    assert!(Arc::ptr_eq(&got, node(&b)));
    assert_eq!(mismatched.len(), 1);
    assert!(Arc::ptr_eq(&mismatched[0], node(&a)));
    let got = shard.get_or_insert(0, ka, &mut Vec::new(), || panic!("must hit"));
    assert!(Arc::ptr_eq(&got, node(&a)));
}

/// A dead entry and a live one for the same contents, as while a node's last
/// reference is being dropped and the contents are interned again: lookups
/// skip the dead entry, and removal takes the entry of the node removed.
#[test]
fn dead_and_live_entries_coexist() {
    let (f, x) = (sym("coexist_f"), sym("coexist_x"));
    let a = Expr::apply(&f, std::slice::from_ref(&x));
    let (n, h) = (node(&a), a.structural_hash());
    let key = Key::Apply {
        head: &f,
        args: a.args(),
    };
    for dead_first in [true, false] {
        let mut shard = intern::Shard::new_for_test();
        for dead in [dead_first, !dead_first] {
            shard.insert_raw(h, if dead { Weak::new() } else { Arc::downgrade(n) });
        }
        let got = shard.get_or_insert(h, key, &mut Vec::new(), || panic!("must hit"));
        assert!(Arc::ptr_eq(&got, n));
        assert!(shard.remove(n));
        let left = shard.targets(h);
        assert_eq!(left.len(), 1);
        assert_ne!(left[0], Arc::as_ptr(n), "removed the dead entry instead");
        assert!(!shard.remove(n));
    }
}

/// The model's symbols: user symbols and builtins with each kind of flag.
const SYMS: [&str; 8] = [
    "prop_0", "prop_1", "prop_2", "Plus", "Sin", "Pi", "Blank", "Pattern",
];

/// A plain model of an expression.
#[derive(Clone, PartialEq, Debug)]
enum Tree {
    /// An index into [`SYMS`].
    Sym(usize),
    Int(i64),
    /// A number atom that is not an `i64`: a larger integer or a rational.
    Num(Rational),
    /// A complex number with a nonzero imaginary part.
    Cplx(Rational, Rational),
    Str(String),
    Apply(Box<Tree>, Vec<Tree>),
}

/// Integers just outside `i64` and non-integer rationals, small and large.
fn big_number() -> impl Strategy<Value = Rational> {
    let edge = (any::<bool>(), 0i128..3).prop_map(|(negative, k)| {
        let x = if negative {
            i128::from(i64::MIN) - 1 - k
        } else {
            i128::from(i64::MAX) + 1 + k
        };
        Rational::from(x)
    });
    let fraction = (-50i64..50, 2i64..9, 0u32..3).prop_filter_map("an integer", |(n, d, scale)| {
        let n = Integer::from(n) * Integer::from(i128::from(i64::MAX).pow(scale));
        let q = Rational::new(n, Integer::from(d));
        (!q.is_integer()).then_some(q)
    });
    prop_oneof![edge, fraction]
}

fn tree() -> impl Strategy<Value = Tree> {
    let leaf = prop_oneof![
        (0..SYMS.len()).prop_map(Tree::Sym),
        prop_oneof![
            Just(0i64),
            Just(1),
            Just(-1),
            Just(i64::MIN),
            Just(i64::MAX),
            -3i64..3
        ]
        .prop_map(Tree::Int),
        big_number().prop_map(Tree::Num),
        (
            prop_oneof![Just(Rational::ZERO), big_number()],
            big_number()
        )
            .prop_map(|(re, im)| Tree::Cplx(re, im)),
        "[ab]{0,2}".prop_map(Tree::Str),
    ];
    leaf.prop_recursive(4, 32, 4, |inner| {
        (inner.clone(), prop::collection::vec(inner, 0..4))
            .prop_map(|(h, a)| Tree::Apply(Box::new(h), a))
    })
}

fn build(t: &Tree) -> Expr {
    match t {
        Tree::Sym(i) => sym(SYMS[*i]),
        Tree::Int(x) => Expr::from(*x),
        Tree::Num(q) => Expr::from(q),
        Tree::Cplx(re, im) => Expr::from(Complex::new(re.clone(), im.clone())),
        Tree::Str(s) => Expr::string(s),
        Tree::Apply(h, a) => Expr::apply(&build(h), &a.iter().map(build).collect::<Vec<_>>()),
    }
}

fn unbuild(e: &Expr) -> Tree {
    match e.view() {
        View::Symbol(s) => Tree::Sym(SYMS.iter().position(|&n| n == s.name()).unwrap()),
        View::Int(x) => Tree::Int(x),
        View::Integer(x) => Tree::Num(Rational::from(x)),
        View::Rational(x) => Tree::Num(x.clone()),
        View::Complex(z) => Tree::Cplx(z.re().clone(), z.im().clone()),
        View::String(s) => Tree::Str(s.to_owned()),
        View::Apply { head, args } => {
            Tree::Apply(Box::new(unbuild(head)), args.iter().map(unbuild).collect())
        }
    }
}

fn model_depth_size(t: &Tree) -> (u16, u32) {
    match t {
        Tree::Apply(h, a) => {
            let (mut d, mut s) = model_depth_size(h);
            for x in a {
                let (dx, sx) = model_depth_size(x);
                d = d.max(dx);
                s += sx;
            }
            (d + 1, s + 1)
        }
        _ => (1, 1),
    }
}

fn model_numeric(t: &Tree) -> bool {
    match t {
        Tree::Int(_) | Tree::Num(_) | Tree::Cplx(..) => true,
        Tree::Sym(i) => SYMS[*i] == "Pi",
        Tree::Apply(h, a) => {
            matches!(**h, Tree::Sym(i) if ["Plus", "Sin"].contains(&SYMS[i]))
                && !a.is_empty()
                && a.iter().all(model_numeric)
        }
        Tree::Str(_) => false,
    }
}

fn model_has_pattern(t: &Tree) -> bool {
    match t {
        Tree::Apply(h, a) => {
            matches!(**h, Tree::Sym(i) if ["Blank", "Pattern"].contains(&SYMS[i]))
                || model_has_pattern(h)
                || a.iter().any(model_has_pattern)
        }
        _ => false,
    }
}

fn symbols(t: &Tree, out: &mut Vec<usize>) {
    match t {
        Tree::Sym(i) => out.push(*i),
        Tree::Apply(h, a) => {
            symbols(h, out);
            for x in a {
                symbols(x, out);
            }
        }
        _ => {}
    }
}

proptest! {
    #[test]
    fn equality_is_structural(a in tree(), b in tree()) {
        let (ea, eb) = (build(&a), build(&b));
        prop_assert_eq!(a == b, ea == eb);
        prop_assert_eq!(ea.clone(), build(&a));
        if a == b {
            prop_assert_eq!(ea.structural_hash(), eb.structural_hash());
        }
    }

    #[test]
    fn view_round_trips(t in tree()) {
        prop_assert_eq!(unbuild(&build(&t)), t);
    }

    #[test]
    fn metadata_matches_model(t in tree()) {
        let e = build(&t);
        prop_assert_eq!((e.depth(), e.size()), model_depth_size(&t));
        let mut syms = Vec::new();
        symbols(&t, &mut syms);
        for i in syms {
            let s = Symbol::new(SYMS[i]);
            prop_assert!(e.may_contain(s), "{:?} lacks {:?}", e, s);
        }
        prop_assert_eq!(e.is_numeric(), model_numeric(&t));
        prop_assert_eq!(e.has_pattern(), model_has_pattern(&t));
    }
}

/// Each number has one form: inline exactly when it is an integer that fits
/// an `i64`, whichever constructor builds it.
#[test]
fn numbers_have_one_form() {
    let inline = |e: &Expr| matches!(e.0, Repr::Int(_));
    for x in [i64::MIN, i64::MIN + 1, -1, 0, 1, i64::MAX - 1, i64::MAX] {
        let from_i64 = Expr::from(x);
        let from_integer = Expr::from(Integer::from(x));
        let from_rational = Expr::from(Rational::from(x));
        assert!(inline(&from_i64) && inline(&from_integer) && inline(&from_rational));
        assert!(from_i64 == from_integer && from_i64 == from_rational);
        assert_eq!(from_i64.to_integer(), Some(Integer::from(x)));
        assert_eq!(from_i64.view(), View::Int(x));
    }
    let edges = [
        i128::from(i64::MIN) - 1,
        i128::from(i64::MAX) + 1,
        i128::from(u64::MAX),
        i128::MIN,
    ];
    for x in edges {
        let n = Integer::from(x);
        let a = Expr::from(&n);
        let b = Expr::from(Rational::new(
            Integer::from(x) * Integer::from(3),
            Integer::from(3),
        ));
        assert!(!inline(&a) && a == b);
        assert!(Arc::ptr_eq(node(&a), node(&b)));
        assert_eq!(a.view(), View::Integer(&n));
        assert_eq!(a.to_integer(), Some(n.clone()));
        assert_eq!(a.to_rational(), Some(Rational::from(n)));
        assert!(a.is_number() && a.is_numeric() && a.is_atom());
    }
    let half = Rational::new(Integer::from(1), Integer::from(2));
    let a = Expr::from(&half);
    assert_eq!(
        a,
        Expr::from(Rational::new(Integer::from(-7), Integer::from(-14)))
    );
    assert_eq!(a.view(), View::Rational(&half));
    assert_eq!(a.to_integer(), None);
    assert_eq!(a.to_rational(), Some(half));
    assert!(a.is_number() && a.is_numeric() && a.as_i64().is_none());
    // Different kinds with related payloads stay apart.
    let big = Expr::from(Integer::from(1i128 << 64));
    assert_ne!(
        big,
        Expr::from(Rational::new(Integer::from(1i128 << 64), Integer::from(3)))
    );
    assert!(!sym("Pi").is_number() && sym("Pi").is_numeric());
    assert!(!Expr::string("1").is_number());
}

/// A complex number is a node only when it is not real; a real value is the
/// rational or integer atom, so each number keeps one form (CAS-51).
#[test]
fn complex_numbers_have_one_form() {
    let q = |n: i64, d: i64| Rational::new(n.into(), d.into());
    let z = |re: Rational, im: Rational| Complex::new(re, im);
    // Real values, whatever builds them.
    for re in [
        q(0, 1),
        q(5, 1),
        q(-7, 2),
        Rational::from(Integer::from(1i128 << 64)),
    ] {
        let e = Expr::from(z(re.clone(), Rational::ZERO));
        assert_eq!(e, Expr::from(&re));
        assert_eq!(e.to_complex(), Some(Complex::from(re)));
    }
    assert_eq!(Expr::from(Complex::ZERO), Expr::from(0));
    // Non-real values: one shared node each, apart from every real.
    let a = Expr::from(z(q(2, 1), q(-3, 4)));
    let b = Expr::from(z(q(4, 2), q(6, -8)));
    assert!(Arc::ptr_eq(node(&a), node(&b)));
    assert_eq!(a.view(), View::Complex(&z(q(2, 1), q(-3, 4))));
    assert_eq!(a.to_complex(), Some(z(q(2, 1), q(-3, 4))));
    assert_eq!((a.to_rational(), a.to_integer()), (None, None));
    assert!(a.is_number() && a.is_numeric() && a.is_atom());
    assert_ne!(a, Expr::from(z(q(-3, 4), q(2, 1))));
    assert_ne!(Expr::from(Complex::I), Expr::from(1));
    assert_eq!(crate::build::i(), Expr::from(Complex::I));
    // The symbol `I` evaluates to the number; it is not the number.
    let unit = Expr::from(Symbol::I);
    assert_ne!(unit, crate::build::i());
    assert!(!unit.is_number() && unit.is_numeric());
    assert_eq!(-crate::build::i(), Expr::from(-Complex::I));
}
