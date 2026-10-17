#![allow(clippy::many_single_char_names, reason = "x, y, e mirror the notation")]
#![allow(
    clippy::manual_assert_eq,
    reason = "`assert_eq!` would print a large or shared expression on failure"
)]

use std::collections::HashSet;

use proptest::prelude::*;

use crate::build::symbols;
use crate::{Expr, Symbol};

fn f(args: &[Expr]) -> Expr {
    Expr::apply(&Expr::symbol("f"), args)
}

/// `f[e, e]` nested `n` deep over `x`: a tree of `2^n` leaves with `n + 2`
/// distinct subexpressions (the levels, `x` and `f`).
fn doubling(n: usize) -> Expr {
    (0..n).fold(Expr::symbol("x"), |e, _| f(&[e.clone(), e]))
}

#[test]
fn parts() {
    let [x, y] = symbols(["x", "y"]);
    let e = f(&[x.clone(), Expr::apply(&f(&[]), std::slice::from_ref(&y))]);
    assert_eq!(e.part(&[]), Some(&e));
    assert_eq!(e.part(&[0]), Some(&Expr::symbol("f")));
    assert_eq!(e.part(&[1]), Some(&x));
    assert_eq!(e.part(&[2, 0]), Some(&f(&[])));
    assert_eq!(e.part(&[2, 1]), Some(&y));
    assert_eq!(e.part(&[3]), None);
    assert_eq!(e.part(&[1, 0]), None);
    assert_eq!(&e[2][1], &y);
    assert_eq!(&e[0], &Expr::symbol("f"));
}

#[test]
#[should_panic(expected = "no part 3")]
fn index_out_of_range() {
    let _ = &f(&[Expr::symbol("x")])[3];
}

#[test]
fn replace_part_shares_the_rest() {
    let [x, y, z] = symbols(["x", "y", "z"]);
    let left = f(&[x.clone(), y.clone()]);
    let e = f(&[left.clone(), f(std::slice::from_ref(&y))]);
    let r = e.replace_part(&[2, 1], z.clone()).unwrap();
    assert_eq!(format!("{r:?}"), "f[f[x, y], f[z]]");
    assert_eq!(r.args()[0], left); // the same node
    let g = Expr::symbol("g");
    assert_eq!(
        format!("{:?}", e.replace_part(&[0], g).unwrap()),
        "g[f[x, y], f[y]]"
    );
    assert_eq!(e.replace_part(&[], z.clone()), Some(z.clone()));
    assert_eq!(e.replace_part(&[3], z.clone()), None);
    assert_eq!(e.replace_part(&[1, 1, 1], z), None);
}

#[test]
fn positions_in_pre_order() {
    let [x] = symbols(["x"]);
    let e = Expr::apply(
        &f(std::slice::from_ref(&x)),
        &[x.clone(), f(std::slice::from_ref(&x))],
    );
    assert_eq!(e.positions(|s| *s == x), [vec![0, 1], vec![1], vec![2, 1]]);
    assert_eq!(
        e.positions(|s| s.head().is_some()),
        [vec![], vec![0], vec![2]]
    );
    let sym_f = Expr::symbol("f");
    assert_eq!(e.positions(|s| *s == sym_f), [vec![0, 0], vec![2, 0]]);
}

#[test]
fn symbols_once_in_canonical_order() {
    let [x, b, a_upper] = symbols(["x", "b", "A"]);
    let e = (&x + &b) * Expr::call(Symbol::SIN, [&a_upper]) + &x;
    let names: Vec<&str> = e.symbols().iter().map(|s| s.name()).collect();
    assert_eq!(names, ["A", "b", "Plus", "Sin", "Times", "x"]);
}

#[test]
fn substitution_is_simultaneous_and_top_down() {
    let [x, y] = symbols(["x", "y"]);
    let e = f(&[x.clone(), y.clone(), f(std::slice::from_ref(&x))]);
    let swapped = e.substitute(&[(x.clone(), y.clone()), (y.clone(), x.clone())]);
    assert_eq!(format!("{swapped:?}"), "f[y, x, f[y]]");
    // Results are not substituted into, and the first rule wins.
    let r = e.substitute(&[
        (x.clone(), f(std::slice::from_ref(&x))),
        (x.clone(), y.clone()),
    ]);
    assert_eq!(format!("{r:?}"), "f[f[x], y, f[f[x]]]");
    // Whole subexpressions, heads and the root.
    assert_eq!(
        format!(
            "{:?}",
            e.substitute(&[(f(std::slice::from_ref(&x)), y.clone())])
        ),
        "f[x, y, y]"
    );
    assert_eq!(
        format!(
            "{:?}",
            e.substitute(&[(Expr::symbol("f"), Expr::symbol("g"))])
        ),
        "g[x, y, g[x]]"
    );
    assert_eq!(e.substitute(&[(e.clone(), x.clone())]), x);
    // Nothing to do: the same handle.
    assert_eq!(e.substitute(&[]), e);
    assert_eq!(e.substitute(&[(Expr::symbol("z"), x)]), e);
}

#[test]
fn map_args_keeps_the_handle_when_unchanged() {
    let [x, y] = symbols(["x", "y"]);
    let e = f(&[x.clone(), y.clone()]);
    assert_eq!(e.map_args(Clone::clone), e);
    assert_eq!(
        format!("{:?}", e.map_args(|a| f(std::slice::from_ref(a)))),
        "f[f[x], f[y]]"
    );
    assert_eq!(x.map_args(|_| y.clone()), x);
}

/// The done-when of CAS-25: on a DAG with `2^n` occurrences, the closures are
/// called once per distinct subexpression and the result keeps the sharing.
#[test]
fn shared_dags_take_linear_time() {
    let n = 100;
    let e = doubling(n);
    let (x, y) = (Expr::symbol("x"), Expr::symbol("y"));

    let mut calls = 0;
    let r = e.replace(|s| {
        calls += 1;
        // Fail fast rather than run for 2^n calls.
        assert!(calls <= n + 2, "called once per occurrence");
        (*s == x).then(|| y.clone())
    });
    assert_eq!(calls, n + 2);
    // Never `assert_eq!` a DAG: printing it is exponential.
    assert!(r.args()[0] == r.args()[1]);
    assert_eq!(r.part(&[1; 100]), Some(&y));

    let s = e.substitute(&[(x.clone(), y.clone())]);
    assert!(s == r);

    let mut calls = 0;
    let m = e.map_bottom_up(|s| {
        calls += 1;
        assert!(calls <= n + 2, "called once per occurrence");
        s.clone()
    });
    assert_eq!(calls, n + 2);
    assert!(m == e);

    assert_eq!(e.distinct().take(n + 3).count(), n + 2);
    assert!(e.contains(&x) && e.free_of(&y));
    assert_eq!(e.symbols(), [Symbol::new("f"), Symbol::new("x")]);
}

#[test]
fn deep_expressions_use_constant_stack() {
    let n = 100_000;
    let x = Expr::symbol("x");
    let mut e = x.clone();
    for _ in 0..n {
        e = f(&[e]);
    }
    let y = Expr::symbol("y");
    assert_eq!(e.pre_order().count(), 2 * n + 1);
    assert_eq!(e.post_order().count(), 2 * n + 1);
    assert_eq!(e.distinct().count(), n + 2);
    assert_eq!(e.positions(|s| *s == x), [vec![1; n]]);
    assert!(e.contains(&x));
    let r = e.substitute(&[(x.clone(), y.clone())]);
    assert_eq!(r.part(&vec![1; n]), Some(&y));
    // `assert!`, not `assert_eq!`: a wrong result can take quadratic space
    // to print.
    assert!(e.replace(|_| None) == e);
    assert!(e.map_bottom_up(Clone::clone) == e);
    assert!(e.replace_part(&vec![1; n], y.clone()) == Some(r));
}

/// Trees over a few symbols and heads, with shared subtrees likely.
fn tree() -> impl Strategy<Value = Expr> {
    let leaf = prop_oneof![
        prop::sample::select(vec!["x", "y", "z", "f"]).prop_map(Expr::symbol),
        (-1i64..2).prop_map(Expr::from),
    ];
    leaf.prop_recursive(4, 32, 3, |inner| {
        (inner.clone(), prop::collection::vec(inner, 0..3)).prop_map(|(h, a)| Expr::apply(&h, &a))
    })
}

/// Children: the head, then the arguments.
fn children(e: &Expr) -> Vec<&Expr> {
    e.head().into_iter().chain(e.args()).collect()
}

fn model_pre<'a>(e: &'a Expr, out: &mut Vec<&'a Expr>) {
    out.push(e);
    for c in children(e) {
        model_pre(c, out);
    }
}

fn model_post<'a>(e: &'a Expr, out: &mut Vec<&'a Expr>) {
    for c in children(e) {
        model_post(c, out);
    }
    out.push(e);
}

fn model_positions(
    e: &Expr,
    at: &mut Vec<usize>,
    pred: &dyn Fn(&Expr) -> bool,
    out: &mut Vec<Vec<usize>>,
) {
    if pred(e) {
        out.push(at.clone());
    }
    for (i, c) in children(e).into_iter().enumerate() {
        at.push(i);
        model_positions(c, at, pred, out);
        at.pop();
    }
}

fn model_replace(e: &Expr, f: &dyn Fn(&Expr) -> Option<Expr>) -> Expr {
    if let Some(r) = f(e) {
        return r;
    }
    match e.head() {
        Some(h) => {
            let args: Vec<Expr> = e.args().iter().map(|a| model_replace(a, f)).collect();
            Expr::apply(&model_replace(h, f), &args)
        }
        None => e.clone(),
    }
}

fn model_bottom_up(e: &Expr, f: &dyn Fn(&Expr) -> Expr) -> Expr {
    let rebuilt = match e.head() {
        Some(h) => {
            let args: Vec<Expr> = e.args().iter().map(|a| model_bottom_up(a, f)).collect();
            Expr::apply(&model_bottom_up(h, f), &args)
        }
        None => e.clone(),
    };
    f(&rebuilt)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(1000))]

    #[test]
    fn orders_match_the_tree(e in tree()) {
        let (mut pre, mut post) = (Vec::new(), Vec::new());
        model_pre(&e, &mut pre);
        model_post(&e, &mut post);
        prop_assert_eq!(e.pre_order().collect::<Vec<_>>(), pre.clone());
        prop_assert_eq!(e.post_order().collect::<Vec<_>>(), post.clone());

        // Distinct: every subexpression once, parts before their parents.
        let distinct: Vec<&Expr> = e.distinct().collect();
        let set: HashSet<&Expr> = distinct.iter().copied().collect();
        prop_assert_eq!(set.len(), distinct.len());
        prop_assert_eq!(&set, &pre.iter().copied().collect::<HashSet<_>>());
        for (i, d) in distinct.iter().enumerate() {
            for c in children(d) {
                prop_assert!(distinct[..i].contains(&c));
            }
        }
        prop_assert_eq!(distinct.last(), Some(&&e));
    }

    #[test]
    fn positions_and_parts_agree(e in tree(), t in tree()) {
        let pred = |s: &Expr| s.args().len() == 1 || s.as_i64().is_some();
        let mut model = Vec::new();
        model_positions(&e, &mut Vec::new(), &pred, &mut model);
        let found = e.positions(pred);
        prop_assert_eq!(&found, &model);
        for p in &found {
            prop_assert!(pred(e.part(p).unwrap()));
            let r = e.replace_part(p, t.clone()).unwrap();
            prop_assert_eq!(r.part(p), Some(&t));
        }
    }

    #[test]
    fn queries_match_the_tree(e in tree(), sub in tree()) {
        let mut pre = Vec::new();
        model_pre(&e, &mut pre);
        prop_assert_eq!(e.contains(&sub), pre.contains(&&sub));
        for s in pre.iter().copied() {
            prop_assert!(e.contains(s));
        }
        let mut symbols: Vec<Symbol> = pre.iter().filter_map(|s| s.as_symbol()).collect();
        symbols.sort_by_key(|&s| Expr::from(s));
        symbols.dedup();
        prop_assert_eq!(e.symbols(), symbols);
    }

    #[test]
    fn rewriting_matches_the_recursive_model(e in tree(), from in tree(), to in tree()) {
        let x = Expr::symbol("x");
        let rule = |s: &Expr| {
            if *s == from {
                Some(to.clone())
            } else if *s == x {
                Some(f(&[]))
            } else {
                None
            }
        };
        prop_assert_eq!(e.replace(rule), model_replace(&e, &rule));
        let rules = [(from.clone(), to.clone()), (x.clone(), f(&[]))];
        prop_assert_eq!(e.substitute(&rules), model_replace(&e, &rule));
        prop_assert_eq!(&e.replace(|_| None), &e);

        let wrap = |s: &Expr| if s.as_i64().is_some() { f(std::slice::from_ref(s)) } else { s.clone() };
        prop_assert_eq!(e.map_bottom_up(wrap), model_bottom_up(&e, &wrap));
    }
}
