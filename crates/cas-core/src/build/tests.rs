use proptest::prelude::*;

use super::{cos, exp, list, log, sin, symbols, tan};
use crate::{Expr, Symbol};

fn full(e: &Expr) -> String {
    format!("{e:?}")
}

#[test]
fn operators_on_every_operand_type() {
    let [x, y] = symbols(["x", "y"]);
    let s = Symbol::new("x");
    let expected = "Plus[2, x]";
    assert_eq!(full(&(x.clone() + 2)), expected);
    assert_eq!(full(&(&x + 2)), expected);
    assert_eq!(full(&(s + 2)), expected);
    assert_eq!(full(&(2 + x.clone())), expected);
    assert_eq!(full(&(2 + &x)), expected);
    assert_eq!(full(&(2 + s)), expected);
    assert_eq!(full(&(&x + &y)), "Plus[x, y]");
    assert_eq!(full(&(s * &y)), "Times[x, y]");
    assert_eq!(full(&(2 - s)), "Plus[2, Times[-1, x]]");
    assert_eq!(full(&(s - &y)), "Plus[x, Times[-1, y]]");
    assert_eq!(full(&(1 / &x)), "Times[1, Power[x, -1]]");
    assert_eq!(full(&s.pow(&y)), "Power[x, y]");
}

#[test]
fn assignment_operators() {
    let [x, y] = symbols(["x", "y"]);
    let mut e = x.clone();
    e += &y;
    e *= 3;
    assert_eq!(full(&e), "Times[3, Plus[x, y]]");
    e -= 1;
    e /= &x;
    assert_eq!(
        full(&e),
        "Times[Power[x, -1], Plus[-1, Times[3, Plus[x, y]]]]"
    );
}

#[test]
fn negation() {
    let [x] = symbols(["x"]);
    assert_eq!(-Expr::from(3), Expr::from(-3));
    assert_eq!(
        full(&-Expr::from(i64::MIN)),
        format!("Times[{}, -1]", i64::MIN)
    );
    assert_eq!(full(&-&x), "Times[-1, x]");
    assert_eq!(full(&-Symbol::new("x")), "Times[-1, x]");
    // Not arithmetic: the coefficient is not folded.
    assert_eq!(full(&-(2 * &x)), "Times[-1, 2, x]");
}

#[test]
fn nothing_is_evaluated() {
    let [x] = symbols(["x"]);
    assert_eq!(full(&(&x + &x)), "Plus[x, x]");
    assert_eq!(full(&(Expr::from(1) + 2)), "Plus[1, 2]");
    assert_eq!(full(&(&x - &x)), "Plus[Times[-1, x], x]");
    assert_eq!(full(&x.pow(1)), "Power[x, 1]");
    // `Power` is not `Flat`.
    assert_eq!(full(&x.pow(2).pow(3)), "Power[Power[x, 2], 3]");
}

#[test]
fn only_builtin_heads_are_normalised() {
    let [x, y] = symbols(["x", "y"]);
    let f = Symbol::new("f");
    let nested = Expr::call(f, [f.call([&y]), x.clone()]);
    assert_eq!(full(&nested), "f[f[y], x]");
    // `List` is a builtin without `Flat` or `Orderless`.
    assert_eq!(full(&list([&y, &x])), "List[y, x]");
    assert_eq!(full(&list([list([&y])])), "List[List[y]]");
    // Only `Orderless` sorts: `Power` is a `Listable` builtin.
    assert_eq!(full(&Expr::call(Symbol::POWER, [&y, &x])), "Power[y, x]");
    // A compound head is never normalised.
    let g = Expr::apply(&Symbol::PLUS.into(), &[]);
    assert_eq!(full(&Expr::call(&g, [&y, &x])), "Plus[][y, x]");
}

#[test]
fn flattening_at_any_depth() {
    let [a, b, c, d] = symbols(["a", "b", "c", "d"]);
    let plus = |args: &[Expr]| Expr::apply(&Symbol::PLUS.into(), args);
    let raw = plus(&[
        d.clone(),
        plus(&[plus(&[c.clone(), plus(&[])]), b.clone()]),
        Expr::apply(&Symbol::TIMES.into(), std::slice::from_ref(&a)),
    ]);
    assert_eq!(
        full(&Expr::call(Symbol::PLUS, raw.args())),
        "Plus[Times[a], b, c, d]"
    );
}

#[test]
fn deep_flattening_uses_constant_stack() {
    let [x] = symbols(["x"]);
    let plus: Expr = Symbol::PLUS.into();
    let mut e = x.clone();
    for i in 0..100_000 {
        e = Expr::apply(&plus, &[e, Expr::from(i)]);
    }
    let flat = Expr::call(Symbol::PLUS, [e]);
    assert_eq!(flat.args().len(), 100_001);
    assert_eq!(flat.args()[0], Expr::from(0));
    assert_eq!(flat.args()[100_000], x);
}

#[test]
fn sums_and_products_of_iterators() {
    let [x, y] = symbols(["x", "y"]);
    let none: [Expr; 0] = [];
    assert_eq!(none.iter().sum::<Expr>(), Expr::from(0));
    assert_eq!(none.iter().product::<Expr>(), Expr::from(1));
    assert_eq!([x.clone()].into_iter().sum::<Expr>(), x);
    let terms = [y.pow(2), &x * &y, x.pow(2)];
    assert_eq!(
        full(&terms.iter().sum()),
        "Plus[Power[x, 2], Times[x, y], Power[y, 2]]"
    );
    assert_eq!(
        full(&[&x + 1, &y + 2].into_iter().product()),
        "Times[Plus[1, x], Plus[2, y]]"
    );
}

#[test]
fn builtin_functions() {
    let [x] = symbols(["x"]);
    let e = [exp(&x), log(&x), sin(&x), cos(&x), tan(&x)];
    assert_eq!(
        full(&list(e)),
        "List[Exp[x], Log[x], Sin[x], Cos[x], Tan[x]]"
    );
    assert_eq!(full(&sin(2)), "Sin[2]");
}

/// A tree over a few names; heads include `Plus` and `Times` so that
/// flattening has nested applications to splice.
fn tree() -> impl Strategy<Value = Expr> {
    let leaf = prop_oneof![
        prop::sample::select(vec!["a", "b", "c"]).prop_map(Expr::symbol),
        (-2i64..3).prop_map(Expr::from),
    ];
    leaf.prop_recursive(4, 32, 4, |inner| {
        let head = prop::sample::select(vec![Symbol::PLUS, Symbol::TIMES, Symbol::new("f")]);
        (head, prop::collection::vec(inner, 0..4)).prop_map(|(h, a)| Expr::apply(&h.into(), &a))
    })
}

/// The flattened arguments of `head`, unsorted, by plain recursion.
fn model_flatten(head: &Expr, args: &[Expr], out: &mut Vec<Expr>) {
    for a in args {
        if a.head() == Some(head) {
            model_flatten(head, a.args(), out);
        } else {
            out.push(a.clone());
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(1000))]

    #[test]
    fn plus_and_times_flatten_and_sort(args in prop::collection::vec(tree(), 0..6)) {
        for head in [Symbol::PLUS, Symbol::TIMES] {
            let h = Expr::from(head);
            let mut model = Vec::new();
            model_flatten(&h, &args, &mut model);
            model.sort();
            let built = Expr::call(head, &args);
            prop_assert_eq!(built.head(), Some(&h));
            prop_assert_eq!(built.args(), &model[..]);
            // Idempotent, and independent of the argument order.
            prop_assert_eq!(&Expr::call(head, built.args()), &built);
            let reversed: Vec<Expr> = args.iter().rev().cloned().collect();
            prop_assert_eq!(&Expr::call(head, reversed), &built);
        }
    }

    #[test]
    fn operators_agree_with_call(a in tree(), b in tree()) {
        prop_assert_eq!(&a + &b, Expr::call(Symbol::PLUS, [&a, &b]));
        prop_assert_eq!(&a * &b, Expr::call(Symbol::TIMES, [&a, &b]));
        prop_assert_eq!(&a + &b, &b + &a);
        prop_assert_eq!((&a + &b) + &a, [&a, &b, &a].into_iter().sum::<Expr>());
    }
}
