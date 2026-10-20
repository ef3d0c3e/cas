//! Tests of the builtins' native functions, through evaluation.

use std::sync::atomic::{AtomicUsize, Ordering};

use cas_num::{Complex, Integer, Rational};
use proptest::prelude::*;

use crate::{Context, Expr, Symbol, View};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, result)` in full form, each evaluated in a new context.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        match ctx.eval(&read(input)) {
            Ok(found) if format!("{found:?}") == *expected => {}
            found => wrong.push(format!("{input}: {found:?}, expected {expected}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn sums() {
    check(&[
        ("Plus[1, 2, 3]", "6"),
        ("Plus[1/2, 1/3]", "5/6"),
        ("Plus[1, Complex[0, 1]]", "Complex[1, 1]"),
        ("Plus[Complex[1, 1], Complex[0, -1]]", "1"),
        ("Plus[x, 1, 2]", "Plus[3, x]"),
        ("Plus[x, x]", "Times[2, x]"),
        ("Plus[x, Times[2, x], y]", "Plus[Times[3, x], y]"),
        ("Plus[x, Times[-1, x]]", "0"),
        ("Plus[Times[x, y], Times[2, x, y]]", "Times[3, x, y]"),
        ("Plus[Times[1/2, x], Times[1/2, x]]", "x"),
        (
            "Plus[Times[Complex[0, 1], x], Times[Complex[0, -1], x], y]",
            "y",
        ),
        ("Plus[0, x]", "x"),
        ("Plus[x]", "x"),
        ("Plus[]", "0"),
        ("Plus[f[x], f[x], g[x]]", "Plus[Times[2, f[x]], g[x]]"),
        (
            "Plus[Power[2, 1/2], Power[8, 1/2]]",
            "Times[3, Power[2, 1/2]]",
        ),
        ("Plus[List[1, 2], List[3, 4]]", "List[4, 6]"),
    ]);
}

#[test]
fn products() {
    check(&[
        ("Times[2, 3]", "6"),
        ("Times[2, x, 3]", "Times[6, x]"),
        ("Times[1/2, 2]", "1"),
        ("Times[x, x]", "Power[x, 2]"),
        ("Times[x, Power[x, -1]]", "1"),
        ("Times[x, Power[x, 2], y]", "Times[Power[x, 3], y]"),
        ("Times[Power[x, a], Power[x, b]]", "Power[x, Plus[a, b]]"),
        ("Times[0, x]", "0"),
        ("Times[1, x]", "x"),
        ("Times[x]", "x"),
        ("Times[]", "1"),
        ("Times[Complex[0, 1], Complex[0, 1]]", "-1"),
        ("Times[2, Times[x, 3]]", "Times[6, x]"),
        // Powers of numbers with one exponent multiply.
        ("Times[Power[2, 1/2], Power[3, 1/2]]", "Power[6, 1/2]"),
        ("Times[Power[2, 1/2], Power[2, 1/2]]", "2"),
        (
            "Times[Power[2, 1/2], Power[6, 1/2]]",
            "Times[2, Power[3, 1/2]]",
        ),
        (
            "Times[Power[2, 1/3], Power[3, 1/2]]",
            "Times[Power[2, 1/3], Power[3, 1/2]]",
        ),
        ("Times[2, Power[2, 1/2]]", "Times[2, Power[2, 1/2]]"),
        ("Times[Power[6, 1/2], Power[2, -1/2]]", "Power[3, 1/2]"),
        // Lists are left to the products of linear algebra (R-10).
        // A scalar scales a list (CAS-62).
        ("Times[2, 3, List[1, 2]]", "List[6, 12]"),
        ("Times[0, List[1, 2]]", "List[0, 0]"),
    ]);
}

#[test]
fn powers() {
    check(&[
        ("Power[2, 10]", "1024"),
        ("Power[2, -2]", "1/4"),
        ("Power[2/3, 3]", "8/27"),
        ("Power[-2, 3]", "-8"),
        ("Power[Complex[0, 1], 2]", "-1"),
        ("Power[Complex[1, 1], -1]", "Complex[1/2, -1/2]"),
        ("Power[x, 0]", "1"),
        ("Power[x, 1]", "x"),
        ("Power[1, x]", "1"),
        ("Power[Power[x, 2], 3]", "Power[x, 6]"),
        ("Power[Power[x, 1/2], 2]", "x"),
        // In the Real domain, the default (CAS-58): `x` is real.
        ("Power[Power[x, 2], 1/2]", "Abs[x]"),
        ("Power[Power[x, a], -1]", "Power[x, Times[-1, a]]"),
        ("Power[Times[2, x], 2]", "Times[4, Power[x, 2]]"),
        (
            "Power[Times[x, y], -1]",
            "Times[Power[x, -1], Power[y, -1]]",
        ),
        ("Power[Times[x, y], 1/2]", "Power[Times[x, y], 1/2]"),
        ("Power[x, y]", "Power[x, y]"),
        // Too large to hold: left unevaluated.
        ("Power[2, 100000000]", "Power[2, 100000000]"),
        ("Sqrt[4]", "2"),
        ("Sqrt[x]", "Power[x, 1/2]"),
        ("Sqrt[-1]", "Complex[0, 1]"),
        ("I", "Complex[0, 1]"),
        ("Times[I, I]", "-1"),
    ]);
}

#[test]
fn radicals() {
    check(&[
        ("Power[12, 1/2]", "Times[2, Power[3, 1/2]]"),
        ("Power[8, 2/3]", "4"),
        ("Power[72, 1/3]", "Times[2, Power[3, 2/3]]"),
        ("Power[9, 1/3]", "Power[3, 2/3]"),
        ("Power[2, 1/2]", "Power[2, 1/2]"),
        ("Power[2, 2/3]", "Power[2, 2/3]"),
        ("Power[12, 2/3]", "Times[2, Power[2, 1/3], Power[3, 2/3]]"),
        // The integer part of the exponent comes out.
        ("Power[2, 3/2]", "Times[2, Power[2, 1/2]]"),
        ("Power[2, -1/2]", "Times[1/2, Power[2, 1/2]]"),
        ("Power[12, 3/2]", "Times[24, Power[3, 1/2]]"),
        // A denominator goes into the numerator.
        ("Power[2/3, 1/2]", "Times[1/3, Power[6, 1/2]]"),
        ("Power[3/4, 1/2]", "Times[1/2, Power[3, 1/2]]"),
        ("Power[1/4, 1/3]", "Times[1/2, Power[2, 1/3]]"),
        (
            "Power[5/3, -2/3]",
            "Times[1/5, Power[3, 2/3], Power[5, 1/3]]",
        ),
        // Square roots of negatives are imaginary; in the Real domain, odd
        // roots are real (O7, CAS-58).
        ("Power[-4, 1/2]", "Complex[0, 2]"),
        ("Power[-4, 3/2]", "Complex[0, -8]"),
        ("Power[-1, 1/2]", "Complex[0, 1]"),
        ("Power[-2, 1/2]", "Times[Complex[0, 1], Power[2, 1/2]]"),
        ("Power[-1/4, 1/2]", "Complex[0, 1/2]"),
        ("Power[-8, 1/3]", "-2"),
        ("Power[Power[2, 1/2], 1/2]", "Power[2, 1/4]"),
        ("Power[0, 1/2]", "0"),
        ("Power[1, 1/3]", "1"),
        // A prime above 2^16, squared: found as a perfect power.
        ("Power[4295098369, 1/2]", "65537"),
        ("Power[Complex[0, 1], 1/2]", "Power[Complex[0, 1], 1/2]"),
    ]);
}

#[test]
fn infinities_and_singular_cases() {
    check(&[
        ("DirectedInfinity[1]", "Infinity"),
        ("DirectedInfinity[]", "ComplexInfinity"),
        ("DirectedInfinity[0]", "ComplexInfinity"),
        ("DirectedInfinity[-3]", "DirectedInfinity[-1]"),
        (
            "DirectedInfinity[Complex[3, 4]]",
            "DirectedInfinity[Complex[3/5, 4/5]]",
        ),
        (
            "DirectedInfinity[Complex[1, 1]]",
            "DirectedInfinity[Complex[1, 1]]",
        ),
        ("DirectedInfinity[x]", "DirectedInfinity[x]"),
        ("Plus[Infinity, 1]", "Infinity"),
        ("Plus[Infinity, x]", "Infinity"),
        ("Plus[Infinity, Infinity]", "Infinity"),
        ("Plus[Infinity, DirectedInfinity[-1]]", "Indeterminate"),
        ("Plus[Infinity, Times[-1, Infinity]]", "Indeterminate"),
        ("Plus[ComplexInfinity, 1]", "ComplexInfinity"),
        ("Plus[ComplexInfinity, ComplexInfinity]", "Indeterminate"),
        ("Plus[ComplexInfinity, Infinity]", "Indeterminate"),
        ("Plus[Indeterminate, 1]", "Indeterminate"),
        ("Times[2, Infinity]", "Infinity"),
        ("Times[-2, Infinity]", "DirectedInfinity[-1]"),
        (
            "Times[Complex[0, 1], Infinity]",
            "DirectedInfinity[Complex[0, 1]]",
        ),
        ("Times[Infinity, Infinity]", "Infinity"),
        (
            "Times[DirectedInfinity[-1], DirectedInfinity[-1]]",
            "Infinity",
        ),
        ("Times[x, Infinity]", "Times[Infinity, x]"),
        ("Times[-2, x, Infinity]", "Times[x, DirectedInfinity[-1]]"),
        ("Times[ComplexInfinity, x]", "ComplexInfinity"),
        ("Times[ComplexInfinity, Infinity]", "ComplexInfinity"),
        ("Times[Indeterminate, x]", "Indeterminate"),
        // The singular cases.
        ("Times[0, Infinity]", "Indeterminate"),
        ("Times[0, ComplexInfinity]", "Indeterminate"),
        ("Power[0, 0]", "Indeterminate"),
        ("Power[0, -1]", "ComplexInfinity"),
        ("Power[0, -3/2]", "ComplexInfinity"),
        ("Times[1, Power[0, -1]]", "ComplexInfinity"),
        ("Power[Infinity, 0]", "Indeterminate"),
        ("Times[Infinity, Power[Infinity, -1]]", "Indeterminate"),
        ("Plus[Infinity, Times[-1, Infinity]]", "Indeterminate"),
        // Powers of infinities.
        ("Power[Infinity, 2]", "Infinity"),
        ("Power[Infinity, 1/2]", "Infinity"),
        ("Power[Infinity, -1]", "0"),
        ("Power[DirectedInfinity[-1], 3]", "DirectedInfinity[-1]"),
        ("Power[DirectedInfinity[-1], 2]", "Infinity"),
        ("Power[ComplexInfinity, 2]", "ComplexInfinity"),
        ("Power[ComplexInfinity, -2]", "0"),
        ("Power[Indeterminate, 2]", "Indeterminate"),
        ("Power[2, Indeterminate]", "Indeterminate"),
        ("Power[x, Infinity]", "Power[x, Infinity]"),
    ]);
}

#[test]
fn comparisons_and_predicates() {
    check(&[
        ("Equal[1, 1]", "True"),
        ("Equal[1, 2]", "False"),
        ("Equal[x, x]", "True"),
        ("Equal[x, y]", "Equal[x, y]"),
        (r#"Equal["a", "b"]"#, "False"),
        ("Equal[True, False]", "False"),
        ("Equal[1, 1, 1]", "True"),
        ("Equal[1, 1, 2]", "False"),
        ("Equal[x, x, y]", "Equal[x, x, y]"),
        ("Equal[Infinity, DirectedInfinity[-1]]", "False"),
        ("Equal[Plus[1, 1], 2]", "True"),
        ("Unequal[1, 2]", "True"),
        ("Unequal[2, 2]", "False"),
        ("Unequal[x, y]", "Unequal[x, y]"),
        ("Less[1, 2]", "True"),
        ("Less[2, 1]", "False"),
        ("Less[1, 2, 3]", "True"),
        ("Less[1, 3, 2]", "False"),
        ("Less[1, x]", "Less[1, x]"),
        ("Less[1, Infinity]", "True"),
        ("Less[DirectedInfinity[-1], -100]", "True"),
        ("Less[Infinity, Infinity]", "False"),
        ("LessEqual[1, 1]", "True"),
        ("Greater[1/2, 1/3]", "True"),
        ("GreaterEqual[1, 2]", "False"),
        ("Less[Complex[0, 1], 2]", "Less[Complex[0, 1], 2]"),
        ("SameQ[x, x]", "True"),
        ("SameQ[x, y]", "False"),
        ("SameQ[1, 1, 1]", "True"),
        ("FreeQ[f[x, y], x]", "False"),
        ("FreeQ[f[x, y], z]", "True"),
        ("FreeQ[f[x], f]", "False"),
        ("IntegerQ[2]", "True"),
        ("IntegerQ[1/2]", "False"),
        ("IntegerQ[x]", "False"),
        ("NumberQ[Complex[0, 1]]", "True"),
        ("NumberQ[Pi]", "False"),
    ]);
}

#[test]
fn logic_and_conditionals() {
    check(&[
        ("And[True, True]", "True"),
        ("And[True, False]", "False"),
        ("And[x, True]", "x"),
        ("And[x, y]", "And[x, y]"),
        ("And[]", "True"),
        ("And[False, x]", "False"),
        ("And[x, False]", "False"),
        ("And[Less[1, 2], Less[2, 3]]", "True"),
        ("Or[False, x]", "x"),
        ("Or[True, x]", "True"),
        ("Or[]", "False"),
        ("Or[x, y]", "Or[x, y]"),
        ("Not[True]", "False"),
        ("Not[False]", "True"),
        ("Not[Not[x]]", "x"),
        ("Not[x]", "Not[x]"),
        ("If[True, a, b]", "a"),
        ("If[False, a, b]", "b"),
        ("If[False, a]", "Null"),
        ("If[x, a, b]", "If[x, a, b]"),
        ("If[Less[1, 2], Plus[1, 1], b]", "2"),
    ]);
}

#[test]
fn short_circuits_skip_what_they_do_not_need() {
    static CALLS: AtomicUsize = AtomicUsize::new(0);
    let counted = crate::eval::tests::define_native("logic_counted", |_, _| {
        CALLS.fetch_add(1, Ordering::Relaxed);
        Ok(None)
    });
    let ctx = Context::new();
    let counted = Expr::apply(&counted.into(), &[]);
    let and = Expr::apply(
        &Symbol::AND.into(),
        &[Symbol::FALSE.into(), counted.clone()],
    );
    let or = Expr::apply(&Symbol::OR.into(), &[Symbol::TRUE.into(), counted.clone()]);
    let branch = Expr::apply(
        &Symbol::IF.into(),
        &[Symbol::TRUE.into(), Expr::from(1), counted.clone()],
    );
    for e in [and, or, branch] {
        ctx.eval(&e).unwrap();
    }
    assert_eq!(CALLS.load(Ordering::Relaxed), 0);
    let and = Expr::apply(&Symbol::AND.into(), &[Symbol::TRUE.into(), counted]);
    ctx.eval(&and).unwrap();
    assert_eq!(CALLS.load(Ordering::Relaxed), 1);
}

#[test]
fn structure() {
    check(&[
        ("Length[List[1, 2, 3]]", "3"),
        ("Length[x]", "0"),
        ("Length[f[]]", "0"),
        ("Part[List[a, b, c], 1]", "a"),
        ("Part[List[a, b, c], -1]", "c"),
        ("Part[f[a, b], 0]", "f"),
        ("Part[List[List[1, 2], List[3, 4]], 2, 1]", "3"),
        ("Part[List[a], 5]", "Part[List[a], 5]"),
        ("Part[List[a], -2]", "Part[List[a], -2]"),
        ("Part[x, 1]", "Part[x, 1]"),
        ("Part[x, 0]", "Part[x, 0]"),
        ("Part[List[a, b], x]", "Part[List[a, b], x]"),
        ("Map[f, List[a, b]]", "List[f[a], f[b]]"),
        ("Map[f, x]", "x"),
        ("Map[f, g[a]]", "g[f[a]]"),
        ("Apply[f, List[a, b]]", "f[a, b]"),
        ("Apply[Plus, List[1, 2]]", "3"),
        ("Apply[f, x]", "x"),
    ]);
}

/// Approximate value of an expression made of numbers, sums, products and
/// powers, as `(re, im)`.
fn approx(e: &Expr) -> Option<(f64, f64)> {
    let real = |q: &Rational| {
        q.to_string().parse::<f64>().ok().or_else(|| {
            let (n, d) = (q.numer().to_string(), q.denom().to_string());
            Some(n.parse::<f64>().ok()? / d.parse::<f64>().ok()?)
        })
    };
    match e.view() {
        View::Complex(z) => Some((real(z.re())?, real(z.im())?)),
        View::Apply { head, args } => {
            let values: Option<Vec<(f64, f64)>> = args.iter().map(approx).collect();
            let values = values?;
            match head.as_symbol()? {
                Symbol::PLUS => Some(
                    values
                        .iter()
                        .fold((0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1)),
                ),
                Symbol::TIMES => Some(values.iter().fold((1.0, 0.0), |a, b| {
                    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
                })),
                Symbol::POWER => {
                    let [(b, bi), (n, ni)] = values.as_slice() else {
                        return None;
                    };
                    (*bi == 0.0 && *ni == 0.0 && *b > 0.0).then(|| (b.powf(*n), 0.0))
                }
                _ => None,
            }
        }
        _ => e.to_rational().and_then(|q| real(&q)).map(|x| (x, 0.0)),
    }
}

fn close(a: (f64, f64), b: (f64, f64)) -> bool {
    let scale = 1.0 + b.0.abs() + b.1.abs();
    (a.0 - b.0).abs() <= 1e-9 * scale && (a.1 - b.1).abs() <= 1e-9 * scale
}

fn rational() -> impl Strategy<Value = Rational> {
    (-50i64..50, 1i64..20).prop_map(|(n, d)| Rational::new(n.into(), d.into()))
}

fn number() -> impl Strategy<Value = Expr> {
    prop_oneof![
        rational().prop_map(Expr::from),
        (rational(), rational()).prop_map(|(re, im)| Expr::from(Complex::new(re, im))),
        any::<i64>().prop_map(Expr::from),
    ]
}

fn value(e: &Expr) -> Complex {
    e.to_complex()
        .unwrap_or_else(|| panic!("{e:?} is not a number"))
}

/// Terms for sums and products: numbers, symbols, products and powers.
fn term() -> impl Strategy<Value = Expr> {
    let atom = prop_oneof![
        rational().prop_map(Expr::from),
        prop::sample::select(vec!["x", "y", "z"]).prop_map(|n| Expr::from(Symbol::new(n))),
        Just(Expr::from(Symbol::INFINITY)),
    ];
    atom.prop_recursive(2, 12, 3, |inner| {
        prop_oneof![
            prop::collection::vec(inner.clone(), 2..4)
                .prop_map(|a| Expr::apply(&Symbol::TIMES.into(), &a)),
            (inner, -2i64..3).prop_map(|(b, n)| Expr::apply(&Symbol::POWER.into(), &[b, n.into()])),
        ]
    })
}

proptest! {
    /// Sums, products and integer powers of numbers agree with cas-num.
    #[test]
    fn arithmetic_matches_cas_num(a in number(), b in number(), c in number(), n in -6i64..7) {
        let ctx = Context::new();
        let sum = ctx.eval(&Expr::apply(&Symbol::PLUS.into(), &[a.clone(), b.clone(), c.clone()])).unwrap();
        prop_assert_eq!(value(&sum), &(&value(&a) + &value(&b)) + &value(&c));
        let product = ctx.eval(&Expr::apply(&Symbol::TIMES.into(), &[a.clone(), b.clone()])).unwrap();
        prop_assert_eq!(value(&product), &value(&a) * &value(&b));
        // `0^0` and `0^-n` are not numbers.
        if !value(&a).is_zero() || n > 0 {
            let power = ctx.eval(&Expr::apply(&Symbol::POWER.into(), &[a.clone(), n.into()])).unwrap();
            prop_assert_eq!(value(&power), value(&a).pow(n));
        }
    }

    /// Rational powers of rationals evaluate to an equal value, in a form
    /// that evaluates to itself.
    #[test]
    fn radicals_keep_their_value(
        b in (1i64..200, 1i64..50).prop_map(|(n, d)| Rational::new(n.into(), d.into())),
        p in -7i64..8,
        q in 2i64..7,
    ) {
        let ctx = Context::new();
        let r = Rational::new(p.into(), q.into());
        let e = Expr::apply(&Symbol::POWER.into(), &[Expr::from(b.clone()), Expr::from(r.clone())]);
        let result = ctx.eval(&e).unwrap();
        prop_assert!(ctx.eval(&result).unwrap() == result, "{:?} gave {:?}", e, result);
        let expected = approx(&e).unwrap();
        let found = approx(&result).unwrap_or_else(|| panic!("{result:?}"));
        prop_assert!(close(found, expected), "{:?} gave {:?}: {:?} against {:?}", e, result, found, expected);
    }

    /// The order of the arguments of `Plus` and `Times` changes nothing.
    #[test]
    fn argument_order_changes_nothing(
        args in prop::collection::vec(term(), 0..5),
        head in prop::sample::select(vec![Symbol::PLUS, Symbol::TIMES]),
    ) {
        let ctx = Context::new();
        let mut reversed = args.clone();
        reversed.reverse();
        let forward = ctx.eval(&Expr::apply(&head.into(), &args));
        let backward = ctx.eval(&Expr::apply(&head.into(), &reversed));
        prop_assert!(forward == backward, "{:?} against {:?}", forward, backward);
    }

    /// Integer bases: radicals of integers are exact when they are perfect
    /// powers.
    #[test]
    fn perfect_powers_come_out(root in 1u64..1000, k in 2u32..5) {
        let ctx = Context::new();
        let n = Integer::from(root).pow(u64::from(k));
        let e = Expr::apply(
            &Symbol::POWER.into(),
            &[Expr::from(n), Expr::from(Rational::new(1.into(), i64::from(k).into()))],
        );
        prop_assert_eq!(ctx.eval(&e).unwrap(), Expr::from(i64::try_from(root).unwrap()));
    }
}
