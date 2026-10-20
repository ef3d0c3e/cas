//! Tests of `D` as system rules (CAS-59): tables, the acceptance test,
//! random polynomials against the closed form, and the product, quotient,
//! power and chain rules against forward differentiation at random rational
//! points.

use proptest::prelude::*;

use crate::trace::{Origin, StepKind};
use crate::{Context, EvalOptions, Expr, Rational, Symbol, TraceLevel};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn apply(head: Symbol, args: &[Expr]) -> Expr {
    Expr::apply(&head.into(), args)
}

fn x() -> Expr {
    Symbol::new("x").into()
}

/// `D[f, x]`, evaluated.
fn d(ctx: &Context, f: &Expr) -> Expr {
    ctx.eval(&apply(Symbol::D, &[f.clone(), x()]))
        .unwrap_or_else(|e| panic!("D[{f:?}, x]: {e}"))
}

#[test]
fn derivatives() {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in [
        // The acceptance test (CAS-DOC-4 §8): no simplification asked for.
        (
            "D[Plus[Power[x, 2], Times[3, x]], x]",
            "Plus[3, Times[2, x]]",
        ),
        ("D[y, x]", "0"),
        ("D[5, x]", "0"),
        ("D[x, x]", "1"),
        ("D[Times[x, y], x]", "y"),
        ("D[Times[x, Sin[x]], x]", "Plus[Times[x, Cos[x]], Sin[x]]"),
        (
            "D[Times[x, Power[Plus[1, x], -1]], x]",
            "Plus[Times[-1, x, Power[Plus[1, x], -2]], Power[Plus[1, x], -1]]",
        ),
        ("D[Power[x, -1], x]", "Times[-1, Power[x, -2]]"),
        ("D[Sqrt[x], x]", "Times[1/2, Power[x, -1/2]]"),
        ("D[Power[x, x], x]", "Times[Power[x, x], Plus[1, Log[x]]]"),
        ("D[Power[E, x], x]", "Power[E, x]"),
        (
            "D[Power[E, Times[3, x]], x]",
            "Times[3, Power[E, Times[3, x]]]",
        ),
        // The table, and the chain rule.
        ("D[Sin[x], x]", "Cos[x]"),
        ("D[Sin[Power[x, 2]], x]", "Times[2, x, Cos[Power[x, 2]]]"),
        ("D[Cos[x], x]", "Times[-1, Sin[x]]"),
        ("D[Tan[x], x]", "Power[Cos[x], -2]"),
        ("D[Exp[Times[2, x]], x]", "Times[2, Exp[Times[2, x]]]"),
        ("D[Log[x], x]", "Power[x, -1]"),
        ("D[Sin[Cos[x]], x]", "Times[-1, Cos[Cos[x]], Sin[x]]"),
        // Functions the table does not know.
        ("D[f[x], x]", "Derivative[1][f][x]"),
        (
            "D[f[Power[x, 2]], x]",
            "Times[2, x, Derivative[1][f][Power[x, 2]]]",
        ),
        ("D[Derivative[1][f][x], x]", "Derivative[2][f][x]"),
        ("D[Abs[x], x]", "Derivative[1][Abs][x]"),
        // Higher and partial derivatives, lists.
        ("D[Power[x, 3], List[x, 0]]", "Power[x, 3]"),
        ("D[Power[x, 3], List[x, 2]]", "Times[6, x]"),
        ("D[Power[x, 3], List[x, 3]]", "6"),
        (
            "D[Times[Power[x, 2], Power[y, 3]], x, y]",
            "Times[6, x, Power[y, 2]]",
        ),
        ("D[List[x, Power[x, 2]], x]", "List[1, Times[2, x]]"),
        // A variable that is not a symbol, or an order that is not one.
        ("D[Power[x, 2], 2]", "D[Power[x, 2], 2]"),
        ("D[Power[x, 3], List[x, -1]]", "D[Power[x, 3], List[x, -1]]"),
    ] {
        match ctx.eval(&read(input)) {
            Ok(found) if found == read(expected) => {}
            found => wrong.push(format!("{input}: {found:?}, expected {expected}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The rules have names, and the trace shows them with where they are
/// written (CAS-DOC-4 §9).
#[test]
fn a_trace_names_the_rules() {
    let ctx = Context::new();
    let (_, report) = ctx.eval_with(
        &read("D[Sin[Power[x, 2]], x]"),
        &EvalOptions {
            trace: TraceLevel::Rules,
        },
    );
    let names: Vec<&str> = report
        .trace
        .steps
        .iter()
        .filter_map(|s| match s.kind {
            StepKind::Rule(id) => match id.origin {
                Origin::System { name, location } => {
                    assert!(location.file().ends_with("calculus.rs"), "{location}");
                    Some(name)
                }
                Origin::User => None,
            },
            _ => None,
        })
        .collect();
    assert_eq!(names, ["D/Sin", "D/power", "D/variable"]);
}

/// A rational `p/q` with small `p` and `q`.
fn small() -> impl Strategy<Value = Rational> {
    (-6i64..=6, 1i64..=4).prop_map(|(p, q)| Rational::new(p.into(), q.into()))
}

/// Rational functions of `x`: sums, products and integer powers of `x` and
/// small rationals, unevaluated.
fn rational_function() -> impl Strategy<Value = Expr> {
    let leaf = prop_oneof![
        2 => Just(x()),
        1 => small().prop_map(Expr::from),
    ];
    leaf.prop_recursive(4, 24, 3, |inner| {
        prop_oneof![
            prop::collection::vec(inner.clone(), 2..4).prop_map(|a| apply(Symbol::PLUS, &a)),
            prop::collection::vec(inner.clone(), 2..4).prop_map(|a| apply(Symbol::TIMES, &a)),
            (inner, -3i64..=3).prop_map(|(b, k)| apply(Symbol::POWER, &[b, k.into()])),
        ]
    })
}

/// `e` and its derivative at `x = r`, by forward differentiation over exact
/// rationals: what the rules are checked against, sharing nothing with them.
/// `None` at a pole.
#[allow(
    clippy::many_single_char_names,
    reason = "values and derivatives, as written in mathematics"
)]
fn forward(e: &Expr, r: &Rational) -> Option<(Rational, Rational)> {
    if *e == x() {
        return Some((r.clone(), Rational::from(1)));
    }
    if let Some(c) = e.to_rational() {
        return Some((c, Rational::from(0)));
    }
    let args = e.args();
    match e.head().and_then(Expr::as_symbol)? {
        Symbol::PLUS => {
            args.iter()
                .try_fold((Rational::from(0), Rational::from(0)), |(v, d), a| {
                    let (va, da) = forward(a, r)?;
                    Some((&v + &va, &d + &da))
                })
        }
        Symbol::TIMES => {
            args.iter()
                .try_fold((Rational::from(1), Rational::from(0)), |(v, d), a| {
                    let (va, da) = forward(a, r)?;
                    Some((&v * &va, &(&d * &va) + &(&v * &da)))
                })
        }
        Symbol::POWER => {
            let [b, k] = args else { return None };
            let k = k.as_i64()?;
            let (v, dv) = forward(b, r)?;
            if v.is_zero() && k < 1 {
                return None;
            }
            let derivative = &(&Rational::from(k) * &v.pow(k - 1)) * &dv;
            Some((v.pow(k), derivative))
        }
        _ => None,
    }
}

/// `e` with `x` at `r`, evaluated.
fn at(ctx: &Context, e: &Expr, r: &Rational) -> Expr {
    ctx.eval(&e.substitute(&[(x(), Expr::from(r.clone()))]))
        .unwrap_or_else(|err| panic!("{e:?} at {r:?}: {err}"))
}

/// Whether `e` has no infinity or `Indeterminate` in it.
fn finite(e: &Expr) -> bool {
    [
        Symbol::COMPLEX_INFINITY,
        Symbol::INDETERMINATE,
        Symbol::DIRECTED_INFINITY,
        Symbol::INFINITY,
    ]
    .into_iter()
    .all(|s| !e.contains(&s.into()))
}

/// A function of the table, at an evaluated value, and its derivative there.
fn table(ctx: &Context, g: Symbol, v: &Expr) -> (Expr, Expr) {
    let read_at = |text: &str| {
        ctx.eval(&read(text).substitute(&[(read("v"), v.clone())]))
            .unwrap()
    };
    let derivative = match g {
        Symbol::SIN => "Cos[v]",
        Symbol::COS => "Times[-1, Sin[v]]",
        Symbol::TAN => "Power[Cos[v], -2]",
        Symbol::EXP => "Exp[v]",
        _ => "Power[v, -1]",
    };
    (
        ctx.eval(&apply(g, std::slice::from_ref(v))).unwrap(),
        read_at(derivative),
    )
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(400))]

    /// `D` of a polynomial `Σ c_i x^i` is `Σ i c_i x^(i - 1)`.
    #[test]
    fn polynomials_have_the_closed_form(coefficients in prop::collection::vec(small(), 0..7)) {
        let ctx = Context::new();
        let term = |c: &Rational, i: usize| {
            apply(Symbol::TIMES, &[c.clone().into(), apply(Symbol::POWER, &[x(), Expr::from(i64::try_from(i).unwrap())])])
        };
        let terms: Vec<Expr> = coefficients.iter().enumerate().map(|(i, c)| term(c, i)).collect();
        let closed: Vec<Expr> = coefficients
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, c)| term(&(c * &Rational::from(i64::try_from(i).unwrap())), i - 1))
            .collect();
        let found = d(&ctx, &apply(Symbol::PLUS, &terms));
        prop_assert_eq!(found, ctx.eval(&apply(Symbol::PLUS, &closed)).unwrap());
    }

    /// The product, quotient, power and chain rules on rational functions:
    /// the derivative at a rational point is what forward differentiation
    /// gives there.
    #[test]
    fn rational_functions_agree_with_forward_differentiation(
        f in rational_function(),
        r in small(),
    ) {
        let ctx = Context::new();
        let derivative = d(&ctx, &f);
        if let Some((_, expected)) = forward(&f, &r) {
            let found = at(&ctx, &derivative, &r);
            // The evaluated derivative may have a pole there that `f` as
            // written does not, as `(x^2)^(-1) x` written out.
            if finite(&found) {
                prop_assert_eq!(&found, &Expr::from(expected), "{:?} at {:?}: {:?}", f, r, derivative);
            }
        }
    }

    /// The chain rule through the table: `D` of `g_k[…g_1[p]…]` at a rational
    /// point is the product of the derivatives along the chain there, a
    /// rational times functions at rationals, whose canonical form is one.
    #[test]
    fn the_chain_rule_agrees_with_forward_differentiation(
        p in rational_function(),
        chain in prop::collection::vec(
            prop::sample::select(vec![Symbol::SIN, Symbol::COS, Symbol::TAN, Symbol::EXP, Symbol::LOG]),
            1..4,
        ),
        r in small(),
    ) {
        let ctx = Context::new();
        let Some((v, dv)) = forward(&p, &r) else {
            return Ok(());
        };
        let f = chain.iter().fold(p, |inner, &g| apply(g, &[inner]));
        let (mut value, mut derivative) = (Expr::from(v), Expr::from(dv));
        for &g in &chain {
            let (gv, dg) = table(&ctx, g, &value);
            derivative = ctx.eval(&apply(Symbol::TIMES, &[dg, derivative])).unwrap();
            value = gv;
        }
        let found = at(&ctx, &d(&ctx, &f), &r);
        // A pole along the chain, as `Log` at 0, tells nothing.
        if finite(&found) && finite(&derivative) {
            prop_assert_eq!(found, derivative, "{:?} at {:?}", f, r);
        }
    }
}
