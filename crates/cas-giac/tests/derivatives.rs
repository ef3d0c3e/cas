//! Numerical derivatives against giac (CAS-595): `numdiff` exactly, after
//! `Expand` and `Together` of the difference, lists element by element;
//! giac's aide example of decimals against the closed form and giac's
//! value; `nDeriv` by value at three points within `10⁻¹⁰` (C0006), since
//! giac answers doubles where this CAS answers exactly.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};
use proptest::prelude::*;

/// The shared library, or `None` with a note: the comparisons skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: comparisons skipped");
    }
    giac
}

/// Whether `a − b` is 0 once expanded and put together, lists element by
/// element.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) || b.head() == Some(&list) {
        return a.head() == b.head()
            && a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(x, y)| same(ctx, x, y));
    }
    let difference = Expr::call(
        Symbol::TOGETHER,
        [Expr::call(Symbol::EXPAND, [a.clone() - b.clone()])],
    );
    ctx.eval(&difference) == Ok(Expr::from(0))
}

/// Whether `|a − b| < tolerance · max(1, |b|)`, numerically.
fn near(ctx: &Context, a: &Expr, b: &Expr, tolerance: Expr) -> bool {
    let gap = Expr::call(
        Symbol::ABS,
        [Expr::call(Symbol::N, [a.clone() - b.clone()])],
    );
    let scale = Expr::call(
        Symbol::MAX,
        [
            Expr::from(1),
            Expr::call(Symbol::ABS, [Expr::call(Symbol::N, [b.clone()])]),
        ],
    );
    let claim = Expr::call(Symbol::LESS, [gap, tolerance * scale]);
    ctx.eval(&claim).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// `numdiff` of exact and symbolic points: orders, several at once, an
/// order past the degree.
#[test]
fn numdiff() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "numdiff([1,3,4],[3,0,1],2)",
        "numdiff([0,1,2,3],[0,1,8,27],3/2,2)",
        "numdiff([0,1,2],[0,1,4],1)",
        "numdiff([0,1,2,3],[0,1,8,27],1,3)",
        "numdiff([1,2],[3,4],1,0)",
        "numdiff([0,1/2,1],[1,2,5],1/4)",
        "numdiff([0,1,2],[0,1,4],1,1,2,5)",
        "numdiff([0,1,2],[a,b,c],t)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's aide example: `f = sin(x) exp(−x)` sampled at decimals, `f″` at
/// `1/π` from the points. Ours is within `4·10⁻⁷` of the closed form
/// `−2 cos(1/π) exp(−1/π)` (giac's own error is `3.91·10⁻⁷`), and within
/// `10⁻⁹` of giac's value.
#[test]
fn numdiff_of_decimals() {
    let Some(giac) = giac() else {
        return;
    };
    let input = "numdiff([0,0.1,0.2,0.4,0.5,0.7,0.8,1],\
        [sin(0)*exp(-0),sin(0.1)*exp(-0.1),sin(0.2)*exp(-0.2),sin(0.4)*exp(-0.4),\
        sin(0.5)*exp(-0.5),sin(0.7)*exp(-0.7),sin(0.8)*exp(-0.8),sin(1)*exp(-1)],1/pi,2)";
    let ctx = Context::new();
    let closed = ours("-2*cos(1/pi)*exp(-1/pi)").unwrap();
    let ours = ours(input).unwrap();
    let theirs = giac.value(input).unwrap();
    assert!(
        near(&ctx, &ours, &closed, Expr::from(4) * Expr::from(10).pow(-7)),
        "ours {ours:?}, f'' {closed:?}"
    );
    assert!(
        near(&ctx, &ours, &theirs, Expr::from(10).pow(-9)),
        "ours {ours:?}, giac's {theirs:?}"
    );
}

/// `nDeriv`: giac's doubles against ours, exact where `h` is, by value at
/// three points of `x`.
#[test]
fn n_deriv() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let x = Expr::symbol("x");
    let mut wrong = Vec::new();
    for input in [
        "nDeriv(x^2,x)",
        "nDeriv(x^2,x,0.1)",
        "nDeriv(x^3,x,0.5)",
        "nDeriv(sin(x),x,1/10)",
    ] {
        let (Ok(theirs), Ok(ours)) = (giac.value(input), ours(input)) else {
            wrong.push(format!("{input}: no answer"));
            continue;
        };
        for point in [Expr::from(3) / 10, Expr::from(17) / 10, Expr::from(29) / 10] {
            let at = |e: &Expr| e.substitute(&[(x.clone(), point.clone())]);
            if !near(&ctx, &at(&ours), &at(&theirs), Expr::from(10).pow(-10)) {
                wrong.push(format!(
                    "{input} at {point:?}: ours {ours:?}, giac's {theirs:?}"
                ));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `nDeriv` of an unknown function and step: giac's
/// `(f(x+h)-(f(x-h)))*0.5/h`, its `0.5` read as `1/2`, the same as ours.
#[test]
fn n_deriv_symbolic() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let input = "nDeriv(f(x),x,h)";
    let half = ours("0.5").unwrap();
    let theirs = giac
        .value(input)
        .unwrap()
        .substitute(&[(half, Expr::from(1) / Expr::from(2))]);
    let ours = ours(input).unwrap();
    assert!(
        same(&ctx, &ours, &theirs),
        "ours {ours:?}, giac's {theirs:?}"
    );
}

/// Where giac errs or answers otherwise and this CAS leaves the call:
/// giac's answer as caseval prints it, whitespace collapsed (C0007,
/// C0047). A repeated abscissa is giac's `undef`; a variable not a symbol
/// gives giac's `0.0` (it substitutes into nothing).
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        (
            "numdiff([1,2],[3,4],1,1/2)",
            "GIAC_ERROR: Expected a nonnegative integer. Error: Bad Argument Type",
            Symbol::NUMERICAL_DERIVATIVE,
        ),
        (
            "numdiff([1,2],[3,4],1,-1)",
            "GIAC_ERROR: Expected a nonnegative integer. Error: Bad Argument Type",
            Symbol::NUMERICAL_DERIVATIVE,
        ),
        (
            "numdiff([1,2,3],[3,4],1)",
            "GIAC_ERROR: Lists X and Y must have equal sizes. Error: Bad Argument Value",
            Symbol::NUMERICAL_DERIVATIVE,
        ),
        (
            "numdiff([1,2],[3,4])",
            "GIAC_ERROR: Too few input arguments. Error: Bad Argument Value",
            Symbol::NUMERICAL_DERIVATIVE,
        ),
        (
            "numdiff([1,1],[3,4],1)",
            "undef",
            Symbol::NUMERICAL_DERIVATIVE,
        ),
        ("nDeriv(x^2,2,0.1)", "0.0", Symbol::N_DERIV),
    ];
    let collapsed = |s: &str| s.split_whitespace().collect::<Vec<_>>().join(" ");
    let mut wrong = Vec::new();
    for (input, answer, head) in cases {
        let printed = giac.eval(input).unwrap();
        if collapsed(&printed) != collapsed(answer) {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        let stays = ours(input)
            .ok()
            .and_then(|e| e.head().and_then(Expr::as_symbol));
        if stays != Some(head) {
            wrong.push(format!("{input}: ours {:?} does not stay", ours(input)));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `p/q` in the language, an integer when `q` divides `p`.
fn ratio(p: i64, q: i64) -> String {
    if p % q == 0 {
        (p / q).to_string()
    } else {
        format!("{p}/{q}")
    }
}

/// A small rational.
fn rational() -> impl Strategy<Value = String> {
    (-9i64..=9, 1i64..=5).prop_map(|(p, q)| ratio(p, q))
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// `numdiff` of 1 to 6 points with distinct abscissae at a rational
    /// point, of an order from 0 to the number of points or none: exactly,
    /// as [`same`].
    #[test]
    fn random_numdiff(
        (xs, ys) in prop::collection::btree_set(-20i64..=20, 1..=6).prop_flat_map(|xs| {
            let n = xs.len();
            (Just(xs.into_iter().collect::<Vec<_>>()), prop::collection::vec(rational(), n))
        }),
        d in 1i64..=4,
        at in rational(),
        order in any::<prop::sample::Index>(),
        given in any::<bool>(),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let abscissae: Vec<String> = xs.iter().map(|&x| ratio(x, d)).collect();
        let order = if given {
            format!(",{}", order.index(xs.len() + 1))
        } else {
            String::new()
        };
        let input = format!(
            "numdiff([{}],[{}],{at}{order})",
            abscissae.join(","),
            ys.join(","),
        );
        let ctx = Context::new();
        match (giac.value(&input), ours(&input)) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, ours) => {
                return Err(TestCaseError::fail(format!(
                    "{input}: giac {theirs:?}, ours {ours:?}"
                )));
            }
        }
    }

    /// `nDeriv` of a polynomial of degree ≤ 4, perhaps over `x^2+k` (no
    /// real pole), with a rational step: giac's doubles against ours by
    /// value at three points within `10⁻¹⁰`, as in [`n_deriv`].
    #[test]
    fn random_n_deriv(
        c in prop::collection::vec(-5i64..=5, 1..=5),
        k in 0i64..=4,
        h in (1i64..=3, 1i64..=9),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let terms: Vec<String> = (0..).zip(&c).map(|(e, a)| format!("({a})*x^{e}")).collect();
        let mut f = terms.join("+");
        if k > 0 {
            f = format!("({f})/(x^2+{k})");
        }
        let input = format!("nDeriv({f},x,{})", ratio(h.0, h.1));
        let (Ok(theirs), Ok(ours)) = (giac.value(&input), ours(&input)) else {
            return Err(TestCaseError::fail(format!("{input}: no answer")));
        };
        let ctx = Context::new();
        let x = Expr::symbol("x");
        for point in [Expr::from(3) / 10, Expr::from(17) / 10, Expr::from(29) / 10] {
            let at = |e: &Expr| e.substitute(&[(x.clone(), point.clone())]);
            if !near(&ctx, &at(&ours), &at(&theirs), Expr::from(10).pow(-10)) {
                return Err(TestCaseError::fail(format!(
                    "{input} at {point:?}: giac {theirs:?}, ours {ours:?}"
                )));
            }
        }
    }
}
