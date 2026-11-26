//! Sampling and intervals against giac (CAS-593): `linspace`,
//! `linear_interpolate`, `interval2center` and `center2interval`. giac
//! answers doubles where this CAS answers exactly, so lists and ranges are
//! compared element by element, numbers within `10⁻¹⁰` relatively.

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

/// Whether `a` is within `10⁻¹⁰ max(1, |b|)` of `b`, lists and ranges
/// element by element.
fn near(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let compound = |e: &Expr| {
        matches!(
            e.head().and_then(Expr::as_symbol),
            Some(Symbol::LIST | Symbol::RANGE)
        )
    };
    if compound(a) || compound(b) {
        return a.head() == b.head()
            && a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(x, y)| near(ctx, x, y));
    }
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
    let claim = Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-10) * scale]);
    ctx.eval(&claim).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// giac's answers, exact or doubles; `center2interval` bracketed, since
/// it answers a sequence (C0023).
#[test]
fn answers() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "linspace(0,1,11)",
        "linspace(1,2,3)",
        "linspace(2,1,3)",
        "linspace(0.0,1.0,3)",
        "linear_interpolate([[1,2,6,9],[3,4,6,7]],1,9,1)",
        "linear_interpolate([[1,2,6,9],[3,4,6,7]],1,9,3)",
        "linear_interpolate([[0,2,3],[1,5,4]],0,3,1/2)",
        "linear_interpolate([[0,1],[0,1]],0,1,0.25)",
        "interval2center(2..5)",
        "interval2center(3)",
        "interval2center(1.5..2)",
        "interval2center([1..3,5..7])",
        "[center2interval([2,5,9],1)]",
        "[center2interval([2,5,8])]",
        "[center2interval([1,2,4])]",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if near(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Where giac errs and this CAS leaves the call: giac's error as caseval
/// prints it, whitespace collapsed (C0007), and the head that stays.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        (
            "linspace(1,2,1)",
            "GIAC_ERROR: linspace(1,2,1) Error: Invalid dimension",
            Symbol::LIN_SPACE,
        ),
        (
            "linear_interpolate([[0,2,3],[1,5,4]],-1,4,1)",
            "GIAC_ERROR: Values out of range Error: Bad Argument Value",
            Symbol::LINEAR_INTERPOLATE,
        ),
        (
            "center2interval([3])",
            "GIAC_ERROR: center2interval([3]) Error: Bad Argument Value",
            Symbol::CENTER_TO_INTERVAL,
        ),
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

/// giac's answer to `input` near ours by [`near`]; else a failure naming
/// the input and both answers.
fn compare(giac: &Giac, input: &str) -> Result<(), TestCaseError> {
    let ctx = Context::new();
    match (giac.value(input), ours(input)) {
        (Ok(theirs), Ok(ours)) if near(&ctx, &ours, &theirs) => Ok(()),
        (theirs, ours) => Err(TestCaseError::fail(format!(
            "{input}: giac {theirs:?}, ours {ours:?}"
        ))),
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// `linspace(a, b, n)` of small rationals, `n ≥ 2` (giac errs on
    /// `n = 1`: departure row linspace(1,2,1)).
    #[test]
    fn random_linspace(a in rational(), b in rational(), n in 2u32..=7) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        compare(giac, &format!("linspace({a},{b},{n})"))?;
    }

    /// `linear_interpolate` of 2 to 6 points with distinct sorted
    /// abscissae, from one abscissa to a later one. The abscissae are
    /// nonnegative: giac misreads negative ones (probe:
    /// `linear_interpolate([[-3,-2,-1],[2,3,9]],-3,-1,1/2)` errs "Values out
    /// of range", `linear_interpolate([[-1,-1/4,0,1/2,5/3],[1/2,7,1,-1/2,9/5]],0,0,3)`
    /// answers `[[0.0],[4.5]]` for the point `(0,1)`). Abscissae and steps
    /// are dyadic, so giac's doubles reach the last point exactly.
    #[test]
    fn random_linear_interpolate(
        (xs, ys) in prop::collection::btree_set(0i64..=40, 2..=6).prop_flat_map(|xs| {
            let n = xs.len();
            (Just(xs.into_iter().collect::<Vec<_>>()), prop::collection::vec(rational(), n))
        }),
        d in prop::sample::select(&[1i64, 2, 4][..]),
        bounds in (any::<prop::sample::Index>(), any::<prop::sample::Index>()),
        step in (1i64..=4, prop::sample::select(&[1i64, 2, 4][..])),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let (i, j) = (bounds.0.index(xs.len()), bounds.1.index(xs.len()));
        let (i, j) = (i.min(j), i.max(j));
        let abscissae: Vec<String> = xs.iter().map(|&x| ratio(x, d)).collect();
        let input = format!(
            "linear_interpolate([[{}],[{}]],{},{},{})",
            abscissae.join(","),
            ys.join(","),
            abscissae[i],
            abscissae[j],
            ratio(step.0, step.1),
        );
        compare(giac, &input)?;
    }
}
