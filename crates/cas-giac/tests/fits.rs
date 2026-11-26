//! B-splines and spline fits against giac (CAS-596): `bspline` exactly,
//! after `Expand` and `Together` of the difference element by element
//! (C0032), `piecewise` branch by branch; `fitspline` at points, ours exact
//! against the closed form and within `10⁻¹⁰` of giac's doubles (C0006);
//! with a variable, by value at three points of each piece within `10⁻⁹`.
//! `fitpoly` (CAS-597) coefficient by coefficient within
//! `10⁻⁹ max(1, max |c|)` of giac's (C0006), a coefficient giac leaves out
//! taken as 0; data rows also exactly against the closed form.

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

/// Whether `|a − b| < tolerance · max(1, |b|)`, numerically, lists element
/// by element.
fn near(ctx: &Context, a: &Expr, b: &Expr, tolerance: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) || b.head() == Some(&list) {
        return a.head() == b.head()
            && a.args().len() == b.args().len()
            && a.args()
                .iter()
                .zip(b.args())
                .all(|(x, y)| near(ctx, x, y, tolerance));
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
    let claim = Expr::call(Symbol::LESS, [gap, tolerance.clone() * scale]);
    ctx.eval(&claim).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// `bspline`: degrees, knots `x`, `x=a..b` and `x=[…]`, complex, symbolic
/// and matrix control points, giac's collapse of one column.
#[test]
fn bspline() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "bspline([2*i,1-2*i,5-i,6-4*i,8+2*i,5+i],x=0..15)",
        "bspline([1,2,3,4],x)",
        "bspline([1,2,3,4,5],x=0..1,2)",
        "bspline([0,1,0],x,1)",
        "bspline([1,2,3,4,5,6],x=0..1)",
        "bspline([[1,2],[3,4],[5,6]],x,1)",
        "bspline([1,2,3],x=[0,1/2,1],1)",
        "bspline([a,b,c],x,1)",
        "bspline([[1],[2],[3]],x,1)",
        "bspline([[1],[2]],x,1)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Whether two conditions in `x` agree at points around the knots.
fn same_condition(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let x = Expr::symbol("x");
    [-1, 0, 1, 2, 3, 4, 5, 6, 7]
        .into_iter()
        .map(|k| Expr::from(k) / Expr::from(6))
        .all(|point| {
            let at = |e: &Expr| ctx.eval(&e.substitute(&[(x.clone(), point.clone())]));
            at(a).is_ok() && at(a) == at(b)
        })
}

/// `bspline(…, piecewise)`: one `Piecewise` per column, branch by branch,
/// the conditions by their truth at points (giac's `1>x` is our `x<1`).
#[test]
fn bspline_piecewise() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "bspline([1,2,3,4],x,3,piecewise)",
        "bspline([1,2,3,4,5],x=0..1,2,piecewise)",
    ] {
        let (Ok(theirs), Ok(ours)) = (giac.value(input), ours(input)) else {
            wrong.push(format!("{input}: no answer"));
            continue;
        };
        let agree = ours.args().len() == theirs.args().len()
            && ours.args().iter().zip(theirs.args()).all(|(o, t)| {
                o.head() == Some(&Expr::from(Symbol::PIECEWISE))
                    && o.args().len() == t.args().len()
                    && o.args().chunks(2).zip(t.args().chunks(2)).all(|(o, t)| {
                        same_condition(&ctx, &o[0], &t[0]) && same(&ctx, &o[1], &t[1])
                    })
            });
        if !agree {
            wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `fitspline` at points: ours exact equals the closed form (the exact
/// least-squares fit) and lies within `10⁻¹⁰` of giac's doubles.
#[test]
fn fitspline_points() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let tolerance = Expr::from(10).pow(-10);
    let six = "[[0,0],[1,1],[2,0],[3,1],[4,0],[5,1]]";
    let seven = "tran([[0,1,2,3,4,5,6],[0,1,0,1,0,1,0]])";
    let rows = [
        (format!("fitspline({six},[1/2,5/2])"), "[61/126,1/2]"),
        (
            format!("fitspline({seven},[1/2,3/2])"),
            "[1461/2212,1569/2212]",
        ),
        (format!("fitspline({seven},[3])"), "[169/553]"),
        (
            format!("fitspline({six},[1/2,5/2],breakpoint=[0,2,5])"),
            "[2053/3276,61/156]",
        ),
        (
            format!("fitspline({six},[1/2,5/2],degree=2)"),
            "[135/334,1/2]",
        ),
        (format!("fitspline({six},[1/2,5/2],1)"), "[61/126,1/2]"),
        (format!("fitspline({six},[1/2,5/2],3)"), "[229/170,1/2]"),
        (
            "fitspline(tran([[0,1,2,3,4,5,6],[0,1,0,1,0,1,0],[1,2,3,4,5,6,8]]),[1/2,5/2])"
                .to_owned(),
            "[[1461/2212,106861/70784],[783/2212,247775/70784]]",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, closed) in rows {
        let closed = ours(closed).unwrap();
        let (Ok(theirs), Ok(ours)) = (giac.value(&input), ours(&input)) else {
            wrong.push(format!("{input}: no answer"));
            continue;
        };
        if ours != closed {
            wrong.push(format!("{input}: ours {ours:?}, closed form {closed:?}"));
        }
        if !near(&ctx, &ours, &theirs, &tolerance) {
            wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `fitspline` with a variable: each piece by value at three points of its
/// interval, `n` equal intervals of `[min s, max s]`, within `10⁻⁹`.
#[test]
fn fitspline_variable() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let x = Expr::symbol("x");
    let tolerance = Expr::from(10).pow(-9);
    let mut wrong = Vec::new();
    for (input, top) in [
        ("fitspline(tran([[0,1,2,3,4,5,6],[0,1,0,1,0,1,0]]),x)", 6),
        ("fitspline([[0,0],[1,1],[2,0],[3,1],[4,0],[5,1]],x)", 5),
    ] {
        let (Ok(theirs), Ok(ours)) = (giac.value(input), ours(input)) else {
            wrong.push(format!("{input}: no answer"));
            continue;
        };
        let pieces = ours.args().len();
        if pieces == 0 || theirs.args().len() != pieces {
            wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}"));
            continue;
        }
        let n = i64::try_from(pieces).unwrap();
        for (j, (o, t)) in (0..).zip(ours.args().iter().zip(theirs.args())) {
            for k in 1..=3 {
                let point = Expr::from(top) * Expr::from(4 * j + k) / Expr::from(4 * n);
                let at = |e: &Expr| e.substitute(&[(x.clone(), point.clone())]);
                if !near(&ctx, &at(o), &at(t), &tolerance) {
                    wrong.push(format!("{input} at {point:?}: ours {o:?}, giac's {t:?}"));
                }
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Where giac answers otherwise: giac's answer as caseval prints it,
/// whitespace collapsed (C0007), and ours. Points outside the data: giac's
/// last polynomial continued, ours stays. `derive=1`: giac's rows of one,
/// ours flat as without it. A variable with `breakpoint=`: giac builds on
/// equal intervals, ours on the breakpoints, so ours at a point of each
/// piece is the points form's. Knots too few: giac's error, ours stays.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let collapsed = |s: &str| s.split_whitespace().collect::<Vec<_>>().join(" ");
    let six = "[[0,0],[1,1],[2,0],[3,1],[4,0],[5,1]]";
    let mut wrong = Vec::new();
    let mut printed = |input: &str, answer: &str| {
        let text = giac.eval(input).unwrap();
        if collapsed(&text) != collapsed(answer) {
            wrong.push(format!("{input}: giac prints {text}, recorded {answer}"));
        }
    };
    for point in ["7", "-1"] {
        printed(&format!("fitspline({six},[{point}])"), "[0.936507936508]");
    }
    printed(
        &format!("fitspline({six},[1/2,5/2],derive=1)"),
        "[[0.600529100529],[-0.28835978836]]",
    );
    printed(
        &format!("fitspline({six},x,breakpoint=[0,1,5])"),
        "[0.141257142857*x^3-0.71619047619*x^2+0.4*x*(-0.500952380952*x+1.31428571429)\
         +1.05142857143*x,0.0297142857143*x^3-0.0365714285715*x^2\
         +0.4*x*(-0.56*x+0.638095238095)-0.4*x*(-0.500952380952*x+1.31428571429)\
         +0.4*x*(-0.0495238095238*x+0.247619047619)-1.0*(-0.56*x+0.638095238095)\
         +2.0*(-0.500952380952*x+1.31428571429)-1.0*(-0.0495238095238*x+0.247619047619)\
         +0.0495238095238*x]",
    );
    printed(
        "bspline([0,1,0,1,0],x=[0,1/2,1],2)",
        "GIAC_ERROR: Expected a variable. Error: Bad Argument Value",
    );
    let head = |e: &Result<Expr, _>| {
        e.as_ref()
            .ok()
            .and_then(|e| e.head().and_then(Expr::as_symbol))
    };
    for (input, stays) in [
        (format!("fitspline({six},[7])"), Symbol::SPLINE_FIT),
        (format!("fitspline({six},[-1])"), Symbol::SPLINE_FIT),
        (
            "bspline([0,1,0,1,0],x=[0,1/2,1],2)".to_owned(),
            Symbol::B_SPLINE,
        ),
    ] {
        let found = ours(&input);
        if head(&found) != Some(stays) {
            wrong.push(format!("{input}: ours {found:?} does not stay"));
        }
    }
    let derived = ours(&format!("fitspline({six},[1/2,5/2],derive=1)"));
    let expected = ours("[227/378,-109/378]").unwrap();
    if derived.as_ref() != Ok(&expected) {
        wrong.push(format!("derive=1: ours {derived:?}, expected {expected:?}"));
    }
    let x = Expr::symbol("x");
    let pieces = ours(&format!("fitspline({six},x,breakpoint=[0,1,5])")).unwrap();
    let points = ours(&format!("fitspline({six},[1/2,3],breakpoint=[0,1,5])")).unwrap();
    let at = |e: &Expr, v: Expr| e.substitute(&[(x.clone(), v)]);
    let sampled = match pieces.args() {
        [first, second] => Expr::call(
            Symbol::LIST,
            [
                at(first, Expr::from(1) / Expr::from(2)),
                at(second, Expr::from(3)),
            ],
        ),
        _ => pieces.clone(),
    };
    if !same(&ctx, &sampled, &points) {
        wrong.push(format!("breakpoint=: ours {pieces:?}, points {points:?}"));
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The coefficients of `p` in `x`, the constant first; `None` when it is
/// no polynomial.
fn coefficients(ctx: &Context, p: &Expr) -> Option<Vec<Expr>> {
    let list = ctx
        .eval(&Expr::call(
            Symbol::COEFFICIENT_LIST,
            [p.clone(), Expr::symbol("x")],
        ))
        .ok()?;
    (list.head() == Some(&Expr::from(Symbol::LIST)))
        .then(|| list.args().iter().rev().cloned().collect())
}

/// Whether the polynomials `a` and `b` in `x` agree coefficient by
/// coefficient within `10⁻⁹ max(1, max |c|)`, `c` over `b`'s, a coefficient
/// one leaves out taken as 0.
fn near_polynomial(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let (Some(a), Some(b)) = (coefficients(ctx, a), coefficients(ctx, b)) else {
        return false;
    };
    let size = |c: &Expr| Expr::call(Symbol::ABS, [Expr::call(Symbol::N, [c.clone()])]);
    let scale = Expr::call(
        Symbol::MAX,
        std::iter::once(Expr::from(1)).chain(b.iter().map(size)),
    );
    let tolerance = Expr::from(10).pow(-9) * scale;
    let zero = Expr::from(0);
    (0..a.len().max(b.len())).all(|k| {
        let gap = size(&(a.get(k).unwrap_or(&zero).clone() - b.get(k).unwrap_or(&zero).clone()));
        let claim = Expr::call(Symbol::LESS, [gap, tolerance.clone()]);
        ctx.eval(&claim).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
    })
}

/// `fitpoly` of functions: Chebyshev's sums, the degree by giac's rule,
/// `degree=` and `limit=`; Runge's function, where ours keeps degree 15
/// with odd coefficients holding 0 and giac prints degree 14.
#[test]
fn fitpoly_functions() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "fitpoly(exp(x),x=0..1)",
        "fitpoly(x*cos(x),x=0..pi/2)",
        "fitpoly(sin(x),x=0..pi)",
        "fitpoly(exp(x),x=0..1,degree=2)",
        "fitpoly(exp(x),x=0..1,limit=2)",
        "fitpoly(x^2,x=0..1)",
        "fitpoly(1/(1+25*x^2),x=-1..1)",
    ] {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if near_polynomial(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `fitpoly` of data with a variable: ours equals the closed form (exact
/// data exactly, decimals within the tolerance) and is near giac's.
#[test]
fn fitpoly_data() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let eight = "[0,1,2,3,4,5,6,7],[1,3,2,5,4,7,6,9]";
    let squares = "[0,1,2,3,4,5,6,7,8,9],[0,1,4,9,16,25,36,49,64,81.5]";
    let quadratic = "89/88*x^2-3/40*x+7/110";
    let rows = [
        (
            "fitpoly([0,1,2,3],[1,2,5,10],x)".to_owned(),
            "x^2+1".to_owned(),
            true,
        ),
        (
            "fitpoly([[0,1],[1,2],[2,5],[3,10]],x)".to_owned(),
            "x^2+1".to_owned(),
            true,
        ),
        (
            "fitpoly([0,1,2,3,4],[1,2,5,10,17.5],x)".to_owned(),
            "x^4/48-x^3/8+59*x^2/48-x/8+1".to_owned(),
            false,
        ),
        (
            format!("fitpoly({eight},x)"),
            format!("expand(lagrange({eight},x))"),
            true,
        ),
        (
            format!("fitpoly({eight},x,degree=2)"),
            "x^2/24+121*x/168+11/8".to_owned(),
            true,
        ),
        (format!("fitpoly({squares},x)"), quadratic.to_owned(), false),
        (
            format!("fitpoly({squares},x,len=3)"),
            quadratic.to_owned(),
            false,
        ),
        (format!("fitpoly({squares},x=2..7)"), "x^2".to_owned(), true),
    ];
    let mut wrong = Vec::new();
    for (input, closed, exact) in rows {
        let closed = ours(&closed).unwrap();
        let (Ok(theirs), Ok(ours)) = (giac.value(&input), ours(&input)) else {
            wrong.push(format!("{input}: no answer"));
            continue;
        };
        let equal = if exact {
            ours == closed
        } else {
            near_polynomial(&ctx, &ours, &closed)
        };
        if !equal {
            wrong.push(format!("{input}: ours {ours:?}, closed form {closed:?}"));
        }
        if !near_polynomial(&ctx, &ours, &theirs) {
            wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `fitpoly` where giac answers otherwise, giac's answer as caseval prints
/// it (C0007): without a variable giac's `poly1[…]`, ours the coefficient
/// list; `cos` with `1e-8`, giac's degree 15 (its `gaussquad` of the
/// weighted norm is off), ours degree 5 near the prototype's; `abs` on
/// `[-1, 1]`, giac's coefficients off ours by up to `1.3·10⁻⁶` (its
/// doubles in the degree-15 expansion), ours near the prototype's.
#[test]
fn fitpoly_departures() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let collapsed = |s: &str| s.split_whitespace().collect::<Vec<_>>().join(" ");
    let aide = "fitpoly([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],\
                [0.25,1.12,1.89,2.63,3.19,3.47,3.53,3.18,2.64,1.71,0.09])";
    let mut wrong = Vec::new();
    let mut printed = |input: &str, answer: &str| {
        let text = giac.eval(input).unwrap();
        if collapsed(&text) != collapsed(answer) {
            wrong.push(format!("{input}: giac prints {text}, recorded {answer}"));
        }
    };
    printed("fitpoly([0,1,2,3],[1,2,5,10])", "poly1[1,0,1.0]");
    printed(
        aide,
        "poly1[-8.20707070707,-0.800116550117,8.89355089355,0.244755244755]",
    );
    printed(
        "fitpoly(cos(x),x=0..2,1e-8)",
        "0.108784834278*x^15-1.63417223592*x^14+11.0469352195*x^13-44.3891631896*x^12\
         +117.971975405*x^11-218.316495173*x^10+288.43624426*x^9-274.432928612*x^8\
         +187.305412639*x^7-90.2555210984*x^6+29.7964454465*x^5-6.37333689193*x^4\
         +0.830922155851*x^3-0.556296581254*x^2+0.00150307735257*x+0.999993355198",
    );
    printed(
        "fitpoly(abs(x),x=-1..1)",
        "2.31497324421e-07*x^15+52.3976513753*x^14-8.68114966579e-07*x^13\
         -201.356816381*x^12+1.30217244987e-06*x^11+312.576792355*x^10\
         -9.94715065872e-07*x^9-251.427068953*x^8+4.06928890584e-07*x^7\
         +112.172944357*x^6-8.54550670226e-08*x^5-28.1207867728*x^4\
         +7.9125062058e-09*x^3+4.71660667691*x^2-2.11942130512e-10*x+0.0434154361851",
    );
    let list = |text: &str| ours(text).unwrap();
    let pairs = [
        ("fitpoly([0,1,2,3],[1,2,5,10])", "[1,0,1]"),
        (
            aide,
            "[-8.20707070707,-0.800116550117,8.89355089355,0.244755244755]",
        ),
    ];
    let tolerance = Expr::from(10).pow(-10);
    for (input, expected) in pairs {
        let found = ours(input);
        if !found
            .as_ref()
            .is_ok_and(|f| near(&ctx, f, &list(expected), &tolerance))
        {
            wrong.push(format!("{input}: ours {found:?}, expected {expected}"));
        }
    }
    for (input, prototype) in [
        (
            "fitpoly(cos(x),x=0..2,1e-8)",
            "-0.00672524426254*x^5+0.0550363612025*x^4-0.0127904889208*x^3\
             -0.494334869161*x^2-0.000931225201068*x+1.00002504404",
        ),
        (
            "fitpoly(abs(x),x=-1..1)",
            "52.3976511367*x^14-201.356815562*x^12+312.576791255*x^10\
             -251.427068223*x^8+112.17294411*x^6-28.1207867333*x^4\
             +4.71660667457*x^2+0.0434154362073",
        ),
    ] {
        let found = ours(input);
        let degree = found
            .as_ref()
            .ok()
            .and_then(|f| coefficients(&ctx, f))
            .map(|c| c.len());
        let expected = if input.contains("cos") { 6 } else { 16 };
        if degree != Some(expected)
            || !found
                .as_ref()
                .is_ok_and(|f| near_polynomial(&ctx, f, &list(prototype)))
        {
            wrong.push(format!(
                "{input}: ours {found:?}, the prototype's {prototype}"
            ));
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

    /// `fitspline` of 5 or 6 points with distinct sorted abscissae (giac
    /// wants more samples than control points: 4 points err) at 1 to 3
    /// points inside the data (outside: departure rows): ours within
    /// `10⁻¹⁰` of giac's doubles, as in [`fitspline_points`].
    #[test]
    fn random_fitspline(
        (xs, ys) in prop::collection::btree_set(-20i64..=20, 5..=6).prop_flat_map(|xs| {
            let n = xs.len();
            (Just(xs.into_iter().collect::<Vec<_>>()), prop::collection::vec(rational(), n))
        }),
        d in 1i64..=4,
        eighths in prop::collection::vec(0i64..=8, 1..=3),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let data: Vec<String> = xs
            .iter()
            .zip(&ys)
            .map(|(&x, y)| format!("[{},{y}]", ratio(x, d)))
            .collect();
        let (first, last) = (xs[0], xs[xs.len() - 1]);
        let points: Vec<String> = eighths
            .iter()
            .map(|&k| ratio(8 * first + (last - first) * k, 8 * d))
            .collect();
        let input = format!("fitspline([{}],[{}])", data.join(","), points.join(","));
        let (Ok(theirs), Ok(ours)) = (giac.value(&input), ours(&input)) else {
            return Err(TestCaseError::fail(format!("{input}: no answer")));
        };
        let ctx = Context::new();
        if !near(&ctx, &ours, &theirs, &Expr::from(10).pow(-10)) {
            return Err(TestCaseError::fail(format!(
                "{input}: giac {theirs:?}, ours {ours:?}"
            )));
        }
    }

    /// `bspline` of more rational control points than the degree (giac's
    /// "Invalid dimension" otherwise), knots `x` or `x=a..b`: exactly, as
    /// [`same`].
    #[test]
    fn random_bspline(
        (degree, control) in (1usize..=3).prop_flat_map(|k| {
            (Just(k), prop::collection::vec(rational(), k + 1..=6))
        }),
        range in prop::option::of((-9i64..=9, 1i64..=9, 1i64..=4)),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let knots = match range {
            Some((a, width, d)) => format!("x={}..{}", ratio(a, d), ratio(a + width, d)),
            None => "x".to_owned(),
        };
        let input = format!("bspline([{}],{knots},{degree})", control.join(","));
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
}
