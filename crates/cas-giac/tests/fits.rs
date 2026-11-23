//! B-splines and spline fits against giac (CAS-596): `bspline` exactly,
//! after `Expand` and `Together` of the difference element by element
//! (C0032), `piecewise` branch by branch; `fitspline` at points, ours exact
//! against the closed form and within `10⁻¹⁰` of giac's doubles (C0006);
//! with a variable, by value at three points of each piece within `10⁻⁹`.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};

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
