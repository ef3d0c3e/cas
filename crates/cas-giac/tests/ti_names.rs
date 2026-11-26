//! Xcas's `type` and its `DOM_` answers (CAS-621), `iPart`, `fPart`, `frac`,
//! `trunc`, `SortA`, `SortD`, `sortperm` and `RandSeed` (CAS-622),
//! `make_symbol`, `symbol_array`, `LSQ`, `LQ`, `SVD` and `SVL` (CAS-623)
//! against giac: giac's answers read back are this CAS's value of the same input,
//! exact answers against giac's decimals
//! within their error (C0006, C0032); departures, giac's answer as caseval
//! prints it, whitespace collapsed (C0047), and this CAS's own; and random
//! rationals and integer lists a run to `type`, the parts, `trunc` and the
//! sorts.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};
use proptest::prelude::*;

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// Each input: giac's answer, read back, is this CAS's value.
fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| giac.agrees(input).err().map(|e| e.to_string()))
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Each input whose giac answer is a decimal: giac's answer, read back, is
/// this CAS's value within giac's error ([`near`]; C0006, C0032), ours
/// exact where the input is (R67).
fn check_near(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(mine)) if near(&ctx, &theirs, &mine) => {}
            (theirs, mine) => wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `type`, compared through `==` since giac prints its answer as the
/// type's name (probe1, probe2, probe6).
#[test]
fn types() {
    check(&[
        "type(2)==DOM_INT",
        "type(2^100)==DOM_INT",
        "type(1/2)==DOM_RAT",
        "type(2.0)==DOM_FLOAT",
        "type(1+i)==DOM_COMPLEX",
        "type(1.0+i)==DOM_COMPLEX",
        "type(x)==DOM_IDENT",
        "type(pi)==DOM_IDENT",
        "type(infinity)==DOM_IDENT",
        "type(undef)==DOM_IDENT",
        "type(\"a\")==DOM_STRING",
        "type([1,2])==DOM_LIST",
        "type(sin)==DOM_FUNC",
        "type(x->x^2)==DOM_FUNC",
        "type(sqrt(2))==DOM_SYMBOLIC",
        "type(x=1)==DOM_SYMBOLIC",
    ]);
}

/// `iPart`, `fPart` and `frac`; a sequence bracketed (C0023; probe2,
/// probe6); giac's decimal answers (every `iPart`, T8.c4 probe) within its
/// error.
#[test]
fn parts() {
    check(&[
        "fPart(1/2)",
        "fPart(-1/2)",
        "fPart(7/2)",
        "fPart(sqrt(2))",
        "frac(7/3)",
    ]);
    check_near(&[
        "iPart(4.3)",
        "iPart(sqrt(2))",
        "iPart(-4.3)",
        "iPart(7/2)",
        "iPart(-7/2)",
        "[iPart(4.3,sqrt(2))]",
        "iPart([1.5,-2.5])",
        "fPart(1.2)",
        "fPart(-1.2)",
        "fPart(-2.5)",
        "fPart([3.4,sqrt(2)])",
    ]);
}

/// `trunc`, to an integer or to `n` decimals (probe2, probe6); giac's
/// decimal answers (every `trunc(x, n)`, T8.c4 probe) within its error.
#[test]
fn truncation() {
    check(&[
        "trunc(4.3)",
        "trunc(-4.7)",
        "trunc(7/2)",
        "trunc(-7/2)",
        "trunc([4.3333,sqrt(2)])",
    ]);
    check_near(&[
        "trunc(sqrt(2),3)",
        "trunc(pi,4)",
        "trunc(1/3,2)",
        "trunc(2/3,1)",
        "trunc(12345.678,-2)",
        "trunc([4.3333,sqrt(2)],2)",
        "trunc(sqrt(2)+i*sqrt(5),4)",
    ]);
}

/// `SortA`, `sorta`, `SortD`, `sortd` and `sortperm` (CAS-622): lists,
/// sequences, matrices and sequences of rows by their first row, ties
/// (`SortD` reverses them, as giac), symbols and strings; `RandSeed`, `srand`
/// (c2probe1, c2probe2).
#[test]
fn sorting() {
    check(&[
        "SortA([3,4,2])",
        "SortA(3,4,2)",
        "SortA([[3,4,2],[6,4,5]])",
        "SortA([3,4,2],[6,4,5])",
        "SortD([3,4,2])",
        "SortD([[3,4,2],[6,4,5]])",
        "sorta([3,1,2],[6,5,4])",
        "sortd(3,4,2)",
        "SortA([x,1,y])",
        "SortD([x,1,y])",
        "SortD([\"a\",\"c\",\"b\"])",
        "SortD([1/2,0.3,1])",
        "SortA([2,1,2,1])",
        "SortA([[2,1,2],[3,2,1]])",
        "SortD([[2,1,2],[3,2,1]])",
        "SortA([[x,3,1],[1,2,3]])",
        "SortA([[3,1],[6,5],[9,8]])",
        "SortA([5])",
        "sortperm([30,25,40,10,20])",
        "sortperm([3,1,2,1])",
        "sortperm([\"b\",\"a\",\"c\"])",
        "sortperm([5])",
        "sortperm([\"ab\",\"cd\",\"ad\",\"bc\",\"de\"],[3,4,1,0,2])",
        "sortperm([3,1,2],[2,0,1])",
        "RandSeed(321)",
        "RandSeed(5)==srand(5)",
    ]);
}

/// Departures, giac's answer as caseval prints it (probe1, probe2, probe6)
/// and ours: giac's types are integers, printed as the type's name, and
/// `true` is one; giac writes `fPart` and `trunc` of a symbol with `floor`,
/// where ours stays; giac's `trunc(x, n)` floors; `make_symbol` of a
/// non-identifier and `symbol_array` of a template that would not give
/// identifiers stay; `QR` and `LU` storing into variables are giac's
/// (CAS-623, no code); giac's `SVL([[1,4],[4,1]])` is not sorted, ours
/// ascends.
#[test]
fn departures() {
    let cases = [
        ("DOM_INT+0", "2", "DOM_INT"),
        ("type(2)", "integer", "DOM_INT"),
        ("type(true)", "integer", "DOM_IDENT"),
        ("fPart(x)", "x-sign(x)*floor(x*sign(x))", "fPart(x)"),
        ("trunc(x)", "floor(x)", "trunc(x)"),
        ("trunc(-sqrt(2),3)", "-1.415", "-707/500"),
        ("make_symbol(\"x y\")", "undef", "make_symbol(\"x y\")"),
        ("make_symbol(\"1a\")", "a", "make_symbol(\"1a\")"),
        (
            "symbol_array(\"%a\",2)",
            "[0a,1a]",
            "symbol_array(\"%a\",2)",
        ),
        (
            "QR([[1,2],[3,4]],Q,R)",
            "[[sqrt(10),7/5*sqrt(10)],[0,sqrt(10)/5]]",
            "QR([[1,2],[3,4]],Q,R)",
        ),
        (
            "LU([[1,2],[3,4]],L,U,P)",
            "[[1,0],[0,1]]",
            "LU([[1,2],[3,4]],L,U,P)",
        ),
        ("SVL([[1,4],[4,1]])", "[5.0,3.0]", "[3,5]"),
    ];
    let collapsed = |s: &str| s.split_whitespace().collect::<Vec<_>>().join(" ");
    let mut wrong = Vec::new();
    for (input, answer, mine) in cases {
        if let Some(giac) = giac() {
            let printed = giac.eval(input).unwrap_or_default();
            if collapsed(&printed) != collapsed(answer) {
                wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
            }
        }
        let found = ours(input);
        if found.as_ref().ok() != ours(mine).ok().as_ref() {
            wrong.push(format!("{input}: ours {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `make_symbol`, `symbol_array` and `LSQ`, `lsq` (CAS-623): giac's
/// answers read back, exact (probe2, probe6, T8c3 probe).
#[test]
fn symbols_and_least_squares() {
    check(&[
        "make_symbol(\"abc\")",
        "symbol_array(\"x\",5)",
        "symbol_array(\"a\",2,3)",
        "symbol_array(\"a%%\",2,3)",
        "LSQ([[1,2],[3,4]],[5,11])",
        "LSQ([[1,2],[3,4],[5,6]],[1,2,4])",
        "lsq([[1,2],[3,4],[5,6]],[[1],[2],[4]])",
        "LSQ([[1,2],[3,4]],[[5,-1],[11,-1]])",
    ]);
}

/// The entries of nested lists, in order.
fn leaves(e: &Expr) -> Vec<Expr> {
    if e.head() == Some(&Expr::from(Symbol::LIST)) {
        e.args().iter().flat_map(leaves).collect()
    } else {
        vec![e.clone()]
    }
}

/// Whether `a` and `b` agree entry by entry within a relative 10⁻⁹, giac's
/// 12 digits through a product (C0006).
fn near(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let numeric = |e: &Expr| {
        ctx.eval(&Expr::call(Symbol::N, [e.clone(), Expr::from(30)]))
            .unwrap()
    };
    let (t, o) = (leaves(&numeric(a)), leaves(&numeric(b)));
    t.len() == o.len()
        && t.iter().zip(&o).all(|(a, b)| {
            let gap = Expr::call(
                Symbol::ABS,
                [Expr::call(
                    Symbol::PLUS,
                    [b.clone(), Expr::from(-1) * a.clone()],
                )],
            );
            let scale = Expr::call(
                Symbol::MAX,
                [Expr::from(1), Expr::call(Symbol::ABS, [a.clone()])],
            );
            let bound = scale * Expr::from(10).pow(-9);
            ctx.eval(&Expr::call(Symbol::LESS, [gap, bound]))
                .ok()
                .and_then(|e| e.as_symbol())
                == Some(Symbol::TRUE)
        })
}

/// The matrix product of `a` and `b` (`Times` of two matrices).
fn times(ctx: &Context, a: &Expr, b: &Expr) -> Expr {
    ctx.eval(&Expr::apply(&Symbol::TIMES.into(), &[a.clone(), b.clone()]))
        .unwrap()
}

/// `SVL` (`svl`'s alias) near giac's; `SVD` and `LQ` by their products,
/// `U diag(s) Vᴴ` and `L Q`, giac's against ours, `SVD`'s `s` ascending as
/// giac's, `LQ`'s `P` giac's identity (probe2, probe6).
#[test]
fn decompositions() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "SVL([[1,2],[3,4]])",
        "SVL([[3,0],[4,5]])",
        "SVL([[1,2],[3,4],[5,6]])",
    ] {
        let (Ok(theirs), Ok(mine)) = (giac.value(input), ours(input)) else {
            wrong.push(format!("{input}: not read"));
            continue;
        };
        if !near(&ctx, &theirs, &mine) {
            wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}"));
        }
    }
    let svd_product = |d: &Expr| -> Option<Expr> {
        let [u, s, v] = d.args() else {
            return None;
        };
        let sigma = Expr::call(Symbol::DIAGONAL_MATRIX, [s.clone()]);
        let vh = Expr::call(
            Symbol::CONJUGATE,
            [Expr::call(Symbol::TRANSPOSE, [v.clone()])],
        );
        Some(times(&ctx, &times(&ctx, u, &sigma), &vh))
    };
    for (a, ascending) in [
        ("[[1,2],[3,4]]", true),
        ("[[3,0],[4,5]]", true),
        ("[[1,4],[4,1]]", false),
    ] {
        let input = format!("SVD({a})");
        let (Ok(theirs), Ok(mine)) = (giac.value(&input), ours(&input)) else {
            wrong.push(format!("{input}: not read"));
            continue;
        };
        let products = svd_product(&theirs).zip(svd_product(&mine));
        let s_agree = !ascending
            || theirs
                .args()
                .get(1)
                .zip(mine.args().get(1))
                .is_some_and(|(t, m)| near(&ctx, t, m));
        if !products.is_some_and(|(t, m)| near(&ctx, &t, &m)) || !s_agree {
            wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}"));
        }
    }
    for a in [
        "[[1,2],[3,4]]",
        "[[6,12,18],[5,14,31],[3,8,18]]",
        "[[4,0,0],[8,-4,3]]",
    ] {
        let input = format!("LQ({a})");
        let (Ok(theirs), Ok(mine)) = (giac.value(&input), ours(&input)) else {
            wrong.push(format!("{input}: not read"));
            continue;
        };
        let ([tl, tq, tp], [ml, mq, mp]) = (theirs.args(), mine.args()) else {
            wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}"));
            continue;
        };
        if !near(&ctx, &times(&ctx, tl, tq), &times(&ctx, ml, mq)) || tp != mp {
            wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// A random rational `p/q` in the language, `q` from 1.
fn rational() -> impl Strategy<Value = String> {
    (-60i64..=60, 1i64..=9).prop_map(|(p, q)| format!("{p}/{q}"))
}

/// giac's answer to `input`, read back, is ours (`exact`) or ours within
/// giac's error ([`near`]); else a failure naming both.
fn sweep_check(giac: &Giac, input: &str, exact: bool) -> Result<(), TestCaseError> {
    if exact {
        return giac
            .agrees(input)
            .map_err(|e| TestCaseError::fail(e.to_string()));
    }
    match (giac.value(input), ours(input)) {
        (Ok(theirs), Ok(mine)) if near(&Context::new(), &theirs, &mine) => Ok(()),
        (theirs, mine) => Err(TestCaseError::fail(format!(
            "{input}: giac {theirs:?}, ours {mine:?}"
        ))),
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// `type` of a random rational against each of `DOM_INT`, `DOM_RAT`
    /// and `DOM_FLOAT`, through `==` as in [`types`].
    #[test]
    fn random_type(
        r in rational(),
        dom in prop::sample::select(&["DOM_INT", "DOM_RAT", "DOM_FLOAT"][..]),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        sweep_check(giac, &format!("type({r})=={dom}"), true)?;
    }

    /// `fPart`, `frac` and `trunc` of a random rational, exact; `iPart`,
    /// giac's decimal, near; `trunc(r, n)` near, for `r ≥ 0` and `n` in
    /// `0..=3` only: giac floors a negative (departure row
    /// trunc(-sqrt(2),3)), and its negative `n`, significant digits
    /// through doubles, gives trunc(150,-2) = 140.0 and errs on
    /// trunc(0,-1) (caseval, T8.c2).
    #[test]
    fn random_parts(
        r in rational(),
        name in prop::sample::select(&["fPart", "frac", "trunc", "iPart"][..]),
        n in prop::option::of(0i64..=3),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        match n {
            Some(n) if !r.starts_with('-') => {
                sweep_check(giac, &format!("trunc({r},{n})"), false)?;
            }
            _ => sweep_check(giac, &format!("{name}({r})"), name != "iPart")?,
        }
    }

    /// `SortA`, `SortD` and `sortperm` of a random list of integers.
    #[test]
    fn random_sorting(
        values in prop::collection::vec(-9i64..=9, 1..=7),
        name in prop::sample::select(&["SortA", "SortD", "sortperm"][..]),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let items: Vec<String> = values.iter().map(i64::to_string).collect();
        sweep_check(giac, &format!("{name}([{}])", items.join(",")), true)?;
    }
}
