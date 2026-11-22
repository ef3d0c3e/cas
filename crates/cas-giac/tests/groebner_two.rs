//! Ideal membership, giac's order of the variables, syzygies and the
//! rational univariate representation against giac (CAS-541): `in_ideal`
//! exactly, `gbasis_vars` exactly, `syzygy`'s generators as sets up to
//! constant factors, `gbasis`'s `rur` element by element. giac's knobs of
//! its modular engine are departures (CAS-542).

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

/// Whether `a` is a non-zero constant multiple of `b`: their quotient, put
/// together, a number.
fn proportional(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let quotient = Expr::call(Symbol::TOGETHER, [a.clone() / b.clone()]);
    ctx.eval(&quotient)
        .is_ok_and(|q| q.is_number() && q != Expr::from(0))
}

/// Whether two answers are one basis: lists of as many elements, each of
/// one a constant multiple of one of the other's.
fn same_basis(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    a.head() == Some(&list)
        && b.head() == Some(&list)
        && a.args().len() == b.args().len()
        && a.args()
            .iter()
            .all(|p| b.args().iter().any(|q| proportional(ctx, p, q)))
}

/// Whether two answers are one expression, their difference expanded 0.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let difference = Expr::call(Symbol::EXPAND, [a.clone() - b.clone()]);
    ctx.eval(&difference).is_ok_and(|d| d == Expr::from(0))
}

/// Whether two answers are the same expression.
fn equal(_: &Context, a: &Expr, b: &Expr) -> bool {
    a == b
}

/// Whether two `rur` answers are one: lists of as many elements, the first
/// two equal (`rur` and the form), the polynomials the same.
fn same_rur(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    a.head() == Some(&list)
        && b.head() == Some(&list)
        && a.args().len() == b.args().len()
        && a.args().len() > 2
        && a.args()[..2] == b.args()[..2]
        && a.args()[2..]
            .iter()
            .zip(&b.args()[2..])
            .all(|(p, q)| same(ctx, p, q))
}

/// giac's answers and this CAS's to `inputs`, alike by `alike`.
fn check(inputs: &[&str], alike: fn(&Context, &Expr, &Expr) -> bool) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if alike(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `in_ideal`: 1 or 0, a list for a list (CAS-541).
#[test]
fn memberships() {
    check(
        &[
            "in_ideal(x^3-1,[x-1],[x])",
            "in_ideal(x^3-1,[x-2],[x])",
            "in_ideal(x*y,[x^2,y],[x,y])",
            "in_ideal([x*y,x+1],[x^2,y],[x,y])",
            "in_ideal((x+y)^2,[y^2,x^2+2*x*y],[x,y])",
            "in_ideal(x+y,[y^2,x^2+2*x*y],[x,y])",
            "in_ideal(x+y,[y^2,x^2+2*x*y],[x,y],plex)",
            "in_ideal(x^3-1,[x-1])",
        ],
        equal,
    );
}

/// `gbasis_vars`: first appearance, stably sorted by weight (CAS-541).
#[test]
fn variables() {
    check(
        &[
            "gbasis_vars([x^2-y,z])",
            "gbasis_vars([x^2+y^2-1,x-y])",
            "gbasis_vars([x*y,x^3-y^2])",
            "gbasis_vars([x*y*z-1,x^2+t])",
            "gbasis_vars([3,x])",
            "gbasis_vars([])",
            "gbasis_vars(x^2+y)",
        ],
        equal,
    );
}

/// `syzygy`'s generators, `A0`, `A1`, … in both (CAS-541); giac's
/// variables there are `[y,x]` for `[x,y]`, the generators the same.
#[test]
fn syzygies() {
    check(
        &["syzygy([x,y])[0]", "syzygy([x*y,x^2-y^3])[0]"],
        same_basis,
    );
}

/// `gbasis`'s `rur` (CAS-541): the form, then `f`, `f'` and the numerators.
#[test]
fn representations() {
    check(
        &[
            "gbasis([x^2+y^2-1,x-y],[x,y],rur)",
            "gbasis([x^2-2,y^2-3],[x,y],rur)",
            "gbasis([x^2+y^2+z^2-1,x-y,y-z],[x,y,z],rur)",
            "gbasis([x^2-2,y-x],[x,y],rur)",
            "gbasis([x*y-1,x^2+y^2-4],[x,y],rur)",
            "gbasis([x^3-y,y^2-x],[x,y],rur)",
            "gbasis([x^2,y],[x,y],rur)",
            "gbasis([x^2*(x-1),y],[x,y],rur)",
            "gbasis([x^2-1,(y-1)^2],[x,y],rur)",
            "gbasis([x-1,y-2],[x,y],rur)",
            "gbasis([x^2-3*x+2,y-x^2],[x,y],rur)",
            "gbasis([x^2+y^2-1,x-y],[x,y],plex,rur)",
        ],
        same_rur,
    );
    // No representation: one variable, not zero-dimensional, the basis.
    check(
        &[
            "gbasis([x^2+y^2-1],[x,y],rur)",
            "gbasis([x^2-y^2,x-y],[x,y],rur)",
            "gbasis([x^3-2],[x],rur)",
        ],
        same_basis,
    );
}

/// Where this CAS departs from giac on purpose: giac's answer as it
/// prints it, and this CAS's as decided.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        // giac reduces by the generators as given; ours by their reduced
        // basis: x² = x(x + y)/2 + x(x − y)/2.
        ("in_ideal(x^2,[x+y,x-y],[x,y])", "0", "1"),
        // Ties: giac's first appearance is in its order of the input, ours
        // in the evaluated one.
        ("gbasis_vars([a*x+b,y^3-a])", "[x,b,a,y]", "[b,x,a,y]"),
        // Symbols only: no kernel `sin(x)`.
        (
            "gbasis_vars([sin(x)+y])",
            "[sin(x),y]",
            "gbasis_vars([sin(x)+y])",
        ),
        // giac 2.0.0 gives [1] whenever the variables are given.
        (
            "syzygy([x,y],[x,y])",
            "[1],[A0,A1,x,y],[A0,A1],[x,y]",
            "[[A1*x-A0*y,A1^2,A0*A1,A0^2],[A0,A1,x,y],[A0,A1],[x,y]]",
        ),
        // giac's generators miss the syzygy (x + 1)A0 − A1.
        (
            "syzygy([x-1,x^2-1])[0]",
            "[A0*x^2-A1*x-A0+A1,A0^2,A0*A1,A1^2]",
            "[A0-A1+A0*x,A1^2,A0*A1,A0^2]",
        ),
        // giac's four syzygies are multiples of the two minimal ones,
        // A1 x − A0 y and A2 x − A1 y, which it misses.
        (
            "syzygy([x^2,x*y,y^2])[0]",
            "[A1*x*y^2-A0*y^3,A1*x^2-A0*x*y,A2*x^2-A0*y^2,A2*x*y-A1*y^2,A0^2,A0*A1,A1^2,A0*A2,A1*A2,A2^2]",
            "[A2*x-A1*y,A1*x-A0*y,A2^2,A1*A2,A0*A2,A1^2,A0*A1,A0^2]",
        ),
        // Another separating form: ours the first of x + k y that is.
        (
            "gbasis([x^2-1,y^2-1],[x,y],rur)",
            "[rur,-6*x-4*y,-x^4+104*x^2-400,-4*x^3+208*x,24*x^2-480,16*x^2+320]",
            "[rur,x+2*y,-x^4+10*x^2-9,-4*x^3+20*x,-4*x^2-12,-8*x^2+24]",
        ),
        // Knobs of giac's modular F4 engine; ours is Buchberger over ℚ
        // (CAS-542): the calls stay.
        ("gbasis_max_pairs(4096)", "32768", "gbasis_max_pairs(4096)"),
        ("gbasis_simult_primes(16)", "20", "gbasis_simult_primes(16)"),
        (
            "gbasis_reinject(0.1,0.05)",
            "[0.2,0.125]",
            "gbasis_reinject(0.1,0.05)",
        ),
        ("rur_gbasis(0)", "-1", "rur_gbasis(0)"),
        ("rur_certify(1)", "0", "rur_certify(1)"),
    ];
    let mut wrong = Vec::new();
    for (input, answer, decided) in cases {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        if ours(input) != ours(decided) {
            wrong.push(format!("{input}: ours is not {decided}"));
        }
        if giac.agrees(input).is_ok() {
            wrong.push(format!("{input} agrees now"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
