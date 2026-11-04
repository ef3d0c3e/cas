//! Polynomial tools II against giac (CAS-DOC-48 §5, CAS-PLAN-42): giac's
//! answers read back are this CAS's by value (their difference over one
//! denominator is 0); `proot`'s decimals within 10⁻¹¹ as a list; `abcuv` and `pade` by
//! their identities where giac's normalization differs (giac's `u` is not
//! always reduced modulo `b`).

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// The entries of nested lists, in order.
fn leaves(e: &Expr) -> Vec<Expr> {
    if e.head() == Some(&Expr::from(Symbol::LIST)) {
        e.args().iter().flat_map(leaves).collect()
    } else {
        vec![e.clone()]
    }
}

/// Whether two answers are one value, entry by entry: the difference
/// put over one denominator is 0 (giac's `pade` scales numerator and
/// denominator), or within 10⁻¹¹ by balls.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let (a, b) = (leaves(a), leaves(b));
    a.len() == b.len()
        && a.iter().zip(&b).all(|(x, y)| {
            let difference = Expr::call(Symbol::TOGETHER, [x.clone() - y.clone()]);
            if ctx.eval(&difference).is_ok_and(|d| d == Expr::from(0)) {
                return true;
            }
            let numeric = |e: &Expr| {
                ctx.eval(&Expr::call(Symbol::N, [e.clone(), Expr::from(30)]))
                    .unwrap()
            };
            let gap = Expr::call(Symbol::ABS, [numeric(x) - numeric(y)]);
            let small = Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-11)]);
            ctx.eval(&small).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
        })
}

fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(t), Ok(o)) if same(&ctx, &o, &t) => {}
            (t, o) => wrong.push(format!("{input}: ours {o:?}, giac's {t:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Roots and coefficients (CAS-482).
#[test]
fn roots_and_coefficients() {
    check(&[
        "proot([1,0,-2])",
        "proot([1,0,1])",
        "proot(x^3-1)",
        "proot(x^4-5*x^2+4)",
        "pcoeff([1,2,3])",
        "pcoeff([1,-1,2,-2])",
        "peval([1,2,3],2)",
        "horner(x^3+2*x+1,2)",
        "horner([1,-3,0,5],-1/2)",
        "ptayl([1,0,2,1],2)",
        "ptayl(x^3+2*x+1,2)",
        "froot(x^2*(x-3)^3)",
        "froot((x-1)^2/(x+2))",
        "e2r(x^2+2*x)",
        "r2e([1,2,0],x)",
        "ldegree(x^3+x^2,x)",
    ]);
}

/// Bézout, Chinese remainders, series (CAS-483).
#[test]
fn bezout_and_series() {
    check(&[
        "egcd(x^2+1,x+1,x)",
        "gcdex(x^3-1,x^2-1,x)",
        "abcuv(x^2+1,x+1,x)",
        "chinrem([2,x],[3,x+1])",
        "chinrem([x+1,x^2+1],[x+2,x^2+2])",
        "divpc(1,1-x,5)",
        "divpc(x+1,x^2+1,4)",
        "pade(exp(x),x,4,2)",
        "pade(exp(x),x,5,3)",
    ]);
}

/// Factoring and forms (CAS-484).
#[test]
fn factoring_and_forms() {
    check(&[
        "cfactor(x^2+1)",
        "cfactor(x^4+4)",
        "canonical_form(2*x^2+4*x+1,x)",
        "canonical_form(x^2-6*x+2,x)",
        "genpoly(15,4,x)",
        "genpoly(1000,10,x)",
        "berlekamp_massey([0,1,1,2,3,5,8,13])",
        "berlekamp_massey([1,2,4,8])",
    ]);
}
