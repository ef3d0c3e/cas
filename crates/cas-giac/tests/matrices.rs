//! Matrix decompositions and functions against giac (CAS-DOC-30 §5,
//! CAS-PLAN-24): entries compared after `Together` of the difference.
//! `jordan` only where giac's change of basis is this CAS's; `exp` of a
//! matrix is `expm` here, giac's `exp`. giac's `lu`, `qr` and `jordan`
//! return sequences, read in brackets as lists.

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

/// Whether two answers agree once put together, lists element by element.
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

#[test]
fn matrices() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, theirs_input) in [
        ("lu([[1,2],[3,4]])", "[lu([[1,2],[3,4]])]"),
        ("lu([[0,1],[1,1]])", "[lu([[0,1],[1,1]])]"),
        ("qr([[3,0],[4,5]])", "[qr([[3,0],[4,5]])]"),
        ("qr([[1,1],[1,0],[0,1]])", "[qr([[1,1],[1,0],[0,1]])]"),
        ("cholesky([[4,2],[2,3]])", "cholesky([[4,2],[2,3]])"),
        ("cholesky([[2,1],[1,2]])", "cholesky([[2,1],[1,2]])"),
        ("jordan([[1,1],[0,1]])", "[jordan([[1,1],[0,1]])]"),
        ("pmin([[1,1],[0,1]],x)", "pmin([[1,1],[0,1]],x)"),
        ("expm([[0,1],[0,0]])", "exp([[0,1],[0,0]])"),
        ("expm([[1,0],[0,2]])", "exp([[1,0],[0,2]])"),
        ("matpow([[1,1],[0,1]],n)", "matpow([[1,1],[0,1]],n)"),
        ("matpow([[2,1],[0,2]],n)", "matpow([[2,1],[0,2]],n)"),
        ("diag([1,2,3])", "diag([1,2,3])"),
        ("diag([[1,2],[3,4]])", "diag([[1,2],[3,4]])"),
        ("hilbert(3)", "hilbert(3)"),
        ("vandermonde([1,2,3])", "vandermonde([1,2,3])"),
        ("l2norm([3,4])", "l2norm([3,4])"),
        (
            "frobenius_norm([[1,2],[3,4]])",
            "frobenius_norm([[1,2],[3,4]])",
        ),
    ] {
        match (giac.value(theirs_input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
