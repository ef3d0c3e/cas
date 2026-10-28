//! Statistics against giac (CAS-DOC-22 §6, CAS-PLAN-16): giac's answers and
//! this CAS's, element by element for lists, one value when their difference
//! expands to 0 or, for giac's decimals, is within 10⁻²⁰.

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

/// Whether two numbers agree to 10⁻²⁰, real and imaginary parts apart.
fn close(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let gap = Expr::call(Symbol::N, [a.clone() - b.clone(), Expr::from(30)]);
    [Symbol::RE, Symbol::IM].into_iter().all(|part| {
        let small = Expr::call(
            Symbol::LESS,
            [
                Expr::call(Symbol::ABS, [Expr::call(part, [gap.clone()])]),
                Expr::from(10).pow(-20),
            ],
        );
        ctx.eval(&small).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
    })
}

/// Whether two answers are one value: lists element by element, otherwise
/// a difference that expands to 0, or numbers within 10⁻²⁰.
fn same_value(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) || b.head() == Some(&list) {
        return a.head() == b.head()
            && a.args().len() == b.args().len()
            && a.args()
                .iter()
                .zip(b.args())
                .all(|(x, y)| same_value(ctx, x, y));
    }
    let difference = Expr::call(Symbol::EXPAND, [a.clone() - b.clone()]);
    ctx.eval(&difference).is_ok_and(|d| d == Expr::from(0)) || close(ctx, a, b)
}

/// giac's answers and this CAS's to `inputs`, one value each.
fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same_value(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Means, with frequencies, of columns, of symbols (CAS-294).
#[test]
fn means() {
    check(&[
        "mean([1,2,3,4])",
        "mean([1,2,3,4],[1,1,1,5])",
        "mean([[1,2],[3,4],[5,7]])",
        "mean([a,b,c])",
        "mean([1/2,1/3])",
        "mean([sqrt(2),1])",
        "mean([x,y],[p,1-p])",
        "mean([0.5,1.5])",
    ]);
}

/// Population variances and deviations, and the sample deviation
/// (CAS-294).
#[test]
fn deviations() {
    check(&[
        "variance([1,2,3,4])",
        "variance([1,2,3,4],[1,1,1,5])",
        "variance([[1,2],[3,4],[5,7]])",
        "variance([a,b])",
        "variance([a,b,c])",
        "variance([1/2,1/3])",
        "variance([sqrt(2),1])",
        "variance([1,2,3],[1/2,1/3,1/6])",
        "stddev([1,2,3,4])",
        "stddev([1,2],[1,3])",
        "stddev([[1,2],[3,4],[5,7]])",
        "stddev([5])",
        "stddevp([1,2,3,4])",
        "stddevp([1,2,3,4],[1,1,1,5])",
        "stddevp([[1,2],[3,4],[5,7]])",
        "stddevp([1,2,3],[2,2,2])",
    ]);
}
