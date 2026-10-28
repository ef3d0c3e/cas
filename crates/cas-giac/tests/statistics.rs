//! Statistics against giac (CAS-DOC-22 §6, CAS-PLAN-16): giac's answers and
//! this CAS's, element by element for lists, one value when their difference
//! expands to 0, when `==` holds (giac's decimals are balls of 53 bits, which
//! `==` takes as equal to the values they overlap, D3), or within 10⁻²⁰.
//! giac prints 12 digits of its decimals, so only those that are exact
//! compare: its medians of irrational values do not.

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
/// a difference that expands to 0, equal by `==`, or numbers within 10⁻²⁰.
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
    let equal = Expr::call(Symbol::EQUAL, [a.clone(), b.clone()]);
    ctx.eval(&difference).is_ok_and(|d| d == Expr::from(0))
        || ctx.eval(&equal).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
        || close(ctx, a, b)
}

/// giac's answers, made by `adjust`, and this CAS's to `inputs`, one value
/// each.
fn check_adjusted(inputs: &[&str], adjust: fn(Expr) -> Expr) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input).map(adjust), ours(input)) {
            (Ok(theirs), Ok(ours)) if same_value(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's answers and this CAS's to `inputs`, one value each.
fn check(inputs: &[&str]) {
    check_adjusted(inputs, |e| e);
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

/// Medians where giac's lower middle value is ours: odd counts, or half the
/// total not reached exactly (CAS-295).
#[test]
fn medians() {
    check(&[
        "median([3,1,2])",
        "median([1,2,3,5,10])",
        "median([1,2,3,4],[1,1,1,5])",
        "median([[1,2],[3,4],[5,7]])",
    ]);
}

/// Quantiles, the first and third quartiles (CAS-295).
#[test]
fn quantiles() {
    check(&[
        "quantile([1,2,3,4,5,6,7,8,9,10],0.25)",
        "quantile([1,2,3,4,5,6,7,8,9,10],1/4)",
        "quantile([1,2,3,4,5,6,7,8,9,10],0.95)",
        "quantile([3,1,2],[1,2,1],0.25)",
        "quantile([1,2,3,4,5,6,7,8],1/3)",
        "quantile([sqrt(2),1,pi],1/3)",
        "quartile1([1,2,3,4,5,6,7,8,9,10])",
        "quartile3([1,2,3,4,5,6,7,8,9,10])",
        "quartile3([1,5,2,8,4])",
        "quartile1([1,2,3,4,5,6,7,8],[1,1,1,1,1,1,1,1])",
    ]);
}

/// Quartiles of odd counts, where giac's median is ours; giac answers a
/// column matrix, this CAS a flat list (CAS-295).
#[test]
fn quartiles() {
    check_adjusted(
        &["quartiles([3,1,2])", "quartiles([1,2,3,4,5,6,7,8,9])"],
        |column| {
            let entries = column.args().iter().map(|row| match row.args() {
                [x] => x.clone(),
                _ => row.clone(),
            });
            Expr::call(Symbol::LIST, entries.collect::<Vec<_>>())
        },
    );
}

/// Population covariances and correlations (CAS-296).
#[test]
fn two_variables() {
    check(&[
        "covariance([1,2,3,4],[2,4,5,9])",
        "covariance([1,2,3,4],[2,4,5,9],[1,1,2,1])",
        "covariance([[1,2],[2,4],[3,5],[4,9]])",
        "covariance([1,2,3],[a,b,c])",
        "correlation([1,2,3,4],[2,4,5,9])",
        "correlation([1,2,3,4],[2,4,5,9],[1,1,2,1])",
        "correlation([[1,2],[2,4],[3,5],[4,9]])",
        "correlation([1,2,3],[2,4,6])",
    ]);
}
