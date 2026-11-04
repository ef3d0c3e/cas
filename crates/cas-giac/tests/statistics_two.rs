//! Statistics II against giac (CAS-DOC-50 §4, CAS-PLAN-44): giac's exact
//! answers read back are this CAS's; its decimals (`frequencies`, the ends
//! of `classes`) entry by entry within 10⁻¹¹, `a .. b` read as `Range[a,
//! b]`.

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

/// The leaves of nested lists and ranges, in order.
fn leaves(e: &Expr) -> Vec<Expr> {
    let nested = [Symbol::LIST, Symbol::RANGE];
    if e.head()
        .and_then(Expr::as_symbol)
        .is_some_and(|h| nested.contains(&h))
    {
        e.args().iter().flat_map(leaves).collect()
    } else {
        vec![e.clone()]
    }
}

fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        let (Ok(theirs), Ok(mine)) = (giac.value(input), ours(input)) else {
            wrong.push(format!("{input}: not read"));
            continue;
        };
        let (t, o) = (leaves(&theirs), leaves(&mine));
        let same = t.len() == o.len()
            && t.iter().zip(&o).all(|(a, b)| {
                if a == b {
                    return true;
                }
                let numeric = |e: &Expr| {
                    ctx.eval(&Expr::call(Symbol::N, [e.clone(), Expr::from(30)]))
                        .unwrap()
                };
                let gap = Expr::call(Symbol::ABS, [numeric(a) - numeric(b)]);
                let small = Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-11)]);
                ctx.eval(&small).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
            });
        if !same {
            wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Frequencies and classes (CAS-491).
#[test]
fn frequencies_and_classes() {
    check(&[
        "frequencies([1,2,2,3,3,3])",
        "frequencies([3,1,2,2])",
        "classes([1,2,3,4,5,6,7],0,3)",
        "classes([0.5,1.5,2.5,2.7],0,1)",
        "classes([1,2,3,4,5,6,7])",
    ]);
}

/// Moments, windows and the other names (CAS-492).
#[test]
fn moments_and_names() {
    check(&[
        "covariance_correlation([1,2,3],[2,4,7])",
        "covariance_correlation([[1,2],[2,4],[3,7]])",
        "rms([3,4])",
        "rms([1,2,3,4])",
        "moving_average([1,2,3,4,5],3)",
        "moving_average([1,2,3,4,5,6],2)",
        "moyenne([1,2,4])",
        "ecart_type([1,2,4])",
        "ecart_type_population([1,2,4])",
        "stdDev([1,2,4])",
    ]);
}
