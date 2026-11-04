//! Matrix decompositions III against giac (CAS-DOC-51 §7, CAS-PLAN-45).
//! giac's `schur` and `svd` are numeric: its singular values within a
//! relative 10⁻¹¹ of this CAS's exact ones, and its Schur form's diagonal,
//! sorted, of this CAS's (giac's vectors differ by signs and order).
//! `gauss_seidel_linsolve` by value; `jacobi_linsolve` is left out, giac's
//! answer not being the solution. `smith`'s `D` up to units (ours monic),
//! and `pcar_hessenberg` read back.

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

/// `N` of `e` at 30 digits.
fn numeric(ctx: &Context, e: Expr) -> Expr {
    ctx.eval(&Expr::call(Symbol::N, [e, Expr::from(30)]))
        .unwrap()
}

/// Whether `a` and `b` agree entry by entry within a relative 10⁻¹¹.
fn near(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let (t, o) = (
        leaves(&numeric(ctx, a.clone())),
        leaves(&numeric(ctx, b.clone())),
    );
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
            let bound = scale * Expr::from(10).pow(-11);
            ctx.eval(&Expr::call(Symbol::LESS, [gap, bound]))
                .ok()
                .and_then(|e| e.as_symbol())
                == Some(Symbol::TRUE)
        })
}

/// The diagonal of the matrix `m`, sorted, as numbers.
fn sorted_diagonal(ctx: &Context, m: &Expr) -> Expr {
    let d: Vec<Expr> = m
        .args()
        .iter()
        .enumerate()
        .map(|(i, row)| row.args()[i].clone())
        .collect();
    let d = numeric(ctx, Expr::call(Symbol::LIST, d));
    ctx.eval(&Expr::call(Symbol::SORT, [d])).unwrap()
}

/// Schur forms and singular values (CAS-496, CAS-497).
#[test]
fn schur_and_svd() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for a in [
        "[[1,2],[3,4]]",
        "[[1,2],[0,3]]",
        "[[2,1],[1,3]]",
        "[[2,0,0],[1,3,0],[4,5,6]]",
        "[[2,1,1],[1,2,1],[1,1,2]]",
    ] {
        let theirs = giac.value(&format!("schur({a})")).unwrap();
        let mine = ours(&format!("schur({a})")).unwrap();
        let (Some(t), Some(m)) = (theirs.args().get(1), mine.args().get(1)) else {
            wrong.push(format!("schur({a}): ours {mine:?}, giac's {theirs:?}"));
            continue;
        };
        if !near(&ctx, &sorted_diagonal(&ctx, t), &sorted_diagonal(&ctx, m)) {
            wrong.push(format!("schur({a}): ours {mine:?}, giac's {theirs:?}"));
        }
    }
    // Singular values with closed forms.
    for a in [
        "[[1,2],[3,4]]",
        "[[3,0],[4,5]]",
        "[[2,1],[1,3]]",
        "[[1,2],[2,4]]",
        "[[2,1,1],[1,2,1],[1,1,2]]",
    ] {
        let theirs = giac.value(&format!("[svd({a})]")).unwrap();
        let mine = ours(&format!("svd({a})")).unwrap();
        let (Some(t), Some(m)) = (theirs.args().get(1), mine.args().get(1)) else {
            wrong.push(format!("svd({a}): ours {mine:?}, giac's {theirs:?}"));
            continue;
        };
        if !near(&ctx, t, m) {
            wrong.push(format!("svd({a}): ours {mine:?}, giac's {theirs:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Gauss–Seidel, with and without relaxation (CAS-499).
#[test]
fn iterative_solutions() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "gauss_seidel_linsolve([[4,1],[1,3]],[1,2],1e-12)",
        "gauss_seidel_linsolve([[10,-1,2],[-1,11,-1],[2,-1,10]],[16,-14,23],1e-12)",
        "gauss_seidel_linsolve(1.2,[[4,1],[1,3]],[1,2],1e-12)",
    ] {
        let (Ok(theirs), Ok(mine)) = (giac.value(input), ours(input)) else {
            wrong.push(format!("{input}: not read"));
            continue;
        };
        if !near(&ctx, &theirs, &mine) {
            wrong.push(format!("{input}: ours {mine:?}, giac's {theirs:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Smith's form over ℚ\[x\] up to units, and `pcar_hessenberg` (CAS-500).
#[test]
fn smith_and_names() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for a in ["[[x,1],[0,x]]", "[[x-1,0],[0,x-2]]", "[[x^2,x],[x,1+x]]"] {
        let theirs = giac.value(&format!("[smith({a})]")).unwrap();
        let mine = ours(&format!("smith({a})")).unwrap();
        let (Some(t), Some(m)) = (theirs.args().get(1), mine.args().get(1)) else {
            wrong.push(format!("smith({a}): ours {mine:?}, giac's {theirs:?}"));
            continue;
        };
        // Each of giac's diagonal entries is a nonzero constant times ours.
        let units = leaves(t).iter().zip(leaves(m)).all(|(x, y)| {
            if *x == Expr::from(0) || y == Expr::from(0) {
                return *x == y;
            }
            let q = ctx
                .eval(&Expr::call(Symbol::TOGETHER, [x.clone() / y.clone()]))
                .unwrap();
            q.to_rational().is_some()
        });
        if !units {
            wrong.push(format!("smith({a}): ours {mine:?}, giac's {theirs:?}"));
        }
    }
    for input in [
        "pcar_hessenberg([[1,2],[3,4]])",
        "pcar_hessenberg([[2,1,0],[1,3,1],[0,1,4]])",
    ] {
        if let Err(e) = giac.agrees(input) {
            wrong.push(e.to_string());
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
