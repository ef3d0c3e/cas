//! Linear algebra against giac (CAS-DOC-12, CAS-PLAN-6): giac's answers,
//! its matrices `matrix[…]` read back as lists, are this CAS's.

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

/// The parameters' values where entries with parameters are compared.
const PARAMETERS: [(&str, i64, i64); 5] = [
    ("a", 3, 7),
    ("b", 5, 11),
    ("c", -2, 13),
    ("d", 7, 5),
    ("x", 2, 9),
];

/// Whether `proposition` evaluates to `True`.
fn holds(ctx: &Context, proposition: &Expr) -> bool {
    ctx.eval(proposition).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Whether two values are one: lists part by part, entries the same, equal
/// by `==`, or equal by `==` with the parameters at their values.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) && b.head() == Some(&list) {
        return a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(a, b)| same(ctx, a, b));
    }
    let equal = |a: &Expr, b: &Expr| holds(ctx, &Expr::call(Symbol::EQUAL, [a.clone(), b.clone()]));
    if a == b || equal(a, b) {
        return true;
    }
    let rules: Vec<(Expr, Expr)> = PARAMETERS
        .iter()
        .map(|&(s, p, q)| (Expr::symbol(s), Expr::from(p) / Expr::from(q)))
        .collect();
    let at = |e: &Expr| ctx.eval(&e.substitute(&rules)).ok();
    match (at(a), at(b)) {
        (Some(a), Some(b)) => a == b || equal(&a, &b),
        _ => false,
    }
}

/// giac's answers and this CAS's, the same by [`same`].
fn check_values(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The basics of CAS-218.
#[test]
fn basics() {
    check(&[
        "tran([[1,2,3],[4,5,6]])",
        "transpose([[1,2],[3,4]])",
        "tran([1,2,3])",
        "idn(3)",
        "identity(2)",
        "trace([[1,2],[3,4]])",
        "trace([[1,2,3],[4,5,6],[7,8,9]])",
        "dot([1,2,3],[4,5,6])",
        "dot([i,0],[1,0])",
        "dot([1,0],[i,0])",
        "dot([1+i,2],[3,4*i])",
        "dot([[1,2],[3,4]],[[1,0],[0,1]])",
        "cross([1,0,0],[0,1,0])",
        "cross([1,2],[3,4])",
        "cross([1,2,3],[4,5,6])",
    ]);
}

/// Determinants and inverses (CAS-219).
#[test]
fn determinants_and_inverses() {
    check_values(&[
        "det([[1,2],[3,4]])",
        "det([[a,b],[c,d]])",
        "det([[1,2,3],[4,5,6],[7,8,10]])",
        "det([[x,1],[1,x]])",
        "det([[2]])",
        "inv([[1,2],[3,4]])",
        "inv([[a,b],[c,d]])",
        "inv([[1,2,3],[0,1,4],[5,6,0]])",
        "inv([[2]])",
        "[[1,2],[3,4]]^-1",
        "[[1,2],[3,4]]^-2",
    ]);
}

/// Row reductions and linear systems (CAS-220).
#[test]
fn row_reductions() {
    check_values(&[
        "rref([[1,2,3],[4,5,6]])",
        "rref([[a,b],[c,d]])",
        "rref([[1,2,3],[4,5,6],[7,8,9]])",
        "rank([[1,2],[2,4]])",
        "rank([[1,2,3],[4,5,6],[7,8,9]])",
        "rank([[a,b],[c,d]])",
        "ker([[1,2,3],[2,4,6]])",
        "ker([[1,0],[0,1]])",
        "ker([[1,2,3],[4,5,6],[7,8,9]])",
        "ker([[a,b],[2*a,2*b]])",
        "image([[1,2,3],[2,4,6]])",
        "image([[1,0,1],[0,1,1]])",
        "linsolve([x+y=1,x-y=3],[x,y])",
        "linsolve([x+y=1,2*x+2*y=2],[x,y])",
        "linsolve([x+y=1,x+y=2],[x,y])",
        "linsolve([[1,1],[1,-1]],[1,3])",
        "linsolve([[1,2],[3,4]],[5,6])",
        "linsolve([[1,1],[1,1]],[1,2])",
        "linsolve([[a,1],[1,-1]],[1,b])",
    ]);
}

/// Characteristic polynomials and eigenvalues (CAS-221): giac's
/// eigenvalues, a sequence, as a list in any order.
#[test]
fn eigenvalues() {
    check_values(&[
        "charpoly([[1,2],[3,4]],x)",
        "charpoly([[1,2,3],[4,5,6],[7,8,10]],x)",
        "charpoly([[1,2],[3,4]],t)",
    ]);
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, matrix) in [
        ("eigenvals([[2,1],[1,2]])", "[[2,1],[1,2]]"),
        ("eigenvals([[0,-1],[1,0]])", "[[0,-1],[1,0]]"),
        ("eigenvals([[1,2],[3,4]])", "[[1,2],[3,4]]"),
        ("eigenvals([[1,1],[0,1]])", "[[1,1],[0,1]]"),
        (
            "eigenvals([[2,0,0],[0,2,0],[0,0,3]])",
            "[[2,0,0],[0,2,0],[0,0,3]]",
        ),
        (
            "eigenvals([[1,2,0],[2,1,0],[0,0,5]])",
            "[[1,2,0],[2,1,0],[0,0,5]]",
        ),
    ] {
        // giac's sequence, bracketed into a list.
        let theirs = giac.value(&format!("[{input}]"));
        let ours = ours(&format!("eigenvals({matrix})"));
        let (Ok(theirs), Ok(ours)) = (theirs, ours) else {
            wrong.push(format!("{input}: unreadable"));
            continue;
        };
        let mut left: Vec<&Expr> = theirs.args().iter().collect();
        let all = ours.args().len() == left.len()
            && ours.args().iter().all(|v| {
                left.iter()
                    .position(|t| same(&ctx, v, t))
                    .map(|k| left.swap_remove(k))
                    .is_some()
            });
        if !all {
            wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The columns of a matrix.
fn columns(m: &Expr) -> Vec<Vec<Expr>> {
    let n = m.args().first().map_or(0, |r| r.args().len());
    (0..n)
        .map(|j| m.args().iter().map(|r| r.args()[j].clone()).collect())
        .collect()
}

/// Whether two vectors are parallel: each of their 2 × 2 minors 0 by
/// value, within 10⁻⁹ of the vectors' scale.
fn parallel(ctx: &Context, u: &[Expr], v: &[Expr]) -> bool {
    let n = u.len();
    if n != v.len() {
        return false;
    }
    let minors: Vec<Expr> = (0..n)
        .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
        .map(|(i, j)| u[i].clone() * v[j].clone() - u[j].clone() * v[i].clone())
        .collect();
    let scale = Expr::call(
        Symbol::PLUS,
        u.iter()
            .chain(v)
            .map(|x| Expr::call(Symbol::ABS, [x.clone()]).pow(2))
            .collect::<Vec<_>>(),
    );
    minors.iter().all(|minor| {
        let small = Expr::call(
            Symbol::LESS,
            [
                Expr::call(Symbol::ABS, [minor.clone()]),
                scale.clone() / Expr::from(1_000_000_000_i64),
            ],
        );
        holds(ctx, &Expr::call(Symbol::N, [small]))
            || holds(
                ctx,
                &Expr::call(Symbol::EQUAL, [minor.clone(), Expr::from(0)]),
            )
    })
}

/// Eigenvectors (CAS-222): each of our columns parallel to one of giac's.
#[test]
fn eigenvectors() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in [
        "egv([[1,2],[2,1]])",
        "eigenvects([[2,1],[1,2]])",
        "eigenvects([[0,-1],[1,0]])",
        "eigenvects([[1,2],[3,4]])",
        "eigenvects([[2,1,0],[1,2,0],[0,0,3]])",
        "eigenvects([[2,0,0],[0,2,0],[0,0,3]])",
    ] {
        let (Ok(theirs), Ok(ours)) = (giac.value(input), ours(input)) else {
            wrong.push(format!("{input}: unreadable"));
            continue;
        };
        let mut left = columns(&theirs);
        let ours_columns = columns(&ours);
        let all = ours_columns.len() == left.len()
            && ours_columns.iter().all(|c| {
                left.iter()
                    .position(|t| parallel(&ctx, c, t))
                    .map(|k| left.swap_remove(k))
                    .is_some()
            });
        if !all {
            wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
