//! Number theory II against giac (CAS-DOC-39, CAS-PLAN-33): giac's answers
//! to `iquorem`, `smod`, `is_pseudoprime`, `cyclotomic`, `znorder`,
//! `znprimroot`, `iabcuv`, `pa2b2`, `isolve`, `dfc` and `dfc2f`, read back,
//! are this CAS's by value, lists part by part. giac's parameters `_Z0`
//! and `_N0` read as `n_0` and its residues `g % n` as `g`; its `isolve`
//! answers of another form, some with redundant parameters, are compared
//! as sets, and its `dfc` of irrationals by the numbers they give back.

use cas_core::{Context, Expr, Symbol, View};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// Whether `proposition` evaluates to `True`.
fn holds(ctx: &Context, proposition: &Expr) -> bool {
    ctx.eval(proposition).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// Whether two values are one: lists and equations part by part, the
/// same, or numbers within 10⁻¹⁰.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    if a == b {
        return true;
    }
    if a.head().is_some()
        && a.head() == b.head()
        && a.head()
            .and_then(Expr::as_symbol)
            .is_some_and(|h| matches!(h, Symbol::LIST | Symbol::EQUATION))
    {
        return a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(a, b)| same(ctx, a, b));
    }
    let Ok(difference) = ctx.eval(&Expr::call(Symbol::N, [a.clone() - b.clone()])) else {
        return false;
    };
    let parts = match difference.view() {
        View::ComplexBall(z) => vec![Expr::from(z.re().clone()), Expr::from(z.im().clone())],
        _ => vec![difference],
    };
    let tolerance = Expr::from(1) / Expr::from(10_000_000_000_i64);
    parts.into_iter().all(|p| {
        let gap = Expr::call(Symbol::ABS, [p]);
        holds(ctx, &Expr::call(Symbol::LESS, [gap, tolerance.clone()]))
    })
}

/// giac's answer to `input` in the language: its lists, its parameters
/// `_Z0` and `_N0` as `n_0`, and a residue `g % n` as `g`.
fn text(giac: &Giac, input: &str) -> String {
    let answer = giac
        .eval(input)
        .unwrap_or_else(|e| panic!("{input}: {e:?}"));
    let answer = answer
        .replace("list[", "[")
        .replace("_Z", "n_")
        .replace("_N", "n_");
    match answer.split_once(" % ") {
        Some((g, _)) => g.to_owned(),
        None => answer,
    }
}

/// This CAS's value of `text`, which must read.
fn value(text: &str) -> Expr {
    ours(text).unwrap_or_else(|e| panic!("{text}: {e:?}"))
}

/// giac's answers and this CAS's, the same by [`same`]; giac's read
/// through [`text`] when `raw`.
fn compare(inputs: &[&str], raw: bool) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        let theirs = if raw {
            ours(&text(giac, input))
        } else {
            giac.value(input)
        };
        match (theirs, ours(input)) {
            (Ok(theirs), Ok(ours)) if same(&ctx, &ours, &theirs) => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

fn check(inputs: &[&str]) {
    compare(inputs, false);
}

/// Division, pseudoprimes and cyclotomic polynomials (CAS-422).
#[test]
fn residues() {
    check(&[
        "iquorem(17,-5)",
        "iquorem(-17,5)",
        "iquorem(0,7)",
        "smod(8,5)",
        "smod(-7,4)",
        "smod(7,-4)",
        "smod(10,4)",
        "smod(-5,10)",
        "smod(5,10)",
        "is_pseudoprime(561)",
        "is_pseudoprime(97)",
        "is_pseudoprime(2)",
        "is_pseudoprime(1)",
        "is_pseudoprime(-7)",
        "cyclotomic(1)",
        "cyclotomic(2)",
        "cyclotomic(12)",
        "cyclotomic(105)",
        "cyclotomic(6,x)",
        "cyclotomic(15,x)",
    ]);
}

/// Orders and primitive roots (CAS-423); giac's `znprimroot` answers a
/// residue.
#[test]
fn orders() {
    check(&[
        "znorder(3,7)",
        "znorder(-2,9)",
        "znorder(4,6)",
        "znorder(10,1)",
        "znorder(0,1)",
        "znprimroot(2)",
        "znprimroot(4)",
        "znprimroot(8)",
        "znprimroot(12)",
    ]);
    compare(
        &[
            "znprimroot(1)",
            "znprimroot(7)",
            "znprimroot(18)",
            "znprimroot(25)",
        ],
        true,
    );
}

/// Diophantine equations (CAS-425): `iabcuv`, `pa2b2`, and `isolve` of
/// products, sums of two squares, systems without solutions and one
/// linear equation in two variables, giac's form.
#[test]
fn diophantine() {
    check(&[
        "iabcuv(21,28,14)",
        "pa2b2(97)",
        "pa2b2(13)",
        "pa2b2(2)",
        "isolve(x*y=6)",
        "isolve(x*y=-6)",
        "isolve(x*y=36)",
        "isolve(x*y=1)",
        "isolve(x^2+y^2=25)",
        "isolve(x^2+y^2=50)",
        "isolve(x^2+y^2=0)",
        "isolve(x^2+y^2=3)",
        "isolve(2*x+4*y=3)",
        "isolve([x+2*y+3*z=1,4*x+5*y+6*z=2])",
    ]);
    compare(
        &[
            "isolve(3*x+5*y=7)",
            "isolve(x-2*y=1)",
            "isolve(x+2*y=1)",
            "isolve(-x+2*y=1)",
            "isolve(4*x-6*y=2)",
            "isolve(-4*x-6*y=2)",
            "isolve(7*x+3*y=1)",
            "isolve(x+y=5)",
            "isolve(y+2*x=1)",
            "isolve([x+y+z=1,x-y=0])",
        ],
        true,
    );
}

/// The right sides of an answer `{x = …, …}`.
fn sides(answer: &Expr) -> Vec<Expr> {
    answer
        .args()
        .iter()
        .map(|eq| eq.args()[1].clone())
        .collect()
}

/// The parameters `n_0, …, n_7` that occur in `values`, in order.
fn parameters(values: &[Expr]) -> Vec<Expr> {
    (0..8)
        .map(|k| Expr::from(Symbol::new(&format!("n_{k}"))))
        .take_while(|p| values.iter().any(|v| v.contains(p)))
        .collect()
}

/// `values` at `t`, the parameters' values, expanded.
fn at(ctx: &Context, values: &[Expr], parameters: &[Expr], t: &[i64]) -> Vec<Expr> {
    let rules = Expr::call(
        Symbol::LIST,
        parameters
            .iter()
            .zip(t)
            .map(|(p, v)| Expr::call(Symbol::RULE, [p.clone(), Expr::from(*v)])),
    );
    values
        .iter()
        .map(|v| {
            let e = Expr::call(
                Symbol::EXPAND,
                [Expr::call(Symbol::REPLACE_ALL, [v.clone(), rules.clone()])],
            );
            ctx.eval(&e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
        })
        .collect()
}

/// Pell's equations (CAS-425): giac's general solution and this CAS's,
/// the same at `n_0` from 0 to 3.
#[test]
fn pell() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    for d in [2, 3, 13, 61] {
        let input = format!("isolve(x^2-{d}*y^2=1)");
        let theirs = sides(&value(&text(giac, &input)));
        let ours = sides(&value(&input));
        let n = [Expr::from(Symbol::new("n_0"))];
        for k in 0..4 {
            assert_eq!(
                at(&ctx, &theirs, &n, &[k]),
                at(&ctx, &ours, &n, &[k]),
                "{input} at {k}: ours {ours:?}, giac's {theirs:?}"
            );
        }
    }
}

/// An answer as `x₀ + Σ tₖ cₖ`: `x₀` and the columns `cₖ`, integers.
fn affine(ctx: &Context, answer: &Expr) -> (Vec<i64>, Vec<Vec<i64>>) {
    let values = sides(answer);
    let parameters = parameters(&values);
    let integers = |t: &[i64]| -> Vec<i64> {
        at(ctx, &values, &parameters, t)
            .iter()
            .map(|v| v.as_i64().unwrap_or_else(|| panic!("{answer:?}: {v:?}")))
            .collect()
    };
    let k = parameters.len();
    let x0 = integers(&vec![0; k]);
    let columns = (0..k)
        .map(|j| {
            let mut t = vec![0; k];
            t[j] = 1;
            integers(&t).iter().zip(&x0).map(|(a, b)| a - b).collect()
        })
        .collect();
    (x0, columns)
}

/// Whether `target` is an integer combination of at most two independent
/// `columns`, by Cramer's rule on a nonsingular minor.
fn in_lattice(columns: &[Vec<i64>], target: &[i64]) -> bool {
    let rows = target.len();
    let weights: Vec<i64> = match columns {
        [] => Vec::new(),
        [first] => {
            let Some(row) = (0..rows).find(|&row| first[row] != 0) else {
                return false;
            };
            if target[row] % first[row] != 0 {
                return false;
            }
            vec![target[row] / first[row]]
        }
        [first, second] => {
            let det = |top: usize, low: usize| first[top] * second[low] - second[top] * first[low];
            let minor = (0..rows)
                .flat_map(|top| (top + 1..rows).map(move |low| (top, low)))
                .find(|&(top, low)| det(top, low) != 0);
            let Some((top, low)) = minor else {
                return false;
            };
            let d = det(top, low);
            let one = target[top] * second[low] - second[top] * target[low];
            let two = first[top] * target[low] - target[top] * first[low];
            if one % d != 0 || two % d != 0 {
                return false;
            }
            vec![one / d, two / d]
        }
        _ => panic!("more than two parameters: {columns:?}"),
    };
    (0..rows).all(|row| {
        let sum: i64 = columns.iter().zip(&weights).map(|(c, w)| c[row] * w).sum();
        sum == target[row]
    })
}

/// Linear equations whose giac answers take another form, some with
/// redundant parameters (CAS-425): the same sets, as this CAS's answer
/// solves the equations identically, and giac's particular solution and
/// columns lie in its lattice.
#[test]
fn linear_sets() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let cases = [
        "x+2*y+3*z=1",
        "2*x+3*y+5*z=1",
        "6*x+10*y+15*z=1",
        "2*x+4*y+6*z=2",
        "3*x+2*y+z=1",
        "[x+y=1,x+y+z=3]",
    ];
    for equations in cases {
        let input = format!("isolve({equations})");
        let theirs = value(&text(giac, &input));
        let ours = value(&input);
        let names =
            |a: &Expr| -> Vec<Expr> { a.args().iter().map(|e| e.args()[0].clone()).collect() };
        assert_eq!(names(&theirs), names(&ours), "{input}");
        // Each equation, its sides' difference with this CAS's answer
        // substituted, expands to 0.
        let rules = Expr::call(
            Symbol::LIST,
            ours.args()
                .iter()
                .map(|eq| Expr::call(Symbol::RULE, [eq.args()[0].clone(), eq.args()[1].clone()])),
        );
        let system = value(&format!("[{}]", equations.trim_matches(['[', ']'])));
        for eq in system.args() {
            let [l, r] = eq.args() else {
                panic!("{eq:?}");
            };
            let gap = Expr::call(
                Symbol::EXPAND,
                [Expr::call(
                    Symbol::REPLACE_ALL,
                    [l.clone() - r.clone(), rules.clone()],
                )],
            );
            let gap = ctx.eval(&gap).unwrap();
            assert_eq!(gap, Expr::from(0), "{input}: ours {ours:?} in {eq:?}");
        }
        let (x0, columns) = affine(&ctx, &ours);
        let (y0, others) = affine(&ctx, &theirs);
        let offset: Vec<i64> = y0.iter().zip(&x0).map(|(a, b)| a - b).collect();
        assert!(
            in_lattice(&columns, &offset),
            "{input}: giac's {theirs:?} starts outside ours {ours:?}"
        );
        for c in &others {
            assert!(
                in_lattice(&columns, c),
                "{input}: giac's {theirs:?} reaches outside ours {ours:?}"
            );
        }
    }
}

/// Continued fractions (CAS-424): of rationals the same quotients; giac's
/// of irrationals take their own forms, a period started late or a last
/// complete quotient, so the numbers they give back, and the quotients
/// before giac's last.
#[test]
fn fractions() {
    check(&[
        "dfc(355/113)",
        "dfc(-7/3)",
        "dfc2f([3,7,16])",
        "dfc2f([1,2,[2]])",
        "dfc2f([2,[1,1,1,4]])",
        "dfc2f([[1]])",
    ]);
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    for (input, x) in [
        ("dfc(sqrt(2),5)", "sqrt(2)"),
        ("dfc(sqrt(7),10)", "sqrt(7)"),
        ("dfc((1+sqrt(5))/2,4)", "(1+sqrt(5))/2"),
    ] {
        let periodic = text(giac, input);
        let theirs = value(&format!("dfc2f({periodic})"));
        let ours = value(&format!("dfc2f(dfc({x}))"));
        let x = value(x);
        assert_eq!(theirs, x, "{input}: giac's {periodic}");
        assert_eq!(ours, x, "{input}: ours");
    }
    for (input, n) in [("dfc(pi,6)", 6), ("dfc(exp(1),8)", 8)] {
        let theirs = giac.value(input).unwrap();
        let ours = value(input);
        assert!(
            same(
                &ctx,
                &ours,
                &Expr::call(Symbol::LIST, theirs.args()[..n].to_vec())
            ),
            "{input}: ours {ours:?}, giac's {theirs:?}"
        );
    }
}
