//! Differential equations against giac (CAS-DOC-16 §9, CAS-PLAN-10). Each
//! equation is a template in `Y`, standing for the unknown applied: this
//! CAS's general solution, printed into it, makes it 0 at a point where the
//! constants have values too, with a constant per order; where both choose
//! the constants alike, giac's solution is ours by value there.

use cas_core::{Context, Expr, Symbol, View};
use cas_giac::{Giac, ours};
use cas_lang::{Names, print};

/// The shared library, or `None` with a note: the comparisons skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: comparisons skipped");
    }
    giac
}

/// The symbols' values where answers are compared, the constants among
/// them.
const POINT: [(&str, i64, i64); 6] = [
    ("x", 2, 9),
    ("t", 5, 7),
    ("a", 3, 4),
    ("c_0", 7, 5),
    ("c_1", -2, 3),
    ("c_2", 5, 11),
];

/// Whether `proposition` evaluates to `True`.
fn holds(ctx: &Context, proposition: &Expr) -> bool {
    ctx.eval(proposition).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// The parts of a number, a complex ball's real and imaginary ones apart.
fn parts(e: &Expr) -> Vec<Expr> {
    match e.view() {
        View::ComplexBall(z) => vec![Expr::from(z.re().clone()), Expr::from(z.im().clone())],
        _ => vec![e.clone()],
    }
}

/// Whether two numbers agree: each part of `a − b` within `10⁻¹⁰ (1 + |b|)`.
fn close(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let Ok(difference) = ctx.eval(&(a.clone() - b.clone())) else {
        return false;
    };
    let scale = Expr::from(1) + Expr::call(Symbol::ABS, [b.clone()]);
    let bound = scale / Expr::from(10_000_000_000_i64);
    parts(&difference).into_iter().all(|p| {
        let gap = Expr::call(Symbol::ABS, [p]);
        holds(ctx, &Expr::call(Symbol::LESS, [gap, bound.clone()]))
    })
}

/// `e` approximated at [`POINT`].
fn at(ctx: &Context, e: &Expr) -> Option<Expr> {
    let rules: Vec<(Expr, Expr)> = POINT
        .iter()
        .map(|&(s, p, q)| (Expr::symbol(s), Expr::from(p) / Expr::from(q)))
        .collect();
    ctx.eval(&Expr::call(Symbol::N, [e.substitute(&rules)]))
        .ok()
}

/// Whether two answers are one: lists part by part, values close at
/// [`POINT`].
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) && b.head() == Some(&list) {
        return a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(a, b)| same(ctx, a, b));
    }
    match (at(ctx, a), at(ctx, b)) {
        (Some(a), Some(b)) => close(ctx, &a, &b),
        _ => false,
    }
}

/// An equation: its template in `Y` (`diff(Y, x) - Y` for `y′ = y`), its
/// variable, its order, and whether giac's constants are ours.
struct Case(&'static str, &'static str, usize, bool);

/// Each case's general solution satisfies its equation, with its constants;
/// those alike are giac's by value.
fn check(cases: &[Case]) {
    let ctx = Context::new();
    let giac = giac();
    let mut wrong = Vec::new();
    for &Case(template, v, order, alike) in cases {
        let unknown = format!("y({v})");
        let input = format!(
            "desolve({} = 0, {unknown})",
            template.replace('Y', &unknown)
        );
        let Ok(solution) = ours(&input) else {
            wrong.push(format!("{input}: unreadable"));
            continue;
        };
        // A nonlinear equation's answer is a list of branches.
        let branches = if solution.head() == Some(&Expr::from(Symbol::LIST)) {
            solution.args().to_vec()
        } else {
            vec![solution.clone()]
        };
        let satisfies = |branch: &Expr| {
            print(branch, Names::english()).is_ok_and(|text| {
                let residual = template.replace('Y', &format!("({text})"));
                ours(&residual)
                    .ok()
                    .and_then(|r| at(&ctx, &r))
                    .is_some_and(|r| close(&ctx, &r, &Expr::from(0)))
            })
        };
        let constants =
            |branch: &Expr| (0..order).all(|k| branch.contains(&Expr::symbol(&format!("c_{k}"))));
        if branches.is_empty() || !branches.iter().all(|b| satisfies(b) && constants(b)) {
            wrong.push(format!("{input}: {solution:?}"));
            continue;
        }
        if let (true, Some(giac)) = (alike, giac) {
            match giac.value(&input) {
                Ok(theirs) if same(&ctx, &solution, &theirs) => {}
                theirs => wrong.push(format!("{input}: ours {solution:?}, giac's {theirs:?}")),
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// First-order linear equations (CAS-259). giac writes `y′ + y = 1/x` with
/// `Ei`; its integral stays here.
#[test]
fn first_order_linear() {
    check(&[
        Case("diff(Y, x) - Y", "x", 1, true),
        Case("diff(Y, x) + Y - x", "x", 1, true),
        Case("diff(Y, x) - 2*Y/x", "x", 1, true),
        Case("diff(Y, x) - Y/x - x", "x", 1, true),
        Case("2*x*Y + (x^2 + 1)*diff(Y, x)", "x", 1, true),
        Case("diff(Y, t) + 2*t*Y", "t", 1, true),
        Case("diff(Y, x) + 2*x*Y - x", "x", 1, true),
        Case("diff(Y, x) - sin(x)", "x", 1, true),
        Case("diff(Y, x) - a*Y", "x", 1, true),
        Case("diff(Y, x) + Y - 1/x", "x", 1, false),
    ]);
}

/// Separable, Bernoulli, homogeneous and exact equations (CAS-260), lists
/// of branches. Alike where giac's constant is ours, `H(x) − G(y) = c_0`;
/// giac's others are parametrised otherwise.
#[test]
fn first_order_nonlinear() {
    check(&[
        Case("diff(Y, x) - Y^2", "x", 1, true),
        Case("diff(Y, x) - x*Y^2", "x", 1, true),
        Case("diff(Y, x) - x/Y", "x", 1, true),
        Case("diff(Y, x) - 1 - Y^2", "x", 1, true),
        Case("diff(Y, x) + Y - Y^2", "x", 1, false),
        Case("diff(Y, x) - sqrt(Y)", "x", 1, false),
        Case("diff(Y, x) - x*exp(Y)", "x", 1, true),
        Case("diff(Y, x) + Y/x - x^2*Y^3", "x", 1, true),
        Case("diff(Y, x) - (x + Y)/x", "x", 1, true),
        Case("diff(Y, x) - (Y^2 + x^2)/(x*Y)", "x", 1, false),
        Case("(3*x^2 + Y) + (x + 2*Y)*diff(Y, x)", "x", 1, false),
    ]);
}

/// Linear equations with constant coefficients (CAS-261). Alike where giac
/// orders the basis as the roots are here; it reverses a double root's and
/// some real ones', and writes `y⁗ = y` with combined constants.
#[test]
fn constant_coefficients() {
    check(&[
        Case("diff(Y, x, x) + Y", "x", 2, true),
        Case("diff(Y, x, x) + 2*diff(Y, x) + Y", "x", 2, false),
        Case("diff(Y, x, x) + diff(Y, x) + Y", "x", 2, true),
        Case("diff(Y, x, x, x) - Y", "x", 3, true),
        Case("diff(Y, x, x) - Y", "x", 2, false),
        Case("diff(Y, x, x, x, x) - Y", "x", 4, false),
        Case("diff(Y, x, x) + a^2*Y", "x", 2, false),
        Case("diff(Y, x, x) - a*Y", "x", 2, true),
        Case("diff(Y, x, x) + diff(Y, x)", "x", 2, true),
    ]);
}

/// Right sides (CAS-262). Alike where giac's particular solution has no
/// homogeneous part and its basis is ordered as ours.
#[test]
fn right_sides() {
    check(&[
        Case("diff(Y, x, x) - 3*diff(Y, x) + 2*Y - exp(x)", "x", 2, false),
        Case("diff(Y, x, x) + 4*Y - sin(x)", "x", 2, true),
        Case("diff(Y, x, x) - x", "x", 2, false),
        Case("diff(Y, x, x) - Y - x^2", "x", 2, false),
        Case("diff(Y, x, x) + Y - cos(x)", "x", 2, false),
        Case("diff(Y, x, x) + Y - x*exp(x)", "x", 2, true),
        Case("diff(Y, x, x) - 2*diff(Y, x) + Y - exp(x)", "x", 2, false),
        Case("diff(Y, x, x) + Y - sin(x) - x", "x", 2, false),
        Case("diff(Y, x, x) + Y - tan(x)", "x", 2, false),
    ]);
}
