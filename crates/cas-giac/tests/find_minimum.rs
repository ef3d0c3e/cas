//! `find_minimum` against giac (CAS-608): closed forms first (C0006), our
//! ball holding the closed form (at 60 digits) with its radius at most
//! `eps`, giac's double within `10⁻⁶` of it (`10⁻⁴` for the `eps = 10⁻⁵`
//! row); departures, where `f` does not compile and ours stays
//! `FindMinimum[…]`, with giac's answer as caseval prints it (C0007).

use cas_core::{Context, Expr, Symbol, View};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the comparisons skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: comparisons skipped");
    }
    giac
}

/// Whether `|a − b| < bound`, `a − b` at 60 digits.
fn within(ctx: &Context, a: &Expr, b: &Expr, bound: &Expr) -> bool {
    let gap = Expr::call(
        Symbol::ABS,
        [Expr::call(
            Symbol::N,
            [a.clone() - b.clone(), Expr::from(60)],
        )],
    );
    let claim = Expr::call(Symbol::LESS, [gap, bound.clone()]);
    ctx.eval(&claim).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// The midpoint and radius of an answer exactly, and the radius as a
/// double; an exact number has radius 0.
fn parts(answer: &Expr) -> Option<(Expr, Expr, f64)> {
    if let View::Ball(b) = answer.view() {
        let mid = b.mid().checked_to_rational()?;
        let rad = b.rad().to_float()?;
        let exact = rad.checked_to_rational()?;
        return Some((Expr::from(mid), Expr::from(exact), rad.to_f64()));
    }
    answer
        .to_rational()
        .map(|q| (Expr::from(q), Expr::from(0), 0.0))
}

/// giac's probe rows (probe3.out) that compile, with their closed forms.
#[test]
fn closed_forms() {
    let ctx = Context::new();
    let root = "sqrt(2)*cos(acos(-sqrt(2)/4)/3-4*pi/3)";
    let mut wrong = Vec::new();
    for (input, closed, eps, giac_bound) in [
        ("find_minimum(cosh(x),-1,1)", "0", 1e-8, 6),
        ("find_minimum(cosh(x),-1,1,1e-5,20)", "0", 1e-5, 4),
        ("find_minimum(x^2-2*x,0,3)", "1", 1e-8, 6),
        ("find_minimum((x-1/3)^2,0,1)", "1/3", 1e-8, 6),
        ("find_minimum(x->x^4-3*x^2+x,-3,3)", root, 1e-8, 6),
        ("find_minimum(sin(x),0,10)", "3*pi/2", 1e-8, 6),
        ("find_minimum(abs(x-1/2),0,1)", "1/2", 1e-8, 6),
    ] {
        let Ok(c) = ours(closed) else {
            wrong.push(format!("{input}: closed form {closed} unread"));
            continue;
        };
        let found = ours(input).ok();
        match found.as_ref().and_then(parts) {
            Some((mid, radius, rad)) => {
                // The ball holds c: |mid − c| ≤ radius, up to 10⁻⁴⁰.
                let bound = Expr::from(10).pow(-40) + radius;
                if !(rad <= eps && within(&ctx, &mid, &c, &bound)) {
                    wrong.push(format!("{input}: ours {found:?}, closed form {closed}"));
                }
            }
            None => wrong.push(format!("{input}: ours {found:?}")),
        }
        if let Some(giac) = giac() {
            let theirs = giac.value(input);
            let tol = Expr::from(10).pow(-giac_bound);
            if !theirs.as_ref().is_ok_and(|t| within(&ctx, t, &c, &tol)) {
                wrong.push(format!("{input}: giac {theirs:?}, closed form {closed}"));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departures, giac's answers as caseval prints them (C0007): Airy, Bessel
/// and Gamma do not compile to ball operations, and ours stays
/// `FindMinimum[…]`.
#[test]
fn departures() {
    let mut wrong = Vec::new();
    for (input, answer) in [
        (
            "find_minimum(Airy_Ai(x+sin(x))+cos(x^2),-5,5)",
            "-3.96344565209",
        ),
        (
            "find_minimum(besselJ(x,2)/Gamma(x+1)+(x+1)^sin(x),0,10)",
            "4.82843934915",
        ),
    ] {
        if let Some(giac) = giac() {
            let text = giac.eval(input).unwrap_or_default();
            if text != answer {
                wrong.push(format!("{input}: giac prints {text}, recorded {answer}"));
            }
        }
        let found = ours(input);
        if found.as_ref().ok().and_then(Expr::head) != Some(&Expr::from(Symbol::FIND_MINIMUM)) {
            wrong.push(format!("{input}: ours {found:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
