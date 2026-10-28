//! Residues and Fourier coefficients against giac (CAS-DOC-25 §8,
//! CAS-PLAN-19): values compared as numbers, at sample points of the
//! parameters, since the two write radicals and exponentials of `iπ`
//! apart (`(-4-4*i)*sqrt(2)/32` against `−e^(iπ/4)/4`).

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

/// `e` with each symbol of `values` put in.
fn at(e: &Expr, values: &[(&str, &str)]) -> Expr {
    let rules: Vec<(Expr, Expr)> = values
        .iter()
        .map(|(name, value)| (Expr::symbol(name), ours(value).expect("a number")))
        .collect();
    e.substitute(&rules)
}

/// Whether two answers agree at the sample points of the parameters.
fn same_value(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    [("1/3", "5/7"), ("-2", "3/2")]
        .iter()
        .all(|&(a_value, b_value)| {
            let values = [("a", a_value), ("b", b_value)];
            close(ctx, &at(a, &values), &at(b, &values))
        })
}

/// giac's answers and this CAS's to `inputs`, as values.
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

/// Residues at poles (CAS-315). Not at `inf`, where giac's sign is the
/// opposite of Mathematica's, which this CAS follows; not at branch
/// points, where giac answers 0 and this CAS stays.
#[test]
fn residues() {
    check(&[
        "residue(1/(z^2+1),z,i)",
        "residue(1/(z^2+1),z,-i)",
        "residue(exp(z)/z^3,z,0)",
        "residue(1/sin(z),z,0)",
        "residue(1/(z^2*(z-1)),z,0)",
        "residue(cot(z),z,pi)",
        "residue(tan(z),z,pi/2)",
        "residue(1/((z-a)*(z-b)),z,a)",
        "residue(1/(z^4+1),z,exp(i*pi/4))",
        "residue(1/(z^2+1)^2,z,i)",
        "residue(exp(z)/(z-1)^2,z,1)",
        "residue(1/(exp(z)-1),z,2*i*pi)",
        "residue(1/(z*(z-1)),z,1)",
        "residue(cos(z)/z^3,z,0)",
        "residue(z/(exp(z)-1)^2,z,0)",
        "residue(1/z,z=0)",
        "residue(z^2,z,1)",
        "residue(sin(z)/z,z,0)",
        "residue(exp(a*z)/(z^2+1),z,i)",
    ]);
}
