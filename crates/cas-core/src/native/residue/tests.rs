//! Tests of `Residue` (CAS-DOC-25 §4, CAS-315): giac's residues at poles,
//! the forms of the call, analytic points, the point at infinity, and what
//! stays.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)`, both evaluated: within `10⁻³⁰` of each other, as
/// numbers, so that `1/(2i)` and `−i/2` agree.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let difference = Expr::call(Symbol::N, [found.clone() - read(expected), Expr::from(40)]);
        let gap = Expr::call(
            Symbol::LESS,
            [
                Expr::call(Symbol::ABS, [difference]),
                Expr::from(10).pow(-30),
            ],
        );
        if ctx.eval(&gap).unwrap().as_symbol() != Some(Symbol::TRUE) {
            wrong.push(format!("{input}: {found:?}, expected {expected}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's `residue` at poles: `1/(z^2+1)` at `i` (`1/(2*i)`), `exp(z)/z^3`
/// at 0 (`1/2`), `1/sin(z)` at 0 (1), `1/(z^2*(z-1))` at 0 (−1),
/// `cot(z)` at `π` (1), `1/((z-a)*(z-b))` at `a` (`1/(a-b)`),
/// `1/(z^2+1)^2` at `i` (`-i/4`), `exp(z)/(z-1)^2` at 1 (`e`),
/// `z/(exp(z)-1)^2` at 0 (1), `cos(z)/z^3` at 0 (`-1/2`), `1/(exp(z)-1)`
/// at `2iπ` (1, by CAS-313), and `1/(z^4+1)` at `exp(iπ/4)`
/// (`(-4-4*i)*sqrt(2)/32`).
#[test]
fn poles() {
    check(&[
        (
            "Residue[Power[Plus[1, Power[z, 2]], -1], z, I]",
            "Times[Complex[0, -1/2]]",
        ),
        ("Residue[Times[Exp[z], Power[z, -3]], z, 0]", "1/2"),
        ("Residue[Power[Sin[z], -1], z, 0]", "1"),
        (
            "Residue[Power[Times[Power[z, 2], Plus[-1, z]], -1], z, 0]",
            "-1",
        ),
        ("Residue[Cot[z], z, Pi]", "1"),
        (
            "Residue[Power[Plus[1, Power[z, 2]], -2], z, I]",
            "Complex[0, -1/4]",
        ),
        ("Residue[Times[Exp[z], Power[Plus[-1, z], -2]], z, 1]", "E"),
        ("Residue[Times[z, Power[Plus[-1, Exp[z]], -2]], z, 0]", "1"),
        ("Residue[Times[Cos[z], Power[z, -3]], z, 0]", "-1/2"),
        (
            "Residue[Power[Plus[-1, Exp[z]], -1], z, Times[Complex[0, 2], Pi]]",
            "1",
        ),
        (
            "Residue[Power[Plus[1, Power[z, 4]], -1], z, Exp[Times[Complex[0, 1/4], Pi]]]",
            "Times[Complex[-1, -1], Power[2, 1/2], 1/8]",
        ),
    ]);
    let ctx = Context::new();
    let found = ctx
        .eval(&read(
            "Residue[Power[Times[Plus[z, Times[-1, a]], Plus[z, Times[-1, b]]], -1], z, a]",
        ))
        .unwrap();
    let expected = ctx.eval(&read("Power[Plus[a, Times[-1, b]], -1]")).unwrap();
    assert_eq!(
        ctx.eval(&Expr::call(Symbol::TOGETHER, [found.clone() - expected]))
            .unwrap(),
        Expr::from(0),
        "{found:?}"
    );
}

/// The forms of the call, analytic points and removable singularities
/// (0, as giac's `residue(z^2,z,1)` and `residue(sin(z)/z,z,0)`), and the
/// point at infinity with Mathematica's sign.
#[test]
fn forms_and_points() {
    check(&[
        ("Residue[Power[z, -1], Equation[z, 0]]", "1"),
        ("Residue[Power[z, -1], List[z, 0]]", "1"),
        ("Residue[Power[z, -1], Rule[z, 0]]", "1"),
        ("Residue[Power[z, 2], z, 1]", "0"),
        ("Residue[Times[Sin[z], Power[z, -1]], z, 0]", "0"),
        ("Residue[c, z, 0]", "0"),
        ("Residue[Power[z, -1], z, Infinity]", "-1"),
        ("Residue[Power[z, -1], z, ComplexInfinity]", "-1"),
        ("Residue[Power[z, -2], z, Infinity]", "0"),
        (
            "Residue[Times[z, Power[Plus[1, Power[z, 2]], -1]], z, Infinity]",
            "-1",
        ),
        ("Residue[Power[z, 3], z, Infinity]", "0"),
        (
            "Residue[Times[Power[z, 2], Power[Plus[1, z], -1]], z, Infinity]",
            "-1",
        ),
    ]);
}

/// The residues of a rational function at its poles and at infinity add
/// up to 0.
#[test]
fn residue_theorem() {
    let ctx = Context::new();
    let f = "Times[Plus[1, Times[2, z], Power[z, 3]], Power[Times[Power[z, 2], Plus[-1, z], Plus[2, z]], -1]]";
    let total = read(&format!(
        "Plus[Residue[{f}, z, 0], Residue[{f}, z, 1], Residue[{f}, z, -2], Residue[{f}, z, Infinity]]"
    ));
    assert_eq!(ctx.eval(&total).unwrap(), Expr::from(0));
}

/// What stays: an essential singularity (giac refuses `exp(1/z)`), branch
/// points (giac's 0 for `sqrt(z)/z`), an unknown function, a variable that
/// is no symbol, and a hidden zero in the denominator (CAS-314).
#[test]
fn what_stays() {
    let ctx = Context::new();
    for input in [
        "Residue[Exp[Power[z, -1]], z, 0]",
        "Residue[Power[z, -1/2], z, 0]",
        "Residue[Times[Log[z], Power[z, -1]], z, 0]",
        "Residue[Log[z], z, 0]",
        "Residue[f[z], z, 0]",
        "Residue[Power[z, -1], Plus[1, z], 0]",
        "Residue[Power[Plus[z, Power[Sin[1], 2], Power[Cos[1], 2], -1], -1], z, 0]",
    ] {
        let e = ctx.eval(&read(input)).unwrap();
        assert_eq!(e.head(), read(input).head(), "{input}: {e:?}");
    }
}
