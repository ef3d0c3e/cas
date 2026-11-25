//! `latex` and `TeX` against giac (CAS-615): the text of an expression's
//! answer in LaTeX, ours printed by the port of giac's `gen2markup`
//! (`cas_lang::markup`). Rows: giac's answer, its quotes removed, is ours,
//! whitespace collapsed (C0047). Departures: where our answer has another
//! shape or order than giac's evaluator leaves (constants first in sums,
//! our product order, relations kept as written), ours asserted and giac's
//! text in a comment; the printer is the same.

use cas_core::Expr;
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// The text of a string answer: giac's quotes removed and its doubled
/// inner quotes undone.
fn unquoted(answer: &str) -> Option<String> {
    let inner = answer.trim().strip_prefix('"')?.strip_suffix('"')?;
    Some(inner.replace("\"\"", "\""))
}

fn collapsed(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ")
}

/// Each input answers a string: giac's text is ours.
fn check_latex(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for input in inputs {
        let answer = giac.eval(input).unwrap_or_default();
        let mine = ours(input);
        let ok = match (unquoted(&answer), mine.as_ref().ok().and_then(Expr::as_str)) {
            (Some(theirs), Some(mine)) => collapsed(&theirs) == collapsed(mine),
            _ => false,
        };
        if !ok {
            wrong.push(format!("{input}: giac {answer}, ours {mine:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The plan's rows (R1.6 probes) and the lowering rows of
/// `cas_lang::markup`'s tests that print as giac prints its own answer;
/// `TeX` is `latex`; a `rootof` stays and prints as giac's.
#[test]
fn answers_as_giac_writes_them() {
    check_latex(&[
        "latex(1/2)",
        "latex(-1/2)",
        "latex(x^2)",
        "TeX(x^2)",
        "latex(sqrt(x))",
        "latex(sqrt(2)/2)",
        "latex(x/2)",
        "latex(2/x)",
        "latex(x*y)",
        "latex(2*x)",
        "latex(-x)",
        "latex(x-y)",
        "latex(-x-y)",
        "latex(1/x)",
        "latex(x^(-2))",
        "latex(sin(x))",
        "latex(cos(x)^2)",
        "latex(sin(x^2))",
        "latex(exp(x))",
        "latex(exp(2*x))",
        "latex(ln(x))",
        "latex(abs(x))",
        "latex(floor(x))",
        "latex(asin(x))",
        "latex(pi)",
        "latex(infinity)",
        "latex(i)",
        "latex(1+2*i)",
        "latex(3.14)",
        "latex(1.5e-7)",
        "latex(-2.5)",
        "latex(3.0)",
        "latex(0.001)",
        "latex(1e-5)",
        "latex(1.23e20)",
        "latex(1e100)",
        "latex(12345678901.0)",
        "latex(1.0/3)",
        "latex([1,2,3])",
        "latex([[1,2],[3,4]])",
        "latex(x=1)",
        "latex(x^2+y^2=1)",
        "latex(f(x))",
        "latex(f(x,y))",
        "latex(diff(f(x),x))",
        "latex(int(f(x),x))",
        r#"latex("abc")"#,
        "latex(alpha)",
        "latex(abc)",
        "latex(x_1)",
        "latex(a_n)",
        "latex(n!)",
        "latex(x^y^z)",
        "latex(a/(b*c))",
        "latex(2/3*x)",
        "latex(x*(-y))",
        "latex(1..2)",
        "latex(true)",
        "latex(undef)",
        "latex(rootof([1,0],[1,0,1,1]))",
    ]);
}

/// More rows of `giac_latex.tsv` (giac's own `latex()` texts): numbers,
/// complex numbers, fractions, strings, names, matrices, functions, big
/// operators, relations.
#[test]
fn more_of_giacs_table() {
    check_latex(&[
        "latex(a*b*c)",
        "latex(log(x))",
        "latex(e)",
        "latex(cos(x))",
        "latex(tan(x))",
        "latex(atan(x))",
        "latex(sinh(x))",
        "latex(2^x)",
        "latex(theta^2)",
        "latex(-2*x)",
        "latex(3*(-2))",
        "latex(123456789012345678901234567890)",
        "latex(0.5)",
        "latex(6.02e23)",
        "latex(1e-10)",
        "latex(-1.5e-7)",
        "latex(3*i)",
        "latex(-i)",
        "latex(2-3*i)",
        "latex(-3/4)",
        "latex(3/4*x)",
        "latex(x/(2*y))",
        "latex(1/(2*x))",
        "latex(-x/3)",
        "latex(y/x)",
        "latex(2/(3*x))",
        r#"latex("a_b")"#,
        r#"latex("50%")"#,
        "latex(beta)",
        "latex(omega)",
        "latex(abc_1)",
        "latex(x12)",
        "latex(alpha1)",
        "latex(x_1_2)",
        "latex([])",
        "latex([[1,-2],[x,y^2]])",
        "latex(x^(1/3))",
        "latex(sin(x)^3)",
        "latex(ln(x)^2)",
        "latex(exp(-x))",
        "latex(2*exp(x))",
        "latex(floor(x/2))",
        "latex(sin(2*x))",
        "latex(3*sin(x))",
        "latex(f(x)/2)",
        "latex((2*n)!)",
        "latex(diff(f(x),x,x))",
        "latex(int(f(x),x,0,1))",
        "latex(sum(f(k),k,1,n))",
        "latex(x>=1)",
        "latex(a..b)",
        "latex(false)",
        "latex(pi/2)",
        "latex(x*y*z^2)",
        "latex(e*x)",
        "latex(-x^2)",
        "latex(x^2*y^3)",
    ]);
}

/// Departures: ours, with giac's text in the comment. giac's evaluator
/// leaves another shape or order (sums ordered with the variable first,
/// `x<y` turned into `y>x`, `x!=y` decided, `ceil` and `binomial`
/// rewritten); its printer and ours write the same shape alike.
#[test]
fn departures() {
    let cases: &[(&str, &str)] = &[
        // Our sums put constants first: giac \sqrt{x+1}.
        ("latex(sqrt(x+1))", r"\sqrt{1+x}"),
        // giac \frac{x^{4}-1}{x^{2}+3}.
        ("latex((x^4-1)/(x^2+3))", r"\frac{-1+x^{4}}{3+x^{2}}"),
        // Ascending degree: giac 2 x^{2}-3 x+1.
        ("latex(2*x^2-3*x+1)", r"1-3 x+2 x^{2}"),
        // giac \frac1{x+1}.
        ("latex(1/(x+1))", r"\frac1{1+x}"),
        // giac \left(x+1\right)^{2}.
        ("latex((x+1)^2)", r"\left(1+x\right)^{2}"),
        // Our factor order: giac \left(x+1\right) \left(x-1\right).
        ("latex((x+1)*(x-1))", r"\left(-1+x\right) \left(1+x\right)"),
        // Our products put cos before sin: giac \sin x\cdot \cos x.
        ("latex(sin(x)*cos(x))", r"\cos x\cdot \sin x"),
        // Ours keeps Times[-1, Plus[1, x]] (C0029): giac -x-1.
        ("latex(-(x+1))", r"-\left(1+x\right)"),
        // giac's evaluator turns x<y into y>x: giac y>x.
        ("latex(x<y)", r"x<y"),
        // giac y\geq x.
        ("latex(x<=y)", r"x\leq y"),
        // giac decides x!=y: \mathrm{true}.
        ("latex(x!=y)", r"x\neq y"),
        // giac rewrites ceil: -\left\lfloor -x\right\rfloor.
        ("latex(ceil(x))", r"\left\lceil x\right\rceil "),
        // giac rewrites binomial into factorials:
        // \frac{n\mathclose{!}}{k\mathclose{!}\cdot \left(n-k\right)\mathclose{!}}.
        ("latex(binomial(n,k))", r"\binom{n}{k}"),
        // giac \left(x^{\frac{1}{3}}\right)^{2}.
        ("latex(x^(2/3))", r"x^{\frac{2}{3}}"),
        // Our product order: giac x \sqrt{2}.
        ("latex(x*sqrt(2))", r"\sqrt{2} x"),
        // Our closed form, factored: giac \frac{-n-1+\left(n+1\right)^{2}}{2}.
        ("latex(sum(k,k,1,n))", r"\frac{1}{2} n \left(1+n\right)"),
        // Our diff(e, x, 2) differentiates by x, then by 2 (held):
        // giac \frac{\mathrm{d}^{2}}{\mathrm{d}x^{2}}f\left(x\right).
        (
            "latex(diff(f(x),x,2))",
            r"\mathrm{derive}\left(\frac{\mathrm{d}}{\mathrm{d}x}f\left(x\right),2\right)",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, mine) in cases {
        let found = ours(input);
        if found.as_ref().ok().and_then(Expr::as_str) != Some(*mine) {
            wrong.push(format!("{input}: ours {found:?}, recorded {mine:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
