//! Xcas's calculus helpers against giac (CAS-574, CAS-576 to CAS-580):
//! `ibpdv` and `ibpu` (one step, the bounds, the finish), `preval` and
//! `int` of one argument, `function_diff` (`fonction_derivee`) applied,
//! `revert`, `sum_riemann` and `wz_certificate`, `euler_lagrange`,
//! `jacobi_equation` and `conjugate_equation` (CAS-575), giac's answers read
//! back and compared with ours (`Giac::agrees`, or by value through
//! [`compare`]); and the departures, each with giac's answer as libgiac
//! prints it.

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

fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| giac.agrees(input).err().map(|e| format!("{e}")))
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's answer to `input`, written for this CAS: its second derivatives
/// `diff(f,v,2)` as `diff(f,v,v)` (ours would differentiate in 2), and its
/// constants `K_n`, numbered by a counter global to giac's session, as
/// `K_0`, `K_1`, … in the order they appear, as ours number them per call.
fn giac_text(giac: &Giac, input: &str) -> String {
    let mut text = giac.eval(input).unwrap();
    for v in ["x", "t"] {
        text = text.replace(&format!(",{v},2)"), &format!(",{v},{v})"));
    }
    let mut seen: Vec<String> = Vec::new();
    let mut out = String::new();
    let mut rest = text.as_str();
    while let Some(at) = rest.find("K_") {
        out.push_str(&rest[..at]);
        let digits = rest[at + 2..]
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(rest.len() - at - 2);
        let name = &rest[at..at + 2 + digits];
        let k = seen.iter().position(|s| s == name).unwrap_or_else(|| {
            seen.push(name.to_owned());
            seen.len() - 1
        });
        out.push_str("K_");
        out.push_str(&k.to_string());
        rest = &rest[at + 2 + digits..];
    }
    out.push_str(rest);
    out
}

fn zero(e: Expr) -> bool {
    Context::new().eval(&Expr::call(Symbol::SIMPLIFY, [e])).ok() == Some(Expr::from(0))
}

fn minus(a: &Expr, b: &Expr) -> Expr {
    Expr::call(
        Symbol::PLUS,
        [
            a.clone(),
            Expr::call(Symbol::TIMES, [Expr::from(-1), b.clone()]),
        ],
    )
}

/// Whether `a` and `b` are the same answer: lists of the same length
/// elementwise; equations `l = r` by `l − r` equal, or opposite when both
/// right sides are 0 (`E = 0` and `−E = 0` are one equation; giac's sign
/// varies with its `radsimp`); other expressions by a 0 difference.
fn equivalent(a: &Expr, b: &Expr) -> bool {
    let head = |e: &Expr| e.head().and_then(Expr::as_symbol);
    match (head(a), head(b)) {
        (Some(Symbol::LIST), Some(Symbol::LIST)) => {
            a.args().len() == b.args().len()
                && a.args().iter().zip(b.args()).all(|(a, b)| equivalent(a, b))
        }
        (Some(Symbol::EQUATION), Some(Symbol::EQUATION)) => {
            let ([la, ra], [lb, rb]) = (a.args(), b.args()) else {
                return false;
            };
            let (da, db) = (minus(la, ra), minus(lb, rb));
            let homogeneous = *ra == Expr::from(0) && *rb == Expr::from(0);
            zero(minus(&da, &db)) || homogeneous && zero(Expr::call(Symbol::PLUS, [da, db]))
        }
        _ => zero(minus(a, b)),
    }
}

/// giac's answer to `theirs` against ours to `input` (`input` writes
/// giac's `y'` as `diff(y(t),t)`).
fn compare(giac: &Giac, theirs: &str, input: &str) -> Result<(), String> {
    let text = giac_text(giac, theirs);
    match (ours(&text), ours(input)) {
        (Ok(t), Ok(o)) if equivalent(&o, &t) => Ok(()),
        (t, o) => Err(format!("{input}: ours {o:?}, giac's {text} ({t:?})")),
    }
}

/// `ibpdv` and `ibpu` (CAS-574): giac's examples.
#[test]
fn integration_by_parts() {
    check(&[
        "ibpdv(ln(x),x)",
        "ibpdv(ln(x),x,x,1,3)",
        "ibpdv(x*ln(x),x^2/2)",
        "ibpdv([x*ln(x),-1],0)",
        "ibpdv(ibpdv(ln(x),x,x,2,3),0,x,2,3)",
        "ibpdv(ln(y),y,y)",
        "ibpdv(x*exp(x),exp(x))",
        "ibpu(ln(x),ln(x))",
        "ibpu(ln(x),ln(x),x,1,3)",
        "ibpu(x*ln(x),ln(x))",
        "ibpu([x*ln(x),-1],0)",
        "ibpu(ibpu(ln(x),ln(x),x,2,3),0,x,2,3)",
        "ibpu(x^2,1)",
        "ibpu(x*exp(x),x)",
        "ibpu(exp(x)*sin(x),sin(x))",
    ]);
}

/// `preval` (CAS-578).
#[test]
fn values_at_bounds() {
    check(&[
        "preval(x^2-2,2,3)",
        "preval(y^2-2,2,3,y)",
        "preval(int(y,y),0,1,y)",
        "preval(x^2,x=1,2)",
        "preval(sin(x)/x,0,pi)",
        "preval(1/x,-1,1)",
        "preval(int(x),0,1)",
        "int(x)",
        "int(sin(x))",
        "int(y)",
    ]);
}

/// `function_diff` (CAS-577), applied: giac's functions print with a
/// hidden parameter `` ` x` `` that the language cannot read.
#[test]
fn derivatives_of_functions() {
    check(&[
        "function_diff(sin)(x)",
        "fonction_derivee(sin)(x)",
        "function_diff(sin+cos)(t)",
        "function_diff(2*exp)(t)",
        "function_diff(ln)(x)",
        "function_diff(ln)([x,y])",
        "function_diff(sin)([x,y])",
        "function_diff(unapply(y^3,y))(2)",
        "function_diff(f)",
    ]);
}

/// giac's fourth `wz_certificate` example, its summand over its sum, and
/// ours.
const WZ_FOURTH: &str = "wz_certificate((a+k)!*(b+k)!*(c+n-k-a-b-1)!/(k+1)!/(n-k)!/(c+k)!,\
    (a-1)!*(b-1)!*(c-a-b-1)!*(c+n-a)!*(c+n-b)!/(c-a-1)!/(c-b-1)!/(n+1)!/(c+n)!,n,k)";
const WZ_FOURTH_F: &str = "((a+k)!*(b+k)!*(c+n-k-a-b-1)!/(k+1)!/(n-k)!/(c+k)!)/\
    ((a-1)!*(b-1)!*(c-a-b-1)!*(c+n-a)!*(c+n-b)!/(c-a-1)!/(c-b-1)!/(n+1)!/(c+n)!)";
const WZ_FOURTH_OURS: &str = "-(a+k)*(b+k)/((a-c-n-1)*(b-c-n-1))";

/// `revert` (CAS-579): polynomials to their degree, other functions to the
/// default order 5 or the one given, Puiseux inverses.
#[test]
fn reversion() {
    check(&[
        "revert(x+x^2+x^4)",
        "revert(x-x^3/6)",
        "revert(1+x)",
        "revert(x)",
        "revert(2*x+x^2)",
        "revert(a*x+x^2)",
        "revert(x+x^7)",
        "revert(x+x^2,x,8)",
        "revert(x+x^2,x,3.5)",
        "revert(t+t^2,t)",
        "revert(sin(x))",
        "revert(tan(x))",
        "revert(exp(x)-1)",
        "revert(ln(1+x))",
        "revert(x^2)",
        "revert(x^2+x^3)",
    ]);
}

/// `sum_riemann` (CAS-576): giac's examples and what is not a Riemann sum.
#[test]
fn riemann_sums() {
    check(&[
        "sum_riemann(1/(n+k),[n,k])",
        "sum_riemann(n/(n^2+k^2),[n,k])",
        "sum_riemann(1/n*sin(k*pi/n),[n,k])",
        "sum_riemann(k^2/n^3,[n,k])",
        "sum_riemann(k/n^2,[n,k])",
        "sum_riemann(k^3/n^4,[n,k])",
        "sum_riemann(1/n,[n,k])",
        "sum_riemann(k,[n,k])",
        "sum_riemann(2/(n+k),[n,k])",
        "sum_riemann(n/(n+k)^2,[n,k])",
        "sum_riemann(1/(2*n+k),[n,k])",
        "sum_riemann(exp(k/n)/n,[n,k])",
        "sum_riemann(1/(n+k^2),[n,k])",
    ]);
}

/// `wz_certificate` (CAS-580): giac's examples, the argument forms, and 0
/// when there is none. giac's factored answers with a sign inside the
/// numerator are compared as differences.
#[test]
fn wz_certificates() {
    check(&[
        "wz_certificate((-1)^k*comb(n,k)/comb(k+a,k),a/(n+a),n,k)",
        "wz_certificate((-1)^k*comb(n,k)*comb(2k,k)*4^(n-k),comb(2n,n),n,k)",
        "normal(wz_certificate((-1)^(n-k)*comb(2n,k)^2,comb(2n,n),n,k)\
         -(-k^2+6*k*n+5*k-10*n^2-17*n-7)/(2*(k-2*n-2)^2))",
        "wz_certificate(comb(n,k)^2,comb(2n,n),n,k)",
        "wz_certificate(comb(n,k),2^n,n,k)",
        "wz_certificate(comb(n,k)/2^n,n,k)",
        "wz_certificate(comb(n,k)/2^n)",
        "wz_certificate(comb(n,k),2^n)",
        "wz_certificate(comb(m,j),2^m,m,j)",
        "wz_certificate(1/k,n,k)",
        "wz_certificate(k,1,n,k)",
        "wz_certificate(1,1,n,k)",
    ]);
}

/// `euler_lagrange` (CAS-575): giac's examples and the argument forms,
/// compared by value ([`compare`]): giac's `diff(f,t,2)` does
/// not read back, and its constants `K_n` count across its calls. The
/// first input is giac's, the second ours where the language lacks `y'`.
#[test]
fn euler_lagrange_equations() {
    let Some(giac) = giac() else {
        return;
    };
    let same = |input: &str| (input.to_owned(), input.to_owned());
    let cases = [
        same("euler_lagrange(sqrt(1+diff(y(x),x)^2),x,y)"),
        (
            "euler_lagrange(x^2*y'^2+y^2)".to_owned(),
            "euler_lagrange(x^2*diff(y(x),x)^2+y^2)".to_owned(),
        ),
        same("euler_lagrange(x^2*diff(y(x),x)^2+y(x)^2)"),
        same("euler_lagrange(x^2*diff(y(x),x)^2+y(x)^2,y(x))"),
        same("euler_lagrange(diff(y(t),t)^2-y(t)^2,t,y)"),
        same("euler_lagrange(diff(y(t),t)^2-y(t)^2,y(t))"),
        same("euler_lagrange(diff(y(t),t)^2-y(t)^2)"),
        (
            "euler_lagrange(sqrt(x'^2+y'^2),[x(t),y(t)])".to_owned(),
            "euler_lagrange(sqrt(diff(x(t),t)^2+diff(y(t),t)^2),[x(t),y(t)])".to_owned(),
        ),
        same("euler_lagrange(sqrt(diff(u(t),t)^2+diff(v(t),t)^2),t,[u,v])"),
        same("euler_lagrange(diff(y(x),x)^2+x*y(x))"),
        same("euler_lagrange(diff(y(x),x)^2+y(x)*sin(x))"),
        same("euler_lagrange(sqrt(1+diff(y(x),x)^2)/y(x))"),
        same("euler_lagrange(y(x)*sqrt(1+diff(y(x),x)^2))"),
        same("euler_lagrange(x*diff(y(x),x)^2)"),
        same("euler_lagrange(y(x)^2)"),
        same("euler_lagrange(exp(x)*(diff(y(x),x)^2+y(x)^2))"),
    ];
    let wrong: Vec<String> = cases
        .iter()
        .filter_map(|(theirs, input)| compare(giac, theirs, input).err())
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `jacobi_equation` (CAS-575): the equation by value, the solution up to
/// its constant (`DSolve`'s name may differ from giac's: their ratio is
/// free of `t`); with `L_y'y' = 0`, `true` or `false`.
#[test]
fn jacobi_equations() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        (
            "jacobi_equation(-1/2*y'(t)^2+y(t)^2/2,t,y,sin(t),h,0)",
            "jacobi_equation(-1/2*diff(y(t),t)^2+y(t)^2/2,t,y,sin(t),h,0)",
        ),
        (
            "jacobi_equation(-1/2*diff(y(t),t)^2+y(t)^2/2,y(t),sin(t),h,0)",
            "jacobi_equation(-1/2*diff(y(t),t)^2+y(t)^2/2,y(t),sin(t),h,0)",
        ),
        (
            "jacobi_equation(diff(y(t),t)^2,t,y,t,h,0)",
            "jacobi_equation(diff(y(t),t)^2,t,y,t,h,0)",
        ),
        (
            "jacobi_equation(diff(y(t),t)^2-y(t)^2,t,y,sin(t),h,0)",
            "jacobi_equation(diff(y(t),t)^2-y(t)^2,t,y,sin(t),h,0)",
        ),
    ];
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (theirs, input) in cases {
        let text = giac_text(giac, theirs);
        // giac's sequence `equation,solution`, split at its one comma
        // outside parentheses.
        let mut depth = 0;
        let comma = text.char_indices().find_map(|(i, c)| {
            match c {
                '(' | '[' => depth += 1,
                ')' | ']' => depth -= 1,
                ',' if depth == 0 => return Some(i),
                _ => {}
            }
            None
        });
        let found = ours(input).ok();
        let parts = found.as_ref().map(Expr::args);
        let ok = match (comma, parts) {
            (Some(i), Some([equation, solution])) => {
                let (eq, sol) = (ours(&text[..i]), ours(&text[i + 1..]));
                match (eq, sol) {
                    (Ok(eq), Ok(sol)) => {
                        let ratio = Expr::call(
                            Symbol::D,
                            [
                                Expr::call(
                                    Symbol::TIMES,
                                    [
                                        solution.clone(),
                                        Expr::call(Symbol::POWER, [sol, Expr::from(-1)]),
                                    ],
                                ),
                                Expr::symbol("t"),
                            ],
                        );
                        let ratio = Expr::call(Symbol::SIMPLIFY, [ratio]);
                        equivalent(equation, &eq) && ctx.eval(&ratio).ok() == Some(Expr::from(0))
                    }
                    _ => false,
                }
            }
            _ => false,
        };
        if !ok {
            wrong.push(format!("{input}: ours {found:?}, giac's {text}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
    check(&[
        "jacobi_equation(y(t)^2,t,y,t,h,0)",
        "jacobi_equation(y(t)*diff(y(t),t),t,y,t,h,0)",
    ]);
}

/// `conjugate_equation` (CAS-575): giac's example and others.
#[test]
fn conjugate_points() {
    check(&[
        "conjugate_equation(A*sin(t)+B*cos(t),[A,B],[1,0],t,0)",
        "conjugate_equation(A*exp(t)+B*exp(-t),[A,B],[1,1],t,0)",
        "conjugate_equation(A*t+B,[A,B],[1,0],t,0)",
        "conjugate_equation(A*sinh(t)+B*cosh(t),[A,B],[1,0],t,0)",
        "conjugate_equation(A*sin(t)+B*cos(t),[A,B],[2,3],t,1)",
    ]);
}

/// Where this CAS departs from giac: the input, giac's answer as it prints
/// it, ours as decided.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        // giac's `id` is the identity; ours has none, so `id` is an
        // unknown function.
        (
            "function_diff(sin+id)",
            "(` x`)->cos(` x`)+1",
            "function_diff(sin+id)",
        ),
        // A function of one parameter: ours ignores the arguments past it,
        // giac multiplies.
        ("function_diff(ln)(x,y)", "1/x/y", "1/x"),
        // giac's hidden parameter; ours is `x`.
        ("function_diff(sin)", "(` x`)->cos(` x`)", "(x)->cos(x)"),
        // giac's series of the product x*exp(x) to order 5 has a term
        // more; ours stops at x^5.
        (
            "revert(x*exp(x))",
            "x-x^2+3/2*x^3-8/3*x^4+125/24*x^5-54/5*x^6",
            "x-x^2+3/2*x^3-8/3*x^4+125/24*x^5",
        ),
        // giac takes sqrt(n^2*(1+x^2)) as n*sqrt(1+x^2); ours keeps the
        // factor with both n and x.
        (
            "sum_riemann(1/sqrt(n^2+k^2),[n,k])",
            "-ln(sqrt(2)-1)",
            "\"Probably not a Riemann sum\"",
        ),
        // giac's sum with the parameters a, b, c finds no antidifference
        // and answers 0; ours is a certificate, checked below.
        (WZ_FOURTH, "0", WZ_FOURTH_OURS),
        // giac's d/dt L_u' differentiates only the unknown of its equation
        // (u stays a constant in v's); ours in every unknown: 2v - u' = 0.
        (
            "euler_lagrange(diff(u(t),t)^2+u(t)*diff(v(t),t)+v(t)^2,t,[u,v])",
            "[-2*diff(u(t),t,2)+diff(v(t),t)=0,2*v(t)=0]",
            "[diff(v(t),t)-2*diff(u(t),t,t)=0,2*v(t)-diff(u(t),t)=0]",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, answer, decided) in cases {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        match (ours(input), ours(decided)) {
            (Ok(found), Ok(decided)) if found == decided => {}
            (found, _) => wrong.push(format!("{input}: ours {found:?}, not {decided}")),
        }
        if giac.agrees(input).is_ok() {
            wrong.push(format!("{input} agrees now"));
        }
    }
    // giac's `@` composition and `@@` powers: the language has neither.
    for (input, answer) in [
        (
            "function_diff(sq@sin+id)",
            "(` x`)->2*cos(` x`)*sin(` x`)+1",
        ),
        ("(function_diff @@3)(ln)('x')", "2/x^3"),
    ] {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        if ours(input).is_ok() {
            wrong.push(format!("{input}: ours reads it now"));
        }
    }
    // Our certificate R of the fourth: F(n+1,k) - F(n,k) = G(n,k+1) - G(n,k)
    // with G = R F(n,k-1), by giac's simplify, and F(n+1,k) != F(n,k).
    let f = WZ_FOURTH_F;
    let g = format!("({WZ_FOURTH_OURS})*subst({f},k=k-1)");
    let residual = format!("simplify((subst({f},n=n+1)-{f}-(subst({g},k=k+1)-{g}))/{f})");
    let printed = giac.eval(&residual).unwrap();
    if printed != "0" {
        wrong.push(format!("the fourth certificate: giac's residual {printed}"));
    }
    if giac
        .eval(&format!("simplify(subst({f},n=n+1)-{f})"))
        .unwrap()
        == "0"
    {
        wrong.push("the fourth: F(n+1,k) = F(n,k)".to_owned());
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
