//! Xcas's expression commands against giac, second part (CAS-568 to
//! CAS-571): `lname`, `indets`, `lvar`, `algvar`, `reorder`, `part`,
//! `eval`, `evalm`, `evalb`, `epsilon2zero`, `isinf`, `isnan`, `compare`,
//! `equal2diff`, `equal2list`, `list2exp` and `exp2list`, and `linstep`,
//! `simplifyDirac` and `periodic` (CAS-570, CAS-572, CAS-573), giac's
//! answers read back and compared with ours (`Giac::agrees`), lists of
//! kernels and identifiers as sets where giac's order of operands is not
//! ours; and the departures, each with giac's answer as libgiac prints it.

use cas_core::{Expr, Symbol};
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

/// `e` with the entries of each innermost list sorted: lists compared as
/// sets, the order of the rows of `algvar` kept.
fn sorted(e: &Expr) -> Expr {
    let is_list = |e: &Expr| e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST);
    if !is_list(e) {
        return e.clone();
    }
    let mut args: Vec<Expr> = e.args().iter().map(sorted).collect();
    if !args.iter().any(is_list) {
        args.sort_by_key(|x| format!("{x:?}"));
    }
    Expr::call(Symbol::LIST, args)
}

/// giac's answer and ours, lists compared as sets.
fn check_sets(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for input in inputs {
        let theirs = giac.value(input).unwrap();
        let ours = ours(input).unwrap();
        if sorted(&theirs) != sorted(&ours) {
            wrong.push(format!("{input}: giac {theirs:?}, ours {ours:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `lname`, `indets`, `lvar`, `algvar` and `reorder` (CAS-568).
#[test]
fn identifiers_and_kernels() {
    check(&[
        "lname(exp(x)*2*sin(y))",
        "indets(x*y+z)",
        "lname(x*y+z)",
        "lname(f(x,y))",
        "lvar(exp(x)*2*sin(y))",
        "lvar(sqrt(x)+x^(1/3)*y)",
        "lvar(x^2+y)",
        "lvar(1/(x+1))",
        "lvar(f(x,y))",
        "lvar([x,sin(y)])",
        "lvar(x^n)",
        "lvar(2^x)",
        "algvar(sqrt(x)+y)",
        "algvar(sqrt(x)+sqrt(y)+z)",
        "algvar(x^(1/3)+y)",
        "algvar(x+y)",
        "algvar(sqrt(x+sqrt(y)))",
        "algvar(x^(1/7)+y)",
        "algvar(x^(1/4)+y)",
        "algvar(x^(1/3)+y^(1/4))",
        "reorder(x^2+2*x+y^2,[y,x])",
        "reorder(a*x+b*y+c,[y,x])",
        "reorder(sin(x)+cos(y),[y,x])",
        "reorder((x+1)^2,[x])",
        "reorder(1/(x+1)+y,[y,x])",
        "reorder(3,[x])",
    ]);
    // giac's operands in its order, ours canonical: the same sets.
    check_sets(&[
        "lname(sin(a*x)+b^2)",
        "lname(pi*x+e^y)",
        "lvar(exp(x)*2*sin(y)+ln(x))",
        "lvar(sin(x)^2+cos(x))",
        "lvar(exp(2*x)+exp(x))",
    ]);
}

/// `part` (CAS-570), where its operands are giac's.
#[test]
fn parts() {
    check(&[
        "part(exp(x)*sin(x),1)",
        "part(x^2,1)",
        "part(x^2,2)",
        "part(sin(x),1)",
        "part(x+y+z,3)",
    ]);
    // A sum's operands in giac's order and in ours.
    check_sets(&["[part(sin(x)+cos(x),1),part(sin(x)+cos(x),2)]"]);
}

/// `eval`, `evalm`, `evalb`, `epsilon2zero`, `isinf`, `isnan` and
/// `compare` (CAS-569).
#[test]
fn evaluations_and_tests() {
    check(&[
        "eval(2*sin(pi))",
        "eval(x+x)",
        "eval(x+x,1)",
        "evalm([[1,2],[3,4]]^2)",
        "evalm(x)",
        "evalb(1<2)",
        "evalb(2>3)",
        "evalb(a==2)",
        "evalb(x==x)",
        "evalb(1=1)",
        "evalb(a<b)",
        "evalb(pi>3)",
        "evalb(sqrt(2)>1)",
        "evalb(sqrt(2)+pi>a)",
        "evalb(1<2 and 3<4)",
        "evalb(0)",
        "evalb(2)",
        "evalb([1<2,2<1])",
        "evalb(x+1==1+x)",
        "evalb(x<x+1)",
        "epsilon2zero(1e-13+x+5)",
        "epsilon2zero(1e-9*x+3.0)",
        "epsilon2zero([1e-13,2.5,1e-11*x])",
        "epsilon2zero(1e-13)",
        "epsilon2zero(9e-13)",
        "epsilon2zero(1.1e-12)",
        "epsilon2zero(1e-13+2e-11*i)",
        "epsilon2zero(x+1e-13*i)",
        "epsilon2zero(1e-5+x,1e-4)",
        "epsilon2zero(1e-5+x,1e-6)",
        "isinf(-inf)",
        "isinf(inf)",
        "isinf(+infinity)",
        "isinf(infinity)",
        "isinf(1)",
        "isinf(x)",
        "isinf(undef)",
        "isnan(undef)",
        "isnan(0)",
        "isnan(inf)",
        "isnan(x)",
        "isnan([undef,1])",
        "isnan([1,undef])",
        "compare(1,2)",
        "compare(2,1)",
        "compare(1,1)",
        "compare(-1,1)",
        "compare(1,-1)",
        "compare(1.0,2)",
        "compare(2,1.0)",
        "compare(1.5,1)",
        "compare(1.5,2.5)",
        "compare(\"ab\",\"cd\")",
        "compare(\"ab\",\"ab\")",
        "compare(1,\"ab\")",
        "compare(x,y)",
        "compare(y,x)",
        "compare(1,x)",
        "compare(x,1)",
        "compare(0,x)",
        "compare(x,0)",
        "compare(1/2,1)",
        "compare(1/2,1/3)",
        "compare(2^40,3)",
        "compare(3,2^40)",
        "compare(i,1)",
        "compare(1+i,2)",
        "compare(1+2*i,2+i)",
        "compare([1,2],[1,3])",
        "compare([1,2],[3])",
        "compare(x^2,x)",
        "compare(sin(x),x)",
        "compare(sin(x),sin(y))",
        "compare(x+1,x-1)",
        "compare(x+y,x)",
    ]);
}

/// `equal2diff`, `equal2list` (bracketed: giac prints a sequence, C0023),
/// `list2exp` and `exp2list` (CAS-571).
#[test]
fn equations_and_lists() {
    check(&[
        "equal2diff(x=2)",
        "equal2diff(x=y+1)",
        "equal2diff([x=1,y=2])",
        "[equal2list(x=2)]",
        "[equal2list(x+1=y)]",
        "equal2list(x)",
        "equal2list([x=1,y=2])",
        "list2exp([[3,9],[- 1,1]],[x, y])",
        "list2exp(solve(x^2-1=0,x),x)",
        "list2exp([1,2],x)",
        "list2exp([1],x)",
        "list2exp([[1,2]],[x,y])",
        "list2exp([[1,2],[3,4]],[x,y])",
        "list2exp([1,2],[x,y])",
        "list2exp([],x)",
        "exp2list((x=2) or (x=0))",
        "exp2list((x=3 and y=9) or (x=-1 and y=1) )",
        "exp2list(x=2)",
        "exp2list((x=1) or (x=2) or (x=3))",
        "exp2list([x=1,x=2])",
        "exp2list([x=1,y=2])",
        "exp2list(0)",
        "exp2list(x)",
    ]);
}

/// `linstep`, `simplifyDirac` and `periodic` (CAS-570, CAS-572, CAS-573).
/// giac's sums of steps and deltas read back in our canonical order; its
/// coefficients `series(c,x,z,1)` of a split delta read back as `c`.
#[test]
fn steps_deltas_periodic() {
    check(&[
        "linstep(Heaviside(x-2)*Heaviside(3-x),x)",
        "linstep(piecewise(x<-1,0,x<1,x+1,2-(x-1)^2),x)",
        "linstep(abs(abs(x)-abs(x-1)))",
        "linstep(Heaviside(x),x)",
        "linstep(x*Heaviside(x-1),x)",
        "linstep(piecewise(x<0,-x,x),x)",
        "linstep(Heaviside(x-1)*Heaviside(x-2),x)",
        "linstep(abs(x),x)",
        "linstep(Heaviside(x^2-1),x)",
        "linstep(abs(x^2-1),x)",
        "linstep(sign(x-1)*x,x)",
        "linstep(Heaviside(2*x-1)*Heaviside(x),x)",
        "linstep(Heaviside(x)^2,x)",
        "linstep(Heaviside(t-1)*t,t)",
        "linstep(abs(x)*abs(x-1),x)",
        "simplifyDirac(Dirac(a*x),x)",
        "simplifyDirac(Dirac(x,1)*x)",
        "simplifyDirac(Dirac(x^2-x-6))",
        "simplifyDirac(Dirac(2*x),x)",
        "simplifyDirac(Dirac(x)*x)",
        "simplifyDirac(Dirac(x-1)*x^2)",
        "simplifyDirac(Dirac(x)*sin(x))",
        "simplifyDirac(Dirac(x,1)*x,x)",
        "simplifyDirac(Dirac(-x))",
        "simplifyDirac(Dirac(x^2-1),x)",
        "simplifyDirac(Dirac(x-1,2)*x^3,x)",
        "simplifyDirac(Dirac(2*x-2,1)*x^2,x)",
        "simplifyDirac(Dirac(x^2+1),x)",
        "simplifyDirac(Dirac(x)*cos(x)+sin(x)^2+cos(x)^2,x)",
        "simplifyDirac(x^2,x)",
        "simplifyDirac(Dirac(t-1)*t^2,t)",
        "simplifyDirac(Dirac(x)/x,x)",
        "simplifyDirac(diff(Dirac(x),x)*x,x)",
        "periodic(x^2,x=-1..1)",
        "periodic(x,x=0..1)",
        "periodic(x^2,x,-1,1)",
        "periodic(sin(x),x=0..pi)",
        "periodic(3,x=0..1)",
        "periodic(y,x=0..1)",
        "periodic(x*y,x,0,2)",
    ]);
}

/// giac's answers this CAS departs from, as libgiac prints them (C0015),
/// with ours as decided.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        // Our sum is cos(x) + exp(x)*sin(x), canonical order: its first
        // operand is cos(x).
        ("part(part(exp(x)*sin(x)+cos(x),1),1)", "exp(x)", "x"),
        // Part 0 is the head itself here, giac's its name as a string.
        ("part(x^2,0)", "\"^\"", "pow"),
        ("part(sin(x),0)", "\"sin\"", "sin"),
        // giac's f(x,y) applies `of` to f and the sequence x,y.
        ("part(f(x,y),2)", "x,y", "y"),
        // An index out of range or an atom leaves Part, as Part takes lists.
        (
            "part(x+y,5)",
            "GIAC_ERROR: part(x+y,5) \n Error: Bad Argument Value",
            "part(x+y,5)",
        ),
        (
            "part([1,2,3],2)",
            "GIAC_ERROR: part([1,2,3],2) \n Error: Bad Argument Value",
            "2",
        ),
        ("part(x,1)", "x", "part(x,1)"),
        // The language has no `equal`: an unknown function here.
        ("equal2diff(equal(x,2))", "x-2", "equal(x,2)"),
        // Two heads: giac orders its operators by address, ours by name
        // (CAS-569).
        ("compare(sin(x),cos(x))", "1", "0"),
        ("compare(cos(x),sin(x))", "0", "1"),
        // giac's step of a downward quadratic is minus the box (CAS-570).
        (
            "linstep(Heaviside(1-x^2),x)",
            "-Heaviside(x+1)+Heaviside(x-1)",
            "Heaviside(x+1)-Heaviside(x-1)",
        ),
        // A period not known positive: giac errs, ours stays (CAS-573).
        (
            "periodic(x,x=a..b)",
            "GIAC_ERROR: periodic(x,x=(a .. b)) \n Error: Bad Argument Type",
            "periodic(x,x=a..b)",
        ),
        (
            "periodic(x,x=1..0)",
            "GIAC_ERROR: periodic(x,x=(1 .. 0)) \n Error: Bad Argument Type",
            "periodic(x,x=1..0)",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, answer, decided) in cases {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        if ours(input).unwrap() != ours(decided).unwrap() {
            wrong.push(format!("{input}: ours is not {decided}"));
        }
        if giac.agrees(input).is_ok() {
            wrong.push(format!("{input} agrees now"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
