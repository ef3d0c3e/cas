//! Xcas's `kovacicsols` against giac (CAS-584): giac's answers read back
//! and compared with ours (`Giac::agrees`); case 2's radicals, where the
//! forms differ, by the ratio of each solution to giac's, a constant; and
//! the departures, each with giac's answer as libgiac prints it. giac quotes
//! its arguments, so its equations in `y''` do not read here and
//! `diff(y(x),x,x)` does not read there: the rows give the coefficients
//! `[a,b,c]`, written without `/` (giac answers 0 to a list holding one).

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

/// giac's examples (its aide and the probes of CAS-584) and the equations
/// with no Liouvillian solution: cases 1 and 2, `[]`, and giac's `0`.
#[test]
fn agreements() {
    let Some(giac) = giac() else {
        return;
    };
    let inputs = [
        // Case 2, x^(1/4) e^(±2√x).
        "kovacicsols([16*x^2,0,-(16*x-3)])",
        // Case 1.
        "kovacicsols([36*(x-1)^2,0,-27])",
        "kovacicsols([x^2,x^2,-6])",
        "kovacicsols([1,0,19-x^2])",
        "kovacicsols([(1+x^2)^2,0,3])",
        "kovacicsols([1,0,3-x^2])",
        "kovacicsols([1,0,3-x^2],x)",
        "kovacicsols([1,0,-(x^2-21)])",
        // giac's x^3/x and x^5/x^2.
        "kovacicsols([t^2,0,-2],t)",
        "kovacicsols([x^2,0,-2])",
        "kovacicsols([x^2,0,-6])",
        // No solution.
        "kovacicsols([1,0,-x])",
        "kovacicsols([1,0,-(x^3+1)])",
        "kovacicsols([1,0,-1])",
        "kovacicsols([x^4,0,-(x^2+1)])",
        // giac's return false.
        "kovacicsols([1,0,1])",
        "kovacicsols([1,0,0])",
        "kovacicsols([t^2,t,t^2-1/4],t)",
    ];
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| giac.agrees(input).err().map(|e| format!("{e}")))
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Case 2 with radicals: giac's fourth roots of `−t (t − 1) (2√(t² − t) +
/// 2t − 1)` and ours of `(t² − t)` and `(4t − 2 + 4√(t² − t))` differ by
/// a constant factor, the same at `t = 3` and `t = 5`.
#[test]
fn radicals_up_to_constants() {
    let Some(giac) = giac() else {
        return;
    };
    let input = "kovacicsols([16*(t-1)^2*t^2,0,3*(t^2-t+1)],t)";
    let mine = ours(input).unwrap();
    let theirs = giac.value(input).unwrap();
    assert_eq!(mine.args().len(), 2, "{mine:?}");
    assert_eq!(theirs.args().len(), 2, "{theirs:?}");
    let ctx = Context::new();
    let t = Expr::from(Symbol::new("t"));
    let at = |ratio: &Expr, point: i64| {
        let rule = Expr::call(Symbol::RULE, [t.clone(), Expr::from(point)]);
        let value = Expr::call(Symbol::REPLACE_ALL, [ratio.clone(), rule]);
        ctx.eval(&Expr::call(Symbol::N, [value, Expr::from(30)]))
            .unwrap()
    };
    for (o, g) in mine.args().iter().zip(theirs.args()) {
        let ratio = Expr::call(Symbol::TIMES, [o.clone(), g.clone().pow(-1)]);
        let (r3, r5) = (at(&ratio, 3), at(&ratio, 5));
        let gap = Expr::call(
            Symbol::ABS,
            [Expr::call(
                Symbol::PLUS,
                [r3.clone(), Expr::from(-1) * r5.clone()],
            )],
        );
        let small = Expr::call(
            Symbol::LESS,
            [
                gap,
                Expr::call(Symbol::POWER, [Expr::from(10), Expr::from(-20)]),
            ],
        );
        assert_eq!(
            ctx.eval(&small).unwrap(),
            Expr::from(Symbol::TRUE),
            "{o:?} over {g:?}: {r3:?} at 3, {r5:?} at 5"
        );
    }
}

/// Where giac's answer is not ours: giac's text, and ours as decided.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        // Case 3 with a y′ term: giac reads the shifted polynomial's
        // coefficients highest first as lowest first (desolve.cc:2417-2422),
        // so its polynomial is the reciprocal one, whose roots are not
        // y′/y (unit test case_three_riccati).
        (
            "kovacicsols([48*t*(t+1)*(5*t-4),8*(25*t+16)*(t-2),-(5*t+68)],t)",
            "omega_^4*(135*t^4-616*t^3-144*t^2+3072*t-4096)/20736-omega_^3*t*(t+1)*\
             (23*t^2-92*t+128)/54-omega_^2*t^2*(t+1)*(15*t^3-80*t^2+80*t+256)/24+\
             omega_*2*t^3*(t-4)*(t+1)^2*(5*t+8)/3-t^4*(t+1)^2*(t+4)*(5*t+4)",
            "-omega_^4*t^4*(1 + t)^2*(4 + t)*(4 + 5*t) + 2*omega_^3*t^3*(-4 + t)*\
             (1 + t)^2*(8 + 5*t)/3 - omega_*t*(1 + t)*(128 - 92*t + 23*t^2)/54 - \
             omega_^2*t^2*(1 + t)*(256 + 80*t - 80*t^2 + 15*t^3)/24 + \
             (-4096 + 3072*t - 144*t^2 - 616*t^3 + 135*t^4)/20736",
        ),
        // Poles ±√2: giac finds none, ours case 2's two solutions.
        (
            "kovacicsols([(x^2-2)^2,0,-1])",
            "[]",
            "[exp(6^(1/2)*ln(-2^(1/2)/(2^(1/2) + x) + x/(2^(1/2) + x))/4)*(-2 + x^2)^(1/2), \
             exp(-6^(1/2)*ln(-2^(1/2)/(2^(1/2) + x) + x/(2^(1/2) + x))/4)*(-2 + x^2)^(1/2)]",
        ),
        // A pole of x³ − 2, past this port: the call stays.
        (
            "kovacicsols([x^3-2,0,-1])",
            "[]",
            "kovacicsols([x^3-2,0,-1])",
        ),
        (
            "kovacicsols([x^3-2,3*x^2,0])",
            "[1]",
            "kovacicsols([x^3-2,3*x^2,0])",
        ),
        // giac quotes its arguments: a `/` in the list reads as a name.
        ("kovacicsols([1,0,-2/x^2])", "0", "[x^2]"),
        // An equation in diff(y(x),x,x) is giac's error, ours solves it.
        (
            "kovacicsols(diff(y(x),x,x)=(x^2-3)*y(x))",
            "GIAC_ERROR: kovacicsols(diff(y(x),x,x)=((x^2-3)*y(x))) \n Error: Bad Argument Value",
            "[x*exp(-x^2/2)]",
        ),
        (
            "kovacicsols(diff(y(x),x,x)=(x^2-3)*y(x),x,y)",
            "GIAC_ERROR: kovacicsols(diff(y(x),x,x)=((x^2-3)*y(x)),x,y) \n Error: Bad Argument Value",
            "[x*exp(-x^2/2)]",
        ),
        (
            "kovacicsols(diff(y(x),x,x)=(x^2-3)*y(x),y(x))",
            "GIAC_ERROR: kovacicsols(diff(y(x),x,x)=((x^2-3)*y(x)),y(x)) \n Error: Bad Argument Value",
            "[x*exp(-x^2/2)]",
        ),
        // giac's errors: the call stays, with the warning.
        (
            "kovacicsols(diff(y(x),x)+y(x)=0)",
            "GIAC_ERROR: kovacicsols(diff(y(x),x)+y(x)=0) \n Error: Bad Argument Value",
            "kovacicsols(diff(y(x),x)+y(x)=0)",
        ),
        (
            "kovacicsols([0,1,1])",
            "GIAC_ERROR: kovacicsols([0,1,1]) \n Error: Bad Argument Value",
            "kovacicsols([0,1,1])",
        ),
        (
            "kovacicsols([x^2,0,-2,1])",
            "GIAC_ERROR: kovacicsols([x^2,0,-2,1]) \n Error: Bad Argument Value",
            "kovacicsols([x^2,0,-2,1])",
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
