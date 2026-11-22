//! Real roots against giac (CAS-538, CAS-539): `VAS`, `VAS_positive`,
//! `posubLMQ`, `poslbdLMQ`, `sturmseq`, `sturm` with bounds,
//! `rationalroot` and `crationalroot` read back and compared by value;
//! `complexroot` against closed forms within `10⁻¹⁰` (giac's `[a..b]`
//! intervals do not read back), with giac's count of roots; and the
//! departures, each with giac's answer.

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

/// `VAS` and `VAS_positive`, giac's examples and more.
const ISOLATION: &[&str] = &[
    "VAS(x^3-7*x+7)",
    "VAS(-x^3+7*x-7)",
    "VAS(x^4-2)",
    "VAS(x^2)",
    "VAS(x^3)",
    "VAS((x-1)^2*(x-2))",
    "VAS(x^2+1)",
    "VAS(x^5-3*x^2+x-1)",
    "VAS(x^4-10*x^2+1)",
    "VAS(1/2*x^2-1)",
    "VAS(x^2-2,x)",
    "VAS(x-100)",
    "VAS(x^2-300*x+20000)",
    "VAS(x-1/3)",
    "VAS((3*x-1)*(x-5))",
    "VAS(x^3-x)",
    "VAS((x^2-2)*(x-1)^2)",
    "VAS(x^2-1/100)",
    "VAS(x^5-x)",
    "VAS(x^4-5*x^2+4)",
    "VAS(x^3-3*x^2+2*x)",
    "VAS(x^6-14*x^4+49*x^2-36)",
    "VAS(2*x^3-3*x^2-11*x+6)",
    "VAS(5)",
    "VAS_positive(x^3-7*x+7)",
    "VAS_positive(x^2-3*x+2)",
    "VAS_positive(-x^2+3*x-2)",
    "VAS_positive(x^4-10*x^2+1)",
    "VAS_positive(x^3-3*x^2+2*x)",
    "VAS_positive(x^2-x)",
    "VAS_positive(x^2+1)",
];

/// `posubLMQ` and `poslbdLMQ`.
const BOUNDS: &[&str] = &[
    "posubLMQ(x^3-7*x+7)",
    "posubLMQ(x^3-7*x+7,x)",
    "posubLMQ(-x^3+7*x-7)",
    "posubLMQ(x^2-2)",
    "posubLMQ(x^2+1)",
    "posubLMQ(x^3+x+1)",
    "posubLMQ(7*x^3-7*x^2+1)",
    "posubLMQ(2*x^2-1)",
    "posubLMQ(x^2-4)",
    "posubLMQ(x^2-3)",
    "posubLMQ(x^2-5)",
    "posubLMQ(x-2)",
    "posubLMQ(x-5/2)",
    "posubLMQ(x^3-8)",
    "posubLMQ(x^2-100)",
    "posubLMQ(x^4-10*x^2+1)",
    "posubLMQ(x^2-3*x+2)",
    "posubLMQ(x^5-3*x^2+x-1)",
    "posubLMQ(x^5+3*x^4-2*x^3+x-7)",
    "posubLMQ(x^2)",
    "posubLMQ(x^2-x)",
    "posubLMQ(-x^2-1)",
    "posubLMQ(1/2*x-3)",
    "posubLMQ(5)",
    "poslbdLMQ(x^3-7*x+7)",
    "poslbdLMQ(x^2-2)",
    "poslbdLMQ(x^3+x+1)",
    "poslbdLMQ(x^4-10*x^2+1)",
    "poslbdLMQ(x-2)",
    "poslbdLMQ(x-3)",
    "poslbdLMQ(x-100)",
    "poslbdLMQ(x^2-3*x+2)",
    "poslbdLMQ(x^5-3*x^2+x-1)",
    "poslbdLMQ(x^5+3*x^4-2*x^3+x-7)",
    "poslbdLMQ(x^2)",
    "poslbdLMQ(x^2-x)",
    "poslbdLMQ(-x+5)",
    "poslbdLMQ(1/3*x^2-7*x+1/5)",
    "poslbdLMQ(5)",
];

/// `sturmseq` and `sturm` without bounds: signs, multiplicities, fractions.
const SEQUENCES: &[&str] = &[
    "sturmseq(x^3-7*x+7)",
    "sturmseq(x^3-7*x+7,x)",
    "sturm(x^3-7*x+7,x)",
    "sturm(y^3-7*y+7,y)",
    "sturmseq(-x^3+7*x-7)",
    "sturmseq(2*x^2-4)",
    "sturmseq(x^4+1)",
    "sturmseq(x^5-3*x^2+x-1)",
    "sturmseq(x^4-2*x^2+3)",
    "sturmseq(x^6+x^3-2*x+5)",
    "sturmseq(3*x^3-x+1/2)",
    "sturmseq((x^2-2)^2)",
    "sturmseq((x-1)^2*(x+2))",
    "sturmseq((x-1)^2*(x^2-2))",
    "sturmseq((x-1)^3*(x+2)^2*(x-3))",
    "sturmseq(-(x-1)^2*(x+2))",
    "sturmseq(-(x-1)*(x+2)^2)",
    "sturmseq(3*(x-1)^3*(2*x+1))",
    "sturmseq(-2*(x-1)*(x+2)^3*(x-3))",
    "sturmseq(-(x-1)*(x+2)^3)",
    "sturmseq(-(x-1)^3*(x+2))",
    "sturmseq(-(x-1)*(x+2)^2*(x-3)^3)",
    "sturmseq(-(x-1)^2*(x+2)^3)",
    "sturmseq(-(x-1)^5*(x+2)^3)",
    "sturmseq(-(x-1)*(x+2)^5*(x-3)^3)",
    "sturmseq(-(x-1)^2*(x+2)^4)",
    "sturmseq(-x^2*(x-1))",
    "sturmseq((x-1)/(x+2))",
    "sturmseq((x^2-1)/(2*x^2+4))",
    "sturmseq(-(x-1)/(x+2)^2)",
    "sturmseq((x-1)/(-x^2-2))",
    "sturmseq((x-1)^2/(3*x+2))",
    "sturmseq(5)",
    "sturmseq(-5)",
    "sturmseq(0)",
    "sturmseq(x)",
    "sturmseq(-3*x+6)",
];

/// `sturm` with bounds: roots and poles in `(a, b]`, or the sign at `a`.
const COUNTS: &[&str] = &[
    "sturm(x^3-7*x+7,x,-5,5)",
    "sturm(x^3-7*x+7,x,-5,1)",
    "sturm(x^3-7*x+7,x,1,5)",
    "sturm(x^3-7*x+7,x,1/2,3/2)",
    "sturm((x-1)^2*(x+2),x,-5,5)",
    "sturm((x-1)/(x+2),x,-5,5)",
    "sturm(x^2-2,x,0,2)",
    "sturm(x^2-1,x,1,2)",
    "sturm(x^2-1,x,-1,1)",
    "sturm(x^2-1,x,1,-1)",
    "sturm(x^2-1,x,1,1)",
    "sturm(x^2-1,x,2,2)",
    "sturm(x^2-1,x,0,0)",
    "sturm(-x^2-1,x,0,0)",
    "sturm(x^2,x,-1,1)",
    "sturm(x^2+1,x,-1,1)",
    "sturm(-x^2-1,x,-1,1)",
    "sturm(x^2-4,x,-1,1)",
    "sturm(5,x,0,1)",
    "sturm(-5,x,0,1)",
];

/// `rationalroot` and `crationalroot`.
const RATIONAL: &[&str] = &[
    "rationalroot(2*x^2-3*x+1)",
    "rationalroot((x-1)^2*(2*x+3))",
    "rationalroot((x-1)*(x+1))",
    "rationalroot(x^2-2)",
    "rationalroot(4*x^4-1)",
    "rationalroot(x^3)",
    "rationalroot(0)",
    "rationalroot(3)",
    "crationalroot(x^2+1)",
    "crationalroot((x^2+1)*(2*x-1))",
    "crationalroot((x^2+4)*(x^2+1))",
    "crationalroot((2*x-1)*(x+3)*(x^2+9))",
    "crationalroot(4*x^4-1)",
    "crationalroot(x^2+x+1)",
    "crationalroot(x^2+x+1/4)",
    "crationalroot(x^2-2*i*x-1)",
];

/// giac's exact answers, read back, are ours.
#[test]
fn exact_answers() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for input in [ISOLATION, BOUNDS, SEQUENCES, COUNTS, RATIONAL].concat() {
        match (giac.value(input), ours(input)) {
            (Ok(theirs), Ok(ours)) if theirs == ours => {}
            (theirs, ours) => wrong.push(format!("{input}: ours {ours:?}, giac's {theirs:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The number of entries of a list giac prints: the commas at depth 1.
fn entries(printed: &str) -> usize {
    let inner = printed.trim();
    let inner = inner
        .strip_prefix('[')
        .and_then(|s| s.strip_suffix(']'))
        .unwrap_or("");
    if inner.trim().is_empty() {
        return 0;
    }
    let mut depth = 0_i32;
    let mut count = 1;
    for c in inner.chars() {
        match c {
            '[' | '(' => depth += 1,
            ']' | ')' => depth -= 1,
            ',' if depth == 0 => count += 1,
            _ => {}
        }
    }
    count
}

/// Whether `a` and `b` agree to `10⁻¹⁰` in this CAS, real and imaginary
/// parts apart: `Abs` of a complex ball that holds 0 is indeterminate.
fn close(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let gap = Expr::call(Symbol::N, [a.clone() - b.clone(), Expr::from(30)]);
    [Symbol::RE, Symbol::IM].into_iter().all(|part| {
        let small = Expr::call(
            Symbol::LESS,
            [
                Expr::call(Symbol::ABS, [Expr::call(part, [gap.clone()])]),
                Expr::from(10).pow(-10),
            ],
        );
        ctx.eval(&small).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
    })
}

/// `complexroot`: each closed form, with its multiplicity, within `10⁻¹⁰`
/// of one of our roots, as many roots as giac lists.
#[test]
fn complex_roots() {
    let ctx = Context::new();
    let i3 = ["-1/2-sqrt(3)/2*i", "-1/2+sqrt(3)/2*i"];
    let cases: [(&str, &[(&str, i64)]); 7] = [
        (
            "complexroot(x^3-1,1/10^12)",
            &[(i3[0], 1), (i3[1], 1), ("1", 1)],
        ),
        ("complexroot(x^2+1,1/10^12)", &[("-i", 1), ("i", 1)]),
        (
            "complexroot((x-1)^2*(x^2+1),1/10^12)",
            &[("-i", 1), ("i", 1), ("1", 2)],
        ),
        (
            "complexroot(x^4-2,1/10^12)",
            &[
                ("-2^(1/4)", 1),
                ("-2^(1/4)*i", 1),
                ("2^(1/4)*i", 1),
                ("2^(1/4)", 1),
            ],
        ),
        (
            "complexroot(x^2+1,1/10^12,-2-2*i,2+2*i)",
            &[("-i", 1), ("i", 1)],
        ),
        ("complexroot(x^2+1,1/10^12,0,2+2*i)", &[("i", 1)]),
        ("complexroot((x-1)^2*(x+2),1/10^12,0,2+2*i)", &[("1", 2)]),
    ];
    let mut wrong = Vec::new();
    for (input, closed) in cases {
        let Ok(found) = ours(input) else {
            wrong.push(format!("{input}: unreadable"));
            continue;
        };
        let pairs: Vec<(Expr, Expr)> = found
            .args()
            .iter()
            .filter_map(|p| match p.args() {
                [z, m] => Some((z.clone(), m.clone())),
                _ => None,
            })
            .collect();
        let all_found = closed.iter().all(|(form, m)| {
            let form = ours(form).expect("a closed form");
            pairs
                .iter()
                .any(|(z, k)| *k == Expr::from(*m) && close(&ctx, z, &form))
        });
        if pairs.len() != closed.len() || !all_found {
            wrong.push(format!("{input}: ours {found:?}"));
        }
        if let Some(giac) = giac() {
            let printed = giac.eval(input).unwrap_or_default();
            if entries(&printed) != closed.len() {
                wrong.push(format!("{input}: giac prints {printed}"));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Where this CAS departs from giac on purpose: giac's answer as it
/// prints it, and this CAS's as decided.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        // giac's VAS_positive assumes a square-free polynomial and finds
        // a root in [0, 1] here; ours isolates the square-free part's.
        (
            "VAS_positive((x-1)^2*(x-2))",
            "[[0,1],1,[1,4]]",
            "[1,[1,4]]",
        ),
        // The sign of the primitive part goes on the factor of highest odd
        // multiplicity; giac, its square-free factorisation's way, keeps
        // it apart here when the highest multiplicity is even.
        (
            "sturmseq(-(x-1)^2*(x+2)^3*(x-5)^4)",
            "[1,[-1,0],[[1,2],1],1]",
            "[1,[[-1,-2],-1],1]",
        ),
        // The order of the rational roots in a multiplicity: giac's comes
        // from its search for linear factors; ours: the real ones
        // descending, then the others by real and imaginary part.
        (
            "rationalroot((x-2)*(x-3)*(x+1)*(2*x-1)*(3*x+1))",
            "[-1,1/2,-1/3,2,3]",
            "[3,2,1/2,-1/3,-1]",
        ),
        (
            "rationalroot((x+2)*(x-3)*(x+1)*(2*x+1)*(3*x-1)*(x-5))",
            "[3,5,-1/2,1/3,-2,-1]",
            "[5,3,1/3,-1/2,-1,-2]",
        ),
        (
            "rationalroot((x-1)^2*(x+1)*(x-3)^2*(x+5)^3)",
            "[-1,1,3,-5]",
            "[-1,3,1,-5]",
        ),
        (
            "crationalroot((x-1)*(x+1)*(x^2+4)*(x^2-2*x+2))",
            "[-1,2*i,-2*i,1-i,1,1+i]",
            "[1,-1,-2*i,2*i,1-i,1+i]",
        ),
    ];
    let mut wrong = Vec::new();
    for (input, answer, decided) in cases {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        if ours(input) != ours(decided) {
            wrong.push(format!("{input}: ours is not {decided}"));
        }
        if giac.agrees(input).is_ok() {
            wrong.push(format!("{input} agrees now"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
