//! Polynomial tools III against giac (CAS-540, CAS-543): `rootof`,
//! `fcoeff`, `homogeneize`, `total_degree`, `factor_xn`, `bezoutian` and
//! `is_irreducible`, giac's answers read back and compared with ours
//! (`Giac::agrees`: the same expression, or `==`); and the departures,
//! each with giac's answer as it prints it; and random inputs a run to
//! `fcoeff`, `total_degree`, `homogeneize` and `is_irreducible`.

use std::collections::BTreeMap;

use cas_giac::{Giac, ours};
use proptest::prelude::*;

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

/// `rootof` stays a call; polynomials become giac's lists (CAS-540).
#[test]
fn roots_of() {
    check(&[
        "rootof([1,0],[1,0,10,0,1])",
        "rootof([1,2],[1,0,-2])",
        "rootof(x+2,x^2-2)",
        "rootof(y+2,y^2-2)",
        "rootof(x^4+x+1)",
        "rootof([],[1,0,-2])",
        "rootof([3],[1,0,-2])",
        "rootof([1,0],[1,0,-2])+1",
        "rootof([1,0],[1,0,-2])^2",
    ]);
}

/// `fcoeff`, `homogeneize`, `total_degree` and `factor_xn` (CAS-540).
#[test]
fn forms() {
    check(&[
        "fcoeff([1,2,-1,1])",
        "fcoeff([1,2,0,1,3,-1])",
        "fcoeff([1,2,0,1,3,-1],y)",
        "fcoeff([])",
        "fcoeff([1/2,1])",
        "fcoeff([1,0])",
        "fcoeff([2,3])",
        "fcoeff([sqrt(2),1])",
        "homogeneize(x^2-y,z)",
        "homogeneize(x^2-1)",
        "homogeneize(x*y-1,y)",
        "homogeneize(x^3+x*y+1,z)",
        "homogeneize(x^2+y^3)",
        "homogeneize(x^2-1,z^2)",
        "homogeneize(5)",
        "homogeneize(x^2+x+1/2)",
        "homogeneize(1/x+1)",
        "total_degree(x^2*y+y,[x,y])",
        "total_degree(x^3*y+x*y,[x,y])",
        "total_degree(x^2*y*z+y,[x,y])",
        "total_degree(x^2*y+y)",
        "total_degree(x^2*y+y,x)",
        "total_degree(x^2*y+y,[z])",
        "total_degree(0,[x,y])",
        "total_degree(7,[x,y])",
        "factor_xn(x^3+x^2)",
        "factor_xn(x^4-1)",
        "factor_xn(x^4+12*x^3+54*x^2+108*x+81)",
        "factor_xn(2*x+1)",
        "factor_xn(x^2*y+y)",
        "factor_xn(y^2+y)",
        "factor_xn(5)",
        "factor_xn(0)",
    ]);
}

/// `bezoutian` (CAS-543), in `x`, the other symbols coefficients.
#[test]
fn bezoutians() {
    check(&[
        "bezoutian(5x^3-x,3x^2+1)",
        "bezoutian(x^2+1,x+2)",
        "bezoutian([1,0,1],[1,2])",
        "bezoutian(x^3+1,x^3+x)",
        "bezoutian(x+1,x+2)",
        "bezoutian(3,x+2)",
        "bezoutian(x^2+y,x+y)",
        "bezoutian(x^2+1,y+2)",
        "bezoutian(y^2+x,y+x)",
        "bezoutian(a*x^2+1,x+a)",
        "bezoutian(x^3-2*x+1,x^2-1)",
    ]);
}

/// `is_irreducible` (CAS-543), over ℚ.
#[test]
fn irreducibility() {
    check(&[
        "is_irreducible(x^2+1)",
        "is_irreducible(x^2-1)",
        "is_irreducible(x^2-2)",
        "is_irreducible(x^4+4)",
        "is_irreducible(4*x^2-4)",
        "is_irreducible(x^2)",
        "is_irreducible(x)",
        "is_irreducible(2*x+2)",
        "is_irreducible(x^2*y^2-1)",
        "is_irreducible(y^2-1)",
        "is_irreducible(x^2+y^2)",
        "is_irreducible(x^2*y+1)",
        "is_irreducible(x^2+1,x)",
        "is_irreducible(5)",
        "is_irreducible(0)",
    ]);
}

/// Where this CAS departs from giac on purpose: giac's answer as it
/// prints it (caseval's text for an error, cookbook C0007), and this
/// CAS's as decided.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let cases = [
        // No algebraic-number type (R1.2-003): rootof stays a call.
        (
            "normal(1/rootof([1,0],[1,0,10,0,1]))",
            "rootof([[-1,0,-10,0],[1,0,10,0,1]])",
            "1/rootof([1,0],[1,0,10,0,1])",
        ),
        (
            "evalf(rootof([1,0],[1,0,-2]))",
            "1.41421356237",
            "rootof([1,0],[1,0,-2])",
        ),
        // One list is the minimal polynomial; giac's caseval gives undef
        // (its command line infinity).
        ("rootof([1,0])", "undef", "rootof([1,0],[1,0])"),
        // giac errs on a symbolic minimal polynomial and on one not monic;
        // ours takes the first and leaves the second.
        (
            "rootof([1,0],x^2-2)",
            "GIAC_ERROR: rootof minimal polynomial must be unitary Error: Bad Argument Value",
            "rootof([1,0],[1,0,-2])",
        ),
        (
            "rootof([1,0],[2,0,-4])",
            "GIAC_ERROR: rootof minimal polynomial must be unitary Error: Bad Argument Value",
            "rootof([1,0],[2,0,-4])",
        ),
        // No polynomial in the variables: the call stays.
        ("factor_xn(1/x+1)", "1+x^-1", "factor_xn(1/x+1)"),
        (
            "total_degree(1/x+y,[x,y])",
            "1",
            "total_degree(1/x+y,[x,y])",
        ),
        (
            "is_irreducible(1/x)",
            "GIAC_ERROR: is_irreducible(1/x) \n Error: Bad Argument Type",
            "is_irreducible(1/x)",
        ),
        // giac's sequence of the answer and the extra argument; a list.
        ("factor_xn(x^2+x,y)", "x^2*(1+x^-1),y", "factor_xn(x^2+x,y)"),
        (
            "homogeneize([x^2-1,x],z)",
            "[x^2-1,x]",
            "homogeneize([x^2-1,x],z)",
        ),
        // giac rejects a variable; ours leaves the call.
        (
            "bezoutian(x^2+1,x+2,x)",
            "GIAC_ERROR: bezoutian(x^2+1,x+2,x) \n Error: Bad Argument Value",
            "bezoutian(x^2+1,x+2,x)",
        ),
        // x y has two factors: giac's 1 is wrong.
        ("is_irreducible(x*y)", "1", "0"),
        ("is_irreducible(x*y+x)", "1", "0"),
        // 2*(y+1)*(x^2*y-x^2-1), the sweep's shrunk input (2.1 T7): giac
        // misses a factor in y alone.
        ("is_irreducible((-2)+(-2)*y+(-2)*x^2+2*x^2*y^2)", "1", "0"),
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
    // giac's `mod` suffix: the language does not read it.
    let modular = "is_irreducible(x^3+x+1 mod 7)";
    if giac.eval(modular).unwrap() != "1" {
        wrong.push(format!("{modular}: giac's answer changed"));
    }
    if ours(modular).is_ok() {
        wrong.push(format!("{modular}: the language reads it now"));
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The variables of the random polynomials, in order.
const VARIABLES: [&str; 3] = ["x", "y", "z"];

/// The terms of a polynomial in the first `n` of [`VARIABLES`]: exponent
/// vector to coefficient.
type Terms = BTreeMap<Vec<u32>, i64>;

/// Up to `len` terms with exponents below 3 and coefficients in `-5..=5`,
/// some of them 0.
fn terms(n: usize, len: usize) -> impl Strategy<Value = Terms> {
    prop::collection::btree_map(prop::collection::vec(0u32..3, n), -5i64..=5, 0..=len)
}

/// As [`terms`], at least one and none 0, so a nonzero polynomial.
fn nonzero(n: usize, len: usize) -> impl Strategy<Value = Terms> {
    prop::collection::btree_map(
        prop::collection::vec(0u32..3, n),
        prop_oneof![-5i64..=-1, 1i64..=5],
        1..=len,
    )
}

/// The polynomial as Xcas text, in parentheses.
fn text(terms: &Terms) -> String {
    let parts: Vec<String> = terms
        .iter()
        .filter(|&(_, &c)| c != 0)
        .map(|(exps, c)| {
            let mut factors = vec![if *c < 0 {
                format!("({c})")
            } else {
                c.to_string()
            }];
            for (v, &k) in VARIABLES.iter().zip(exps) {
                match k {
                    0 => {}
                    1 => factors.push((*v).to_owned()),
                    _ => factors.push(format!("{v}^{k}")),
                }
            }
            factors.join("*")
        })
        .collect();
    if parts.is_empty() {
        "(0)".to_owned()
    } else {
        format!("({})", parts.join("+"))
    }
}

/// Whether the polynomial in `x`, `y` has `x` and, read as a polynomial in
/// `x`, a content of positive degree in `y`: the gcd over ℚ[y] of its
/// coefficients, by primitive pseudo-remainders on `i128`.
fn y_content(p: &Terms) -> bool {
    let mut rows: BTreeMap<u32, Vec<i128>> = BTreeMap::new();
    for (e, &c) in p.iter().filter(|&(_, &c)| c != 0) {
        rows.entry(e[0]).or_insert_with(|| vec![0; 3])[e[1] as usize] += i128::from(c);
    }
    if rows.keys().all(|&k| k == 0) {
        return false;
    }
    let mut rows = rows.into_values();
    let first = rows.next().unwrap_or_default();
    let g = rows.fold(first, gcd_y);
    g.iter().rposition(|&c| c != 0).is_some_and(|d| d > 0)
}

/// `a` with trailing zeros dropped and its integer content divided out.
fn primitive(mut a: Vec<i128>) -> Vec<i128> {
    while a.last() == Some(&0) {
        a.pop();
    }
    let g = a.iter().fold(0i128, |g, &c| {
        let (mut g, mut c) = (g, c.abs());
        while c != 0 {
            (g, c) = (c, g % c);
        }
        g
    });
    if g > 1 {
        for c in &mut a {
            *c /= g;
        }
    }
    a
}

/// The gcd in ℚ[y] of `a` and `b`, up to a constant (empty for 0).
fn gcd_y(a: Vec<i128>, b: Vec<i128>) -> Vec<i128> {
    let (mut a, mut b) = (primitive(a), primitive(b));
    while !b.is_empty() {
        while a.len() >= b.len() && !a.is_empty() {
            let (la, lb) = (a[a.len() - 1], b[b.len() - 1]);
            let shift = a.len() - b.len();
            for c in &mut a {
                *c *= lb;
            }
            for (i, &c) in b.iter().enumerate() {
                a[i + shift] -= la * c;
            }
            a = primitive(a);
        }
        (a, b) = (b, a);
    }
    a
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// `fcoeff` of up to four distinct roots `k/2`, `k` in `-6..=6`, each
    /// with a multiplicity in `-2..=3`, perhaps in `y`.
    #[test]
    fn random_fcoeff(
        roots in prop::collection::btree_map(-6i64..=6, -2i64..=3, 0..=4),
        in_y in any::<bool>(),
    ) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let list: Vec<String> = roots
            .iter()
            .map(|(k, m)| if k % 2 == 0 { format!("{},{m}", k / 2) } else { format!("{k}/2,{m}") })
            .collect();
        let input = format!("fcoeff([{}]{})", list.join(","), if in_y { ",y" } else { "" });
        giac.agrees(&input).map_err(|e| TestCaseError::fail(e.to_string()))?;
    }

    /// `total_degree` of a polynomial in `x`, `y`, `z` in some of them.
    #[test]
    fn random_total_degree(p in terms(3, 4), vars in 0usize..4) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let vars = ["[x,y]", "[x]", "[x,y,z]", "[y,z]"][vars];
        let input = format!("total_degree({},{vars})", text(&p));
        giac.agrees(&input).map_err(|e| TestCaseError::fail(e.to_string()))?;
    }

    /// `homogeneize` of a nonzero polynomial in `x`, `y`, with `z` or
    /// giac's default `t`.
    #[test]
    fn random_homogeneize(p in nonzero(2, 4), named in any::<bool>()) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let input = format!("homogeneize({}{})", text(&p), if named { ",z" } else { "" });
        giac.agrees(&input).map_err(|e| TestCaseError::fail(e.to_string()))?;
    }

    /// `is_irreducible` of a nonzero polynomial in `x`, `y`, over ℚ.
    #[test]
    fn random_is_irreducible(p in nonzero(2, 4)) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        // A factor in y alone, as in x*y, x*y+x and (y+1)*(x^2*y-x^2-1):
        // giac's 1 there is wrong (README departure).
        if y_content(&p) {
            return Ok(());
        }
        let input = format!("is_irreducible({})", text(&p));
        giac.agrees(&input).map_err(|e| TestCaseError::fail(e.to_string()))?;
    }
}
