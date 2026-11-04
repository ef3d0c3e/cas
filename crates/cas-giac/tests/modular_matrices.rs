//! Modular arithmetic II against giac (CAS-DOC-44, CAS-PLAN-38): giac's
//! linear algebra of matrices of `%` objects, read back, is this CAS's;
//! its `ker` and `solve` answers, partly plain integers (`solve((x^3-1) %
//! 13, x)` is `[-4,1,3]`), are compared modulo `p`, and its `eigenvals`,
//! in its own order, as sorted lists. `ichinrem` and `powmod` read back.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// giac's answers to `inputs`, read back, are this CAS's.
fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| giac.agrees(input).err().map(|e| e.to_string()))
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's answers and this CAS's, both taken modulo `p`, entry by entry,
/// and sorted when `sorted`: to each pair of inputs, giac's and this CAS's.
fn check_modulo(p: i64, inputs: &[(&str, &str)], sorted: bool) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let modulo = |e: Expr| {
        let mut found = ctx
            .eval(&Expr::call(Symbol::MODULAR, [e, Expr::from(p)]))
            .unwrap();
        if sorted {
            let mut items = found.args().to_vec();
            items.sort_by_key(|x| format!("{x:?}"));
            found = Expr::call(Symbol::LIST, items);
        }
        found
    };
    let mut wrong = Vec::new();
    for (theirs, mine) in inputs {
        match (giac.value(theirs), ours(mine)) {
            (Ok(t), Ok(o)) if modulo(t.clone()) == modulo(o.clone()) => {}
            (t, o) => wrong.push(format!("{mine}: ours {o:?}, giac's {t:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Determinants, inverses and characteristic polynomials (CAS-457).
#[test]
fn determinants_and_inverses() {
    check(&[
        "det([[1,2],[3,4]] % 7)",
        "det([[2,3,1],[4,1,5],[6,2,2]] % 101)",
        "det([[1,2],[3,4]] % 6)",
        "det([[3,4],[5,6]] % 12)",
        "inv([[1,2],[3,4]] % 7)",
        "inv([[2,3,1],[4,1,5],[6,2,2]] % 11)",
        "charpoly([[1,2],[3,4]] % 7, x)",
        "charpoly([[2,3,1],[4,1,5],[6,2,2]] % 7, x)",
        "charpoly([[3,4],[5,6]] % 12, x)",
    ]);
}

/// Reduced forms and linear systems over 𝔽ₚ (CAS-458).
#[test]
fn reduced_forms() {
    check(&[
        "rref([[1,2,3],[4,5,6],[7,8,10]] % 5)",
        "rref([[1,2],[2,4]] % 7)",
        "rref([[1,1,0,1],[0,1,1,1],[1,0,1,0]] % 2)",
        "rank([[1,2,3],[4,5,6],[7,8,9]] % 3)",
        "rank([[1,1,0,1],[0,1,1,1],[1,0,1,0]] % 2)",
        "image([[1,2,3],[2,4,7]] % 5)",
        "linsolve([x+2*y=1,3*x+4*y=2] % 7,[x,y])",
        "linsolve([x+2*y=1 % 7,3*x+4*y=2],[x,y])",
        "linsolve([x+y=1,x+y=2] % 5,[x,y])",
        "linsolve([x+y=1,2*x+2*y=2] % 5,[x,y])",
        "linsolve([x+y+z=1,x-y=2,y+2*z=0] % 13,[x,y,z])",
    ]);
    // giac writes the −1 of a kernel vector plainly, or modulo 2 as 1 % 2.
    check_modulo(
        5,
        &[("ker([[1,2,3],[2,4,6]] % 5)", "ker([[1,2,3],[2,4,6]] % 5)")],
        false,
    );
    check_modulo(
        7,
        &[("ker([[1,2],[2,4]] % 7)", "ker([[1,2],[2,4]] % 7)")],
        false,
    );
    check_modulo(
        2,
        &[(
            "ker([[1,1,0,1],[0,1,1,1],[1,0,1,0]] % 2)",
            "ker([[1,1,0,1],[0,1,1,1],[1,0,1,0]] % 2)",
        )],
        false,
    );
}

/// Roots modulo a prime and eigenvalues over 𝔽ₚ (CAS-459).
#[test]
fn roots_and_eigenvalues() {
    check_modulo(
        7,
        &[
            ("solve((x^2-2) % 7, x)", "solve((x^2-2) % 7, x)"),
            ("solve(x^2 % 7 = 2 % 7, x)", "solve(x^2 % 7 = 2 % 7, x)"),
        ],
        false,
    );
    check_modulo(
        13,
        &[("solve((x^3-1) % 13, x)", "solve((x^3-1) % 13, x)")],
        false,
    );
    check_modulo(
        5,
        &[("solve((x^5-x) % 5, x)", "solve((x^5-x) % 5, x)")],
        false,
    );
    check(&["solve((x^2+1) % 3, x)"]);
    check_modulo(
        7,
        &[(
            "[eigenvals([[2,1],[1,2]] % 7)]",
            "eigenvals([[2,1],[1,2]] % 7)",
        )],
        true,
    );
    check_modulo(
        11,
        &[(
            "[eigenvals([[2,0,0],[1,3,0],[4,5,6]] % 11)]",
            "eigenvals([[2,0,0],[1,3,0],[4,5,6]] % 11)",
        )],
        true,
    );
    check_modulo(
        5,
        &[(
            "[eigenvals([[1,1],[0,1]] % 5)]",
            "eigenvals([[1,1],[0,1]] % 5)",
        )],
        true,
    );
}

/// Chinese remainders and modular powers of polynomials (CAS-460).
#[test]
fn remainders_and_powers() {
    check(&[
        "ichinrem(2 % 3, 3 % 5)",
        "ichinrem([2 % 3, 3 % 5, 1 % 7])",
        "ichinrem(1 % 4, 3 % 6)",
        "ichinrem([3 % 7, 4 % 9, 5 % 11])",
        "powmod(x+1,10,3,x^2+1,x)",
        "powmod(x,100,2,x^3+x+1,x)",
        "powmod(3*x+1,7,7,x^2+3,x)",
        "powmod(x^3+2*x,50,5,x^4+x+2,x)",
    ]);
}
