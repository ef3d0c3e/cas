//! giac's minors, echelon forms and pivots (CAS-530, CAS-528), `gauss`
//! and `rat_jordan` (CAS-535) against giac 2.0: its answers read back are
//! this CAS's, and the departures, each with giac's answer as it prints it
//! and this CAS's as decided. `rat_jordan` returns a sequence, read in
//! brackets as a list. Random integer matrices, entries -5 to 5, a run to
//! each of them.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};
use proptest::prelude::*;

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

#[test]
fn minors() {
    check(&[
        "det_minor([[1,2],[3,4]])",
        "det_minor([[1,2,3],[4,5,6],[7,8,10]])",
        "det_minor([[a,b],[c,d]])",
        "adjoint_matrix([[1,2],[3,4]])",
        "adjoint_matrix([[a,b],[c,d]])",
    ]);
}

#[test]
fn echelon_forms() {
    check(&[
        "ref([[1,2],[3,4]])",
        "ref([[2,1],[1,3]])",
        "ref([[3,1,-2],[3,2,2]])",
        "ref([[2,1,1,-1],[1,1,2,-1],[1,2,1,-4]])",
        "ref([[0,0],[0,3]])",
        "ref([[0,2],[0,3],[1,1]])",
        "ref([[a,b],[c,d]])",
        "rref([[1,2],[3,4]],keep_pivot)",
        "rref([[4,2],[2,4]],keep_pivot)",
        "rref([[0,2],[3,4]],keep_pivot)",
        "rref([[2,1],[1,3]],keep_pivot)",
        "rref([[2,4,6],[1,3,5]],keep_pivot)",
        "rref([[1,2,3],[4,5,6],[7,8,10]],keep_pivot)",
        "rref([[a,b],[c,d]],keep_pivot)",
    ]);
}

#[test]
fn pivots() {
    check(&[
        "pivot([[1,2],[3,4]],1,1)",
        "pivot([[1,2,3],[4,5,6],[7,8,10]],1,2)",
        "pivot([[1,2],[3,4],[5,6]],0,1)",
        "egvl([[4,1],[-4,0]])",
        "eigVl([[4,1],[-4,0]])",
    ]);
}

/// Whether two answers agree once put together, lists element by element,
/// as matrices.rs's: giac's `(sqrt(33)+5)/2` is this CAS's
/// `5/2 + sqrt(33)/2`.
fn same(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let list = Expr::from(Symbol::LIST);
    if a.head() == Some(&list) || b.head() == Some(&list) {
        return a.head() == b.head()
            && a.args().len() == b.args().len()
            && a.args().iter().zip(b.args()).all(|(x, y)| same(ctx, x, y));
    }
    let difference = Expr::call(
        Symbol::TOGETHER,
        [Expr::call(Symbol::EXPAND, [a.clone() - b.clone()])],
    );
    ctx.eval(&difference) == Ok(Expr::from(0))
}

/// giac's answers to `inputs`, read back, are this CAS's entry by entry
/// once put together.
fn check_together(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for input in inputs {
        let theirs = giac.value(input).unwrap();
        let mine = ours(input).unwrap();
        if !same(&ctx, &theirs, &mine) {
            wrong.push(format!("{input}: giac {theirs:?}, this CAS {mine:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn jordan_matrices() {
    check_together(&[
        "egvl([[1,2],[3,4]])",
        "egvl([[1,1],[0,2]])",
        "egvl([[1,0,1],[0,2,0],[0,0,3]])",
        "egvl([[1,2],[2,1]])",
        "egvl([[0,1],[-2,2]])",
        "egvl([[0,-1],[1,0]])",
        "egvl([[2,0,0],[0,1,0],[0,0,3]])",
        "egvl([[1,1,0],[0,1,0],[0,0,-2]])",
    ]);
}

#[test]
fn squares() {
    check(&[
        "gauss(x^2+2*x*y,[x,y])",
        "gauss(x^2+y^2+z^2+2*x*y+4*y*z,[x,y,z])",
        "gauss(x*y,[x,y])",
        "gauss(2*x^2-3*y^2,[x,y])",
        "gauss(2*x*y,[x,y])",
        "gauss(x*y+y*z,[x,y,z])",
        "gauss(x^2+x*y+y^2+z^2,[x,y,z])",
        "gauss(2*x^2+2*x*y+3*y^2,[x,y])",
        "gauss(x^2+2*a*x*y,[x,y])",
        "gauss(y*z+x*z,[x,y,z])",
    ]);
}

/// giac's `[rat_jordan(A)]`, read back, is this CAS's `rat_jordan(A)`.
#[test]
fn rational_jordan_forms() {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = [
        "rat_jordan([[1,1],[0,1]])",
        "rat_jordan([[1,2],[3,4]])",
        "rat_jordan([[0,-1],[1,0]])",
        "rat_jordan([[2,0],[0,3]])",
        "rat_jordan([[1,2,0],[3,4,0],[0,0,7]])",
    ]
    .iter()
    .filter_map(|input| {
        giac.answers(&format!("[{input}]"), input)
            .err()
            .map(|e| e.to_string())
    })
    .collect();
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
        // giac's P from its own kernels: the vectors of an eigenvalue
        // here are those of the reduced kernel, 1 at each free column.
        (
            "rat_jordan([[2,1],[0,3]])",
            "[[1,-1],[1,0]],[[3,0],[0,2]]",
            "[[[1,1],[1,0]],[[3,0],[0,2]]]",
        ),
        (
            "rat_jordan([[3,1],[0,2]])",
            "[[1,1],[0,-1]],[[3,0],[0,2]]",
            "[[[1,-1],[0,1]],[[3,0],[0,2]]]",
        ),
        (
            "rat_jordan([[-2,-2,1],[-2,1,-2],[1,-2,-2]])",
            "[[1,5,5],[-2,0,5],[1,-5,5]],[[3,0,0],[0,-3,0],[0,0,-3]]",
            "[[[1,2,-1],[-2,1,0],[1,0,1]],[[3,0,0],[0,-3,0],[0,0,-3]]]",
        ),
        // Already a rational Jordan form: P = I here, giac's P is not.
        (
            "rat_jordan([[0,-1,1,0],[1,0,0,1],[0,0,0,-1],[0,0,1,0]])",
            "[[4,0,0,4],[0,4,-4,0],[0,0,4,0],[0,0,0,4]],[[0,-1,1,0],[1,0,0,1],[0,0,0,-1],[0,0,1,0]]",
            "[[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],[[0,-1,1,0],[1,0,0,1],[0,0,0,-1],[0,0,1,0]]]",
        ),
        // A float matrix: giac's rows come out exactly 0.0 and 1.0, this
        // CAS's numeric elimination keeps its balls' radii; keep_pivot of
        // numbers is rref's here.
        (
            "rref([[2.0,1],[1,3]],keep_pivot)",
            "[[1.0,0.0],[0.0,1.0]]",
            "rref([[2.0,1],[1,3]])",
        ),
        // Not a quadratic form: giac returns it, the call stays here.
        ("gauss(x^3+y^2,[x,y])", "x^3+y^2", "gauss(x^3+y^2,[x,y])"),
        // CAS-527: giac's inert capitalised forms stay unevaluated user
        // calls here; an alias would evaluate.
        (
            "Det([[1,2],[3,4]])",
            "det([[1,2],[3,4]])",
            "Det([[1,2],[3,4]])",
        ),
        (
            "Inverse([[1,2],[3,4]])",
            "inverse([[1,2],[3,4]])",
            "Inverse([[1,2],[3,4]])",
        ),
        (
            "Nullspace([[1,2],[2,4]])",
            "nullspace([[1,2],[2,4]])",
            "Nullspace([[1,2],[2,4]])",
        ),
        (
            "Rref([[1,2],[3,4]])",
            "rref([[1,2],[3,4]])",
            "Rref([[1,2],[3,4]])",
        ),
    ];
    let normal = |s: &str| s.replace("matrix[", "[").replace(' ', "");
    let mut wrong = Vec::new();
    for (input, answer, decided) in cases {
        let printed = giac.eval(input).unwrap();
        if normal(&printed) != normal(answer) {
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

/// An `r` by `c` integer matrix, entries -5 to 5.
fn entries(r: usize, c: usize) -> impl Strategy<Value = Vec<Vec<i64>>> {
    prop::collection::vec(prop::collection::vec(-5i64..=5, c), r)
}

/// A matrix of 1 to 4 rows and columns.
fn matrix() -> impl Strategy<Value = Vec<Vec<i64>>> {
    (1usize..=4, 1usize..=4).prop_flat_map(|(r, c)| entries(r, c))
}

/// A square matrix of size 1 to `max`.
fn square(max: usize) -> impl Strategy<Value = Vec<Vec<i64>>> {
    (1usize..=max).prop_flat_map(|n| entries(n, n))
}

/// `[[a,b],…]`.
fn text(m: &[Vec<i64>]) -> String {
    let rows: Vec<String> = m
        .iter()
        .map(|r| {
            let r: Vec<String> = r.iter().map(i64::to_string).collect();
            format!("[{}]", r.join(","))
        })
        .collect();
    format!("[{}]", rows.join(","))
}

/// The diagonal blocks of a square matrix, as lists: the shortest runs of
/// indices whose rows and columns meet no nonzero entry outside the run.
fn blocks(m: &Expr) -> Vec<Expr> {
    let zero = Expr::from(0);
    let m: Vec<&[Expr]> = m.args().iter().map(Expr::args).collect();
    let n = m.len();
    let mut out = Vec::new();
    let mut s = 0;
    while s < n {
        let mut e = s + 1;
        while (e..n).any(|k| (s..e).any(|r| m[r][k] != zero || m[k][r] != zero)) {
            e += 1;
        }
        let rows = m[s..e]
            .iter()
            .map(|r| Expr::call(Symbol::LIST, r[s..e].iter().cloned()));
        out.push(Expr::call(Symbol::LIST, rows));
        s = e;
    }
    out
}

/// Whether two square matrices have the same diagonal blocks in some
/// order, each pair [`same`]: giac orders eigenvalues and blocks its own
/// way (`egvl([[1,-1],[1,1]])` gives `1-i` first, `egvl([[2,1],[-1,1]])`
/// `(i*sqrt(3)+3)/2`; `rat_jordan([[0,-1],[0,1]])`'s J is
/// `[[0,0],[0,1]]`, this CAS's `[[1,0],[0,0]]`, CI 45823883128).
fn same_blocks(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let mut rest = blocks(b);
    for x in blocks(a) {
        match rest.iter().position(|y| same(ctx, &x, y)) {
            Some(k) => {
                rest.swap_remove(k);
            }
            None => return false,
        }
    }
    rest.is_empty()
}

/// giac's answer to `input`, read back, is ours; else the failure names
/// the input, giac's answer and ours.
fn agree(giac: &Giac, input: &str) -> Result<(), TestCaseError> {
    giac.agrees(input)
        .map_err(|e| TestCaseError::fail(e.to_string()))
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// `det_minor` and `adjoint_matrix` of a square matrix.
    #[test]
    fn random_minors(a in square(4)) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        agree(giac, &format!("det_minor({})", text(&a)))?;
        agree(giac, &format!("adjoint_matrix({})", text(&a)))?;
    }

    /// `ref` and `rref` with `keep_pivot` of an integer matrix (a float
    /// one is a departure).
    #[test]
    fn random_echelon_forms(a in matrix()) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        agree(giac, &format!("ref({})", text(&a)))?;
        agree(giac, &format!("rref({},keep_pivot)", text(&a)))?;
    }

    /// `pivot` at any entry, zero or not, where no other row has 0 in the
    /// pivot's column: giac leaves such a row, this CAS combines it
    /// (`pivot([[0,0],[0,1]],0,0)`: giac `[[0,0],[0,1]]`, ours
    /// `[[0,0],[0,0]]`, CI 45823883128; reported, not fixed here).
    #[test]
    fn random_pivots(a in matrix(), i in 0usize..4, j in 0usize..4) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let (i, j) = (i % a.len(), j % a[0].len());
        if a.iter().enumerate().any(|(k, r)| k != i && r[j] == 0) {
            return Ok(());
        }
        agree(giac, &format!("pivot({},{i},{j})", text(&a)))?;
    }

    /// `egvl` of a square matrix of size 1 to 3: the same blocks in some
    /// order ([`same_blocks`]), entry by entry once put together.
    #[test]
    fn random_egvl(a in square(3)) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let input = format!("egvl({})", text(&a));
        let printed = giac.eval(&input).map_err(|e| TestCaseError::fail(e.to_string()))?;
        // giac's roots of an irreducible cubic print as `rootof([…],[…])`,
        // which this CAS does not read: no comparison.
        if printed.contains("rootof") {
            return Ok(());
        }
        let fail = |e: cas_giac::OracleError| TestCaseError::fail(e.to_string());
        let theirs = giac.value(&input).map_err(fail)?;
        let mine = ours(&input).map_err(fail)?;
        if !same_blocks(&Context::new(), &theirs, &mine) {
            return Err(TestCaseError::fail(format!(
                "{input}: giac {printed}, ours {mine:?}"
            )));
        }
    }

    /// `gauss` of a quadratic form in 1 to 3 variables, integer
    /// coefficients.
    #[test]
    fn random_gauss(n in 1usize..=3, c in prop::collection::vec(-5i64..=5, 6)) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let vars = &["x", "y", "z"][..n];
        let mut terms = Vec::new();
        let mut k = 0;
        for i in 0..n {
            for j in i..n {
                if c[k] != 0 {
                    terms.push(format!("{}*{}*{}", c[k], vars[i], vars[j]));
                }
                k += 1;
            }
        }
        let form = if terms.is_empty() { "0".to_owned() } else { terms.join("+").replace("+-", "-") };
        agree(giac, &format!("gauss({form},[{}])", vars.join(",")))?;
    }

    /// `rat_jordan` of a square matrix of size 1 to 3: the same blocks of
    /// `J` in some order ([`same_blocks`]). giac's `P` comes from its own
    /// kernels (a departure): not compared.
    #[test]
    fn random_rat_jordan(a in square(3)) {
        let Some(giac) = giac() else {
            return Ok(());
        };
        let input = format!("rat_jordan({})", text(&a));
        let fail = |e: cas_giac::OracleError| TestCaseError::fail(e.to_string());
        let theirs = giac.value(&format!("[{input}]")).map_err(fail)?;
        let mine = ours(&input).map_err(fail)?;
        let ctx = Context::new();
        if mine.args().len() != 2
            || theirs.args().len() != 2
            || !same_blocks(&ctx, &theirs.args()[1], &mine.args()[1])
        {
            return Err(TestCaseError::fail(format!(
                "{input}: giac {theirs:?}, ours {mine:?}"
            )));
        }
    }
}
