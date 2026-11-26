//! giac's reshaping, bordering, triangles, quadratic forms, changes of
//! basis and `basis` (CAS-639) against giac 2.0: its answers read back are
//! this CAS's, and the departures, each with giac's answer as it prints it
//! and this CAS's as decided.

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

#[test]
fn reshapes() {
    check(&[
        "REDIM([[1,2,3],[3,-9,6],[4,5,6]],1,2)",
        "REDIM([[1,2,3],[3,-9,6],[4,5,6]],4,4)",
        "REDIM([[1,2,3],[3,-9,6],[4,5,6]],[2,5])",
        "REDIM([1,2,3],5)",
        "REDIM([1,2,3],2)",
        "redim([[1,2],[3,4]],3,1)",
        "REDIM([[1,2],[3,4]],1,3)",
        "REDIM([[1,2],[3,4]],2)",
        "REDIM([[1,2],[3,4]],5)",
        "REDIM([[1,2],[3,4]],[3])",
        "REDIM([1,2],[3])",
        "REDIM([1,2],2,2)",
        "REDIM([1,2],[2,2])",
        "REDIM([1,2,3,4],[2,2])",
        "REDIM([1,2,3,4,5],2,3)",
        "REDIM([1,2],0)",
        "REDIM([1,2],-1)",
        "REDIM([[1,2],[3,4]],0)",
        "REDIM([],2)",
        "border([[1,2,3,4],[4,5,6,8],[7,8,9,10]],[1,3,5])",
        "border([[1,2,3],[4,5,6],[7,8,9]],[1,0,1])",
        "border([[1,2],[3,4]],[[1],[2]])",
        "border([[1,2],[3,4]],[[5,7],[6,8]])",
        "border([[1,2]],[x])",
        "semi_augment([[68,-21],[56,59],[1,2]],[[68,-21],[56,59]])",
        "semi_augment([[1,2]],[[3,4],[5,6]])",
        "semi_augment([1,2],[3,4])",
        "semi_augment([[1,2]],[[3,4,5]])",
        "semi_augment([[1,2]],5)",
        "semi_augment([[1,2]],[3,4])",
        "semi_augment([1,2],5)",
        "semi_augment(5,[1,2])",
        "semi_augment([],[1])",
        "newMat(2,3)",
        "newMat(1,1)",
        "newMat(0,2)",
        "newMat(0,0)",
        "newMat(0,-1)",
    ]);
}

#[test]
fn triangles_and_transposes() {
    check(&[
        "upper([[1,2,3],[4,5,6],[7,8,9]])",
        "upper([[1,2,3],[4,5,6]])",
        "upper([[1,2],[3,4],[5,6]])",
        "upper([[1]])",
        "upper([[1,2]])",
        "lower([[1,2,3],[4,5,6],[7,8,9]])",
        "lower([[1,2,3],[4,5,6]])",
        "lower([[1,2],[3,4],[5,6]])",
        "upper(\"hello\")",
        "upper(\"Hello World 1é\")",
        "upper(\"\")",
        "upper(\"aBc\")",
        "lower(\"HELLO\")",
        "lower(\"ÉA b\")",
        "lower(\"aBc\")",
        "trn([[1,2+i],[3,4]])",
        "trn([[1,2],[3,4],[5,6]])",
        "trn([[1,2],[3,4]])",
        "trn([1,2+i])",
        "trn([[i,1],[2,3-2*i]])",
    ]);
}

#[test]
fn solutions_and_forms() {
    check(&[
        "simult([[3,1],[3,2]],[[-2],[2]])",
        "simult([[3,1],[3,2]],[[-2,1],[2,-1]])",
        "simult([[1,2,3],[0,1,4],[5,6,0]],[[1],[0],[0]])",
        "simult([[1,2],[3,4]],[[1,0],[0,1]])",
        "simult([[2,0],[0,4]],[[1],[1]])",
        "changebase([[1,2],[1,3]],[[1,1],[0,1]])",
        "changebase([[1,2],[1,3]],[[1,0],[1,1]])",
        "changebase([[a,b],[c,d]],[[1,1],[0,1]])",
        "changebase([[1,2],[3,4]],[[2,0],[0,2]])",
        "a2q([[1,2],[4,4]],[x,y])",
        "a2q([[1,3],[3,4]],[x,y])",
        "a2q([[1,0,0],[0,2,0],[0,0,3]],[x,y,z])",
        "a2q([[a,b],[b,c]],[x,y])",
        "a2q([[1,2],[3,4]],[u,v])",
        "a2q([[1,2],[3,4]],x)",
        "q2a(x^2+2*x*y+2*y^2,[x,y])",
        "q2a(x*y+y*z,[x,y,z])",
        "q2a(3*x^2,[x])",
        "q2a(a*x^2+b*x*y,[x,y])",
        "q2a(x^2-y^2,[x,y])",
        "q2a(2*x*y+z^2,[x,y,z])",
        "q2a(0,[x,y])",
        "q2a(x^2,x)",
        "basis([[1,0],[0,1]])",
        "basis([[0,0],[0,0]])",
        "basis([[1,2],[2,4]])",
        "basis([[0,1],[1,0]])",
        "basis([[2,4],[1,3]])",
        "basis([[0,1,2],[0,2,5],[3,0,0]])",
        "basis([1,2])",
        "basis([[1,2],[3,4]],x)",
        "charpoly([[1,2],[3,4]],fadeev)",
        "charpoly([[1,2],[3,4]],x,fadeev)",
        "charpoly([[1,2],[3,4]],y,fadeev)",
        "pcar([[1,2],[3,4]],fadeev)",
        "pcar([[2,1],[1,2]],fadeev)",
        "convert([[1,2],[3,4]],array)",
        "convert([1,2,3],array)",
        "convert(x+1,array)",
        "convert(5,array)",
    ]);
}

/// The vectors of `e`, a list or a set, in reduced row echelon form.
fn reduced(e: &Expr) -> Option<Expr> {
    let vectors = Expr::apply(&Symbol::LIST.into(), e.args());
    let rref = Expr::apply(&Symbol::ROW_REDUCE.into(), &[vectors]);
    Context::new().eval(&rref).ok()
}

/// Where giac's rows are fraction-free multiples: the same span, giac's
/// answer and this CAS's having one reduced row echelon form.
#[test]
fn same_spans() {
    let Some(giac) = giac() else {
        return;
    };
    let inputs = [
        "basis([[1,2,3],[4,5,6],[7,8,9],[10,11,12]])",
        "basis([[2,1,0],[1,3,1],[0,1,4]])",
        "basis([[1,2,3],[2,4,6],[1,1,1]])",
        "basis([[a,1],[1,a]])",
        "basis([[1,2],[3,4],[5,6]])",
    ];
    let mut wrong = Vec::new();
    for input in inputs {
        let theirs = giac.value(input).ok().and_then(|e| reduced(&e));
        let mine = ours(input).ok().and_then(|e| reduced(&e));
        if theirs.is_none() || theirs != mine {
            wrong.push(format!("{input}: giac's span {theirs:?}, ours {mine:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departures of CAS-639: giac errs (its error text by `rror`, or no
/// answer), takes `|r|` and `|c|`, or takes symbols as real; this CAS's as
/// decided, the call staying with a warning where giac errs.
fn departed(cases: &[(&str, Option<&str>, &str)]) {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for &(input, answer, decided) in cases {
        let printed = giac.eval(input);
        let recorded = match (answer, &printed) {
            (None, Ok(p)) => p.contains("rror"),
            (None, Err(_)) => true,
            (Some(a), Ok(p)) => p.replace("matrix[", "[").replace(' ', "") == a,
            (Some(_), Err(_)) => false,
        };
        if !recorded {
            wrong.push(format!(
                "{input}: giac prints {printed:?}, recorded {answer:?}"
            ));
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

/// Departures of `REDIM`, `border`, `semi_augment` and `simult`.
#[test]
fn reshape_departures() {
    departed(&[
        // Dimensions below 1, or no vector or matrix, where giac errs.
        ("REDIM([[1,2],[3,4]],0,0)", None, "REDIM([[1,2],[3,4]],0,0)"),
        (
            "REDIM([[1,2],[3,4]],-1,2)",
            None,
            "REDIM([[1,2],[3,4]],-1,2)",
        ),
        ("REDIM(5,2)", None, "REDIM(5,2)"),
        ("REDIM([[1,2],[3,4]],2.5)", None, "REDIM([[1,2],[3,4]],2.5)"),
        // Rows that do not fit, or no matrix, where giac errs.
        (
            "border([[1,2],[3,4]],[1,2,3])",
            None,
            "border([[1,2],[3,4]],[1,2,3])",
        ),
        (
            "border([[1,2],[3,4]],[1])",
            None,
            "border([[1,2],[3,4]],[1])",
        ),
        ("border([1,2],[3])", None, "border([1,2],[3])"),
        ("border([[1,2],[3,4]],5)", None, "border([[1,2],[3,4]],5)"),
        // No list stays, where giac answers 0.
        ("semi_augment(5,6)", Some("0"), "semi_augment(5,6)"),
        // Singular, a vector, not square, or rows that do not fit.
        (
            "simult([[1,1],[1,1]],[[1],[2]])",
            None,
            "simult([[1,1],[1,1]],[[1],[2]])",
        ),
        (
            "simult([[3,1],[3,2]],[-2,2])",
            None,
            "simult([[3,1],[3,2]],[-2,2])",
        ),
        (
            "simult([[1,2],[3,4],[5,6]],[[1],[2],[3]])",
            None,
            "simult([[1,2],[3,4],[5,6]],[[1],[2],[3]])",
        ),
        (
            "simult([[1,2],[3,4]],[[1,2],[3,4],[5,6]])",
            None,
            "simult([[1,2],[3,4]],[[1,2],[3,4],[5,6]])",
        ),
    ]);
}

/// Departures of `newMat`, `upper`, `lower`, `a2q`, `q2a`, `changebase`
/// and `trn`.
#[test]
fn form_departures() {
    departed(&[
        // A dimension below 1 stays, where giac takes its absolute value.
        ("newMat(2,-1)", Some("[[0],[0]]"), "newMat(2,-1)"),
        ("newMat(-2,1)", Some("[[0],[0]]"), "newMat(-2,1)"),
        ("newMat(2,0)", Some("[[],[]]"), "newMat(2,0)"),
        ("newMat(1.5,2)", None, "newMat(1.5,2)"),
        // Neither a matrix nor a string, where giac errs.
        ("upper(5)", None, "upper(5)"),
        ("upper([1,2,3])", None, "upper([1,2,3])"),
        ("upper([])", None, "upper([])"),
        ("lower(5)", None, "lower(5)"),
        ("lower([])", None, "lower([])"),
        // Sizes that do not fit, where giac errs.
        (
            "a2q([[1,2],[3,4],[5,6]],[x,y])",
            None,
            "a2q([[1,2],[3,4],[5,6]],[x,y])",
        ),
        (
            "a2q([[1,2],[2,1]],[x,y,z])",
            None,
            "a2q([[1,2],[2,1]],[x,y,z])",
        ),
        ("a2q([[1,2],[3,4]],[x])", None, "a2q([[1,2],[3,4]],[x])"),
        // Not a quadratic form in the variables, where giac errs.
        (
            "q2a(x^2+2*x*y+2*y^2+x,[x,y])",
            None,
            "q2a(x^2+2*x*y+2*y^2+x,[x,y])",
        ),
        ("q2a(x^3,[x])", None, "q2a(x^3,[x])"),
        ("q2a(x^2+y^2,[x])", None, "q2a(x^2+y^2,[x])"),
        ("q2a(sin(x),[x])", None, "q2a(sin(x),[x])"),
        ("q2a(x*y,[x])", None, "q2a(x*y,[x])"),
        // Singular, or sizes that do not fit, where giac errs.
        (
            "changebase([[1,2],[1,3]],[[1,2],[2,4]])",
            None,
            "changebase([[1,2],[1,3]],[[1,2],[2,4]])",
        ),
        (
            "changebase([[1,2],[1,3]],[[1,1,1],[0,1,1]])",
            None,
            "changebase([[1,2],[1,3]],[[1,1,1],[0,1,1]])",
        ),
        (
            "changebase([[1,2],[1,3]],[[1,2,3],[4,5,6],[7,8,10]])",
            None,
            "changebase([[1,2],[1,3]],[[1,2,3],[4,5,6],[7,8,10]])",
        ),
        // A symbol's conjugate stays, where giac takes symbols as real.
        (
            "trn([[a,b],[c,d]])",
            Some("[[a,c],[b,d]]"),
            "[[conj(a),conj(c)],[conj(b),conj(d)]]",
        ),
        // A scalar stays as written, where giac holds `conj(tran(...))`.
        ("trn(2+i)", Some("conj(tran(2+i))"), "trn(2+i)"),
    ]);
}
