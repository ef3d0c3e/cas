use cas_core::{Complex, Expr, Integer, Rational, Symbol};
use proptest::prelude::*;

use super::{MAX_PRINT_DEPTH, PrintErrorKind as K, fold_arithmetic, print};
use crate::lower::lower;
use crate::names::Names;
use crate::parser::parse;

fn sym(name: &str) -> Expr {
    Expr::from(Symbol::new(name))
}

fn int(n: i64) -> Expr {
    Expr::from(n)
}

fn ratio(n: i64, d: i64) -> Rational {
    Rational::new(n.into(), d.into())
}

fn rat(n: i64, d: i64) -> Expr {
    Expr::from(ratio(n, d))
}

/// The complex number `re + im·i`, with the parts as `(numerator,
/// denominator)`.
fn cx(re: (i64, i64), im: (i64, i64)) -> Expr {
    Expr::from(Complex::new(ratio(re.0, re.1), ratio(im.0, im.1)))
}

fn apply(head: &Expr, args: impl AsRef<[Expr]>) -> Expr {
    Expr::apply(head, args.as_ref())
}

fn app(head: Symbol, args: impl AsRef<[Expr]>) -> Expr {
    apply(&head.into(), args)
}

fn plus(args: impl AsRef<[Expr]>) -> Expr {
    app(Symbol::PLUS, args)
}

fn times(args: impl AsRef<[Expr]>) -> Expr {
    app(Symbol::TIMES, args)
}

fn pow(b: &Expr, n: &Expr) -> Expr {
    app(Symbol::POWER, [b.clone(), n.clone()])
}

/// `f[args…]` for a user's `f`.
fn call(f: &str, args: impl AsRef<[Expr]>) -> Expr {
    apply(&sym(f), args)
}

fn lambda(params: &[&str], body: &Expr) -> Expr {
    let params: Vec<Expr> = params.iter().map(|p| sym(p)).collect();
    app(Symbol::FUNCTION, [app(Symbol::LIST, params), body.clone()])
}

fn printed(e: &Expr) -> String {
    print(e, Names::english()).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

/// The expression `text` lowers to with `names`: one statement, without
/// errors or notes.
fn read_with(text: &str, names: &Names) -> Expr {
    let parse = parse(text);
    assert!(parse.errors().is_empty(), "{text}: {:?}", parse.errors());
    let lowered = lower(&parse, names);
    assert!(
        lowered.diagnostics.is_empty(),
        "{text}: {:?}",
        lowered.diagnostics
    );
    assert_eq!(lowered.statements.len(), 1, "{text}");
    lowered.statements[0].expr.clone()
}

/// Whether `text` has a sign right after an arithmetic operator, outside
/// strings: `x + -y`, `2*-3`, `x^-1`, which the printer never writes (it
/// writes `x - y`, `2*(-3)`, `x^(-1)`).
fn sign_after_operator(text: &str) -> bool {
    let mut in_string = false;
    let mut escaped = false;
    let mut after_operator = false;
    for c in text.chars() {
        if in_string {
            match (escaped, c) {
                (true, _) => escaped = false,
                (false, '\\') => escaped = true,
                (false, '"') => in_string = false,
                _ => {}
            }
            continue;
        }
        match c {
            '"' => {
                in_string = true;
                after_operator = false;
            }
            ' ' => {}
            '+' | '-' if after_operator => return true,
            '+' | '-' | '*' | '/' | '^' => after_operator = true,
            _ => after_operator = false,
        }
    }
    false
}

/// Prints `e` with `names` and reads the text back: the same expression up
/// to the arithmetic the text writes out. Returns the text.
fn round_trip_with(e: &Expr, names: &Names) -> String {
    let text = print(e, names).unwrap_or_else(|err| panic!("{e:?}: {err}"));
    assert!(
        !sign_after_operator(&text),
        "{e:?} printed as {text}: a sign after an operator"
    );
    let back = read_with(&text, names);
    assert_eq!(
        fold_arithmetic(&back),
        fold_arithmetic(e),
        "{e:?} printed as {text} reads back as {back:?}"
    );
    text
}

fn round_trip(e: &Expr) -> String {
    round_trip_with(e, Names::english())
}

/// `(expression, text)`: each prints as its text, which reads back as it.
fn check(cases: &[(Expr, &str)]) {
    let mut wrong = Vec::new();
    for (e, expected) in cases {
        let found = round_trip(e);
        if found != *expected {
            wrong.push(format!("{e:?}: found {found}, expected {expected}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `text` reads back as an expression that prints as `text` again: the
/// printer is the inverse of lowering on it.
fn check_text(texts: &[&str]) {
    let mut wrong = Vec::new();
    for text in texts {
        let found = printed(&read_with(text, Names::english()));
        if found != *text {
            wrong.push(format!("{text}: printed back as {found}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn atoms() {
    let big = Expr::from(Integer::from(10).pow(30));
    check(&[
        (sym("x"), "x"),
        (sym("θ"), "θ"),
        (sym("x_1"), "x_1"),
        (sym("_t"), "_t"),
        // Capitalised user names are escaped; the text drops the `$`.
        (sym("$E"), "E"),
        (sym("$Sin"), "Sin"),
        (int(42), "42"),
        (int(-42), "-42"),
        (big.clone(), "1000000000000000000000000000000"),
        (-&big, "-1000000000000000000000000000000"),
        (rat(3, 4), "3/4"),
        (rat(-3, 4), "-3/4"),
        (Expr::string("plain"), r#""plain""#),
        (
            Expr::string("π \"q\" \\ \n\t\r\u{1}"),
            r#""π \"q\" \\ \n\t\u{d}\u{1}""#,
        ),
        (Expr::string(""), r#""""#),
    ]);
}

#[test]
fn complex_numbers() {
    check(&[
        (cx((0, 1), (1, 1)), "i"),
        (cx((0, 1), (-1, 1)), "-i"),
        (cx((0, 1), (3, 1)), "3*i"),
        (cx((0, 1), (-3, 1)), "-3*i"),
        (cx((0, 1), (1, 2)), "i/2"),
        (cx((0, 1), (-3, 4)), "-3*i/4"),
        (cx((2, 1), (3, 1)), "2 + 3*i"),
        (cx((2, 1), (-1, 1)), "2 - i"),
        (cx((-1, 2), (3, 4)), "-1/2 + 3*i/4"),
        (cx((-2, 1), (-3, 1)), "-2 - 3*i"),
    ]);
}

#[test]
fn constants_and_functions_as_values() {
    let minus_infinity = app(Symbol::DIRECTED_INFINITY, [int(-1)]);
    check(&[
        (Symbol::PI.into(), "pi"),
        (Symbol::E.into(), "e"),
        (Symbol::TRUE.into(), "true"),
        (Symbol::FALSE.into(), "false"),
        (Symbol::INDETERMINATE.into(), "undef"),
        (Symbol::COMPLEX_INFINITY.into(), "infinity"),
        (Symbol::INFINITY.into(), "+infinity"),
        (minus_infinity.clone(), "-infinity"),
        (Symbol::SIN.into(), "sin"),
        (Symbol::LOG.into(), "ln"),
        (Symbol::FACTORIAL.into(), "factorial"),
        (Symbol::REPLACE_ALL.into(), "subst"),
        // A signed infinity is an operand with a sign.
        (plus([sym("x"), Symbol::INFINITY.into()]), "x + (+infinity)"),
        (plus([sym("x"), minus_infinity.clone()]), "x + (-infinity)"),
        (times([int(2), minus_infinity]), "2*(-infinity)"),
        (app(Symbol::LIST, [Symbol::INFINITY.into()]), "[+infinity]"),
        // `-infinity` is a whole form: the negated complex infinity is not it.
        (
            times([int(-1), Symbol::COMPLEX_INFINITY.into()]),
            "-(infinity)",
        ),
        (
            times([int(-1), Symbol::COMPLEX_INFINITY.into(), sym("x")]),
            "-(infinity)*x",
        ),
        (
            plus([sym("x"), times([int(-1), Symbol::COMPLEX_INFINITY.into()])]),
            "x - infinity",
        ),
    ]);
}

/// A product whose first factor has its own sign prints as the product with
/// the sign in front: `Times[Times[-1, x], y]` and `Times[-1, x, y]` read
/// back as each other, so they must print alike everywhere (CAS-127).
#[test]
fn a_first_factor_with_its_own_sign() {
    let [x, y, z] = ["x", "y", "z"].map(sym);
    check(&[
        // Written bare: `-x*y` reads as `(-x)*y`, the nested product.
        (times([times([int(-1), x.clone()]), y.clone()]), "-x*y"),
        (times([rat(1, 2), int(-1), rat(1, 2)]), "-1*(1/2)/2"),
        // After an operator, parenthesised whole.
        (
            plus([x.clone(), times([times([int(-1), y.clone()]), z.clone()])]),
            "x + (-y*z)",
        ),
        (
            plus([
                x.clone(),
                times([int(-1), times([int(-1), y.clone()]), pow(&z, &int(-1))]),
            ]),
            "x - (-y/z)",
        ),
        (
            times([x.clone(), times([times([int(-1), y.clone()]), z.clone()])]),
            "x*(-y*z)",
        ),
    ]);
}

#[test]
fn sums_and_products() {
    let [x, y, z] = ["x", "y", "z"].map(sym);
    check(&[
        (plus([x.clone(), times([int(-2), y.clone()])]), "x - 2*y"),
        (plus([x.clone(), times([int(-1), y.clone()])]), "x - y"),
        (plus([x.clone(), int(-3)]), "x - 3"),
        (plus([x.clone(), rat(-1, 2)]), "x - 1/2"),
        (plus([x.clone(), cx((0, 1), (-2, 1))]), "x - 2*i"),
        (plus([x.clone(), cx((1, 1), (-2, 1))]), "x + (1 - 2*i)"),
        (plus([times([int(-1), x.clone()]), y.clone()]), "-x + y"),
        (
            plus([x.clone(), plus([y.clone(), z.clone()])]),
            "x + (y + z)",
        ),
        (
            plus([plus([x.clone(), y.clone()]), z.clone()]),
            "(x + y) + z",
        ),
        (
            plus([x.clone(), times([int(-1), plus([y.clone(), z.clone()])])]),
            "x - (y + z)",
        ),
        (
            plus([x.clone(), times([rat(-1, 2), y.clone(), pow(&z, &int(-1))])]),
            "x - y/(2*z)",
        ),
        (times([rat(1, 2), x.clone()]), "x/2"),
        (
            times([rat(-3, 4), x.clone(), pow(&y, &int(-2))]),
            "-3*x/(4*y^2)",
        ),
        (times([rat(1, 2), pow(&x, &int(-1))]), "1/(2*x)"),
        (pow(&x, &int(-1)), "1/x"),
        (pow(&x, &int(-2)), "1/x^2"),
        (times([x.clone(), pow(&y, &int(-1))]), "x/y"),
        (times([pow(&x, &int(-1)), y.clone()]), "y/x"),
        (times([int(-1), x.clone()]), "-x"),
        (times([int(-1), x.clone(), y.clone()]), "-x*y"),
        (times([int(-2), x.clone()]), "-2*x"),
        (times([int(-1), plus([x.clone(), y.clone()])]), "-(x + y)"),
        (times([int(-1), pow(&x, &int(2))]), "-x^2"),
        (times([int(-1), times([int(-1), x.clone()])]), "-(-x)"),
        (times([int(-1), pow(&x, &int(-1))]), "-1/x"),
        (times([x.clone(), int(-2)]), "x*(-2)"),
        (times([x.clone(), times([y.clone(), z.clone()])]), "x*(y*z)"),
        (times([x.clone(), rat(1, 2)]), "x*(1/2)"),
        (
            times([
                plus([x.clone(), y.clone()]),
                plus([x.clone(), times([int(-1), y.clone()])]),
            ]),
            "(x + y)*(x - y)",
        ),
        (times([cx((0, 1), (2, 1)), x.clone()]), "2*i*x"),
        (times([cx((0, 1), (-1, 1)), x.clone()]), "-i*x"),
        (times([cx((0, 1), (3, 4)), x.clone()]), "3*i*x/4"),
        (times([cx((2, 1), (3, 1)), x.clone()]), "(2 + 3*i)*x"),
        (times([int(1), x.clone()]), "1*x"),
        (times([int(2), int(3)]), "2*3"),
        (times([int(0), x.clone()]), "0*x"),
        // `-3` would read back as the number: the 1 is written.
        (times([int(-1), int(3)]), "-1*3"),
        (times([int(-1), int(0)]), "-1*0"),
        (times([int(-1), int(-3)]), "-1*(-3)"),
        (times([int(-1), cx((0, 1), (1, 1))]), "-1*i"),
        (times([int(-1), int(3), x.clone()]), "-1*3*x"),
        (plus([x.clone(), times([int(-1), int(3)])]), "x - 3"),
        (times([int(2), pow(&int(3), &int(-1))]), "2/3"),
        (times([x.clone(), pow(&int(-2), &int(-1))]), "x/(-2)"),
        (
            times([x.clone(), pow(&plus([y.clone(), z.clone()]), &int(-1))]),
            "x/(y + z)",
        ),
        (
            times([x.clone(), pow(&times([y.clone(), z.clone()]), &int(-1))]),
            "x/(y*z)",
        ),
        // A zero in a denominator absorbs nothing: found by fuzzing.
        (
            times([
                rat(101, 100),
                pow(&Symbol::E.into(), &int(-1)),
                pow(&times([int(0), z.clone()]), &int(-1)),
            ]),
            "101/(100*e*(0*z))",
        ),
        (pow(&rat(1, 2), &int(-1)), "1/(1/2)"),
        (
            times([y.clone(), pow(&pow(&x, &int(-1)), &int(-1))]),
            "y/(1/x)",
        ),
    ]);
}

#[test]
fn powers() {
    let [x, y, n] = ["x", "y", "n"].map(sym);
    check(&[
        (pow(&x, &int(2)), "x^2"),
        (pow(&plus([x.clone(), int(1)]), &int(2)), "(x + 1)^2"),
        (pow(&int(-2), &x), "(-2)^x"),
        (pow(&x, &rat(1, 2)), "x^(1/2)"),
        (pow(&x, &rat(-1, 2)), "x^(-1/2)"),
        (pow(&x, &times([int(-1), y.clone()])), "x^(-y)"),
        (pow(&int(2), &pow(&int(3), &int(2))), "2^3^2"),
        (pow(&pow(&int(2), &int(3)), &int(2)), "(2^3)^2"),
        (pow(&rat(1, 2), &x), "(1/2)^x"),
        (pow(&x, &cx((0, 1), (1, 1))), "x^i"),
        (pow(&cx((0, 1), (1, 1)), &int(2)), "i^2"),
        (pow(&cx((0, 1), (2, 1)), &int(2)), "(2*i)^2"),
        (pow(&Symbol::SIN.into(), &int(2)), "(sin)^2"),
        (pow(&app(Symbol::SIN, [x.clone()]), &int(2)), "sin(x)^2"),
        (pow(&Symbol::PI.into(), &int(2)), "pi^2"),
        (pow(&app(Symbol::FACTORIAL, [n.clone()]), &int(2)), "n!^2"),
        (pow(&app(Symbol::LIST, [x.clone()]), &int(2)), "[x]^2"),
        (app(Symbol::ELEMENTWISE_POWER, [x.clone(), int(2)]), "x.^2"),
        (
            pow(&x, &app(Symbol::ELEMENTWISE_POWER, [y.clone(), int(2)])),
            "x^y.^2",
        ),
        (pow(&x, &pow(&y, &int(-1))), "x^(1/y)"),
    ]);
}

#[test]
fn elementwise_and_factorials() {
    let [a, b, c, n] = ["a", "b", "c", "n"].map(sym);
    let et = |args: &[Expr]| app(Symbol::ELEMENTWISE_TIMES, args);
    let ed = |args: &[Expr]| app(Symbol::ELEMENTWISE_DIVIDE, args);
    let fact = |e: &Expr| app(Symbol::FACTORIAL, [e.clone()]);
    check(&[
        (et(&[a.clone(), b.clone()]), "a.*b"),
        (et(&[a.clone(), b.clone(), c.clone()]), "a.*b.*c"),
        (ed(&[a.clone(), b.clone()]), "a./b"),
        (ed(&[ed(&[a.clone(), b.clone()]), c.clone()]), "(a./b)./c"),
        (et(&[times([a.clone(), b.clone()]), c.clone()]), "(a*b).*c"),
        (times([et(&[a.clone(), b.clone()]), c.clone()]), "(a.*b)*c"),
        (et(&[a.clone(), times([int(-1), b.clone()])]), "a.*(-b)"),
        (et(&[times([int(-1), a.clone()]), b.clone()]), "-a.*b"),
        (et(&[int(2), a.clone()]), "2.*a"),
        (fact(&n), "n!"),
        (fact(&plus([n.clone(), int(1)])), "(n + 1)!"),
        (fact(&fact(&n)), "n!!"),
        (fact(&int(-3)), "(-3)!"),
        (
            app(Symbol::FACTORIAL, [a.clone(), b.clone()]),
            "factorial(a, b)",
        ),
        (app(Symbol::EQUAL, [fact(&n), int(3)]), "n! == 3"),
    ]);
}

#[test]
#[allow(
    clippy::many_single_char_names,
    reason = "symbols named as in mathematics"
)]
fn relations_and_logic() {
    let [a, b, c, p, q, r, x] = ["a", "b", "c", "p", "q", "r", "x"].map(sym);
    let rel = |s: Symbol, l: &Expr, r: &Expr| app(s, &[l.clone(), r.clone()]);
    let and = |args: &[Expr]| app(Symbol::AND, args);
    let or = |args: &[Expr]| app(Symbol::OR, args);
    let not = |e: &Expr| app(Symbol::NOT, [e.clone()]);
    check(&[
        (rel(Symbol::EQUATION, &x, &int(1)), "x = 1"),
        (rel(Symbol::EQUAL, &x, &int(-1)), "x == -1"),
        (rel(Symbol::UNEQUAL, &x, &int(1)), "x != 1"),
        (rel(Symbol::LESS, &x, &int(1)), "x < 1"),
        (rel(Symbol::LESS_EQUAL, &x, &int(1)), "x <= 1"),
        (rel(Symbol::GREATER, &x, &int(1)), "x > 1"),
        (rel(Symbol::GREATER_EQUAL, &x, &int(1)), "x >= 1"),
        (
            rel(
                Symbol::ELEMENT,
                &x,
                &app(Symbol::FINITE_SET, [int(1), int(2)]),
            ),
            "x in {1, 2}",
        ),
        (
            and(&[
                rel(Symbol::LESS, &int(0), &x),
                rel(Symbol::LESS_EQUAL, &x, &int(1)),
            ]),
            "0 < x <= 1",
        ),
        (
            and(&[
                rel(Symbol::LESS, &a, &b),
                rel(Symbol::LESS, &b, &c),
                rel(Symbol::ELEMENT, &c, &p),
            ]),
            "a < b < c in p",
        ),
        (
            and(&[
                rel(Symbol::LESS, &int(0), &x),
                rel(Symbol::LESS, &a, &int(1)),
            ]),
            "0 < x and a < 1",
        ),
        // An equation never chains.
        (
            and(&[rel(Symbol::EQUATION, &a, &b), rel(Symbol::EQUATION, &b, &c)]),
            "a = b and b = c",
        ),
        (
            and(&[
                and(&[rel(Symbol::LESS, &a, &b), rel(Symbol::LESS, &b, &c)]),
                p.clone(),
            ]),
            "a < b < c and p",
        ),
        (
            rel(Symbol::LESS, &rel(Symbol::LESS, &a, &b), &c),
            "(a < b) < c",
        ),
        (
            rel(Symbol::EQUATION, &rel(Symbol::EQUATION, &a, &b), &c),
            "(a = b) = c",
        ),
        (
            rel(Symbol::EQUATION, &a, &app(Symbol::RANGE, [int(1), int(2)])),
            "a = 1..2",
        ),
        (and(&[p.clone(), q.clone()]), "p and q"),
        (
            or(&[and(&[p.clone(), q.clone()]), r.clone()]),
            "p and q or r",
        ),
        (
            and(&[or(&[p.clone(), q.clone()]), r.clone()]),
            "(p or q) and r",
        ),
        (
            and(&[and(&[p.clone(), q.clone()]), r.clone()]),
            "(p and q) and r",
        ),
        (
            and(&[p.clone(), and(&[q.clone(), r.clone()])]),
            "p and (q and r)",
        ),
        (not(&and(&[p.clone(), q.clone()])), "not (p and q)"),
        (not(&rel(Symbol::LESS, &a, &b)), "not a < b"),
        (and(&[not(&p), q.clone()]), "not p and q"),
        (not(&not(&p)), "not not p"),
        (rel(Symbol::EQUAL, &not(&p), &q), "(not p) == q"),
        (
            rel(Symbol::LESS, &x, &plus([a.clone(), b.clone()])),
            "x < a + b",
        ),
    ]);
}

#[test]
fn brackets_ranges_and_indexes() {
    let [m, n, v, x] = ["m", "n", "v", "x"].map(sym);
    let range = |a: &Expr, b: &Expr| app(Symbol::RANGE, [a.clone(), b.clone()]);
    let at = |args: Vec<Expr>| app(Symbol::AT, args);
    check(&[
        (range(&int(0), &n), "0..n"),
        (range(&int(1), &plus([n.clone(), int(-1)])), "1 .. n - 1"),
        (range(&int(-1), &int(1)), "-1..1"),
        (range(&rat(1, 2), &int(3)), "1/2..3"),
        (range(&range(&int(0), &int(1)), &int(2)), "(0..1)..2"),
        (app(Symbol::LIST, []), "[]"),
        (app(Symbol::LIST, [int(1), x.clone()]), "[1, x]"),
        (app(Symbol::FINITE_SET, []), "{}"),
        (app(Symbol::FINITE_SET, [x.clone(), int(1)]), "{x, 1}"),
        (at(vec![v.clone(), int(0)]), "v[0]"),
        (at(vec![v.clone(), int(-1)]), "v[-1]"),
        (at(vec![m.clone(), int(0), int(1)]), "m[0, 1]"),
        (at(vec![v.clone(), range(&int(1), &int(2))]), "v[1..2]"),
        (at(vec![at(vec![m.clone(), int(0)]), int(1)]), "m[0][1]"),
        (at(vec![plus([m.clone(), v.clone()]), int(0)]), "(m + v)[0]"),
        (at(vec![v.clone()]), "v[]"),
        (at(vec![app(Symbol::LIST, [x.clone()]), int(0)]), "[x][0]"),
        (at(vec![int(-2), int(0)]), "(-2)[0]"),
    ]);
}

#[test]
#[allow(
    clippy::many_single_char_names,
    reason = "symbols named as in mathematics"
)]
fn calls() {
    let [a, b, f, u, x, y] = ["a", "b", "f", "u", "x", "y"].map(sym);
    let sin = |e: &Expr| app(Symbol::SIN, [e.clone()]);
    let subst = |args: &[Expr]| app(Symbol::REPLACE_ALL, args);
    let pi_of = |args: Vec<Expr>| apply(&Symbol::PI.into(), args);
    check(&[
        (call("f", [x.clone()]), "f(x)"),
        (call("f", []), "f()"),
        (call("f", [x.clone(), y.clone()]), "f(x, y)"),
        (call("$F", [x.clone()]), "F(x)"),
        (sin(&x), "sin(x)"),
        (app(Symbol::SQRT, [x.clone()]), "sqrt(x)"),
        (app(Symbol::HOLD, [x.clone()]), "hold(x)"),
        // The call-or-product rule would read these as products (§3.1).
        (call("f", [f.clone()]), "call(f, f)"),
        (call("f", [plus([f.clone(), int(1)])]), "call(f, f + 1)"),
        (call("x", [plus([x.clone(), int(1)])]), "call(x, x + 1)"),
        (call("f", [call("g", [f.clone()])]), "call(f, g(f))"),
        (
            call("f", [app(Symbol::AT, [f.clone(), int(0)])]),
            "call(f, f[0])",
        ),
        (call("f", [lambda(&["f"], &f)]), "call(f, f -> f)"),
        (call("f", [call("f", [f.clone()])]), "call(f, call(f, f))"),
        (sin(&Symbol::SIN.into()), "call(sin, sin)"),
        (pi_of(vec![x.clone()]), "call(pi, x)"),
        (apply(&cx((0, 1), (1, 1)), [x.clone()]), "call(i, x)"),
        (
            apply(&Symbol::INFINITY.into(), [x.clone()]),
            "call(+infinity, x)",
        ),
        (
            apply(&Symbol::INFINITY.into(), [x.clone(), y.clone()]),
            "call(+infinity, x, y)",
        ),
        // Only a single argument can read as a product.
        (call("f", [f.clone(), x.clone()]), "f(f, x)"),
        (call("f", [call("f", [x.clone()])]), "f(f(x))"),
        (
            call("f", [apply(&call("f", [a.clone()]), [b.clone()])]),
            "f(f(a)(b))",
        ),
        (pi_of(vec![x.clone(), y.clone()]), "pi(x, y)"),
        (pi_of(vec![]), "pi()"),
        // Heads that cannot stand before parentheses.
        (apply(&int(2), [x.clone()]), "call(2, x)"),
        (
            apply(&plus([a.clone(), b.clone()]), [x.clone()]),
            "call(a + b, x)",
        ),
        (apply(&Expr::string("s"), []), "call(\"s\")"),
        (apply(&app(Symbol::LIST, []), [x.clone()]), "call([], x)"),
        (
            apply(&app(Symbol::FACTORIAL, [x.clone()]), [y.clone()]),
            "call(x!, y)",
        ),
        (
            apply(&app(Symbol::IF, [a.clone(), f.clone()]), [x.clone()]),
            "call(if a then f end, x)",
        ),
        // Heads that can.
        (apply(&call("f", [a.clone()]), [b.clone()]), "f(a)(b)"),
        (
            apply(&app(Symbol::AT, [sym("v"), int(0)]), [x.clone()]),
            "v[0](x)",
        ),
        (
            apply(&lambda(&["x"], &pow(&x, &int(2))), &[int(3)]),
            "(x -> x^2)(3)",
        ),
        (
            apply(&pi_of(vec![x.clone()]), [y.clone()]),
            "call(pi, x)(y)",
        ),
        (apply(&call("f", [f.clone()]), [y.clone()]), "call(f, f)(y)"),
        // `subst` reads an equation or a list after its first argument as
        // rules; parentheses keep them values.
        (
            subst(&[u.clone(), app(Symbol::EQUATION, [x.clone(), int(2)])]),
            "subst(u, (x = 2))",
        ),
        (
            subst(&[
                u.clone(),
                app(Symbol::LIST, [app(Symbol::EQUATION, [x.clone(), int(2)])]),
            ]),
            "subst(u, ([x = 2]))",
        ),
        (
            subst(&[app(Symbol::EQUATION, [x.clone(), int(2)]), y.clone()]),
            "subst(x = 2, y)",
        ),
        (subst(&[u.clone(), y.clone()]), "subst(u, y)"),
    ]);
}

/// The derivatives of unknown functions, which have no name (CAS-DOC-7
/// L13): printed as the `diff` or `subst` that evaluates back to them
/// (CAS-258), not read back by lowering alone.
#[test]
fn derivatives() {
    use cas_core::Context;

    let [f, x] = ["f", "x"].map(sym);
    let derivative = |n: i64, h: &Expr, a: &Expr| {
        let d = app(Symbol::DERIVATIVE, [int(n)]);
        apply(&apply(&d, [h.clone()]), [a.clone()])
    };
    let cases = [
        (derivative(1, &f, &x), "diff(f(x), x)"),
        (derivative(2, &f, &x), "diff(f(x), x, x)"),
        (derivative(3, &f, &sym("t")), "diff(f(t), t, t, t)"),
        (derivative(5, &f, &x), "diff(f(x), [x, 5])"),
        (derivative(1, &f, &int(0)), "subst(diff(f(x), x), x = 0)"),
        (
            derivative(2, &f, &plus([int(1), x.clone()])),
            "subst(diff(f(t), t, t), t = 1 + x)",
        ),
        // The variable is neither the function nor in the argument.
        (derivative(1, &x, &int(0)), "subst(diff(x(t), t), t = 0)"),
        (
            derivative(1, &f, &call("g", [x.clone(), sym("t")])),
            "subst(diff(f(u), u), u = g(x, t))",
        ),
        (
            derivative(1, &f, &Expr::from(Symbol::PI)),
            "subst(diff(f(x), x), x = pi)",
        ),
        (
            times([int(2), derivative(1, &f, &times([int(2), x.clone()]))]),
            "2*subst(diff(f(t), t), t = 2*x)",
        ),
    ];
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (e, text) in &cases {
        let found = printed(e);
        let back = ctx.eval(&read_with(&found, Names::english()));
        if found != *text || back.as_ref() != Ok(e) {
            wrong.push(format!("{e:?}: {found}, evaluated back {back:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
    // A bare derivative, a function and not a value, has no form.
    let bare = apply(&app(Symbol::DERIVATIVE, [int(1)]), [f]);
    assert!(print(&bare, Names::english()).is_err());
}

#[test]
fn lambdas_and_conditionals() {
    let [c, d, x, y] = ["c", "d", "x", "y"].map(sym);
    let null: Expr = Symbol::NULL.into();
    let if_ = |args: &[Expr]| app(Symbol::IF, args);
    check(&[
        (lambda(&["x"], &pow(&x, &int(2))), "x -> x^2"),
        (
            lambda(&["x", "y"], &plus([x.clone(), y.clone()])),
            "(x, y) -> x + y",
        ),
        (lambda(&[], &int(1)), "() -> 1"),
        (lambda(&["$X"], &sym("$X")), "X -> X"),
        (
            lambda(&["x"], &lambda(&["y"], &plus([x.clone(), y.clone()]))),
            "x -> y -> x + y",
        ),
        (plus([lambda(&["x"], &x), int(1)]), "(x -> x) + 1"),
        (call("f", [lambda(&["x"], &x), int(1)]), "f(x -> x, 1)"),
        (app(Symbol::LIST, [lambda(&["x"], &x)]), "[x -> x]"),
        (if_(&[c.clone(), int(1)]), "if c then 1 end"),
        (if_(&[c.clone(), int(1), int(2)]), "if c then 1 else 2 end"),
        (
            if_(&[c.clone(), int(1), if_(&[d.clone(), int(2), int(3)])]),
            "if c then 1 elif d then 2 else 3 end",
        ),
        (
            if_(&[c.clone(), int(1), if_(&[d.clone(), int(2)])]),
            "if c then 1 elif d then 2 end",
        ),
        (
            if_(&[c.clone(), null.clone(), int(2)]),
            "if c then else 2 end",
        ),
        (
            if_(&[c.clone(), int(1), null.clone()]),
            "if c then 1 else end",
        ),
        (if_(&[c.clone(), null]), "if c then end"),
        (
            if_(&[c.clone(), if_(&[d.clone(), int(1)]), int(2)]),
            "if c then if d then 1 end else 2 end",
        ),
        (
            plus([if_(&[c.clone(), int(1), int(2)]), int(1)]),
            "(if c then 1 else 2 end) + 1",
        ),
        (
            lambda(&["x"], &if_(&[c.clone(), x.clone(), int(0)])),
            "x -> if c then x else 0 end",
        ),
        (
            if_(&[lambda(&["x"], &x), lambda(&["y"], &y)]),
            "if x -> x then y -> y end",
        ),
    ]);
}

/// Text that the printer gives back unchanged: lowering builds these
/// structures exactly.
#[test]
fn text_reads_back_and_prints_the_same() {
    check_text(&[
        "x - 2*y",
        "x - y/(2*z)",
        "-x^2 + 2*x - 1",
        "1/x + 1/(x*y)",
        "(x + 1)*(x - 1)",
        "2^3^2",
        "sin(x)^2 + cos(x)^2",
        "f(a)(b)",
        "call(x, x + 1)",
        "0 < x <= 1",
        "x -> y -> x + y",
        "[1, 2, 3][0]",
        "{1, 2}",
        "n! == 3",
        "if x < 0 then -x elif x == 0 then 0 else x end",
        "not (p and q) or r",
        "subst(x^2, (x = 2))",
        "(a.*b)./c",
        "E = m*c^2",
        "x + (-infinity)",
    ]);
}

#[test]
fn locales() {
    let names = Names::with_locale("[builtins]\nSin: sinus\n[constants]\nimaginary: j\n").unwrap();
    let x = sym("x");
    assert_eq!(
        round_trip_with(&app(Symbol::SIN, [x.clone()]), &names),
        "sinus(x)"
    );
    assert_eq!(round_trip_with(&cx((2, 1), (3, 1)), &names), "2 + 3*j");
    // English names stay reserved: a user's `sin` has no name.
    let error = print(&call("sin", [x]), &names).unwrap_err();
    assert_eq!((error.part, error.kind), (sym("sin"), K::NoName));
}

#[test]
fn errors() {
    let [x, y] = ["x", "y"].map(sym);
    let blank = app(Symbol::BLANK, []);
    let cases = [
        // Internal builtins have no name.
        (
            Expr::from(Symbol::PLUS),
            Expr::from(Symbol::PLUS),
            K::NoName,
        ),
        (Symbol::NULL.into(), Symbol::NULL.into(), K::NoName),
        (Symbol::I.into(), Symbol::I.into(), K::NoName),
        // Names that would read back as other symbols.
        (sym("sin"), sym("sin"), K::NoName),
        (sym("pi"), sym("pi"), K::NoName),
        (sym("Foo"), sym("Foo"), K::NoName),
        (sym("$x"), sym("$x"), K::NoName),
        (sym("a b"), sym("a b"), K::NoName),
        (sym("if"), sym("if"), K::NoName),
        (sym(""), sym(""), K::NoName),
        (sym("$"), sym("$"), K::NoName),
        (sym("e\u{301}"), sym("e\u{301}"), K::NoName),
        (call("sin", [x.clone()]), sym("sin"), K::NoName),
        // Operators with another number of arguments, patterns, statements.
        (times([x.clone()]), times([x.clone()]), K::NoSyntax),
        (plus([]), plus([]), K::NoSyntax),
        (
            app(Symbol::POWER, [x.clone(), y.clone(), int(2)]),
            app(Symbol::POWER, [x.clone(), y.clone(), int(2)]),
            K::NoSyntax,
        ),
        (
            app(Symbol::PATTERN, [x.clone(), blank.clone()]),
            app(Symbol::PATTERN, [x.clone(), blank.clone()]),
            K::NoSyntax,
        ),
        (
            app(Symbol::SET, [x.clone(), int(1)]),
            app(Symbol::SET, [x.clone(), int(1)]),
            K::NoSyntax,
        ),
        (
            app(Symbol::DIRECTED_INFINITY, [int(1)]),
            app(Symbol::DIRECTED_INFINITY, [int(1)]),
            K::NoSyntax,
        ),
        (
            app(Symbol::FUNCTION, [x.clone(), x.clone()]),
            app(Symbol::FUNCTION, [x.clone(), x.clone()]),
            K::NoSyntax,
        ),
        (
            app(Symbol::IF, [x.clone()]),
            app(Symbol::IF, [x.clone()]),
            K::NoSyntax,
        ),
        // Inside something else: the part is the innermost.
        (
            call("f", [plus([int(1), times([x.clone()])])]),
            times([x.clone()]),
            K::NoSyntax,
        ),
        (
            app(Symbol::LIST, [Symbol::NULL.into()]),
            Symbol::NULL.into(),
            K::NoName,
        ),
        // Lambda parameters are names of the user's.
        (lambda(&["Pi"], &x), Symbol::PI.into(), K::Parameter),
        (
            app(Symbol::FUNCTION, [app(Symbol::LIST, [int(1)]), x.clone()]),
            int(1),
            K::Parameter,
        ),
        (lambda(&["sin"], &x), sym("sin"), K::Parameter),
    ];
    for (e, part, kind) in cases {
        let error = print(&e, Names::english()).expect_err(&format!("{e:?}"));
        assert_eq!((&error.part, error.kind), (&part, kind), "{e:?}");
        let _ = error.to_string();
    }
    let error = print(&plus([x.clone()]), Names::english()).unwrap_err();
    assert_eq!(
        error.to_string(),
        "the language has no syntax for Plus with 1 argument"
    );
}

/// Nested to the deepest that prints, in the forms that take the most
/// levels of the parser per level of the expression: the text parses. One
/// more level is too deep.
#[test]
fn depth() {
    let leaf = cx((2, 1), (3, 4));
    let shapes: [fn(Expr) -> Expr; 4] = [
        |e| plus([sym("x"), e]),
        |e| plus([sym("x"), times([int(-1), e])]),
        |e| app(Symbol::IF, [sym("x"), e]),
        |e| times([sym("x"), pow(&e, &sym("x"))]),
    ];
    for wrap in shapes {
        let mut e = leaf.clone();
        while usize::from(wrap(e.clone()).depth()) <= MAX_PRINT_DEPTH {
            e = wrap(e);
        }
        assert!(usize::from(e.depth()) >= MAX_PRINT_DEPTH - 1);
        round_trip(&e);
        let deeper = wrap(e);
        if usize::from(deeper.depth()) > MAX_PRINT_DEPTH {
            let error = print(&deeper, Names::english()).unwrap_err();
            assert_eq!((error.part, error.kind), (deeper, K::TooDeep));
        }
    }
}

/// What a text reads back as may nest deeper than what printed it: two
/// divisors, one of them a product, read back as one parenthesised product
/// of products, `j/((z*w)*z)` (the `lang_print` fuzz target, CAS-165). At
/// the limit the expression prints and reads back, and what it reads back
/// as is one level too deep to print.
#[test]
fn reading_back_may_nest_deeper() {
    let [j, z, w] = ["j", "z", "w"].map(sym);
    let quotient = times([j, pow(&times([z.clone(), w]), &int(-1)), pow(&z, &int(-1))]);
    let mut e = quotient;
    while usize::from(app(Symbol::LIST, [e.clone()]).depth()) <= MAX_PRINT_DEPTH {
        e = app(Symbol::LIST, [e]);
    }
    let text = round_trip(&e);
    assert!(text.contains("[j/((z*w)*z)]"), "{text}");
    let back = read_with(&text, Names::english());
    assert_eq!(back.depth(), e.depth() + 1);
    let error = print(&back, Names::english()).unwrap_err();
    assert_eq!(error.kind, K::TooDeep);
}

/// `fold_arithmetic` folds each part once (CAS-176): it folded a power's
/// exponent twice, so a tower of `n` powers took `2^n` folds, and the
/// `lang_print` fuzz target timed out on towers through parentheses and
/// applications (`t^(I(A(e tt^t^t^t^ …`). A tower of 120 round-trips at
/// once.
#[test]
fn towers_fold_once() {
    let start = std::time::Instant::now();
    let mut tower = sym("x");
    for _ in 0..120 {
        tower = pow(&sym("y"), &tower);
    }
    assert_eq!(fold_arithmetic(&tower), tower);
    round_trip(&tower);
    let elapsed = start.elapsed();
    assert!(elapsed < std::time::Duration::from_secs(5), "{elapsed:?}");
}

/// Every form in every operand position of every other: the precedence
/// pairs, with the signs, numbers and heads that need parentheses.
#[test]
fn every_form_in_every_position() {
    let [a, b, c, f] = ["a", "b", "c", "f"].map(sym);
    let fillers = fillers(&a, &b, &c, &f);
    let contexts: Vec<fn(&Expr) -> Expr> = vec![
        |h| h.clone(),
        |h| plus([h.clone(), sym("y")]),
        |h| plus([sym("y"), h.clone()]),
        |h| plus([sym("y"), times([int(-1), h.clone()])]),
        |h| plus([sym("y"), times([int(-2), h.clone()])]),
        |h| times([h.clone(), sym("y")]),
        |h| times([sym("y"), h.clone()]),
        |h| times([int(-1), h.clone()]),
        |h| times([int(-1), h.clone(), sym("y")]),
        |h| times([sym("y"), pow(h, &int(-1))]),
        |h| times([sym("y"), pow(h, &int(-2))]),
        |h| times([h.clone(), pow(&sym("y"), &int(-1))]),
        |h| times([rat(1, 2), h.clone()]),
        |h| times([rat(-3, 4), h.clone(), pow(&sym("y"), &int(-1))]),
        |h| times([cx((0, 1), (2, 1)), h.clone()]),
        |h| pow(h, &sym("y")),
        |h| pow(&sym("y"), h),
        |h| pow(h, &int(-1)),
        |h| app(Symbol::ELEMENTWISE_POWER, [h.clone(), sym("y")]),
        |h| app(Symbol::ELEMENTWISE_POWER, [sym("y"), h.clone()]),
        |h| app(Symbol::ELEMENTWISE_TIMES, [h.clone(), sym("y")]),
        |h| app(Symbol::ELEMENTWISE_TIMES, [sym("y"), h.clone()]),
        |h| app(Symbol::ELEMENTWISE_DIVIDE, [h.clone(), sym("y")]),
        |h| app(Symbol::ELEMENTWISE_DIVIDE, [sym("y"), h.clone()]),
        |h| app(Symbol::FACTORIAL, [h.clone()]),
        |h| app(Symbol::EQUATION, [h.clone(), sym("y")]),
        |h| app(Symbol::EQUATION, [sym("y"), h.clone()]),
        |h| app(Symbol::LESS, [h.clone(), sym("y")]),
        |h| app(Symbol::LESS, [sym("y"), h.clone()]),
        |h| app(Symbol::ELEMENT, [h.clone(), sym("y")]),
        |h| {
            app(
                Symbol::AND,
                &[
                    app(Symbol::LESS, [sym("y"), h.clone()]),
                    app(Symbol::LESS, [h.clone(), sym("z")]),
                ],
            )
        },
        |h| app(Symbol::AND, [h.clone(), sym("p")]),
        |h| app(Symbol::AND, [sym("p"), h.clone()]),
        |h| app(Symbol::OR, [h.clone(), sym("p")]),
        |h| app(Symbol::OR, [sym("p"), h.clone()]),
        |h| app(Symbol::NOT, [h.clone()]),
        |h| app(Symbol::RANGE, [h.clone(), sym("y")]),
        |h| app(Symbol::RANGE, [sym("y"), h.clone()]),
        |h| app(Symbol::LIST, [h.clone()]),
        |h| app(Symbol::FINITE_SET, [h.clone(), sym("y")]),
        |h| app(Symbol::AT, [h.clone(), int(0)]),
        |h| app(Symbol::AT, [sym("v"), h.clone()]),
        |h| call("f", [h.clone()]),
        |h| call("f", [h.clone(), sym("y")]),
        |h| call("g", [h.clone()]),
        |h| apply(h, &[sym("y")]),
        |h| apply(h, &[]),
        |h| app(Symbol::SIN, [h.clone()]),
        |h| apply(&Symbol::PI.into(), [h.clone()]),
        |h| app(Symbol::REPLACE_ALL, [sym("x"), h.clone()]),
        |h| app(Symbol::REPLACE_ALL, [h.clone(), sym("x")]),
        |h| lambda(&["x"], h),
        |h| app(Symbol::IF, [h.clone(), sym("y")]),
        |h| app(Symbol::IF, [sym("c"), h.clone()]),
        |h| app(Symbol::IF, [sym("c"), sym("y"), h.clone()]),
    ];
    let mut count = 0;
    for context in &contexts {
        for filler in &fillers {
            round_trip(&context(filler));
            count += 1;
        }
    }
    assert!(count > 2000, "{count}");
}

/// Every form the printer writes, with simple operands.
fn fillers(a: &Expr, b: &Expr, c: &Expr, f: &Expr) -> Vec<Expr> {
    let [a, b, c, f] = [a, b, c, f].map(Clone::clone);
    vec![
        a.clone(),
        sym("$E"),
        int(0),
        int(2),
        int(-2),
        Expr::from(Integer::from(10).pow(25)),
        rat(1, 2),
        rat(-1, 2),
        cx((0, 1), (1, 1)),
        cx((0, 1), (-1, 1)),
        cx((0, 1), (2, 1)),
        cx((0, 1), (-1, 2)),
        cx((2, 1), (3, 1)),
        cx((-2, 1), (-3, 1)),
        Expr::string("s"),
        Symbol::PI.into(),
        Symbol::TRUE.into(),
        Symbol::COMPLEX_INFINITY.into(),
        Symbol::INFINITY.into(),
        app(Symbol::DIRECTED_INFINITY, [int(-1)]),
        Symbol::SIN.into(),
        plus([a.clone(), b.clone()]),
        plus([a.clone(), times([int(-1), b.clone()])]),
        plus([times([int(-1), a.clone()]), b.clone()]),
        plus([a.clone(), int(-3)]),
        times([a.clone(), b.clone()]),
        times([int(-1), a.clone()]),
        times([int(-2), a.clone()]),
        times([int(-1), a.clone(), b.clone()]),
        times([a.clone(), pow(&b, &int(-1))]),
        times([rat(1, 2), a.clone()]),
        times([cx((0, 1), (-1, 1)), a.clone()]),
        pow(&a, &int(-1)),
        pow(&a, &int(2)),
        pow(&a, &rat(1, 2)),
        app(Symbol::ELEMENTWISE_POWER, [a.clone(), b.clone()]),
        app(Symbol::ELEMENTWISE_TIMES, [a.clone(), b.clone()]),
        app(Symbol::ELEMENTWISE_DIVIDE, [a.clone(), b.clone()]),
        app(Symbol::FACTORIAL, [a.clone()]),
        app(Symbol::EQUATION, [a.clone(), b.clone()]),
        app(Symbol::EQUAL, [a.clone(), b.clone()]),
        app(Symbol::LESS, [a.clone(), b.clone()]),
        app(Symbol::ELEMENT, [a.clone(), b.clone()]),
        app(
            Symbol::AND,
            &[
                app(Symbol::LESS, [a.clone(), b.clone()]),
                app(Symbol::LESS, [b.clone(), c.clone()]),
            ],
        ),
        app(
            Symbol::AND,
            &[
                app(Symbol::LESS, [a.clone(), b.clone()]),
                app(Symbol::LESS, [c.clone(), a.clone()]),
            ],
        ),
        app(Symbol::AND, [a.clone(), b.clone()]),
        app(Symbol::OR, [a.clone(), b.clone()]),
        app(Symbol::NOT, [a.clone()]),
        app(Symbol::RANGE, [a.clone(), b.clone()]),
        app(Symbol::LIST, [a.clone(), b.clone()]),
        app(Symbol::FINITE_SET, [a.clone()]),
        app(Symbol::AT, [a.clone(), int(0)]),
        call("f", [a.clone()]),
        call("f", []),
        call("f", [f.clone()]),
        apply(&call("g", [a.clone()]), [b.clone()]),
        app(Symbol::SIN, [a.clone()]),
        apply(&Symbol::PI.into(), [a.clone()]),
        lambda(&["a"], &a),
        app(Symbol::IF, [a.clone(), b.clone(), c.clone()]),
        app(Symbol::IF, [a.clone(), b.clone()]),
        app(
            Symbol::REPLACE_ALL,
            [a.clone(), app(Symbol::EQUATION, [b.clone(), int(1)])],
        ),
    ]
}

/// Atoms of every kind the printer writes.
fn leaf() -> impl Strategy<Value = Expr> {
    prop_oneof![
        prop::sample::select(vec!["x", "y", "f", "g", "n", "$E", "θ", "x_1"]).prop_map(sym),
        (-20i64..20).prop_map(int),
        (-20i64..20, 2i64..9).prop_map(|(n, d)| rat(n, d)),
        ((-5i64..5, 1i64..4), (-5i64..5, 1i64..4))
            .prop_filter("not real", |(_, im)| im.0 != 0)
            .prop_map(|(re, im)| cx(re, im)),
        prop::sample::select(vec!["", "s", "a\"b", "\\", "\n\u{7}"]).prop_map(Expr::string),
        prop::sample::select(vec![
            Symbol::PI,
            Symbol::E,
            Symbol::TRUE,
            Symbol::INDETERMINATE,
            Symbol::INFINITY,
            Symbol::COMPLEX_INFINITY,
            Symbol::SIN,
            Symbol::LOG,
            Symbol::FACTORIAL,
            Symbol::REPLACE_ALL,
        ])
        .prop_map(Expr::from),
        Just(app(Symbol::DIRECTED_INFINITY, [int(-1)])),
    ]
}

/// A product as evaluation leaves it: a coefficient first, and factors with
/// negative exponents.
fn evaluated_product(inner: BoxedStrategy<Expr>) -> impl Strategy<Value = Expr> {
    (
        (-6i64..6, 1i64..4),
        prop::collection::vec(inner.clone(), 1..3),
        prop::collection::vec((inner, -3i64..-1), 0..2),
    )
        .prop_map(|((n, d), factors, inverses)| {
            let mut all = vec![rat(n, d)];
            all.extend(factors);
            all.extend(inverses.iter().map(|(b, k)| pow(b, &int(*k))));
            times(all)
        })
}

/// `And` of comparisons: chains that share their operands, and comparisons
/// that need not, which are a chain only if they do.
fn comparisons(inner: BoxedStrategy<Expr>) -> impl Strategy<Value = Expr> {
    let ops = prop::sample::select(vec![
        Symbol::LESS,
        Symbol::LESS_EQUAL,
        Symbol::EQUAL,
        Symbol::ELEMENT,
    ]);
    prop_oneof![
        (
            prop::collection::vec(inner.clone(), 3..5),
            prop::collection::vec(ops, 4),
        )
            .prop_map(|(operands, ops)| {
                let tests: Vec<Expr> = operands
                    .windows(2)
                    .zip(ops)
                    .map(|(w, op)| app(op, [w[0].clone(), w[1].clone()]))
                    .collect();
                app(Symbol::AND, tests)
            }),
        prop::collection::vec((inner.clone(), inner), 2..4).prop_map(|pairs| {
            let tests: Vec<Expr> = pairs
                .into_iter()
                .map(|(a, b)| app(Symbol::LESS, [a, b]))
                .collect();
            app(Symbol::AND, tests)
        }),
    ]
}

/// `If` with two or three arguments, some of them `Null`.
fn conditional(inner: BoxedStrategy<Expr>) -> impl Strategy<Value = Expr> {
    let branch = prop_oneof![inner.clone(), Just(Expr::from(Symbol::NULL))];
    (inner, branch.clone(), prop::option::of(branch)).prop_map(|(c, a, b)| {
        let mut args = vec![c, a];
        args.extend(b);
        app(Symbol::IF, args)
    })
}

/// Random expressions over every form the printer writes.
fn expr() -> impl Strategy<Value = Expr> {
    leaf().prop_recursive(5, 48, 4, |inner| {
        let args = |n: std::ops::Range<usize>| prop::collection::vec(inner.clone(), n);
        let pair = (inner.clone(), inner.clone());
        let binary = prop::sample::select(vec![
            Symbol::ELEMENTWISE_POWER,
            Symbol::ELEMENTWISE_DIVIDE,
            Symbol::EQUATION,
            Symbol::EQUAL,
            Symbol::UNEQUAL,
            Symbol::LESS,
            Symbol::GREATER_EQUAL,
            Symbol::ELEMENT,
            Symbol::RANGE,
        ]);
        let variadic = prop::sample::select(vec![
            Symbol::ELEMENTWISE_TIMES,
            Symbol::AND,
            Symbol::OR,
            Symbol::FINITE_SET,
            Symbol::LIST,
        ]);
        let unary = prop::sample::select(vec![
            Symbol::NOT,
            Symbol::FACTORIAL,
            Symbol::SIN,
            Symbol::HOLD,
        ]);
        let params = prop::sample::select(vec![&["x"][..], &["x", "y"][..], &[][..], &["f"][..]]);
        prop_oneof![
            args(2..4).prop_map(plus),
            args(2..4).prop_map(times),
            evaluated_product(inner.clone()),
            pair.clone().prop_map(|(b, n)| pow(&b, &n)),
            (inner.clone(), -3i64..0).prop_map(|(b, k)| pow(&b, &int(k))),
            (binary, pair.clone()).prop_map(|(s, (a, b))| app(s, [a, b])),
            (variadic, args(2..4)).prop_map(|(s, a)| app(s, a)),
            comparisons(inner.clone()),
            (unary, inner.clone()).prop_map(|(s, a)| app(s, [a])),
            args(0..3).prop_map(|a| app(Symbol::LIST, a)),
            (inner.clone(), args(0..3)).prop_map(|(base, mut i)| {
                i.insert(0, base);
                app(Symbol::AT, i)
            }),
            // Calls, with heads that may occur in their arguments.
            (prop::sample::select(vec!["f", "g", "x"]), args(0..3)).prop_map(|(f, a)| call(f, a)),
            (inner.clone(), args(0..3)).prop_map(|(head, a)| apply(&head, a)),
            (params, inner.clone()).prop_map(|(params, body)| lambda(params, &body)),
            conditional(inner.clone()),
            pair.prop_map(|(e, r)| app(Symbol::REPLACE_ALL, [e, r])),
        ]
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2000))]

    /// Any expression of these forms prints, and its text reads back as it.
    #[test]
    fn round_trips(e in expr()) {
        round_trip(&e);
    }

    /// Text built by lowering prints back to text that lowers to the same
    /// structure, and printing that is stable.
    #[test]
    fn printing_lowered_text_is_stable(e in expr()) {
        let text = printed(&e);
        let back = read_with(&text, Names::english());
        let again = printed(&back);
        prop_assert_eq!(printed(&read_with(&again, Names::english())), again);
    }
}

/// Decimal literals, as lowering makes them (CAS-148), print back as
/// themselves.
#[test]
fn decimal_literals_print_back() {
    check_text(&[
        "0.1",
        "1.5",
        "0.0015",
        "2.50",
        "2e10",
        "5e0",
        "1.5e-30",
        "-0.5",
        "x - 0.5",
        "-0.5*x",
        "x - 0.5*y",
        "2*(-0.5)",
        "1.5 + 2.5*i",
        "[0.25, 1e-7]",
    ]);
}

/// A product led by a negative decimal, subtracted in a sum, reads back as
/// `-1` times the product: the same up to arithmetic (the `lang_print` fuzz
/// target, CAS-148).
#[test]
fn a_subtracted_decimal_product_reads_back() {
    let decimal = |digits: i64, exponent: i64| app(Symbol::DECIMAL, [int(digits), int(exponent)]);
    let e = plus([decimal(-9, 4), times([decimal(-9, 3), sym("y")])]);
    assert_eq!(round_trip(&e), "-9e4 - 9e3*y");
    assert_eq!(round_trip(&times([sym("y"), decimal(-5, -1)])), "y*(-0.5)");
}

fn ball_expr(b: cas_num::ball::Ball) -> Expr {
    Expr::from(b)
}

fn decimal_ball(s: &str) -> cas_num::ball::Ball {
    cas_num::ball::Ball::from_decimal(s, 53).unwrap()
}

/// In the input form proper, a ball is the literal that reads back to its
/// midpoint, or `ball(c, r)` when it is wider than that literal (CAS-148).
#[test]
fn balls_print_as_literals() {
    use cas_num::ball::{Ball, ComplexBall, Float};

    let third = Ball::from_rational(&ratio(1, 3), 53);
    let sum = decimal_ball("0.1").add(&decimal_ball("0.2"), 53);
    let tiny = Ball::exact(Float::from(&Integer::from(1)).mul_2exp(&Integer::from(-60)));
    let cases = [
        (ball_expr(decimal_ball("0.1")), "0.1"),
        (ball_expr(decimal_ball("-2.5")), "-2.5"),
        (ball_expr(third), "0.3333333333333333"),
        (ball_expr(tiny), "8.673617379884035e-19"),
        (ball_expr(decimal_ball("6.02e23")), "6.02e+23"),
        (ball_expr(decimal_ball("2.0")), "2.0"),
        (ball_expr(sum), "ball(0.3, 7.64e-17)"),
        (
            Expr::from(ComplexBall::new(decimal_ball("1.5"), decimal_ball("-2.5"))),
            "1.5 - 2.5*i",
        ),
        (times([ball_expr(decimal_ball("-0.5")), sym("x")]), "-0.5*x"),
        (plus([sym("x"), ball_expr(decimal_ball("-0.5"))]), "x - 0.5"),
    ];
    let mut wrong = Vec::new();
    for (e, expected) in cases {
        let found = printed(&e);
        if found != expected {
            wrong.push(format!("{e:?}: found {found}, expected {expected}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Shown for people, a ball is the decimal it shows (D2).
#[test]
fn balls_print_as_shown() {
    use super::{Balls, print_with};

    let shown = |e: &Expr| {
        print_with(
            e,
            Names::english(),
            Balls::Shown {
                working_precision: 53,
            },
        )
        .unwrap()
    };
    let sum = decimal_ball("0.1").add(&decimal_ball("0.2"), 53);
    assert_eq!(shown(&ball_expr(sum.clone())), "0.3");
    let zero = sum.sub(&decimal_ball("0.3"), 53);
    assert!(shown(&ball_expr(zero)).starts_with("ball(0.0, "));
    assert_eq!(
        shown(&ball_expr(cas_num::ball::Ball::from_rational(
            &ratio(1, 3),
            53
        ))),
        "0.333333333333333"
    );
}

/// A ball, its midpoint dyadic, its radius 0 or below its midpoint.
fn ball_strategy() -> impl Strategy<Value = cas_num::ball::Ball> {
    use cas_num::ball::{Ball, Float, Mag};

    (
        any::<i64>(),
        -300i64..300,
        prop_oneof![Just(None), (1i64..120).prop_map(Some)],
    )
        .prop_map(|(m, e, j)| {
            let mid = Float::from(&Integer::from(m)).mul_2exp(&Integer::from(e));
            let rad = match j {
                Some(j) if !mid.is_zero() => {
                    Mag::power_of_two(&(mid.exponent() - &Integer::from(j)))
                }
                _ => Mag::ZERO,
            };
            Ball::new(mid, rad)
        })
}

proptest! {
    /// A ball printed in the input form reads back, through the session's
    /// literals, as a ball of the same midpoint containing it.
    #[test]
    fn balls_read_back(b in ball_strategy(), working in prop::sample::select(vec![24u64, 53, 100, 200])) {
        use super::{Balls, print_with};
        use cas_core::{Context, View};

        let text = print_with(&ball_expr(b.clone()), Names::english(), Balls::Exact { working_precision: working })
            .unwrap_or_else(|e| panic!("{b}: {e}"));
        let mut ctx = Context::new();
        ctx.set_working_precision(working);
        let lowered = crate::literal::decimal_balls(&read_with(&text, Names::english()), working);
        let back = ctx.eval(&lowered).unwrap();
        let View::Ball(back) = back.view() else {
            panic!("{b} printed as {text}, reads back as {back:?}");
        };
        prop_assert_eq!(back.mid(), b.mid(), "{} printed as {}", b, text);
        prop_assert!(b.rad() <= back.rad(), "{} printed as {}", b, text);
    }
}

/// The set operators print infix, with parentheses where the reading
/// would group otherwise (CAS-DOC-32 §3).
#[test]
fn set_operators() {
    let (a, b, c) = (sym("a"), sym("b"), sym("c"));
    let call = |s: Symbol, args: &[Expr]| Expr::call(s, args.to_vec());
    check(&[
        (call(Symbol::UNION, &[a.clone(), b.clone()]), "a union b"),
        (
            call(Symbol::INTERSECTION, &[a.clone(), b.clone(), c.clone()]),
            "a intersect b intersect c",
        ),
        (
            call(Symbol::COMPLEMENT, &[a.clone(), b.clone(), c.clone()]),
            "a minus b minus c",
        ),
        (
            call(
                Symbol::UNION,
                &[a.clone(), call(Symbol::COMPLEMENT, &[b.clone(), c.clone()])],
            ),
            "a union (b minus c)",
        ),
        (
            call(
                Symbol::COMPLEMENT,
                &[call(Symbol::UNION, &[a.clone(), b.clone()]), c.clone()],
            ),
            "(a union b) minus c",
        ),
        (
            call(
                Symbol::UNION,
                &[
                    a.clone(),
                    call(Symbol::INTERSECTION, &[b.clone(), c.clone()]),
                ],
            ),
            "a union b intersect c",
        ),
        (
            call(
                Symbol::ELEMENT,
                &[sym("x"), call(Symbol::UNION, &[a.clone(), b.clone()])],
            ),
            "x in a union b",
        ),
    ]);
}
