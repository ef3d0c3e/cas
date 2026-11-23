//! Tests of the variation and sign tables (CAS-582): giac's tables, its
//! `undef`, and what stays.

use crate::{Context, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// A row: a cell starting with `'` is that string, any other full form.
fn row(cells: &[&str]) -> Expr {
    Expr::call(
        Symbol::LIST,
        cells
            .iter()
            .map(|c| c.strip_prefix('\'').map_or_else(|| read(c), Expr::string))
            .collect::<Vec<_>>(),
    )
}

fn table(rows: &[&[&str]]) -> Expr {
    Expr::call(
        Symbol::LIST,
        rows.iter().map(|r| row(r)).collect::<Vec<_>>(),
    )
}

/// `input` in full form, evaluated: the table `expected`, evaluated.
fn check(input: &str, expected: &Expr) {
    let ctx = Context::new();
    let found = ctx.eval(&read(input)).unwrap();
    let expected = ctx.eval(expected).unwrap();
    assert_eq!(found, expected, "{input}");
}

/// giac's `tabsign(x^2-x,x)`.
#[test]
fn sign_table() {
    check(
        "TabSign[Plus[Power[x, 2], Times[-1, x]], x]",
        &table(&[
            &[
                "x",
                "Times[-1, Infinity]",
                "' ",
                "0",
                "' ",
                "1",
                "' ",
                "Infinity",
            ],
            &[
                "Plus[Power[x, 2], Times[-1, x]]",
                "Infinity",
                "'+",
                "0",
                "'-",
                "0",
                "'+",
                "Infinity",
            ],
        ]),
    );
}

/// giac's `tabsign((x-1)/(x+2),x)`: the two columns at the pole merge.
#[test]
fn sign_table_pole() {
    check(
        "TabSign[Times[Plus[x, -1], Power[Plus[x, 2], -1]], x]",
        &table(&[
            &[
                "x",
                "Times[-1, Infinity]",
                "' ",
                "-2",
                "' ",
                "1",
                "' ",
                "Infinity",
            ],
            &[
                "Times[Plus[x, -1], Power[Plus[x, 2], -1]]",
                "1",
                "'+",
                "'||",
                "'-",
                "0",
                "'+",
                "1",
            ],
        ]),
    );
}

/// giac's `tabvar(x^2,x,-3,5)`, the four rows.
#[test]
fn variation_table() {
    check(
        "TabVar[Power[x, 2], x, -3, 5]",
        &table(&[
            &["x", "-3", "' ", "0", "' ", "5"],
            &["Equation[\"y'\", Times[2, x]]", "-6", "'-", "0", "'+", "10"],
            &["Equation[y, Power[x, 2]]", "9", "'↓", "0", "'↑", "25"],
            &["\"(y')'\"", "2", "'+ (∪)", "2", "'+ (∪)", "2"],
        ]),
    );
}

/// `diff = 1`: no row of `f''`; `x = a..b`.
#[test]
fn options() {
    check(
        "TabVar[Power[x, 2], Equation[x, Range[-3, 5]], Equation[D, 1]]",
        &table(&[
            &["x", "-3", "' ", "0", "' ", "5"],
            &["Equation[\"y'\", Times[2, x]]", "-6", "'-", "0", "'+", "10"],
            &["Equation[y, Power[x, 2]]", "9", "'↓", "0", "'↑", "25"],
        ]),
    );
}

/// giac's `undef`: infinitely many critical points.
#[test]
fn undefined() {
    let ctx = Context::new();
    for input in ["TabVar[Plus[x, Sin[x]], x]", "TabVar[Abs[Plus[x, -1]], x]"] {
        assert_eq!(
            ctx.eval(&read(input)).unwrap(),
            Expr::from(Symbol::INDETERMINATE),
            "{input}"
        );
    }
}

/// Parametric tables, parameters and constants stay.
#[test]
fn stays() {
    let ctx = Context::new();
    for input in [
        "TabVar[List[Sin[t], Cos[t]], t]",
        "TabVar[Times[a, Power[x, 2]], x]",
        "TabSign[3, x]",
    ] {
        let e = read(input);
        assert_eq!(ctx.eval(&e).unwrap(), e, "{input}");
    }
}
