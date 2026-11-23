//! `euler_lagrange` under giac's `assume` (CAS-575): its example with
//! `assume(y>=0)`. giac's assumptions are global to its session, so the
//! example runs alone in this test binary, the `assume` and the `purge` in
//! one input. Ours has no assumption: its first integral keeps
//! `√((1 + y′²)/y)·y`, giac's `√((1 + y′²) y)`, equal for `y > 0`; so the
//! two answers are compared at a point with `y > 0`.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};

/// `e` at `y(t) = 2/3`, `y′(t) = 5/7`, `y″(t) = 3/11` and `K_0 = 0`,
/// numerically.
fn at_point(e: &Expr) -> Expr {
    let read = |text: &str| Expr::parse_full_form(text).unwrap();
    let values = [
        (read("Derivative[2][y][t]"), read("3/11")),
        (read("Derivative[1][y][t]"), read("5/7")),
        (read("y[t]"), read("2/3")),
        (Expr::symbol("$K_0"), Expr::from(0)),
    ];
    let e = e.replace(|part| {
        values
            .iter()
            .find_map(|(from, to)| (part == from).then(|| to.clone()))
    });
    Expr::call(Symbol::N, [e])
}

/// Whether `|e| < 10⁻¹⁰`.
fn small(ctx: &Context, e: Expr) -> bool {
    let bound = Expr::parse_full_form("Power[10, -10]").unwrap();
    let less = Expr::call(Symbol::LESS, [Expr::call(Symbol::ABS, [e]), bound]);
    ctx.eval(&less).ok() == Some(Expr::from(Symbol::TRUE))
}

#[test]
fn under_an_assumption() {
    let Some(giac) = Giac::shared() else {
        eprintln!("libgiac is missing: comparisons skipped");
        return;
    };
    let input = "[assume(y>=0),euler_lagrange(sqrt((1+y'^2)/y),t,y),purge(y)][1]";
    // giac's second derivative diff(y(t),t,2) as ours; its one constant,
    // numbered by a counter global to giac's session, as K_0.
    let text = giac.eval(input).unwrap().replace(",t,2)", ",t,t)");
    let at = text.find("=K_").expect("giac's constant") + 3;
    let digits = text[at..].find(|c: char| !c.is_ascii_digit()).unwrap_or(0);
    let text = format!("{}0{}", &text[..at], &text[at + digits..]);
    let theirs = ours(&text).unwrap();
    let found = ours("euler_lagrange(sqrt((1+diff(y(t),t)^2)/y(t)),t,y)").unwrap();
    let ctx = Context::new();
    // l − r of each equation, at the point.
    let sides = |e: &Expr| -> Vec<Expr> {
        e.args()
            .iter()
            .map(|eq| {
                let [l, r] = eq.args() else {
                    panic!("{eq:?} is not an equation");
                };
                let minus_r = Expr::call(Symbol::TIMES, [Expr::from(-1), r.clone()]);
                ctx.eval(&at_point(&Expr::call(Symbol::PLUS, [l.clone(), minus_r])))
                    .unwrap()
            })
            .collect()
    };
    let (ours_sides, giac_sides) = (sides(&found), sides(&theirs));
    assert!(text.contains("=K_0,"), "giac's {text}");
    assert_eq!(ours_sides.len(), 2, "ours {found:?}");
    assert_eq!(giac_sides.len(), 2, "giac's {text}");
    // The first integral as it is; the equation = 0 up to its sign, which
    // giac's radsimp flips.
    let difference = |a: &Expr, b: &Expr, sign: i64| {
        let b = Expr::call(Symbol::TIMES, [Expr::from(sign), b.clone()]);
        Expr::call(Symbol::PLUS, [a.clone(), b])
    };
    assert!(
        small(&ctx, difference(&ours_sides[0], &giac_sides[0], -1)),
        "first integral: ours {found:?}, giac's {text}"
    );
    assert!(
        small(&ctx, difference(&ours_sides[1], &giac_sides[1], -1))
            || small(&ctx, difference(&ours_sides[1], &giac_sides[1], 1)),
        "equation: ours {found:?}, giac's {text}"
    );
}
