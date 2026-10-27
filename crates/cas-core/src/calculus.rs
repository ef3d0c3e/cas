//! Calculus (CAS-DOC-4 §8, CAS-59): `D` written as system rules (§5.4), the
//! acceptance test of the rewriting design.
//!
//! `D[f, x]` is the derivative of `f` in the symbol `x`. The rules, tried in
//! this order:
//!
//! - what is free of `x` has derivative 0, which the Bloom filter of `f`'s
//!   symbols tells at once (`FreeQ`), and `x` has 1;
//! - `Integrate[f, x]`, an antiderivative that stays, has `f` (CAS-DOC-9);
//! - lists thread, and sums are linear: `D` is mapped over all the terms
//!   at once, not split off one at a time, which would rebuild the rest of
//!   the sum at each;
//! - the product rule, `D[a_ b__, x_]`, through AC matching;
//! - powers: `f^n` with `n` free of `x`, `E^g`, and `f^g` in general;
//! - the chain rule for `Sin`, `Cos`, `Tan`, `ArcSin`, `ArcCos`, `ArcTan`,
//!   the hyperbolic `Sinh` … `ArcTanh` (CAS-DOC-15 §2), `Erf` and `Erfc`
//!   (§8), `Exp` and `Log`, from a table, and `Abs`, whose derivative is
//!   `Abs[f] / f` times `f′`, giac's `abs(x)/x`, so that the logarithms of
//!   absolute values integrals give differentiate back (CAS-DOC-9);
//! - a function the table does not know: `Derivative[1][g][u] D[u, x]`, and
//!   `Derivative[n + 1][g][u]` for the derivative of `Derivative[n][g][u]`;
//! - higher derivatives, `D[f, {x, n}]` for an integer `n >= 0`, and partial
//!   ones, `D[f, x, y, …]`.
//!
//! A variable that is not a symbol leaves `D` as it is. The results go
//! through the evaluator, so `D[x^2 + 3 x, x]` is `3 + 2 x` with no
//! simplification asked for. Each rule has a stable name, `"D/product"`,
//! that a trace shows with the line it is written on (§9).

use crate::rules::SystemRule;

/// The variable of a derivative: a symbol.
const X: &str = "PatternTest[Pattern[x, Blank[]], SymbolQ]";

/// The rule `D[f, x] :> rhs` for a symbol `x`, written at the caller's line.
#[track_caller]
fn d(name: &'static str, f: &str, rhs: &str) -> SystemRule {
    SystemRule::new(name, &format!("D[{f}, {X}]"), rhs)
}

/// The rules of `D`, in the order they are tried.
#[allow(clippy::too_many_lines, reason = "a table, one rule after another")]
pub(crate) fn rules() -> Vec<SystemRule> {
    vec![
        d(
            "D/constant",
            "Pattern[c, Blank[]]",
            "Condition[0, FreeQ[c, x]]",
        ),
        d("D/variable", "Pattern[x, Blank[]]", "1"),
        // The derivative of an antiderivative that stays (CAS-DOC-9).
        d(
            "D/Integrate",
            "Integrate[Pattern[f, Blank[]], Pattern[x, Blank[]]]",
            "f",
        ),
        d(
            "D/list",
            "Pattern[l, Blank[List]]",
            "Map[Function[D[Slot[1], x]], l]",
        ),
        d(
            "D/sum",
            "Plus[Pattern[a, Blank[]], Pattern[b, BlankSequence[]]]",
            "Map[Function[D[Slot[1], x]], Plus[a, b]]",
        ),
        d(
            "D/product",
            "Times[Pattern[a, Blank[]], Pattern[b, BlankSequence[]]]",
            "Plus[Times[D[a, x], b], Times[a, D[Times[b], x]]]",
        ),
        d(
            "D/power",
            "Power[Pattern[f, Blank[]], Pattern[n, Blank[]]]",
            "Condition[Times[n, Power[f, Plus[n, -1]], D[f, x]], FreeQ[n, x]]",
        ),
        d(
            "D/E",
            "Power[E, Pattern[g, Blank[]]]",
            "Times[Power[E, g], D[g, x]]",
        ),
        d(
            "D/exponential",
            "Power[Pattern[f, Blank[]], Pattern[g, Blank[]]]",
            "Times[Power[f, g], Plus[Times[D[g, x], Log[f]], Times[g, D[f, x], Power[f, -1]]]]",
        ),
        d(
            "D/Sin",
            "Sin[Pattern[f, Blank[]]]",
            "Times[Cos[f], D[f, x]]",
        ),
        d(
            "D/Cos",
            "Cos[Pattern[f, Blank[]]]",
            "Times[-1, Sin[f], D[f, x]]",
        ),
        d(
            "D/Tan",
            "Tan[Pattern[f, Blank[]]]",
            "Times[Power[Cos[f], -2], D[f, x]]",
        ),
        d(
            "D/ArcSin",
            "ArcSin[Pattern[f, Blank[]]]",
            "Times[Power[Plus[1, Times[-1, Power[f, 2]]], -1/2], D[f, x]]",
        ),
        d(
            "D/ArcCos",
            "ArcCos[Pattern[f, Blank[]]]",
            "Times[-1, Power[Plus[1, Times[-1, Power[f, 2]]], -1/2], D[f, x]]",
        ),
        d(
            "D/ArcTan",
            "ArcTan[Pattern[f, Blank[]]]",
            "Times[Power[Plus[1, Power[f, 2]], -1], D[f, x]]",
        ),
        // The hyperbolic functions (CAS-DOC-15 §2).
        d(
            "D/Sinh",
            "Sinh[Pattern[f, Blank[]]]",
            "Times[Cosh[f], D[f, x]]",
        ),
        d(
            "D/Cosh",
            "Cosh[Pattern[f, Blank[]]]",
            "Times[Sinh[f], D[f, x]]",
        ),
        d(
            "D/Tanh",
            "Tanh[Pattern[f, Blank[]]]",
            "Times[Power[Cosh[f], -2], D[f, x]]",
        ),
        d(
            "D/ArcSinh",
            "ArcSinh[Pattern[f, Blank[]]]",
            "Times[Power[Plus[1, Power[f, 2]], -1/2], D[f, x]]",
        ),
        d(
            "D/ArcCosh",
            "ArcCosh[Pattern[f, Blank[]]]",
            "Times[Power[Plus[-1, f], -1/2], Power[Plus[1, f], -1/2], D[f, x]]",
        ),
        d(
            "D/ArcTanh",
            "ArcTanh[Pattern[f, Blank[]]]",
            "Times[Power[Plus[1, Times[-1, Power[f, 2]]], -1], D[f, x]]",
        ),
        // The error functions (CAS-DOC-15 §8).
        d(
            "D/Erf",
            "Erf[Pattern[f, Blank[]]]",
            "Times[2, Power[Pi, -1/2], Power[E, Times[-1, Power[f, 2]]], D[f, x]]",
        ),
        d(
            "D/Erfc",
            "Erfc[Pattern[f, Blank[]]]",
            "Times[-2, Power[Pi, -1/2], Power[E, Times[-1, Power[f, 2]]], D[f, x]]",
        ),
        d(
            "D/Exp",
            "Exp[Pattern[f, Blank[]]]",
            "Times[Exp[f], D[f, x]]",
        ),
        d(
            "D/Log",
            "Log[Pattern[f, Blank[]]]",
            "Times[Power[f, -1], D[f, x]]",
        ),
        // abs(f)/f times f', as giac writes it (CAS-DOC-9).
        d(
            "D/Abs",
            "Abs[Pattern[f, Blank[]]]",
            "Times[Abs[f], Power[f, -1], D[f, x]]",
        ),
        d(
            "D/derivative",
            "Derivative[Pattern[n, Blank[]]][Pattern[g, Blank[]]][Pattern[u, Blank[]]]",
            "Times[Derivative[Plus[n, 1]][g][u], D[u, x]]",
        ),
        d(
            "D/function",
            "Pattern[g, Blank[]][Pattern[u, Blank[]]]",
            "Times[Derivative[1][g][u], D[u, x]]",
        ),
        SystemRule::new(
            "D/order",
            &format!("D[Pattern[f, Blank[]], List[{X}, 0]]"),
            "f",
        ),
        SystemRule::new(
            "D/higher",
            &format!(
                "D[Pattern[f, Blank[]], List[{X}, PatternTest[Pattern[n, Blank[]], IntegerQ]]]"
            ),
            "Condition[D[D[f, x], List[x, Plus[n, -1]]], Greater[n, 0]]",
        ),
        SystemRule::new(
            "D/partial",
            "D[Pattern[f, Blank[]], Pattern[x, Blank[]], Pattern[y, BlankSequence[]]]",
            "D[D[f, x], y]",
        ),
    ]
}

#[cfg(test)]
mod tests;
