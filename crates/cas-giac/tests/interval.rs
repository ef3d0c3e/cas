//! `convert(…,interval)` against giac (CAS-592): giac's intervals
//! `[lo..hi]`, ours balls. Each of our balls overlaps giac's interval at
//! the same place, widened by a unit of its last printed digit, and the
//! rest, each interval and ball a symbol, is the same expression.

use cas_core::{Context, Expr, Rational, View};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// A decimal giac prints, `-1.25`, and a unit of its last digit.
fn decimal(text: &str) -> Option<(Rational, Rational)> {
    let (int, frac) = text.split_once('.').unwrap_or((text, ""));
    let den = format!("1{}", "0".repeat(frac.len()));
    let q: Rational = format!("{int}{frac}/{den}").parse().ok()?;
    let unit: Rational = format!("1/{den}").parse().ok()?;
    Some((q, unit))
}

/// giac's text with each interval `[lo..hi]` the symbol `iv<k>`, and the
/// intervals, widened by a unit of the last printed digit.
fn intervals(text: &str) -> Option<(String, Vec<(Rational, Rational)>)> {
    let mut out = String::new();
    let mut found = Vec::new();
    let mut rest = text;
    while let Some(dots) = rest.find("..") {
        let open = rest[..dots].rfind('[')?;
        let close = dots + rest[dots..].find(']')?;
        let (lo, lo_unit) = decimal(&rest[open + 1..dots])?;
        let (hi, hi_unit) = decimal(&rest[dots + 2..close])?;
        out.push_str(&rest[..open]);
        out += "iv";
        out += &found.len().to_string();
        found.push((&lo - &lo_unit, &hi + &hi_unit));
        rest = &rest[close + 1..];
    }
    out.push_str(rest);
    Some((out, found))
}

/// giac's answer to `input` and ours: overlapping intervals and balls in
/// the same places.
fn agrees(giac: &Giac, input: &str) -> Result<(), String> {
    let text = giac.eval(input).map_err(|e| e.to_string())?;
    let (shape, theirs) = intervals(&text).ok_or_else(|| format!("{input}: giac's {text}"))?;
    let shape = ours(&shape).map_err(|e| e.to_string())?;
    let value = ours(input).map_err(|e| e.to_string())?;
    let differs = |why: &str| format!("{input}: {why}: ours {value:?}, giac's {text}");
    let balls: Vec<Expr> = value
        .pre_order()
        .filter(|e| matches!(e.view(), View::Ball(_)))
        .cloned()
        .collect();
    if balls.len() != theirs.len() {
        return Err(differs("not as many numbers"));
    }
    let mut rules = Vec::new();
    for (k, (ball, (lo, hi))) in balls.iter().zip(&theirs).enumerate() {
        let View::Ball(b) = ball.view() else {
            unreachable!()
        };
        let mid = b.mid().to_rational();
        let rad = b
            .rad()
            .to_float()
            .ok_or_else(|| differs("an infinite ball"))?
            .to_rational();
        if &(&mid + &rad) < lo || &(&mid - &rad) > hi {
            return Err(differs(&format!("number {k} outside giac's interval")));
        }
        rules.push((ball.clone(), Expr::symbol(&format!("iv{k}"))));
    }
    let ctx = Context::new();
    let ours_shape = ctx
        .eval(&value.substitute(&rules))
        .map_err(|e| e.to_string())?;
    if ours_shape != shape {
        return Err(differs("not the same expression"));
    }
    Ok(())
}

fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| agrees(giac, input).err())
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// giac's `convert(…,interval)` with and without digits; lists and the
/// numbers of an expression.
#[test]
fn intervals_overlap() {
    check(&[
        "convert(sqrt(2),interval)",
        "convert(sqrt(2),interval,20)",
        "convert(exp(1),interval,30)",
        "convert(1/3,interval)",
        "convert([1/3,pi],interval)",
        "convert(x+1,interval)",
        "convert(sqrt(2)*sin(2)+sqrt(2)*cos(2),interval,20)",
    ]);
}
