//! `export_mathml` against giac (CAS-616): an expression's answer in
//! `MathML`, ours printed by the port of giac's `gen2markup`
//! (`cas_lang::markup`), in giac's three forms: `export_mathml(e, display)`
//! (presentation), `export_mathml(e, content)` and `export_mathml(e)`
//! (both, in `<semantics>`). Rows: giac's answer, its quotes removed, is
//! ours, whitespace collapsed (C0047); in the default form the annotation
//! is each one's input form (giac `x+1`, ours `1 + x`), so it is left out
//! of the comparison and ours checked against our `string`. Departures:
//! the inputs of `latex.rs`'s, ours asserted in presentation and content;
//! giac's texts are in `cas-lang`'s `markup/giac_mathml.tsv`.
//!
//! giac's `mathml(e)` (`cas_lang::mathml`, a port of `mathml.cc`): the same
//! inputs, giac's text, its quotes removed, is ours exactly (a newline as
//! giac's table shows it, ' | ', on both sides, C0047); departures as
//! above, ours asserted (giac's texts in `markup/giac_mathml_doc.tsv`).

use cas_core::Expr;
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// The text of a string answer: giac's quotes removed and its doubled
/// inner quotes undone.
fn unquoted(answer: &str) -> Option<String> {
    let inner = answer.trim().strip_prefix('"')?.strip_suffix('"')?;
    Some(inner.replace("\"\"", "\""))
}

fn collapsed(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ")
}

const OPEN: &str = "<annotation encoding='giac'>";
const CLOSE: &str = "</annotation></semantics></math>";

/// A default-form text without its annotation, and the annotation.
fn split_annotation(text: &str) -> Option<(&str, &str)> {
    let (head, rest) = text.split_once(OPEN)?;
    Some((head, rest.strip_suffix(CLOSE)?))
}

/// giac's `str_to_mml`: trimmed, `&`, `<`, `>` escaped.
fn escaped(s: &str) -> String {
    s.trim_matches([' ', '\n', '\t'])
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// Ours of `input` as a string.
fn ours_text(input: &str) -> Option<String> {
    ours(input)
        .ok()
        .as_ref()
        .and_then(Expr::as_str)
        .map(str::to_owned)
}

/// Each input `x`: giac's `export_mathml(x, display)`,
/// `export_mathml(x, content)` and `export_mathml(x)` are ours, the last
/// one's annotation ours `string(x)`.
fn check_mathml(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for x in inputs {
        for kind in [",display", ",content", ""] {
            let input = format!("export_mathml({x}{kind})");
            let answer = giac.eval(&input).unwrap_or_default();
            let mine = ours_text(&input);
            let ok = match (unquoted(&answer), &mine) {
                (Some(theirs), Some(mine)) if kind.is_empty() => {
                    match (split_annotation(&theirs), split_annotation(mine)) {
                        (Some((theirs, _)), Some((mine, note))) => {
                            collapsed(theirs) == collapsed(mine)
                                && ours_text(&format!("string({x})"))
                                    .is_some_and(|s| escaped(&s) == note)
                        }
                        _ => false,
                    }
                }
                (Some(theirs), Some(mine)) => collapsed(&theirs) == collapsed(mine),
                _ => false,
            };
            if !ok {
                wrong.push(format!("{input}: giac {answer}, ours {mine:?}"));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The rows of `latex.rs`'s `answers_as_giac_writes_them` (giac's `MathML`
/// of each is in `giac_mathml.tsv`).
const ANSWERS: [&str; 59] = [
    "1/2",
    "-1/2",
    "x^2",
    "sqrt(x)",
    "sqrt(2)/2",
    "x/2",
    "2/x",
    "x*y",
    "2*x",
    "-x",
    "x-y",
    "-x-y",
    "1/x",
    "x^(-2)",
    "sin(x)",
    "cos(x)^2",
    "sin(x^2)",
    "exp(x)",
    "exp(2*x)",
    "ln(x)",
    "abs(x)",
    "floor(x)",
    "asin(x)",
    "pi",
    "infinity",
    "i",
    "1+2*i",
    "3.14",
    "1.5e-7",
    "-2.5",
    "3.0",
    "0.001",
    "1e-5",
    "1.23e20",
    "1e100",
    "12345678901.0",
    "1.0/3",
    "[1,2,3]",
    "[[1,2],[3,4]]",
    "x=1",
    "x^2+y^2=1",
    "f(x)",
    "f(x,y)",
    "diff(f(x),x)",
    "int(f(x),x)",
    r#""abc""#,
    "alpha",
    "abc",
    "x_1",
    "a_n",
    "n!",
    "x^y^z",
    "a/(b*c)",
    "2/3*x",
    "x*(-y)",
    "1..2",
    "true",
    "undef",
    "rootof([1,0],[1,0,1,1])",
];

#[test]
fn answers_as_giac_writes_them() {
    check_mathml(&ANSWERS);
}

/// The rows of `latex.rs`'s `more_of_giacs_table` but `diff(f(x),x,x)`,
/// which giac's `MathML` table lacks.
const MORE: [&str; 56] = [
    "a*b*c",
    "log(x)",
    "e",
    "cos(x)",
    "tan(x)",
    "atan(x)",
    "sinh(x)",
    "2^x",
    "theta^2",
    "-2*x",
    "3*(-2)",
    "123456789012345678901234567890",
    "0.5",
    "6.02e23",
    "1e-10",
    "-1.5e-7",
    "3*i",
    "-i",
    "2-3*i",
    "-3/4",
    "3/4*x",
    "x/(2*y)",
    "1/(2*x)",
    "-x/3",
    "y/x",
    "2/(3*x)",
    r#""a_b""#,
    r#""50%""#,
    "beta",
    "omega",
    "abc_1",
    "x12",
    "alpha1",
    "x_1_2",
    "[]",
    "[[1,-2],[x,y^2]]",
    "x^(1/3)",
    "sin(x)^3",
    "ln(x)^2",
    "exp(-x)",
    "2*exp(x)",
    "floor(x/2)",
    "sin(2*x)",
    "3*sin(x)",
    "f(x)/2",
    "(2*n)!",
    "int(f(x),x,0,1)",
    "sum(f(k),k,1,n)",
    "x>=1",
    "a..b",
    "false",
    "pi/2",
    "x*y*z^2",
    "e*x",
    "-x^2",
    "x^2*y^3",
];

#[test]
fn more_of_giacs_table() {
    check_mathml(&MORE);
}

/// Departures: ours in presentation and content, giac's answer in a
/// comment (its `MathML` is in `giac_mathml.tsv`): giac's evaluator
/// leaves another shape or order, its printer and ours write the same shape
/// alike.
#[test]
#[allow(clippy::too_many_lines)]
fn departures() {
    let cases: &[(&str, &str, &str)] = &[
        // Our sums put constants first: giac sqrt(x+1).
        (
            "sqrt(x+1)",
            r"<msqrt><mrow><mn>1</mn><mo>+</mo><mi>x</mi></mrow></msqrt>",
            r"<apply id='id4'><root/><apply id='id3'><plus/><cn id='id1' type='integer'>1</cn><ci id='id2'>x</ci></apply></apply>",
        ),
        // giac (x^4-1)/(x^2+3).
        (
            "(x^4-1)/(x^2+3)",
            r"<mfrac><mrow><mo>&minus;</mo><mn>1</mn><mo>+</mo><msup><mi>x</mi><mn>4</mn></msup></mrow><mrow><mn>3</mn><mo>+</mo><msup><mi>x</mi><mn>2</mn></msup></mrow></mfrac>",
            r"<apply id='id11'><divide/><apply id='id5'><plus/><apply><minus/><cn id='id1' type='integer'>1</cn></apply><apply id='id4'><power/><ci id='id2'>x</ci><cn id='id3' type='integer'>4</cn></apply></apply><apply id='id10'><plus/><cn id='id6' type='integer'>3</cn><apply id='id9'><power/><ci id='id7'>x</ci><cn id='id8' type='integer'>2</cn></apply></apply></apply>",
        ),
        // Ascending degree: giac 2*x^2-3*x+1.
        (
            "2*x^2-3*x+1",
            r"<mrow><mn>1</mn><mo>&minus;</mo><mrow><mn>3</mn><mo>&it;</mo><mi>x</mi></mrow><mo>+</mo><mrow><mn>2</mn><mo>&it;</mo><msup><mi>x</mi><mn>2</mn></msup></mrow></mrow>",
            r"<apply id='id10'><plus/><cn id='id1' type='integer'>1</cn><apply><minus/><apply id='id4'><times/><cn id='id2' type='integer'>3</cn><ci id='id3'>x</ci></apply></apply><apply id='id9'><times/><cn id='id5' type='integer'>2</cn><apply id='id8'><power/><ci id='id6'>x</ci><cn id='id7' type='integer'>2</cn></apply></apply></apply>",
        ),
        // giac 1/(x+1).
        (
            "1/(x+1)",
            r"<mfrac><mn>1</mn><mrow><mn>1</mn><mo>+</mo><mi>x</mi></mrow></mfrac>",
            r"<apply id='id4'><divide/><cn type='integer'>1</cn><apply id='id3'><plus/><cn id='id1' type='integer'>1</cn><ci id='id2'>x</ci></apply></apply>",
        ),
        // giac (x+1)^2.
        (
            "(x+1)^2",
            r"<msup><mfenced><mrow><mn>1</mn><mo>+</mo><mi>x</mi></mrow></mfenced><mn>2</mn></msup>",
            r"<apply id='id5'><power/><apply id='id3'><plus/><cn id='id1' type='integer'>1</cn><ci id='id2'>x</ci></apply><cn id='id4' type='integer'>2</cn></apply>",
        ),
        // Our factor order: giac (x+1)*(x-1).
        (
            "(x+1)*(x-1)",
            r"<mrow><mfenced><mrow><mo>&minus;</mo><mn>1</mn><mo>+</mo><mi>x</mi></mrow></mfenced><mo>&it;</mo><mfenced><mrow><mn>1</mn><mo>+</mo><mi>x</mi></mrow></mfenced></mrow>",
            r"<apply id='id7'><times/><apply id='id3'><plus/><apply><minus/><cn id='id1' type='integer'>1</cn></apply><ci id='id2'>x</ci></apply><apply id='id6'><plus/><cn id='id4' type='integer'>1</cn><ci id='id5'>x</ci></apply></apply>",
        ),
        // Our products put cos before sin: giac sin(x)*cos(x).
        (
            "sin(x)*cos(x)",
            r"<mrow><mrow><mi>cos</mi><mo>&af;</mo><mi>x</mi></mrow><mo>&middot;</mo><mrow><mi>sin</mi><mo>&af;</mo><mi>x</mi></mrow></mrow>",
            r"<apply id='id7'><times/><apply id='id3'><cos id='id1'/><ci id='id2'>x</ci></apply><apply id='id6'><sin id='id4'/><ci id='id5'>x</ci></apply></apply>",
        ),
        // Ours keeps Times[-1, Plus[1, x]] (C0029): giac -x-1.
        (
            "-(x+1)",
            r"<mrow><mo>&minus;</mo><mfenced><mrow><mn>1</mn><mo>+</mo><mi>x</mi></mrow></mfenced></mrow>",
            r"<apply id='id3'><minus/><apply><plus/><cn id='id1' type='integer'>1</cn><ci id='id2'>x</ci></apply></apply>",
        ),
        // giac's evaluator turns x<y into y>x (giac <mi>y</mi><mo>&gt;</mo><mi>x</mi>).
        (
            "x<y",
            r"<mrow><mi>x</mi><mo>&lt;</mo><mi>y</mi></mrow>",
            r"<apply id='id3'><lt/><ci id='id1'>x</ci><ci id='id2'>y</ci></apply>",
        ),
        // giac y>=x.
        (
            "x<=y",
            r"<mrow><mi>x</mi><mo>&leq;</mo><mi>y</mi></mrow>",
            r"<apply id='id3'><leq/><ci id='id1'>x</ci><ci id='id2'>y</ci></apply>",
        ),
        // giac decides x!=y: true.
        (
            "x!=y",
            r"<mrow><mi>x</mi><mo>&#x2260;</mo><mi>y</mi></mrow>",
            r"<apply id='id3'><neq/><ci id='id1'>x</ci><ci id='id2'>y</ci></apply>",
        ),
        // giac rewrites ceil: -floor(-x).
        (
            "ceil(x)",
            r"<mfenced open='&lceil;' close='&rceil;'><mi>x</mi></mfenced>",
            r"<apply id='id2'><ceiling/><ci id='id1'>x</ci></apply>",
        ),
        // giac rewrites binomial into factorials.
        (
            "binomial(n,k)",
            r"<mfenced><mfrac linethickness='0'><mi>n</mi><mi>k</mi></mfrac></mfenced>",
            r"<apply id='id3'><csymbol cd='combinat1'>binomial</csymbol><ci id='id1'>n</ci><ci id='id2'>k</ci></apply>",
        ),
        // giac (x^(1/3))^2.
        (
            "x^(2/3)",
            r"<msup><mi>x</mi><mfrac><mn>2</mn><mn>3</mn></mfrac></msup>",
            r"<apply id='id3'><power/><ci id='id1'>x</ci><cn id='id2' type='rational'>2<sep/>3</cn></apply>",
        ),
        // Our product order: giac x*sqrt(2).
        (
            "x*sqrt(2)",
            r"<mrow><msqrt><mn>2</mn></msqrt><mo>&it;</mo><mi>x</mi></mrow>",
            r"<apply id='id4'><times/><apply id='id2'><root/><cn id='id1' type='integer'>2</cn></apply><ci id='id3'>x</ci></apply>",
        ),
        // Our closed form, factored: giac (-n-1+(n+1)^2)/2.
        (
            "sum(k,k,1,n)",
            r"<mrow><mfrac><mn>1</mn><mn>2</mn></mfrac><mo>&it;</mo><mi>n</mi><mo>&it;</mo><mfenced><mrow><mn>1</mn><mo>+</mo><mi>n</mi></mrow></mfenced></mrow>",
            r"<apply id='id6'><times/><cn id='id1' type='rational'>1<sep/>2</cn><ci id='id2'>n</ci><apply id='id5'><plus/><cn id='id3' type='integer'>1</cn><ci id='id4'>n</ci></apply></apply>",
        ),
        // Our diff(e, x, 2) differentiates by x, then by 2 (held): giac d^2/dx^2 f(x).
        (
            "diff(f(x),x,2)",
            r"<mrow><mi>derive</mi><mfenced><mrow><mfrac><mi mathvariant='normal'>d</mi><mrow><mi>x</mi></mrow></mfrac><mo>&af;</mo><mrow><mi>f</mi><mo>&af;</mo><mo>&af;</mo><mfenced><mi>x</mi></mfenced></mrow></mrow><mn>2</mn></mfenced></mrow>",
            r"<apply id='id13'><ci id='id6' type='function'>derive</ci><apply id='id11'><diff/><bvar><ci id='id10'>x</ci></bvar><apply id='id9'><ci id='id7'>f</ci><ci id='id8'>x</ci></apply></apply><cn id='id12' type='integer'>2</cn></apply>",
        ),
    ];
    let math = |body: &str| {
        format!("<math mode='display' xmlns='http://www.w3.org/1998/Math/MathML'>{body}</math>")
    };
    let mut wrong = Vec::new();
    for (x, display, content) in cases {
        for (kind, body) in [("display", display), ("content", content)] {
            let input = format!("export_mathml({x},{kind})");
            let found = ours_text(&input);
            if found.as_deref() != Some(math(body).as_str()) {
                wrong.push(format!(
                    "{input}: ours {found:?}, recorded {:?}",
                    math(body)
                ));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Each input `x`: giac's `mathml(x)` is ours, character for character.
fn check_document(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let lines = |s: &str| s.replace('\n', " | ");
    let mut wrong = Vec::new();
    for x in inputs {
        let input = format!("mathml({x})");
        let answer = giac.eval(&input).unwrap_or_default();
        let mine = ours_text(&input);
        let ok = matches!((unquoted(&answer), &mine),
            (Some(theirs), Some(mine)) if lines(&theirs) == lines(mine));
        if !ok {
            wrong.push(format!("{input}: giac {answer}, ours {mine:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The rows of [`answers_as_giac_writes_them`] and [`more_of_giacs_table`]
/// in giac's `mathml` document.
#[test]
fn mathml_documents() {
    check_document(&[&ANSWERS[..], &MORE[..]].concat());
}

/// giac's XHTML around a `mathml` formula (mathml.cc :69, :1980, :70).
fn document(body: &str) -> String {
    format!(
        "<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n\
<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.1 plus MathML 2.0//EN\"\n\
\"http://www.w3.org/TR/MathML2/dtd/xhtml-math11-f.dtd\" [\n\
<!ENTITY mathml \"http://www.w3.org/1998/Math/MathML\">\n\
]>\n<html xmlns=\"http://www.w3.org/1999/xhtml\">\n<body>\n\
\n<math mode=\"display\" xmlns=\"http://www.w3.org/1998/Math/MathML\">\n\n\
{body}\n\n</math><br/>\n\n</body>\n</html>\n"
    )
}

/// The departures of [`departures`] in giac's `mathml` document: ours,
/// the formula of the prototype (`.run/1.6/T5/c3/`) on our shapes (C0068).
#[test]
#[allow(clippy::too_many_lines)]
fn mathml_departures() {
    let cases: &[(&str, &str)] = &[
        // Our sums put constants first: giac sqrt(x+1).
        (
            "sqrt(x+1)",
            r"<msqrt><mrow><mn>1</mn><mo>+</mo><mi>x</mi></mrow></msqrt>",
        ),
        // giac (x^4-1)/(x^2+3).
        (
            "(x^4-1)/(x^2+3)",
            r"<mfrac><mrow><mo>-</mo><mn>1</mn><mo>+</mo><msup><mrow> <mi>x</mi></mrow><mrow><mn>4</mn></mrow></msup></mrow><mrow><mn>3</mn><mo>+</mo><msup><mrow> <mi>x</mi></mrow><mrow><mn>2</mn></mrow></msup></mrow></mfrac>",
        ),
        // Ascending degree: giac 2*x^2-3*x+1.
        (
            "2*x^2-3*x+1",
            r"<mn>1</mn><mo>-</mo><mrow><mn>3</mn><mi>x</mi></mrow><mo>+</mo><mn>2</mn><msup><mrow> <mi>x</mi></mrow><mrow><mn>2</mn></mrow></msup>",
        ),
        // giac 1/(x+1).
        (
            "1/(x+1)",
            r"<mfrac><mrow><mn>1</mn></mrow><mrow><mn>1</mn><mo>+</mo><mi>x</mi></mrow></mfrac>",
        ),
        // giac (x+1)^2.
        (
            "(x+1)^2",
            r#"<msup><mfenced open="(" close=")"><mrow>(<mn>1</mn><mo>+</mo><mi>x</mi>)</mrow></mfenced><mrow><mn>2</mn></mrow></msup>"#,
        ),
        // Our factor order: giac (x+1)*(x-1).
        (
            "(x+1)*(x-1)",
            r" <mo>(</mo><mo>-</mo><mn>1</mn><mo>+</mo><mi>x</mi> <mo>)</mo> <mo>(</mo><mn>1</mn><mo>+</mo><mi>x</mi> <mo>)</mo>",
        ),
        // Our products put cos before sin: giac sin(x)*cos(x).
        (
            "sin(x)*cos(x)",
            r"<mi>cos</mi><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow><mi>sin</mi><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow>",
        ),
        // Ours keeps Times[-1, Plus[1, x]] (C0029): giac -x-1.
        (
            "-(x+1)",
            r"<mo>-</mo><mrow><mo>(</mo><mn>1</mn><mo>+</mo><mi>x</mi><mo>)</mo></mrow>",
        ),
        // giac's evaluator turns x<y into y>x (giac <mi>y</mi><mo>&gt;</mo><mi>x</mi>).
        ("x<y", r"<mi>x</mi><mo>&lt;</mo><mi>y</mi>"),
        // giac y>=x.
        ("x<=y", r"<mi>x</mi><mo>&le;</mo><mi>y</mi>"),
        // giac decides x!=y: true.
        ("x!=y", r"<mi>x</mi><mo>≠</mo><mi>y</mi>"),
        // giac rewrites ceil: -floor(-x).
        (
            "ceil(x)",
            r"<mi>ceil</mi><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow>",
        ),
        // giac rewrites binomial into factorials.
        (
            "binomial(n,k)",
            r"<mi>binomial</mi><mrow><mo>(</mo><mi>n</mi>,<mi>k</mi><mo>)</mo></mrow>",
        ),
        // giac (x^(1/3))^2.
        (
            "x^(2/3)",
            r"<msup><mrow> <mi>x</mi></mrow><mrow><mfrac><mrow><mn>2</mn></mrow><mrow><mn>3</mn></mrow></mfrac></mrow></msup>",
        ),
        // Our product order: giac x*sqrt(2).
        (
            "x*sqrt(2)",
            r"<msqrt><mrow><mn>2</mn></mrow></msqrt><mi>x</mi>",
        ),
        // Our closed form, factored: giac (-n-1+(n+1)^2)/2.
        (
            "sum(k,k,1,n)",
            r"<mfrac><mrow><mi>n</mi> <mo>(</mo><mn>1</mn><mo>+</mo><mi>n</mi> <mo>)</mo></mrow><mrow><mn>2</mn></mrow></mfrac>",
        ),
        // Our diff(e, x, 2) differentiates by x, then by 2 (held): giac d^2/dx^2 f(x).
        (
            "diff(f(x),x,2)",
            r"<mfrac><mrow><mo>&part;</mo><mo>(</mo><mfrac><mrow><mo>&part;</mo><mo>(</mo><mi>f</mi><mo>(</mo><mi>x</mi><mo>)</mo><mo>)</mo></mrow><mrow><mo>&part;</mo><mi>x</mi></mrow></mfrac><mo>)</mo></mrow><mrow><mo>&part;</mo><mn>2</mn></mrow></mfrac>",
        ),
    ];
    let mut wrong = Vec::new();
    for (x, body) in cases {
        let input = format!("mathml({x})");
        let found = ours_text(&input);
        if found.as_deref() != Some(document(body).as_str()) {
            wrong.push(format!("{input}: ours {found:?}, recorded {body:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
