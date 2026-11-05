//! Tests of the signal processing builtins (CAS-521 …): exact windows,
//! giac's values, parameters, slices, and what stays.

use crate::{Context, Expr, Symbol, View};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, expected)` in full form, both evaluated: the same expression.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let expected = ctx.eval(&read(expected)).unwrap();
        if found != expected {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The doubles of a list's entries, by `N`.
fn doubles(ctx: &Context, e: &Expr) -> Vec<f64> {
    e.args()
        .iter()
        .map(|x| {
            let v = ctx.eval(&Expr::call(Symbol::N, [x.clone()])).unwrap();
            if let View::Ball(b) = v.view() {
                return b.mid().to_f64();
            }
            let Some(q) = x.to_rational() else {
                panic!("{x:?}: {v:?}");
            };
            let rounded = cas_num::ball::Round::Nearest;
            cas_num::ball::Float::from_rational(&q, 64, rounded)
                .0
                .to_f64()
        })
        .collect()
}

const FIVE: &str = "List[1, 1, 1, 1, 1]";
const SIX: &str = "List[1, 1, 1, 1, 1, 1]";

/// Exact windows of five and six ones.
#[test]
fn exact_windows() {
    let cases = [
        ("HannWindow", FIVE, "List[0, 1/2, 1, 1/2, 0]"),
        ("HammingWindow", FIVE, "List[2/25, 27/50, 1, 27/50, 2/25]"),
        ("BlackmanWindow", FIVE, "List[0, 17/50, 1, 17/50, 0]"),
        (
            "BlackmanHarrisWindow",
            FIVE,
            "List[3/50000, 21747/100000, 1, 21747/100000, 3/50000]",
        ),
        ("BartlettHannWindow", FIVE, "List[0, 1/2, 1, 1/2, 0]"),
        (
            "BohmanWindow",
            FIVE,
            "List[0, Power[Pi, -1], 1, Power[Pi, -1], 0]",
        ),
        (
            "RiemannWindow",
            FIVE,
            "List[0, Times[2, Power[Pi, -1]], 1, Times[2, Power[Pi, -1]], 0]",
        ),
        (
            "CosineWindow",
            FIVE,
            "List[0, Power[2, -1/2], 1, Power[2, -1/2], 0]",
        ),
        ("TriangleWindow", FIVE, "List[1/5, 3/5, 1, 3/5, 1/5]"),
        ("WelchWindow", SIX, "List[0, 16/25, 24/25, 24/25, 16/25, 0]"),
        (
            "ParzenWindow",
            SIX,
            "List[0, 16/125, 101/125, 101/125, 16/125, 0]",
        ),
        ("ParzenWindow", FIVE, "List[0, 1/4, 1, 1/4, 0]"),
        ("TukeyWindow", FIVE, "List[0, 1, 1, 1, 0]"),
    ];
    let cases: Vec<(String, String)> = cases
        .iter()
        .map(|(w, u, expected)| (format!("{w}[{u}]"), (*expected).to_owned()))
        .collect();
    let cases: Vec<(&str, &str)> = cases
        .iter()
        .map(|(a, b)| (a.as_str(), b.as_str()))
        .collect();
    check(&cases);
    check(&[
        // Symbolic samples stay symbolic: giac's [0.0,0.75*b,0.75*c,…].
        (
            "HannWindow[List[a, b, c, d]]",
            "List[0, Times[3/4, b], Times[3/4, c], 0]",
        ),
        // Parameters: the triangle's d, Blackman's a, Tukey's ends, the
        // cosine's power.
        (
            "TriangleWindow[List[1, 1, 1, 1, 1], 1]",
            "List[1/3, 2/3, 1, 2/3, 1/3]",
        ),
        (
            "TriangleWindow[List[1, 1, 1, 1, 1], -1]",
            "List[0, 1/2, 1, 1/2, 0]",
        ),
        (
            "BlackmanWindow[List[1, 1, 1, 1, 1], 1/5]",
            "List[0, 3/10, 1, 3/10, 0]",
        ),
        ("TukeyWindow[List[1, 1, 1, 1, 1], 0]", "List[1, 1, 1, 1, 1]"),
        (
            "TukeyWindow[List[1, 1, 1, 1, 1], 1]",
            "List[0, 1/2, 1, 1/2, 0]",
        ),
        (
            "CosineWindow[List[1, 1, 1, 1, 1], 2]",
            "List[0, 1/2, 1, 1/2, 0]",
        ),
        // A slice, from index 1 to 5: its own window.
        (
            "HannWindow[List[9, 1, 1, 1, 1, 1, 9], Range[1, 5]]",
            "List[0, 1/2, 1, 1/2, 0]",
        ),
        (
            "TriangleWindow[List[9, 1, 1, 1, 1, 1, 9], 1, Range[1, 5]]",
            "List[1/3, 2/3, 1, 2/3, 1/3]",
        ),
    ]);
}

/// The windows with transcendental values, against giac's doubles.
#[test]
fn giac_values() {
    let ctx = Context::new();
    let cases: [(&str, &[f64]); 8] = [
        (
            "GaussianWindow[List[1, 1, 1, 1, 1]]",
            &[
                1.928_749_847_963_918e-22,
                3.726_653_172_078_671e-6,
                1.0,
                3.726_653_172_078_671e-6,
                1.928_749_847_963_918e-22,
            ],
        ),
        (
            "GaussianWindow[List[1, 1, 1, 1, 1], 3/10]",
            &[
                0.003_865_920_139_472_5,
                0.249_352_208_777_296,
                1.0,
                0.249_352_208_777_296,
                0.003_865_920_139_472_5,
            ],
        ),
        (
            "PoissonWindow[List[1, 1, 1, 1, 1, 1]]",
            &[
                0.367_879_441_171,
                0.548_811_636_094,
                0.818_730_753_078,
                0.818_730_753_078,
                0.548_811_636_094,
                0.367_879_441_171,
            ],
        ),
        (
            "HannPoissonWindow[List[1, 1, 1, 1, 1, 1]]",
            &[
                0.0,
                0.189_609_756_915,
                0.740_548_923_068,
                0.740_548_923_068,
                0.189_609_756_915,
                0.0,
            ],
        ),
        (
            "HannPoissonWindow[List[1, 1, 1, 1, 1], 2]",
            &[0.0, 0.183_939_720_586, 1.0, 0.183_939_720_586, 0.0],
        ),
        (
            "TukeyWindow[List[1, 1, 1, 1, 1, 1]]",
            &[0.0, 0.904_508_497_187, 1.0, 1.0, 0.904_508_497_187, 0.0],
        ),
        (
            "BartlettHannWindow[List[1, 1, 1, 1, 1, 1]]",
            &[
                0.0,
                0.358_573_542_138,
                0.879_426_457_862,
                0.879_426_457_862,
                0.358_573_542_138,
                0.0,
            ],
        ),
        (
            "BohmanWindow[List[1, 1, 1, 1, 1, 1]]",
            &[
                0.0,
                0.179_123_893_706,
                0.834_311_452_258,
                0.834_311_452_258,
                0.179_123_893_706,
                0.0,
            ],
        ),
    ];
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let values = doubles(&ctx, &found);
        assert_eq!(values.len(), expected.len(), "{input}: {found:?}");
        for (a, b) in values.iter().zip(expected) {
            assert!(
                (a - b).abs() <= 1e-11 * (1.0 + b.abs()),
                "{input}: {values:?}, expected {expected:?}"
            );
        }
    }
}

/// A decimal parameter makes the weights numbers, as giac's.
#[test]
fn inexact_parameters() {
    let ctx = Context::new();
    let found = ctx
        .eval(&read("TukeyWindow[List[1, 1, 1, 1, 1, 1], Decimal[4, -1]]"))
        .unwrap();
    for x in found.args() {
        assert!(
            matches!(x.view(), View::Ball(_)) || x.to_rational().is_some(),
            "{found:?}"
        );
    }
    let values = doubles(&ctx, &found);
    let expected = [0.0, 1.0, 1.0, 1.0, 1.0, 0.0];
    for (a, b) in values.iter().zip(&expected) {
        assert!((a - b).abs() <= 1e-12, "{values:?}");
    }
}

/// Out of range parameters, too few samples, bad slices: the call stays.
#[test]
fn what_stays() {
    for input in [
        "HannWindow[List[1]]",
        "HannWindow[List[]]",
        "GaussianWindow[List[1, 1, 1], 3/5]",
        "TriangleWindow[List[1, 1, 1], 2]",
        "TukeyWindow[List[1, 1, 1], 2]",
        "CosineWindow[List[1, 1, 1], 0]",
        "PoissonWindow[List[1, 1, 1], -1]",
        "HannWindow[List[1, 1, 1], Range[1, 3]]",
        "HannWindow[List[1, 1, 1], Range[2, 2]]",
        "HammingWindow[List[1, 1, 1], 2]",
        "HannWindow[x]",
    ] {
        check(&[(input, input)]);
    }
}

/// RC filters (CAS-522): giac's values, an impulse decaying by `1 − α` and
/// the high pass of a constant by `1 − α` per step; what stays.
#[test]
fn rc_filters() {
    let ctx = Context::new();
    let cases: [(&str, &[f64]); 5] = [
        (
            "LowPass[List[1, 0, 0, 0, 0], 1000]",
            &[
                1.0,
                0.875_292_019_174,
                0.766_136_118_83,
                0.670_592_830_412,
                0.586_964_552_575,
            ],
        ),
        (
            "LowPass[List[1, 1, 1, 1], 1000, 10000]",
            &[1.0, 1.0, 1.0, 1.0],
        ),
        (
            "HighPass[List[1, 1, 1, 1], 1000, 10000]",
            &[1.0, 0.614_130_454_902, 0.377_156_215_638, 0.231_623_118_279],
        ),
        (
            "HighPass[List[1, 0, 0, 0], 1000, 10000]",
            &[1.0, 0.0, 0.0, 0.0],
        ),
        (
            "LowPass[List[1, 2, 3], 1000]",
            &[1.0, 1.124_707_980_83, 1.358_571_862],
        ),
    ];
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let values = doubles(&ctx, &found);
        assert_eq!(values.len(), expected.len(), "{input}: {found:?}");
        for (a, b) in values.iter().zip(expected) {
            assert!((a - b).abs() <= 1e-11, "{input}: {values:?}");
        }
    }
    check(&[("LowPass[List[], 1000]", "List[]")]);
    for input in [
        "LowPass[List[a, b], 1000]",
        "LowPass[List[1, 2], -1]",
        "HighPass[List[1, 2], 1000, 0]",
        "LowPass[List[1, 2], c]",
        "HighPass[x, 1000]",
    ] {
        check(&[(input, input)]);
    }
}

/// Correlations (CAS-523): giac's examples exactly, symbolic and complex
/// samples, empty lists.
#[test]
fn correlations() {
    check(&[
        (
            "CrossCorrelation[List[1, 2], List[3, 4, 5]]",
            "List[6, 11, 14, 5]",
        ),
        (
            "AutoCorrelation[List[1, -1, 0, 2, 1]]",
            "List[1, 1, -2, 1, 7, 1, -2, 1, 1]",
        ),
        ("CrossCorrelation[List[1, 2, 3], List[1]]", "List[3, 2, 1]"),
        ("CrossCorrelation[List[1], List[1, 2, 3]]", "List[1, 2, 3]"),
        (
            "CrossCorrelation[List[a, b], List[c, d]]",
            "List[Times[b, c], Plus[Times[a, c], Times[b, d]], Times[a, d]]",
        ),
        // The conjugate of the first: r_0 of [1, i] with itself is 2.
        (
            "AutoCorrelation[List[1, Complex[0, 1]]]",
            "List[Complex[0, -1], 2, Complex[0, 1]]",
        ),
        ("CrossCorrelation[List[], List[1]]", "List[]"),
        ("AutoCorrelation[List[2]]", "List[4]"),
        ("AutoCorrelation[x]", "AutoCorrelation[x]"),
    ]);
}
