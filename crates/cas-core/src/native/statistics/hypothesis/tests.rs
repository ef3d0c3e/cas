//! Tests of the hypothesis tests (CAS-DOC-52): giac's answers where its
//! conventions are this CAS's, and the calls that stay.

use crate::{Context, Expr};

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

/// Each input stays as it is.
fn stays(inputs: &[&str]) {
    let ctx = Context::new();
    for input in inputs {
        let e = ctx.eval(&read(input)).unwrap();
        assert_eq!(e, ctx.eval(&read(input)).unwrap());
        assert_eq!(e, read(input), "{input}: {e:?}");
    }
}

/// `ZTest` and `TTest` (CAS-503).
#[test]
fn location_tests() {
    check(&[
        // giac's normalt([10,30],.5,.02,'!=',0.1): 0.
        ("ZTest[List[10, 30], 1/2, 1/50, \"!=\", 1/10]", "0"),
        // giac's normalt([0.48,50],0.5,0.1,'<'): z = −√2, p = 0.0786.
        ("ZTest[List[12/25, 50], 1/2, 1/10, \"<\"]", "1"),
        ("ZTest[List[12/25, 50], 1/2, 1/10, \"<\", 1/10]", "0"),
        // z = 3 √3 / 2 ≈ 2.6: rejected both ways at 5 %, but not "<".
        ("ZTest[List[3, 2, 5], 1, 1, \">\"]", "0"),
        ("ZTest[List[3, 2, 5], 1, 1]", "0"),
        ("ZTest[List[3, 2, 5], 1, 1, \"<\"]", "1"),
        // giac's studentt([1,2,3,4,5],2,'>'): t = √2, 4 df, p = 0.115.
        ("TTest[List[1, 2, 3, 4, 5], 2, \">\"]", "1"),
        ("TTest[List[1, 2, 3, 4, 5], 0, \">\"]", "0"),
        // giac's studentt([0.48,20],0.5,0.1,'<').
        ("TTest[List[12/25, 20], 1/2, 1/10, \"<\"]", "1"),
        // Decimal data: balls.
        ("ZTest[List[0.48, 50], 0.5, 0.1, \"<\"]", "1"),
    ]);
    stays(&[
        // No σ for a Z test.
        "ZTest[List[1, 2, 3], 1]",
        // An alternative that is not one.
        "ZTest[List[1, 2, 3], 1, 1, \"=\"]",
        // A level outside (0, 1).
        "ZTest[List[1, 2, 3], 1, 1, \"<\", 2]",
        // A size that is not a positive integer; a negative σ.
        "ZTest[List[1, 1/2], 1, 1]",
        "ZTest[List[1, 2, 3], 1, -1]",
        // Symbols.
        "TTest[List[a, b, c], 1]",
    ]);
}

/// `ChiSquareTest` (CAS-504): giac's examples, exact.
#[test]
fn chi_square() {
    check(&[
        // giac's 1.10344827586.
        ("ChiSquareTest[List[57, 30], List[3/5, 2/5]]", "32/29"),
        (
            "ChiSquareTest[List[57, 30], List[3/5, 2/5], 1/100]",
            "32/29",
        ),
        // Uniform: giac's 6.57142857143 and 0.0810810810811.
        ("ChiSquareTest[List[1, 1, 1, 5, 2, 4]]", "46/7"),
        ("ChiSquareTest[List[57, 54]]", "3/37"),
        (
            "ChiSquareTest[List[11, 16, 17, 22, 14, 10], List[1/6, 1/6, 1/6, 1/6, 1/6, 1/6]]",
            "32/5",
        ),
        // Two samples: giac's 0.202238009215 and 0.0634482758621.
        (
            "ChiSquareTest[List[17, 15, 12, 15], List[15, 13, 13, 14]]",
            "26644251/131747000",
        ),
        ("ChiSquareTest[List[57, 30], List[3, 2]]", "46/725"),
        // Three: giac's 2.40477818022.
        (
            "ChiSquareTest[List[List[17, 15, 12, 15], List[15, 13, 13, 14], List[1, 2, 3, 4]]]",
            "6008699/2498650",
        ),
        // Raw data, counted: 10 zeros and 13 ones; giac's 0.115942028986.
        (
            "ChiSquareTest[List[1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0], List[2/5, 3/5]]",
            "8/69",
        ),
    ]);
    stays(&[
        "ChiSquareTest[List[-1, 3]]",
        "ChiSquareTest[List[a, 3]]",
        // Three distinct values for two probabilities.
        "ChiSquareTest[List[1, 2, 3], List[1/2, 1/2]]",
        // No expected count above 0.
        "ChiSquareTest[List[0, 0], List[1, 1]]",
    ]);
}

/// The ball of `input`'s value, as a double, and its radius.
fn ball(input: &str) -> (f64, f64) {
    let ctx = Context::new();
    let value = ctx.eval(&read(input)).unwrap();
    match value.view() {
        crate::View::Ball(b) => (
            b.mid().to_f64(),
            b.rad().to_float().map_or(f64::INFINITY, |r| r.to_f64()),
        ),
        _ => panic!("{input}: {value:?}"),
    }
}

/// `KolmogorovCDF` (CAS-505), against both series in doubles: giac's
/// `kolmogorovd(1.0)` is 0.730000328323.
#[test]
fn kolmogorov_distribution() {
    for (x, expected) in [
        ("1", 0.730_000_328_322_645_5),
        ("9999/10000", 0.729_893_118_031_144_5),
        ("10001/10000", 0.730_107_507_743_515_6),
        ("3", 0.999_999_969_540_040_6),
        ("1/5", 5.050_407_338_670_114e-13),
        // giac's 2.2e-16 is its rounding.
        ("1/20", 2.423_167_479_157_9e-213),
    ] {
        let (mid, rad) = ball(&format!("N[KolmogorovCDF[{x}], 30]"));
        assert!(
            ((mid - expected) / expected).abs() < 1e-13,
            "K({x}) = {mid}, expected {expected}"
        );
        assert!(rad <= mid.abs() * 1e-25, "K({x}): radius {rad}");
    }
    check(&[
        ("KolmogorovCDF[0]", "0"),
        ("KolmogorovCDF[-2]", "0"),
        ("KolmogorovCDF[Ball[-1/2, 0]]", "0"),
        ("KolmogorovCDF[1]", "KolmogorovCDF[1]"),
        ("KolmogorovCDF[x]", "KolmogorovCDF[x]"),
    ]);
}

/// `KolmogorovSmirnovTest` (CAS-506): giac's
/// `kolmogorovt([1,2,3],[1.5,2.5,8])` is `["D=",0.333,"K=",0.408,…]`, and
/// against `uniformd(0,2)` `D` is 0.55.
#[test]
fn kolmogorov_smirnov() {
    check(&[
        (
            "KolmogorovSmirnovTest[List[1, 2, 3], List[3/2, 5/2, 8]]",
            "List[\"D=\", 1/3, \"K=\", Times[1/6, Power[6, 1/2]], \"1-kolmogorovd(K)=\", Plus[1, Times[-1, KolmogorovCDF[Times[1/6, Power[6, 1/2]]]]]]",
        ),
        // Ties across the samples step together.
        (
            "Part[KolmogorovSmirnovTest[List[1, 1, 2], List[1, 2, 2]], 2]",
            "1/3",
        ),
        (
            "Part[KolmogorovSmirnovTest[List[1/10, 1/2, 9/10], UniformPDF[0, 2]], 2]",
            "11/20",
        ),
        (
            "KolmogorovSmirnovTest[List[a, 1], List[1, 2]]",
            "KolmogorovSmirnovTest[List[a, 1], List[1, 2]]",
        ),
        (
            "KolmogorovSmirnovTest[List[1, 2], f[5]]",
            "KolmogorovSmirnovTest[List[1, 2], f[5]]",
        ),
    ]);
    // giac's D = 0.539827837277, K = 0.935009241504, 1 − K(K) =
    // 0.346240413235 against normald(0,1).
    let ctx = Context::new();
    let found = ctx
        .eval(&read(
            "KolmogorovSmirnovTest[List[N[1/10], 1/2, 9/10], NormalPDF[0, 1]]",
        ))
        .unwrap();
    let expected = [0.539_827_837_277, 0.935_009_241_504, 0.346_240_413_235];
    for (k, want) in [1, 3, 5].into_iter().zip(expected) {
        let crate::View::Ball(b) = found.args()[k].view() else {
            panic!("{found:?}");
        };
        assert!((b.mid().to_f64() - want).abs() < 1e-11, "{found:?}");
    }
}
