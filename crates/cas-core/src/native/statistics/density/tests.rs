//! Tests of kernel density estimation (CAS-629): `kernel_density`'s three
//! forms against the prototype's direct sums (crates/cas-core/proto/density/proto.py, which
//! reproduces giac's numbers within its FFT noise, C0084), exactness on
//! exact data, each warning, the options not ported and the limits.

use crate::context::Limits;
use crate::trace::Warning;
use crate::{Context, EvalOptions, Expr, Symbol, View};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, text: &str) -> Expr {
    ctx.eval(&read(text))
        .unwrap_or_else(|e| panic!("{text}: {e}"))
}

/// A ball's midpoint in `f64`, or a small integer's value.
fn float(e: &Expr) -> f64 {
    if let View::Ball(b) = e.view() {
        return b.mid().to_f64();
    }
    let i = e
        .as_i64()
        .unwrap_or_else(|| panic!("{e:?}: no small number"));
    i32::try_from(i).map(f64::from).expect("a small integer")
}

/// Whether a ball appears anywhere in `e`.
fn has_ball(e: &Expr) -> bool {
    match e.view() {
        View::Ball(_) | View::ComplexBall(_) => true,
        View::Apply { head, args } => has_ball(head) || args.iter().any(has_ball),
        _ => false,
    }
}

/// `f` at `x = v` in `f64`, through `N`.
fn at(ctx: &Context, f: &Expr, var: &str, v: &str) -> f64 {
    let rule = Expr::call(Symbol::RULE, [Expr::symbol(var), read(v)]);
    let value = Expr::call(
        Symbol::N,
        [Expr::call(Symbol::REPLACE_ALL, [f.clone(), rule])],
    );
    float(&ctx.eval(&value).unwrap())
}

/// Each entry of the list `found` within `1e-12` of the largest expected.
fn close(input: &str, found: &Expr, expected: &[f64]) {
    let scale = expected.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    assert_eq!(
        found.head(),
        Some(&Expr::from(Symbol::LIST)),
        "{input}: {found:?}"
    );
    assert_eq!(found.args().len(), expected.len(), "{input}: {found:?}");
    for (x, v) in found.args().iter().zip(expected) {
        assert!(
            matches!(x.view(), View::Ball(_)) && (float(x) - v).abs() <= 1e-12 * scale,
            "{input}: {found:?}, expected {expected:?}"
        );
    }
}

/// The one warning of `input`, which stays.
fn warning(ctx: &Context, input: &str) -> String {
    let e = read(input);
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    assert_eq!(found.unwrap().head(), e.head(), "{input}");
    match &report.warnings[..] {
        [Warning::Native { message, .. }] => message.to_string(),
        other => panic!("{input}: {other:?}"),
    }
}

const DATA: &str = "List[1, 2, 3, 2]";
/// `[1.0, 2.5, 3.0, 2.0, 4.5]`.
const DECIMALS: &str = "List[0x1p0, 0x1.4p1, 0x1.8p1, 0x1p1, 0x1.2p2]";

/// The exact form: exact on exact data with an exact bandwidth
/// (`(e^(−8(x−1)²) + 2e^(−8(x−2)²) + e^(−8(x−3)²))/√(2π)`), balls for a
/// selected one or decimal data; the variable a bare symbol.
#[test]
fn exact_forms() {
    let ctx = Context::new();
    let f = eval(
        &ctx,
        &format!("KernelDensity[{DATA}, Equation[Bandwidth, 1/4], Exact]"),
    );
    assert!(!has_ball(&f), "{f:?}");
    assert!((at(&ctx, &f, "x", "2") - 0.798_152_221_254_395_1).abs() < 1e-14);
    assert!((at(&ctx, &f, "x", "1/2") - 0.053_990_978_664_953_77).abs() < 1e-14);
    let same = eval(
        &ctx,
        &format!("KernelDensity[{DATA}, Equation[Bandwidth, 1/4], Equation[output, Exact]]"),
    );
    assert_eq!(same, f);
    // kde([1,2,3,2],bandwidth=1/2,exact,t).
    let t = eval(
        &ctx,
        &format!("KernelDensity[{DATA}, Equation[Bandwidth, 1/2], Exact, t]"),
    );
    assert!(!has_ball(&t) && !t.contains(&Expr::symbol("x")), "{t:?}");
    assert!((at(&ctx, &t, "t", "2") - 0.452_933_246_914_620_84).abs() < 1e-14);
    // The plug-in bandwidth, h = 0.8610922659380383.
    let f = eval(&ctx, &format!("KernelDensity[{DATA}, Exact]"));
    assert!(has_ball(&f), "{f:?}");
    for (v, y) in [
        ("1/2", 0.150_372_116_580_086_8),
        ("2", 0.349_673_926_693_619_3),
        ("13/4", 0.195_625_487_494_870_2),
    ] {
        assert!((at(&ctx, &f, "x", v) - y).abs() < 1e-12, "{v}: {f:?}");
    }
    // The rule of thumb, h = 0.599499836007844.
    let f = eval(
        &ctx,
        &format!("KernelDensity[{DATA}, Equation[Bandwidth, GaussSquares], Exact]"),
    );
    assert!((at(&ctx, &f, "x", "2") - 0.415_503_909_590_284_9).abs() < 1e-12);
    // Decimal data, h = 1.2432990560829538.
    let f = eval(&ctx, &format!("KernelDensity[{DECIMALS}, Exact]"));
    assert!((at(&ctx, &f, "x", "1") - 0.160_427_744_221_985_9).abs() < 1e-12);
    assert!((at(&ctx, &f, "x", "11/4") - 0.226_943_481_290_805).abs() < 1e-12);
}

/// The list form against the prototype's direct sums.
#[test]
fn lists() {
    let ctx = Context::new();
    let rows: [(&str, &[f64]); 8] = [
        (
            "Equation[Bandwidth, 1/4], Equation[Bins, 4]",
            &[
                7.446_045_870_629_99e-6,
                0.398_964_618_539_044_6,
                1.196_834_287_250_168_7,
                2.233_813_761_188_997e-5,
            ],
        ),
        (
            "Equation[Bins, 7]",
            &[
                0.241_654_462_671_644_12,
                0.296_920_048_459_782_6,
                0.335_684_871_540_775_64,
                0.349_673_926_693_619_3,
                0.335_684_871_540_775_64,
                0.296_920_048_459_782_54,
                0.241_654_462_671_644_1,
            ],
        ),
        (
            "Equation[Bandwidth, 1/4], Equation[range, List[0, 4]], Equation[Bins, 4]",
            &[
                2.656_301_908_521_912_5e-7,
                0.398_943_077_292_005_23,
                1.196_827_106_834_489_1,
                7.968_905_725_565_738e-7,
            ],
        ),
        (
            "Equation[Bandwidth, 1/4], Range[0, 4], Equation[Bins, 4]",
            &[
                2.656_301_908_521_912_5e-7,
                0.398_943_077_292_005_23,
                1.196_827_106_834_489_1,
                7.968_905_725_565_738e-7,
            ],
        ),
        (
            "Equation[Bandwidth, 1/4], Equation[Bins, 3], Equation[output, list]",
            &[
                0.398_942_280_419_702_14,
                0.797_884_560_821_134_8,
                0.398_942_280_419_702_14,
            ],
        ),
        (
            "Equation[Bandwidth, 1/4], Equation[Bins, 2]",
            &[0.398_942_280_401_432_7, 1.196_826_841_204_298_2],
        ),
        (
            "Equation[Bandwidth, GaussSquares], Equation[Bins, 5]",
            &[
                0.010_923_345_122_723_206,
                0.188_202_068_015_440_7,
                0.354_563_615_826_848_9,
                0.188_202_068_015_440_7,
                0.010_923_345_122_723_206,
            ],
        ),
        (
            "Exact, Equation[Bandwidth, 1/4], Equation[Bins, 4], Equation[output, list]",
            &[
                7.446_045_870_629_99e-6,
                0.398_964_618_539_044_6,
                1.196_834_287_250_168_7,
                2.233_813_761_188_997e-5,
            ],
        ),
    ];
    for (options, expected) in rows {
        let input = format!("KernelDensity[{DATA}, {options}]");
        close(&input, &eval(&ctx, &input), expected);
    }
}

/// The list form on decimal data, and past 1000 data.
#[test]
fn more_lists() {
    let ctx = Context::new();
    // Decimal data, the rule of thumb (`normal`).
    let input =
        format!("KernelDensity[{DECIMALS}, Equation[Bandwidth, Together], Equation[Bins, 5]]");
    close(
        &input,
        &eval(&ctx, &input),
        &[
            0.001_428_804_968_962_359_4,
            0.143_452_354_729_941_75,
            0.420_355_457_240_067_4,
            0.143_452_354_729_941_75,
            0.001_428_804_968_962_359_4,
        ],
    );
    // Past 1000 data, the plug-in rule on the bins (h = 1.7598446249772945).
    let squares: Vec<String> = (1..=1001u64).map(|k| (k * k % 37).to_string()).collect();
    let input = format!(
        "KernelDensity[List[{}], Equation[Bins, 12]]",
        squares.join(", ")
    );
    close(
        &input,
        &eval(&ctx, &input),
        &[
            0.022_962_072_566_460_295,
            0.030_185_773_598_022_856,
            0.023_148_859_706_239_085,
            0.041_063_577_879_167_985,
            0.020_932_590_760_565_157,
            0.016_605_104_458_020_728,
            0.014_447_421_243_884_971,
            0.008_715_586_618_556_098,
            0.041_063_352_332_718_33,
            0.035_319_400_003_778_666,
            0.031_004_027_181_920_538,
            0.016_592_982_492_018_862,
        ],
    );
}

/// The piecewise form: giac's breakpoints, exact on an exact grid, and its
/// linear pieces; 0 outside the grid.
#[test]
fn piecewise_forms() {
    let ctx = Context::new();
    let input = format!(
        "KernelDensity[{DATA}, Equation[Bandwidth, 1/4], Equation[range, List[0, 4]], Equation[Bins, 4], Piecewise]"
    );
    let f = eval(&ctx, &input);
    assert_eq!(f.head(), Some(&Expr::from(Symbol::PIECEWISE)), "{f:?}");
    let args = f.args();
    assert_eq!(args.len(), 9, "{f:?}");
    assert_eq!(args[0], eval(&ctx, "Less[x, 0]"));
    assert_eq!(args[1], Expr::from(0));
    assert_eq!(args[2], eval(&ctx, "Less[x, 4/3]"));
    assert_eq!(args[4], eval(&ctx, "Less[x, 8/3]"));
    assert_eq!(args[6], eval(&ctx, "LessEqual[x, 4]"));
    assert_eq!(args[8], Expr::from(0));
    for (v, y) in [
        ("1/2", 0.149_603_820_003_371_25),
        ("2", 0.797_885_092_063_247_2),
        ("10/3", 0.598_413_951_862_530_5),
        ("5", 0.0),
        ("-1", 0.0),
    ] {
        assert!((at(&ctx, &f, "x", v) - y).abs() < 1e-12, "{v}: {f:?}");
    }
    // t=0..4 names the variable and sets the grid.
    let input = format!(
        "KernelDensity[{DATA}, Equation[t, Range[0, 4]], Equation[Bandwidth, 1/4], Equation[Bins, 4], Equation[output, Piecewise]]"
    );
    let g = eval(&ctx, &input);
    assert_eq!(g.args()[2], eval(&ctx, "Less[t, 4/3]"));
    assert!((at(&ctx, &g, "t", "2") - 0.797_885_092_063_247_2).abs() < 1e-12);
    // interp=1 is the default.
    let input = format!(
        "KernelDensity[{DATA}, Equation[Bandwidth, 1/4], Equation[range, List[0, 4]], Equation[Bins, 4], Equation[InterpolatingPolynomial, 1], Piecewise]"
    );
    assert_eq!(eval(&ctx, &input), f);
}

/// giac's errors stay with its message; options not ported stay without
/// one; the limits.
#[test]
fn stays() {
    let ctx = Context::new();
    for (input, message) in [
        (
            "KernelDensity[5]",
            "Expected a list or a sequence of arguments",
        ),
        ("KernelDensity[List[1], Exact]", "Data too small"),
        (
            "KernelDensity[List[1, 2, a, 2]]",
            "Expected a real constant",
        ),
        (
            "KernelDensity[List[1, 2, 3, 2], Equation[Bandwidth, -1]]",
            "Invalid bandwidth specification",
        ),
        (
            "KernelDensity[List[1, 2, 3, 2], Equation[Bandwidth, foo]]",
            "Invalid bandwidth specification",
        ),
        (
            "KernelDensity[List[1, 2, 3, 2], Equation[Bins, 0]]",
            "Number of bins must be positive",
        ),
        (
            "KernelDensity[List[1, 2, 3, 2], Equation[Bins, 5/2]]",
            "Number of bins must be positive",
        ),
        (
            "KernelDensity[List[1, 2, 3, 2], Equation[range, List[4, 0]]]",
            "Invalid range specification",
        ),
        (
            "KernelDensity[List[1, 2, 3, 2], Equation[range, 3]]",
            "Invalid range specification",
        ),
        (
            "KernelDensity[List[1, 2, 3, 2], Range[4, 0]]",
            "Invalid range specification",
        ),
        (
            "KernelDensity[List[1, 2, 3, 2], Equation[t, 3]]",
            "Invalid range specification",
        ),
        (
            "KernelDensity[List[1, 2, 3, 2], Equation[foo, List[1, 2]]]",
            "Invalid range specification",
        ),
        (
            "KernelDensity[List[1, 2, 3, 2], Equation[output, foo]]",
            "Invalid method specification",
        ),
        (
            "KernelDensity[List[1, 2, 3, 2], Equation[InterpolatingPolynomial, 0]]",
            "Invalid interpolation method specification",
        ),
        (
            "KernelDensity[List[1, 2, 3, 2], Equation[Spline, 0]]",
            "Invalid interpolation method specification",
        ),
        ("KernelDensity[List[1, 2, 3, 2], 3]", "Unrecognized option"),
        (
            "KernelDensity[List[1, 2, 3, 2], Equation[Bandwidth, 1/4], Equation[Bins, 1]]",
            "Invalid dimension",
        ),
    ] {
        assert_eq!(warning(&ctx, input), message, "{input}");
    }
    for input in [
        "KernelDensity[]",
        "KernelDensity[List[1, 2, 3, 2], Equation[InterpolatingPolynomial, 2]]",
        "KernelDensity[List[1, 2, 3, 2], Equation[Spline, 1]]",
        "KernelDensity[List[1, 2, 3, 2], Equation[Evaluate, 2]]",
        // All data equal: giac's plug-in rule divides by 0.
        "KernelDensity[List[2, 2, 2]]",
    ] {
        let e = read(input);
        let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
        assert_eq!(found.unwrap(), e, "{input}");
        assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
    }
    let mut small = Context::new();
    let limits = small.limits();
    small.set_limits(Limits {
        terms: 20,
        ..limits
    });
    let many: Vec<String> = (1..=25).map(|k: u32| k.to_string()).collect();
    for input in [
        // 4 bins times 4 data within, 10 bins not.
        "KernelDensity[List[1, 2, 3, 2], Equation[Bandwidth, 1/4], Equation[Bins, 10]]".to_owned(),
        // The kernel's weights: 4 bins of 2L + 1 = 7 weights (L = ⌊1 + 3⌋).
        "KernelDensity[List[1, 2, 3, 2], Equation[Bandwidth, 1], Equation[range, List[0, 4]], Equation[Bins, 4]]".to_owned(),
        // 25 data in the exact form.
        format!("KernelDensity[List[{}], Exact]", many.join(", ")),
    ] {
        assert_eq!(warning(&small, &input), "too many points", "{input}");
    }
    small.set_limits(Limits { terms: 5, ..limits });
    // 4 distinct values, 6 pairs, for the plug-in rule.
    assert_eq!(
        warning(&small, "KernelDensity[List[1, 2, 3, 4], Exact]"),
        "too many points"
    );
}
