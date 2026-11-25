//! Tests of the wavelet transforms (CAS-634): giac's rows, which the
//! prototype (`.run/1.7/T9/proto_dwt.py`, `logs/T9.c1.md`) reproduces from GSL
//! to giac's printed digits, Haar exactly, round trips, padding, giac's
//! errors with their warnings, and the limits.

#![allow(
    clippy::approx_constant,
    reason = "giac's printed digits of √2 and 1/√2 are expected values"
)]

use crate::context::Limits;
use crate::expr::View;
use crate::trace::Warning;
use crate::{Context, EvalOptions, Expr};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, text: &str) -> Expr {
    ctx.eval(&read(text))
        .unwrap_or_else(|e| panic!("{text}: {e}"))
}

/// The entries of nested lists, in order.
fn leaves(e: &Expr) -> Vec<Expr> {
    if e.head() == Some(&Expr::from(crate::Symbol::LIST)) {
        e.args().iter().flat_map(leaves).collect()
    } else {
        vec![e.clone()]
    }
}

/// The value of `text` as doubles, entry by entry, through `N`.
fn doubles(ctx: &Context, text: &str) -> Vec<f64> {
    leaves(&eval(ctx, &format!("N[{text}]")))
        .iter()
        .map(|x| match x.view() {
            View::Ball(b) => b.mid().to_f64(),
            View::Int(n) => f64::from(i32::try_from(n).unwrap()),
            _ => panic!("{text}: {x:?}"),
        })
        .collect()
}

/// `(input, expected)`: each entry within `10⁻⁹` of giac's, relative to
/// `max(1, |giac's|)`, its 12 printed digits.
fn near(ctx: &Context, cases: &[(String, &[f64])]) {
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = doubles(ctx, input);
        let close = found.len() == expected.len()
            && found
                .iter()
                .zip(*expected)
                .all(|(f, g)| (f - g).abs() <= 1e-9 * g.abs().max(1.0));
        if !close {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `(input, expected)` in full form, both evaluated: the same expression.
fn check(ctx: &Context, cases: &[(String, &str)]) {
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let found = eval(ctx, input);
        let expected = eval(ctx, expected);
        if found != expected {
            wrong.push(format!("{input}: {found:?}, expected {expected:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `input` with its arguments evaluated: the call as it stays.
fn stayed(ctx: &Context, input: &str) -> Expr {
    read(input).map_args(|a| ctx.eval(a).unwrap_or_else(|e| panic!("{input}: {e}")))
}

/// The value of `input` and the messages of its warnings.
fn with_warnings(ctx: &Context, input: &str) -> (Expr, Vec<String>) {
    let (found, report) = ctx.eval_with(&read(input), &EvalOptions::default());
    let messages = report
        .warnings
        .iter()
        .map(|w| match w {
            Warning::Native { message, .. } => message.to_string(),
            other => panic!("{input}: {other:?}"),
        })
        .collect();
    (found.unwrap(), messages)
}

fn dwt(args: &str) -> String {
    format!("DiscreteWaveletTransform[{args}]")
}

fn idwt(args: &str) -> String {
    format!("InverseDiscreteWaveletTransform[{args}]")
}

const V4: &str = "List[1, 2, 3, 4]";
const V8: &str = "List[1, 2, 3, 4, 5, 6, 7, 8]";
const M44: &str =
    "List[List[15, 25, 5, 70], List[21, 33, 12, 7], List[64, 18, 50, 10], List[5, 17, 35, 28]]";
const M4B: &str = "List[List[1, 2, 3, 4], List[4, 5, 6, 7], List[1, 0, 0, 1], List[2, 2, 1, 1]]";
const R2: &str = "Power[2, 1/2]";

/// Haar on exact data: `a + b√2` exactly (the prototype's `Q2` rows).
#[test]
fn exact_haar() {
    let ctx = Context::new();
    let s = |b: &str| format!("Times[{b}, {R2}]");
    check(
        &ctx,
        &[
            (
                dwt(&format!("{V4}, \"haar\"")),
                &format!("List[5, -2, {}, {}]", s("-1/2"), s("-1/2")),
            ),
            (
                dwt(&format!("{V4}, \"haar\", center")),
                &format!("List[5, 0, {}, {}]", s("3/2"), s("-1/2")),
            ),
            (
                dwt(&format!("{V4}, \"haar\", centre")),
                &format!("List[5, 0, {}, {}]", s("3/2"), s("-1/2")),
            ),
            (
                dwt("List[1, 2, 3, 4, 5], \"haar\""),
                &format!(
                    "List[{}, {}, -2, 5/2, {}, {}, {}, 0]",
                    s("15/4"),
                    s("5/4"),
                    s("-1/2"),
                    s("-1/2"),
                    s("5/2")
                ),
            ),
            (
                dwt("List[List[1, 2], List[3, 4], List[5, 6]], \"haar\""),
                &format!(
                    "List[List[21/4, 21/4, {}, 0], List[-1/4, -1/4, {}, 0], List[{}, {}, 0, 0], List[{}, {}, -1/2, 0]]",
                    s("-3/4"),
                    s("-1/4"),
                    s("-1"),
                    s("-1"),
                    s("11/4"),
                    s("11/4")
                ),
            ),
            (
                idwt(&format!("{}, \"haar\"", dwt(&format!("{V4}, \"haar\"")))),
                V4,
            ),
            (
                idwt(&format!(
                    "{}, \"haar\", ColumnSpace",
                    dwt(&format!("{M4B}, \"haar\", ColumnSpace"))
                )),
                M4B,
            ),
            (
                idwt(&format!(
                    "{}, \"haar\", 2",
                    dwt(&format!("{M44}, \"haar\""))
                )),
                M44,
            ),
            (dwt("List[5], \"haar\""), "List[5]"),
            (dwt("List[]"), "List[]"),
            (idwt("List[], \"bspline\""), "List[]"),
        ],
    );
    near(
        &ctx,
        &[(
            idwt(&format!("{M4B}, \"haar\", center")),
            &[
                0.707_106_781_187,
                -3.414_213_562_37,
                -0.292_893_218_813,
                -1.414_213_562_37,
                -4.181_980_515_34,
                8.803_300_858_9,
                0.025_126_265_847_1,
                3.181_980_515_34,
                -3.181_980_515_34,
                4.974_873_734_15,
                -1.803_300_858_9,
                2.181_980_515_34,
                1.414_213_562_37,
                -1.707_106_781_19,
                -0.585_786_437_627,
                -0.707_106_781_187,
            ],
        )],
    );
}

/// giac's decimal rows (p2.out, p1.out), every wavelet family, centred,
/// inverse, 2-D standard and non-standard.
#[test]
#[allow(clippy::too_many_lines, reason = "one row per giac answer")]
fn decimals() {
    let ctx = Context::new();
    let daubechies4: &[f64] = &[
        5.0,
        -1.732_050_807_57,
        2.220_446_049_25e-16,
        -1.414_213_562_37,
    ];
    let cases: Vec<(String, &[f64])> = vec![
        (dwt(V4), daubechies4),
        (dwt(&format!("{V4}, 4")), daubechies4),
        (dwt(&format!("{V4}, \"daubechies\", 4")), daubechies4),
        (
            dwt(&format!("{V4}, \"bspline\"")),
            &[5.0, 2.0, -0.707_106_781_187, -0.707_106_781_187],
        ),
        (
            dwt("List[0x1.8p+0, 2, 3, 4], \"haar\""),
            &[5.25, -1.75, -0.353_553_390_593, -0.707_106_781_187],
        ),
        (
            dwt(&format!("{V8}, 6")),
            &[
                12.727_922_061_4,
                -1.479_388_396_55,
                -4.408_985_905_29,
                2.246_708_245_13,
                0.0,
                0.0,
                -3.793_767_650_89,
                0.965_340_526_142,
            ],
        ),
        (
            dwt("List[0, 1/2, 1/4, 1, 3/4, -1/4, -1/2, 0]"),
            &[
                0.618_718_433_538,
                1.055_475_572_15,
                -0.121_935_118_395,
                -0.207_311_706_132,
                -0.385_905_771_231,
                0.491_631_707_41,
                -0.459_279_326_772,
                -0.176_776_695_297,
            ],
        ),
        (
            dwt("List[1, 2, 3]"),
            &[3.0, 1.0, 1.931_851_652_58, -0.517_638_090_205],
        ),
        (
            dwt("List[1/2, 1/3, 1/4, 1/5]"),
            &[
                0.641_666_666_667,
                0.136_821_535_725,
                -0.026_882_890_650_8,
                0.180_089_359_908,
            ],
        ),
        (
            dwt(&format!("{V4}, center")),
            &[
                5.0,
                1.732_050_807_57,
                -1.414_213_562_37,
                2.220_446_049_25e-16,
            ],
        ),
        (
            dwt(&format!("{V4}, center, ColumnSpace")),
            &[
                5.0,
                1.732_050_807_57,
                -1.414_213_562_37,
                2.220_446_049_25e-16,
            ],
        ),
        (
            dwt(&format!("{V8}, center")),
            &[
                12.727_922_061_4,
                2.967_127_832_99,
                4.830_127_018_92,
                -1.366_025_403_78,
                -2.828_427_124_75,
                2.220_446_049_25e-16,
                4.440_892_098_5e-16,
                4.440_892_098_5e-16,
            ],
        ),
        (
            dwt(&format!("{V8}, \"bspline\", center")),
            &[
                12.727_922_061_4,
                -3.535_533_905_93,
                2.5,
                -2.5,
                4.949_747_468_31,
                -0.707_106_781_187,
                -0.707_106_781_187,
                -0.707_106_781_187,
            ],
        ),
        (
            dwt(&format!("{V8}, \"bspline\", 309")),
            &[
                12.727_922_061_4,
                -9.222_771_649_93,
                -1.630_371_093_75,
                1.630_371_093_75,
                0.0,
                2.220_446_049_25e-16,
                2.220_446_049_25e-16,
                2.828_427_124_75,
            ],
        ),
        (
            dwt(&format!("{V8}, \"bspline\", 105")),
            &[
                12.727_922_061_4,
                5.656_854_249_49,
                -2.0,
                -2.0,
                -0.707_106_781_187,
                -0.707_106_781_187,
                -0.707_106_781_187,
                -0.707_106_781_187,
            ],
        ),
        (
            dwt(&format!("{V8}, \"bspline\", 208, center")),
            &[
                12.727_922_061_4,
                7.919_043_522_12,
                3.200_195_312_5,
                0.799_804_687_5,
                -2.828_427_124_75,
                2.220_446_049_25e-16,
                0.0,
                -4.440_892_098_5e-16,
            ],
        ),
        (
            dwt(&format!("{V8}, \"bspline\", 301")),
            &[
                12.727_922_061_4,
                -5.656_854_249_49,
                -1.0,
                1.0,
                0.0,
                2.220_446_049_25e-16,
                2.220_446_049_25e-16,
                2.828_427_124_75,
            ],
        ),
        (
            dwt(&format!("{V8}, \"daubechies\", 20")),
            &[
                12.727_922_061_4,
                -3.445_470_555_71,
                -1.496_572_443_6,
                4.026_316_147_39,
                -2.688_357_849_88,
                1.674_510_126_43,
                -0.916_791_219_661,
                -0.897_788_181_637,
            ],
        ),
        (
            idwt(&format!("{V8}, \"bspline\", 202, center")),
            &[
                2.578_427_124_75,
                -5.992_640_687_12,
                4.492_640_687_12,
                -1.992_640_687_12,
                5.906_854_249_49,
                -5.992_640_687_12,
                6.821_067_811_87,
                -2.992_640_687_12,
            ],
        ),
        (
            idwt(&format!("{V8}, \"daubechies\", 8, center")),
            &[
                5.559_128_339_54,
                -1.503_625_989_36,
                5.124_201_182_38,
                -6.457_676_990_54,
                2.053_709_545_44,
                -2.479_300_497_91,
                7.061_950_805_86,
                -6.529_959_270_66,
            ],
        ),
        (
            dwt(M44),
            &[
                103.75,
                -7.238_620_667_98,
                -2.324_725_813_08,
                1.971_172_422_48,
                -6.568_747_686_9,
                15.376_760_388_8,
                -21.396_659_667_3,
                -0.404_199_289_747,
                23.244_435_026_5,
                12.385_360_261_5,
                20.773_316_684_9,
                39.776_193_115_5,
                11.757_350_642_2,
                -31.676_003_230_2,
                -35.149_432_726_7,
                -13.900_077_073_7,
            ],
        ),
        (
            dwt("List[List[1, 2], List[3, 4]]"),
            &[5.0, -1.0, -2.0, -1.665_334_536_94e-16],
        ),
        (
            dwt("List[List[1, 2], List[3, 4]], ColumnSpace"),
            &[5.0, -1.0, -2.0, -1.665_334_536_94e-16],
        ),
        (
            dwt("List[List[1, 2, 3], List[4, 5, 6]]"),
            &[
                5.25,
                2.024_519_052_84,
                2.656_296_022_29,
                -0.181_422_288_142,
                5.671_633_369_87,
                2.265_544_456_62,
                2.662_642_020_14,
                0.010_991_590_687_8,
                -2.069_307_577_97,
                -0.812_832_717_365,
                -1.007_772_228_31,
                0.032_291_281_149_7,
                -1.112_672_937_37,
                -0.636_056_022_068,
                -0.016_746_824_526_9,
                -0.507_772_228_311,
            ],
        ),
        (
            dwt(&format!("{M4B}, ColumnSpace")),
            &[
                10.0,
                -1.549_038_105_68,
                -0.091_506_350_946_1,
                -1.774_519_052_84,
                6.196_152_422_71,
                -2.0,
                -0.591_506_350_946,
                0.457_531_754_731,
                -2.792_468_245_27,
                -2.475_480_947_16,
                -0.341_506_350_946,
                0.841_506_350_946,
                -1.890_544_456_62,
                -0.841_506_350_946,
                0.158_493_649_054,
                -0.658_493_649_054,
            ],
        ),
        (
            idwt(&format!("{M4B}, ColumnSpace")),
            &[
                5.137_259_526_42,
                0.152_803_983_042,
                3.911_778_579_26,
                -2.188_702_367_9,
                2.920_753_175_47,
                2.643_829_386_83,
                0.896_234_122_635,
                -3.045_753_175_47,
                3.844_791_281_15,
                -2.719_791_281_15,
                2.106_170_613_17,
                -3.244_310_333_99,
                1.195_272_228_31,
                -6.710_816_684_93,
                1.987_740_473_58,
                -2.887_259_526_42,
            ],
        ),
        (
            dwt(&format!("{M4B}, \"bspline\", ColumnSpace, center")),
            &[
                10.0,
                -0.5,
                1.0,
                -4.163_336_342_34e-17,
                2.0,
                -0.5,
                1.5,
                -0.5,
                -1.0,
                -1.0,
                -2.0,
                1.0,
                4.5,
                5.5,
                1.5,
                -0.5,
            ],
        ),
        (dwt("List[5]"), &[5.0]),
    ];
    near(&ctx, &cases);
}

/// giac's rounding noise (`2.2·10⁻¹⁶`) is about `10⁻³⁰` with GSL's
/// literals read exactly.
#[test]
fn noise() {
    let ctx = Context::new();
    let found = doubles(&ctx, &dwt(V4));
    assert!(found[2].abs() < 1e-25, "{found:?}");
}

/// `idwt(dwt(x))` is `x` within GSL's literals' error, about `10⁻²⁸`.
#[test]
fn round_trips() {
    let ctx = Context::new();
    for (x, options) in [
        (V8, ""),
        (V8, ", \"bspline\", 309, center"),
        (M44, ""),
        (M4B, ", ColumnSpace"),
        (M4B, ", \"daubechies\", 20"),
        (M4B, ", \"bspline\", 105, center, ColumnSpace"),
    ] {
        let back = idwt(&format!("{}{options}", dwt(&format!("{x}{options}"))));
        let found = doubles(&ctx, &back);
        let expected = doubles(&ctx, x);
        assert_eq!(found.len(), expected.len(), "{back}");
        for (f, e) in found.iter().zip(&expected) {
            assert!((f - e).abs() < 1e-12, "{back}: {found:?}");
        }
    }
}

/// Padding: giac's warning, then the answer.
#[test]
fn padding() {
    let ctx = Context::new();
    let (found, warnings) = with_warnings(&ctx, &dwt("List[1, 2, 3]"));
    assert_eq!(found.args().len(), 4);
    assert_eq!(
        warnings,
        ["data size is not a power of 2, padding with zeros up to size 4"]
    );
    let (found, warnings) = with_warnings(&ctx, &dwt("List[List[1, 2, 3], List[4, 5, 6]]"));
    assert_eq!(found.args().len(), 4);
    assert_eq!(
        warnings,
        ["data size is not a power of 2, padding with zeros up to size 4×4"]
    );
    let (_, warnings) = with_warnings(&ctx, &dwt(V8));
    assert!(warnings.is_empty(), "{warnings:?}");
}

/// giac's errors: the call stays with its message.
#[test]
fn stays() {
    let ctx = Context::new();
    let cases = [
        (dwt(&format!("{V4}, 2")), "Wavelet not implemented"),
        (
            dwt(&format!("{V4}, \"haar\", 4")),
            "Wavelet not implemented",
        ),
        (
            dwt(&format!("{V4}, \"bspline\", 104")),
            "Wavelet not implemented",
        ),
        (dwt(&format!("{V4}, 22")), "Wavelet not implemented"),
        (dwt(&format!("{V4}, haar")), "Invalid argument 2"),
        (dwt(&format!("{V4}, center, 1/2")), "Invalid argument 3"),
        (
            dwt(&format!("{V4}, Equation[centre, True]")),
            "Invalid argument 2",
        ),
        (dwt(&format!("{V4}, \"foo\"")), "Unsupported wavelet type"),
        (dwt("List[x, 1, 2, 3]"), "Data must be numeric"),
        (dwt("List[x, 1, 2, 3], \"haar\""), "Data must be numeric"),
        (dwt("List[List[1, 2], List[3]]"), "Data must be numeric"),
        (
            idwt("List[List[1, x], List[3, 4]]"),
            "Failed to convert data to GSL matrix",
        ),
        (
            dwt(&format!("{V4}, \"haar\", \"daubechies\"")),
            "Failed to allocate memory",
        ),
    ];
    for (input, message) in &cases {
        let (found, warnings) = with_warnings(&ctx, input);
        assert_eq!(found, stayed(&ctx, input), "{input}");
        assert_eq!(warnings, [*message], "{input}");
    }
}

/// Malformed calls stay without a warning.
#[test]
fn stays_silently() {
    let ctx = Context::new();
    for input in [
        "DiscreteWaveletTransform[]".to_owned(),
        dwt("5"),
        idwt("x"),
        dwt(&format!("{V4}, 4, 4, 4, 4, 4")),
    ] {
        let (found, warnings) = with_warnings(&ctx, &input);
        assert_eq!(found, stayed(&ctx, &input), "{input}");
        assert!(warnings.is_empty(), "{input}: {warnings:?}");
    }
}

/// Sizes past `Limits::terms` (padded entries times the filter length)
/// stay before anything is built.
#[test]
fn limits() {
    let mut small = Context::new();
    let limits = small.limits();
    small.set_limits(Limits {
        terms: 20,
        ..limits
    });
    for input in [
        dwt(V8),
        dwt(&format!("{M4B}, \"haar\"")),
        dwt("List[1, 2, 3, 4, 5, 6, 7, 8, 9], \"haar\""),
    ] {
        let (found, warnings) = with_warnings(&small, &input);
        assert_eq!(found, stayed(&small, &input), "{input}");
        assert!(warnings.is_empty(), "{input}: {warnings:?}");
    }
    check(
        &small,
        &[(
            dwt(&format!("{V4}, \"haar\"")),
            &format!("List[5, -2, Times[-1/2, {R2}], Times[-1/2, {R2}]]"),
        )],
    );
}
