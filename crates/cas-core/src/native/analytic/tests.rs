//! Tests of the analytic signal builtins (CAS-635): hilbert's list form,
//! instfreq and instphase, emd; giac's values, exact forms, what stays.

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

/// The inputs stay unevaluated.
fn stays(inputs: &[&str]) {
    let ctx = Context::new();
    for input in inputs {
        let e = read(input);
        let args = ctx
            .eval(&Expr::call(Symbol::LIST, e.args().to_vec()))
            .unwrap();
        let expected = Expr::apply(e.head().unwrap(), args.args());
        assert_eq!(ctx.eval(&e).unwrap(), expected, "{input}");
    }
}

/// A real number's double, by `N`.
fn double(ctx: &Context, x: &Expr) -> f64 {
    let v = ctx.eval(&Expr::call(Symbol::N, [x.clone()])).unwrap();
    if let View::Ball(b) = v.view() {
        b.mid().to_f64()
    } else {
        let Some(q) = v.to_rational() else {
            panic!("{x:?}: {v:?}");
        };
        cas_num::ball::Float::from_rational(&q, 64, cas_num::ball::Round::Nearest)
            .0
            .to_f64()
    }
}

/// The parts `(re, im)` of a list's entries, as doubles (C0010, C0013).
fn parts(ctx: &Context, e: &Expr) -> Vec<(f64, f64)> {
    e.args()
        .iter()
        .map(|x| {
            let re = ctx.eval(&Expr::call(Symbol::RE, [x.clone()])).unwrap();
            let im = ctx.eval(&Expr::call(Symbol::IM, [x.clone()])).unwrap();
            (double(ctx, &re), double(ctx, &im))
        })
        .collect()
}

/// `input` evaluates to a list whose entries are `expected` within `tol`.
fn near(input: &str, expected: &[(f64, f64)], tol: f64) {
    let ctx = Context::new();
    let found = ctx.eval(&read(input)).unwrap();
    assert_eq!(
        found.head().and_then(Expr::as_symbol),
        Some(Symbol::LIST),
        "{input}: {found:?}"
    );
    let found = parts(&ctx, &found);
    assert_eq!(found.len(), expected.len(), "{input}: {found:?}");
    for (k, ((a, b), (c, d))) in found.iter().zip(expected).enumerate() {
        assert!(
            (a - c).abs() <= tol && (b - d).abs() <= tol,
            "{input}[{k}]: {a}+{b}i, expected {c}+{d}i"
        );
    }
}

fn real(values: &[f64]) -> Vec<(f64, f64)> {
    values.iter().map(|v| (*v, 0.0)).collect()
}

const L: &str = "List[0x1.6666666666666p-1, 0x1.999999999999ap-3, -0x1.0000000000000p-1, \
                 -0x1.3333333333333p-1, -0x1.999999999999ap-3, 0x1.3333333333333p-2, \
                 0x1.0000000000000p-1, -0x1.999999999999ap-4, -0x1.3333333333333p-2]";

/// giac's `hilbert` on lists; the integer form stays the matrix.
#[test]
fn hilbert_lists() {
    near(
        "HilbertMatrix[List[1, 0, 0, 0]]",
        &[(1.0, 0.0), (0.0, 0.5), (0.0, 0.0), (0.0, -0.5)],
        1e-15,
    );
    near(
        &format!("HilbertMatrix[{L}]"),
        &[
            (0.7, -0.166_169_214_531),
            (0.2, 0.744_636_918_643),
            (-0.5, 0.408_186_837_669),
            (-0.6, -0.127_861_937_19),
            (-0.2, -0.562_185_310_255),
            (0.3, -0.412_752_500_28),
            (0.5, 0.063_091_044_045),
            (-0.1, 0.452_633_687_721),
            (-0.3, -0.399_579_525_821),
        ],
        1e-11,
    );
    let third = 0.577_350_269_19;
    near(
        "HilbertMatrix[List[1, 2, 3]]",
        &[(1.0, third), (2.0, -2.0 * third), (3.0, third)],
        1e-10,
    );
    // The imaginary part is dropped, with giac's warning.
    near(
        "HilbertMatrix[List[1, Complex[0, 1], 0]]",
        &[(1.0, 0.0), (0.0, third), (0.0, -third)],
        1e-10,
    );
    check(&[("HilbertMatrix[2]", "List[List[1, 1/2], List[1/2, 1/3]]")]);
    stays(&["HilbertMatrix[List[]]", "HilbertMatrix[List[x, 1]]"]);
}

/// The list forms, exact on exact samples, giac's doubles on decimals.
#[test]
fn instantaneous_lists() {
    check(&[
        (
            "InstantaneousFrequency[List[1, Complex[0, 1], -1, Complex[0, -1]]]",
            "List[Times[1/2, Pi], Times[1/2, Pi], Times[1/2, Pi]]",
        ),
        (
            "InstantaneousPhase[List[1, Complex[0, 1], -1, Complex[0, -1]]]",
            "List[0, Times[1/2, Pi], Pi, Times[3/2, Pi]]",
        ),
        ("InstantaneousFrequency[List[1, 2]]", "List[0]"),
        ("InstantaneousFrequency[List[1]]", "List[]"),
        ("InstantaneousPhase[List[]]", "List[]"),
    ]);
    near(
        &format!("InstantaneousFrequency[HilbertMatrix[{L}]]"),
        &real(&[
            1.541_472_248_07,
            1.148_545_502_45,
            0.894_607_943_316,
            1.019_041_669_82,
            0.970_297_400_326,
            1.067_809_822_17,
            1.662_714_037_98,
            2.280_150_198_13,
        ]),
        1e-10,
    );
    near(
        &format!("InstantaneousPhase[HilbertMatrix[{L}]]"),
        &real(&[
            -0.233_070_552_139,
            1.308_401_695_93,
            2.456_947_198_38,
            3.351_555_141_7,
            4.370_596_811_51,
            5.340_894_211_84,
            6.408_704_034_01,
            8.071_418_071_99,
            10.351_568_270_1,
        ]),
        1e-9,
    );
    near(
        "InstantaneousPhase[List[0x1p+0, -0x1p+0, 0x1p+0]]",
        &real(&[0.0, std::f64::consts::PI, 2.0 * std::f64::consts::PI]),
        1e-12,
    );
}

/// The expression forms against giac's, by the simplified difference.
#[test]
fn instantaneous_expressions() {
    let ctx = Context::new();
    let cases = [
        (
            "InstantaneousFrequency[Power[E, Times[Complex[0, 1], Power[x, 2]]], x]",
            "Times[2, x]",
        ),
        (
            "InstantaneousPhase[Power[E, Times[Complex[0, 1], Power[x, 2]]], x]",
            "Power[x, 2]",
        ),
        (
            "InstantaneousFrequency[Power[E, Times[Complex[0, 2], x]]]",
            "2",
        ),
        (
            "InstantaneousPhase[Power[E, Times[Complex[0, 2], x]]]",
            "Times[2, x]",
        ),
        (
            "InstantaneousFrequency[Plus[Cos[t], Times[Complex[0, 1], Sin[t]]], t]",
            "1",
        ),
    ];
    for (input, expected) in cases {
        let found = ctx.eval(&read(input)).unwrap();
        let gap = Expr::call(
            Symbol::SIMPLIFY,
            [Expr::call(
                Symbol::PLUS,
                [
                    found.clone(),
                    Expr::call(Symbol::TIMES, [Expr::from(-1), read(expected)]),
                ],
            )],
        );
        assert_eq!(ctx.eval(&gap).unwrap(), Expr::from(0), "{input}: {found:?}");
    }
    stays(&["InstantaneousFrequency[3]", "InstantaneousPhase[x]"]);
}

const DATA: &str = "List[1, 2, 1, 3, 1, 2, 1, 3, 1]";

const IMF: [f64; 9] = [
    0.0,
    0.508_370_535_714,
    -0.875,
    0.978_794_642_857,
    -0.8125,
    0.518_725_198_413,
    -0.263_888_888_889,
    1.881_634_424_6,
    0.0,
];

const RESIDUE: [f64; 9] = [
    1.0,
    1.491_629_464_29,
    1.875,
    2.021_205_357_14,
    1.8125,
    1.481_274_801_59,
    1.263_888_888_89,
    1.118_365_575_4,
    1.0,
];

/// The entries of `e` as doubles.
fn doubles(ctx: &Context, e: &Expr) -> Vec<f64> {
    e.args().iter().map(|x| double(ctx, x)).collect()
}

fn close(found: &[f64], expected: &[f64], what: &str) {
    assert_eq!(found.len(), expected.len(), "{what}: {found:?}");
    for (a, b) in found.iter().zip(expected) {
        assert!((a - b).abs() <= 1e-10, "{what}: {found:?}");
    }
}

/// giac's `emd` answers, the options, and what stays.
#[test]
fn empirical_mode_decompositions() {
    let ctx = Context::new();
    let found = ctx
        .eval(&read(&format!("EmpiricalModeDecomposition[{DATA}]")))
        .unwrap();
    assert_eq!(
        found.head().and_then(Expr::as_symbol),
        Some(Symbol::SEQUENCE),
        "{found:?}"
    );
    let [imfs, residue] = found.args() else {
        panic!("{found:?}");
    };
    assert_eq!(imfs.args().len(), 1, "{found:?}");
    close(&doubles(&ctx, &imfs.args()[0]), &IMF, "imf");
    close(&doubles(&ctx, residue), &RESIDUE, "residue");
    // residue = false: the IMFs alone; the same with threshold, limit and
    // extrema giac's rows leave unchanged.
    for options in [
        "Equation[Residue, False]",
        "Equal[Residue, False], Equation[Threshold, 1/2], Equation[Limit, 1], Equation[Extrema, 6]",
    ] {
        let found = ctx
            .eval(&read(&format!(
                "EmpiricalModeDecomposition[{DATA}, {options}]"
            )))
            .unwrap();
        assert_eq!(
            found.head().and_then(Expr::as_symbol),
            Some(Symbol::LIST),
            "{found:?}"
        );
        assert_eq!(found.args().len(), 1, "{found:?}");
        close(&doubles(&ctx, &found.args()[0]), &IMF, options);
    }
    // No IMF: the data is the residue.
    for input in [
        "EmpiricalModeDecomposition[List[1, 2, 3, 4]]".to_owned(),
        format!("EmpiricalModeDecomposition[{DATA}, Equation[Extrema, 7]]"),
        "EmpiricalModeDecomposition[List[5]]".to_owned(),
    ] {
        let found = ctx.eval(&read(&input)).unwrap();
        let [imfs, residue] = found.args() else {
            panic!("{input}: {found:?}");
        };
        assert!(imfs.args().is_empty(), "{input}: {found:?}");
        let data = read(&input).args()[0].clone();
        close(&doubles(&ctx, residue), &doubles(&ctx, &data), &input);
    }
    // threshold = 0 ends: an answer or the sift bound, never a hang.
    ctx.eval(&read(&format!(
        "EmpiricalModeDecomposition[{DATA}, Equation[Threshold, 0]]"
    )))
    .unwrap();
    stays(&[
        "EmpiricalModeDecomposition[List[]]",
        "EmpiricalModeDecomposition[List[x, 1, 2]]",
        "EmpiricalModeDecomposition[List[1, 2, 1, 3, 1], Equation[Threshold, -1]]",
        "EmpiricalModeDecomposition[List[1, 2, 1, 3, 1], Equation[Residue, 0]]",
        "EmpiricalModeDecomposition[List[1, 2, 1, 3, 1], Equation[Extrema, 1]]",
        "EmpiricalModeDecomposition[List[1, 2, 1, 3, 1], Equation[Limit, 0]]",
        "EmpiricalModeDecomposition[3]",
    ]);
}

/// Two IMFs of giac's 13-sample row, and the residue.
#[test]
fn two_imfs() {
    let ctx = Context::new();
    let found = ctx
        .eval(&read(
            "EmpiricalModeDecomposition[List[0, 1, 0, -1, 0, 1, 0, -1, 0, 1, 0, -1, 0]]",
        ))
        .unwrap();
    let [imfs, residue] = found.args() else {
        panic!("{found:?}");
    };
    assert_eq!(imfs.args().len(), 2, "{found:?}");
    close(
        &doubles(&ctx, &imfs.args()[1]),
        &[
            0.0,
            0.205_319_250_132,
            0.041_027_771_578_8,
            -0.312_767_397_451,
            -0.590_114_443_035,
            -0.665_691_321_33,
            -0.490_138_026_619,
            -0.182_006_428_017,
            0.111_264_243_165,
            0.294_597_181_678,
            0.308_870_902_19,
            0.188_615_913_522,
            0.0,
        ],
        "imf 2",
    );
    close(
        &doubles(&ctx, residue),
        &[
            0.0,
            0.278_599_007_121,
            0.447_196_540_561,
            0.497_305_828_184,
            0.463_419_369_83,
            0.384_723_157_766,
            0.295_082_776_942,
            0.207_082_183_045,
            0.127_984_924_447,
            0.065_054_549_519_7,
            0.024_504_418_254_5,
            0.006_153_996_594_49,
            0.0,
        ],
        "residue",
    );
}

/// `input` evaluates to a matrix whose rows are `expected` within `tol`.
fn near_rows(input: &str, expected: &[&[(f64, f64)]], tol: f64) {
    let ctx = Context::new();
    let found = ctx.eval(&read(input)).unwrap();
    assert_eq!(found.args().len(), expected.len(), "{input}: {found:?}");
    for (k, (row, want)) in found.args().iter().zip(expected).enumerate() {
        let got = parts(&ctx, row);
        assert_eq!(got.len(), want.len(), "{input}[{k}]: {got:?}");
        for ((a, b), (c, d)) in got.iter().zip(*want) {
            assert!(
                (a - c).abs() <= tol && (b - d).abs() <= tol,
                "{input}[{k}]: {got:?}, expected {want:?}"
            );
        }
    }
}

const EIGHT: &str = "List[1, 2, 3, 4, 5, 6, 7, 8]";

/// giac's `stft` and `istft` (`proto_stft.py`, `giac/p6.out`), and what stays.
#[test]
fn short_time_fourier_transforms() {
    let tail = (-0.16, 0.16);
    near_rows(
        &format!("ShortTimeFourierTransform[{EIGHT}, 2]"),
        &[
            &[(0.8, 0.0), tail],
            &[(1.44, 0.0), tail],
            &[(2.08, 0.0), tail],
        ],
        1e-12,
    );
    near_rows(
        &format!("ShortTimeFourierTransform[{EIGHT}, 4, HannWindow]"),
        &[&[
            (13.5, 0.0),
            (-9.053_300_858_899_11, -1.189_339_828_220_18),
            (1.5, -1.5),
            (1.553_300_858_899_11, 3.310_660_171_779_82),
        ]],
        1e-12,
    );
    // Seven samples padded to eight; the imaginary part is dropped.
    near_rows(
        "ShortTimeFourierTransform[List[1, 2, 3, 4, 5, 6, 7], 2]",
        &[
            &[(0.8, 0.0), tail],
            &[(1.44, 0.0), tail],
            &[(1.44, 0.0), (-0.16, -0.48)],
        ],
        1e-12,
    );
    near_rows(
        "ShortTimeFourierTransform[List[1, Complex[2, 1], 3, 4], 2]",
        &[&[(0.8, 0.0), tail]],
        1e-12,
    );
    check(&[
        (&format!("ShortTimeFourierTransform[{EIGHT}]"), "List[]"),
        (&format!("ShortTimeFourierTransform[{EIGHT}, 8]"), "List[]"),
    ]);
    near(
        "InverseShortTimeFourierTransform[List[List[1, 2], List[3, 4]]]",
        &real(&[1.5, 0.0, 3.0, 0.0, -0.5, 0.0]),
        1e-12,
    );
    near(
        "InverseShortTimeFourierTransform[List[List[1, 2, 3]]]",
        &real(&[2.0, -1.0 / 6.0, -0.5, 1.0 / 3.0, -0.5, -1.0 / 6.0]),
        1e-12,
    );
    let s = 0.220_970_869_120_796;
    near(
        "InverseShortTimeFourierTransform[List[List[1/2, -1, 2, 1/4]]]",
        &real(&[0.4375, -s, -0.375, s, 0.8125, s, -0.375, -s]),
        1e-12,
    );
    stays(&[
        &format!("ShortTimeFourierTransform[{EIGHT}, 3]"),
        &format!("ShortTimeFourierTransform[{EIGHT}, 1]"),
        &format!("ShortTimeFourierTransform[{EIGHT}, 0]"),
        &format!("ShortTimeFourierTransform[{EIGHT}, 2, 4]"),
        &format!("ShortTimeFourierTransform[{EIGHT}, 2, Sin]"),
        &format!("ShortTimeFourierTransform[{EIGHT}, 2, HannWindow, 1]"),
        "ShortTimeFourierTransform[List[x, 1, 2, 3], 2]",
        "ShortTimeFourierTransform[List[]]",
        "ShortTimeFourierTransform[3]",
        "InverseShortTimeFourierTransform[List[List[1, 2], List[3]]]",
        "InverseShortTimeFourierTransform[List[List[1, x], List[3, 4]]]",
        "InverseShortTimeFourierTransform[List[1, 2]]",
    ]);
}

/// A matrix cell: row, column, value.
type Cell = (usize, usize, f64);

/// The nonzero cells of a matrix, and its size.
fn cells(input: &str) -> ((usize, usize), Vec<Cell>) {
    let ctx = Context::new();
    let found = ctx.eval(&read(input)).unwrap();
    let rows = found.args();
    let n = rows.first().map_or(0, |r| r.args().len());
    let mut out = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        assert_eq!(row.args().len(), n, "{input}: {found:?}");
        for (j, v) in doubles(&ctx, row).into_iter().enumerate() {
            if v != 0.0 {
                out.push((i, j, v));
            }
        }
    }
    ((rows.len(), n), out)
}

fn same_cells(input: &str, size: (usize, usize), expected: &[Cell]) {
    let (found_size, found) = cells(input);
    assert_eq!(found_size, size, "{input}");
    assert_eq!(found.len(), expected.len(), "{input}: {found:?}");
    for ((i, j, v), (k, l, w)) in found.iter().zip(expected) {
        assert!(
            i == k && j == l && (v - w).abs() <= 1e-12,
            "{input}: {found:?}"
        );
    }
}

const MATRIX: &str = "Equation[output, Matrix]";

/// giac's `hht` matrices (`proto_hht.py`, `giac/p7.out`), and what stays.
#[test]
fn hilbert_huang_transforms() {
    // giac's top of the range lands one bin low: 98 of 100.
    same_cells(
        &format!("HilbertHuangTransform[List[List[1, -1, 1, -1]], {MATRIX}]"),
        (100, 4),
        &[(0, 0, 1.0), (98, 1, 1.0), (98, 2, 1.0), (98, 3, 1.0)],
    );
    same_cells(
        &format!(
            "HilbertHuangTransform[List[List[1, 0, -1, 0]], {MATRIX}, Equation[SampleRate, 8], \
             Equation[Frequencies, Range[0, 4]]]"
        ),
        (100, 4),
        &[(0, 0, 1.0), (49, 1, 1.0), (49, 2, 1.0), (49, 3, 1.0)],
    );
    // Ours honours bins=n (giac's 161 rows).
    same_cells(
        &format!(
            "HilbertHuangTransform[List[List[1, -1, 1, -1]], {MATRIX}, Equation[Bins, 5], \
             Equation[Frequencies, Range[0, 4]]]"
        ),
        (5, 4),
        &[(0, 0, 1.0), (3, 1, 1.0), (3, 2, 1.0), (3, 3, 1.0)],
    );
    // A range that is not positive, or not numeric: zeros.
    same_cells(
        &format!(
            "HilbertHuangTransform[List[List[1, 0, -1, 0]], {MATRIX}, \
             Equation[Frequencies, Range[a, b]]]"
        ),
        (100, 4),
        &[],
    );
    same_cells(
        &format!(
            "HilbertHuangTransform[List[List[1, 0, -1, 0]], {MATRIX}, \
             Equation[Frequencies, Range[2, 2]]]"
        ),
        (100, 4),
        &[],
    );
    let ones = "List[List[1, 0, -1, 0]]";
    stays(&[
        "HilbertHuangTransform[3]",
        &format!("HilbertHuangTransform[{ones}]"),
        &format!("HilbertHuangTransform[{ones}, Equation[output, foo]]"),
        &format!("HilbertHuangTransform[{ones}, Log, {MATRIX}]"),
        &format!("HilbertHuangTransform[{ones}, {MATRIX}, Equation[SampleRate, 0]]"),
        &format!("HilbertHuangTransform[{ones}, {MATRIX}, Equation[Bins, 0]]"),
        &format!("HilbertHuangTransform[{ones}, {MATRIX}, Equation[foo, 3]]"),
        &format!(
            "HilbertHuangTransform[{ones}, {MATRIX}, Equation[Log, 3], Equation[SampleRate, 8]]"
        ),
        &format!("HilbertHuangTransform[List[1, 2, 3], {MATRIX}]"),
    ]);
}
