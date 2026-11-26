//! Tests of clustering (CAS-630): `cluster`'s labels and `kmeans`'
//! partitions (up to relabelling) against the prototype's answers
//! (crates/cas-core/proto/clustering/proto.py, R1.7-012), the seed-dependent ones among the
//! partitions it found over 200 seeds, a run repeated under `SeedRandom`,
//! each warning, and the limits; `output=part`, strings' Levenshtein
//! distance and `cluster`'s silhouette (1.8/T5, CAS-646: giac's answers in
//! crates/cas-core/proto/giac/probe_clustering2.out, the prototype's `t5`);
//! `kmeans`' `index=` (T5.c3: `probe_clustering3.out`, the prototype's `t5c3`),
//! its seed-independent values and counts.

use crate::context::Limits;
use crate::trace::Warning;
use crate::{Context, EvalOptions, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

fn eval(ctx: &Context, text: &str) -> Expr {
    ctx.eval(&read(text))
        .unwrap_or_else(|e| panic!("{text}: {e}"))
}

/// giac's D6, its decimals as the oracle reads them.
const D6: &str = "N[List[List[11/10, 3/2], List[16/5, 5/2], List[3/2, 9/10], \
                  List[7/2, 27/10], List[4/5, 11/10], List[31/10, 11/5]]]";
const E6: &str = "List[List[0, 0], List[0, 1], List[5, 5], List[5, 6], List[10, 0], List[10, 1]]";
const P8: &str = "List[List[0, 0], List[1, 0], List[0, 2], List[4, 4], List[5, 5], List[4, 6], \
                  List[9, 1], List[10, 0]]";

/// The labels of an answer.
fn labels(input: &str, found: &Expr) -> Vec<i64> {
    found
        .args()
        .iter()
        .map(|l| l.as_i64().unwrap_or_else(|| panic!("{input}: {found:?}")))
        .collect()
}

/// Labels renumbered by first appearance.
fn canonical(labels: &[i64]) -> Vec<i64> {
    let mut seen: Vec<i64> = Vec::new();
    labels
        .iter()
        .map(|l| {
            let i = seen.iter().position(|s| s == l).unwrap_or_else(|| {
                seen.push(*l);
                seen.len() - 1
            });
            i64::try_from(i).expect("small")
        })
        .collect()
}

/// `cluster`'s labels as they are (probe4, probe8 and the prototype's
/// extra rows).
#[test]
fn cluster_labels() {
    let ctx = Context::new();
    let alike = [0, 1, 0, 1, 0, 1];
    let mut rows: Vec<(String, Vec<i64>)> = vec![
        (format!("Cluster[{D6}]"), alike.to_vec()),
        (format!("Cluster[{D6}, 2]"), alike.to_vec()),
        (format!("Cluster[{D6}, Equation[Count, 2]]"), alike.to_vec()),
        (format!("Cluster[{D6}, 3]"), vec![0, 1, 2, 1, 0, 1]),
        (
            format!("Cluster[{D6}, Equation[Count, 3]]"),
            vec![0, 1, 2, 1, 0, 1],
        ),
        (
            format!("Cluster[{D6}, Equation[Count, Range[2, 4]]]"),
            alike.to_vec(),
        ),
        (
            format!("Cluster[{D6}, Equation[output, list]]"),
            alike.to_vec(),
        ),
        (format!("Cluster[{E6}]"), vec![0, 0, 1, 1, 2, 2]),
        (format!("Cluster[{E6}, 2]"), vec![0, 0, 0, 0, 1, 1]),
        (
            format!("Cluster[{E6}, Equation[Type, \"single\"], 3]"),
            vec![0, 0, 1, 1, 2, 2],
        ),
        (
            "Cluster[List[1, 2, 10, 11, 20, 21]]".to_owned(),
            vec![0, 0, 1, 1, 2, 2],
        ),
        (
            "Cluster[List[1, 2, 10, 11, 20, 21], 3]".to_owned(),
            vec![0, 0, 1, 1, 2, 2],
        ),
        ("Cluster[List[1/2, 1, 3, 7/2]]".to_owned(), vec![0, 0, 1, 1]),
        ("Cluster[List[List[1, 2]]]".to_owned(), vec![0]),
        (
            "Cluster[List[List[1, 2], List[3, 4]]]".to_owned(),
            vec![0, 1],
        ),
        (format!("Cluster[{P8}]"), vec![0, 0, 0, 1, 1, 1, 2, 2]),
        (
            format!("Cluster[{P8}, Equation[Type, \"ward\"]]"),
            vec![0, 0, 0, 1, 1, 1, 2, 2],
        ),
        (
            format!("Cluster[{P8}, Equation[Type, \"ward\"], Equation[Count, 3]]"),
            vec![0, 0, 0, 1, 1, 1, 2, 2],
        ),
        (
            format!("Cluster[{P8}, Equation[Type, \"complete\"]]"),
            vec![0, 0, 0, 0, 0, 0, 1, 1],
        ),
        (
            format!("Cluster[{P8}, Equation[Type, \"average\"], Equation[Count, 3]]"),
            vec![0, 0, 0, 1, 1, 1, 2, 2],
        ),
        // |p − q| as the distance: heights 1, 1, 1, 8, 9, the gap of 7.
        (
            "Cluster[List[1, 2, 10, 11, 20, 21], \
             Equation[distance, Function[List[p, q], Abs[Plus[p, Times[-1, q]]]]]]"
                .to_owned(),
            vec![0, 0, 1, 1, 2, 2],
        ),
    ];
    for t in ["complete", "single", "average", "ward", "weighted"] {
        rows.push((
            format!("Cluster[{D6}, Equation[Type, \"{t}\"]]"),
            alike.to_vec(),
        ));
        for count in ["3", "Range[2, 4]"] {
            rows.push((
                format!("Cluster[{E6}, Equation[Type, \"{t}\"], Equation[Count, {count}]]"),
                vec![0, 0, 1, 1, 2, 2],
            ));
        }
    }
    for (input, expected) in rows {
        let found = eval(&ctx, &input);
        assert_eq!(labels(&input, &found), expected, "{input}");
    }
    assert_eq!(
        eval(&ctx, &format!("Cluster[{D6}, Equation[output, Count]]")),
        Expr::from(2)
    );
}

/// `kmeans`' seed-independent partitions, up to relabelling.
#[test]
fn kmeans_partitions() {
    let ctx = Context::new();
    let four = "List[List[0, 0], List[0, 1], List[5, 5], List[5, 6]]";
    for seed in [0, 1, 42] {
        eval(&ctx, &format!("SeedRandom[{seed}]"));
        for (input, expected) in [
            (format!("KMeans[{D6}]"), vec![0, 1, 0, 1, 0, 1]),
            (format!("KMeans[{D6}, 2]"), vec![0, 1, 0, 1, 0, 1]),
            (
                format!("KMeans[{D6}, Equation[output, list]]"),
                vec![0, 1, 0, 1, 0, 1],
            ),
            (format!("KMeans[{four}, 2]"), vec![0, 0, 1, 1]),
            (format!("KMeans[Matrix[{four}], 2]"), vec![0, 0, 1, 1]),
            (
                "KMeans[List[List[0, 0], List[0, 1], List[5, 5], List[5, 6], List[1/2, 1/3]], 2]"
                    .to_owned(),
                vec![0, 0, 1, 1, 0],
            ),
            (
                "KMeans[List[List[1], List[2], List[10], List[11]], 2]".to_owned(),
                vec![0, 0, 1, 1],
            ),
            ("KMeans[List[1, 2, 10, 11], 2]".to_owned(), vec![0, 0, 1, 1]),
        ] {
            let found = eval(&ctx, &input);
            assert_eq!(canonical(&labels(&input, &found)), expected, "{input}");
        }
        for input in [
            format!("KMeans[{E6}, 3, Equation[output, Count]]"),
            format!("KMeans[{P8}, Equation[output, Count]]"),
        ] {
            assert_eq!(eval(&ctx, &input), Expr::from(3), "{input}");
        }
    }
}

/// The seed-dependent rows: each answer among the prototype's partitions
/// over 200 seeds.
#[test]
fn kmeans_local_optima() {
    let ctx = Context::new();
    let rows: [(String, &[&[i64]]); 7] = [
        (
            format!("KMeans[{D6}, 3]"),
            &[&[0, 1, 0, 2, 0, 1], &[0, 1, 2, 1, 0, 1]],
        ),
        (
            "KMeans[List[1, 2, 10, 11, 20, 21], 3]".to_owned(),
            &[&[0, 0, 1, 1, 2, 2], &[0, 0, 0, 0, 1, 2]],
        ),
        (
            format!("KMeans[{E6}]"),
            &[
                &[0, 0, 1, 1, 2, 2],
                &[0, 0, 0, 0, 1, 1],
                &[0, 0, 1, 1, 1, 1],
            ],
        ),
        (
            format!("KMeans[{E6}, Equation[Count, Range[2, 4]]]"),
            &[
                &[0, 0, 1, 1, 2, 2],
                &[0, 0, 0, 0, 1, 1],
                &[0, 0, 1, 1, 1, 1],
            ],
        ),
        (
            format!("KMeans[{E6}, 3]"),
            &[
                &[0, 0, 1, 1, 2, 2],
                &[0, 1, 2, 2, 2, 2],
                &[0, 0, 0, 0, 1, 2],
            ],
        ),
        (
            format!("KMeans[{P8}, 3]"),
            &[&[0, 0, 0, 1, 1, 1, 2, 2], &[0, 0, 0, 0, 0, 0, 1, 2]],
        ),
        // min_k = max_k = 2 overrides the 5 given.
        (
            "KMeans[List[List[1, 2], List[3, 4], List[5, 6]], 5]".to_owned(),
            &[&[0, 1, 1], &[0, 0, 1]],
        ),
    ];
    for seed in 0..8 {
        eval(&ctx, &format!("SeedRandom[{seed}]"));
        for (input, partitions) in &rows {
            let found = canonical(&labels(input, &eval(&ctx, input)));
            assert!(
                partitions.contains(&found.as_slice()),
                "{input}, seed {seed}: {found:?}"
            );
        }
    }
}

/// `SeedRandom` repeats a run.
#[test]
fn seeded_runs_repeat() {
    let ctx = Context::new();
    for seed in [1, 7, 42] {
        for call in [format!("KMeans[{E6}, 3]"), format!("KMeans[{P8}]")] {
            let twice = eval(
                &ctx,
                &format!("List[SeedRandom[{seed}], {call}, SeedRandom[{seed}], {call}]"),
            );
            assert_eq!(twice.args()[1], twice.args()[3], "{call}, seed {seed}");
        }
    }
}

/// `k ≤ 1` or `k ≥ m`: giac's `[]`.
#[test]
fn kmeans_sizes() {
    let ctx = Context::new();
    for input in [format!("KMeans[{E6}, 1]"), format!("KMeans[{E6}, 10]")] {
        assert_eq!(eval(&ctx, &input), read("List[]"), "{input}");
    }
}

/// `input` with its arguments evaluated: the call as it stays.
fn stayed(ctx: &Context, input: &str) -> Expr {
    read(input).map_args(|a| ctx.eval(a).unwrap_or_else(|e| panic!("{input}: {e}")))
}

/// The one warning of `input`, which stays.
fn warning(ctx: &Context, input: &str) -> String {
    let e = read(input);
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    assert_eq!(found.unwrap(), stayed(ctx, input), "{input}");
    match &report.warnings[..] {
        [Warning::Native { message, .. }] => message.to_string(),
        other => panic!("{input}: {other:?}"),
    }
}

/// giac's errors stay with its message.
#[test]
fn stays() {
    let ctx = Context::new();
    for (input, message) in [
        ("KMeans[List[]]".to_owned(), "No data found"),
        ("Cluster[List[]]".to_owned(), "No data found"),
        ("KMeans[List[a, 1, 2]]".to_owned(), "Data must be numeric"),
        ("Cluster[List[a, 1, 2]]".to_owned(), "Bad distance function"),
        (
            format!(
                "Cluster[{E6}, Equation[Type, \"complete\"], Equation[distance, Function[List[p, q], x]]]"
            ),
            "Bad distance function",
        ),
        ("KMeans[List[1, 2]]".to_owned(), "Invalid k-range"),
        (
            "KMeans[List[List[1, 2], List[3, 4]], 2]".to_owned(),
            "Invalid k-range",
        ),
        (
            format!("Cluster[{E6}, Equation[CountLess, 5], Equation[CountGreater, 3]]"),
            "Invalid k-range",
        ),
        (
            format!("KMeans[{E6}, Equation[output, foo]]"),
            "Invalid output specification",
        ),
        (
            format!("Cluster[{E6}, Equation[output, foo]]"),
            "Invalid output specification",
        ),
        (
            format!("KMeans[{E6}, Equation[output, FirstPosition]]"),
            "No index specified for output",
        ),
        (
            format!("Cluster[{E6}, Equation[output, FirstPosition]]"),
            "No index specified for output",
        ),
        (
            format!("Cluster[{E6}, Equation[Type, \"centroid\"]]"),
            "Unknown linkage type",
        ),
        (
            format!("Cluster[{E6}, Equation[Type, complete]]"),
            "Expected a string",
        ),
        (
            format!("Cluster[{E6}, Equation[distance, 3]]"),
            "Invalid distance function specification",
        ),
        (format!("Cluster[{E6}, foo]"), "Invalid option"),
        (format!("Cluster[{E6}, Equation[foo, 1]]"), "Invalid option"),
        (format!("KMeans[{E6}, Range[2, 4]]"), "Invalid option"),
        (format!("KMeans[{E6}, 0]"), "Expected a positive integer"),
        (
            format!("Cluster[{E6}, Equation[Count, 0]]"),
            "Expected a positive integer",
        ),
        (
            format!("KMeans[{E6}, Equation[Limit, 0]]"),
            "Expected a positive integer",
        ),
        (
            format!("Cluster[{E6}, Equation[Count, Range[3, 2]]]"),
            "Expected a range of positive integers",
        ),
        (
            format!("Cluster[{E6}, Equation[Count, x]]"),
            "Expected an integer or range",
        ),
    ] {
        assert_eq!(warning(&ctx, &input), message, "{input}");
    }
}

/// Graphics and indices, another shape of data: the call stays without a
/// warning.
#[test]
fn stays_silently() {
    let ctx = Context::new();
    for input in [
        format!("KMeans[{E6}, Equation[output, plot]]"),
        format!("Cluster[{E6}, Equation[output, tree]]"),
        format!("Cluster[{E6}, Equation[display, 1]]"),
        "Cluster[List[List[1, 2], List[3]]]".to_owned(),
        "KMeans[List[List[1, 2], 3, 4]]".to_owned(),
        "KMeans[x]".to_owned(),
        "Cluster[]".to_owned(),
    ] {
        let e = read(&input);
        let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
        assert_eq!(found.unwrap(), stayed(&ctx, &input), "{input}");
        assert!(report.warnings.is_empty(), "{input}: {:?}", report.warnings);
    }
}

/// Past `Limits::terms`: "too many points".
#[test]
fn limits() {
    let mut small = Context::new();
    let limits = small.limits();
    small.set_limits(Limits {
        terms: 20,
        ..limits
    });
    for input in [format!("KMeans[{E6}, 3]"), format!("Cluster[{P8}]")] {
        assert_eq!(warning(&small, &input), "too many points", "{input}");
    }
}

const E12: &str = "List[List[0, 0], List[0, 1], List[1, 0], List[10, 10], List[10, 11], \
                   List[11, 10], List[20, 0], List[20, 1], List[21, 0], List[0, 20], \
                   List[1, 20], List[0, 21]]";
const WORDS: &str = "List[\"hello\", \"hallo\", \"world\", \"word\", \"help\"]";

/// `output=part`: the clusters by label, each the data items in order.
#[test]
fn parts() {
    let ctx = Context::new();
    for (input, expected) in [
        (
            format!("Cluster[{E6}, Equation[output, Part]]"),
            "List[List[List[0, 0], List[0, 1]], List[List[5, 5], List[5, 6]], \
             List[List[10, 0], List[10, 1]]]",
        ),
        (
            "Cluster[List[1, 2, 10, 11, 20, 21], Equation[output, Part]]".to_owned(),
            "List[List[1, 2], List[10, 11], List[20, 21]]",
        ),
        (
            format!("Cluster[{WORDS}, Equation[Type, \"complete\"], Equation[output, Part]]"),
            "List[List[\"hello\", \"hallo\"], List[\"world\", \"word\"], List[\"help\"]]",
        ),
    ] {
        assert_eq!(eval(&ctx, &input), read(expected), "{input}");
    }
    // kmeans' labels depend on the draw: the clusters as a set.
    let input = "KMeans[List[List[0, 0], List[0, 1], List[5, 5], List[5, 6]], 2, \
                 Equation[output, Part]]";
    let found = eval(&ctx, input);
    let mut groups = found.args().to_vec();
    groups.sort();
    let mut expected = read("List[List[List[0, 0], List[0, 1]], List[List[5, 5], List[5, 6]]]")
        .args()
        .to_vec();
    expected.sort();
    assert_eq!(groups, expected, "{input}: {found:?}");
}

/// Strings: giac's Levenshtein distance, on bytes.
#[test]
fn strings() {
    assert_eq!(super::levenshtein("é".as_bytes(), b"e"), 2);
    assert_eq!(super::levenshtein(b"kitten", b"sitting"), 3);
    let ctx = Context::new();
    for (input, expected) in [
        (format!("Cluster[{WORDS}]"), vec![0, 0, 1, 1, 0]),
        (
            format!("Cluster[{WORDS}, Equation[Type, \"ward\"], Equation[Count, 2]]"),
            vec![0, 0, 1, 1, 0],
        ),
        (
            "Cluster[List[\"kitten\", \"sitting\", \"mitten\", \"abc\", \"abd\"], 3]".to_owned(),
            vec![0, 1, 0, 2, 2],
        ),
        (
            "Cluster[List[\"aaaa\", \"aaab\", \"bbbb\", \"bbba\", \"cccc\"]]".to_owned(),
            vec![0, 0, 1, 1, 2],
        ),
        ("Cluster[List[\"x\"]]".to_owned(), vec![0]),
    ] {
        let found = eval(&ctx, &input);
        assert_eq!(labels(&input, &found), expected, "{input}");
    }
    assert_eq!(
        eval(&ctx, &format!("Cluster[{WORDS}, Equation[output, Count]]")),
        Expr::from(2)
    );
}

/// The values of an `output=index` answer: `[k, value]` rows.
fn index_values(input: &str, found: &Expr) -> Vec<(i64, f64)> {
    found
        .args()
        .iter()
        .map(|row| match row.args() {
            [k, v] => match (k.as_i64(), v.view()) {
                (Some(k), crate::expr::View::Ball(b)) => (k, b.mid().to_f64()),
                _ => panic!("{input}: {found:?}"),
            },
            _ => panic!("{input}: {found:?}"),
        })
        .collect()
}

const SIL: &str = "Equation[FirstPosition, \"silhouette\"]";

/// `cluster`'s silhouette labels (giac's rows, and the prototype's where
/// giac aborts: ours starts at `N − 1` clusters).
#[test]
fn silhouette_labels() {
    let ctx = Context::new();
    let sil = SIL;
    for (input, expected) in [
        (
            format!("Cluster[{E12}, {sil}, Equation[Count, Range[2, 6]]]"),
            vec![0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 4, 5],
        ),
        (
            format!("Cluster[{E12}, {sil}, Equation[Count, Range[2, 3]]]"),
            vec![0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3],
        ),
        (
            format!("Cluster[{E12}, Equation[FirstPosition, List[\"dunn\", \"mcclain-rao\"]]]"),
            vec![0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3],
        ),
        (
            format!("Cluster[{P8}, {sil}, Equation[CountGreater, 7]]"),
            vec![0, 0, 1, 2, 2, 2, 3, 3],
        ),
        (
            format!(
                "Cluster[{P8}, {sil}, Equation[CountGreater, 7], Equation[Type, \"complete\"]]"
            ),
            vec![0, 0, 0, 0, 0, 0, 1, 1],
        ),
        (
            "Cluster[List[0, 0, 0, 5], Equation[FirstPosition, True], Equation[CountGreater, 3]]"
                .to_owned(),
            vec![0, 0, 0, 1],
        ),
        // Departures: giac aborts (max_k 16 ≥ N).
        (
            format!("Cluster[{E12}, {sil}]"),
            vec![0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 4, 5],
        ),
        (
            format!("Cluster[{E12}, {sil}, Equation[Type, \"ward\"]]"),
            vec![0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0],
        ),
    ] {
        let found = eval(&ctx, &input);
        assert_eq!(labels(&input, &found), expected, "{input}");
    }
    assert_eq!(
        eval(
            &ctx,
            &format!("Cluster[{P8}, {sil}, Equation[CountGreater, 7], Equation[output, Count]]")
        ),
        Expr::from(4)
    );
}

/// `cluster`'s `output=index`: giac's silhouette values.
#[test]
fn silhouette_values() {
    let ctx = Context::new();
    let sil = SIL;
    let e12 = [
        (11, 0.25),
        (10, 0.992_808_471_454_880_3),
        (9, 0.993_180_830_751_805_6),
        (8, 0.993_296_570_002_618_4),
        (7, 0.993_296_570_002_618_4),
        (6, 0.993_313_996_053_753_3),
        (5, 0.993_313_996_053_753_3),
        (4, 0.993_313_996_053_753_3),
    ];
    for (input, expected) in [
        (
            format!(
                "Cluster[{E12}, {sil}, Equation[CountGreater, 11], Equation[output, FirstPosition]]"
            ),
            e12.to_vec(),
        ),
        (
            format!("Cluster[{E12}, {sil}, Equation[output, FirstPosition]]"),
            e12.to_vec(),
        ),
        (
            format!(
                "Cluster[{P8}, {sil}, Equation[CountGreater, 7], Equation[output, FirstPosition]]"
            ),
            vec![
                (7, 0.775),
                (6, 0.775),
                (5, 0.897_916_666_666_666_7),
                (4, 0.956_896_551_724_137_9),
                (3, 0.956_896_551_724_137_9),
            ],
        ),
        (
            format!(
                "Cluster[{WORDS}, {sil}, Equation[CountGreater, 4], Equation[output, FirstPosition]]"
            ),
            vec![
                (4, 0.583_333_333_333_333_3),
                (3, 0.75),
                (2, 0.767_857_142_857_142_8),
            ],
        ),
        (
            "Cluster[List[0, 0, 0, 5], Equation[FirstPosition, True], Equation[CountGreater, 3], \
             Equation[output, FirstPosition]]"
                .to_owned(),
            vec![(3, f64::MIN_POSITIVE), (2, 1.0)],
        ),
    ] {
        let found = index_values(&input, &eval(&ctx, &input));
        assert_eq!(found.len(), expected.len(), "{input}: {found:?}");
        for ((k, v), (ek, ev)) in found.iter().zip(&expected) {
            assert!(k == ek && (v - ev).abs() <= 1e-15, "{input}: {found:?}");
        }
    }
}

/// `index=` and `output=index` errors stay with giac's message.
#[test]
fn index_errors() {
    let ctx = Context::new();
    let sil = SIL;
    for (input, message) in [
        (
            format!("Cluster[{E12}, Equation[FirstPosition, \"foo\"]]"),
            "Unknown index name",
        ),
        (
            format!("Cluster[{E12}, Equation[FirstPosition, List[\"a\", \"b\"]]]"),
            "Unknown index name",
        ),
        (
            format!("Cluster[{E12}, Equation[FirstPosition, 1]]"),
            "Expected a string or list of strings",
        ),
        (
            format!("Cluster[{E12}, Equation[FirstPosition, List[]]]"),
            "Expected a string or list of strings",
        ),
        (
            format!(
                "Cluster[{E12}, Equation[FirstPosition, False], Equation[output, FirstPosition]]"
            ),
            "No index specified for output",
        ),
        (
            format!(
                "Cluster[{E12}, Equation[FirstPosition, \"mcclain-rao\"], Equation[output, FirstPosition]]"
            ),
            "Unable to transpose",
        ),
        (
            format!("Cluster[{E12}, {sil}, 4, Equation[output, FirstPosition]]"),
            "Unable to transpose",
        ),
        (
            format!(
                "Cluster[{E12}, {sil}, Equation[Count, Range[2, 3]], Equation[output, FirstPosition]]"
            ),
            "Unable to transpose",
        ),
    ] {
        assert_eq!(warning(&ctx, &input), message, "{input}");
    }
}

/// The silhouette's steps and strings' lengths count against
/// `Limits::terms`.
#[test]
fn index_and_string_limits() {
    let mut small = Context::new();
    let limits = small.limits();
    small.set_limits(Limits {
        terms: 20,
        ..limits
    });
    for input in [
        format!("Cluster[{E6}, Equation[FirstPosition, \"silhouette\"]]"),
        "Cluster[List[\"abcdefghij\", \"b\", \"c\"]]".to_owned(),
    ] {
        assert_eq!(warning(&small, &input), "too many points", "{input}");
    }
}

/// Two well-separated groups of three: `kmeans`' indices are
/// seed-independent there (the prototype over 200 seeds, giac's probes).
const S6: &str =
    "List[List[0, 0], List[0, 1], List[1, 0], List[10, 10], List[10, 11], List[11, 10]]";

/// `kmeans`' `output=index` rows: a value, `None` for "not defined",
/// `ComplexInfinity` (giac's `infinity`) as an infinity.
fn kmeans_values(input: &str, found: &Expr) -> Vec<(i64, Option<f64>)> {
    found
        .args()
        .iter()
        .map(|row| match row.args() {
            [k, v] => {
                let k = k.as_i64().unwrap_or_else(|| panic!("{input}: {found:?}"));
                let v = match v.view() {
                    crate::expr::View::Ball(b) => Some(b.mid().to_f64()),
                    _ if v.as_str() == Some("not defined") => None,
                    _ if v.as_symbol() == Some(Symbol::COMPLEX_INFINITY) => Some(f64::INFINITY),
                    _ => panic!("{input}: {found:?}"),
                };
                (k, v)
            }
            _ => panic!("{input}: {found:?}"),
        })
        .collect()
}

fn kmeans_index(index: &str, rest: &str) -> String {
    format!("KMeans[{S6}, Equation[FirstPosition, {index}]{rest}]")
}

/// `kmeans`' `output=index` (the prototype's values; giac prints them to
/// 12 digits): the first index's, "not defined", infinities.
#[test]
fn kmeans_index_values() {
    let ctx = Context::new();
    let out = ", Equation[output, FirstPosition]";
    let calinski = [
        (2, Some(450.000_000_000_000_06)),
        (3, Some(246.136_363_636_363_65)),
        (4, Some(201.111_111_111_111_1)),
        (5, Some(151.083_333_333_333_3)),
    ];
    let ray_turi = [
        (2, Some(0.002_222_222_222_222_222)),
        (3, Some(0.244_444_444_444_444_4)),
        (4, Some(0.133_333_333_333_333_3)),
        (5, Some(0.083_333_333_333_333_33)),
    ];
    for (input, expected) in [
        (kmeans_index("\"calinski-harabasz\"", out), &calinski[..]),
        (kmeans_index("True", out), &calinski[..]),
        (kmeans_index("\"ray-turi\"", out), &ray_turi[..]),
        // The first index in giac's order, not the list's.
        (
            kmeans_index("List[\"pbm\", \"ray-turi\"]", out),
            &ray_turi[..],
        ),
        (
            kmeans_index("\"ball-hall\"", out),
            &[
                (2, Some(0.444_444_444_444_444_5)),
                (3, Some(0.231_481_481_481_481_48)),
                (4, Some(0.125)),
                (5, Some(0.05)),
            ][..],
        ),
        (
            kmeans_index("\"banfeld-raftery\"", out),
            &[
                (2, Some(-4.865_581_297_297_972)),
                (3, None),
                (4, None),
                (5, None),
            ][..],
        ),
        (
            kmeans_index("\"all\"", out),
            &[
                (2, Some(-4.865_581_297_297_972)),
                (3, None),
                (4, None),
                (5, None),
            ][..],
        ),
    ] {
        let found = kmeans_values(&input, &eval(&ctx, &input));
        assert_eq!(found.len(), expected.len(), "{input}: {found:?}");
        for ((k, v), (ek, ev)) in found.iter().zip(expected) {
            let near = match (v, ev) {
                (Some(v), Some(e)) => (v - e).abs() <= 1e-12 * e.abs().max(1.0),
                (None, None) => true,
                _ => false,
            };
            assert!(k == ek && near, "{input}: {found:?}");
        }
    }
    // Seed-dependent at k = 4; k = 2, 3 and 5 the same in every run: an
    // undefined determinant ratio, a singular W's infinite trace of W⁻¹B.
    for (index, first, last) in [
        (
            "\"det\"",
            [338.249_999_999_999_94, 4.005_921_539_600_296],
            None,
        ),
        (
            "\"trace-WiB\"",
            [281.25, 418.857_142_857_142_95],
            Some(f64::INFINITY),
        ),
    ] {
        let input = kmeans_index(index, out);
        let found = kmeans_values(&input, &eval(&ctx, &input));
        assert_eq!(found.len(), 4, "{input}: {found:?}");
        for ((k, v), e) in found.iter().zip(first) {
            assert!(
                v.is_some_and(|v| (v - e).abs() <= 1e-12 * e) && (*k == 2 || *k == 3),
                "{input}: {found:?}"
            );
        }
        assert_eq!(found[3], (5, last), "{input}");
    }
}

/// `kmeans`' `index=`: the suggested number of clusters and its labels; the
/// last of `k` and `index=` wins; no suggestion answers no labels.
#[test]
fn kmeans_index_counts() {
    let ctx = Context::new();
    let count = ", Equation[output, Count]";
    for (input, k) in [
        (kmeans_index("\"calinski-harabasz\"", count), 2),
        (kmeans_index("\"pbm\"", count), 5),
        (kmeans_index("\"ball-hall\"", count), 3),
        (kmeans_index("\"trace-W\"", count), 4),
        (kmeans_index("\"all\"", count), 2),
        (
            format!("KMeans[{S6}, 2, Equation[FirstPosition, \"pbm\"]{count}]"),
            5,
        ),
        (kmeans_index("\"pbm\"", &format!(", 3{count}")), 3),
        (kmeans_index("False", count), 2),
        // An elbow needs three k: none suggested.
        (
            kmeans_index(
                "\"ball-hall\"",
                &format!(", Equation[Count, Range[2, 3]]{count}"),
            ),
            0,
        ),
    ] {
        assert_eq!(eval(&ctx, &input), Expr::from(k), "{input}");
    }
    let input = kmeans_index("\"calinski-harabasz\"", "");
    assert_eq!(
        canonical(&labels(&input, &eval(&ctx, &input))),
        [0, 0, 0, 1, 1, 1],
        "{input}"
    );
    let empty = read("List[]");
    for rest in [
        ", Equation[Count, Range[2, 3]]",
        ", Equation[Count, Range[2, 3]], Equation[output, Part]",
        // k fixed by a range of one: no index values.
        ", Equation[Count, Range[2, 2]], Equation[output, FirstPosition]",
        // k = 1: AS 136's fault 3.
        ", Equation[Count, Range[1, 3]]",
    ] {
        let input = kmeans_index("\"ball-hall\"", rest);
        assert_eq!(eval(&ctx, &input), empty, "{input}");
    }
}

/// `kmeans`' `index=` errors stay with giac's message.
#[test]
fn kmeans_index_errors() {
    let ctx = Context::new();
    let out = ", Equation[output, FirstPosition]";
    for (input, message) in [
        (kmeans_index("\"foo\"", ""), "Unknown index name"),
        (kmeans_index("\"Ball-hall\"", ""), "Unknown index name"),
        (kmeans_index("List[\"foo\"]", ""), "Unknown index name"),
        (
            kmeans_index("5", ""),
            "Expected a string or list of strings",
        ),
        (
            kmeans_index("List[]", ""),
            "Expected a string or list of strings",
        ),
        (kmeans_index("False", out), "No index specified for output"),
        (
            kmeans_index("\"pbm\"", &format!(", Equation[Count, 2]{out}")),
            "No index specified for output",
        ),
    ] {
        assert_eq!(warning(&ctx, &input), message, "{input}");
    }
}
