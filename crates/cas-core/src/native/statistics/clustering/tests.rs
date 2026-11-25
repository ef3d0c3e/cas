//! Tests of clustering (CAS-630): `cluster`'s labels and `kmeans`'
//! partitions (up to relabelling) against the prototype's answers
//! (.run/1.7/T6/proto.py, R1.7-012), the seed-dependent ones among the
//! partitions it found over 200 seeds, a run repeated under `SeedRandom`,
//! each warning, and the limits.

use crate::context::Limits;
use crate::trace::Warning;
use crate::{Context, EvalOptions, Expr};

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
        format!("KMeans[{E6}, Equation[output, Part]]"),
        format!("KMeans[{E6}, Equation[FirstPosition, True]]"),
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
