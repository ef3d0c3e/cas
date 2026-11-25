//! Xcas's clustering against giac (CAS-630): `cluster` and `kmeans`
//! (giac's answers in .run/1.7/giac/probe4.out and probe8.out; the
//! prototype .run/1.7/T6/proto.py, R1.7-012). `cluster` is deterministic:
//! labels compared as they are. `kmeans` starts from random centres and
//! giac's `_rand` stream is not ours: partitions compared up to relabelling
//! (labels renumbered by first appearance, C0087), the seed-independent rows
//! exactly, the seed-dependent ones as departures, each answer among the
//! partitions the prototype found over 200 seeds (giac's among them).
//! Invalid calls stay with giac's message as the warning (C0022); graphics
//! (`output=plot`, `output=tree`) stay without one. Seeded rows (C0014):
//! ours repeats under `srand`.

use cas_core::{Context, EvalOptions, Expr, Warning};
use cas_giac::{Giac, ours};
use cas_lang::{Hook, Names, decimal_balls, lower, parse};

const D6: &str = "[[1.1,1.5],[3.2,2.5],[1.5,0.9],[3.5,2.7],[0.8,1.1],[3.1,2.2]]";
const E6: &str = "[[0,0],[0,1],[5,5],[5,6],[10,0],[10,1]]";

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// The labels of a list of integers.
fn labels(e: &Expr) -> Option<Vec<i64>> {
    e.args().iter().map(Expr::as_i64).collect()
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

/// `cluster`'s labels (probe4, probe8), giac's and ours, as they are.
#[test]
fn cluster_rows() {
    let rows = [
        (format!("cluster({D6})"), "[0,1,0,1,0,1]"),
        (format!("cluster({D6},2)"), "[0,1,0,1,0,1]"),
        (format!("cluster({D6},count=2)"), "[0,1,0,1,0,1]"),
        (format!("cluster({D6},3)"), "[0,1,2,1,0,1]"),
        (format!("cluster({D6},count=3)"), "[0,1,2,1,0,1]"),
        (format!("cluster({D6},count=2..4)"), "[0,1,0,1,0,1]"),
        (format!("cluster({D6},type=\"complete\")"), "[0,1,0,1,0,1]"),
        (format!("cluster({D6},type=\"single\")"), "[0,1,0,1,0,1]"),
        (format!("cluster({D6},type=\"average\")"), "[0,1,0,1,0,1]"),
        (format!("cluster({D6},type=\"ward\")"), "[0,1,0,1,0,1]"),
        (format!("cluster({D6},type=\"weighted\")"), "[0,1,0,1,0,1]"),
        (format!("cluster({D6},output=list)"), "[0,1,0,1,0,1]"),
        (format!("cluster({D6},output=count)"), "2"),
        (format!("cluster({E6})"), "[0,0,1,1,2,2]"),
        (format!("cluster({E6},2)"), "[0,0,0,0,1,1]"),
        (format!("cluster({E6},type=\"single\",3)"), "[0,0,1,1,2,2]"),
        ("cluster([1,2,10,11,20,21])".to_owned(), "[0,0,1,1,2,2]"),
        ("cluster([1,2,10,11,20,21],3)".to_owned(), "[0,0,1,1,2,2]"),
        ("cluster([1/2,1,3,7/2])".to_owned(), "[0,0,1,1]"),
        ("cluster([[1,2]])".to_owned(), "[0]"),
        ("cluster([[1,2],[3,4]])".to_owned(), "[0,1]"),
    ];
    let giac = giac();
    let mut wrong = Vec::new();
    for (input, answer) in &rows {
        if let Some(giac) = giac
            && let Err(e) = giac.answers(input, answer)
        {
            wrong.push(format!("{input}: giac {e:?}"));
        }
        let (mine, expected) = (ours(input), ours(answer));
        if !matches!((&mine, &expected), (Ok(m), Ok(x)) if m == x) {
            wrong.push(format!("{input}: ours {mine:?}, wanted {answer}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// The last value of giac's answer to `input`, which may begin with
/// `srand(s);` (C0014).
fn giac_last(giac: &Giac, input: &str) -> Option<Expr> {
    let answer = giac.eval(input).ok()?;
    let mut depth = 0i32;
    let mut start = 0;
    for (i, c) in answer.char_indices() {
        match c {
            '(' | '[' => depth += 1,
            ')' | ']' => depth -= 1,
            ',' if depth == 0 => start = i + 1,
            _ => {}
        }
    }
    ours(&answer[start..]).ok()
}

/// Ours of `input` after `srand(seed)` in the same input, or of `input`.
fn ours_seeded(seed: Option<u32>, input: &str) -> Option<Expr> {
    match seed {
        Some(s) => ours(&format!("[srand({s}),{input}]"))
            .ok()
            .and_then(|e| e.args().get(1).cloned()),
        None => ours(input).ok(),
    }
}

/// `input` after `srand(seed);` for giac.
fn with_seed(seed: Option<u32>, input: &str) -> String {
    seed.map_or_else(|| input.to_owned(), |s| format!("srand({s});{input}"))
}

/// `kmeans`' seed-independent rows (R1.7-012: one partition over 200
/// prototype seeds): giac's partition and ours, up to relabelling.
#[test]
fn kmeans_partitions() {
    let four = "[[0,0],[0,1],[5,5],[5,6]]";
    let rows = [
        (None, format!("kmeans({D6})"), "[1,0,1,0,1,0]"),
        (Some(42), format!("kmeans({D6})"), "[1,0,1,0,1,0]"),
        (Some(1), format!("kmeans({D6})"), "[1,0,1,0,1,0]"),
        (Some(42), format!("kmeans({D6},2)"), "[1,0,1,0,1,0]"),
        (
            Some(42),
            format!("kmeans({D6},output=list)"),
            "[1,0,1,0,1,0]",
        ),
        (
            Some(42),
            "kmeans([[1],[2],[10],[11]],2)".to_owned(),
            "[1,1,0,0]",
        ),
        (None, format!("kmeans(matrix({four}),2)"), "[1,1,0,0]"),
        (None, format!("kmeans({four},2)"), "[1,1,0,0]"),
        (
            None,
            "kmeans([[0,0],[0,1],[5,5],[5,6],[1/2,1/3]],2)".to_owned(),
            "[1,1,0,0,1]",
        ),
    ];
    let giac = giac();
    let mut wrong = Vec::new();
    for (seed, input, answer) in &rows {
        let expected = ours(answer)
            .ok()
            .and_then(|e| labels(&e))
            .map(|l| canonical(&l));
        let partition = |e: Option<Expr>| e.and_then(|e| labels(&e)).map(|l| canonical(&l));
        if let Some(giac) = giac {
            let theirs = partition(giac_last(giac, &with_seed(*seed, input)));
            if theirs != expected {
                wrong.push(format!("{input}: giac {theirs:?}, recorded {answer}"));
            }
        }
        let mine = partition(ours_seeded(*seed, input));
        if mine.is_none() || mine != expected {
            wrong.push(format!(
                "{input}: ours {mine:?}, wanted {answer} up to relabelling"
            ));
        }
    }
    let count = format!("kmeans({E6},3,output=count)");
    if let Some(giac) = giac
        && let Err(e) = giac.answers(&count, "3")
    {
        wrong.push(format!("{count}: giac {e:?}"));
    }
    if ours(&count).ok() != ours("3").ok() {
        wrong.push(format!("{count}: ours {:?}", ours(&count)));
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// A seed, an input, the partitions allowed, whether giac is asked.
type Row = (Option<u32>, String, &'static [&'static [i64]], bool);

/// Departures, `kmeans`' seed-dependent rows (R1.7-012): giac's and ours
/// each among the prototype's partitions over 200 seeds (giac's recorded
/// answer in the comment). `kmeans([1,2,10,11,20,21],3)`: giac crashes
/// (`mcols` of scalar data, optimization.cc:10655), not asked; ours clusters
/// the 1-D points.
#[test]
fn kmeans_departures() {
    const E: [&[i64]; 3] = [
        &[0, 0, 1, 1, 2, 2],
        &[0, 1, 2, 2, 2, 2],
        &[0, 0, 0, 0, 1, 2],
    ];
    const AUTO: [&[i64]; 2] = [&[0, 0, 1, 1, 2, 2], &[0, 0, 0, 0, 1, 1]];
    let rows: [Row; 9] = [
        // giac [1,0,1,2,1,0].
        (
            Some(42),
            format!("kmeans({D6},3)"),
            &[&[0, 1, 0, 2, 0, 1], &[0, 1, 2, 1, 0, 1]],
            true,
        ),
        (
            Some(42),
            "kmeans([1,2,10,11,20,21],3)".to_owned(),
            &[&[0, 0, 1, 1, 2, 2], &[0, 0, 0, 0, 1, 2]],
            false,
        ),
        // giac [1,1,0,0,2,2].
        (None, format!("kmeans({E6})"), &AUTO, true),
        // giac [1,1,0,0,2,2].
        (None, format!("kmeans({E6},count=2..4)"), &AUTO, true),
        // giac [2,2,1,1,0,0].
        (None, format!("kmeans({E6},3)"), &E, true),
        // giac [2,2,0,0,1,1].
        (Some(7), format!("kmeans({E6},3)"), &E, true),
        // giac [2,2,1,1,0,0].
        (None, format!("kmeans({E6},count=3)"), &E, true),
        // giac [2,2,1,1,0,0].
        (None, format!("kmeans({E6},3,output=list)"), &E, true),
        (Some(1), format!("kmeans({E6},3)"), &E, true),
    ];
    let giac = giac();
    let mut wrong = Vec::new();
    let partition = |e: Option<Expr>| e.and_then(|e| labels(&e)).map(|l| canonical(&l));
    for (seed, input, partitions, ask_giac) in &rows {
        let among = |p: &Option<Vec<i64>>| {
            p.as_ref()
                .is_some_and(|p| partitions.contains(&p.as_slice()))
        };
        if *ask_giac && let Some(giac) = giac {
            let theirs = partition(giac_last(giac, &with_seed(*seed, input)));
            if !among(&theirs) {
                wrong.push(format!("{input}: giac {theirs:?}"));
            }
        }
        let mine = partition(ours_seeded(*seed, input));
        if !among(&mine) {
            wrong.push(format!("{input}: ours {mine:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Seeded rows: `srand(s)` repeats ours.
#[test]
fn seeded_rows() {
    let mut wrong = Vec::new();
    for s in [1, 7, 42, 2024] {
        let input = format!("[srand({s}),kmeans({E6},3),srand({s}),kmeans({E6},3)]");
        match ours(&input) {
            Ok(e) if e.args().len() == 4 && e.args()[1] == e.args()[3] => {}
            other => wrong.push(format!("{input}: {other:?}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `text` read as [`ours`] reads it, and its value with the warnings of
/// its evaluation.
fn evaluated(text: &str) -> Option<(Expr, Expr, Vec<Warning>)> {
    let parsed = parse(text);
    if !parsed.errors().is_empty() {
        return None;
    }
    let lowered = lower(&parsed, Names::english());
    let [statement] = &lowered.statements[..] else {
        return None;
    };
    let e = decimal_balls(&statement.expr, cas_core::DEFAULT_WORKING_PRECISION);
    let mut ctx = Context::new();
    ctx.set_language(Some(Hook::english()));
    let (found, report) = ctx.eval_with(&e, &EvalOptions::default());
    Some((e, found.ok()?, report.warnings))
}

/// Whether ours of `input` stays, its head the input's, with a warning of
/// `expected` (or none for `None`).
fn stays_with(input: &str, expected: Option<&str>) -> bool {
    evaluated(input).is_some_and(|(e, found, warnings)| {
        found.head() == e.head()
            && match expected {
                Some(expected) => warnings
                    .iter()
                    .any(|w| matches!(w, Warning::Native { message, .. } if message == expected)),
                None => warnings.is_empty(),
            }
    })
}

/// Whitespace collapsed (C0077).
fn collapsed(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ")
}

/// giac's errors (probe4, probe8) stay with its message, among them the
/// agreeing stays rows `kmeans(D6,2..4)` (a range only as `count=`),
/// `kmeans([[1,2],[3,4]],2)` and `cluster(D6,type=complete)` (C0056).
#[test]
fn invalid_calls_stay() {
    let value = |m: &str| format!("GIAC_ERROR: {m}. Error: Bad Argument Value");
    let kind = |m: &str| format!("GIAC_ERROR: {m}. Error: Bad Argument Type");
    let output = "Invalid output specification";
    let mut rows = vec![
        (
            format!("cluster({D6},type=\"centroid\")"),
            value("Unknown linkage type"),
        ),
        (
            format!("cluster({D6},type=\"median\")"),
            value("Unknown linkage type"),
        ),
        (
            format!("cluster({D6},distance=\"manhattan\")"),
            kind("Invalid distance function specification"),
        ),
        (
            format!("cluster({D6},output=index,2)"),
            value("No index specified for output"),
        ),
        (
            format!("cluster({D6},output=index)"),
            value("No index specified for output"),
        ),
        ("cluster([x,y])".to_owned(), value("Bad distance function")),
        ("cluster([])".to_owned(), value("No data found")),
        (format!("kmeans({E6},2..4)"), value("Invalid option")),
        (format!("kmeans({D6},2..4)"), value("Invalid option")),
        (format!("kmeans({D6},2,index)"), value("Invalid option")),
        (format!("kmeans({D6},2,count)"), value("Invalid option")),
        (
            format!("kmeans({E6},3,output=index,1)"),
            value("No index specified for output"),
        ),
        (
            format!("kmeans({D6},output=index)"),
            value("No index specified for output"),
        ),
        (
            "kmeans([[1,2],[3,4]],2)".to_owned(),
            value("Invalid k-range"),
        ),
        (
            "kmeans([[1,2],[3,4]],3)".to_owned(),
            value("Invalid k-range"),
        ),
        (
            "kmeans([[1,2],[3,4]],1)".to_owned(),
            value("Invalid k-range"),
        ),
    ];
    for t in ["complete", "single", "average", "ward"] {
        rows.push((format!("cluster({D6},type={t})"), kind("Expected a string")));
    }
    for o in [
        "centroids",
        "centers",
        "tree",
        "centroid",
        "means",
        "clusters",
    ] {
        rows.push((format!("kmeans({E6},3,output={o})"), value(output)));
    }
    for o in ["centroids", "all"] {
        rows.push((format!("kmeans({D6},output={o})"), value(output)));
    }
    let giac = giac();
    let mut wrong = Vec::new();
    for (input, error) in &rows {
        if let Some(giac) = giac {
            let printed = giac.eval(input).unwrap_or_default();
            if collapsed(&printed) != collapsed(error) {
                wrong.push(format!("{input}: giac prints {printed}, recorded {error}"));
            }
        }
        let message = error
            .strip_prefix("GIAC_ERROR: ")
            .and_then(|m| m.split_once(". Error:"))
            .map_or("", |(m, _)| m);
        if !stays_with(input, Some(message)) {
            wrong.push(format!(
                "{input}: ours {:?}, wanted the call with {message:?}",
                evaluated(input)
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Graphics, not ported: giac draws, ours stays without a warning.
#[test]
fn graphics_stay() {
    let giac = giac();
    let mut wrong = Vec::new();
    for (input, drawing) in [
        (
            format!("kmeans({E6},3,output=plot)"),
            "[title=\"Number of clusters: 3\"",
        ),
        (
            format!("cluster({D6},output=plot)"),
            "[title=\"Number of clusters: 2\"",
        ),
        (format!("cluster({D6},output=tree)"), "Graphic_object"),
    ] {
        if let Some(giac) = giac {
            let printed = giac.eval(&input).unwrap_or_default();
            if !printed.starts_with(drawing) {
                wrong.push(format!("{input}: giac prints {printed}"));
            }
        }
        if !stays_with(&input, None) {
            wrong.push(format!("{input}: ours {:?}", evaluated(&input)));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
