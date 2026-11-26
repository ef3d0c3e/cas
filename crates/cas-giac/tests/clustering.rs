//! Xcas's clustering against giac (CAS-630): `cluster` and `kmeans`
//! (giac's answers in crates/cas-core/proto/giac/probe4.cas and probe8.cas; the
//! prototype crates/cas-core/proto/clustering/proto.py, R1.7-012). `cluster` is deterministic:
//! labels compared as they are. `kmeans` starts from random centres and
//! giac's `_rand` stream is not ours: partitions compared up to relabelling
//! (labels renumbered by first appearance, C0087), the seed-independent rows
//! exactly, the seed-dependent ones as departures, each answer among the
//! partitions the prototype found over 200 seeds (giac's among them).
//! Invalid calls stay with giac's message as the warning (C0022); graphics
//! (`output=plot`, `output=tree`) stay without one. Seeded rows (C0014):
//! ours repeats under `srand`. 1.8/T5 (CAS-646, giac's answers in
//! crates/cas-core/proto/giac/probe_clustering2.out, the prototype's `t5`):
//! `output=part`, lists of strings, `cluster`'s `index=`; the silhouette's
//! values within 1e-9 (C0100). giac aborts on the silhouette with `max_k ≥
//! N` (optimization.cc:11325): those rows are departures, giac not asked
//! (C0021). T5.c3 (CAS-646, giac's answers in
//! crates/cas-core/proto/giac/probe_clustering3.out, the prototype's `t5c3`):
//! `kmeans`' `index=`, its seed-independent values and counts as rows,
//! the seed-dependent ones as departures (each `[k, value]` among the
//! prototype's over 200 seeds, giac's among them).

use cas_core::{Context, EvalOptions, Expr, Symbol, View, Warning};
use cas_giac::{Giac, ours};
use cas_lang::{Hook, Names, decimal_balls, lower, parse};

const D6: &str = "[[1.1,1.5],[3.2,2.5],[1.5,0.9],[3.5,2.7],[0.8,1.1],[3.1,2.2]]";
const E6: &str = "[[0,0],[0,1],[5,5],[5,6],[10,0],[10,1]]";
const E12: &str =
    "[[0,0],[0,1],[1,0],[10,10],[10,11],[11,10],[20,0],[20,1],[21,0],[0,20],[1,20],[0,21]]";
const P8: &str = "[[0,0],[1,0],[0,2],[4,4],[5,5],[4,6],[9,1],[10,0]]";
const WORDS: &str = r#"["hello","hallo","world","word","help"]"#;
/// Two well-separated groups of three: `kmeans`' indices seed-independent.
const S6: &str = "[[0,0],[0,1],[1,0],[10,10],[10,11],[11,10]]";

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
    // E6 is symmetric under x -> 10 - x: each 2-cluster answer has its
    // mirror ([0,0,1,1,1,1], seen on CI 45606968339 at seed 3).
    const AUTO: [&[i64]; 3] = [
        &[0, 0, 1, 1, 2, 2],
        &[0, 0, 0, 0, 1, 1],
        &[0, 0, 1, 1, 1, 1],
    ];
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

/// giac's answer and ours, the same expression: `output=part`, strings'
/// Levenshtein distance, and the labels of `index=` (none of these
/// silhouette rows has `max_k ≥ N`).
fn answer_rows(rows: &[(String, &str)]) {
    let giac = giac();
    let mut wrong = Vec::new();
    for (input, answer) in rows {
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

/// `cluster`'s `output=part`, and `kmeans`' where one partition is reached.
#[test]
fn parts() {
    let groups4 = "[[[0,0],[0,1],[1,0]],[[10,10],[10,11],[11,10]],[[20,0],[20,1],[21,0]],\
                   [[0,20],[1,20],[0,21]]]";
    answer_rows(&[
        (
            format!("cluster({D6},output=part)"),
            "[[[1.1,1.5],[1.5,0.9],[0.8,1.1]],[[3.2,2.5],[3.5,2.7],[3.1,2.2]]]",
        ),
        (
            format!("cluster({D6},count=3,output=part)"),
            "[[[1.1,1.5],[0.8,1.1]],[[3.2,2.5],[3.5,2.7],[3.1,2.2]],[[1.5,0.9]]]",
        ),
        (
            format!("cluster({D6},2,output=part,count=3)"),
            "[[[1.1,1.5],[0.8,1.1]],[[3.2,2.5],[3.5,2.7],[3.1,2.2]],[[1.5,0.9]]]",
        ),
        (format!("cluster({E12},count=4,output=part)"), groups4),
        (format!("cluster({E12},4,output=part)"), groups4),
        (
            "cluster([1,2,10,11,20,21],output=part)".to_owned(),
            "[[1,2],[10,11],[20,21]]",
        ),
        (
            format!("cluster({WORDS},output=part)"),
            r#"[["hello","hallo","help"],["world","word"]]"#,
        ),
        (
            format!("cluster({WORDS},type=\"complete\",output=part)"),
            r#"[["hello","hallo"],["world","word"],["help"]]"#,
        ),
        (
            format!("cluster({WORDS},type=\"average\",count=3,output=part)"),
            r#"[["hello","hallo"],["world","word"],["help"]]"#,
        ),
    ]);
    // kmeans: one partition over 200 prototype seeds (D6, 2), or giac's
    // among the prototype's two (D6, 3; ours a departure, C0044).
    let giac = giac();
    let mut wrong = Vec::new();
    let data = ours(D6).expect("D6");
    let partition = |e: Option<Expr>| {
        let e = e?;
        let labels: Option<Vec<i64>> = data
            .args()
            .iter()
            .map(|x| {
                e.args()
                    .iter()
                    .position(|g| g.args().contains(x))
                    .and_then(|i| i64::try_from(i).ok())
            })
            .collect();
        labels.map(|l| canonical(&l))
    };
    let two: &[&[i64]] = &[&[0, 1, 0, 1, 0, 1]];
    let three: &[&[i64]] = &[&[0, 1, 2, 1, 0, 1], &[0, 1, 0, 2, 0, 1]];
    for (seed, input, allowed) in [
        (None, format!("kmeans({D6},2,output=part)"), two),
        (Some(42), format!("kmeans({D6},2,output=part)"), two),
        (None, format!("kmeans({D6},3,output=part)"), three),
    ] {
        let among =
            |p: &Option<Vec<i64>>| p.as_ref().is_some_and(|p| allowed.contains(&p.as_slice()));
        if let Some(giac) = giac {
            let theirs = partition(giac_last(giac, &with_seed(seed, &input)));
            if !among(&theirs) {
                wrong.push(format!("{input}: giac {theirs:?}"));
            }
        }
        let mine = partition(ours_seeded(seed, &input));
        if !among(&mine) {
            wrong.push(format!("{input}: ours {mine:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Lists of strings: giac's Levenshtein distance, on bytes.
#[test]
fn strings() {
    let kitten = r#"["kitten","sitting","mitten","abc","abd"]"#;
    let prefixes = r#"["a","ab","abc","abcd","xyz"]"#;
    let letters = r#"["aaaa","aaab","bbbb","bbba","cccc"]"#;
    answer_rows(&[
        (format!("cluster({WORDS})"), "[0,0,1,1,0]"),
        (format!("cluster({WORDS},count=2)"), "[0,0,1,1,0]"),
        (format!("cluster({WORDS},output=count)"), "2"),
        (
            format!("cluster({WORDS},type=\"ward\",count=2)"),
            "[0,0,1,1,0]",
        ),
        (
            r#"cluster(["abc","abd","xyz","xyw"],count=2,output=list)"#.to_owned(),
            "[0,0,1,1]",
        ),
        (format!("cluster({kitten},count=2)"), "[0,0,0,1,1]"),
        (format!("cluster({kitten},count=3)"), "[0,1,0,2,2]"),
        (format!("cluster({prefixes},count=2)"), "[0,0,0,0,1]"),
        (
            format!("cluster({prefixes},count=2,type=\"complete\")"),
            "[0,0,0,0,1]",
        ),
        (format!("cluster({letters})"), "[0,0,1,1,2]"),
        (format!("cluster({letters},output=count)"), "3"),
        (r#"cluster(["x"])"#.to_owned(), "[0]"),
        (r#"cluster(["aé","ae","zz","zy"])"#.to_owned(), "[0,1,2,3]"),
    ]);
}

/// `index=`'s labels and counts where giac does not abort.
#[test]
fn index_labels() {
    let sil = "index=\"silhouette\"";
    let split = "[0,0,0,1,1,1,2,2,2,3,4,5]";
    let four = "[0,0,0,1,1,1,2,2,2,3,3,3]";
    answer_rows(&[
        (format!("cluster({E12},{sil},count=2..6)"), split),
        (format!("cluster({E12},{sil},count=2..6,output=count)"), "6"),
        (
            format!("cluster({E12},{sil},count=2..6,output=part)"),
            "[[[0,0],[0,1],[1,0]],[[10,10],[10,11],[11,10]],[[20,0],[20,1],[21,0]],\
             [[0,20]],[[1,20]],[[0,21]]]",
        ),
        (format!("cluster({E12},{sil},count_sup=11)"), split),
        (format!("cluster({E12},{sil},count=2..3)"), four),
        (format!("cluster({E12},4,{sil})"), four),
        (
            format!("cluster({E12},{sil},count=11..11)"),
            "[0,0,1,2,3,4,5,6,7,8,9,10]",
        ),
        (format!("cluster({E12},index=false)"), four),
        (
            format!("cluster({E12},index=[\"dunn\",\"mcclain-rao\"],count=2..6)"),
            four,
        ),
        (format!("cluster({D6},index=\"dunn\")"), "[0,1,0,1,0,1]"),
        (format!("cluster({D6},{sil},count_sup=5)"), "[0,1,0,1,0,1]"),
        (
            format!("cluster({P8},{sil},count_sup=7)"),
            "[0,0,1,2,2,2,3,3]",
        ),
        (format!("cluster({P8},{sil},count_sup=7,output=count)"), "4"),
        (
            format!("cluster({P8},{sil},count_sup=7,type=\"complete\")"),
            "[0,0,0,0,0,0,1,1]",
        ),
        (format!("cluster({WORDS},{sil},count_sup=4)"), "[0,0,1,1,0]"),
        (
            "cluster([0,0,0,5],index=\"silhouette\",count_sup=3)".to_owned(),
            "[0,0,0,1]",
        ),
        (r#"cluster(["x"],index="silhouette")"#.to_owned(), "[0]"),
    ]);
}

/// The `[k, value]` rows of an `output=index` answer.
fn index_rows(e: &Expr) -> Option<Vec<(i64, f64)>> {
    e.args()
        .iter()
        .map(|row| match row.args() {
            [k, v] => match v.view() {
                View::Ball(b) => Some((k.as_i64()?, b.mid().to_f64())),
                _ => Some((k.as_i64()?, f64::from(i32::try_from(v.as_i64()?).ok()?))),
            },
            _ => None,
        })
        .collect()
}

/// Rows equal in `k` and within 1e-9 in value (giac prints 12 digits).
fn near(a: &[(i64, f64)], b: &[(i64, f64)]) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b)
            .all(|((k, v), (l, w))| k == l && (v - w).abs() <= 1e-9 * v.abs().max(1.0))
}

const E12_VALUES: [(i64, f64); 8] = [
    (11, 0.25),
    (10, 0.992_808_471_454_880_3),
    (9, 0.993_180_830_751_805_6),
    (8, 0.993_296_570_002_618_4),
    (7, 0.993_296_570_002_618_4),
    (6, 0.993_313_996_053_753_3),
    (5, 0.993_313_996_053_753_3),
    (4, 0.993_313_996_053_753_3),
];

/// `output=index`: giac's silhouette values and ours, within 1e-9.
#[test]
fn silhouette_values() {
    let sil = "index=\"silhouette\"";
    let giac = giac();
    let mut wrong = Vec::new();
    let mut rows = vec![
        format!("cluster({E12},{sil},count=2..6,output=index)"),
        format!("cluster({E12},{sil},count_sup=11,output=index)"),
        format!("cluster({E12},{sil},count_sup=20,output=index)"),
        format!("cluster({E12},{sil},count_inf=3,count_sup=11,output=index)"),
        format!("cluster({E12},index=\"all\",count=2..6,output=index)"),
        format!("cluster({E12},index=true,count=2..6,output=index)"),
        format!("cluster({D6},{sil},count_sup=5,output=index)"),
        format!("cluster({D6},{sil},count_sup=5,type=\"ward\",output=index)"),
        format!("cluster({P8},{sil},count_sup=7,output=index)"),
        format!("cluster({P8},{sil},count_sup=7,type=\"complete\",output=index)"),
        format!("cluster({WORDS},{sil},count_sup=4,output=index)"),
        "cluster([0,0,0,5],index=\"silhouette\",count_sup=3,output=index)".to_owned(),
    ];
    for t in ["complete", "ward", "average", "weighted"] {
        rows.push(format!(
            "cluster({E12},{sil},count=2..6,type=\"{t}\",output=index)"
        ));
    }
    for input in &rows {
        let mine = ours(input).ok().as_ref().and_then(index_rows);
        if mine.as_ref().is_none_or(Vec::is_empty) {
            wrong.push(format!("{input}: ours {mine:?}"));
        }
        if let Some(giac) = giac {
            let theirs = giac.value(input).ok().as_ref().and_then(index_rows);
            if !matches!((&mine, &theirs), (Some(m), Some(t)) if near(m, t)) {
                wrong.push(format!("{input}: ours {mine:?}, giac {theirs:?}"));
            }
        }
    }
    let first = ours(&rows[1]).ok().as_ref().and_then(index_rows);
    if !first.as_ref().is_some_and(|v| near(v, &E12_VALUES)) {
        wrong.push(format!("{}: ours {first:?}", rows[1]));
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departures: the silhouette with `max_k ≥ N` (the default 16) aborts
/// libgiac (optimization.cc:11325, `probe_clustering2.out` `<stderr>`), giac
/// not asked; ours starts at `N − 1` clusters (the prototype's answers).
#[test]
fn silhouette_departures() {
    let sil = "index=\"silhouette\"";
    let mut wrong = Vec::new();
    for (input, answer) in [
        (format!("cluster({E12},{sil})"), "[0,0,0,1,1,1,2,2,2,3,4,5]"),
        (
            format!("cluster({E12},{sil},type=\"ward\")"),
            "[0,0,0,0,0,0,1,1,1,0,0,0]",
        ),
        (format!("cluster({D6},{sil})"), "[0,1,0,1,0,1]"),
        (format!("cluster({D6},index=true)"), "[0,1,0,1,0,1]"),
    ] {
        let (mine, expected) = (ours(&input), ours(answer));
        if !matches!((&mine, &expected), (Ok(m), Ok(x)) if m == x) {
            wrong.push(format!("{input}: ours {mine:?}, wanted {answer}"));
        }
    }
    let d6: &[(i64, f64)] = &[
        (5, 0.493_433_395_872_420_5),
        (4, 0.960_414_356_130_622_3),
        (3, 0.960_535_421_629_749_8),
        (2, 0.965_962_580_045_557_6),
    ];
    let words: &[(i64, f64)] = &[
        (4, 0.583_333_333_333_333_3),
        (3, 0.75),
        (2, 0.767_857_142_857_142_8),
    ];
    for (input, values) in [
        (
            format!("cluster({E12},{sil},output=index)"),
            &E12_VALUES[..],
        ),
        (format!("cluster({D6},{sil},output=index)"), d6),
        (
            format!("cluster({D6},index=[\"silhouette\",\"dunn\"],output=index)"),
            d6,
        ),
        (format!("cluster({D6},index=\"all\",output=index)"), d6),
        (format!("cluster({WORDS},{sil},output=index)"), words),
    ] {
        let mine = ours(&input).ok().as_ref().and_then(index_rows);
        if !mine.as_ref().is_some_and(|m| near(m, values)) {
            wrong.push(format!("{input}: ours {mine:?}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `index=`'s and strings' errors stay with giac's message.
#[test]
fn index_errors_stay() {
    let value = |m: &str| format!("GIAC_ERROR: {m}. Error: Bad Argument Value");
    let kind = |m: &str| format!("GIAC_ERROR: {m}. Error: Bad Argument Type");
    let strings = "Expected a string or list of strings";
    let unknown = "Unknown index name";
    let mut rows = vec![
        (
            format!("cluster({D6},index=[\"a\",\"b\",\"c\",\"d\",\"e\",\"f\"])"),
            value(unknown),
        ),
        (
            format!("cluster({D6},index=[\"a\",\"b\",\"c\",\"d\",\"e\",\"f\"],output=part)"),
            value(unknown),
        ),
        (format!("cluster({D6},index=\"abc\")"), value(unknown)),
        (format!("cluster({D6},index=[\"a\",\"b\"])"), value(unknown)),
        (format!("cluster({E12},index=\"foo\")"), value(unknown)),
        (
            format!("cluster({E12},index=\"Silhouette\")"),
            value(unknown),
        ),
        (
            format!("cluster({D6},index=[1,2,3,4,5,6],output=part)"),
            kind(strings),
        ),
        (format!("cluster({E12},index=1)"), kind(strings)),
        (
            format!("cluster({E12},index=[\"silhouette\",1])"),
            kind(strings),
        ),
        (format!("cluster({E12},index=[])"), kind(strings)),
        (
            format!("cluster({E12},index=[],output=index)"),
            kind(strings),
        ),
        (
            format!("cluster({E12},index=false,output=index)"),
            value("No index specified for output"),
        ),
        (
            r#"kmeans(["hello","hallo","world","word"],2)"#.to_owned(),
            value("Data must be numeric"),
        ),
        (
            r#"kmeans(["a","b","c"],output=part)"#.to_owned(),
            value("Data must be numeric"),
        ),
        (
            r#"cluster([["a","b"],["a","c"]])"#.to_owned(),
            value("Bad distance function"),
        ),
    ];
    for name in ["silhouette", "mcclain", "dunn", "foo"] {
        rows.push((format!("cluster({D6},index={name})"), kind(strings)));
    }
    stay_rows(&rows);
}

/// `output=index` with no silhouette values: giac's `mtran` fails.
#[test]
fn transpose_errors_stay() {
    let transpose = "GIAC_ERROR: Unable to transpose Error: Bad Argument Value";
    let sil = "index=\"silhouette\"";
    let rows = vec![
        (
            format!("cluster({D6},index=\"mcclain-rao\",output=index)"),
            transpose.to_owned(),
        ),
        (
            format!("cluster({E12},index=\"mcclain-rao\",count=2..3,output=index)"),
            transpose.to_owned(),
        ),
        (
            format!("cluster({E12},index=[\"dunn\",\"mcclain-rao\"],count=2..6,output=index)"),
            transpose.to_owned(),
        ),
        (
            format!("cluster({E12},{sil},count=2..3,output=index)"),
            transpose.to_owned(),
        ),
        (
            format!("cluster({E12},4,{sil},output=index)"),
            transpose.to_owned(),
        ),
        (
            format!("cluster({E12},{sil},count=11..11,output=index)"),
            transpose.to_owned(),
        ),
        (
            format!("cluster({E12},index=\"dunn\",output=index,count=2)"),
            transpose.to_owned(),
        ),
        (
            r#"cluster(["x"],index="silhouette",output=index)"#.to_owned(),
            transpose.to_owned(),
        ),
        (
            r#"cluster(["a","a","a"],index="silhouette",count_sup=2,output=index)"#.to_owned(),
            transpose.to_owned(),
        ),
    ];
    stay_rows(&rows);
}

/// giac prints each recorded error; ours stays with its message.
fn stay_rows(rows: &[(String, String)]) {
    let giac = giac();
    let mut wrong = Vec::new();
    for (input, error) in rows {
        if let Some(giac) = giac {
            let printed = giac.eval(input).unwrap_or_default();
            if collapsed(&printed) != collapsed(error) {
                wrong.push(format!("{input}: giac prints {printed}, recorded {error}"));
            }
        }
        let message = error
            .strip_prefix("GIAC_ERROR: ")
            .and_then(|m| m.split_once(" Error:"))
            .map_or("", |(m, _)| m.strip_suffix('.').unwrap_or(m));
        if !stays_with(input, Some(message)) {
            wrong.push(format!(
                "{input}: ours {:?}, wanted the call with {message:?}",
                evaluated(input)
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `kmeans`' `output=index` rows: a value, `None` for "not defined",
/// `ComplexInfinity` (giac's `infinity`) as an infinity.
fn kmeans_rows(e: &Expr) -> Option<Vec<(i64, Option<f64>)>> {
    e.args()
        .iter()
        .map(|row| match row.args() {
            [k, v] => {
                let v = match v.view() {
                    View::Ball(b) => Some(b.mid().to_f64()),
                    _ if v.as_str() == Some("not defined") => None,
                    _ if v.as_symbol() == Some(Symbol::COMPLEX_INFINITY) => Some(f64::INFINITY),
                    _ => return None,
                };
                Some((k.as_i64()?, v))
            }
            _ => None,
        })
        .collect()
}

/// A value within 1e-9 (giac prints 12 digits), the same infinity, or both
/// undefined.
fn near_value(a: Option<f64>, b: Option<f64>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => a.total_cmp(&b).is_eq() || (a - b).abs() <= 1e-9 * a.abs().max(1.0),
        (None, None) => true,
        _ => false,
    }
}

/// `kmeans`' seed-independent `output=index` rows (the prototype's one
/// answer over 200 seeds): giac's, ours and the recorded one alike.
#[test]
fn kmeans_index_values() {
    let not = "\"not defined\"";
    let calinski = "[[2,450.0],[3,246.136363636],[4,201.111111111],[5,151.083333333]]";
    let ray = "[[2,0.00222222222222],[3,0.244444444444],[4,0.133333333333],[5,0.0833333333333]]";
    let banfeld = format!("[[2,-4.8655812973],[3,{not}],[4,{not}],[5,{not}]]");
    let trace = "[[2,2.66666666667],[3,1.83333333333],[4,1.0],[5,0.5]]";
    let rows = [
        (format!("kmeans({S6},index=\"calinski-harabasz\",output=index)"), calinski.to_owned()),
        (format!("kmeans({S6},index=true,output=index)"), calinski.to_owned()),
        (
            "kmeans([[0.0,0],[0,1],[1,0],[10,10],[10,11],[11,10]],index=\"calinski-harabasz\",output=index)"
                .to_owned(),
            calinski.to_owned(),
        ),
        (format!("kmeans({S6},index=\"ray-turi\",output=index)"), ray.to_owned()),
        (format!("kmeans({S6},index=[\"pbm\",\"ray-turi\"],output=index)"), ray.to_owned()),
        (format!("kmeans({S6},index=\"banfeld-raftery\",output=index)"), banfeld.clone()),
        (format!("kmeans({S6},index=\"all\",output=index)"), banfeld),
        (
            format!("kmeans({S6},index=\"scott-symons\",output=index)"),
            format!("[[2,-11.4572550293],[3,{not}],[4,{not}],[5,{not}]]"),
        ),
        (
            format!("kmeans({S6},index=\"ratkowsky-lance\",output=index)"),
            "[[2,0.9955849674],[3,0.996966765383],[4,0.998346650839],[5,0.999173667398]]".to_owned(),
        ),
        (
            format!("kmeans({S6},index=\"log-ss\",output=index)"),
            "[[2,4.72295322164],[3,5.10042059797],[4,5.70932265503],[5,6.40412592209]]".to_owned(),
        ),
        (
            format!("kmeans({S6},index=\"ball-hall\",output=index)"),
            "[[2,0.444444444444],[3,0.231481481481],[4,0.125],[5,0.05]]".to_owned(),
        ),
        (format!("kmeans({S6},index=\"trace-W\",output=index)"), trace.to_owned()),
        (format!("kmeans(matrix({S6}),index=\"trace-W\",output=index)"), trace.to_owned()),
        (
            "kmeans([[0,0],[0,1],[1/2,0],[10,10],[10,11],[21/2,10]],index=\"trace-W\",output=index)"
                .to_owned(),
            "[[2,1.66666666667],[3,0.958333333333],[4,0.25],[5,0.125]]".to_owned(),
        ),
    ];
    let giac = giac();
    let mut wrong = Vec::new();
    let same = |a: &Option<Vec<(i64, Option<f64>)>>, b: &Option<Vec<(i64, Option<f64>)>>| {
        matches!((a, b), (Some(a), Some(b)) if a.len() == b.len()
            && a.iter().zip(b).all(|((k, v), (l, w))| k == l && near_value(*v, *w)))
    };
    for (input, answer) in &rows {
        let expected = ours(answer).ok().as_ref().and_then(kmeans_rows);
        let mine = ours(input).ok().as_ref().and_then(kmeans_rows);
        if expected.is_none() || !same(&mine, &expected) {
            wrong.push(format!("{input}: ours {mine:?}, wanted {answer}"));
        }
        if let Some(giac) = giac {
            let theirs = giac.value(input).ok().as_ref().and_then(kmeans_rows);
            if !same(&theirs, &expected) {
                wrong.push(format!("{input}: giac {theirs:?}, recorded {answer}"));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `kmeans`' seed-independent counts and labels with `index=`: the `k` the
/// indices suggest (none for an elbow over two `k`: no labels, count 0);
/// the last of `k`, `count=k` and `index=` wins.
#[test]
fn kmeans_index_counts() {
    let count = |index: &str, rest: &str| format!("kmeans({S6},index={index}{rest},output=count)");
    let mut rows = vec![
        (count("\"calinski-harabasz\"", ""), "2"),
        (count("\"davies-bouldin\"", ""), "2"),
        (count("\"pbm\"", ""), "5"),
        (count("\"ratkowsky-lance\"", ""), "5"),
        (count("\"ksq-detW\"", ""), "4"),
        (count("\"log-ss\"", ""), "4"),
        (count("\"ball-hall\"", ""), "3"),
        (count("\"trace-W\"", ""), "4"),
        (count("\"all\"", ""), "2"),
        (count("false", ""), "2"),
        (count("\"pbm\"", ",3"), "3"),
        (format!("kmeans({S6},2,index=\"pbm\",output=count)"), "5"),
        (count("\"ball-hall\"", ",count=2..3"), "0"),
        (format!("kmeans({S6},index=\"ball-hall\",count=2..3)"), "[]"),
        (
            format!("kmeans({S6},index=\"ball-hall\",count=2..3,output=part)"),
            "[]",
        ),
        (
            format!("kmeans({S6},count=2..2,index=\"pbm\",output=index)"),
            "[]",
        ),
        (format!("kmeans({S6},count=1..3,index=\"pbm\")"), "[]"),
    ];
    for (index, k) in [
        ("davies-bouldin", "4"),
        ("ratkowsky-lance", "4"),
        ("det", "3"),
        ("log-det", "3"),
        ("ksq-detW", "3"),
        ("log-ss", "3"),
        ("ball-hall", "3"),
        ("trace-W", "3"),
        ("trace-WiB", "3"),
        ("all", "3"),
    ] {
        rows.push((
            format!("kmeans({E12},count=2..4,index=\"{index}\",output=count)"),
            k,
        ));
    }
    for index in ["det", "ksq-detW", "ball-hall", "trace-WiB"] {
        rows.push((
            format!("kmeans({E12},count=3..4,index=\"{index}\",output=count)"),
            "0",
        ));
    }
    answer_rows(&rows);
    let giac = giac();
    let mut wrong = Vec::new();
    let partition = |e: Option<Expr>| e.and_then(|e| labels(&e)).map(|l| canonical(&l));
    for index in [
        "\"banfeld-raftery\"",
        "\"ray-turi\"",
        "\"scott-symons\"",
        "\"calinski-harabasz\"",
        "true",
        "\"all\"",
    ] {
        let input = format!("kmeans({S6},index={index})");
        let mine = partition(ours(&input).ok());
        if mine.as_deref() != Some(&[0, 0, 0, 1, 1, 1][..]) {
            wrong.push(format!("{input}: ours {mine:?}"));
        }
        if let Some(giac) = giac {
            let theirs = partition(giac.value(&input).ok());
            if theirs != mine {
                wrong.push(format!("{input}: giac {theirs:?}"));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Departures, `kmeans`' seed-dependent `index=` rows: each `[k, value]`
/// of giac's answer and ours among the prototype's at that `k` over 200
/// seeds (giac's recorded answers in `probe_clustering3.out`).
#[test]
fn kmeans_index_departures() {
    type Allowed = &'static [(i64, &'static [Option<f64>])];
    let det: Allowed = &[
        (2, &[Some(338.25)]),
        (3, &[Some(4.005_921_539_6)]),
        (4, &[None, Some(2_405.333_333_33)]),
        (5, &[None]),
    ];
    let rows: [(&str, Allowed); 6] = [
        ("det", det),
        ("log-det", det),
        (
            "ksq-detW",
            &[
                (2, &[Some(7.111_111_111_11)]),
                (3, &[Some(7.0)]),
                (4, &[Some(0.0), Some(4.0)]),
                (5, &[Some(0.0)]),
            ],
        ),
        (
            "trace-WiB",
            &[
                (2, &[Some(281.25)]),
                (3, &[Some(418.857_142_857)]),
                (4, &[Some(f64::INFINITY), Some(603.333_333_333)]),
                (5, &[Some(f64::INFINITY)]),
            ],
        ),
        (
            "davies-bouldin",
            &[
                (2, &[Some(0.062_853_936_105_5)]),
                (3, &[Some(0.165_571_716_426), Some(0.165_299_149_551)]),
                (4, &[Some(0.223_606_797_75)]),
                (5, &[Some(0.099_956_096_710_8), Some(0.100_131_499_76)]),
            ],
        ),
        (
            "pbm",
            &[
                (2, &[Some(644_112.5)]),
                (3, &[Some(616_181.899_806), Some(627_538.153_25)]),
                (4, &[Some(1_209_500.138_89), Some(1_203_774.694_44)]),
                (5, &[Some(3_239_227.448_89), Some(3_096_320.355_56)]),
            ],
        ),
    ];
    let giac = giac();
    let mut wrong = Vec::new();
    for (index, allowed) in rows {
        let input = format!("kmeans({S6},index=\"{index}\",output=index)");
        let among = |r: &Option<Vec<(i64, Option<f64>)>>| {
            r.as_ref().is_some_and(|r| {
                r.len() == allowed.len()
                    && r.iter()
                        .zip(allowed)
                        .all(|((k, v), (l, vs))| k == l && vs.iter().any(|w| near_value(*v, *w)))
            })
        };
        let mine = ours(&input).ok().as_ref().and_then(kmeans_rows);
        if !among(&mine) {
            wrong.push(format!("{input}: ours {mine:?}"));
        }
        if let Some(giac) = giac {
            let theirs = giac.value(&input).ok().as_ref().and_then(kmeans_rows);
            if !among(&theirs) {
                wrong.push(format!("{input}: giac {theirs:?}"));
            }
        }
    }
    for (index, counts) in [("trace-WiB", ["3", "4"]), ("det", ["0", "3"])] {
        let input = format!("kmeans({S6},index=\"{index}\",output=count)");
        let allowed: Vec<Option<Expr>> = counts.iter().map(|c| ours(c).ok()).collect();
        let mine = ours(&input).ok();
        if !allowed.contains(&mine) {
            wrong.push(format!("{input}: ours {mine:?}"));
        }
        if let Some(giac) = giac {
            let theirs = giac.value(&input).ok();
            if !allowed.contains(&theirs) {
                wrong.push(format!("{input}: giac {theirs:?}"));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `kmeans`' `index=` errors stay with giac's message.
#[test]
fn kmeans_index_errors_stay() {
    let value = |m: &str| format!("GIAC_ERROR: {m}. Error: Bad Argument Value");
    let kind = |m: &str| format!("GIAC_ERROR: {m}. Error: Bad Argument Type");
    let strings = "Expected a string or list of strings";
    let none = "No index specified for output";
    stay_rows(&[
        (
            format!("kmeans({S6},index=[\"foo\"])"),
            value("Unknown index name"),
        ),
        (
            format!("kmeans({S6},index=\"Ball-hall\")"),
            value("Unknown index name"),
        ),
        (
            format!("kmeans({E12},index=\"foo\")"),
            value("Unknown index name"),
        ),
        (format!("kmeans({S6},index=[])"), kind(strings)),
        (format!("kmeans({S6},index=5)"), kind(strings)),
        (
            format!("kmeans({S6},index=\"pbm\",count=2,output=index)"),
            value(none),
        ),
        (
            format!("kmeans({S6},index=\"pbm\",2,output=index)"),
            value(none),
        ),
        (
            format!("kmeans({E12},index=false,output=index)"),
            value(none),
        ),
    ]);
}
