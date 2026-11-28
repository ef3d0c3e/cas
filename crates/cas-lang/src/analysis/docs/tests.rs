use super::{OVERLAY, Source, docs, kept, records};
use crate::names::{Names, Target};

/// The examples hidden at load, pinned from CI's first run.
const HIDDEN: usize = 0;

/// The first 20 hidden examples, in the order of the files.
const FIRST_HIDDEN: [&str; 0] = [];

#[test]
fn every_builtin_and_the_imaginary_unit_have_a_record() {
    let mut missing: Vec<&str> = Names::english()
        .names()
        .filter(|(_, t)| matches!(t, Target::Builtin(_) | Target::Imaginary))
        .filter(|(_, t)| docs().get(t).is_none())
        .map(|(name, _)| name)
        .collect();
    missing.sort_unstable();
    assert!(missing.is_empty(), "no record: {missing:?}");
}

#[test]
fn overlay_keys_are_named_targets() {
    let overlay = records(OVERLAY, "names/en_docs_overlay.txt", true);
    assert!(!overlay.is_empty());
    for (line, target, _) in &overlay {
        assert!(
            Names::english().spelling(target).is_some(),
            "overlay line {line}: {target:?} has no English name"
        );
        assert!(docs().get(target).is_some_and(|d| d.source != Source::Giac));
    }
}

#[test]
fn overlay_examples_are_kept() {
    for (line, target, fields) in records(OVERLAY, "names/en_docs_overlay.txt", true) {
        let doc = docs().get(&target).expect("loaded");
        for example in fields.examples.unwrap_or_default() {
            assert!(
                kept(&example, Names::english()),
                "overlay record of line {line}: {example}"
            );
            assert!(doc.examples.contains(&example), "{example}");
        }
    }
}

#[test]
fn see_keys_resolve_to_records() {
    for (target, doc) in docs().iter() {
        for see in &doc.see {
            assert!(docs().get(see).is_some(), "{target:?} sees {see:?}");
        }
    }
}

#[test]
fn hidden_examples_are_pinned() {
    let d = docs();
    let first: Vec<&str> = d.hidden.iter().take(20).map(String::as_str).collect();
    assert!(
        d.hidden_examples() == HIDDEN && first == FIRST_HIDDEN,
        "hidden {}, the first 20: {first:#?}",
        d.hidden_examples()
    );
}
