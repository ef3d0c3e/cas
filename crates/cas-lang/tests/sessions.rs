//! Scripted sessions (CAS-121): each `sessions/*.cas` runs in a new
//! session, and its transcript must be the `.out` file next to it. Run with
//! `BLESS=1` to write the transcripts instead, then read the diff.
//!
//! A transcript shows each statement (`>`, and `.` for its next lines), its
//! messages with the line, column and text they are about, and its result
//! after `=`.
#![allow(missing_docs, reason = "test binary")]

use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use cas_lang::{Level, Output, Session};

/// The line and column, from 1, of the byte `at` of `text`.
fn position(text: &str, at: usize) -> (usize, usize) {
    let before = &text[..at];
    let line = before.matches('\n').count() + 1;
    let column = before.rsplit('\n').next().unwrap_or("").chars().count() + 1;
    (line, column)
}

fn transcript(text: &str, outputs: &[Output]) -> String {
    let mut out = String::new();
    for o in outputs {
        for (i, line) in text[o.range.clone()].lines().enumerate() {
            let mark = if i == 0 { '>' } else { '.' };
            writeln!(out, "{mark} {line}").unwrap();
        }
        for m in &o.messages {
            let level = match m.level {
                Level::Error => "error",
                Level::Warning => "warning",
                Level::Note => "note",
            };
            let (line, column) = position(text, m.range.start);
            let about = &text[m.range.clone()];
            writeln!(out, "{level} at {line}:{column} `{about}`: {}", m.text).unwrap();
            if let Some(fix) = &m.fix {
                writeln!(out, "  fix: `{}`", fix.replacement).unwrap();
            }
        }
        if let Some(result) = &o.result {
            writeln!(out, "= {}", result.text()).unwrap();
        }
    }
    out
}

#[test]
fn sessions() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/sessions");
    let bless = std::env::var_os("BLESS").is_some();
    let mut scripts: Vec<_> = fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|e| e == "cas"))
        .collect();
    scripts.sort();
    assert!(!scripts.is_empty());
    let mut wrong = Vec::new();
    for script in scripts {
        let text = fs::read_to_string(&script).unwrap();
        let found = transcript(&text, &Session::new().run(&text));
        let expected_path = script.with_extension("out");
        if bless {
            fs::write(&expected_path, &found).unwrap();
            continue;
        }
        let expected = fs::read_to_string(&expected_path).unwrap_or_default();
        if found != expected {
            wrong.push(format!(
                "{}:\n--- expected\n{expected}--- found\n{found}",
                script.display()
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
