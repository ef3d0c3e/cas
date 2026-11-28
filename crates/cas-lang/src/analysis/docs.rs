//! The documentation of the builtins, by target (CAS-DOC-60 §6): the table
//! generated from giac's `aide_cas` (`names/en_docs.txt`, never hand-edited)
//! with the hand-edited overlay (`names/en_docs_overlay.txt`) merged over it.
//! A record of the overlay for a target the table lacks is new
//! ([`Source::Ours`]); one for a target it has replaces each field kind it
//! holds and adds its notes ([`Source::GiacAmended`]). At load, an example
//! is kept only when it parses with no error and lowers with
//! [`Names::english`] with no error diagnostic; the others are counted
//! ([`Docs::hidden_examples`]).

use std::collections::{HashMap, HashSet};
use std::sync::LazyLock;

use cas_core::Symbol;

use crate::lower::lower;
use crate::names::{Names, Target};
use crate::parser::parse;

/// The table generated from `aide_cas`.
const GENERATED: &str = include_str!("../../names/en_docs.txt");

/// The hand-edited overlay.
const OVERLAY: &str = include_str!("../../names/en_docs_overlay.txt");

/// The documentation of the builtins.
#[derive(Debug)]
pub struct Docs {
    docs: HashMap<Target, Doc>,
    /// The examples dropped at load, in the order of the files.
    hidden: Vec<String>,
}

/// The documentation of one target.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Doc {
    /// Where the record comes from.
    pub source: Source,
    /// The names of the `aide_cas` entry, its credit; empty for ours.
    pub giac: Vec<String>,
    /// The arguments, as `aide_cas` writes them.
    pub signature: Option<String>,
    /// One sentence on what it does.
    pub summary: Option<String>,
    /// Related targets.
    pub see: Vec<Target>,
    /// Example inputs that our parser and lowering accept.
    pub examples: Vec<String>,
    /// Remarks of ours, from the overlay.
    pub notes: Vec<String>,
    /// The departures from giac recorded for it.
    pub departures: Vec<Departure>,
}

/// A departure from giac, as `crates/cas-giac/README.md` (`## Tests`) records it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Departure {
    /// The README's line.
    pub line: u32,
    /// The README's text.
    pub text: String,
}

/// Where a [`Doc`] comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Source {
    /// The generated table, from `aide_cas`.
    Giac,
    /// The generated table with fields of the overlay.
    GiacAmended,
    /// The overlay alone.
    Ours,
}

impl Docs {
    /// The documentation of `target`.
    #[must_use]
    pub fn get(&self, target: &Target) -> Option<&Doc> {
        self.docs.get(target)
    }

    /// Every target with its documentation, in no particular order.
    pub fn iter(&self) -> impl Iterator<Item = (&Target, &Doc)> {
        self.docs.iter()
    }

    /// The number of examples dropped at load.
    #[must_use]
    pub fn hidden_examples(&self) -> usize {
        self.hidden.len()
    }
}

/// The English documentation, loaded once.
///
/// # Panics
///
/// If `names/en_docs.txt` or `names/en_docs_overlay.txt` does not load,
/// which the tests rule out.
#[must_use]
pub fn docs() -> &'static Docs {
    static DOCS: LazyLock<Docs> = LazyLock::new(load);
    &DOCS
}

/// The fields of a record as a file gives them; `None` for a kind it does
/// not hold.
#[derive(Default)]
struct Fields {
    giac: Option<Vec<String>>,
    signature: Option<String>,
    summary: Option<String>,
    see: Option<Vec<Target>>,
    examples: Option<Vec<String>>,
    notes: Vec<String>,
    departures: Option<Vec<Departure>>,
}

impl Fields {
    fn into_doc(self, source: Source) -> Doc {
        Doc {
            source,
            giac: self.giac.unwrap_or_default(),
            signature: self.signature,
            summary: self.summary,
            see: self.see.unwrap_or_default(),
            examples: self.examples.unwrap_or_default(),
            notes: self.notes,
            departures: self.departures.unwrap_or_default(),
        }
    }

    /// Replaces each field kind of `doc` that `self` holds and adds the notes.
    fn amend(self, doc: &mut Doc) {
        doc.source = Source::GiacAmended;
        if let Some(giac) = self.giac {
            doc.giac = giac;
        }
        if self.signature.is_some() {
            doc.signature = self.signature;
        }
        if self.summary.is_some() {
            doc.summary = self.summary;
        }
        if let Some(see) = self.see {
            doc.see = see;
        }
        if let Some(examples) = self.examples {
            doc.examples = examples;
        }
        doc.notes.extend(self.notes);
        if let Some(departures) = self.departures {
            doc.departures = departures;
        }
    }
}

/// The target a record key names: a builtin by internal name, or
/// `imaginary`.
fn resolve(key: &str, builtins: &HashMap<&str, Symbol>) -> Option<Target> {
    if key == "imaginary" {
        Some(Target::Imaginary)
    } else {
        builtins.get(key).map(|&s| Target::Builtin(s))
    }
}

/// The records of `text`, in order, each with its header's line; `file`
/// names it in panics. Only an overlay may hold `note:`.
fn records(text: &str, file: &str, overlay: bool) -> Vec<(usize, Target, Fields)> {
    let builtins: HashMap<&str, Symbol> = Symbol::builtins().map(|s| (s.name(), s)).collect();
    let mut records: Vec<(usize, Target, Fields)> = Vec::new();
    for (index, line) in text.lines().enumerate() {
        let number = index + 1;
        let fail = |message: String| -> ! { panic!("{file}:{number}: {message}") };
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some(key) = line.strip_prefix('[').and_then(|l| l.strip_suffix(']')) {
            let target =
                resolve(key, &builtins).unwrap_or_else(|| fail(format!("`{key}` is not a target")));
            records.push((number, target, Fields::default()));
            continue;
        }
        let Some((_, _, fields)) = records.last_mut() else {
            fail("a field before any record".into());
        };
        let Some((kind, value)) = line.split_once(':') else {
            fail("no `kind:`".into());
        };
        let value = value.strip_prefix(' ').unwrap_or(value).to_owned();
        let single = |slot: &mut Option<String>| {
            if slot.replace(value.clone()).is_some() {
                fail(format!("a second `{kind}:`"));
            }
        };
        match kind {
            "giac" => {
                let names = value.split_whitespace().map(str::to_owned).collect();
                if fields.giac.replace(names).is_some() {
                    fail("a second `giac:`".into());
                }
            }
            "sig" => single(&mut fields.signature),
            "doc" => single(&mut fields.summary),
            "see" => {
                let see = value
                    .split_whitespace()
                    .map(|k| {
                        resolve(k, &builtins)
                            .unwrap_or_else(|| fail(format!("see `{k}` is not a target")))
                    })
                    .collect();
                if fields.see.replace(see).is_some() {
                    fail("a second `see:`".into());
                }
            }
            "ex" => fields.examples.get_or_insert_with(Vec::new).push(value),
            "note" if overlay => fields.notes.push(value),
            "dep" => {
                let departure = value
                    .split_once(' ')
                    .and_then(|(l, t)| {
                        Some(Departure {
                            line: l.parse().ok()?,
                            text: t.to_owned(),
                        })
                    })
                    .unwrap_or_else(|| fail("`dep:` is not `<line> <text>`".into()));
                fields
                    .departures
                    .get_or_insert_with(Vec::new)
                    .push(departure);
            }
            _ => fail(format!("unknown field `{kind}:`")),
        }
    }
    records
}

/// Whether `example` parses with no error and lowers with no error
/// diagnostic.
fn kept(example: &str, names: &Names) -> bool {
    let parsed = parse(example);
    parsed.errors().is_empty() && !lower(&parsed, names).has_errors()
}

fn load() -> Docs {
    let mut docs: HashMap<Target, Doc> = HashMap::new();
    let mut order = Vec::new();
    for (line, target, fields) in records(GENERATED, "names/en_docs.txt", false) {
        if docs
            .insert(target.clone(), fields.into_doc(Source::Giac))
            .is_some()
        {
            panic!("names/en_docs.txt:{line}: a second record for {target:?}");
        }
        order.push(target);
    }
    let mut seen = HashSet::new();
    for (line, target, fields) in records(OVERLAY, "names/en_docs_overlay.txt", true) {
        assert!(
            seen.insert(target.clone()),
            "names/en_docs_overlay.txt:{line}: a second record for {target:?}"
        );
        if let Some(doc) = docs.get_mut(&target) {
            fields.amend(doc);
        } else {
            docs.insert(target.clone(), fields.into_doc(Source::Ours));
            order.push(target);
        }
    }
    let names = Names::english();
    let mut hidden = Vec::new();
    for target in &order {
        if let Some(doc) = docs.get_mut(target) {
            doc.examples.retain(|example| {
                let keep = kept(example, names);
                if !keep {
                    hidden.push(example.clone());
                }
                keep
            });
        }
    }
    Docs { docs, hidden }
}

#[cfg(test)]
mod tests;
