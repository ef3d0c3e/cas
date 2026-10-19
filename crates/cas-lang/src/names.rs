//! Surface names: the name table and its locales (CAS-DOC-7 §5).
//!
//! The table maps the names the language reserves (`sin`, `pi`, `diff`, …) to
//! what they stand for, and gives each target the spelling printers use.
//! English is the base, in `names/en.txt`, whose header describes the format.
//! A locale is an overlay in the same format that adds names and may change
//! spellings, so translating the language means writing a data file.
//!
//! Every other name belongs to the user. One that starts with an ASCII
//! capital becomes a symbol with a `$` in front ([`escape_name`]). Every
//! builtin's internal name starts with an ASCII capital, so `E := m*c^2`
//! defines `$E` and never touches the builtin `E`, whatever builtins are
//! added later.

use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::sync::LazyLock;

use cas_core::Symbol;

use crate::lexer::lex;
use crate::literal::name;
use crate::syntax::SyntaxKind;

/// The English names, the base of every locale.
const ENGLISH: &str = include_str!("../names/en.txt");

/// What a surface name stands for.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Target {
    /// A builtin of cas-core.
    Builtin(Symbol),
    /// A builtin that the design pages name but cas-core does not define yet,
    /// by its internal name. Its names are reserved already; lowering reports
    /// that it is not available.
    Planned(Box<str>),
    /// The imaginary unit, which lowers to a complex number (CAS-51).
    Imaginary,
    /// A form that lowering or the session handles itself.
    Syntax(SyntaxForm),
}

/// The forms of [`Target::Syntax`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SyntaxForm {
    /// `call(f, a, …)`, always a call (CAS-DOC-7 §3.1).
    Call,
    /// `assume(…)`, run by the session, which keeps the child context.
    Assume,
    /// `purge(…)`, run by the session, which clears definitions.
    Purge,
}

impl SyntaxForm {
    fn parse(word: &str) -> Option<SyntaxForm> {
        Some(match word {
            "call" => SyntaxForm::Call,
            "assume" => SyntaxForm::Assume,
            "purge" => SyntaxForm::Purge,
            _ => return None,
        })
    }
}

/// A type written after `::` in a pattern (CAS-DOC-7 §6).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum TypeName {
    /// Anything: marks a pattern variable without constraining it.
    Any,
    /// An integer.
    Integer,
    /// A rational number, integers included.
    Rational,
    /// A number of any kind.
    Number,
    /// A real value, as far as the assumptions tell.
    Real,
    /// A complex value, as far as the assumptions tell.
    Complex,
    /// A list.
    List,
    /// A matrix: a list of rows of equal length.
    Matrix,
    /// A string.
    String,
    /// A symbol.
    Symbol,
}

impl TypeName {
    fn parse(word: &str) -> Option<TypeName> {
        Some(match word {
            "any" => TypeName::Any,
            "integer" => TypeName::Integer,
            "rational" => TypeName::Rational,
            "number" => TypeName::Number,
            "real" => TypeName::Real,
            "complex" => TypeName::Complex,
            "list" => TypeName::List,
            "matrix" => TypeName::Matrix,
            "string" => TypeName::String,
            "symbol" => TypeName::Symbol,
            _ => return None,
        })
    }
}

/// A name table: the English names, plus those of at most one locale.
///
/// ```
/// use cas_core::Symbol;
/// use cas_lang::{Names, Target};
///
/// let names = Names::english();
/// assert_eq!(names.target("ln"), Some(&Target::Builtin(Symbol::LOG)));
/// assert_eq!(names.target("log"), Some(&Target::Builtin(Symbol::LOG)));
/// assert_eq!(names.spelling(&Target::Builtin(Symbol::LOG)), Some("ln"));
/// assert!(names.is_reserved("pi") && !names.is_reserved("x"));
///
/// let french = Names::with_locale("[builtins]\nExpand: developper\n").unwrap();
/// assert_eq!(french.target("developper"), french.target("expand"));
/// assert_eq!(french.spelling(french.target("expand").unwrap()), Some("developper"));
/// ```
#[derive(Clone, Debug)]
pub struct Names {
    targets: HashMap<Box<str>, Target>,
    spellings: HashMap<Target, Box<str>>,
    types: HashMap<Box<str>, TypeName>,
    type_spellings: HashMap<TypeName, Box<str>>,
    internal: Vec<Symbol>,
}

impl Names {
    /// The English table, `names/en.txt`.
    ///
    /// # Panics
    ///
    /// If `names/en.txt` does not load, which the tests rule out.
    #[must_use]
    pub fn english() -> &'static Names {
        static ENGLISH_NAMES: LazyLock<Names> = LazyLock::new(|| {
            let mut names = Names::empty();
            if let Err(e) = names.load(ENGLISH, false) {
                panic!("names/en.txt: {e}");
            }
            names
        });
        &ENGLISH_NAMES
    }

    /// The English table with the overlay of a locale, in the format of
    /// `names/en.txt`.
    ///
    /// # Errors
    ///
    /// A [`NamesError`] with its line for a malformed overlay, a target the
    /// English table does not name, or a name that already stands for
    /// something else.
    pub fn with_locale(overlay: &str) -> Result<Names, NamesError> {
        let mut names = Names::english().clone();
        names.load(overlay, true)?;
        Ok(names)
    }

    /// What `name` stands for, if the table reserves it. `name` is as
    /// [`name`](crate::name) gives it: `pi`, not `π`.
    #[must_use]
    pub fn target(&self, name: &str) -> Option<&Target> {
        self.targets.get(name)
    }

    /// Whether the table reserves `name`. A reserved name cannot be assigned,
    /// defined or bound as a parameter, a loop variable or a local
    /// (CAS-DOC-7 L7).
    #[must_use]
    pub fn is_reserved(&self, name: &str) -> bool {
        self.targets.contains_key(name)
    }

    /// The spelling printers use for `target`; `None` for a target without
    /// a name, such as an internal builtin.
    #[must_use]
    pub fn spelling(&self, target: &Target) -> Option<&str> {
        self.spellings.get(target).map(AsRef::as_ref)
    }

    /// The type `name` stands for after `::`. Type names are a namespace of
    /// their own: `list` is a type and may also name a function or variable.
    #[must_use]
    pub fn type_name(&self, name: &str) -> Option<TypeName> {
        self.types.get(name).copied()
    }

    /// The spelling printers use for `ty`.
    #[must_use]
    pub fn type_spelling(&self, ty: TypeName) -> Option<&str> {
        self.type_spellings.get(&ty).map(AsRef::as_ref)
    }

    /// Every reserved name with its target, in no particular order.
    pub fn names(&self) -> impl Iterator<Item = (&str, &Target)> {
        self.targets
            .iter()
            .map(|(name, target)| (name.as_ref(), target))
    }

    /// Every type name with its type, in no particular order.
    pub fn type_names(&self) -> impl Iterator<Item = (&str, TypeName)> {
        self.types.iter().map(|(name, &ty)| (name.as_ref(), ty))
    }

    /// The builtins without a surface name, each reached through syntax or
    /// not at all.
    #[must_use]
    pub fn internal(&self) -> &[Symbol] {
        &self.internal
    }

    fn empty() -> Names {
        Names {
            targets: HashMap::new(),
            spellings: HashMap::new(),
            types: HashMap::new(),
            type_spellings: HashMap::new(),
            internal: Vec::new(),
        }
    }

    /// Adds the entries of `text`, the base table or a locale's overlay.
    fn load(&mut self, text: &str, locale: bool) -> Result<(), NamesError> {
        let builtins: HashMap<&str, Symbol> = Symbol::builtins().map(|s| (s.name(), s)).collect();
        let mut section = None;
        let mut names_seen = HashSet::new();
        let mut targets_seen = HashSet::new();
        let mut types_seen = HashSet::new();
        for (index, raw) in text.lines().enumerate() {
            let fail = |kind| NamesError {
                line: index + 1,
                kind,
            };
            let line = raw.split('#').next().unwrap_or_default().trim();
            if line.is_empty() {
                continue;
            }
            if line.starts_with('[') {
                let found = line
                    .strip_prefix('[')
                    .and_then(|l| l.strip_suffix(']'))
                    .and_then(Section::parse)
                    .ok_or_else(|| fail(NamesErrorKind::UnknownSection(line.into())))?;
                if locale && !found.in_locale() {
                    return Err(fail(NamesErrorKind::SectionNotInLocale(line.into())));
                }
                section = Some(found);
                continue;
            }
            let Some(section) = section else {
                return Err(fail(NamesErrorKind::OutsideSection));
            };
            if section == Section::Internal {
                for word in line.split_whitespace() {
                    let &symbol = builtins
                        .get(word)
                        .ok_or_else(|| fail(NamesErrorKind::UnknownTarget(word.into())))?;
                    if self.internal.contains(&symbol) {
                        return Err(fail(NamesErrorKind::Repeated(word.into())));
                    }
                    self.internal.push(symbol);
                }
                continue;
            }
            let (word, list) = line
                .split_once(':')
                .ok_or_else(|| fail(NamesErrorKind::MissingColon))?;
            let word = word.trim();
            let names: Vec<&str> = list.split_whitespace().collect();
            let Some(&first) = names.first() else {
                return Err(fail(NamesErrorKind::NoNames));
            };
            if let Some(bad) = names.iter().find(|n| !is_name(n)) {
                return Err(fail(NamesErrorKind::NotAName((*bad).into())));
            }
            if section == Section::Types {
                let ty = TypeName::parse(word)
                    .ok_or_else(|| fail(NamesErrorKind::UnknownTarget(word.into())))?;
                if !types_seen.insert(ty) {
                    return Err(fail(NamesErrorKind::Repeated(word.into())));
                }
                for &n in &names {
                    if !names_seen.insert((true, n)) {
                        return Err(fail(NamesErrorKind::Repeated(n.into())));
                    }
                    match self.types.get(n) {
                        Some(&other) if other != ty => {
                            return Err(fail(NamesErrorKind::Conflict(n.into())));
                        }
                        _ => {
                            self.types.insert(n.into(), ty);
                        }
                    }
                }
                self.type_spellings.insert(ty, first.into());
                continue;
            }
            let target = self
                .resolve(section, word, locale, &builtins)
                .map_err(fail)?;
            if !targets_seen.insert(target.clone()) {
                return Err(fail(NamesErrorKind::Repeated(word.into())));
            }
            for &n in &names {
                if !names_seen.insert((false, n)) {
                    return Err(fail(NamesErrorKind::Repeated(n.into())));
                }
                match self.targets.get(n) {
                    Some(other) if *other != target => {
                        return Err(fail(NamesErrorKind::Conflict(n.into())));
                    }
                    _ => {
                        self.targets.insert(n.into(), target.clone());
                    }
                }
            }
            self.spellings.insert(target, first.into());
        }
        Ok(())
    }

    /// The target `word` names in `section`.
    fn resolve(
        &self,
        section: Section,
        word: &str,
        locale: bool,
        builtins: &HashMap<&str, Symbol>,
    ) -> Result<Target, NamesErrorKind> {
        let unknown = || NamesErrorKind::UnknownTarget(word.into());
        match section {
            // A locale names any builtin the base names, defined or planned.
            Section::Builtins if locale => {
                let target = match builtins.get(word) {
                    Some(&symbol) => Target::Builtin(symbol),
                    None => Target::Planned(word.into()),
                };
                if self.spellings.contains_key(&target) {
                    Ok(target)
                } else {
                    Err(unknown())
                }
            }
            Section::Builtins => builtins
                .get(word)
                .map(|&symbol| Target::Builtin(symbol))
                .ok_or_else(unknown),
            Section::Planned if builtins.contains_key(word) => {
                Err(NamesErrorKind::AlreadyBuiltin(word.into()))
            }
            Section::Planned if is_internal_name(word) => Ok(Target::Planned(word.into())),
            Section::Constants if word == "imaginary" => Ok(Target::Imaginary),
            Section::Syntax => SyntaxForm::parse(word)
                .map(Target::Syntax)
                .ok_or_else(unknown),
            _ => Err(unknown()),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Section {
    Builtins,
    Planned,
    Constants,
    Syntax,
    Types,
    Internal,
}

impl Section {
    fn parse(name: &str) -> Option<Section> {
        Some(match name {
            "builtins" => Section::Builtins,
            "planned" => Section::Planned,
            "constants" => Section::Constants,
            "syntax" => Section::Syntax,
            "types" => Section::Types,
            "internal" => Section::Internal,
            _ => return None,
        })
    }

    /// Whether a locale's overlay may use the section.
    fn in_locale(self) -> bool {
        !matches!(self, Section::Planned | Section::Internal)
    }
}

/// Whether `text` is one identifier of the language, in normal form C.
fn is_name(text: &str) -> bool {
    let lexed = lex(text);
    lexed.len() == 1
        && lexed.kind(0) == SyntaxKind::Ident
        && lexed.errors().is_empty()
        && name(text) == text
}

/// Whether `word` follows the rule for builtins' internal names: an ASCII
/// capital, then ASCII letters and digits.
fn is_internal_name(word: &str) -> bool {
    word.starts_with(|c: char| c.is_ascii_uppercase())
        && word.chars().all(|c| c.is_ascii_alphanumeric())
}

/// Why a name table does not load.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NamesError {
    /// The line, from 1.
    pub line: usize,
    /// What is wrong.
    pub kind: NamesErrorKind,
}

/// What is wrong in a [`NamesError`].
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum NamesErrorKind {
    /// A `[section]` that does not exist.
    UnknownSection(String),
    /// `[planned]` or `[internal]` in a locale's overlay.
    SectionNotInLocale(String),
    /// An entry before any `[section]`.
    OutsideSection,
    /// An entry without the `:` after its target.
    MissingColon,
    /// A target with no name after its `:`.
    NoNames,
    /// A target the section does not know.
    UnknownTarget(String),
    /// A `[planned]` entry for a builtin cas-core defines: it belongs in
    /// `[builtins]`.
    AlreadyBuiltin(String),
    /// A name that is not one identifier in normal form C, or a keyword.
    NotAName(String),
    /// A name that already stands for another target.
    Conflict(String),
    /// A name or target given twice in one file.
    Repeated(String),
}

impl fmt::Display for NamesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "line {}: ", self.line)?;
        match &self.kind {
            NamesErrorKind::UnknownSection(s) => write!(f, "unknown section {s}"),
            NamesErrorKind::SectionNotInLocale(s) => {
                write!(f, "a locale cannot have the section {s}")
            }
            NamesErrorKind::OutsideSection => f.write_str("entry before any [section]"),
            NamesErrorKind::MissingColon => f.write_str("expected `Target: name …`"),
            NamesErrorKind::NoNames => f.write_str("no name after the `:`"),
            NamesErrorKind::UnknownTarget(t) => write!(f, "unknown target {t}"),
            NamesErrorKind::AlreadyBuiltin(t) => {
                write!(f, "{t} is a builtin now: move it to [builtins]")
            }
            NamesErrorKind::NotAName(n) => write!(f, "{n} is not a name"),
            NamesErrorKind::Conflict(n) => write!(f, "{n} already stands for something else"),
            NamesErrorKind::Repeated(n) => write!(f, "{n} is given twice"),
        }
    }
}

impl std::error::Error for NamesError {}

/// The symbol name a user's name lowers to: `$` in front when it starts with
/// an ASCII capital, so it cannot be a builtin's internal name.
///
/// ```
/// use cas_lang::{escape_name, unescape_name};
///
/// assert_eq!(escape_name("E"), "$E");
/// assert_eq!(escape_name("x"), "x");
/// assert_eq!(escape_name("Δ"), "Δ");
/// assert_eq!(unescape_name("$E"), "E");
/// ```
#[must_use]
pub fn escape_name(name: &str) -> Cow<'_, str> {
    if name.starts_with(|c: char| c.is_ascii_uppercase()) {
        Cow::Owned(format!("${name}"))
    } else {
        Cow::Borrowed(name)
    }
}

/// The user's name for a symbol name that [`escape_name`] gave.
#[must_use]
pub fn unescape_name(symbol_name: &str) -> &str {
    symbol_name.strip_prefix('$').unwrap_or(symbol_name)
}

#[cfg(test)]
mod tests;
