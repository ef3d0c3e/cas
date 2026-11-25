//! giac's `gen2markup` (markup.cc :1233-3627), its LaTeX and MathML
//! (presentation and content) arms: a [`Gen`] to a [`MarkupBlock`]. Each
//! function cites the lines it ports; the prototypes
//! `.run/1.6/T4/proto/markup.py` (1.6/T4.c1, LaTeX) and
//! `.run/1.6/T5/proto/markup.py` (1.6/T5.c1, MathML) were checked against
//! giac 2.0.0.21's `latex()` on 301 inputs and `export_mathml` on 893
//! (301 inputs in its three forms), and these are their steps.
//!
//! One run writes the texts its [`Flags`] select. MathML content elements
//! take ids from one counter, post-order (children first), and the
//! presentation tag of the same node refers to it (`xref`): giac counts
//! only in content (`++idc`), so a presentation-only run has no ids.

use std::cell::Cell;

use cas_core::Integer;

use super::ident::{
    func2markup, func2mml, idnt2markup, idnt2mml, is_elemfunc, is_greek_letter, my_isalpha,
    string2markup, string2mml,
};
use super::mml;
use super::tree::{Gen, Sub};

/// Operator priorities (:65-80): 0 is an atom, larger binds weaker.
pub(super) mod priority {
    pub(in super::super) const APPLY: u8 = 1;
    pub(in super::super) const UNARY: u8 = 2;
    pub(in super::super) const EXP: u8 = 3;
    pub(in super::super) const MUL: u8 = 4;
    pub(in super::super) const ADD: u8 = 5;
    pub(in super::super) const SET: u8 = 6;
    pub(in super::super) const INEQ: u8 = 7;
    pub(in super::super) const EQ: u8 = 8;
    pub(in super::super) const COMP: u8 = 9;
    pub(in super::super) const ASSGN: u8 = 15;
}

/// Block kinds (`_MLBLOCK_*`, :83-101), bits of [`MarkupBlock::kind`].
pub(super) mod kind {
    pub(in super::super) const GENERAL: u32 = 0;
    pub(in super::super) const NUMERIC_EXACT: u32 = 1;
    pub(in super::super) const NUMERIC_APPROX: u32 = 2;
    pub(in super::super) const SUBTYPE_IDNT: u32 = 4;
    pub(in super::super) const FRACTION: u32 = 8;
    pub(in super::super) const ROOT: u32 = 16;
    pub(in super::super) const NEG: u32 = 32;
    pub(in super::super) const SUMPROD: u32 = 64;
    pub(in super::super) const FACTORIAL: u32 = 128;
    pub(in super::super) const MATRIX: u32 = 256;
    pub(in super::super) const ELEMAPP: u32 = 512;
    pub(in super::super) const FUNC: u32 = 1024;
    pub(in super::super) const POWER: u32 = 2048;
    pub(in super::super) const DERIVATIVE: u32 = 4096;
    pub(in super::super) const HAS_SUBSCRIPT: u32 = 8192;
    pub(in super::super) const LEADING_DIGIT: u32 = 16384;
    pub(in super::super) const IDNT_NAME: u32 = 32768;
}

/// The flags of a conversion (`_MARKUP_*`, :103-115).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Flags(u32);

impl Flags {
    pub(super) const LATEX: Flags = Flags(1);
    pub(super) const TOPLEVEL: Flags = Flags(2);
    pub(super) const MATHML_PRESENTATION: Flags = Flags(8);
    pub(super) const MATHML_CONTENT: Flags = Flags(16);
    pub(super) const ELEMPOW: Flags = Flags(32);
    pub(super) const FACTOR: Flags = Flags(64);

    fn has(self, f: Flags) -> bool {
        self.0 & f.0 != 0
    }

    fn tex(self) -> bool {
        self.has(Flags::LATEX)
    }

    fn pres(self) -> bool {
        self.has(Flags::MATHML_PRESENTATION)
    }

    fn cont(self) -> bool {
        self.has(Flags::MATHML_CONTENT)
    }

    fn without(self, f: Flags) -> Flags {
        Flags(self.0 & !f.0)
    }
}

impl std::ops::BitOr for Flags {
    type Output = Flags;

    fn bitor(self, f: Flags) -> Flags {
        Flags(self.0 | f.0)
    }
}

/// giac's `MarkupBlock` (:179-192), its LaTeX and MathML fields.
#[derive(Clone, Debug)]
pub(super) struct MarkupBlock {
    /// How weakly the text binds: 0 for an atom, see [`priority`].
    pub(super) priority: u8,
    /// [`kind`] bits.
    pub(super) kind: u32,
    /// A minus sign that belongs in front, not in the texts yet.
    pub(super) neg: bool,
    pub(super) latex: String,
    /// MathML presentation.
    pub(super) markup: String,
    /// MathML content.
    pub(super) content: String,
    /// Where an elementary function's application splits (`\sin |x`), in
    /// the markup text after its first tag (a byte length): 0 in a run
    /// without presentation, -1 for none.
    pub(super) split_pos: i32,
    /// The same split point in `latex`, a byte offset.
    pub(super) split_pos_tex: usize,
    /// An implicit application: `\sin x`.
    pub(super) appl: bool,
}

impl MarkupBlock {
    fn new() -> MarkupBlock {
        MarkupBlock {
            priority: 0,
            kind: kind::GENERAL,
            neg: false,
            latex: String::new(),
            markup: String::new(),
            content: String::new(),
            split_pos: -1,
            split_pos_tex: 0,
            appl: false,
        }
    }

    fn ctype(&self, t: u32) -> bool {
        self.kind & t != 0
    }

    /// The bytes of its texts, LaTeX, presentation and content together.
    fn len(&self) -> usize {
        self.latex.len() + self.markup.len() + self.content.len()
    }

    /// Another block's texts appended to its own.
    fn push(&mut self, other: &MarkupBlock) {
        self.latex.push_str(&other.latex);
        self.markup.push_str(&other.markup);
        self.content.push_str(&other.content);
    }
}

/// The text grew past the length allowed.
#[derive(Debug)]
pub(super) struct TooLong;

type Block = Result<MarkupBlock, TooLong>;

/// :958: `\left(…\right)`, `<mfenced>…</mfenced>`, an atom now.
fn parenthesize(ml: &mut MarkupBlock, flags: Flags) {
    if flags.tex() {
        ml.latex = format!("\\left({}\\right)", ml.latex);
    }
    if flags.pres() {
        ml.markup = format!("<mfenced>{}</mfenced>", ml.markup);
    }
    ml.priority = 0;
}

/// :972: the pending minus written in front; a product now, of kind `NEG`.
/// MathML: `<apply><minus/>…</apply>` and `<mrow><mo>&minus;</mo>…</mrow>`,
/// the operand's id moved onto them.
pub(super) fn prepend_minus(ml: &mut MarkupBlock, flags: Flags) {
    minus(ml, flags, true);
}

/// :972 with `mrow` false, as a sum's terms have it: the content in an
/// `<apply>` without id, the presentation sign bare before the term.
fn prepend_minus_term(ml: &mut MarkupBlock, flags: Flags) {
    minus(ml, flags, false);
}

fn minus(ml: &mut MarkupBlock, flags: Flags, mrow: bool) {
    if !ml.neg {
        return;
    }
    if ml.priority > priority::MUL {
        parenthesize(ml, flags);
    }
    if flags.tex() {
        ml.latex.insert(0, '-');
    }
    if flags.cont() {
        ml.content = if mrow {
            let id = mml::extract_id(&mut ml.content, true);
            let key = if id.is_empty() { "" } else { "id" };
            mml::tag("apply", &format!("<minus/>{}", ml.content), 0, &[(key, &id)])
        } else {
            format!("<apply><minus/>{}</apply>", ml.content)
        };
    }
    if flags.pres() {
        ml.markup = if mrow {
            let id = mml::extract_id(&mut ml.markup, false);
            let key = if id.is_empty() { "" } else { "xref" };
            mml::tag("mrow", &format!("{}{}", mml::MINUS, ml.markup), 0, &[(key, &id)])
        } else {
            format!("{}{}", mml::MINUS, ml.markup)
        };
    }
    ml.neg = false;
    ml.priority = priority::MUL;
    ml.kind = kind::NEG;
}

/// :936: the operands of an n-ary operator, nested ones of the same
/// operator flattened.
fn flatten_operands(g: &Gen) -> Vec<Gen> {
    let Gen::Symb(op, arg) = g else {
        return vec![g.clone()];
    };
    let Gen::Vect(items, _) = &**arg else {
        return vec![(**arg).clone()];
    };
    let mut ops = Vec::new();
    for x in items {
        if matches!(x, Gen::Symb(o, _) if o == op) {
            ops.extend(flatten_operands(x));
        } else {
            ops.push(x.clone());
        }
    }
    ops
}

/// :1053: `derive(…derive(f, k1)…, kn)` of an identifier `f`: `f` and the
/// orders `k + 1`, innermost first.
fn is_partialdiff(g: &Gen) -> Option<(&Gen, Vec<Gen>)> {
    let mut df = g;
    let mut deg = Vec::new();
    while let Gen::Symb(op, arg) = df
        && op == "derive"
        && let Gen::Vect(a, _) = &**arg
    {
        let [f, Gen::Int(k)] = a.as_slice() else {
            return None;
        };
        if k.is_negative() {
            return None;
        }
        deg.push(Gen::Int(k + &Integer::from(1)));
        df = f;
    }
    if !matches!(df, Gen::Idnt(_)) || deg.is_empty() {
        return None;
    }
    deg.reverse();
    Some((df, deg))
}

/// :1069: the variables of `derive(e, vars…)` with their counts, repeats
/// merged. `false` when one is not a variable.
fn get_derive_vars(args: &[Gen], vars: &mut Vec<(Gen, u64)>) -> bool {
    if let [_, v @ Gen::Idnt(_), Gen::Int(k)] = args
        && let Some(k) = k
            .to_i64()
            .and_then(|k| u64::try_from(k).ok())
            .filter(|&k| k > 0)
    {
        match vars.last_mut() {
            Some((last, n)) if last == v => *n += k,
            _ => vars.push((v.clone(), k)),
        }
        return true;
    }
    for g in args.iter().skip(1) {
        if g.is_sym(&["dollar"]) {
            let a = g.args();
            let (Some(h @ Gen::Idnt(_)), Some(Gen::Int(m))) = (a.first(), a.last()) else {
                return false;
            };
            let Some(m) = m
                .to_i64()
                .and_then(|m| u64::try_from(m).ok())
                .filter(|&m| m > 0)
            else {
                return false;
            };
            if let Some((last, n)) = vars.last_mut()
                && last == h
            {
                *n += m;
            }
            // giac pushes even after merging (:1086-1087).
            vars.push((h.clone(), m));
        } else if let Gen::Idnt(_) = g {
            match vars.last_mut() {
                Some((last, n)) if last == g => *n += 1,
                _ => vars.push((g.clone(), 1)),
            }
        } else if let Gen::Vect(items, _) = g {
            // The recursion skips the first item (giac quirk).
            if !get_derive_vars(items, vars) {
                return false;
            }
        } else {
            return false;
        }
    }
    true
}

/// The argument of `inv(…)`.
fn inv_arg(x: &Gen) -> Option<&Gen> {
    match x {
        Gen::Symb(op, arg) if op == "inv" => Some(arg),
        _ => None,
    }
}

/// giac's `ckmatrix`: rows of one length > 0, no vector inside.
fn ckmatrix(items: &[Gen]) -> bool {
    let rows: Vec<&Vec<Gen>> = items
        .iter()
        .filter_map(|r| match r {
            Gen::Vect(row, _) => Some(row),
            _ => None,
        })
        .collect();
    let Some(first) = rows.first() else {
        return false;
    };
    let n = first.len();
    rows.len() == items.len()
        && n > 0
        && rows
            .iter()
            .all(|r| r.len() == n && !r.iter().any(|x| matches!(x, Gen::Vect(..))))
}

/// The header of giac's `<math>` (:3745).
const MATH: &str = "<math mode='display' xmlns='http://www.w3.org/1998/Math/MathML'>";

/// The converter, with the longest text allowed: every block is checked as
/// it is built (its LaTeX, presentation and content together), so a long
/// text is refused before it is all made.
pub(super) struct Engine {
    max_len: usize,
    /// giac's `idc`: the last content id given.
    idc: Cell<u32>,
}

impl Engine {
    pub(super) fn new(max_len: usize) -> Engine {
        Engine {
            max_len,
            idc: Cell::new(0),
        }
    }

    fn check(&self, n: usize) -> Result<(), TooLong> {
        if n > self.max_len {
            Err(TooLong)
        } else {
            Ok(())
        }
    }

    /// giac's `++idc`: the next content id.
    fn next_id(&self) -> u32 {
        self.idc.set(self.idc.get() + 1);
        self.idc.get()
    }

    /// giac's `idc` as it is: the id of the node just made.
    fn id(&self) -> u32 {
        self.idc.get()
    }

    /// `g` at the top level with `flags`, its minus written; ids from 1.
    fn top(&self, g: &Gen, flags: Flags) -> Block {
        self.idc.set(0);
        let flags = flags | Flags::TOPLEVEL;
        let mut ml = self.gen2markup(g, flags)?;
        prepend_minus(&mut ml, flags);
        self.check(ml.len())?;
        Ok(ml)
    }

    /// :3629 `export_latex`.
    pub(super) fn export_latex(&self, g: &Gen) -> Result<String, TooLong> {
        Ok(self.top(g, Flags::LATEX | Flags::ELEMPOW)?.latex)
    }

    /// :3765 `export_mathml`: presentation with `xref`s, content with ids,
    /// and `printed`, the input form, as the annotation (escaped).
    pub(super) fn export_mathml(&self, g: &Gen, printed: &str) -> Result<String, TooLong> {
        let ml = self.top(
            g,
            Flags::MATHML_PRESENTATION | Flags::MATHML_CONTENT | Flags::ELEMPOW,
        )?;
        self.check(ml.len() + printed.len())?;
        self.whole(format!(
            "{MATH}<semantics>{}<annotation-xml encoding='MathML-Content'>{}</annotation-xml>\
             <annotation encoding='giac'>{}</annotation></semantics></math>",
            ml.markup,
            ml.content,
            mml::str_to_mml(printed)
        ))
    }

    /// :3755 `export_mathml_presentation`.
    pub(super) fn export_mathml_presentation(&self, g: &Gen) -> Result<String, TooLong> {
        let ml = self.top(g, Flags::MATHML_PRESENTATION | Flags::ELEMPOW)?;
        self.whole(format!("{MATH}{}</math>", ml.markup))
    }

    /// :3747 `export_mathml_content` (no `ELEMPOW`: content has no use
    /// for it).
    pub(super) fn export_mathml_content(&self, g: &Gen) -> Result<String, TooLong> {
        let ml = self.top(g, Flags::MATHML_CONTENT)?;
        self.whole(format!("{MATH}{}</math>", ml.content))
    }

    /// The whole text, checked against the cap.
    fn whole(&self, text: String) -> Result<String, TooLong> {
        self.check(text.len())?;
        Ok(text)
    }

    /// :1233 `gen2markup`, its texts checked against the cap.
    pub(super) fn gen2markup(&self, g: &Gen, flags: Flags) -> Block {
        let ml = self.markup(g, flags)?;
        self.check(ml.len())?;
        Ok(ml)
    }

    /// A number's MathML (:1327-1330, :1378-1399): `<cn type='…'>s</cn>`
    /// with the next id, `<mn>s</mn>`.
    fn number(&self, ml: &mut MarkupBlock, flags: Flags, s: &str, ty: &str) {
        if flags.cont() {
            ml.content = mml::tag("cn", s, self.next_id(), &[("type", ty)]);
        }
        if flags.pres() {
            ml.markup = mml::tag("mn", s, self.id(), &[]);
        }
    }

    /// The content of most operators: `<apply id=…>head a b</apply>` with
    /// the next id, made before the presentation takes that id.
    fn wrap(&self, ml: &mut MarkupBlock, flags: Flags, head: &str, a: &str, b: &str) {
        if flags.cont() {
            ml.content = mml::tag("apply", &format!("{head}{a}{b}"), self.next_id(), &[]);
        }
    }

    /// :2196-2201, :2233-2238: a square root, `<root/>`, `<msqrt>`,
    /// `\sqrt{…}`.
    fn root(&self, ml: &mut MarkupBlock, flags: Flags, tmp: &MarkupBlock) {
        self.wrap(ml, flags, "<root/>", &tmp.content, "");
        if flags.pres() {
            ml.markup = mml::tag("msqrt", &tmp.markup, self.id(), &[]);
        }
        if flags.tex() {
            ml.latex = format!("\\sqrt{{{}}}", tmp.latex);
        }
    }

    /// A child's block, its minus written.
    fn minus_prepended(&self, g: &Gen, flags: Flags) -> Block {
        let mut ml = self.gen2markup(g, flags)?;
        prepend_minus(&mut ml, flags);
        Ok(ml)
    }

    fn markup(&self, g: &Gen, flags_orig: Flags) -> Block {
        let toplevel = flags_orig.has(Flags::TOPLEVEL);
        let isfactor = flags_orig.has(Flags::FACTOR);
        let flags = flags_orig.without(Flags::TOPLEVEL).without(Flags::FACTOR);
        let mut ml = MarkupBlock::new();
        match g {
            // :1308-1321: content `<true>`, unclosed, as giac writes it.
            Gen::Bool(b) => {
                let name = if *b { "true" } else { "false" };
                if flags.cont() {
                    ml.content = mml::insert_id(&format!("<{name}>"), self.next_id(), true);
                }
                if flags.pres() {
                    ml.markup = idnt2mml(name, self.id()).0;
                }
                if flags.tex() {
                    ml.latex = idnt2markup(name).0;
                }
                ml.kind = kind::SUBTYPE_IDNT | kind::LEADING_DIGIT;
            }
            // :1323-1337.
            Gen::Int(n) => {
                ml.kind = kind::NUMERIC_EXACT | kind::LEADING_DIGIT;
                ml.neg = n.is_negative();
                let s = n.abs().to_string();
                self.number(&mut ml, flags, &s, "integer");
                if flags.tex() {
                    ml.latex = s;
                }
            }
            // :1362-1407.
            Gen::Dbl { text, neg } => {
                ml.kind = kind::NUMERIC_APPROX | kind::LEADING_DIGIT;
                ml.neg = *neg;
                self.number(&mut ml, flags, text, "double");
                if let Some((mant, minus, ex)) = double_parts(text) {
                    ml.priority = priority::MUL;
                    if flags.pres() {
                        let minus = if minus { mml::MINUS } else { "" };
                        ml.markup = mml::tag(
                            "mrow",
                            &format!(
                                "<mn>{mant}</mn><mo>&times;</mo><msup><mn>10</mn><mrow>\
                                 {minus}<mn>{ex}</mn></mrow></msup>"
                            ),
                            self.id(),
                            &[],
                        );
                    }
                    if flags.tex() {
                        let minus = if minus { "-" } else { "" };
                        ml.latex = format!("{mant}\\times10^{{{minus}{ex}}}");
                    }
                } else if flags.tex() {
                    ml.latex.clone_from(text);
                }
            }
            // :1408-1462.
            Gen::Cplx(re, im) => return self.cplx(re, im, flags),
            // :1463-1480.
            Gen::Frac(num, den) => {
                if let Gen::Cplx(_, im) = &**num
                    && !im.is_zero()
                {
                    let division = Gen::symb_n("division", vec![(**num).clone(), (**den).clone()]);
                    return self.gen2markup(&division, flags);
                }
                ml.kind = kind::FRACTION;
                ml.priority = priority::MUL;
                if let (Gen::Int(n), Gen::Int(d)) = (&**num, &**den) {
                    ml.neg = n.is_negative() != d.is_negative();
                    let (n, d) = (n.abs(), d.abs());
                    if flags.cont() {
                        let id = self.next_id();
                        ml.content = mml::tag(
                            "cn",
                            &format!("{n}<sep/>{d}"),
                            id,
                            &[("type", "rational")],
                        );
                    }
                    if flags.pres() {
                        let frac = format!("<mn>{n}</mn><mn>{d}</mn>");
                        ml.markup = mml::tag("mfrac", &frac, self.id(), &[]);
                    }
                    if flags.tex() {
                        ml.latex = format!("\\frac{{{n}}}{{{d}}}");
                    }
                }
            }
            // :1481-1485, :408-468.
            Gen::Str(s) => {
                if flags.cont() {
                    ml.content = mml::str_to_mml_tag(s, "cs", self.next_id());
                }
                if flags.pres() {
                    ml.markup = string2mml(s, self.id());
                }
                if flags.tex() {
                    ml.latex = string2markup(s);
                }
            }
            // :1647-1657.
            Gen::Func(name) => {
                if flags.cont() {
                    ml.content = func2mml(name, self.next_id(), true);
                }
                if flags.pres() {
                    ml.markup = func2mml(name, self.id(), false);
                }
                if flags.tex() {
                    ml.latex = func2markup(name);
                }
                ml.kind = kind::FUNC;
            }
            // :1486-1542.
            Gen::Idnt(name) => {
                let constant = match name.as_str() {
                    "infinity" | "unsigned_inf" => Some(("<infinity/>", "&infin;", "\\infty ")),
                    "pi" => Some(("<pi/>", "&pi;", "\\pi ")),
                    "euler_gamma" => Some(("<eulergamma/>", "&gamma;", "\\gamma ")),
                    _ => None,
                };
                if let Some((c, p, t)) = constant {
                    if flags.cont() {
                        ml.content = mml::insert_id(c, self.next_id(), true);
                    }
                    if flags.pres() {
                        ml.markup = mml::tag("mi", p, self.id(), &[]);
                    }
                    if flags.tex() {
                        ml.latex.push_str(t);
                    }
                    return Ok(ml);
                }
                if flags.cont() {
                    ml.content = mml::tag("ci", name, self.next_id(), &[]);
                }
                let (mut has_subscript, mut isname) = (false, false);
                if flags.pres() {
                    (ml.markup, isname) = idnt2mml(name, self.id());
                    // giac looks for `<msub>`, so an `<msub xref=…>` is
                    // none (:1527).
                    has_subscript = ml.markup.contains("<msub>");
                }
                if flags.tex() {
                    (ml.latex, isname) = idnt2markup(name);
                    has_subscript = ml.latex.contains("_{");
                }
                if has_subscript {
                    ml.kind |= kind::HAS_SUBSCRIPT;
                }
                if isname {
                    ml.kind |= kind::IDNT_NAME;
                }
            }
            // :1658-1773.
            Gen::Vect(items, sub) => return self.vect(items, *sub, toplevel, flags),
            Gen::Symb(..) => return self.symb(g, flags, isfactor),
        }
        Ok(ml)
    }

    /// :1408-1462: `\mathrm{i}`, `2 \mathrm{i}`, `1+2 \mathrm{i}`.
    fn cplx(&self, re: &Gen, im: &Gen, flags: Flags) -> Block {
        if im.is_zero() {
            return self.gen2markup(re, flags);
        }
        let mut ml = MarkupBlock::new();
        let isone = im.is_one() || im.is_minus_one();
        if re.is_zero() {
            if isone {
                // :1414-1417: the id as it is, no new one (giac).
                if flags.cont() {
                    ml.content = mml::insert_id("<imaginaryi/>", self.id(), true);
                }
                if flags.pres() {
                    ml.markup = mml::insert_id(mml::I, self.id(), false);
                }
                if flags.tex() {
                    ml.latex.push_str("\\mathrm{i}");
                }
                ml.neg = im.is_minus_one();
            } else {
                let mut tmp = self.gen2markup(im, flags | Flags::FACTOR)?;
                ml.neg = tmp.neg;
                ml.priority = priority::MUL;
                if tmp.priority >= ml.priority {
                    parenthesize(&mut tmp, flags);
                }
                // :1430: giac puts the presentation of the factor in the
                // content (empty in a content-only run).
                if flags.cont() {
                    let times = format!("<times/><imaginaryi/>{}", tmp.markup);
                    ml.content = mml::tag("apply", &times, self.next_id(), &[]);
                }
                if flags.pres() {
                    let i = format!("{}{}{}", tmp.markup, mml::ITIMES, mml::I);
                    ml.markup = mml::tag("mrow", &i, self.id(), &[]);
                }
                if flags.tex() {
                    ml.latex = format!("{} \\mathrm{{i}}", tmp.latex);
                }
            }
            return Ok(ml);
        }
        ml.priority = priority::ADD;
        let mut left = self.minus_prepended(re, flags)?;
        if left.priority >= ml.priority {
            parenthesize(&mut left, flags);
        }
        let mut right = self.gen2markup(im, flags)?;
        if right.priority >= priority::MUL {
            parenthesize(&mut right, flags);
        }
        if flags.cont() {
            let sign = if right.neg { "<minus/>" } else { "<plus/>" };
            let im = if isone {
                "<imaginaryi/>".to_owned()
            } else {
                let times = format!("<times/><imaginaryi/>{}", right.content);
                mml::tag("apply", &times, 0, &[])
            };
            let c = format!("{sign}{}{im}", left.content);
            ml.content = mml::tag("apply", &c, self.next_id(), &[]);
        }
        if flags.pres() {
            let sign = if right.neg { mml::MINUS } else { mml::PLUS };
            let im = if isone {
                String::new()
            } else {
                format!("{}{}", right.markup, mml::ITIMES)
            };
            let m = format!("{}{sign}{im}{}", left.markup, mml::I);
            ml.markup = mml::tag("mrow", &m, self.id(), &[]);
        }
        if flags.tex() {
            ml.latex = left.latex;
            ml.latex.push(if right.neg { '-' } else { '+' });
            if !isone {
                ml.latex.push_str(&right.latex);
                ml.latex.push(' ');
            }
            ml.latex.push_str("\\mathrm{i}");
        }
        Ok(ml)
    }

    /// :1658-1773: a matrix as an array, else the items joined by commas
    /// in the brackets of the subtype.
    fn vect(&self, items: &[Gen], sub: Sub, toplevel: bool, flags: Flags) -> Block {
        let mut ml = MarkupBlock::new();
        if sub == Sub::List && ckmatrix(items) {
            ml.kind = kind::MATRIX;
            let mut cols = 0;
            let mut rows = Vec::new();
            for r in items {
                let Gen::Vect(row, _) = r else {
                    unreachable!("checked by ckmatrix");
                };
                cols = row.len();
                let mut cells = Vec::new();
                ml.content.push_str("<matrixrow>");
                ml.markup.push_str("<mtr>");
                for x in row {
                    let tmp = self.minus_prepended(x, flags)?;
                    ml.content.push_str(&tmp.content);
                    if flags.pres() {
                        ml.markup.push_str(&format!("<mtd>{}</mtd>", tmp.markup));
                    }
                    cells.push(tmp.latex);
                }
                ml.content.push_str("</matrixrow>");
                ml.markup.push_str("</mtr>");
                rows.push(cells.join("&"));
                self.check(ml.len() + rows.last().expect("pushed").len())?;
            }
            if flags.cont() {
                ml.content = mml::tag("matrix", &ml.content, self.next_id(), &[]);
            } else {
                ml.content.clear();
            }
            if flags.pres() {
                let table = format!("<mtable>{}</mtable>", ml.markup);
                let fences = [("open", "["), ("close", "]")];
                ml.markup = mml::tag("mfenced", &table, self.id(), &fences);
            } else {
                ml.markup.clear();
            }
            if flags.tex() {
                ml.latex = format!(
                    "\\left[\\begin{{array}}{{{}}}{}\\end{{array}}\\right]",
                    "c".repeat(cols),
                    rows.join("\\\\")
                );
            }
            return Ok(ml);
        }
        if sub == Sub::Set && items.is_empty() {
            if flags.cont() {
                ml.content = mml::insert_id("<emptyset/>", self.next_id(), true);
            }
            if flags.pres() {
                ml.markup = mml::tag("mi", "&empty;", self.id(), &[]);
            }
            if flags.tex() {
                ml.latex.push_str("\\emptyset ");
            }
            return Ok(ml);
        }
        let (ld, rd) = match sub {
            Sub::Seq if toplevel => ("", ""),
            Sub::Seq => ("\\left(", "\\right)"),
            Sub::Set => ("\\left\\{", "\\right\\}"),
            Sub::List => ("\\left[", "\\right]"),
        };
        ml.latex.push_str(ld);
        for (k, x) in items.iter().enumerate() {
            if k > 0 {
                ml.latex.push(',');
            }
            ml.push(&self.minus_prepended(x, flags)?);
            self.check(ml.len())?;
        }
        ml.latex.push_str(rd);
        if !flags.tex() {
            ml.latex.clear();
        }
        // :1720-1762: the items run together in a fence of the subtype's
        // brackets (none for a top-level sequence).
        let (element, open, close) = match sub {
            Sub::Seq if toplevel => ("list", "", ""),
            Sub::Seq => ("list", "(", ")"),
            Sub::Set => ("set", "{", "}"),
            Sub::List => ("vector", "[", "]"),
        };
        if flags.cont() {
            ml.content = mml::tag(element, &ml.content, self.next_id(), &[]);
        }
        if flags.pres() {
            let fences = [("open", open), ("close", close)];
            ml.markup = mml::tag("mfenced", &ml.markup, self.id(), &fences);
        }
        Ok(ml)
    }

    /// `_SYMB` (:1774-3611), operator by operator in giac's order.
    #[allow(clippy::too_many_lines)]
    fn symb(&self, g: &Gen, flags: Flags, isfactor: bool) -> Block {
        let Gen::Symb(op, arg) = g else {
            unreachable!("a symbolic");
        };
        let op = op.as_str();
        let arg = &**arg;
        let args = g.args();
        let vectarg = matches!(arg, Gen::Vect(..));
        let isbinary = vectarg && args.len() == 2;
        let mut ml = MarkupBlock::new();
        match op {
            // :1777: the minus stays pending.
            "neg" => {
                let mut ml = self.gen2markup(arg, flags)?;
                ml.neg = !ml.neg;
                return Ok(ml);
            }
            // :1793-1847.
            "plus" => {
                if vectarg {
                    ml.priority = priority::ADD;
                    for (k, x) in flatten_operands(g).iter().enumerate() {
                        let mut tmp = self.gen2markup(x, flags)?;
                        if !tmp.neg && k > 0 {
                            if flags.pres() {
                                ml.markup.push_str(mml::PLUS);
                            }
                            if flags.tex() {
                                ml.latex.push('+');
                            }
                        }
                        prepend_minus_term(&mut tmp, flags);
                        if tmp.priority >= ml.priority {
                            parenthesize(&mut tmp, flags);
                        }
                        ml.push(&tmp);
                        self.check(ml.len())?;
                    }
                    if flags.cont() {
                        let c = format!("<plus/>{}", ml.content);
                        ml.content = mml::tag("apply", &c, self.next_id(), &[]);
                    }
                    if flags.pres() {
                        ml.markup = mml::tag("mrow", &ml.markup, self.id(), &[]);
                    }
                } else {
                    // :1827: unary plus, the id moved onto the new tags.
                    ml = self.gen2markup(arg, flags)?;
                    if !ml.neg && ml.priority == 0 {
                        if flags.cont() {
                            let id = mml::extract_id(&mut ml.content, true);
                            let c = format!("<plus/>{}", ml.content);
                            ml.content = mml::tag("apply", &c, 0, &[("id", &id)]);
                        }
                        if flags.pres() {
                            let id = mml::extract_id(&mut ml.markup, false);
                            let m = format!("{}{}", mml::PLUS, ml.markup);
                            ml.markup = mml::tag("mrow", &m, 0, &[("xref", &id)]);
                        }
                        if flags.tex() {
                            ml.latex.insert(0, '+');
                        }
                        ml.priority = priority::ADD;
                    }
                }
                return Ok(ml);
            }
            // :1848.
            "pointminus" if isbinary => {
                ml.priority = priority::ADD;
                let (mut left, mut right) = self.leftright(args, None, flags)?;
                if left.priority >= ml.priority {
                    parenthesize(&mut left, flags);
                }
                if flags.cont() {
                    let op = if right.neg { ".+" } else { ".-" };
                    let c = format!("<ci>{op}</ci>{}{}", left.content, right.content);
                    ml.content = mml::tag("apply", &c, self.next_id(), &[]);
                }
                if right.priority >= ml.priority {
                    parenthesize(&mut right, flags);
                }
                if flags.pres() {
                    let sign = if right.neg { "&oplus;" } else { "&ominus;" };
                    let m = format!("{}<mo>{sign}</mo>{}", left.markup, right.markup);
                    ml.markup = mml::tag("mrow", &m, self.id(), &[]);
                }
                if flags.tex() {
                    let sign = if right.neg { "\\oplus " } else { "\\ominus " };
                    ml.latex = format!("{}{sign}{}", left.latex, right.latex);
                }
                return Ok(ml);
            }
            // :1869-2032.
            "prod" | "ampersand_times" if vectarg => return self.prod(g, flags, ml),
            // :2033.
            "pointprod" if vectarg => {
                ml.priority = priority::MUL;
                for (k, x) in flatten_operands(g).iter().enumerate() {
                    let mut tmp = self.minus_prepended(x, flags)?;
                    if tmp.priority >= ml.priority {
                        parenthesize(&mut tmp, flags);
                    }
                    if k > 0 {
                        if flags.pres() {
                            ml.markup.push_str("<mo>&odot;</mo>");
                        }
                        if flags.tex() {
                            ml.latex.push_str("\\odot ");
                        }
                    }
                    ml.push(&tmp);
                    self.check(ml.len())?;
                }
                if flags.cont() {
                    let c = format!("<ci>.*</ci>{}", ml.content);
                    ml.content = mml::tag("apply", &c, self.next_id(), &[]);
                }
                if flags.pres() {
                    ml.markup = mml::tag("mrow", &ml.markup, self.id(), &[]);
                }
                return Ok(ml);
            }
            // :2082: `\frac1{X}` for any X.
            "inv" => {
                let tmp = self.gen2markup(arg, flags)?;
                ml.neg = tmp.neg;
                ml.kind = kind::FRACTION;
                ml.priority = priority::MUL;
                let one = "<cn type='integer'>1</cn>";
                self.wrap(&mut ml, flags, "<divide/>", one, &tmp.content);
                if flags.pres() {
                    let m = format!("<mn>1</mn>{}", tmp.markup);
                    ml.markup = mml::tag("mfrac", &m, self.id(), &[]);
                }
                if flags.tex() {
                    ml.latex = format!("\\frac1{{{}}}", tmp.latex);
                }
                return Ok(ml);
            }
            // :2097.
            "division" | "rdiv" if isbinary => {
                ml.kind = kind::FRACTION;
                ml.priority = priority::MUL;
                let num = self.gen2markup(&args[0], flags)?;
                let den = self.gen2markup(&args[1], flags)?;
                ml.neg = num.neg != den.neg;
                self.wrap(&mut ml, flags, "<divide/>", &num.content, &den.content);
                if flags.pres() {
                    let m = format!("{}{}", num.markup, den.markup);
                    ml.markup = mml::tag("mfrac", &m, self.id(), &[]);
                }
                if flags.tex() {
                    ml.latex = format!("\\frac{{{}}}{{{}}}", num.latex, den.latex);
                }
                return Ok(ml);
            }
            // :2130.
            "conj" => {
                let mut tmp = self.gen2markup(arg, flags)?;
                if isfactor && tmp.priority > priority::MUL {
                    parenthesize(&mut tmp, flags);
                }
                if tmp.ctype(kind::LEADING_DIGIT) {
                    ml.kind = kind::LEADING_DIGIT;
                }
                ml.priority = priority::UNARY;
                ml.neg = tmp.neg;
                self.wrap(&mut ml, flags, "<conjugate/>", &tmp.content, "");
                if flags.pres() {
                    let m = format!("{}<mo>&#xaf;</mo>", tmp.markup);
                    ml.markup = mml::tag("mover", &m, self.id(), &[]);
                }
                if flags.tex() {
                    ml.latex = format!("\\overline{{{}}}", tmp.latex);
                }
                return Ok(ml);
            }
            // :2148.
            "exp" | "EXP" => {
                if arg.is_one() {
                    if flags.cont() {
                        ml.content = mml::insert_id("<exponentiale/>", self.next_id(), true);
                    }
                    if flags.pres() {
                        ml.markup = mml::insert_id(mml::E, self.id(), false);
                    }
                    if flags.tex() {
                        ml.latex.push_str("\\mathrm{e}");
                    }
                } else {
                    let tmp = self.minus_prepended(arg, flags)?;
                    ml.priority = priority::EXP;
                    ml.kind = kind::POWER;
                    self.wrap(&mut ml, flags, "<exp/>", &tmp.content, "");
                    if flags.pres() {
                        let m = format!("{}{}", mml::E, tmp.markup);
                        ml.markup = mml::tag("msup", &m, self.id(), &[]);
                    }
                    if flags.tex() {
                        ml.latex = format!("\\mathrm{{e}}^{{{}}}", tmp.latex);
                    }
                }
                return Ok(ml);
            }
            // :2191.
            "sqrt" => {
                ml.priority = priority::UNARY;
                ml.kind = kind::ROOT;
                let tmp = self.minus_prepended(arg, flags)?;
                self.root(&mut ml, flags, &tmp);
                return Ok(ml);
            }
            // :2224-2280.
            "pow" | "matpow" | "pointpow" if isbinary => return self.pow(op, args, flags, ml),
            // :2281: abs drops the child's pending minus (giac quirk).
            "abs" | "floor" | "ceil" => {
                let mut tmp = self.gen2markup(arg, flags)?;
                if op != "abs" {
                    prepend_minus(&mut tmp, flags);
                }
                let (element, ld_mml, rd_mml, ld, rd) = match op {
                    "abs" => ("<abs/>", "|", "|", "|", "|"),
                    "floor" => (
                        "<floor/>",
                        "&lfloor;",
                        "&rfloor;",
                        "\\lfloor ",
                        "\\rfloor ",
                    ),
                    _ => ("<ceiling/>", "&lceil;", "&rceil;", "\\lceil ", "\\rceil "),
                };
                self.wrap(&mut ml, flags, element, &tmp.content, "");
                if flags.pres() {
                    let fences = [("open", ld_mml), ("close", rd_mml)];
                    ml.markup = mml::tag("mfenced", &tmp.markup, self.id(), &fences);
                }
                if flags.tex() {
                    ml.latex = format!("\\left{ld}{}\\right{rd}", tmp.latex);
                }
                return Ok(ml);
            }
            // :2375-2487 (an empty sequence, which giac never has, goes to
            // the general arm).
            "integrate" | "int" | "Int" | "sum" | "somme" | "add" | "product" | "mul"
                if !(vectarg && args.is_empty()) =>
            {
                return self.bigop(g, flags, ml);
            }
            // :2488-2548.
            "limit" | "limite" if vectarg && args.len() > 1 => return self.limit(args, flags, ml),
            // :2560-2727; giac falls through when it is no derivative it
            // writes, with the fields set so far.
            "derive" | "deriver" if self.derive(g, flags, &mut ml)? => return Ok(ml),
            _ => {}
        }
        match op {
            // :2728: sto's arguments are (value, variable).
            "same" | "different" | "sto" | "array_sto" | "equal" if isbinary => {
                let sto = matches!(op, "sto" | "array_sto");
                ml.priority = if sto {
                    priority::ASSGN
                } else if op == "equal" {
                    priority::EQ
                } else {
                    priority::COMP
                };
                let (mut left, mut right) = self.leftright(args, Some(ml.priority), flags)?;
                if sto {
                    std::mem::swap(&mut left, &mut right);
                }
                let (csymb, opm, opt) = match op {
                    "same" => ("<ci>==</ci>", "&#x2a75;", "=="),
                    "different" => ("<neq/>", "&#x2260;", "\\neq "),
                    "equal" => ("<eq/>", "=", "="),
                    "sto" => ("<ci>:=</ci>", "&#x2254;", ":="),
                    _ => ("<ci>=&lt;</ci>", "&larr;", "\\leftarrow"),
                };
                self.infix(&mut ml, flags, &left, &right, csymb, opm, opt);
                return Ok(ml);
            }
            // :2770.
            "inferieur_strict" | "inferieur_egal" | "superieur_strict" | "superieur_egal"
                if isbinary =>
            {
                ml.priority = priority::INEQ;
                let (left, right) = self.leftright(args, Some(ml.priority), flags)?;
                let (csymb, opm, opt) = match op {
                    "inferieur_strict" => ("<lt/>", "&lt;", "<"),
                    "inferieur_egal" => ("<leq/>", "&leq;", "\\leq "),
                    "superieur_strict" => ("<gt/>", "&gt;", ">"),
                    _ => ("<geq/>", "&geq;", "\\geq "),
                };
                self.infix(&mut ml, flags, &left, &right, csymb, opm, opt);
                return Ok(ml);
            }
            // :2942: content `<interval>`, not an `<apply>`.
            "interval" if isbinary => {
                ml.priority = priority::SET;
                let (left, right) = self.leftright(args, Some(ml.priority), flags)?;
                if flags.cont() {
                    let c = format!("{}{}", left.content, right.content);
                    ml.content = mml::tag("interval", &c, self.next_id(), &[]);
                }
                if flags.pres() {
                    let m = format!("{}<mo>&#x2025;</mo>{}", left.markup, right.markup);
                    ml.markup = mml::tag("mrow", &m, self.id(), &[]);
                }
                if flags.tex() {
                    ml.latex = format!("{}\\ldots {}", left.latex, right.latex);
                }
                return Ok(ml);
            }
            // :3317.
            "factorial" => {
                ml.priority = priority::UNARY;
                ml.kind = kind::FACTORIAL;
                let mut tmp = self.minus_prepended(arg, flags)?;
                if tmp.priority >= ml.priority {
                    parenthesize(&mut tmp, flags);
                }
                self.wrap(&mut ml, flags, "<factorial/>", &tmp.content, "");
                if flags.pres() {
                    let m = format!("{}<mo>!</mo>", tmp.markup);
                    ml.markup = mml::tag("mrow", &m, self.id(), &[]);
                }
                if flags.tex() {
                    ml.latex = format!("{}\\mathclose{{!}}", tmp.latex);
                }
                return Ok(ml);
            }
            // :3334.
            "binomial" | "nCr" | "comb" if isbinary => {
                let left = self.minus_prepended(&args[0], flags)?;
                let right = self.minus_prepended(&args[1], flags)?;
                let head = mml::csymbol("binomial", "combinat1");
                self.wrap(&mut ml, flags, &head, &left.content, &right.content);
                if flags.pres() {
                    let m = format!(
                        "<mfrac linethickness='0'>{}{}</mfrac>",
                        left.markup, right.markup
                    );
                    ml.markup = mml::tag("mfenced", &m, self.id(), &[]);
                }
                if flags.tex() {
                    ml.latex = format!("\\binom{{{}}}{{{}}}", left.latex, right.latex);
                }
                return Ok(ml);
            }
            // :3434.
            "at" if isbinary => return self.at(args, flags, ml),
            _ => {}
        }
        self.general(g, flags, ml)
    }

    /// :1039 `get_leftright`: both sides, minus written, wrapped when they
    /// bind as weakly as `priority` or weaker.
    fn leftright(
        &self,
        args: &[Gen],
        priority: Option<u8>,
        flags: Flags,
    ) -> Result<(MarkupBlock, MarkupBlock), TooLong> {
        let mut left = self.minus_prepended(&args[0], flags)?;
        if priority.is_some_and(|p| left.priority >= p) {
            parenthesize(&mut left, flags);
        }
        let mut right = self.minus_prepended(&args[args.len() - 1], flags)?;
        if priority.is_some_and(|p| right.priority >= p) {
            parenthesize(&mut right, flags);
        }
        Ok((left, right))
    }

    /// The relations' texts (:2760-2767, :2797-2804): content
    /// `<apply>csymb left right</apply>`, presentation
    /// `<mrow>left<mo>opm</mo>right</mrow>`, LaTeX `left opt right`.
    #[allow(clippy::too_many_arguments)]
    fn infix(
        &self,
        ml: &mut MarkupBlock,
        flags: Flags,
        left: &MarkupBlock,
        right: &MarkupBlock,
        csymb: &str,
        opm: &str,
        opt: &str,
    ) {
        self.wrap(ml, flags, csymb, &left.content, &right.content);
        if flags.pres() {
            let m = format!("{}<mo>{opm}</mo>{}", left.markup, right.markup);
            ml.markup = mml::tag("mrow", &m, self.id(), &[]);
        }
        if flags.tex() {
            ml.latex = format!("{}{opt}{}", left.latex, right.latex);
        }
    }

    /// :2224-2280: `\sqrt{x}` for an exponent 1/2, the power on the name of
    /// an elementary function (`\cos ^{2}x`, ELEMPOW), else `x^{e}`.
    fn pow(&self, op: &str, args: &[Gen], flags: Flags, mut ml: MarkupBlock) -> Block {
        ml.priority = priority::EXP;
        let mut left = self.minus_prepended(&args[0], flags)?;
        let e = &args[1];
        // is_one(2*exponent) :2230: a number, not a symbolic 1/2.
        if op == "pow" && e.equals(1, 2) {
            ml.kind = kind::ROOT;
            ml.priority = priority::UNARY;
            self.root(&mut ml, flags, &left);
            return Ok(ml);
        }
        let right = self.minus_prepended(e, flags)?;
        let head = if op == "pointpow" {
            "<ci>.^</ci>"
        } else {
            "<power/>"
        };
        self.wrap(&mut ml, flags, head, &left.content, &right.content);
        // :2248.
        if flags.has(Flags::ELEMPOW)
            && op == "pow"
            && right.priority == 0
            && !(e.is_one() || e.is_minus_one())
            && left.ctype(kind::ELEMAPP)
            && left.split_pos >= 0
        {
            ml.kind = kind::ELEMAPP;
            ml.appl = left.appl;
            // :2254-2257: `<msup>` around the name, after the outer tag:
            // `<mrow><msup><mi>cos</mi><mn>2</mn></msup><mo>&af;</mo>…`.
            if flags.pres() {
                let pos = left.markup.find('>').expect("a tag") + 1;
                let end = pos + usize::try_from(left.split_pos).expect("not negative");
                ml.markup = format!(
                    "{}<msup>{}{}</msup>{}",
                    &left.markup[..pos],
                    &left.markup[pos..end],
                    right.markup,
                    &left.markup[end..]
                );
            }
            if flags.tex() {
                let k = left.split_pos_tex;
                ml.latex = format!(
                    "{}^{{{}}}{}",
                    &left.latex[..k],
                    right.latex,
                    &left.latex[k..]
                );
            }
        } else {
            ml.kind = kind::POWER;
            if left.ctype(kind::LEADING_DIGIT) {
                ml.kind |= kind::LEADING_DIGIT;
            }
            if left.appl || left.priority >= ml.priority {
                parenthesize(&mut left, flags);
            }
            if flags.pres() {
                let m = format!("{}{}", left.markup, right.markup);
                ml.markup = mml::tag("msup", &m, self.id(), &[]);
            }
            if flags.tex() {
                ml.latex = format!("{}^{{{}}}", left.latex, right.latex);
            }
        }
        Ok(ml)
    }

    /// :1869-2032: a product, the numbers and constants first, a leading
    /// integer fraction gathered, factors `inv(…)` in a denominator,
    /// `\cdot ` or a space between factors.
    #[allow(clippy::too_many_lines)]
    fn prod(&self, g: &Gen, flags: Flags, mut ml: MarkupBlock) -> Block {
        let mut rest = flatten_operands(g);
        let (mut ncnst, mut cnsts) = (Vec::new(), Vec::new());
        // From the end: integers and fractions, then pi and euler_gamma, so
        // both come out reversed; a pure imaginary c (im != 1) gives its
        // imaginary part to the constants and is replaced by i in place.
        for i in (0..rest.len()).rev() {
            match &rest[i] {
                Gen::Int(_) | Gen::Frac(..) => ncnst.push(rest.remove(i)),
                Gen::Idnt(n) if n == "pi" || n == "euler_gamma" => cnsts.push(rest.remove(i)),
                Gen::Cplx(re, im) if re.is_zero() && !im.is_one() => {
                    if im.is_minus_one() {
                        ml.neg = !ml.neg;
                    } else {
                        let im = (**im).clone();
                        cnsts.push(im);
                    }
                    rest[i] = Gen::i();
                }
                _ => {}
            }
        }
        let mut args: Vec<Gen> = ncnst.into_iter().chain(cnsts).chain(rest).collect();
        ml.priority = priority::MUL;
        // :1897-1937: a leading fraction when the one inv(…) has an integer
        // argument and the numerator has several factors.
        let (mut num_count, mut den_count) = (0usize, 0usize);
        let mut hasleadingfrac = true;
        let (mut np, mut dp): (Option<usize>, Option<usize>) = (None, None);
        for (k, x) in args.iter().enumerate() {
            if inv_arg(x).is_some() {
                den_count += 1;
                if !hasleadingfrac {
                    continue;
                }
                if matches!(inv_arg(x), Some(Gen::Int(_))) && dp.is_none() {
                    dp = Some(k);
                } else {
                    dp = None;
                    hasleadingfrac = false;
                }
            } else {
                num_count += 1;
                if hasleadingfrac && np.is_none() && matches!(x, Gen::Int(_)) {
                    np = Some(k);
                }
            }
        }
        if num_count == 1 {
            hasleadingfrac = false;
        }
        if let (true, Some(dp)) = (hasleadingfrac, dp) {
            let di = inv_arg(&args.remove(dp)).expect("an inv").clone();
            let np = np.map(|np| if np > dp { np - 1 } else { np });
            // giac's fraction(n, d): not reduced.
            let ni = match np {
                Some(np) => args.remove(np),
                None => Gen::int(1),
            };
            args.insert(0, Gen::Frac(Box::new(ni), Box::new(di)));
            num_count = args.iter().filter(|x| inv_arg(x).is_none()).count();
            den_count = args.len() - num_count;
        }
        let (mut pnum, mut pden) = (MarkupBlock::new(), MarkupBlock::new());
        // The numerator's and denominator's texts so far (`latex`: LaTeX).
        let (mut num, mut den) = (MarkupBlock::new(), MarkupBlock::new());
        let (mut nc, mut dc, mut neg_count) = (0usize, 0usize, 0usize);
        for x in &args {
            let isinv = inv_arg(x).is_some();
            let several = if isinv { den_count > 1 } else { num_count > 1 };
            let flg = if several {
                flags | Flags::FACTOR
            } else {
                flags
            };
            let fc = inv_arg(x).unwrap_or(x);
            let mut tmp = self.gen2markup(fc, flg)?;
            if tmp.neg {
                neg_count += 1;
            }
            if isinv {
                dc += 1;
            } else {
                nc += 1;
            }
            let not_last = if isinv {
                dc < den_count
            } else {
                nc < num_count
            };
            if several
                && (tmp.priority > ml.priority
                    || (tmp.priority == ml.priority && !matches!(fc, Gen::Cplx(..))))
                && !tmp.ctype(kind::FRACTION)
                && (!tmp.ctype(kind::SUMPROD) || not_last)
            {
                parenthesize(&mut tmp, flags);
            }
            let prev = if isinv { &pden } else { &pnum };
            let is_cdot = tmp.ctype(kind::LEADING_DIGIT)
                || tmp.ctype(kind::SUBTYPE_IDNT)
                || tmp.ctype(kind::FACTORIAL)
                || tmp.ctype(kind::IDNT_NAME)
                || prev.ctype(kind::IDNT_NAME)
                || prev.ctype(kind::SUBTYPE_IDNT)
                || (tmp.ctype(kind::FRACTION) && prev.ctype(kind::FRACTION))
                || (prev.ctype(kind::ELEMAPP) && prev.appl);
            let side = if isinv { &mut den } else { &mut num };
            if flags.pres() && !side.markup.is_empty() {
                side.markup
                    .push_str(if is_cdot { mml::CDOT } else { mml::ITIMES });
            }
            if flags.tex() && !side.latex.is_empty() {
                side.latex.push_str(if is_cdot { "\\cdot " } else { " " });
            }
            side.push(&tmp);
            self.check(num.len() + den.len())?;
            if isinv {
                pden = tmp;
            } else {
                pnum = tmp;
            }
        }
        if neg_count % 2 == 1 {
            ml.neg = !ml.neg;
        }
        if num.content.is_empty() {
            num.content.push_str("<cn type='integer'>1</cn>");
        }
        if num.markup.is_empty() {
            num.markup.push_str("<mn>1</mn>");
        }
        if num.latex.is_empty() {
            num.latex.push('1');
        }
        // :1993-2018: `<times/>` or `<divide/>` over all the factors (giac
        // puts every factor of the numerator and denominator as an
        // argument of `<divide/>`); `<mfrac>` with an `<mrow>` for a side
        // of several factors.
        if flags.cont() {
            let head = if den.content.is_empty() {
                "<times/>"
            } else {
                "<divide/>"
            };
            self.wrap(&mut ml, flags, head, &num.content, &den.content);
        }
        if flags.pres() {
            ml.markup = if den.markup.is_empty() {
                mml::tag("mrow", &num.markup, self.id(), &[])
            } else {
                let num = if nc > 1 {
                    mml::tag("mrow", &num.markup, 0, &[])
                } else {
                    num.markup
                };
                let den = if dc > 1 {
                    mml::tag("mrow", &den.markup, 0, &[])
                } else {
                    den.markup.clone()
                };
                mml::tag("mfrac", &(num + &den), self.id(), &[])
            };
        }
        if flags.tex() {
            ml.latex = if den.latex.is_empty() {
                num.latex
            } else {
                format!("\\frac{{{}}}{{{}}}", num.latex, den.latex)
            };
        }
        if !(den.latex.is_empty() && den.markup.is_empty() && den.content.is_empty()) {
            ml.kind = kind::FRACTION;
        }
        Ok(ml)
    }

    /// :2375-2487: `\int _{a}^{b}f\,\mathrm{d}x`, `\sum _{k=1}^{n}…`,
    /// `\prod _{k}…`; content `<apply><int/><bvar>…</bvar><lowlimit>…`,
    /// presentation `<munderover><mo>&int;</mo>…`.
    #[allow(clippy::too_many_lines)]
    fn bigop(&self, g: &Gen, flags: Flags, mut ml: MarkupBlock) -> Block {
        let Gen::Symb(op, arg) = g else {
            unreachable!("a symbolic");
        };
        let isint = matches!(op.as_str(), "integrate" | "int" | "Int");
        let isprod = matches!(op.as_str(), "product" | "mul");
        let (mut big, mut big_tex, head) = if isint {
            ("<mo>&int;</mo>".to_owned(), "\\int ".to_owned(), "<int/>")
        } else if isprod {
            ("<mo>&prod;</mo>".to_owned(), "\\prod ".to_owned(), "<product/>")
        } else {
            ("<mo>&sum;</mo>".to_owned(), "\\sum ".to_owned(), "<sum/>")
        };
        ml.priority = priority::MUL;
        ml.kind = kind::SUMPROD;
        let mut var = MarkupBlock::new();
        let (mut lb, mut ub) = (None, None);
        let mut e;
        if let Gen::Vect(a, _) = &**arg {
            e = self.gen2markup(&a[0], flags)?;
            if a.len() > 1 {
                let bounds = match &a[1] {
                    eq @ Gen::Symb(o, _) if o == "equal" => match eq.args() {
                        [lh, iv] if iv.is_sym(&["interval"]) => Some((lh, iv.args())),
                        _ => None,
                    },
                    _ => None,
                };
                if let Some((lh, [low, .., up])) = bounds {
                    lb = Some(self.gen2markup(low, flags)?);
                    ub = Some(self.gen2markup(up, flags)?);
                    var = self.gen2markup(lh, flags)?;
                } else {
                    var = self.gen2markup(&a[1], flags)?;
                }
                prepend_minus(&mut var, flags);
                if isint {
                    if var.priority != 0 {
                        parenthesize(&mut var, flags);
                    }
                    if flags.pres() {
                        var.markup = format!("{}{}{}", mml::ITIMES, mml::D, var.markup);
                    }
                    if flags.tex() {
                        var.latex = format!("\\,\\mathrm{{d}}{}", var.latex);
                    }
                }
            }
            if lb.is_none() && a.len() > 2 {
                lb = Some(self.gen2markup(&a[2], flags)?);
            }
            if ub.is_none() && a.len() > 3 {
                ub = Some(self.gen2markup(&a[3], flags)?);
            }
            if let Some(lb) = &mut lb {
                prepend_minus(lb, flags);
            }
            if let Some(ub) = &mut ub {
                prepend_minus(ub, flags);
            }
            if !isint {
                lb = Some(match lb {
                    Some(mut lb) => {
                        if flags.pres() {
                            let m = format!("{}<mo>=</mo>{}", var.markup, lb.markup);
                            lb.markup = mml::tag("mrow", &m, 0, &[]);
                        }
                        if flags.tex() {
                            lb.latex = format!("{}={}", var.latex, lb.latex);
                        }
                        lb
                    }
                    None => var.clone(),
                });
            }
            match (&lb, &ub) {
                (Some(lb), None) => {
                    big = mml::tag("munder", &(big + &lb.markup), 0, &[]);
                    big_tex = format!("{big_tex}_{{{}}}", lb.latex);
                }
                (Some(lb), Some(ub)) => {
                    let m = format!("{big}{}{}", lb.markup, ub.markup);
                    big = mml::tag("munderover", &m, 0, &[]);
                    big_tex = format!("{big_tex}_{{{}}}^{{{}}}", lb.latex, ub.latex);
                }
                _ => {}
            }
        } else {
            e = self.gen2markup(arg, flags)?;
        }
        if isprod {
            prepend_minus(&mut e, flags);
        } else {
            ml.neg = e.neg;
        }
        if (isprod && !e.ctype(kind::FRACTION) && e.priority >= ml.priority)
            || (!isprod && (e.ctype(kind::NEG) || e.priority > ml.priority))
        {
            parenthesize(&mut e, flags);
        }
        if flags.cont() {
            let mut s = String::from(head);
            if !var.content.is_empty() {
                s.push_str(&mml::tag("bvar", &var.content, 0, &[]));
            }
            match (&lb, &ub) {
                (Some(lb), None) => {
                    let t = if isint {
                        "domainofapplication"
                    } else {
                        "condition"
                    };
                    s.push_str(&mml::tag(t, &lb.content, 0, &[]));
                }
                (Some(lb), Some(ub)) => {
                    s.push_str(&mml::tag("lowlimit", &lb.content, 0, &[]));
                    s.push_str(&mml::tag("uplimit", &ub.content, 0, &[]));
                }
                _ => {}
            }
            s.push_str(&e.content);
            ml.content = mml::tag("apply", &s, self.next_id(), &[]);
        }
        if flags.pres() {
            let v = if isint { var.markup.as_str() } else { "" };
            let m = format!("{big}{}{v}", e.markup);
            ml.markup = mml::tag("mrow", &m, self.id(), &[]);
        }
        if flags.tex() {
            ml.latex = big_tex + &e.latex;
            if isint {
                ml.latex.push_str(&var.latex);
            }
        }
        Ok(ml)
    }

    /// :2488-2548: `\lim_{x\to a}`, `^{+}`/`^{-}` for a direction ±1;
    /// content `<apply><limit/><bvar>…</bvar><condition><tendsto…>`,
    /// presentation `<munder><mi>lim</mi>…`.
    fn limit(&self, a: &[Gen], flags: Flags, mut ml: MarkupBlock) -> Block {
        ml.priority = priority::MUL;
        ml.kind = kind::SUMPROD;
        let mut dirn = 0i64;
        let mut e = self.gen2markup(&a[0], flags)?;
        ml.neg = e.neg;
        if e.priority > ml.priority {
            parenthesize(&mut e, flags);
        } else if flags.pres() {
            e.markup.insert_str(0, mml::APPLY);
        }
        let (mut var, mut dest) = (MarkupBlock::new(), MarkupBlock::new());
        let direction = |g: &Gen| match g {
            Gen::Int(n) => n.to_i64().unwrap_or(0),
            _ => 0,
        };
        if a[1].is_sym(&["equal"]) {
            let eq = a[1].args();
            var = self.gen2markup(&eq[0], flags)?;
            dest = self.gen2markup(&eq[eq.len() - 1], flags)?;
            if a.len() == 3 {
                dirn = direction(&a[2]);
            }
        } else if a.len() >= 3 {
            var = self.gen2markup(&a[1], flags)?;
            dest = self.gen2markup(&a[2], flags)?;
            if a.len() == 4 {
                dirn = direction(&a[3]);
            }
        }
        prepend_minus(&mut var, flags);
        prepend_minus(&mut dest, flags);
        if flags.cont() {
            let mut s = String::from("<limit/>");
            if !var.content.is_empty() {
                s.push_str(&mml::tag("bvar", &var.content, 0, &[]));
                if !dest.content.is_empty() {
                    let (k, v) = match dirn.signum() {
                        0 => ("", ""),
                        -1 => ("type", "below"),
                        _ => ("type", "above"),
                    };
                    let c = format!("{}{}", var.content, dest.content);
                    let t = mml::tag("tendsto", &c, 0, &[(k, v)]);
                    s.push_str(&mml::tag("condition", &t, 0, &[]));
                }
            }
            s.push_str(&e.content);
            ml.content = mml::tag("apply", &s, self.next_id(), &[]);
        }
        if flags.pres() {
            let to = match dirn.signum() {
                0 => dest.markup.clone(),
                -1 => format!("<msup>{}{}</msup>", dest.markup, mml::MINUS),
                _ => format!("<msup>{}{}</msup>", dest.markup, mml::PLUS),
            };
            let m = format!(
                "<munder><mi>lim</mi><mrow>{}<mo>&rarr;</mo>{to}</mrow></munder>{}",
                var.markup, e.markup
            );
            ml.markup = mml::tag("mrow", &m, self.id(), &[]);
        }
        if flags.tex() {
            let side = match dirn.signum() {
                0 => "",
                -1 => "^{-}",
                _ => "^{+}",
            };
            ml.latex = format!(
                "\\lim_{{{}\\to {}{side}}}{}",
                var.latex, dest.latex, e.latex
            );
        }
        Ok(ml)
    }

    /// :2560-2727: partial derivatives of a function `f^{\left(1,2\right)}`,
    /// `f'`, `f''`, `f'''`, `f^{(4)}`, and `derive(e, vars…)` as a quotient
    /// of differentials (content `<diff/>`, `<partialdiff/>`, `nthdiff`).
    /// `false` when giac falls through to the general arm, with `ml`'s
    /// fields set so far.
    #[allow(clippy::too_many_lines)]
    fn derive(&self, g: &Gen, flags: Flags, ml: &mut MarkupBlock) -> Result<bool, TooLong> {
        if let Some((f, deg)) = is_partialdiff(g) {
            ml.priority = priority::EXP;
            ml.kind = kind::POWER | kind::DERIVATIVE;
            let left = self.gen2markup(f, flags)?;
            let right = self.gen2markup(&Gen::Vect(deg, Sub::Seq), flags)?;
            if flags.tex() {
                ml.latex = format!("{}^{{{}}}", left.latex, right.latex);
            }
            let head = "<partialdiff/>";
            self.wrap(ml, flags, head, &right.content, &left.content);
            if flags.pres() {
                let m = format!("{}{}", left.markup, right.markup);
                ml.markup = mml::tag("msup", &m, self.id(), &[]);
            }
            return Ok(true);
        }
        let Gen::Symb(_, arg) = g else {
            unreachable!("a symbolic");
        };
        let mut arg = &**arg;
        let vectarg = matches!(arg, Gen::Vect(..));
        let mut n = 0u64;
        while let Gen::Symb(op, inner) = arg
            && matches!(op.as_str(), "derive" | "deriver")
        {
            arg = inner;
            n += 1;
        }
        if n > 0 && vectarg {
            return Ok(false);
        }
        ml.kind = kind::DERIVATIVE;
        let items: &[Gen] = match arg {
            Gen::Vect(items, _) if vectarg => items,
            _ => std::slice::from_ref(arg),
        };
        let Some(first) = items.first() else {
            return Ok(false);
        };
        let mut tmp = self.gen2markup(first, flags)?;
        let simp = tmp.priority == 0;
        if tmp.priority > priority::APPLY {
            parenthesize(&mut tmp, flags);
        }
        ml.neg = tmp.neg;
        if n > 0 {
            ml.priority = priority::EXP;
            ml.kind |= kind::POWER;
            let k = (n + 1).to_string();
            let head = format!("{}<cn>{k}</cn>", mml::csymbol("nthdiff", "calculus1"));
            self.wrap(ml, flags, &head, &tmp.content, "");
            if flags.pres() {
                let pm = match n + 1 {
                    2 => "<mo>&Prime;</mo>".to_owned(),
                    3 => "<mo>&tprime;</mo>".to_owned(),
                    _ => mml::tag("mfenced", &format!("<mn>{k}</mn>"), 0, &[]),
                };
                let m = format!("{}{pm}", tmp.markup);
                ml.markup = mml::tag("msup", &m, self.id(), &[]);
            }
            if flags.tex() {
                ml.latex = tmp.latex;
                match n + 1 {
                    2 => ml.latex.push_str("''"),
                    3 => ml.latex.push_str("'''"),
                    _ => {
                        ml.latex.push_str("^{(");
                        ml.latex.push_str(&k);
                        ml.latex.push_str(")}");
                    }
                }
            }
            return Ok(true);
        }
        if !vectarg {
            ml.priority = priority::EXP;
            self.wrap(ml, flags, "<diff/>", &tmp.content, "");
            if flags.pres() {
                let m = format!("{}<mo>&prime;</mo>", tmp.markup);
                ml.markup = mml::tag("msup", &m, self.id(), &[]);
            }
            if flags.tex() {
                ml.latex = tmp.latex + "'";
            }
            return Ok(true);
        }
        let mut vars = Vec::new();
        if items.len() > 1 && get_derive_vars(items, &mut vars) {
            ml.priority = priority::MUL;
            ml.kind |= kind::FRACTION;
            let d = mml::D;
            if let [(v, 1)] = vars.as_slice() {
                let var = self.gen2markup(v, flags)?;
                let bvar = format!("<diff/><bvar>{}</bvar>", var.content);
                self.wrap(ml, flags, &bvar, &tmp.content, "");
                if flags.pres() {
                    ml.markup = if simp {
                        let m = format!(
                            "<mrow>{d}{}</mrow><mrow>{d}{}</mrow>",
                            tmp.markup, var.markup
                        );
                        mml::tag("mfrac", &m, self.id(), &[])
                    } else {
                        let v = mml::tag("mrow", &var.markup, 0, &[]);
                        let f = mml::tag("mfrac", &format!("{d}{v}"), 0, &[]);
                        let m = format!("{f}{}{}", mml::APPLY, tmp.markup);
                        mml::tag("mrow", &m, self.id(), &[])
                    };
                }
                if flags.tex() {
                    ml.latex = if simp {
                        format!(
                            "\\frac{{\\mathrm{{d}}{}}}{{\\mathrm{{d}}{}}}",
                            tmp.latex, var.latex
                        )
                    } else {
                        format!(
                            "\\frac{{\\mathrm{{d}}}}{{\\mathrm{{d}}{}}}{}",
                            var.latex, tmp.latex
                        )
                    };
                }
                return Ok(true);
            }
            let nv = vars.len();
            let (ds, ds_tex) = if nv > 1 {
                ("<mo>&part;</mo>", "\\partial ")
            } else {
                (d, "\\mathrm{d}")
            };
            let mut parts = Vec::new();
            for (c, (v, k)) in vars.iter().enumerate() {
                let var = self.gen2markup(v, flags)?;
                if flags.cont() {
                    let degree = if *k == 1 {
                        String::new()
                    } else {
                        format!("<degree><cn>{k}</cn></degree>")
                    };
                    let b = format!("{}{degree}", var.content);
                    ml.content.push_str(&mml::tag("bvar", &b, self.id(), &[]));
                }
                if flags.pres() {
                    if *k == 1 {
                        ml.markup.push_str(ds);
                        ml.markup.push_str(&var.markup);
                    } else {
                        let m = format!("<mrow>{ds}{}</mrow><mn>{k}</mn>", var.markup);
                        ml.markup.push_str(&mml::tag("msup", &m, 0, &[]));
                    }
                    if c + 1 < nv {
                        ml.markup.push_str(mml::ITIMES);
                    }
                }
                let power = if *k == 1 {
                    String::new()
                } else {
                    format!("^{{{k}}}")
                };
                parts.push(format!("{ds_tex}{}{power}", var.latex));
                self.check(ml.len())?;
            }
            let body = parts.join(" ");
            let ns: u64 = vars.iter().map(|(_, k)| k).sum();
            if flags.cont() {
                let head = if nv > 1 { "<partialdiff/>" } else { "<diff/>" };
                let c = format!(
                    "{head}{}<degree><cn type='integer'>{ns}</cn></degree>{}",
                    ml.content, tmp.content
                );
                ml.content = mml::tag("apply", &c, self.next_id(), &[]);
            }
            if flags.pres() {
                let sup = format!("<msup>{ds}<mn>{ns}</mn></msup>");
                ml.markup = if simp {
                    let m = format!(
                        "<mrow>{sup}{}</mrow><mrow>{}</mrow>",
                        tmp.markup, ml.markup
                    );
                    mml::tag("mfrac", &m, self.id(), &[])
                } else {
                    let f = mml::tag("mfrac", &format!("{sup}<mrow>{}</mrow>", ml.markup), 0, &[]);
                    let m = format!("{f}{}{}", mml::APPLY, tmp.markup);
                    mml::tag("mrow", &m, self.id(), &[])
                };
            }
            if flags.tex() {
                ml.latex = if simp {
                    format!("\\frac{{{ds_tex}^{{{ns}}}{}}}{{{body}}}", tmp.latex)
                } else {
                    format!("\\frac{{{ds_tex}^{{{ns}}}}}{{{body}}}{}", tmp.latex)
                };
            }
            return Ok(true);
        }
        Ok(false)
    }

    /// :3434: `x_{1}`, `a_{i j}`, `a_{i,10}`; content `<selector/>`,
    /// presentation `<msub>` (indices joined by `&ic;` or `,`).
    fn at(&self, a: &[Gen], flags: Flags, mut ml: MarkupBlock) -> Block {
        ml.priority = priority::APPLY;
        ml.kind = kind::HAS_SUBSCRIPT;
        let mut tmp = self.minus_prepended(&a[0], flags)?;
        if tmp.priority > ml.priority || tmp.ctype(kind::HAS_SUBSCRIPT) {
            parenthesize(&mut tmp, flags);
        }
        let sub = match &a[1] {
            Gen::Vect(idx, Sub::Seq) => {
                let one_char = |x: &Gen| match x {
                    Gen::Int(n) => n.to_i64().is_some_and(|n| (1..=9).contains(&n)),
                    Gen::Idnt(name) => {
                        let mut chars = name.chars();
                        matches!((chars.next(), chars.next()), (Some(c), None) if my_isalpha(c))
                            || is_greek_letter(name)
                    }
                    _ => false,
                };
                let (sep, sept) = if idx.iter().all(one_char) {
                    (mml::ICOMMA, " ")
                } else {
                    ("<mo>,</mo>", ",")
                };
                let mut sub = MarkupBlock::new();
                let mut subs = Vec::new();
                for (k, x) in idx.iter().enumerate() {
                    let s = self.minus_prepended(x, flags)?;
                    sub.content.push_str(&s.content);
                    if k > 0 {
                        sub.markup.push_str(sep);
                    }
                    sub.markup.push_str(&s.markup);
                    subs.push(s.latex);
                    self.check(sub.len())?;
                }
                sub.markup = mml::tag("mrow", &sub.markup, 0, &[]);
                sub.latex = subs.join(sept);
                sub
            }
            x => self.minus_prepended(x, flags)?,
        };
        self.wrap(&mut ml, flags, "<selector/>", &tmp.content, &sub.content);
        if flags.pres() {
            let m = format!("{}{}", tmp.markup, sub.markup);
            ml.markup = mml::tag("msub", &m, self.id(), &[]);
        }
        if flags.tex() {
            ml.latex = format!("{}_{{{}}}", tmp.latex, sub.latex);
        }
        Ok(ml)
    }

    /// :3532-3611: `f\left(x\right)` for an unknown `f` (`of`), a named
    /// function by `func2markup`; one argument of an elementary function
    /// without parentheses when it is an atom (`\sin x`). Several
    /// arguments are not minus-prepended (giac bug: `f(-x,y)` prints
    /// `f\left(x,y\right)`). MathML: content `<apply>head args</apply>`,
    /// presentation `<mrow>head<mfenced>args</mfenced></mrow>` or
    /// `<mrow>head<mo>&af;</mo>arg</mrow>`.
    fn general(&self, g: &Gen, flags: Flags, mut ml: MarkupBlock) -> Block {
        let Gen::Symb(op, arg) = g else {
            unreachable!("a symbolic");
        };
        ml.priority = priority::APPLY;
        let gp = if op == "of"
            && let [f, gp] = g.args()
        {
            let mut tmp = self.minus_prepended(f, flags)?;
            if !tmp.ctype(kind::DERIVATIVE) && tmp.priority >= ml.priority {
                parenthesize(&mut tmp, flags);
            }
            ml.content = tmp.content;
            if flags.pres() {
                ml.markup = tmp.markup + mml::APPLY;
            }
            ml.latex = tmp.latex;
            gp
        } else {
            if flags.cont() {
                ml.content = func2mml(op, self.next_id(), true);
            }
            if flags.pres() {
                ml.markup = func2mml(op, self.id(), false);
            }
            if flags.tex() {
                ml.latex = func2markup(op);
            }
            &**arg
        };
        if let Gen::Vect(items, Sub::Seq) = gp {
            let mut args = MarkupBlock::new();
            let mut texs = Vec::new();
            for x in items {
                let tmp = self.gen2markup(x, flags)?;
                args.content.push_str(&tmp.content);
                args.markup.push_str(&tmp.markup);
                texs.push(tmp.latex);
                self.check(ml.len() + args.len())?;
            }
            if flags.cont() {
                let c = format!("{}{}", ml.content, args.content);
                ml.content = mml::tag("apply", &c, self.next_id(), &[]);
            }
            if flags.pres() {
                let m = format!("{}{}", ml.markup, mml::tag("mfenced", &args.markup, 0, &[]));
                ml.markup = mml::tag("mrow", &m, self.id(), &[]);
            }
            if flags.tex() {
                ml.latex = format!("{}\\left({}\\right)", ml.latex, texs.join(","));
            }
        } else {
            let mut tmp = self.minus_prepended(gp, flags)?;
            if op != "of" && is_elemfunc(op) {
                ml.kind = kind::ELEMAPP;
                ml.split_pos = i32::try_from(ml.markup.len()).expect("a short text");
                ml.split_pos_tex = ml.latex.len();
                ml.appl = tmp.priority == 0
                    && !tmp.ctype(kind::NUMERIC_EXACT)
                    && !tmp.ctype(kind::NUMERIC_APPROX)
                    && !tmp.ctype(kind::SUBTYPE_IDNT)
                    && !tmp.ctype(kind::IDNT_NAME);
            }
            if !ml.appl {
                parenthesize(&mut tmp, flags);
            }
            if flags.cont() {
                let c = format!("{}{}", ml.content, tmp.content);
                ml.content = mml::tag("apply", &c, self.next_id(), &[]);
            }
            if flags.pres() {
                let m = format!("{}{}{}", ml.markup, mml::APPLY, tmp.markup);
                ml.markup = mml::tag("mrow", &m, self.id(), &[]);
            }
            if flags.tex() {
                ml.latex.push_str(&tmp.latex);
            }
        }
        Ok(ml)
    }
}

/// :1380-1388: a double's text with an exponent (nonzero mantissa and
/// exponent) as mantissa (with a point), sign and exponent digits, written
/// `mant\times10^{ex}` of priority MUL: `1e-10` is `1.0\times10^{-10}`,
/// `1.23e+20` `1.23\times10^{20}`. `None` for a text written as it is.
fn double_parts(s: &str) -> Option<(String, bool, &str)> {
    use super::ident::atof;
    let cpos = s.find('e')?;
    let mant = &s[..cpos];
    let after = s.get(cpos + 2..).unwrap_or("");
    if atof(mant) == 0.0 || atof(after) == 0.0 {
        return None;
    }
    let ex = after.trim_start_matches('0');
    let point = if mant.contains('.') { "" } else { ".0" };
    Some((
        format!("{mant}{point}"),
        s[cpos + 1..].starts_with('-'),
        ex,
    ))
}
