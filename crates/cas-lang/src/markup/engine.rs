//! giac's `gen2markup` (markup.cc :1233-3627), its LaTeX arms: a [`Gen`]
//! to a [`MarkupBlock`]. Each function cites the lines it ports; the
//! prototype `.run/1.6/T4/proto/markup.py` (1.6/T4.c1) was checked against
//! giac 2.0.0.21's `latex()` on 301 inputs, and these are its steps.

use cas_core::Integer;

use super::ident::{
    func2markup, idnt2markup, is_elemfunc, is_greek_letter, my_isalpha, string2markup,
};
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
    pub(in super::super) const POWER: u32 = 2048;
    pub(in super::super) const DERIVATIVE: u32 = 4096;
    pub(in super::super) const HAS_SUBSCRIPT: u32 = 8192;
    pub(in super::super) const LEADING_DIGIT: u32 = 16384;
    pub(in super::super) const IDNT_NAME: u32 = 32768;
}

/// The flags of a conversion (`_MARKUP_*`, :103-115). `MathML` adds its own.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Flags(u32);

impl Flags {
    pub(super) const LATEX: Flags = Flags(1);
    pub(super) const TOPLEVEL: Flags = Flags(2);
    pub(super) const ELEMPOW: Flags = Flags(32);
    pub(super) const FACTOR: Flags = Flags(64);

    fn has(self, f: Flags) -> bool {
        self.0 & f.0 != 0
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

/// giac's `MarkupBlock` (:179-192), its LaTeX fields.
#[derive(Clone, Debug)]
pub(super) struct MarkupBlock {
    /// How weakly the text binds: 0 for an atom, see [`priority`].
    pub(super) priority: u8,
    /// [`kind`] bits.
    pub(super) kind: u32,
    /// A minus sign that belongs in front, not in the text yet.
    pub(super) neg: bool,
    pub(super) latex: String,
    /// Where an elementary function's application splits (`\sin |x`), in
    /// the markup text: 0 in a LaTeX-only run, -1 for none.
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
            split_pos: -1,
            split_pos_tex: 0,
            appl: false,
        }
    }

    fn ctype(&self, t: u32) -> bool {
        self.kind & t != 0
    }
}

/// The text grew past the length allowed.
#[derive(Debug)]
pub(super) struct TooLong;

type Block = Result<MarkupBlock, TooLong>;

/// :958: `\left(…\right)`, an atom now.
fn parenthesize(ml: &mut MarkupBlock) {
    ml.latex = format!("\\left({}\\right)", ml.latex);
    ml.priority = 0;
}

/// :972: the pending minus written in front; a product now, of kind `NEG`.
pub(super) fn prepend_minus(ml: &mut MarkupBlock) {
    if !ml.neg {
        return;
    }
    if ml.priority > priority::MUL {
        parenthesize(ml);
    }
    ml.latex.insert(0, '-');
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

/// The converter, with the longest text allowed: every block is checked as
/// it is built, so a long text is refused before it is all made.
pub(super) struct Engine {
    pub(super) max_len: usize,
}

impl Engine {
    fn check(&self, s: &str) -> Result<(), TooLong> {
        if s.len() > self.max_len {
            Err(TooLong)
        } else {
            Ok(())
        }
    }

    /// :3629 `export_latex`: `g` at the top level, its minus written.
    pub(super) fn export_latex(&self, g: &Gen) -> Result<String, TooLong> {
        let mut ml = self.gen2markup(g, Flags::TOPLEVEL | Flags::LATEX | Flags::ELEMPOW)?;
        prepend_minus(&mut ml);
        self.check(&ml.latex)?;
        Ok(ml.latex)
    }

    /// :1233 `gen2markup`, its text checked against the cap.
    pub(super) fn gen2markup(&self, g: &Gen, flags: Flags) -> Block {
        let ml = self.markup(g, flags)?;
        self.check(&ml.latex)?;
        Ok(ml)
    }

    /// A child's block, its minus written.
    fn minus_prepended(&self, g: &Gen, flags: Flags) -> Block {
        let mut ml = self.gen2markup(g, flags)?;
        prepend_minus(&mut ml);
        Ok(ml)
    }

    fn markup(&self, g: &Gen, flags_orig: Flags) -> Block {
        let toplevel = flags_orig.has(Flags::TOPLEVEL);
        let isfactor = flags_orig.has(Flags::FACTOR);
        let flags = flags_orig.without(Flags::TOPLEVEL).without(Flags::FACTOR);
        let mut ml = MarkupBlock::new();
        match g {
            // :1308-1321.
            Gen::Bool(b) => {
                ml.latex = idnt2markup(if *b { "true" } else { "false" }).0;
                ml.kind = kind::SUBTYPE_IDNT | kind::LEADING_DIGIT;
            }
            // :1323-1337.
            Gen::Int(n) => {
                ml.kind = kind::NUMERIC_EXACT | kind::LEADING_DIGIT;
                ml.neg = n.is_negative();
                ml.latex = n.abs().to_string();
            }
            // :1362-1407.
            Gen::Dbl { text, neg } => {
                ml.kind = kind::NUMERIC_APPROX | kind::LEADING_DIGIT;
                ml.neg = *neg;
                ml.latex = double_latex(text, &mut ml.priority);
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
                    ml.latex = format!("\\frac{{{}}}{{{}}}", n.abs(), d.abs());
                }
            }
            // :1481-1485.
            Gen::Str(s) => ml.latex = string2markup(s),
            // :1486-1542.
            Gen::Idnt(name) => match name.as_str() {
                "infinity" | "unsigned_inf" => ml.latex.push_str("\\infty "),
                "pi" => ml.latex.push_str("\\pi "),
                "euler_gamma" => ml.latex.push_str("\\gamma "),
                _ => {
                    let (latex, isname) = idnt2markup(name);
                    if latex.contains("_{") {
                        ml.kind |= kind::HAS_SUBSCRIPT;
                    }
                    if isname {
                        ml.kind |= kind::IDNT_NAME;
                    }
                    ml.latex = latex;
                }
            },
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
                ml.latex.push_str("\\mathrm{i}");
                ml.neg = im.is_minus_one();
            } else {
                let mut tmp = self.gen2markup(im, flags | Flags::FACTOR)?;
                ml.neg = tmp.neg;
                ml.priority = priority::MUL;
                if tmp.priority >= ml.priority {
                    parenthesize(&mut tmp);
                }
                ml.latex = format!("{} \\mathrm{{i}}", tmp.latex);
            }
            return Ok(ml);
        }
        ml.priority = priority::ADD;
        let mut left = self.minus_prepended(re, flags)?;
        if left.priority >= ml.priority {
            parenthesize(&mut left);
        }
        let mut right = self.gen2markup(im, flags)?;
        if right.priority >= priority::MUL {
            parenthesize(&mut right);
        }
        ml.latex = left.latex;
        ml.latex.push(if right.neg { '-' } else { '+' });
        if !isone {
            ml.latex.push_str(&right.latex);
            ml.latex.push(' ');
        }
        ml.latex.push_str("\\mathrm{i}");
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
                for x in row {
                    cells.push(self.minus_prepended(x, flags)?.latex);
                }
                rows.push(cells.join("&"));
                self.check(rows.last().expect("pushed"))?;
            }
            ml.latex = format!(
                "\\left[\\begin{{array}}{{{}}}{}\\end{{array}}\\right]",
                "c".repeat(cols),
                rows.join("\\\\")
            );
            return Ok(ml);
        }
        if sub == Sub::Set && items.is_empty() {
            ml.latex.push_str("\\emptyset ");
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
            ml.latex.push_str(&self.minus_prepended(x, flags)?.latex);
            self.check(&ml.latex)?;
        }
        ml.latex.push_str(rd);
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
                            ml.latex.push('+');
                        }
                        prepend_minus(&mut tmp);
                        if tmp.priority >= ml.priority {
                            parenthesize(&mut tmp);
                        }
                        ml.latex.push_str(&tmp.latex);
                        self.check(&ml.latex)?;
                    }
                } else {
                    // :1827: unary plus.
                    ml = self.gen2markup(arg, flags)?;
                    if !ml.neg && ml.priority == 0 {
                        ml.latex.insert(0, '+');
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
                    parenthesize(&mut left);
                }
                if right.priority >= ml.priority {
                    parenthesize(&mut right);
                }
                let sign = if right.neg { "\\oplus " } else { "\\ominus " };
                ml.latex = format!("{}{sign}{}", left.latex, right.latex);
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
                        parenthesize(&mut tmp);
                    }
                    if k > 0 {
                        ml.latex.push_str("\\odot ");
                    }
                    ml.latex.push_str(&tmp.latex);
                    self.check(&ml.latex)?;
                }
                return Ok(ml);
            }
            // :2082: `\frac1{X}` for any X.
            "inv" => {
                let tmp = self.gen2markup(arg, flags)?;
                ml.neg = tmp.neg;
                ml.kind = kind::FRACTION;
                ml.priority = priority::MUL;
                ml.latex = format!("\\frac1{{{}}}", tmp.latex);
                return Ok(ml);
            }
            // :2097.
            "division" | "rdiv" if isbinary => {
                ml.kind = kind::FRACTION;
                ml.priority = priority::MUL;
                let num = self.gen2markup(&args[0], flags)?;
                let den = self.gen2markup(&args[1], flags)?;
                ml.neg = num.neg != den.neg;
                ml.latex = format!("\\frac{{{}}}{{{}}}", num.latex, den.latex);
                return Ok(ml);
            }
            // :2130.
            "conj" => {
                let mut tmp = self.gen2markup(arg, flags)?;
                if isfactor && tmp.priority > priority::MUL {
                    parenthesize(&mut tmp);
                }
                if tmp.ctype(kind::LEADING_DIGIT) {
                    ml.kind = kind::LEADING_DIGIT;
                }
                ml.priority = priority::UNARY;
                ml.neg = tmp.neg;
                ml.latex = format!("\\overline{{{}}}", tmp.latex);
                return Ok(ml);
            }
            // :2148.
            "exp" | "EXP" => {
                if arg.is_one() {
                    ml.latex.push_str("\\mathrm{e}");
                } else {
                    let tmp = self.minus_prepended(arg, flags)?;
                    ml.priority = priority::EXP;
                    ml.kind = kind::POWER;
                    ml.latex = format!("\\mathrm{{e}}^{{{}}}", tmp.latex);
                }
                return Ok(ml);
            }
            // :2191.
            "sqrt" => {
                ml.priority = priority::UNARY;
                ml.kind = kind::ROOT;
                let tmp = self.minus_prepended(arg, flags)?;
                ml.latex = format!("\\sqrt{{{}}}", tmp.latex);
                return Ok(ml);
            }
            // :2224-2280.
            "pow" | "matpow" | "pointpow" if isbinary => return self.pow(op, args, flags, ml),
            // :2281: abs drops the child's pending minus (giac quirk).
            "abs" | "floor" | "ceil" => {
                let mut tmp = self.gen2markup(arg, flags)?;
                if op != "abs" {
                    prepend_minus(&mut tmp);
                }
                let (ld, rd) = match op {
                    "abs" => ("|", "|"),
                    "floor" => ("\\lfloor ", "\\rfloor "),
                    _ => ("\\lceil ", "\\rceil "),
                };
                ml.latex = format!("\\left{ld}{}\\right{rd}", tmp.latex);
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
                let opt = match op {
                    "same" => "==",
                    "different" => "\\neq ",
                    "equal" => "=",
                    "sto" => ":=",
                    _ => "\\leftarrow",
                };
                ml.latex = format!("{}{opt}{}", left.latex, right.latex);
                return Ok(ml);
            }
            // :2770.
            "inferieur_strict" | "inferieur_egal" | "superieur_strict" | "superieur_egal"
                if isbinary =>
            {
                ml.priority = priority::INEQ;
                let (left, right) = self.leftright(args, Some(ml.priority), flags)?;
                let opt = match op {
                    "inferieur_strict" => "<",
                    "inferieur_egal" => "\\leq ",
                    "superieur_strict" => ">",
                    _ => "\\geq ",
                };
                ml.latex = format!("{}{opt}{}", left.latex, right.latex);
                return Ok(ml);
            }
            // :2942.
            "interval" if isbinary => {
                ml.priority = priority::SET;
                let (left, right) = self.leftright(args, Some(ml.priority), flags)?;
                ml.latex = format!("{}\\ldots {}", left.latex, right.latex);
                return Ok(ml);
            }
            // :3317.
            "factorial" => {
                ml.priority = priority::UNARY;
                ml.kind = kind::FACTORIAL;
                let mut tmp = self.minus_prepended(arg, flags)?;
                if tmp.priority >= ml.priority {
                    parenthesize(&mut tmp);
                }
                ml.latex = format!("{}\\mathclose{{!}}", tmp.latex);
                return Ok(ml);
            }
            // :3334.
            "binomial" | "nCr" | "comb" if isbinary => {
                let left = self.minus_prepended(&args[0], flags)?;
                let right = self.minus_prepended(&args[1], flags)?;
                ml.latex = format!("\\binom{{{}}}{{{}}}", left.latex, right.latex);
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
            parenthesize(&mut left);
        }
        let mut right = self.minus_prepended(&args[args.len() - 1], flags)?;
        if priority.is_some_and(|p| right.priority >= p) {
            parenthesize(&mut right);
        }
        Ok((left, right))
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
            ml.latex = format!("\\sqrt{{{}}}", left.latex);
            return Ok(ml);
        }
        let right = self.minus_prepended(e, flags)?;
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
            let k = left.split_pos_tex;
            ml.latex = format!(
                "{}^{{{}}}{}",
                &left.latex[..k],
                right.latex,
                &left.latex[k..]
            );
        } else {
            ml.kind = kind::POWER;
            if left.ctype(kind::LEADING_DIGIT) {
                ml.kind |= kind::LEADING_DIGIT;
            }
            if left.appl || left.priority >= ml.priority {
                parenthesize(&mut left);
            }
            ml.latex = format!("{}^{{{}}}", left.latex, right.latex);
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
        let (mut numt, mut dent) = (String::new(), String::new());
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
                parenthesize(&mut tmp);
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
            let sign = if is_cdot { "\\cdot " } else { " " };
            let side = if isinv { &mut dent } else { &mut numt };
            if !side.is_empty() {
                side.push_str(sign);
            }
            side.push_str(&tmp.latex);
            self.check(side)?;
            if isinv {
                pden = tmp;
            } else {
                pnum = tmp;
            }
        }
        if neg_count % 2 == 1 {
            ml.neg = !ml.neg;
        }
        if numt.is_empty() {
            numt.push('1');
        }
        if dent.is_empty() {
            ml.latex = numt;
        } else {
            ml.kind = kind::FRACTION;
            ml.latex = format!("\\frac{{{numt}}}{{{dent}}}");
        }
        Ok(ml)
    }

    /// :2375-2487: `\int _{a}^{b}f\,\mathrm{d}x`, `\sum _{k=1}^{n}…`,
    /// `\prod _{k}…`.
    fn bigop(&self, g: &Gen, flags: Flags, mut ml: MarkupBlock) -> Block {
        let Gen::Symb(op, arg) = g else {
            unreachable!("a symbolic");
        };
        let isint = matches!(op.as_str(), "integrate" | "int" | "Int");
        let isprod = matches!(op.as_str(), "product" | "mul");
        let mut big_tex = String::from(if isint {
            "\\int "
        } else if isprod {
            "\\prod "
        } else {
            "\\sum "
        });
        ml.priority = priority::MUL;
        ml.kind = kind::SUMPROD;
        let mut var = MarkupBlock::new();
        let mut e;
        if let Gen::Vect(a, _) = &**arg {
            e = self.gen2markup(&a[0], flags)?;
            let (mut lb, mut ub) = (None, None);
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
                prepend_minus(&mut var);
                if isint {
                    if var.priority != 0 {
                        parenthesize(&mut var);
                    }
                    var.latex = format!("\\,\\mathrm{{d}}{}", var.latex);
                }
            }
            if lb.is_none() && a.len() > 2 {
                lb = Some(self.gen2markup(&a[2], flags)?);
            }
            if ub.is_none() && a.len() > 3 {
                ub = Some(self.gen2markup(&a[3], flags)?);
            }
            if let Some(lb) = &mut lb {
                prepend_minus(lb);
            }
            if let Some(ub) = &mut ub {
                prepend_minus(ub);
            }
            if !isint {
                lb = Some(match lb {
                    Some(mut lb) => {
                        lb.latex = format!("{}={}", var.latex, lb.latex);
                        lb
                    }
                    None => var.clone(),
                });
            }
            match (lb, ub) {
                (Some(lb), None) => big_tex = format!("{big_tex}_{{{}}}", lb.latex),
                (Some(lb), Some(ub)) => {
                    big_tex = format!("{big_tex}_{{{}}}^{{{}}}", lb.latex, ub.latex);
                }
                _ => {}
            }
        } else {
            e = self.gen2markup(arg, flags)?;
        }
        if isprod {
            prepend_minus(&mut e);
        } else {
            ml.neg = e.neg;
        }
        if (isprod && !e.ctype(kind::FRACTION) && e.priority >= ml.priority)
            || (!isprod && (e.ctype(kind::NEG) || e.priority > ml.priority))
        {
            parenthesize(&mut e);
        }
        ml.latex = big_tex + &e.latex;
        if isint {
            ml.latex.push_str(&var.latex);
        }
        Ok(ml)
    }

    /// :2488-2548: `\lim_{x\to a}`, `^{+}`/`^{-}` for a direction ±1.
    fn limit(&self, a: &[Gen], flags: Flags, mut ml: MarkupBlock) -> Block {
        ml.priority = priority::MUL;
        ml.kind = kind::SUMPROD;
        let mut dirn = 0i64;
        let mut e = self.gen2markup(&a[0], flags)?;
        ml.neg = e.neg;
        if e.priority > ml.priority {
            parenthesize(&mut e);
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
        prepend_minus(&mut var);
        prepend_minus(&mut dest);
        let side = match dirn.signum() {
            0 => "",
            -1 => "^{-}",
            _ => "^{+}",
        };
        ml.latex = format!(
            "\\lim_{{{}\\to {}{side}}}{}",
            var.latex, dest.latex, e.latex
        );
        Ok(ml)
    }

    /// :2560-2727: partial derivatives of a function `f^{\left(1,2\right)}`,
    /// `f'`, `f''`, `f'''`, `f^{(4)}`, and `derive(e, vars…)` as a quotient
    /// of differentials. `false` when giac falls through to the general
    /// arm, with `ml`'s fields set so far.
    #[allow(clippy::too_many_lines)]
    fn derive(&self, g: &Gen, flags: Flags, ml: &mut MarkupBlock) -> Result<bool, TooLong> {
        if let Some((f, deg)) = is_partialdiff(g) {
            ml.priority = priority::EXP;
            ml.kind = kind::POWER | kind::DERIVATIVE;
            let left = self.gen2markup(f, flags)?;
            let right = self.gen2markup(&Gen::Vect(deg, Sub::Seq), flags)?;
            ml.latex = format!("{}^{{{}}}", left.latex, right.latex);
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
            parenthesize(&mut tmp);
        }
        ml.neg = tmp.neg;
        if n > 0 {
            ml.priority = priority::EXP;
            ml.kind |= kind::POWER;
            ml.latex = tmp.latex;
            match n + 1 {
                2 => ml.latex.push_str("''"),
                3 => ml.latex.push_str("'''"),
                k => {
                    ml.latex.push_str("^{(");
                    ml.latex.push_str(&k.to_string());
                    ml.latex.push_str(")}");
                }
            }
            return Ok(true);
        }
        if !vectarg {
            ml.priority = priority::EXP;
            ml.latex = tmp.latex + "'";
            return Ok(true);
        }
        let mut vars = Vec::new();
        if items.len() > 1 && get_derive_vars(items, &mut vars) {
            ml.priority = priority::MUL;
            ml.kind |= kind::FRACTION;
            if let [(v, 1)] = vars.as_slice() {
                let var = self.gen2markup(v, flags)?;
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
                return Ok(true);
            }
            let ds = if vars.len() > 1 {
                "\\partial "
            } else {
                "\\mathrm{d}"
            };
            let mut parts = Vec::new();
            for (v, k) in &vars {
                let var = self.gen2markup(v, flags)?;
                let power = if *k == 1 {
                    String::new()
                } else {
                    format!("^{{{k}}}")
                };
                parts.push(format!("{ds}{}{power}", var.latex));
            }
            let body = parts.join(" ");
            let ns: u64 = vars.iter().map(|(_, k)| k).sum();
            ml.latex = if simp {
                format!("\\frac{{{ds}^{{{ns}}}{}}}{{{body}}}", tmp.latex)
            } else {
                format!("\\frac{{{ds}^{{{ns}}}}}{{{body}}}{}", tmp.latex)
            };
            return Ok(true);
        }
        Ok(false)
    }

    /// :3434: `x_{1}`, `a_{i j}`, `a_{i,10}`.
    fn at(&self, a: &[Gen], flags: Flags, mut ml: MarkupBlock) -> Block {
        ml.priority = priority::APPLY;
        ml.kind = kind::HAS_SUBSCRIPT;
        let mut tmp = self.minus_prepended(&a[0], flags)?;
        if tmp.priority > ml.priority || tmp.ctype(kind::HAS_SUBSCRIPT) {
            parenthesize(&mut tmp);
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
                let sept = if idx.iter().all(one_char) { " " } else { "," };
                let mut subs = Vec::new();
                for x in idx {
                    subs.push(self.minus_prepended(x, flags)?.latex);
                }
                subs.join(sept)
            }
            x => self.minus_prepended(x, flags)?.latex,
        };
        ml.latex = format!("{}_{{{sub}}}", tmp.latex);
        Ok(ml)
    }

    /// :3532-3611: `f\left(x\right)` for an unknown `f` (`of`), a named
    /// function by `func2markup`; one argument of an elementary function
    /// without parentheses when it is an atom (`\sin x`). Several
    /// arguments are not minus-prepended (giac bug: `f(-x,y)` prints
    /// `f\left(x,y\right)`).
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
                parenthesize(&mut tmp);
            }
            ml.latex = tmp.latex;
            gp
        } else {
            ml.latex = func2markup(op);
            &**arg
        };
        if let Gen::Vect(items, Sub::Seq) = gp {
            ml.latex.push_str("\\left(");
            for (k, x) in items.iter().enumerate() {
                if k > 0 {
                    ml.latex.push(',');
                }
                ml.latex.push_str(&self.gen2markup(x, flags)?.latex);
                self.check(&ml.latex)?;
            }
            ml.latex.push_str("\\right)");
        } else {
            let mut tmp = self.minus_prepended(gp, flags)?;
            if op != "of" && is_elemfunc(op) {
                ml.kind = kind::ELEMAPP;
                ml.split_pos = 0;
                ml.split_pos_tex = ml.latex.len();
                ml.appl = tmp.priority == 0
                    && !tmp.ctype(kind::NUMERIC_EXACT)
                    && !tmp.ctype(kind::NUMERIC_APPROX)
                    && !tmp.ctype(kind::SUBTYPE_IDNT)
                    && !tmp.ctype(kind::IDNT_NAME);
            }
            if !ml.appl {
                parenthesize(&mut tmp);
            }
            ml.latex.push_str(&tmp.latex);
        }
        Ok(ml)
    }
}

/// :1362-1407: a double's text; with an exponent (nonzero mantissa and
/// exponent), `mant\times10^{ex}` of priority MUL: `1e-10` is
/// `1.0\times10^{-10}`, `1.23e+20` `1.23\times10^{20}`.
fn double_latex(s: &str, priority: &mut u8) -> String {
    use super::ident::atof;
    let Some(cpos) = s.find('e') else {
        return s.to_owned();
    };
    let mant = &s[..cpos];
    let after = s.get(cpos + 2..).unwrap_or("");
    if atof(mant) == 0.0 || atof(after) == 0.0 {
        return s.to_owned();
    }
    *priority = priority::MUL;
    let ex = after.trim_start_matches('0');
    let point = if mant.contains('.') { "" } else { ".0" };
    let minus = if s[cpos + 1..].starts_with('-') {
        "-"
    } else {
        ""
    };
    format!("{mant}{point}\\times10^{{{minus}{ex}}}")
}
