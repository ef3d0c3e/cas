//! giac's print of a `gen` (`gen::print`, gen.cc :14569; `add_print`,
//! symbolic.cc :175-663; `printsommetasoperator`, unary.cc :599), for the
//! texts `mathml.cc` takes printed, and giac's tests on numbers and heads
//! the converter shares (`need_parenthesis`, usual.cc :8218). Heads with a
//! print of their own beyond relations, intervals, `of`, `at`,
//! `factorial`, `sqrt` and `==` print as `f(args)`.

use std::cmp::Ordering;

use crate::markup::ident::my_isalpha;
use crate::markup::tree::{Gen, Sub};

/// giac's `gen` types, in their order (gen.h), as far as `mathml.cc` and
/// the print compare them.
pub(super) mod ty {
    pub(in crate::mathml) const INT: u8 = 0;
    pub(in crate::mathml) const DOUBLE: u8 = 1;
    pub(in crate::mathml) const ZINT: u8 = 2;
    pub(in crate::mathml) const REAL: u8 = 3;
    pub(in crate::mathml) const CPLX: u8 = 4;
    pub(in crate::mathml) const IDNT: u8 = 6;
    pub(in crate::mathml) const VECT: u8 = 7;
    pub(in crate::mathml) const SYMB: u8 = 8;
    pub(in crate::mathml) const FRAC: u8 = 10;
    pub(in crate::mathml) const STRNG: u8 = 12;
    pub(in crate::mathml) const FUNC: u8 = 13;
}

/// The heads giac registers with a `printsommet` (their own print), by
/// giac name, sorted: the `define_unary_function_eval*` of giac's sources
/// with one (`.run/1.6/T5/c3/printsommet.py`). `need_parenthesis` reads it.
const PRINTSOMMET: &[&str] = &[
    "Ans",
    "Archive",
    "Circle",
    "ClrDraw",
    "ClrGraph",
    "ClrIO",
    "CopyVar",
    "CyclePic",
    "DelFold",
    "DelVar",
    "Dialog",
    "Digits",
    "DispG",
    "DispHome",
    "DrawFunc",
    "DrawInv",
    "DrawParm",
    "DrawPol",
    "DrawSlp",
    "DropDown",
    "DrwCtour",
    "Exec",
    "Get",
    "GetCalc",
    "GetFold",
    "Graph",
    "HAngle",
    "HComplex",
    "HDigits",
    "HFormat",
    "HLanguage",
    "IFTE",
    "Input",
    "InputStr",
    "Line",
    "LineHorz",
    "LineTan",
    "LineVert",
    "MOD",
    "NOP",
    "NTHROOT",
    "NewFold",
    "NewPic",
    "Output",
    "PIECEWISE",
    "POLYFORM",
    "Pause",
    "PopUp",
    "Popup",
    "Prompt",
    "PtOff",
    "PtOn",
    "PtText",
    "PxlOff",
    "PxlOn",
    "RANDOM",
    "RPN_CASE",
    "RPN_FOR",
    "RPN_LOCAL",
    "RPN_UNTIL",
    "RPN_WHILE",
    "RandSeed",
    "RclPic",
    "Request",
    "RplcPic",
    "SetFold",
    "StoPic",
    "Text",
    "Title",
    "Unarchiv",
    "ZoomRcl",
    "ZoomSto",
    "alg",
    "all_trig_solutions",
    "ampersand_times",
    "and",
    "angle_radian",
    "approx_mode",
    "at",
    "avance",
    "backquote",
    "backward",
    "baisse_crayon",
    "binary_minus",
    "bloc",
    "break",
    "cache_tortue",
    "cap",
    "case",
    "cell",
    "check_type",
    "clearscreen",
    "comment",
    "complex_mode",
    "complex_variables",
    "compose",
    "composepow",
    "continue",
    "crayon",
    "curve",
    "deSolve",
    "debut_enregistrement",
    "decrement",
    "del",
    "derive",
    "dessine_tortue",
    "deuxpoints",
    "different",
    "disque",
    "disque_centre",
    "div",
    "divcrement",
    "division",
    "dollar",
    "double_deux_points",
    "droite_tangente",
    "ecris",
    "efface",
    "efface_logo",
    "epsilon",
    "equal",
    "equal2",
    "erase",
    "et",
    "eval_level",
    "factorial",
    "feuille",
    "fin_enregistrement",
    "for",
    "gbasis_max_pairs",
    "gbasis_reinject",
    "gbasis_simult_primes",
    "goto",
    "hasard",
    "ifte",
    "increasing_power",
    "increment",
    "inferieur_egal",
    "inferieur_strict",
    "integer_format",
    "intersect",
    "interval",
    "inverser",
    "iquo",
    "irem",
    "jump",
    "keep_algext",
    "label",
    "leftopen_interval",
    "leftrightopen_interval",
    "leve_crayon",
    "lis",
    "lis_phrase",
    "local",
    "maple_if",
    "maple_lib",
    "maple_op",
    "maple_subs",
    "maple_subsop",
    "minus",
    "montre_tortue",
    "multcrement",
    "nodisp",
    "normalmod",
    "not",
    "ntl_on",
    "of",
    "ou",
    "oufr",
    "pas_de_cote",
    "piecewise",
    "plus",
    "pnt",
    "pointdivision",
    "pointminus",
    "pointplus",
    "pointpow",
    "pointprod",
    "polar_complex",
    "polygone_rempli",
    "position",
    "pour",
    "pourcent",
    "pow",
    "print",
    "proba_epsilon",
    "prod",
    "prog_eval_level",
    "program",
    "quote",
    "quote_pow",
    "ramene",
    "rectangle_plein",
    "recule",
    "repete",
    "reset",
    "return",
    "rightopen_interval",
    "rond",
    "rootof",
    "rpn",
    "rpn_prog",
    "rur_separate_max_tries",
    "same",
    "saute",
    "sauve",
    "scientific_format",
    "si",
    "sialorssinon",
    "signe",
    "size",
    "skip",
    "speed",
    "sqrt",
    "sto",
    "struct_dot",
    "subs",
    "subsop",
    "superieur_egal",
    "superieur_strict",
    "symmetric_difference",
    "tangente",
    "tantque",
    "threads",
    "ti_semi",
    "tilocal",
    "tourne_droite",
    "tourne_gauche",
    "towards",
    "triangle_plein",
    "try_catch",
    "turtle_stack",
    "union",
    "unit",
    "vers",
    "virgule",
    "when",
    "with_sqrt",
    "xor",
];

/// Whether giac's `op` has a print of its own.
pub(super) fn has_printsommet(op: &str) -> bool {
    PRINTSOMMET.binary_search(&op).is_ok()
}

/// giac's type of `g`.
pub(super) fn gtype(g: &Gen) -> u8 {
    match g {
        Gen::Int(n) if n.to_i64().is_some_and(|v| i32::try_from(v).is_ok()) => ty::INT,
        Gen::Int(_) => ty::ZINT,
        Gen::Bool(_) => ty::INT,
        Gen::Dbl { .. } => ty::DOUBLE,
        Gen::Cplx(..) => ty::CPLX,
        Gen::Idnt(_) => ty::IDNT,
        Gen::Vect(..) => ty::VECT,
        Gen::Symb(..) => ty::SYMB,
        Gen::Frac(..) => ty::FRAC,
        Gen::Str(_) => ty::STRNG,
        Gen::Func(_) => ty::FUNC,
    }
}

/// The sign of a real number (an integer, a boolean, a fraction, a
/// double), `None` for anything else.
pub(super) fn sign(g: &Gen) -> Option<Ordering> {
    match g {
        Gen::Int(n) => Some(if n.is_negative() {
            Ordering::Less
        } else if n.is_zero() {
            Ordering::Equal
        } else {
            Ordering::Greater
        }),
        Gen::Bool(b) => Some(if *b {
            Ordering::Greater
        } else {
            Ordering::Equal
        }),
        Gen::Frac(a, _) => match &**a {
            Gen::Int(_) | Gen::Bool(_) => sign(a),
            _ => None,
        },
        Gen::Dbl { text, neg } => {
            if text.parse::<f64>().is_ok_and(|v| v == 0.0) {
                Some(Ordering::Equal)
            } else if *neg {
                Some(Ordering::Less)
            } else {
                Some(Ordering::Greater)
            }
        }
        _ => None,
    }
}

/// giac's `is_positive` as the printers use it: a number not below 0, or
/// anything that is not a real number.
pub(super) fn is_positive(g: &Gen) -> bool {
    sign(g) != Some(Ordering::Less)
}

/// A real number below 0.
pub(super) fn is_negative(g: &Gen) -> bool {
    sign(g) == Some(Ordering::Less)
}

/// Whether `g` is the number `p/q` (a boolean is 0 or 1).
pub(super) fn equals(g: &Gen, p: i64, q: i64) -> bool {
    match g {
        Gen::Bool(b) => q == 1 && p == i64::from(*b),
        _ => g.equals(p, q),
    }
}

/// giac's `is_zero`, a complex number's parts both 0.
pub(super) fn is_zero(g: &Gen) -> bool {
    match g {
        Gen::Cplx(a, b) => is_zero(a) && is_zero(b),
        _ => equals(g, 0, 1),
    }
}

pub(super) fn is_one(g: &Gen) -> bool {
    equals(g, 1, 1)
}

/// An integer of giac's (`_INT_`, `_ZINT`; a boolean is an `_INT_`).
pub(super) fn is_int(g: &Gen) -> bool {
    matches!(g, Gen::Int(_) | Gen::Bool(_))
}

/// `-g`, folded into a number.
pub(super) fn neg_num(g: &Gen) -> Gen {
    match g {
        Gen::Int(n) => Gen::Int(-n),
        Gen::Bool(b) => Gen::int(-i64::from(*b)),
        Gen::Dbl { text, neg } => Gen::Dbl {
            text: text.clone(),
            neg: !neg,
        },
        Gen::Frac(a, b) => Gen::Frac(Box::new(neg_num(a)), b.clone()),
        Gen::Cplx(a, b) => Gen::Cplx(Box::new(neg_num(a)), Box::new(neg_num(b))),
        _ => Gen::symb("neg", g.clone()),
    }
}

/// giac's `need_parenthesis` of a head (usual.cc :8233-8242).
pub(super) fn need_par_head(op: &str) -> bool {
    match op {
        "pow" | "division" | "prod" => false,
        "neg" | "inv" | "minus" | "and" | "et" | "ou" | "oufr" | "xor" | "same" | "equal"
        | "equal2" | "superieur_egal" | "superieur_strict" | "inferieur_egal"
        | "inferieur_strict" => true,
        _ => has_printsommet(op),
    }
}

/// giac's `need_parenthesis` (usual.cc :8218).
pub(super) fn need_parenthesis(g: &Gen) -> bool {
    match g {
        Gen::Int(_) | Gen::Bool(_) => is_negative(g),
        Gen::Cplx(a, b) => {
            let (ra, ib) = (neg_num(a), neg_num(b));
            if is_zero(&ra) {
                sign(&ib) == Some(Ordering::Greater)
            } else if is_zero(&ib) {
                sign(&ra) == Some(Ordering::Greater)
            } else {
                true
            }
        }
        Gen::Frac(..) => true,
        Gen::Symb(op, _) | Gen::Func(op) => need_par_head(op),
        _ => false,
    }
}

/// A head's printed name (its `_xxx_s`), where it is not giac's internal
/// name.
pub(super) fn pname(op: &str) -> &str {
    match op {
        "plus" => "+",
        "prod" => "*",
        "pow" => "^",
        "neg" => "-",
        "equal" => "=",
        "same" => "==",
        "different" => "!=",
        "interval" => "..",
        "inferieur_strict" => "<",
        "inferieur_egal" => "<=",
        "superieur_strict" => ">",
        "superieur_egal" => ">=",
        "division" => "/",
        "derive" => "diff",
        "ou" => "or",
        "et" => "and",
        "pointminus" => ".-",
        "equal2" => "%=",
        _ => op,
    }
}

/// An `inequation` for the product's print: a relation.
fn is_inequation(g: &Gen) -> bool {
    g.is_sym(&[
        "superieur_egal",
        "superieur_strict",
        "inferieur_strict",
        "inferieur_egal",
        "different",
        "same",
        "equal",
    ])
}

/// `g.print()` (gen.cc :14569) of a non-symbolic `g`.
pub(super) fn print_atom(g: &Gen) -> String {
    match g {
        Gen::Int(n) => n.to_string(),
        Gen::Bool(b) => (if *b { "true" } else { "false" }).to_owned(),
        Gen::Dbl { text, neg } => {
            if *neg {
                format!("-{text}")
            } else {
                text.clone()
            }
        }
        Gen::Frac(a, b) => format!("{}/{}", gprint(a), gprint(b)),
        Gen::Cplx(re, im) => print_cplx(re, im),
        Gen::Idnt(s) => s.clone(),
        Gen::Str(s) => format!("\"{s}\""),
        Gen::Func(s) => format!("'{s}'"),
        Gen::Vect(items, sub) => {
            let inner = items.iter().map(gprint).collect::<Vec<_>>().join(",");
            match (sub, items.is_empty()) {
                (Sub::Seq, true) => "seq[]".to_owned(),
                (Sub::Set, true) => "set[ ]".to_owned(),
                (Sub::Seq, false) => inner,
                (Sub::Set, false) => format!("set[{inner}]"),
                (Sub::List, _) => format!("[{inner}]"),
            }
        }
        Gen::Symb(..) => gprint(g),
    }
}

/// A complex number (gen.cc :14606-14638): `i`, `-i`, `2*i`, `1+i`,
/// `1-2*i`.
fn print_cplx(re: &Gen, im: &Gen) -> String {
    if is_zero(im) {
        return gprint(re);
    }
    if is_zero(re) {
        return if is_one(im) {
            "i".to_owned()
        } else if equals(im, -1, 1) {
            "-i".to_owned()
        } else {
            format!("{}*i", gprint(im))
        };
    }
    let r = gprint(re);
    if is_one(im) {
        format!("{r}+i")
    } else if equals(im, -1, 1) {
        format!("{r}-i")
    } else if is_negative(im) {
        format!("{r}-{}*i", gprint(&neg_num(im)))
    } else {
        format!("{r}+{}*i", gprint(im))
    }
}

/// `g` printed.
pub(super) fn gprint(g: &Gen) -> String {
    let mut s = String::new();
    add_print(&mut s, g);
    s
}

/// The last character of `s` replaced by `c`.
fn set_last(s: &mut String, c: char) {
    s.pop();
    s.push(c);
}

/// `add_print` (symbolic.cc :601): `g` appended, a `+-` or `--` it makes
/// folded.
fn add_print(s: &mut String, g: &Gen) {
    match g {
        Gen::Symb(op, f) => return add_print_symbolic(s, op, f),
        // add_print_int (:512).
        Gen::Int(n) if gtype(g) == ty::INT => {
            if n.is_negative() {
                if s.ends_with('+') {
                    set_last(s, '-');
                } else if s.ends_with('-') {
                    if s.len() == 1 {
                        s.clear();
                    } else {
                        set_last(s, '+');
                    }
                } else {
                    s.push('-');
                }
            }
            s.push_str(&n.abs().to_string());
            return;
        }
        Gen::Vect(items, Sub::List) => {
            s.push('[');
            s.push_str(&items.iter().map(gprint).collect::<Vec<_>>().join(","));
            s.push(']');
            return;
        }
        Gen::Frac(a, b) if is_int(a) && is_int(b) => {
            add_print(s, a);
            s.push('/');
            add_print(s, b);
            return;
        }
        _ => {}
    }
    let tmp = print_atom(g);
    if let Some(rest) = tmp.strip_prefix('-') {
        if s.ends_with('+') {
            s.pop();
        } else if s.ends_with('-') {
            set_last(s, '+');
            s.push_str(rest);
            return;
        }
    }
    s.push_str(&tmp);
}

/// `add_print_symbolic` (symbolic.cc :542).
fn add_print_symbolic(s: &mut String, op: &str, f: &Gen) {
    if has_printsommet(op) && !matches!(op, "plus" | "prod" | "pow") {
        s.push_str(&printsommet(op, f));
        return;
    }
    if let Gen::Vect(items, sub) = f
        && items.is_empty()
    {
        s.push_str(pname(op));
        s.push_str(match sub {
            Sub::List => "([])",
            Sub::Set => "(set[ ])",
            Sub::Seq => "(NULL)",
        });
        return;
    }
    match op {
        "prod" => return add_print_prod(s, f),
        "plus" => return add_print_plus(s, f),
        "pow" => return add_print_pow(s, f),
        "neg" => return add_print_neg(s, f),
        "inv" => return add_print_inv(s, f),
        _ => {}
    }
    s.push_str(pname(op));
    s.push('(');
    match f {
        Gen::Vect(items, sub) if !matches!(items[0], Gen::Vect(..)) => {
            let (open, close) = match sub {
                Sub::Seq => ("", ""),
                Sub::Set => ("set[", "]"),
                Sub::List => ("[", "]"),
            };
            s.push_str(open);
            for (k, x) in items.iter().enumerate() {
                if k > 0 {
                    s.push(',');
                }
                add_print(s, x);
            }
            s.push_str(close);
        }
        _ => add_print(s, f),
    }
    s.push(')');
}

/// `add_print_plus` (symbolic.cc :175).
fn add_print_plus(s: &mut String, f: &Gen) {
    let Gen::Vect(v, _) = f else {
        s.push('+');
        return add_print(s, f);
    };
    let Some((first, rest)) = v.split_first() else {
        return;
    };
    if rest.is_empty() {
        s.push('+');
    }
    add_print(s, first);
    for x in rest {
        s.push('+');
        let par = x.is_sym(&[
            "superieur_egal",
            "superieur_strict",
            "inferieur_strict",
            "inferieur_egal",
        ]);
        if par {
            s.push('(');
        }
        add_print(s, x);
        if par {
            s.push(')');
        }
    }
}

/// `add_print_neg` (symbolic.cc :209).
fn add_print_neg(s: &mut String, f: &Gen) {
    if s.ends_with('+') {
        set_last(s, '-');
    } else if s.ends_with('-') {
        set_last(s, '+');
    } else {
        s.push('-');
    }
    let bare = match f {
        Gen::Cplx(..) => false,
        Gen::Symb(op, _) => op == "inv" || op == "prod" || !need_parenthesis(f),
        _ => true,
    };
    if bare {
        add_print(s, f);
    } else {
        s.push('(');
        add_print(s, f);
        s.push(')');
    }
}

/// `add_print_inv` (symbolic.cc :231).
fn add_print_inv(s: &mut String, f: &Gen) {
    let (mut f, mut isneg) = match f {
        Gen::Symb(op, a) if op == "neg" && !matches!(**a, Gen::Vect(..)) => ((**a).clone(), true),
        _ => (f.clone(), false),
    };
    if gtype(&f) < ty::CPLX && sign(&f).is_some_and(Ordering::is_le) {
        f = neg_num(&f);
        isneg = !isneg;
    }
    if isneg && s.ends_with('-') {
        if s.len() == 1 {
            s.clear();
        } else {
            set_last(s, '+');
        }
        isneg = false;
    }
    if isneg && s.ends_with('+') {
        set_last(s, '-');
        isneg = false;
    }
    let bt = matches!(&f, Gen::Symb(op, _) if matches!(op.as_str(), "plus" | "inv" | "prod") || need_par_head(op));
    if !bt && !matches!(f, Gen::Cplx(..)) {
        s.push_str(if isneg { "-1/" } else { "1/" });
        add_print(s, &f);
    } else {
        // giac writes "-1", no "/(" (:261).
        s.push_str(if isneg { "-1" } else { "1/(" });
        add_print(s, &f);
        s.push(')');
    }
}

/// `add_print_prod` (symbolic.cc :268; `print_rewrite_prod_inv` off).
fn add_print_prod(s: &mut String, f: &Gen) {
    let Gen::Vect(v, _) = f else {
        return add_print(s, f);
    };
    for (i, e) in v.iter().enumerate() {
        match e {
            Gen::Symb(op, ff) if op == "inv" => {
                if i == 0 {
                    add_print(s, e);
                    continue;
                }
                let par = match &**ff {
                    Gen::Cplx(..) => true,
                    Gen::Symb(fop, _) => {
                        matches!(fop.as_str(), "plus" | "prod" | "inv")
                            || need_par_head(fop)
                            || (fop == "pow"
                                && ff.args().first().is_some_and(|b| b.is_sym(&["inv"])))
                    }
                    _ => false,
                };
                s.push_str(if par { "/(" } else { "/" });
                add_print(s, ff);
                if par {
                    s.push(')');
                }
            }
            Gen::Symb(op, _) => {
                if i > 0 {
                    s.push('*');
                }
                if op == "plus" || op == "neg" || is_inequation(e) {
                    s.push('(');
                    add_print(s, e);
                    s.push(')');
                } else {
                    add_print(s, e);
                }
            }
            _ => {
                if i > 0 {
                    s.push('*');
                }
                if matches!(e, Gen::Cplx(..)) && need_parenthesis(e) {
                    s.push('(');
                    s.push_str(&print_atom(e));
                    s.push(')');
                } else {
                    add_print(s, e);
                }
            }
        }
    }
}

/// `add_print_pow` (symbolic.cc :344).
fn add_print_pow(s: &mut String, f: &Gen) {
    let Gen::Vect(v, _) = f else {
        s.push_str("pow(");
        add_print(s, f);
        s.push(')');
        return;
    };
    let [arg, pui] = v.as_slice() else {
        s.push_str("pow(");
        add_print(s, f);
        s.push(')');
        return;
    };
    let half = matches!(pui, Gen::Frac(..));
    if half && pui.equals(1, 2) {
        s.push_str("sqrt(");
        add_print(s, arg);
        s.push(')');
        return;
    }
    if half && pui.equals(-1, 2) {
        s.push_str("1/sqrt(");
        add_print(s, arg);
        s.push(')');
        return;
    }
    let pui_par = matches!(pui, Gen::Symb(..) | Gen::Frac(..) | Gen::Cplx(..));
    let plain_head = matches!(arg, Gen::Symb(op, _) if op != "neg" && !has_printsommet(op));
    if matches!(arg, Gen::Idnt(_)) || plain_head {
        let argpar = arg.is_sym(&["inv"]);
        if argpar {
            s.push('(');
        }
        add_print(s, arg);
        if argpar {
            s.push(')');
        }
        s.push('^');
        if pui_par {
            s.push('(');
            add_print(s, pui);
            s.push(')');
        } else {
            add_print(s, pui);
        }
        return;
    }
    let argpar = (gtype(arg) > ty::CPLX || !is_positive(arg)) && !matches!(arg, Gen::Idnt(_));
    if argpar {
        s.push('(');
    }
    add_print(s, arg);
    if argpar {
        s.push(')');
    }
    s.push('^');
    let par = pui_par || matches!(pui, Gen::Vect(_, Sub::Seq));
    if par {
        s.push('(');
    }
    add_print(s, pui);
    if par {
        s.push(')');
    }
}

/// `printsommetasoperator` (unary.cc :599): the arguments with `sop`
/// between them.
fn printsommetasoperator(f: &Gen, sop: &str) -> String {
    let Gen::Vect(v, _) = f else {
        return gprint(f);
    };
    let spaced;
    let sop = if sop
        .chars()
        .next()
        .is_some_and(|c| my_isalpha(c) || c == '%' || c == '.')
    {
        spaced = format!(" {sop} ");
        spaced.as_str()
    } else {
        sop
    };
    let Some((a, rest)) = v.split_first() else {
        return String::new();
    };
    let mut s = if matches!(a, Gen::Frac(..)) && sop != "=" {
        format!("({})", gprint(a))
    } else if sop == "=" || matches!(a, Gen::Idnt(_)) || (gtype(a) <= ty::CPLX && is_positive(a)) {
        gprint(a)
    } else {
        format!("({})", gprint(a))
    };
    for b in rest {
        s.push_str(sop);
        let par = matches!(
            b,
            Gen::Symb(..) | Gen::Frac(..) | Gen::Cplx(..) | Gen::Vect(_, Sub::Seq)
        );
        if par && sop != " . " {
            s.push('(');
            s.push_str(&gprint(b));
            s.push(')');
        } else {
            s.push_str(&gprint(b));
        }
    }
    s
}

/// The print of a head with one of its own (usual.cc `printas*`).
fn printsommet(op: &str, f: &Gen) -> String {
    match op {
        // printassame (:6507).
        "same" => format!("({})", printsommetasoperator(f, "==")),
        // printasdifferent (:5711).
        "different" => printsommetasoperator(f, "<>"),
        "equal" | "interval" | "inferieur_strict" | "inferieur_egal" | "superieur_strict"
        | "superieur_egal" | "minus" | "union" | "intersect" | "xor" | "equal2" => {
            printsommetasoperator(f, pname(op))
        }
        "of" | "at" => match f {
            Gen::Vect(v, _) if v.len() == 2 => {
                let (open, close) = if op == "of" { ("(", ")") } else { ("[", "]") };
                format!("{}{open}{}{close}", gprint(&v[0]), gprint(&v[1]))
            }
            _ => format!("{}({})", pname(op), gprint(f)),
        },
        // printasfactorial (:7721).
        "factorial" => {
            if matches!(f, Gen::Idnt(_)) || (gtype(f) <= ty::DOUBLE && is_positive(f)) {
                format!("{}!", gprint(f))
            } else {
                format!("({})!", gprint(f))
            }
        }
        // printassqrt (:1396).
        _ => format!("{}({})", pname(op), gprint(f)),
    }
}
