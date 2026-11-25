//! `gen2mathml` (mathml.cc :1879-1975) and `symbolic2mathml`
//! (:1571-1773), with the helpers they call, one function per giac
//! function, each citing its lines. Every part's text is checked against
//! the cap as it is made; the tree comes flattened (`flat`).

use super::print::{
    equals, gprint, gtype, is_int, is_one, is_positive, is_zero, need_parenthesis, neg_num, pname,
    print_atom, ty,
};
use crate::markup::tree::{Gen, Sub};

/// The text grew past the length allowed.
#[derive(Debug)]
pub(super) struct TooLong;

type Text = Result<String, TooLong>;

/// The heads with a `texprint` in giac (the `define_unary_function_eval*`
/// of its sources with one): `symbolic2mathml` sends them to
/// `mathml_print` (:1603).
const TEXPRINT: &[&str] = &[
    "MOD",
    "NTHROOT",
    "ampersand_times",
    "and",
    "at",
    "binary_minus",
    "complement",
    "compose",
    "composepow",
    "conj",
    "cross",
    "decrement",
    "derive",
    "deuxpoints",
    "div",
    "divcrement",
    "division",
    "dollar",
    "equal",
    "equal2",
    "et",
    "exp",
    "im",
    "increment",
    "inferieur_egal",
    "inferieur_strict",
    "integrate",
    "intersect",
    "interval",
    "leftopen_interval",
    "leftrightopen_interval",
    "lim",
    "limit",
    "minus",
    "moyal",
    "multcrement",
    "normalmod",
    "of",
    "ou",
    "oufr",
    "pointdivision",
    "pointpow",
    "pointprod",
    "polar_complex",
    "program",
    "quote_pow",
    "re",
    "rightopen_interval",
    "same",
    "sqrt",
    "sto",
    "struct_dot",
    "superieur_egal",
    "superieur_strict",
    "symmetric_difference",
    "ti_semi",
    "tilocal",
    "union",
    "xor",
];

/// `and`/`or` and the set operators, with their sign (:1592-1593).
const ANDOR: [(&str, &str); 7] = [
    ("and", "&and;"),
    ("et", "&and;"),
    ("ou", "&or;"),
    ("xor", "&veebar;"),
    ("union", "&cup;"),
    ("intersect", "&cap;"),
    ("symmetric_difference", "&Delta;"),
];

const TIMES: &str = "<mo>&times;</mo>";

/// The converter, with the longest text allowed.
#[derive(Debug)]
pub(super) struct Writer {
    max_len: usize,
}

impl Writer {
    pub(super) fn new(max_len: usize) -> Writer {
        Writer { max_len }
    }

    /// `s`, refused when longer than the cap.
    pub(super) fn checked(&self, s: String) -> Text {
        self.fits(&s)?;
        Ok(s)
    }

    fn fits(&self, s: &str) -> Result<(), TooLong> {
        if s.len() > self.max_len {
            Err(TooLong)
        } else {
            Ok(())
        }
    }

    /// `gen2mathml(e)` (:1879), its text checked.
    pub(super) fn node(&self, e: &Gen) -> Text {
        let s = self.gen2mathml(e)?;
        self.checked(s)
    }

    /// `gen2mathml` (:1879-1975).
    fn gen2mathml(&self, e: &Gen) -> Text {
        Ok(match e {
            Gen::Symb(op, _) if op == "sum" => self.printassum(e)?,
            // :607 `mathml_printasabs`.
            Gen::Symb(op, f) if op == "abs" => {
                format!(
                    "<mo>&VerticalBar;</mo>{}<mo>&VerticalBar;</mo>",
                    self.node(f)?
                )
            }
            // :1888-1898: `<mn>` and giac's print.
            Gen::Int(_) | Gen::Bool(_) | Gen::Dbl { .. } => format!("<mn>{}</mn>", print_atom(e)),
            Gen::Cplx(re, im) => cplx(re, im),
            // :1920-1927.
            Gen::Idnt(name) => match name.as_str() {
                "infinity" => "<mn>&infin;</mn>".to_owned(),
                "pi" => "<mi>&pi;</mi>".to_owned(),
                "undef" => "<mi>undef</mi>".to_owned(),
                _ => idnt2mathml(name),
            },
            Gen::Symb(..) => self.symbolic(e)?,
            // :1934-1950.
            Gen::Vect(items, sub) if *sub != Sub::Seq && ckmatrix(items) => self.matrix(items)?,
            Gen::Vect(items, sub) => self.vect(items, *sub)?,
            // :1956-1958.
            Gen::Frac(a, b) => format!(
                "<mfrac><mrow>{}</mrow><mrow>{}</mrow></mfrac>",
                self.node(a)?,
                self.node(b)?
            ),
            Gen::Str(s) => string2mathml(s),
            // :1964-1965.
            Gen::Func(_) => format!("<mi>{}</mi>", print_atom(e)),
        })
    }

    /// `matrix2mathml` (:200).
    fn matrix(&self, rows: &[Gen]) -> Text {
        let mut s = String::from("<mrow><mo>(</mo><mtable>");
        for r in rows {
            s.push_str("<mtr>");
            for x in r.args_of_vect() {
                s.push_str("<mtd>");
                s.push_str(&self.node(x)?);
                s.push_str("</mtd>");
                self.fits(&s)?;
            }
            s.push_str("</mtr>");
        }
        s.push_str("</mtable><mo>)</mo></mrow>");
        Ok(s)
    }

    /// `_VECT2mathml` (:221), the brackets `begin_VECT_string` gives
    /// (gen.cc :13218).
    fn vect(&self, items: &[Gen], sub: Sub) -> Text {
        let (open, close) = match sub {
            Sub::Seq => ("", ""),
            Sub::Set => ("set[", "]"),
            Sub::List => ("[", "]"),
        };
        let mut s = String::from(open);
        for (k, x) in items.iter().enumerate() {
            if k > 0 {
                s.push_str("<mo>,</mo>");
            }
            s.push_str(&self.node(x)?);
            self.fits(&s)?;
        }
        s.push_str(close);
        Ok(s)
    }

    /// `prod_vect2mathml` (:274) when `bra`, `prod_vect2mathml_no_bra`
    /// (:321) else: the factors side by side, `&times;` after a number
    /// before a number (`needs_times` :256) and before a positive integer's
    /// power.
    fn prod_vect(&self, v: &[Gen], bra: bool) -> Text {
        match v {
            [] => return Ok("<mn>1</mn>".to_owned()),
            [x] if !bra => return self.node(x),
            _ => {}
        }
        let mut s = String::new();
        let mut prec_int = false;
        for (k, it) in v.iter().enumerate() {
            match it {
                Gen::Symb(op, _) if op == "plus" || op == "neg" => {
                    s.push_str(" <mo>(</mo>");
                    s.push_str(&self.node(it)?);
                    s.push_str(" <mo>)</mo>");
                }
                Gen::Cplx(re, im) if !is_zero(re) && !is_zero(im) => {
                    s.push_str("<mo>(</mo>");
                    s.push_str(&self.node(it)?);
                    s.push_str("<mo>)</mo>");
                }
                _ => {
                    // :292-299 (the `&times;` not doubled), :338-345.
                    if prec_int && pow_times(it) && !(bra && s.ends_with(TIMES)) {
                        s.push_str(TIMES);
                    }
                    s.push_str(&self.node(it)?);
                }
            }
            self.fits(&s)?;
            prec_int = gtype(it) <= ty::REAL;
            if let Some(next) = v.get(k + 1)
                && prec_int
                && gtype(next) <= ty::REAL
            {
                s.push_str(TIMES);
            }
        }
        Ok(s)
    }

    /// `mathml_printsommetasoperator` (:369): the arguments with `sop`
    /// between them.
    fn sommet_as_operator(&self, f: &Gen, sop: &str) -> Text {
        let Gen::Vect(v, _) = f else {
            return Ok(format!("<mi>{}</mi>", gprint(f)));
        };
        let Some((a, rest)) = v.split_first() else {
            return Ok(String::new());
        };
        let eq = sop == "<mo>=</mo>" || sop == "<mo>≈</mo>";
        let mut s = if matches!(a, Gen::Frac(..)) || !need_parenthesis(a) || a.is_sym(&["of"]) || eq
        {
            self.node(a)?
        } else {
            format!("<mo>(</mo>{}<mo>)</mo>", self.node(a)?)
        };
        for b in rest {
            s.push_str(sop);
            if matches!(b, Gen::Symb(..)) && !eq {
                s.push_str("<mo>(</mo>");
                s.push_str(&self.node(b)?);
                s.push_str("<mo>)</mo>");
            } else {
                s.push_str(&self.node(b)?);
            }
            self.fits(&s)?;
        }
        Ok(s)
    }

    /// `mathml_equal2arrow` (:448): `x = a` as `x &rarr; a`.
    fn equal2arrow(&self, g: &Gen) -> Text {
        match g {
            Gen::Symb(op, _) if op == "equal" => match g.args() {
                [a, b, ..] => Ok(format!(
                    "<mrow>{}<mo>&rarr;</mo>{}</mrow>",
                    self.node(a)?,
                    self.node(b)?
                )),
                _ => self.node(g),
            },
            _ => self.node(g),
        }
    }

    /// `mathml_printaslimit` (:457).
    #[allow(clippy::many_single_char_names)]
    fn limit(&self, f: &Gen) -> Text {
        const LIM: &str = "<mo>lim</mo>";
        let Gen::Vect(v, _) = f else {
            return Ok(format!("{LIM}{}", self.node(f)?));
        };
        Ok(match v.as_slice() {
            [a] => format!("{LIM}{}", self.node(a)?),
            [a, x] => format!(
                "<munder>{LIM}{}</munder>{}",
                self.equal2arrow(x)?,
                self.node(a)?
            ),
            [a, x, d] => {
                let side = if is_one(d) {
                    Some("+")
                } else if equals(d, -1, 1) {
                    Some("-")
                } else {
                    None
                };
                match side {
                    Some(side) => format!(
                        "<munder>{LIM}<mrow><msup>{}<mo>{side}</mo></msup></mrow></munder>{}",
                        self.equal2arrow(x)?,
                        self.node(a)?
                    ),
                    None => format!(
                        "<munder>{LIM}{}</munder>{}",
                        self.equal2arrow(x)?,
                        self.node(a)?
                    ),
                }
            }
            _ => LIM.to_owned(),
        })
    }

    /// `mathml_printasintegrate` (:485): bounds and variable printed.
    fn integrate(&self, f: &Gen) -> Text {
        const INT: &str = "<mo>&int;</mo>";
        let Gen::Vect(v, _) = f else {
            return Ok(format!("{INT}{}", self.node(f)?));
        };
        Ok(match v.as_slice() {
            [a] => format!("{INT}{}", self.node(a)?),
            [a, x] => format!("{INT}{}<mi>d</mi>{}", self.node(a)?, self.node(x)?),
            [a, x, lo, hi] => format!(
                "<msubsup>{INT}<mi>{}</mi><mi>{}</mi></msubsup>{}<mi>d</mi><mi>{}</mi>",
                gprint(lo),
                gprint(hi),
                self.node(a)?,
                gprint(x)
            ),
            _ => INT.to_owned(),
        })
    }

    /// `mathml_printasderive` (:509).
    fn derive(&self, f: &Gen) -> Text {
        let v = match f {
            Gen::Vect(v, _) if v.len() >= 2 => v,
            _ => {
                return Ok(
                    if matches!(f, Gen::Idnt(_)) || gtype(f) <= ty::CPLX || f.is_sym(&["derive"]) {
                        format!("{}'", self.node(f)?)
                    } else {
                        format!(
                            "<msup><mrow><mo>(</mo>{}<mo>)</mo></mrow><mi>'</mi></msup>",
                            self.node(f)?
                        )
                    },
                );
            }
        };
        let (lp, rp) = if matches!(v[0], Gen::Symb(..)) {
            ("<mo>(</mo>", "<mo>)</mo>")
        } else {
            ("", "")
        };
        if let Some(n) = v.get(2) {
            let n = self.node(n)?;
            return Ok(format!(
                "<mfrac><mrow><msup><mo>&part;</mo><mrow>{n}</mrow></msup>{lp}{}{rp}</mrow>\
                 <mrow><mo>&part;</mo><msup>{}<mrow>{n}</mrow></msup></mrow></mfrac>",
                self.node(&v[0])?,
                self.node(&v[1])?
            ));
        }
        Ok(format!(
            "<mfrac><mrow><mo>&part;</mo>{lp}{}{rp}</mrow><mrow><mo>&part;</mo>{}</mrow></mfrac>",
            self.node(&v[0])?,
            self.node(&v[1])?
        ))
    }

    /// `mathml_printassto` (:433): `x := value`.
    fn sto(&self, f: &Gen) -> Text {
        let Gen::Vect(v, _) = f else {
            return Ok(format!("<mi>:=</mi><mo>(</mo>{}<mo>)</mo>", self.node(f)?));
        };
        let [value, name] = v.as_slice() else {
            return Ok(format!("<mi>:=</mi><mo>(</mo>{}<mo>)</mo>", self.node(f)?));
        };
        // giac compares the value's type to `_SEQ__VECT` (1), which is
        // `_DOUBLE_`'s number: a double gets brackets.
        let s = format!("{}<mi>:=</mi>", self.node(name)?);
        Ok(if gtype(value) == ty::DOUBLE {
            format!("{s}<mo>(</mo>{}<mo>)</mo>", self.node(value)?)
        } else {
            format!("{s}{}", self.node(value)?)
        })
    }

    /// `mathml_print` (:526-598): the heads with a `texprint`.
    fn mathml_print(&self, e: &Gen) -> Text {
        let Gen::Symb(op, f) = e else {
            unreachable!("a symbolic");
        };
        let sop = match op.as_str() {
            // :530-566: `≈` when the right side is a double.
            "equal" | "equal2" => {
                let approx = match e.args() {
                    [_, rhs, ..] => match rhs {
                        Gen::Dbl { .. } => true,
                        Gen::Cplx(a, b) => {
                            matches!(**a, Gen::Dbl { .. }) || matches!(**b, Gen::Dbl { .. })
                        }
                        _ => false,
                    },
                    _ => false,
                };
                Some(if approx { "<mo>≈</mo>" } else { "<mo>=</mo>" })
            }
            "different" => Some("<mo>≠</mo>"),
            "inferieur_egal" => Some("<mo>&le;</mo>"),
            "superieur_egal" => Some("<mo>&ge;</mo>"),
            "inferieur_strict" => Some("<mo>&lt;</mo>"),
            "superieur_strict" => Some("<mo>&gt;</mo>"),
            _ => None,
        };
        if let Some(sop) = sop {
            return self.sommet_as_operator(f, sop);
        }
        Ok(match op.as_str() {
            // :418-424, :403, :407.
            "re" => format!(
                "<mi>Re</mi><mrow><mo>(</mo>{}<mo>)</mo></mrow>",
                self.node(f)?
            ),
            "im" => format!(
                "<mi>Im</mi><mrow><mo>(</mo>{}<mo>)</mo></mrow>",
                self.node(f)?
            ),
            "exp" => format!("<msup><mi>e</mi><mrow>{}</mrow></msup>", self.node(f)?),
            "sqrt" => format!("<msqrt>{}</msqrt>", self.node(f)?),
            "integrate" => self.integrate(f)?,
            "derive" => self.derive(f)?,
            "sto" => self.sto(f)?,
            "limit" => self.limit(f)?,
            // :594-597: giac's print, `<` and `>` escaped, newlines spaces.
            _ => format!(
                "<mtext>{}</mtext>",
                gprint(e)
                    .replace('<', "&lt;")
                    .replace('>', "&gt;")
                    .replace('\n', " ")
            ),
        })
    }

    /// `mathml_printassum` (:600).
    fn printassum(&self, e: &Gen) -> Text {
        let [body, k, lo, hi] = e.args() else {
            return Ok(format!("<mtext>{}</mtext>", gprint(e)));
        };
        Ok(format!(
            "<munderover><mi>&Sigma;</mi><mrow>{}<mo>=</mo>{}</mrow><mrow>{}</mrow></munderover>{}",
            self.node(k)?,
            self.node(lo)?,
            self.node(hi)?,
            self.node(body)?
        ))
    }

    /// `symbolic2mathml` (:1571-1773).
    #[allow(clippy::too_many_lines)]
    fn symbolic(&self, e: &Gen) -> Text {
        let Gen::Symb(op, f) = e else {
            unreachable!("a symbolic");
        };
        let op = op.as_str();
        let opstring = pname(op);
        // :1575-1582.
        if op == "of" || op == "at" {
            let (a, b) = match e.args() {
                [a, b, ..] => (a, b),
                _ => (f.as_ref(), f.as_ref()),
            };
            return Ok(if op == "of" {
                format!("{}<mo>(</mo>{}<mo>)</mo>", self.node(a)?, self.node(b)?)
            } else {
                format!(
                    "<msub><mrow>{}</mrow><mrow>{}</mrow></msub>",
                    self.node(a)?,
                    self.node(b)?
                )
            });
        }
        // :1583-1589: a program without locals, `x → body`.
        if op == "program"
            && let Gen::Vect(v, _) = &**f
            && v.len() >= 3
            && !v[2].is_sym(&["local"])
        {
            return Ok(format!(
                "{}<mo>→</mo>{}",
                self.node(&v[0])?,
                self.node(&v[2])?
            ));
        }
        if op == "complement" {
            return Ok(format!("{}<mo>&#x2201;</mo>", self.node(f)?));
        }
        let andor = ANDOR.iter().find(|(o, _)| *o == op).map(|(_, m)| *m);
        // :1603-1604.
        if andor.is_none()
            && opstring != "/"
            && (TEXPRINT.binary_search(&op).is_ok() || op == "different")
        {
            return self.mathml_print(e);
        }
        // :1609-1610.
        if op == "member"
            && let [a, b] = e.args()
        {
            return Ok(format!(
                "{}<mtext>∈</mtext>{}",
                self.node(a)?,
                self.node(b)?
            ));
        }
        // :1611-1623.
        if op == "vector" {
            let mut s = String::from("<mover><mrow accent=\"true\">");
            for x in e.args() {
                s.push_str(&self.node(x)?);
                self.fits(&s)?;
            }
            s.push_str("</mrow><mo mathsize=\"100%\">&rarr;</mo></mover>");
            return Ok(s);
        }
        // :1629-1630.
        if let Gen::Vect(v, _) = &**f
            && v.is_empty()
        {
            return Ok(format!("<mi>{opstring}()</mi>"));
        }
        // :1631-1654: one argument (or a first argument that is a vector).
        let v = match &**f {
            Gen::Vect(v, _) if !matches!(v[0], Gen::Vect(..)) => v,
            _ => return self.unary(op, f),
        };
        // :1657-1673.
        if let Some(sign) = andor {
            let mut s = String::new();
            for (i, x) in v.iter().enumerate() {
                let par = matches!(x, Gen::Symb(xop, _) if xop != "complement" && xop != op);
                s.push_str("<mrow>");
                if par {
                    s.push('(');
                }
                s.push_str(&self.node(x)?);
                if par {
                    s.push(')');
                }
                s.push_str("</mrow>");
                if i + 1 != v.len() {
                    s.push_str("<mo> ");
                    s.push_str(sign);
                    s.push_str(" </mo>");
                }
                self.fits(&s)?;
            }
            return Ok(s);
        }
        match op {
            "plus" => self.plus(v),
            "prod" => self.prod(v),
            "pow" => self.pow(v),
            // :1761-1764.
            "division" if v.len() >= 2 => Ok(format!(
                "<mfrac><mrow>{}</mrow><mrow>{}</mrow></mfrac>",
                self.node(&v[0])?,
                self.node(&v[1])?
            )),
            // :1766-1772.
            _ => {
                let mut s = format!("<mi>{opstring}</mi><mrow><mo>(</mo>");
                for (i, x) in v.iter().enumerate() {
                    if i > 0 {
                        s.push(',');
                    }
                    s.push_str(&self.node(x)?);
                    self.fits(&s)?;
                }
                s.push_str("<mo>)</mo></mrow>");
                Ok(s)
            }
        }
    }

    /// :1631-1654: a head on one argument.
    fn unary(&self, op: &str, f: &Gen) -> Text {
        let opstring = pname(op);
        Ok(match op {
            "factorial" if matches!(f, Gen::Symb(..)) => format!(
                "<mrow><mo>(</mo>{}<mo>)</mo></mrow><mtext>!</mtext>",
                self.node(f)?
            ),
            "factorial" => format!("{}<mtext>!</mtext>", self.node(f)?),
            "neg" | "plus" => {
                let bare = match f {
                    Gen::Symb(fop, _) => {
                        (matches!(fop.as_str(), "at" | "of")
                            && matches!(f.args().first(), Some(Gen::Idnt(_))))
                            || fop == "inv"
                            || fop == "pow"
                    }
                    _ => true,
                };
                if bare {
                    format!("<mo>{opstring}</mo>{}", self.node(f)?)
                } else {
                    format!(
                        "<mo>{opstring}</mo><mrow><mo>(</mo>{}<mo>)</mo></mrow>",
                        self.node(f)?
                    )
                }
            }
            "inv" => format!(
                "<mfrac><mrow><mn>1</mn></mrow><mrow>{}</mrow></mfrac>",
                self.node(f)?
            ),
            // A first argument that is a vector: giac reads two.
            "pow" if f.args().len() >= 2 => format!(
                "<msup><mrow>{}</mrow><mrow>{}</mrow></msup>",
                self.node(&f.args()[0])?,
                self.node(&f.args()[1])?
            ),
            _ => format!(
                "<mi>{opstring}</mi><mrow><mo>(</mo>{}<mo>)</mo></mrow>",
                self.node(f)?
            ),
        })
    }

    /// :1674-1697: the terms, `-` for a `neg` or a negative integer, `+`
    /// before the others but one whose text starts with `<mo>-`.
    fn plus(&self, v: &[Gen]) -> Text {
        let mut s = String::new();
        for (i, x) in v.iter().enumerate() {
            match x {
                Gen::Symb(op, a) if op == "neg" => {
                    if a.is_sym(&["plus"]) {
                        s.push_str("<mo>-</mo><mrow><mo>(</mo>");
                        s.push_str(&self.node(a)?);
                        s.push_str("<mo>)</mo></mrow>");
                    } else {
                        s.push_str("<mo>-</mo><mrow>");
                        s.push_str(&self.node(a)?);
                        s.push_str("</mrow>");
                    }
                }
                _ if is_int(x) && !is_positive(x) => {
                    s.push_str("<mo>-</mo><mn>");
                    s.push_str(&gprint(&neg_num(x)));
                    s.push_str("</mn>");
                }
                _ => {
                    let adds = self.node(x)?;
                    if i > 0 && !adds.starts_with("<mo>-") {
                        s.push_str("<mo>+</mo>");
                    }
                    s.push_str(&adds);
                }
            }
            self.fits(&s)?;
        }
        Ok(s)
    }

    /// :1698-1722: the factors, an `inv` one in a denominator; a factor
    /// after a denominator closes the fraction with `&times;`.
    fn prod(&self, v: &[Gen]) -> Text {
        let mut s = String::new();
        let (mut num, mut den) = (Vec::new(), Vec::new());
        for x in v {
            if let Gen::Symb(op, d) = x
                && op == "inv"
            {
                den.push((**d).clone());
                continue;
            }
            if !den.is_empty() {
                s.push_str(&self.fraction(&num, &den)?);
                s.push_str(TIMES);
                self.fits(&s)?;
                num.clear();
                den.clear();
            }
            num.push(x.clone());
        }
        if den.is_empty() {
            s.push_str(&self.prod_vect(&num, true)?);
        } else {
            s.push_str(&self.fraction(&num, &den)?);
        }
        Ok(s)
    }

    fn fraction(&self, num: &[Gen], den: &[Gen]) -> Text {
        Ok(format!(
            "<mfrac><mrow>{}</mrow><mrow>{}</mrow></mfrac>",
            self.prod_vect(num, false)?,
            self.prod_vect(den, false)?
        ))
    }

    /// :1723-1759: `<msqrt>` for the power 1/2, `1/<msqrt>` for -1/2, else
    /// `<msup>` with the base bracketed (`s_bra`) when it is a symbolic of
    /// several arguments, a fraction, a complex or a negative integer.
    fn pow(&self, v: &[Gen]) -> Text {
        let (Some(base), Some(ex)) = (v.first(), v.last()) else {
            unreachable!("not empty");
        };
        if equals(ex, 1, 2) {
            return Ok(format!("<msqrt><mrow>{}</mrow></msqrt>", self.node(base)?));
        }
        if equals(ex, -1, 2) {
            return Ok(format!(
                "<mfrac><mn>1</mn><msqrt>{}</msqrt></mfrac>",
                self.node(base)?
            ));
        }
        let s0 = self.node(base)?;
        let s1 = self.node(v.get(1).unwrap_or(ex))?;
        let bra = match base {
            Gen::Symb(_, mf) => match &**mf {
                Gen::Vect(items, _) => items.len() != 1 || matches!(items[0], Gen::Vect(..)),
                _ => true,
            },
            Gen::Frac(..) => true,
            Gen::Cplx(_, im) => !is_zero(im),
            Gen::Int(_) | Gen::Bool(_) => !is_positive(base),
            _ => false,
        };
        Ok(if bra {
            format!(
                "<msup><mfenced open=\"(\" close=\")\"><mrow>({s0})</mrow></mfenced><mrow>{s1}</mrow></msup>"
            )
        } else {
            format!("<msup><mrow> {s0}</mrow><mrow>{s1}</mrow></msup>")
        })
    }
}

impl Gen {
    /// A vector's items (a matrix row's).
    fn args_of_vect(&self) -> &[Gen] {
        match self {
            Gen::Vect(items, _) => items,
            _ => &[],
        }
    }
}

/// :292-299, :338-345: a power of a positive integer, not a square root,
/// takes `&times;` after a number.
fn pow_times(it: &Gen) -> bool {
    match it {
        Gen::Symb(op, _) if op == "pow" => {
            let a = it.args();
            a.last().is_some_and(|ex| !equals(ex, 1, 2))
                && a.first().is_some_and(|b| is_int(b) && is_positive(b))
        }
        _ => false,
    }
}

/// giac's `ckmatrix`: rows of one length, at least one column, no vector
/// in a cell.
fn ckmatrix(items: &[Gen]) -> bool {
    let Some(Gen::Vect(first, _)) = items.first() else {
        return false;
    };
    !first.is_empty()
        && items.iter().all(|r| {
            matches!(r, Gen::Vect(row, _) if row.len() == first.len()
                && !row.iter().any(|x| matches!(x, Gen::Vect(..))))
        })
}

/// :1899-1919: a complex number by its parts, `<mn>1</mn><mo>+</mo><mn>2</mn><mi>i</mi>`.
fn cplx(re: &Gen, im: &Gen) -> String {
    let mut part_re = String::new();
    let mut part_im = "<mi>i</mi>".to_owned();
    if !is_zero(re) {
        part_re = format!("<mn>{}</mn>", gprint(re));
    }
    if !is_zero(im) {
        if is_positive(im) {
            if !is_zero(re) {
                part_re.push_str("<mo>+</mo>");
            }
        } else {
            part_re.push_str("<mo>-</mo>");
        }
    }
    if !is_one(&neg_num(im)) && !is_one(im) {
        let abs = if is_positive(im) {
            im.clone()
        } else {
            neg_num(im)
        };
        part_im = format!("<mn>{}</mn>{part_im}", gprint(&abs));
    }
    if is_zero(im) {
        part_im.clear();
        if is_zero(re) {
            "<mn>0.0</mn>".clone_into(&mut part_re);
        }
    }
    part_re + &part_im
}

/// `string2mathml` (:111): `&`, `<`, `>` escaped, lines as rows of a
/// one-column `<mtable>`.
fn string2mathml(m: &str) -> String {
    let s = m
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('\n', "</mi></mtd></mtr><mtr><mtd><mi>");
    format!("<mtable columnalign=\"left\"><mtr><mtd><mi>{s}</mi></mtd></mtr></mtable>")
}

/// `idnt2mathml_` (:1808): a greek letter's or a set's name as its entity,
/// a non-letter tail after it kept (`alpha1`, `&alpha;1`).
fn idnt2mathml_(s0: &str) -> String {
    let b = s0.as_bytes();
    let n = b.len();
    if n == 0 {
        return String::new();
    }
    // giac reads bytes as signed chars: a non-ASCII byte is no letter.
    let mut j = n - 1;
    while j >= 2 && !b[j].is_ascii_alphabetic() {
        j -= 1;
    }
    if !s0.is_char_boundary(j + 1) {
        return s0.to_owned();
    }
    let (s, sadd) = s0.split_at(j + 1);
    let entity = match s {
        "mu" | "nu" | "pi" | "xi" | "Xi" | "chi" | "phi" | "Phi" | "eta" | "rho" | "tau"
        | "psi" | "Psi" | "beta" | "zeta" | "alpha" | "delta" | "Delta" | "gamma" | "Gamma"
        | "kappa" | "theta" | "Theta" | "sigma" | "Sigma" | "Omega" | "omega" | "aleph"
        | "nabla" | "lambda" | "Lambda" | "approx" | "epsilon" | "product" => format!("&{s};"),
        "im" => "&Im;".to_owned(),
        "re" => "&Re;".to_owned(),
        "RR" => "ℝ".to_owned(),
        "QQ" => "ℚ".to_owned(),
        "ZZ" => "ℤ".to_owned(),
        "CC" => "ℂ".to_owned(),
        "NN" => "ℕ".to_owned(),
        "HH" => "ℍ".to_owned(),
        "PP" => "ℙ".to_owned(),
        _ => return s0.to_owned(),
    };
    entity + sadd
}

/// `idnt2mathml` (:1870): `x_1` as `<msub>`, at its last `_` not at the
/// end.
fn idnt2mathml(s0: &str) -> String {
    let b = s0.as_bytes();
    let n = b.len();
    for j in (1..n.saturating_sub(1)).rev() {
        if b[j] == b'_' {
            return format!(
                "<msub><mi>{}</mi><mn>{}</mn></msub>",
                idnt2mathml_(&s0[..j]),
                &s0[j + 1..]
            );
        }
    }
    format!("<mi>{}</mi>", idnt2mathml_(s0))
}
