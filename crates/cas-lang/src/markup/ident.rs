//! Names and strings in LaTeX and `MathML`, giac's way: identifiers
//! (markup.cc `idnt2markup` :621-934, the non-unit arms), function names
//! (`func2markup` :1106-1226) and strings (`string2markup` :304-468).

use super::mml;

/// giac's `is_greek_letter` (:141). No `iota`, `upsilon` or `omicron`.
pub(super) fn is_greek_letter(s: &str) -> bool {
    matches!(
        s,
        "mu" | "nu"
            | "pi"
            | "Pi"
            | "xi"
            | "Xi"
            | "chi"
            | "phi"
            | "Phi"
            | "eta"
            | "rho"
            | "tau"
            | "psi"
            | "Psi"
            | "beta"
            | "zeta"
            | "alpha"
            | "delta"
            | "Delta"
            | "gamma"
            | "Gamma"
            | "kappa"
            | "theta"
            | "Theta"
            | "sigma"
            | "Sigma"
            | "Omega"
            | "omega"
            | "lambda"
            | "Lambda"
            | "epsilon"
            | "varepsilon"
    )
}

/// giac's `my_isalpha`: an ASCII letter or any non-ASCII character.
pub(super) fn my_isalpha(c: char) -> bool {
    c.is_ascii_alphabetic() || !c.is_ascii()
}

/// giac's `is_double_letter` (:167): `xx`.
fn is_double_letter(s: &str) -> bool {
    let mut chars = s.chars();
    match (chars.next(), chars.next(), chars.next()) {
        (Some(a), Some(b), None) => my_isalpha(a) && a == b,
        _ => false,
    }
}

/// giac's `is_digits` (:171): true for "".
fn is_digits(s: &str) -> bool {
    s.chars().all(|c| c.is_ascii_digit())
}

fn char_len(s: &str) -> usize {
    s.chars().count()
}

/// C's `atof`: the value of the longest prefix that reads as a decimal
/// (after white space), `inf` or `nan`; 0 when there is none.
pub(super) fn atof(s: &str) -> f64 {
    let s = s.trim_start();
    let b = s.as_bytes();
    let mut i = usize::from(matches!(b.first(), Some(b'+' | b'-')));
    let lower = s[i..].to_ascii_lowercase();
    for word in ["infinity", "inf", "nan"] {
        if lower.starts_with(word) {
            return s[..i + word.len()].parse().unwrap_or(0.0);
        }
    }
    let digits = |i: &mut usize| {
        let start = *i;
        while b.get(*i).is_some_and(u8::is_ascii_digit) {
            *i += 1;
        }
        *i > start
    };
    let mut any = digits(&mut i);
    if b.get(i) == Some(&b'.') {
        let mut j = i + 1;
        if digits(&mut j) || any {
            any = true;
            i = j;
        }
    }
    if !any {
        return 0.0;
    }
    if matches!(b.get(i), Some(b'e' | b'E')) {
        let mut j = i + 1;
        if matches!(b.get(j), Some(b'+' | b'-')) {
            j += 1;
        }
        if digits(&mut j) {
            i = j;
        }
    }
    s[..i].parse().unwrap_or(0.0)
}

/// `phi` and `epsilon` are written `\varphi`, `\varepsilon`.
fn greek(s: &str) -> String {
    if matches!(s, "phi" | "epsilon") {
        format!("\\var{s}")
    } else {
        format!("\\{s}")
    }
}

/// A greek letter as a `MathML` entity: `&alpha;` (no `var` forms, :875).
fn greek_mml(s: &str) -> String {
    format!("&{s};")
}

/// An identifier in LaTeX, and whether it is a name of more than one
/// letter (giac's `isname`), :621-934: `x12` is `x_{12}`, `alpha1`
/// `\alpha_{1}`, `abc_1` `\mathrm{abc}_{1}`, `a_bc` `a_{\mathrm{bc}}`,
/// `x_1_2` `x_{1,2}`, `xx` `\mathbf{x}`, `abc` `\mathrm{abc}`.
pub(super) fn idnt2markup(s_orig: &str) -> (String, bool) {
    idnt(s_orig, None)
}

/// The same identifier in `MathML` presentation (`typ` 2), its outer tag
/// with id `idc`: `<mi>x</mi>`, `<msub><mi>x</mi><mrow>1</mrow></msub>`,
/// `<mi>&alpha;</mi>`, `<mi mathvariant='bold'>x</mi>` for `xx` (a name
/// here, unlike in LaTeX). giac's slips kept: `alpha_` is `&alpha;_`
/// untagged, indices are joined by an unclosed `<mo>,<mo>`.
pub(super) fn idnt2mml(s_orig: &str, idc: u32) -> (String, bool) {
    idnt(s_orig, Some(idc))
}

/// :621-934: LaTeX when `mml` is `None`, else `MathML` with that id.
fn idnt(s_orig: &str, mml: Option<u32>) -> (String, bool) {
    // :845: `greek_`; `phi_` is `varphi` in both.
    if let Some(stem) = s_orig.strip_suffix('_')
        && is_greek_letter(stem)
    {
        let text = match mml {
            None => format!("{}\\_", greek(stem)),
            Some(_) if matches!(stem, "phi" | "epsilon") => format!("&var{stem};_"),
            Some(_) => format!("{}_", greek_mml(stem)),
        };
        return (text, false);
    }
    // Trailing digits are an index after `log`, one letter, a greek letter
    // or a double letter.
    let stem_end = s_orig.trim_end_matches(|c: char| c.is_ascii_digit()).len();
    let mut s = &s_orig[..stem_end];
    let mut ssub = "";
    if stem_end < s_orig.len()
        && (s == "log" || char_len(s) == 1 || is_greek_letter(s) || is_double_letter(s))
    {
        ssub = &s_orig[stem_end..];
    }
    // Else the first `_` splits.
    if ssub.is_empty()
        && let Some(pos) = s_orig.find('_')
    {
        s = &s_orig[..pos];
        let au = &s_orig[pos + 1..];
        if char_len(s) == 1 || is_greek_letter(s) || is_double_letter(s) || is_digits(au) {
            ssub = au;
        }
    }
    if ssub.is_empty() {
        s = s_orig;
    }
    let mut len = char_len(s);
    let mut s = s.to_owned();
    let greek_of = |s: &str| match mml {
        None => greek(s),
        Some(_) => greek_mml(s),
    };
    if is_greek_letter(&s) {
        s = greek_of(&s);
        len = 1;
    }
    let mut ssub = ssub.to_owned();
    let len_sub = if !ssub.is_empty() && is_greek_letter(&ssub) {
        ssub = greek_of(&ssub);
        1
    } else {
        char_len(&ssub)
    };
    // :885: LaTeX `\mathbf{x}`, one letter; MathML a bold name of two.
    let mut bold = false;
    if is_double_letter(&s) {
        let first = s.chars().next().expect("two letters");
        if mml.is_some() {
            s = first.to_string();
            bold = true;
        } else {
            s = format!("\\mathbf{{{first}}}");
            len = 1;
        }
    }
    if len > 3 && s.starts_with("` ") && s.ends_with('`') {
        s = s.chars().skip(2).take(len - 3).collect();
    }
    if let Some(idc) = mml {
        return (idnt_mml(&s, len, bold, &ssub, idc), len > 1);
    }
    let mut ret = if len == 1 {
        s
    } else {
        format!("\\mathrm{{{s}}}")
    };
    if ssub.is_empty() {
        ret = ret.replace('_', "\\_");
    } else {
        let tokens: Vec<&str> = ssub.split('_').collect();
        if tokens
            .iter()
            .all(|t| !t.is_empty() && (*t == "0" || atof(t) != 0.0))
        {
            ssub = tokens.join(",");
        }
        if len_sub == 1 || atof(&ssub) != 0.0 {
            ret = format!("{ret}_{{{ssub}}}");
        } else {
            ret = format!("{ret}_{{\\mathrm{{{ssub}}}}}");
        }
    }
    (ret, len > 1)
}

/// :897-921, `MathML`: `<mi>` (bold for a double letter), the id on it, or
/// on the `<msub>` when there is an index: `<mrow>1<mo>,<mo>2</mrow>` for
/// numbers, else `<mn>` or `<mi>`.
fn idnt_mml(s: &str, len: usize, bold: bool, ssub: &str, idc: u32) -> String {
    let mi_id = if ssub.is_empty() { idc } else { 0 };
    let variant = if bold && len > 1 { "mathvariant" } else { "" };
    let ret = mml::tag("mi", s, mi_id, &[(variant, "bold")]);
    if ssub.is_empty() {
        return ret;
    }
    let tokens: Vec<&str> = ssub.split('_').collect();
    let index = if tokens
        .iter()
        .all(|t| !t.is_empty() && (*t == "0" || atof(t) != 0.0))
    {
        format!("<mrow>{}</mrow>", tokens.join("<mo>,<mo>"))
    } else if ssub == "0" || atof(ssub) != 0.0 {
        format!("<mn>{ssub}</mn>")
    } else {
        format!("<mi>{ssub}</mi>")
    };
    mml::tag("msub", &(ret + &index), idc, &[])
}

/// The elementary functions (:1106): applied to an atom they need no
/// parentheses, and a power of them goes on the name (`\sin ^{2}x`).
pub(super) fn is_elemfunc(op: &str) -> bool {
    const ELEM: [&str; 19] = [
        "ln", "sin", "cos", "sec", "csc", "tan", "cot", "asec", "acsc", "asin", "acos", "atan",
        "acot", "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
    ];
    ELEM.contains(&op) || ELEM.iter().any(|f| f.to_ascii_uppercase() == op)
}

/// The LaTeX of the functions with one (:1128-1217).
fn func_tex(name: &str) -> Option<&'static str> {
    Some(match name {
        "ln" => "\\ln ",
        "sin" => "\\sin ",
        "cos" => "\\cos ",
        "tan" => "\\tan ",
        "cot" => "\\cot ",
        "sinh" => "\\sinh ",
        "cosh" => "\\cosh ",
        "tanh" => "\\tanh ",
        "asin" => "\\arcsin ",
        "acos" => "\\arccos ",
        "atan" => "\\arctan ",
        "acot" => "\\operatorname{arccot}",
        "sec" => "\\sec ",
        "csc" => "\\csc ",
        "asec" => "\\operatorname{arcsec}",
        "acsc" => "\\operatorname{arccsc}",
        "asinh" => "\\operatorname{arsinh}",
        "acosh" => "\\operatorname{arcosh}",
        "atanh" => "\\operatorname{artanh}",
        _ => return None,
    })
}

/// The special functions with a symbol of their own (:1128-1217).
fn func_special(name: &str) -> Option<&'static str> {
    Some(match name {
        "Dirac" => "\\operatorname{\\delta}",
        "Heaviside" => "\\operatorname{\\theta}",
        "Gamma" => "\\Gamma",
        "Beta" => "\\Beta",
        "euler" => "\\phi",
        "Airy_Ai" => "\\operatorname{Ai}",
        "Airy_Bi" => "\\operatorname{Bi}",
        "Psi" => "\\Psi",
        "Zeta" => "\\zeta",
        _ => return None,
    })
}

/// A function's name in LaTeX (:1128): `\sin `, `\operatorname{arccot}`,
/// `\Gamma`, else the identifier's (`\mathrm{erf}`), in `\operatorname`
/// when it is several characters starting with a letter (:1219).
pub(super) fn func2markup(name: &str) -> String {
    let lower = name.to_ascii_lowercase();
    let low = if name.to_ascii_uppercase() == name && func_tex(&lower).is_some() {
        lower.as_str()
    } else {
        name
    };
    if let Some(tex) = func_tex(low) {
        return tex.to_owned();
    }
    let ret = match func_special(name) {
        Some(tex) => tex.to_owned(),
        None => idnt2markup(name).0,
    };
    if char_len(&ret) > 1 && ret.chars().next().is_some_and(my_isalpha) {
        format!("\\operatorname{{{ret}}}")
    } else {
        ret
    }
}

/// The `MathML` names of the functions with one (:1133-1175): content
/// element, presentation name.
fn func_mml(name: &str) -> Option<(&'static str, &'static str)> {
    Some(match name {
        "ln" => ("ln", "ln"),
        "sin" => ("sin", "sin"),
        "cos" => ("cos", "cos"),
        "tan" => ("tan", "tan"),
        "cot" => ("cot", "cot"),
        "sinh" => ("sinh", "sinh"),
        "cosh" => ("cosh", "cosh"),
        "tanh" => ("tanh", "tanh"),
        "asin" => ("arcsin", "arcsin"),
        "acos" => ("arccos", "arccos"),
        "atan" => ("arctan", "arctan"),
        "acot" => ("arccot", "arccot"),
        "sec" => ("sec", "sec"),
        "csc" => ("csc", "csc"),
        "asec" => ("arcsec", "arcsec"),
        "acsc" => ("arccsc", "arccsc"),
        "asinh" => ("arcsinh", "arsinh"),
        "acosh" => ("arccosh", "arcosh"),
        "atanh" => ("arctanh", "artanh"),
        _ => return None,
    })
}

/// A function's name in `MathML` (:1128, `typ` 3 presentation, 2 content),
/// its tag with id `idc`: `<mi>sin</mi>` / `<sin/>`, `&Gamma;` /
/// `<csymbol cd='hypergeo0'>gamma</csymbol>`, else the identifier's
/// (`<mi>erf</mi>`) / `<ci type='function'>erf</ci>`.
pub(super) fn func2mml(name: &str, idc: u32, content: bool) -> String {
    let lower = name.to_ascii_lowercase();
    let low = if name.to_ascii_uppercase() == name && func_mml(&lower).is_some() {
        lower.as_str()
    } else {
        name
    };
    let ret = if let Some((c, p)) = func_mml(low) {
        if content {
            format!("<{c}/>")
        } else {
            format!("<mi>{p}</mi>")
        }
    } else if content {
        let own = match name {
            "id" => "<ident/>",
            "gcd" => "<gcd/>",
            "lcm" => "<lcm/>",
            "irem" | "rem" => "<rem/>",
            "iquo" | "quo" => "<quotient/>",
            "arg" => "<arg/>",
            "max" => "<max/>",
            "min" => "<min/>",
            "det" => "<determinant/>",
            "Gamma" => return mml::insert_id(&mml::csymbol("gamma", "hypergeo0"), idc, true),
            "Beta" => return mml::insert_id(&mml::csymbol("beta", "hypergeo0"), idc, true),
            "euler" => return mml::insert_id(&mml::csymbol("euler", "integer2"), idc, true),
            "Airy_Ai" => return mml::insert_id(&mml::csymbol("Ai", "airy"), idc, true),
            "Airy_Bi" => return mml::insert_id(&mml::csymbol("Bi", "airy"), idc, true),
            _ => return mml::tag("ci", name, idc, &[("type", "function")]),
        };
        own.to_owned()
    } else {
        match name {
            "Dirac" => "<mi>&delta;</mi>",
            "Heaviside" => "<mi>&theta;</mi>",
            "Gamma" => "<mi mathvariant='normal'>&Gamma;</mi>",
            "Beta" => "<mi mathvariant='normal'>&Beta;</mi>",
            "euler" => "<mi mathvariant='normal'>&phi;</mi>",
            "Airy_Ai" => "<mi>Ai</mi>",
            "Airy_Bi" => "<mi>Bi</mi>",
            "Psi" => "<mi mathvariant='normal'>&Psi;</mi>",
            "Zeta" => "<mi mathvariant='normal'>&zeta;</mi>",
            _ => return idnt2mml(name, idc).0,
        }
        .to_owned()
    };
    mml::insert_id(&ret, idc, content)
}

/// giac's `str_to_tex` (:333) for a quoted string: `\ % _ $ { } ^`
/// escaped, in TeX quotes.
fn str_to_tex(s: &str) -> String {
    let mut ret = String::from("``");
    for c in mml::trim_string(s).chars() {
        if "\\%_${}^".contains(c) {
            ret.push('\\');
        }
        ret.push(c);
    }
    ret.push_str("''");
    ret
}

/// A string (:408): `\text{``s''}`; several lines an array of them, blank
/// lines but the last dropped.
pub(super) fn string2markup(s: &str) -> String {
    let lines = string_lines(s);
    let text: Vec<String> = lines
        .iter()
        .map(|l| format!("\\text{{{}}}", str_to_tex(l)))
        .collect();
    let text = text.join("\\\\");
    if lines.len() > 1 {
        format!("\\begin{{array}}{{l}}{text}\\end{{array}}")
    } else {
        text
    }
}

/// :416-426: the lines of a string, blank ones (spaces only) dropped but
/// the last.
fn string_lines(s: &str) -> Vec<&str> {
    let mut lines: Vec<&str> = s.split('\n').collect();
    let last = lines.pop().expect("split gives one piece at least");
    lines.retain(|l| !l.trim_matches(' ').is_empty());
    lines.push(last);
    lines
}

/// A string in `MathML` presentation (:453-464): `<ms>s</ms>` with id
/// `idc`; several lines a left-aligned `<mtable>` of them.
pub(super) fn string2mml(s: &str, idc: u32) -> String {
    match &string_lines(s)[..] {
        [line] => mml::tag("ms", &mml::str_to_mml(line), idc, &[]),
        lines => {
            let mut rows = String::new();
            for l in lines {
                rows.push_str("<mtr><mtd>");
                rows.push_str(&mml::str_to_mml_tag(l, "ms", 0));
                rows.push_str("</mtd></mtr>");
            }
            mml::tag("mtable", &rows, idc, &[("columnalign", "left")])
        }
    }
}
