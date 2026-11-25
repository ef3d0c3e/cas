//! MathML tags as giac writes them (markup.cc :35-47, :206-302, :382-407):
//! an element with its id, ids put into and taken out of a tag, and text
//! escaped for MathML.

/// `<mo>&it;</mo>`, invisible times (:35).
pub(super) const ITIMES: &str = "<mo>&it;</mo>";
/// `<mo>&middot;</mo>` (:36).
pub(super) const CDOT: &str = "<mo>&middot;</mo>";
/// `<mo>+</mo>` (:39).
pub(super) const PLUS: &str = "<mo>+</mo>";
/// `<mo>&minus;</mo>` (:40).
pub(super) const MINUS: &str = "<mo>&minus;</mo>";
/// `<mo>&af;</mo>`, function application (:41).
pub(super) const APPLY: &str = "<mo>&af;</mo>";
/// `<mo>&ic;</mo>`, invisible comma (:42).
pub(super) const ICOMMA: &str = "<mo>&ic;</mo>";
/// The imaginary unit (:43).
pub(super) const I: &str = "<mi mathvariant='normal'>i</mi>";
/// Euler's number (:44).
pub(super) const E: &str = "<mi mathvariant='normal'>e</mi>";
/// The differential's d (:45).
pub(super) const D: &str = "<mi mathvariant='normal'>d</mi>";

/// The content elements (:206-239) that take an `id`, besides those
/// starting with `c`.
const CONTENT_ELEMENTS: [&str; 31] = [
    "cn",
    "ci",
    "cs",
    "csymbol",
    "cerror",
    "share",
    "bind",
    "bvar",
    "semantics",
    "apply",
    "domainofapplication",
    "interval",
    "condition",
    "lowlimit",
    "uplimit",
    "degree",
    "momentabout",
    "logbase",
    "inverse",
    "lambda",
    "piecewise",
    "piece",
    "otherwise",
    "set",
    "list",
    "vector",
    "matrix",
    "matrixrow",
    "declare",
    "reln",
    "fn",
];

/// `mml_tag` (:248): `<tag id='idN' k='v'…>s</tag>`, the id `id` for a
/// content element, `xref` for a presentation one, none when `idc` is 0;
/// attributes with an empty key or value are left out.
pub(super) fn tag(tag: &str, s: &str, idc: u32, attrs: &[(&str, &str)]) -> String {
    let mut r = format!("<{tag}");
    if idc != 0 {
        let key = if tag.starts_with('c') || CONTENT_ELEMENTS.contains(&tag) {
            "id"
        } else {
            "xref"
        };
        r.push_str(&format!(" {key}='id{idc}'"));
    }
    for (k, v) in attrs {
        if !k.is_empty() && !v.is_empty() {
            r.push_str(&format!(" {k}='{v}'"));
        }
    }
    format!("{r}>{s}</{tag}>")
}

/// `mml_csymbol` (:272): `<csymbol cd='cd'>s</csymbol>`.
pub(super) fn csymbol(s: &str, cd: &str) -> String {
    tag("csymbol", s, 0, &[("cd", cd)])
}

/// `extract_id` (:277): the id (`id` in content, `xref` in presentation)
/// of the first tag of `s`, taken out of it; "" when it has none.
pub(super) fn extract_id(s: &mut String, content: bool) -> String {
    let attr = if content { " id='" } else { " xref='" };
    let (Some(start), Some(gt)) = (s.find(attr), s.find('>')) else {
        return String::new();
    };
    if start > gt {
        return String::new();
    }
    let from = start + attr.len();
    let end = from + s[from..].find('\'').expect("a closed attribute");
    let id = s[from..end].to_owned();
    s.replace_range(start..=end, "");
    id
}

/// `insert_id` (:293): the id put into the first tag of `s`, before its
/// `>` or `/>`; `s` as it is when `idc` is 0.
pub(super) fn insert_id(s: &str, idc: u32, content: bool) -> String {
    if idc == 0 {
        return s.to_owned();
    }
    let Some(mut pos) = s.find('>').filter(|&p| p > 0) else {
        return s.to_owned();
    };
    if s.as_bytes()[pos - 1] == b'/' {
        pos -= 1;
    }
    let key = if content { "id" } else { "xref" };
    format!("{} {key}='id{idc}'{}", &s[..pos], &s[pos..])
}

/// giac's `trim_string` (:304): leading spaces, tabs and line breaks and
/// trailing ones dropped.
pub(super) fn trim_string(s: &str) -> &str {
    s.trim_start_matches([' ', '\n', '\t'])
        .trim_end_matches([' ', '\n', '\t'])
}

/// `str_to_mml` (:382) without a tag: `s` trimmed, `&`, `<`, `>` escaped.
pub(super) fn str_to_mml(s: &str) -> String {
    trim_string(s)
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// `str_to_mml` (:382) in a tag with id `idc`.
pub(super) fn str_to_mml_tag(s: &str, t: &str, idc: u32) -> String {
    tag(t, &str_to_mml(s), idc, &[])
}
