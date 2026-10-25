//! giac's answers in this CAS's grammar: the NumWorks formatting of
//! `nws_caseval` undone, and giac's own notations rewritten where the
//! language differs.

/// The index of the bracket that closes the one opening at `open`.
fn closing(s: &[u8], open: usize) -> Option<usize> {
    let mut depth = 0usize;
    for (i, &c) in s.iter().enumerate().skip(open) {
        match c {
            b'[' | b'(' | b'{' => depth += 1,
            b']' | b')' | b'}' => {
                depth -= 1;
                if depth == 0 {
                    return Some(i);
                }
            }
            _ => {}
        }
    }
    None
}

/// Whether `s` is a decimal as giac prints approximations: `-0.5`,
/// `1.41421356237`, `2.80666336043e+43429`.
fn is_decimal(s: &str) -> bool {
    let s = s.strip_prefix('-').unwrap_or(s);
    let (mantissa, exponent) = match s.split_once('e') {
        Some((m, e)) => (m, Some(e.strip_prefix(['+', '-']).unwrap_or(e))),
        None => (s, None),
    };
    let digits = |t: &str| !t.is_empty() && t.bytes().all(|c| c.is_ascii_digit());
    let (int, frac) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    digits(int) && (frac.is_empty() || digits(frac)) && exponent.is_none_or(digits)
}

/// `nws_caseval`'s answer as `caseval` would print it:
/// - `√(…)` is `sqrt(…)`;
/// - an exact scalar is followed by `=` and its approximation, which goes
///   (lists and matrices have none);
/// - a matrix prints its rows without commas, `[[1,2][3,4]]`;
/// - a vector prints as a matrix of one row, `[[1,2]]`, so a one-row matrix
///   reads back as a vector: the formats cannot tell them apart.
///
/// An equation whose right side is a decimal loses it too.
pub(crate) fn undo_nws(answer: &str) -> String {
    let mut s = answer.replace('√', "sqrt");
    // The approximation after the last `=` outside brackets.
    let bytes = s.as_bytes();
    let mut depth = 0i32;
    let mut last = None;
    for (i, &c) in bytes.iter().enumerate() {
        match c {
            b'[' | b'(' | b'{' => depth += 1,
            b']' | b')' | b'}' => depth -= 1,
            b'=' if depth == 0 => last = Some(i),
            _ => {}
        }
    }
    if let Some(i) = last {
        if i > 0 && is_decimal(&s[i + 1..]) {
            s.truncate(i);
        }
    }
    let s = s.replace("][", "],[");
    // `[[…]]` that is one row: a vector.
    let bytes = s.as_bytes();
    if s.starts_with("[[") && s.ends_with("]]") && closing(bytes, 1) == Some(s.len() - 2) {
        return s[1..s.len() - 1].to_owned();
    }
    s
}

/// giac's answer in the language's grammar (CAS-DOC-7): its lists
/// `list[…]` and coefficient lists `poly1[…]` are `[…]`, its sets `set[…]`
/// are `{…}`, and its signed infinities `+infinity` and `-infinity`, where
/// they stand alone, are `inf` and `-inf` (`infinity` alone is unsigned,
/// the language's).
pub(crate) fn for_language(answer: &str) -> String {
    let mut s = signed_infinities(answer)
        .replace("list[", "[")
        .replace("poly1[", "[");
    while let Some(at) = s.find("set[") {
        let open = at + "set".len();
        let Some(close) = closing(s.as_bytes(), open) else {
            break;
        };
        s = format!("{}{{{}}}{}", &s[..at], &s[open + 1..close], &s[close + 1..]);
    }
    s
}

/// `+infinity` and `-infinity` as `inf` and `-inf` where the sign starts
/// the answer or follows `(`, `[` or `,`.
fn signed_infinities(answer: &str) -> String {
    let mut out = String::with_capacity(answer.len());
    let mut rest = answer;
    let mut previous: Option<char> = None;
    while let Some(c) = rest.chars().next() {
        let starts = previous.is_none_or(|p| matches!(p, '(' | '[' | ','));
        if starts && rest.starts_with("+infinity") {
            out.push_str("inf");
            rest = &rest["+infinity".len()..];
            previous = Some('f');
            continue;
        }
        if starts && rest.starts_with("-infinity") {
            out.push_str("-inf");
            rest = &rest["-infinity".len()..];
            previous = Some('f');
            continue;
        }
        out.push(c);
        previous = Some(c);
        rest = &rest[c.len_utf8()..];
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nws_formatting_is_undone() {
        for (nws, plain) in [
            ("(x-√(2))*(x+√(2))", "(x-sqrt(2))*(x+sqrt(2))"),
            ("√(3)-√(2)=0.317837245196", "sqrt(3)-sqrt(2)"),
            ("-1/2=-0.5", "-1/2"),
            ("1/3=0.333333333333", "1/3"),
            ("exp(1)=2.71828182846", "exp(1)"),
            ("2*x+1", "2*x+1"),
            ("3.14159265359", "3.14159265359"),
            ("[[x,1/3]]", "[x,1/3]"),
            ("[[1,2][3,4]]", "[[1,2],[3,4]]"),
            ("[[1/2,√(2)]]", "[1/2,sqrt(2)]"),
            ("x-1", "x-1"),
            (
                "1267650600228229401496703205376=1.26765060023e+30",
                "1267650600228229401496703205376",
            ),
            (
                "1/1267650600228229401496703205376=7.88860905221e-31",
                "1/1267650600228229401496703205376",
            ),
            ("-√(2)=-1.41421356237", "-sqrt(2)"),
        ] {
            assert_eq!(undo_nws(nws), plain, "{nws}");
        }
    }

    #[test]
    fn giac_notations_are_rewritten() {
        assert_eq!(for_language("list[1,2]"), "[1,2]");
        assert_eq!(for_language("poly1[1,3,0]"), "[1,3,0]");
        assert_eq!(for_language("set[1,2]"), "{1,2}");
        assert_eq!(for_language("[set[1,[2,3]],x]"), "[{1,[2,3]},x]");
        assert_eq!(for_language("x+1"), "x+1");
        assert_eq!(for_language("+infinity"), "inf");
        assert_eq!(for_language("-infinity"), "-inf");
        assert_eq!(for_language("infinity"), "infinity");
        assert_eq!(for_language("list[+infinity,-infinity]"), "[inf,-inf]");
        assert_eq!(for_language("x+infinity"), "x+infinity");
    }
}
