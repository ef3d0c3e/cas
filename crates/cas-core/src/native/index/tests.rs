//! Tests of indexing: tables checked against giac 1.9 where Xcas agrees,
//! and a reference implementation on random lists, matrices and cuboids.

use proptest::prelude::*;

use crate::{Context, EvalError, Expr, Symbol};

fn read(text: &str) -> Expr {
    Expr::parse_full_form(text).unwrap_or_else(|e| panic!("{text:?}: {e}"))
}

/// `(input, result)` in full form, `{l}` and `{m}` standing for a list and
/// a matrix.
fn check(cases: &[(&str, &str)]) {
    let ctx = Context::new();
    let mut wrong = Vec::new();
    for (input, expected) in cases {
        let input = input.replace("{l}", L).replace("{m}", M);
        match ctx.eval(&read(&input)) {
            Ok(found) if format!("{found:?}") == *expected => {}
            found => wrong.push(format!("{input}: {found:?}, expected {expected}")),
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

fn check_errors(cases: &[(&str, EvalError)]) {
    let ctx = Context::new();
    for (input, error) in cases {
        let input = input.replace("{l}", L).replace("{m}", M);
        assert_eq!(ctx.eval(&read(&input)), Err(*error), "{input}");
    }
}

const L: &str = "List[10, 20, 30, 40]";
const M: &str = "List[List[1, 2, 3], List[4, 5, 6], List[7, 8, 9]]";

#[test]
fn elements_and_slices() {
    check(&[
        ("At[{l}, 0]", "10"),
        ("At[{l}, 3]", "40"),
        ("At[{l}, -1]", "40"),
        ("At[{l}, -4]", "10"),
        ("At[{l}, Range[1, 2]]", "List[20, 30]"),
        ("At[{l}, Range[1, -1]]", "List[20, 30, 40]"),
        ("At[{l}, Range[-2, -1]]", "List[30, 40]"),
        ("At[{l}, Range[0, 0]]", "List[10]"),
        ("At[{l}, Range[0, -1]]", "List[10, 20, 30, 40]"),
        // Empty slices, as in Xcas: the end before the start, also just
        // past the end or before the first element.
        ("At[{l}, Range[2, 1]]", "List[]"),
        ("At[{l}, Range[3, 0]]", "List[]"),
        ("At[{l}, Range[4, 3]]", "List[]"),
        ("At[List[10], Range[1, -1]]", "List[]"),
        ("At[List[], Range[0, -1]]", "List[]"),
        // The indices evaluate.
        ("At[{l}, Plus[1, 1]]", "30"),
        ("At[{l}, Range[Plus[0, 1], Times[2, 1]]]", "List[20, 30]"),
        // An index that is not an integer, or a list that is not one, stays:
        // `x[0]` is an indexed name.
        ("At[{l}, k]", "At[List[10, 20, 30, 40], k]"),
        ("At[{l}, 1/2]", "At[List[10, 20, 30, 40], 1/2]"),
        (
            "At[{l}, Range[0, n]]",
            "At[List[10, 20, 30, 40], Range[0, n]]",
        ),
        (
            "At[{l}, Range[0, 1, 2]]",
            "At[List[10, 20, 30, 40], Range[0, 1, 2]]",
        ),
        ("At[x, 0]", "At[x, 0]"),
        ("At[f[1, 2], 0]", "At[f[1, 2], 0]"),
        ("At[5, 0]", "At[5, 0]"),
        ("At[{l}]", "At[List[10, 20, 30, 40]]"),
    ]);
}

#[test]
fn matrices() {
    check(&[
        ("At[{m}, 1]", "List[4, 5, 6]"),
        ("At[{m}, 1, 2]", "6"),
        ("At[At[{m}, 1], 2]", "6"),
        ("At[{m}, -1, -1]", "9"),
        // A column slice, a row slice and a block.
        ("At[{m}, Range[0, 1], 1]", "List[2, 5]"),
        ("At[{m}, Range[0, -1], 0]", "List[1, 4, 7]"),
        ("At[{m}, 1, Range[0, 1]]", "List[4, 5]"),
        (
            "At[{m}, Range[0, 1], Range[0, 1]]",
            "List[List[1, 2], List[4, 5]]",
        ),
        ("At[{m}, Range[1, 0], 0]", "List[]"),
        // Past a level that is not a list, what is left stays.
        ("At[List[x, List[1, 2]], 0, 1]", "At[x, 1]"),
        ("At[{m}, 0, k]", "At[List[1, 2, 3], k]"),
        (
            "At[{m}, Range[0, 1], k]",
            "List[At[List[1, 2, 3], k], At[List[4, 5, 6], k]]",
        ),
        ("At[{m}, 0, 0, 0]", "At[1, 0]"),
    ]);
}

#[test]
fn out_of_range() {
    // Unlike Xcas, a slice past the end is an error, not cut short.
    check_errors(&[
        ("At[{l}, 4]", EvalError::IndexOutOfRange),
        ("At[{l}, -5]", EvalError::IndexOutOfRange),
        ("At[List[], 0]", EvalError::IndexOutOfRange),
        (
            "At[{l}, 1000000000000000000000000000000]",
            EvalError::IndexOutOfRange,
        ),
        (
            "At[{l}, -1000000000000000000000000000000]",
            EvalError::IndexOutOfRange,
        ),
        ("At[{l}, Range[3, 5]]", EvalError::IndexOutOfRange),
        ("At[{l}, Range[5, 4]]", EvalError::IndexOutOfRange),
        ("At[{l}, Range[-5, 0]]", EvalError::IndexOutOfRange),
        ("At[{m}, 3, 0]", EvalError::IndexOutOfRange),
        ("At[{m}, 0, 3]", EvalError::IndexOutOfRange),
        ("At[{m}, Range[0, 1], 3]", EvalError::IndexOutOfRange),
    ]);
}

#[test]
fn replacements() {
    check(&[
        ("ReplaceAt[{l}, 1, 99]", "List[10, 99, 30, 40]"),
        ("ReplaceAt[{l}, -1, 7]", "List[10, 20, 30, 7]"),
        (
            "ReplaceAt[{l}, Range[1, 2], List[a, b]]",
            "List[10, a, b, 40]",
        ),
        (
            "ReplaceAt[{l}, Range[2, 1], List[]]",
            "List[10, 20, 30, 40]",
        ),
        (
            "ReplaceAt[{m}, 1, 2, 0]",
            "List[List[1, 2, 3], List[4, 5, 0], List[7, 8, 9]]",
        ),
        (
            "ReplaceAt[{m}, 0, List[7, 7, 7]]",
            "List[List[7, 7, 7], List[4, 5, 6], List[7, 8, 9]]",
        ),
        (
            "ReplaceAt[{m}, Range[0, 1], 0, List[a, b]]",
            "List[List[a, 2, 3], List[b, 5, 6], List[7, 8, 9]]",
        ),
        (
            "ReplaceAt[{m}, Range[1, 2], Range[1, 2], List[List[a, b], List[c, d]]]",
            "List[List[1, 2, 3], List[4, a, b], List[7, c, d]]",
        ),
        // A list replaces one element whole.
        (
            "ReplaceAt[{l}, 0, List[1, 2]]",
            "List[List[1, 2], 20, 30, 40]",
        ),
        // What cannot be replaced stays.
        (
            "ReplaceAt[{l}, k, 1]",
            "ReplaceAt[List[10, 20, 30, 40], k, 1]",
        ),
        ("ReplaceAt[x, 0, 1]", "ReplaceAt[x, 0, 1]"),
        (
            "ReplaceAt[List[x, y], 0, 0, 1]",
            "ReplaceAt[List[x, y], 0, 0, 1]",
        ),
        ("ReplaceAt[{l}, 1]", "ReplaceAt[List[10, 20, 30, 40], 1]"),
    ]);
    // Unlike Xcas: no extension past the end, and no scalar spread over a
    // slice (CAS-DOC-4 O8).
    check_errors(&[
        ("ReplaceAt[{l}, 4, 1]", EvalError::IndexOutOfRange),
        ("ReplaceAt[{l}, 6, 1]", EvalError::IndexOutOfRange),
        ("ReplaceAt[{m}, 0, 3, 1]", EvalError::IndexOutOfRange),
        ("ReplaceAt[{l}, Range[1, 2], 5]", EvalError::Shape),
        (
            "ReplaceAt[{l}, Range[1, 2], List[a, b, c]]",
            EvalError::Shape,
        ),
        (
            "ReplaceAt[{m}, Range[0, 1], Range[0, 1], List[a, b]]",
            EvalError::Shape,
        ),
    ]);
}

/// Each level of indexing takes a level of `Limits::recursion`, so deep
/// indices stop with `RecursionLimit` before the stack runs out.
#[test]
fn deep_indices_count_against_the_recursion_limit() {
    let ctx = Context::new();
    let mut deep = Expr::from(1);
    for _ in 0..200 {
        deep = list(&[deep]);
    }
    let zeros = vec![Expr::from(0); 200];
    let at = Expr::apply(&Symbol::AT.into(), &[&[deep.clone()], &zeros[..]].concat());
    let replace = Expr::apply(
        &Symbol::REPLACE_AT.into(),
        &[&[deep.clone()], &zeros[..], &[Expr::from(2)]].concat(),
    );
    let mut ev = crate::eval::Eval::new(&ctx, 100);
    assert_eq!(super::at(&mut ev, &at), Err(EvalError::RecursionLimit));
    assert_eq!(
        super::replace_at(&mut ev, &replace),
        Err(EvalError::RecursionLimit)
    );
    let mut ev = crate::eval::Eval::new(&ctx, 300);
    assert_eq!(super::at(&mut ev, &at), Ok(Some(Expr::from(1))));
    let mut two = Expr::from(2);
    for _ in 0..200 {
        two = list(&[two]);
    }
    assert_eq!(super::replace_at(&mut ev, &replace), Ok(Some(two)));
}

/// An index as the strategies make it.
#[derive(Clone, Debug)]
enum Ix {
    Int(i64),
    Slice(i64, i64),
    /// A symbol, or a number that is not an integer.
    Other(Expr),
}

impl Ix {
    fn to_expr(&self) -> Expr {
        match self {
            Ix::Int(k) => Expr::from(*k),
            Ix::Slice(a, b) => {
                Expr::apply(&Symbol::RANGE.into(), &[Expr::from(*a), Expr::from(*b)])
            }
            Ix::Other(e) => e.clone(),
        }
    }
}

fn index() -> impl Strategy<Value = Ix> {
    prop_oneof![
        5 => (-6i64..6).prop_map(Ix::Int),
        4 => (-6i64..6, -6i64..6).prop_map(|(a, b)| Ix::Slice(a, b)),
        1 => prop::sample::select(vec![read("k"), read("1/2")]).prop_map(Ix::Other),
    ]
}

fn list(items: &[Expr]) -> Expr {
    Expr::apply(&Symbol::LIST.into(), items)
}

/// A vector, a matrix or a cuboid of small integers and symbols, of up to
/// four elements a side.
fn array() -> impl Strategy<Value = Expr> {
    let entry = prop_oneof![(-9i64..10).prop_map(Expr::from), Just(read("x"))];
    prop::collection::vec(0usize..=4, 1..=3).prop_flat_map(move |sides| {
        let count = sides.iter().product::<usize>();
        prop::collection::vec(entry.clone(), count)
            .prop_map(move |entries| nest(&sides, &mut entries.into_iter()))
    })
}

/// Lists of the sides `sides`, outermost first, of `entries` in order.
fn nest(sides: &[usize], entries: &mut impl Iterator<Item = Expr>) -> Expr {
    match sides.split_first() {
        None => entries.next().expect("one entry for each place"),
        Some((&side, rest)) => list(&(0..side).map(|_| nest(rest, entries)).collect::<Vec<_>>()),
    }
}

/// What the reference makes of an index into `len` elements: the
/// positions it selects, by filtering every position, and whether as a
/// slice; `None` for one that is not an index.
fn positions(ix: &Ix, len: usize) -> Option<Result<(Vec<usize>, bool), EvalError>> {
    let n = i64::try_from(len).unwrap();
    let from_end = |k: i64| if k < 0 { k + n } else { k };
    match *ix {
        Ix::Other(_) => None,
        Ix::Int(k) => {
            let k = from_end(k);
            let all: Vec<usize> = (0..len)
                .filter(|&p| i64::try_from(p).unwrap() == k)
                .collect();
            Some(if all.is_empty() {
                Err(EvalError::IndexOutOfRange)
            } else {
                Ok((all, false))
            })
        }
        Ix::Slice(a, b) => {
            let (a, b) = (from_end(a), from_end(b));
            if a < 0 || a > n || b < -1 || b >= n {
                return Some(Err(EvalError::IndexOutOfRange));
            }
            let all = (0..len)
                .filter(|&p| (a..=b).contains(&i64::try_from(p).unwrap()))
                .collect();
            Some(Ok((all, true)))
        }
    }
}

fn is_list(e: &Expr) -> bool {
    e.head().and_then(Expr::as_symbol) == Some(Symbol::LIST)
}

fn apply(head: Symbol, first: &Expr, ixs: &[Ix], last: Option<&Expr>) -> Expr {
    let mut args = vec![first.clone()];
    args.extend(ixs.iter().map(Ix::to_expr));
    args.extend(last.cloned());
    Expr::apply(&head.into(), &args)
}

/// The reference `At`: `None` when it stays as it is.
fn reference_at(v: &Expr, ixs: &[Ix]) -> Result<Option<Expr>, EvalError> {
    let Some((ix, rest)) = ixs.split_first() else {
        return Ok(Some(v.clone()));
    };
    if !is_list(v) {
        return Ok(None);
    }
    let Some(selected) = positions(ix, v.args().len()) else {
        return Ok(None);
    };
    let (positions, slice) = selected?;
    let mut parts = Vec::new();
    for p in positions {
        let x = &v.args()[p];
        parts.push(reference_at(x, rest)?.unwrap_or_else(|| apply(Symbol::AT, x, rest, None)));
    }
    Ok(Some(if slice { list(&parts) } else { parts.remove(0) }))
}

/// The reference `ReplaceAt`.
fn reference_replace(v: &Expr, ixs: &[Ix], x: &Expr) -> Result<Option<Expr>, EvalError> {
    let Some((ix, rest)) = ixs.split_first() else {
        return Ok(Some(x.clone()));
    };
    if !is_list(v) {
        return Ok(None);
    }
    let Some(selected) = positions(ix, v.args().len()) else {
        return Ok(None);
    };
    let (positions, slice) = selected?;
    let values = if slice {
        if !is_list(x) || x.args().len() != positions.len() {
            return Err(EvalError::Shape);
        }
        x.args().to_vec()
    } else {
        vec![x.clone()]
    };
    let mut elements = v.args().to_vec();
    for (p, value) in positions.into_iter().zip(values) {
        match reference_replace(&elements[p], rest, &value)? {
            Some(new) => elements[p] = new,
            None => return Ok(None),
        }
    }
    Ok(Some(list(&elements)))
}

/// A value with the shape of what `v[ixs]` selects, of new symbols, if the
/// selection is a value.
fn replacement(selected: &Expr, counter: &mut usize) -> Expr {
    if is_list(selected) {
        let parts: Vec<Expr> = selected
            .args()
            .iter()
            .map(|e| replacement(e, counter))
            .collect();
        list(&parts)
    } else {
        *counter += 1;
        Symbol::new(&format!("r{counter}")).into()
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2000))]

    /// `At` agrees with the reference for every index form, also past the
    /// levels of the array.
    #[test]
    fn indexing_agrees_with_the_reference(
        v in array(),
        ixs in prop::collection::vec(index(), 1..=4),
    ) {
        let e = apply(Symbol::AT, &v, &ixs, None);
        let expected = reference_at(&v, &ixs).map(|r| r.unwrap_or_else(|| e.clone()));
        prop_assert_eq!(Context::new().eval(&e), expected);
    }

    /// `ReplaceAt` agrees with the reference, with values of the shape of
    /// the selection and of other shapes; after a replacement, `At` gives
    /// the value back.
    #[test]
    fn replacing_agrees_with_the_reference(
        v in array(),
        ixs in prop::collection::vec(index(), 1..=4),
        scalar in any::<bool>(),
    ) {
        let ctx = Context::new();
        let at = apply(Symbol::AT, &v, &ixs, None);
        let x = match ctx.eval(&at) {
            Ok(selected) if !scalar => replacement(&selected, &mut 0),
            _ => read("z"),
        };
        let e = apply(Symbol::REPLACE_AT, &v, &ixs, Some(&x));
        let expected = reference_replace(&v, &ixs, &x).map(|r| r.unwrap_or_else(|| e.clone()));
        let found = ctx.eval(&e);
        prop_assert_eq!(&found, &expected);
        if let (Ok(replaced), Ok(Some(_))) = (&found, reference_replace(&v, &ixs, &x)) {
            prop_assert_eq!(ctx.eval(&apply(Symbol::AT, replaced, &ixs, None)), Ok(x));
        }
    }
}
