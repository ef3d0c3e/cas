//! Set operations (CAS-DOC-32 §4, CAS-367), one check per input: random
//! finite sets and lists of small integers and a few symbols, and the
//! identities of their operations:
//! - `Union` and `Intersection` are commutative and idempotent;
//! - `|A ∪ B| = |A| + |B| − |A ∩ B|`;
//! - `A ∖ B` is included in `A` and meets `B` nowhere;
//! - De Morgan's laws relative to `U = A ∪ B`: `U ∖ (A ∩ B)` is `(U ∖ A) ∪
//!   (U ∖ B)`, and `U ∖ (A ∪ B)` is empty;
//! - `IncludedQ[A, B]`, where decided, is `Union[A, B] == B` as sets;
//! - `Element[x, A]`, where decided, is `IncludedQ[{x}, A]`.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol};
use libfuzzer_sys::fuzz_target;

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

fn call(head: Symbol, args: Vec<Expr>) -> Expr {
    Expr::call(head, args)
}

/// An element: a small integer, sometimes a symbol.
fn element(u: &mut Unstructured<'_>) -> Result<Expr> {
    Ok(if u.ratio(1u8, 6)? {
        Expr::symbol(["x", "y", "z"][u.int_in_range(0..=2usize)?])
    } else {
        Expr::from(u.int_in_range(-5i64..=5)?)
    })
}

/// A finite set or, sometimes, a list read as one.
fn set(u: &mut Unstructured<'_>) -> Result<Expr> {
    let n = u.int_in_range(0..=6usize)?;
    let mut items = Vec::with_capacity(n);
    for _ in 0..n {
        items.push(element(u)?);
    }
    let head = if u.ratio(1u8, 4)? {
        Symbol::LIST
    } else {
        Symbol::FINITE_SET
    };
    Ok(call(head, items))
}

fn size(e: &Expr) -> usize {
    e.args().len()
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let ctx = Context::new();
    let (a, b) = (set(u)?, set(u)?);
    let as_set = |e: &Expr| eval(&ctx, &call(Symbol::FINITE_SET, e.args().to_vec()));
    let (sa, sb) = (as_set(&a), as_set(&b));
    let union = eval(&ctx, &call(Symbol::UNION, vec![a.clone(), b.clone()]));
    let meet = eval(
        &ctx,
        &call(Symbol::INTERSECTION, vec![a.clone(), b.clone()]),
    );
    let swapped = eval(&ctx, &call(Symbol::UNION, vec![b.clone(), a.clone()]));
    assert_eq!(union, swapped, "Union[{a:?}, {b:?}]");
    let swapped = eval(
        &ctx,
        &call(Symbol::INTERSECTION, vec![b.clone(), a.clone()]),
    );
    assert_eq!(meet, swapped, "Intersection[{a:?}, {b:?}]");
    for (head, name) in [
        (Symbol::UNION, "Union"),
        (Symbol::INTERSECTION, "Intersection"),
    ] {
        let twice = eval(&ctx, &call(head, vec![a.clone(), a.clone()]));
        assert_eq!(twice, sa, "{name}[{a:?}, {a:?}]");
    }
    assert_eq!(
        size(&union) + size(&meet),
        size(&sa) + size(&sb),
        "|{a:?} ∪ {b:?}| = {union:?}, ∩ {meet:?}"
    );
    let difference = eval(&ctx, &call(Symbol::COMPLEMENT, vec![a.clone(), b.clone()]));
    let included = eval(
        &ctx,
        &call(Symbol::INCLUDED_Q, vec![difference.clone(), a.clone()]),
    );
    assert_eq!(included, Expr::from(Symbol::TRUE), "{difference:?} ⊆ {a:?}");
    let apart = eval(
        &ctx,
        &call(Symbol::INTERSECTION, vec![difference.clone(), b.clone()]),
    );
    assert_eq!(size(&apart), 0, "({a:?} ∖ {b:?}) ∩ {b:?} = {apart:?}");
    // De Morgan's laws relative to the union.
    let complement = |x: &Expr| {
        eval(
            &ctx,
            &call(Symbol::COMPLEMENT, vec![union.clone(), x.clone()]),
        )
    };
    let left = complement(&meet);
    let right = eval(
        &ctx,
        &call(Symbol::UNION, vec![complement(&a), complement(&b)]),
    );
    assert_eq!(left, right, "De Morgan for {a:?}, {b:?}");
    assert_eq!(size(&complement(&union)), 0, "U ∖ U for {a:?}, {b:?}");
    // Inclusion and membership, where decided.
    let included = eval(&ctx, &call(Symbol::INCLUDED_Q, vec![a.clone(), b.clone()]));
    if included.as_symbol() == Some(Symbol::TRUE) {
        assert_eq!(union, sb, "{a:?} ⊆ {b:?} but the union is {union:?}");
    }
    if included.as_symbol() == Some(Symbol::FALSE) {
        assert_ne!(union, sb, "{a:?} ⊄ {b:?} but the union is {union:?}");
    }
    let x = element(u)?;
    let member = eval(&ctx, &call(Symbol::ELEMENT, vec![x.clone(), b.clone()]));
    let single = eval(
        &ctx,
        &call(
            Symbol::INCLUDED_Q,
            vec![call(Symbol::FINITE_SET, vec![x.clone()]), b.clone()],
        ),
    );
    if member.as_symbol().is_some() && single.as_symbol().is_some() {
        assert_eq!(member, single, "{x:?} in {b:?}");
    }
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let _ = run(&mut Unstructured::new(data));
});
