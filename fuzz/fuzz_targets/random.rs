//! Random numbers (CAS-DOC-33 §7, CAS-374), one check per input: a seed
//! and a form with its parameters, then
//! - the draw lies in its range or its law's support: reals in `[0, 1)`
//!   and between their ends, integers between theirs, distinct where they
//!   must be, permutations of what they permute, samples at distinct
//!   places, variates in their supports, vectors, matrices and polynomials
//!   of their sizes;
//! - the same seed draws the same again;
//! - the generator's naturals below a bound of up to 4096 bits stay below
//!   it.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_core::{Context, Expr, Symbol, View};
use cas_num::random::Generator;
use cas_num::{Integer, Natural, Rational};
use libfuzzer_sys::fuzz_target;

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn call(head: &str, args: Vec<Expr>) -> Expr {
    Expr::call(Symbol::new(head), args)
}

fn eval(ctx: &Context, e: &Expr) -> Expr {
    ctx.eval(e).unwrap_or_else(|err| panic!("{e:?}: {err}"))
}

/// A number as a rational: itself, or a ball's midpoint.
fn value(e: &Expr) -> Option<Rational> {
    match e.view() {
        View::Ball(b) => b.mid().checked_to_rational(),
        _ => e.to_rational(),
    }
}

fn integer(e: &Expr) -> Option<Integer> {
    e.to_integer()
}

/// Whether `x` lies between `a` and `b`, both included.
fn between(x: &Rational, a: &Rational, b: &Rational) -> bool {
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    lo <= x && x <= hi
}

fn small(u: &mut Unstructured<'_>) -> Result<Rational> {
    Ok(q(u.int_in_range(-12..=12)?, u.int_in_range(1..=4)?))
}

/// A probability `i/d` in `[0, 1]`, or in `(0, 1]` for `positive`.
fn probability(u: &mut Unstructured<'_>, positive: bool) -> Result<Rational> {
    let d = u.int_in_range(1..=8)?;
    let low = i64::from(positive);
    Ok(q(u.int_in_range(low..=d)?, d))
}

fn letters(n: usize) -> Vec<Expr> {
    (0..n)
        .map(|i| Expr::symbol(["a", "b", "c", "d", "e", "f", "g", "h"][i % 8]))
        .collect()
}

/// A form and the check of its draw.
type Check = Box<dyn Fn(&Expr) -> bool>;

/// A variate's form, and its support.
fn variate(u: &mut Unstructured<'_>) -> Result<(Expr, Check)> {
    let natural = |u: &mut Unstructured<'_>| -> Result<i64> { u.int_in_range(1..=6) };
    Ok(match u.int_in_range(0u8..=7)? {
        0 => {
            let n = u.int_in_range(0..=20i64)?;
            let p = probability(u, false)?;
            let form = call("RandomBinomial", vec![Expr::from(n), Expr::from(p)]);
            let check: Check = Box::new(move |x| {
                integer(x).is_some_and(|k| Integer::ZERO <= k && k <= Integer::from(n))
            });
            (form, check)
        }
        1 => {
            let lambda = q(u.int_in_range(1..=40)?, 4);
            let form = call("RandomPoisson", vec![Expr::from(lambda)]);
            (
                form,
                Box::new(|x| integer(x).is_some_and(|k| !k.is_negative())),
            )
        }
        2 => {
            let p = probability(u, true)?;
            let form = call("RandomGeometric", vec![Expr::from(p)]);
            (
                form,
                Box::new(|x| integer(x).is_some_and(|k| k >= Integer::ONE)),
            )
        }
        3 => {
            let lambda = q(u.int_in_range(1..=40)?, 4);
            let form = call("RandomExponential", vec![Expr::from(lambda)]);
            (form, Box::new(at_least_zero))
        }
        4 => {
            let sigma = q(u.int_in_range(1..=40)?, 4);
            let form = call(
                "RandomNormal",
                vec![Expr::from(small(u)?), Expr::from(sigma)],
            );
            (form, Box::new(|x| value(x).is_some()))
        }
        5 => {
            let form = call("RandomChiSquare", vec![Expr::from(natural(u)?)]);
            (form, Box::new(at_least_zero))
        }
        6 => {
            let form = call("RandomStudentT", vec![Expr::from(natural(u)?)]);
            (form, Box::new(|x| value(x).is_some()))
        }
        _ => {
            let args = vec![Expr::from(natural(u)?), Expr::from(natural(u)?)];
            (call("RandomFisher", args), Box::new(at_least_zero))
        }
    })
}

/// A real from 0 on: a ball that reaches 0 or above.
fn at_least_zero(x: &Expr) -> bool {
    match x.view() {
        View::Ball(b) => !b.mid().is_negative() || b.contains_rational(&Rational::ZERO),
        _ => value(x).is_some_and(|v| !v.is_negative()),
    }
}

/// A form of `Random`, `RandomInteger`, the samples and the permutations,
/// and its check.
fn uniform(u: &mut Unstructured<'_>) -> Result<(Expr, Check)> {
    Ok(match u.int_in_range(0u8..=7)? {
        0 => {
            let check: Check = Box::new(|x| {
                matches!(x.view(), View::Ball(b) if b.is_exact())
                    && value(x).is_some_and(|v| !v.is_negative() && v < Rational::ONE)
            });
            (call("Random", vec![]), check)
        }
        1 => {
            let n = Integer::from(u.int_in_range(-(1i64 << 40)..=1i64 << 40)?)
                * Integer::from(u.int_in_range(1i64..=1 << 20)?);
            if n.is_zero() {
                return Ok((call("Random", vec![Expr::from(1)]), Box::new(|_| true)));
            }
            let toward = if n.is_negative() {
                &n + &Integer::ONE
            } else {
                &n - &Integer::ONE
            };
            let (a, b) = (Rational::from(Integer::ZERO), Rational::from(toward));
            let check: Check =
                Box::new(move |x| integer(x).is_some_and(|k| between(&Rational::from(k), &a, &b)));
            (call("Random", vec![Expr::from(n)]), check)
        }
        2 => {
            let (a, b) = (small(u)?, small(u)?);
            let form = call("Random", vec![Expr::from(a.clone()), Expr::from(b.clone())]);
            // The ends' balls round: a little room.
            let slack = q(1, 1 << 20);
            let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
            let (lo, hi) = (&lo - &slack, &hi + &slack);
            let check: Check = Box::new(move |x| value(x).is_some_and(|v| between(&v, &lo, &hi)));
            (form, check)
        }
        3 => {
            let (a, b) = (u.int_in_range(-20..=20i64)?, u.int_in_range(-20..=20i64)?);
            let k = u.int_in_range(0..=12usize)?;
            let width = usize::try_from(a.abs_diff(b)).expect("small") + 1;
            let form = call(
                "Random",
                vec![Expr::from(k as i64), Expr::from(a), Expr::from(b)],
            );
            let check: Check = Box::new(move |x| {
                if k > width {
                    return x.head() == Some(&Expr::from(Symbol::new("Random")));
                }
                let (qa, qb) = (q(a, 1), q(b, 1));
                distinct_within(x, k, |v| {
                    value(v).is_some_and(|v| v.is_integer() && between(&v, &qa, &qb))
                })
            });
            (form, check)
        }
        4 => {
            let a = Integer::from(u.arbitrary::<i64>()?);
            let b = Integer::from(u.arbitrary::<i64>()?)
                * Integer::from(u.int_in_range(1i64..=1 << 30)?);
            let (qa, qb) = (Rational::from(a.clone()), Rational::from(b.clone()));
            let form = call("RandomInteger", vec![Expr::from(a), Expr::from(b)]);
            let check: Check = Box::new(move |x| {
                integer(x).is_some_and(|k| between(&Rational::from(k), &qa, &qb))
            });
            (form, check)
        }
        5 => {
            let n = u.int_in_range(1..=8usize)?;
            let items = Expr::call(Symbol::LIST, letters(n));
            let k = u.int_in_range(0..=9usize)?;
            let form = call("RandomSample", vec![items.clone(), Expr::from(k as i64)]);
            let check: Check = Box::new(move |x| {
                if k > n {
                    return x.head() == Some(&Expr::from(Symbol::new("RandomSample")));
                }
                distinct_within(x, k, |v| items.args().contains(v))
            });
            (form, check)
        }
        6 => {
            let n = u.int_in_range(0..=8usize)?;
            let items = letters(n);
            let form = call(
                "RandomPermutation",
                vec![Expr::call(Symbol::LIST, items.clone())],
            );
            let check: Check = Box::new(move |x| {
                let mut found = x.args().to_vec();
                let mut expected = items.clone();
                found.sort();
                expected.sort();
                found == expected
            });
            (form, check)
        }
        _ => {
            let n = u.int_in_range(0..=20i64)?;
            let form = call("RandomPermutation", vec![Expr::from(n)]);
            let check: Check = Box::new(move |x| {
                let mut found: Vec<Option<Integer>> = x.args().iter().map(integer).collect();
                found.sort();
                let expected: Vec<Option<Integer>> =
                    (0..n).map(|i| Some(Integer::from(i))).collect();
                found == expected
            });
            (form, check)
        }
    })
}

/// Whether `x` is a list of `k` distinct items that each pass `fits`.
fn distinct_within(x: &Expr, k: usize, fits: impl Fn(&Expr) -> bool) -> bool {
    let items = x.args();
    x.head() == Some(&Expr::from(Symbol::LIST))
        && items.len() == k
        && items.iter().all(&fits)
        && items
            .iter()
            .enumerate()
            .all(|(i, v)| !items[..i].contains(v))
}

/// A vector, a matrix or a polynomial of entries in `[0, m)` or the
/// default `[−99, 99]`, and its check.
fn compound(u: &mut Unstructured<'_>) -> Result<(Expr, Check)> {
    let n = u.int_in_range(0..=6i64)?;
    let (spec, low, high) = if u.arbitrary()? {
        let m = u.int_in_range(1..=9i64)?;
        (vec![Expr::from(m)], 0, m - 1)
    } else {
        (Vec::new(), -99, 99)
    };
    let fits = move |v: &Expr| {
        integer(v).is_some_and(|k| Integer::from(low) <= k && k <= Integer::from(high))
    };
    Ok(match u.int_in_range(0u8..=2)? {
        0 => {
            let mut args = vec![Expr::from(n)];
            args.extend(spec);
            let check: Check = Box::new(move |x| {
                x.head() == Some(&Expr::from(Symbol::LIST))
                    && x.args().len() == usize::try_from(n).expect("small")
                    && x.args().iter().all(fits)
            });
            (call("RandomVector", args), check)
        }
        1 => {
            let p = u.int_in_range(1..=4i64)?;
            let mut args = vec![Expr::from(n), Expr::from(p)];
            args.extend(spec);
            let check: Check = Box::new(move |x| {
                x.args().len() == usize::try_from(n).expect("small")
                    && x.args().iter().all(|row| {
                        row.args().len() == usize::try_from(p).expect("small")
                            && row.args().iter().all(fits)
                    })
            });
            (call("RandomMatrix", args), check)
        }
        _ => {
            let mut args = vec![Expr::from(n), Expr::symbol("x")];
            args.extend(spec);
            // Monic of degree n: its value at x = 2^20 is near 2^(20 n).
            let check: Check =
                Box::new(move |x| x.head() != Some(&Expr::from(Symbol::new("RandomPolynomial"))));
            (call("RandomPolynomial", args), check)
        }
    })
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let seed = Integer::from(u.arbitrary::<i64>()?);
    // The generator alone: a natural below a bound of up to 4096 bits.
    let words: Vec<u64> = (0..u.int_in_range(1..=64usize)?)
        .map(|_| u.arbitrary::<u64>())
        .collect::<Result<_>>()?;
    let bound = Natural::from_limbs(&words);
    if !bound.is_zero() {
        let mut g = Generator::seeded(&seed);
        assert!(g.below(&bound) < bound, "a draw past {bound}");
    }
    let (form, check) = match u.int_in_range(0u8..=2)? {
        0 => uniform(u)?,
        1 => variate(u)?,
        _ => compound(u)?,
    };
    let draw = || {
        let ctx = Context::new();
        ctx.seed_random(&seed);
        eval(&ctx, &form)
    };
    let first = draw();
    // A form that stays drew nothing it shows.
    if first.head() != form.head() || form.args().is_empty() {
        assert!(check(&first), "{form:?} drew {first:?}");
    }
    assert_eq!(draw(), first, "{form:?} from the seed {seed} again");
    Ok(())
}

fuzz_target!(|data: &[u8]| {
    let mut u = Unstructured::new(data);
    let _ = run(&mut u);
});
