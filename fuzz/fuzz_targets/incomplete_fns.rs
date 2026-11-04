//! The incomplete gamma and beta functions of real balls (CAS-DOC-46 §5,
//! CAS-472), on balls of `a, b` in `(0, 64]` and `x` in `[0, 128]` (`[0,
//! 1]` for the beta), at up to 256 bits:
//! - `Q(a, x) Γ(a)` meets MPFR's `gamma_inc` bracket at the midpoint, the
//!   ends and an interior point of the box;
//! - `P + Q` holds 1, and `P(a + 1, x) = P(a, x) − xᵃ e^(−x)/Γ(a + 1)`
//!   overlaps;
//! - `I_x(a, b) + I_(1−x)(b, a)` holds 1, and `I_x(a + 1, b) = I_x(a, b) −
//!   xᵃ (1 − x)ᵇ/(a B(a, b))` overlaps;
//! - exact inputs give finite values with `prec − 3` bits.
#![no_main]

use arbitrary::{Result, Unstructured};
use cas_num::Integer;
use cas_num::ball::{Ball, Float, Mag};
use cas_num_fuzz::ball::{Oracle, bracket, meets, points, to_mp, warm_up};
use libfuzzer_sys::fuzz_target;
use rug::Float as MpFloat;

fuzz_target!(|data: &[u8]| {
    warm_up();
    let _ = run(&mut Unstructured::new(data));
});

/// A dyadic `m/2^k` in `(0, 2^top]` (or `[0, …]` when `zero`), and a ball
/// about it: exact, an ulp's radius, or a sizeable one.
fn ball(u: &mut Unstructured<'_>, top: i64, zero: bool, prec: u64) -> Result<Ball> {
    let k = u.int_in_range(0i64..=40)?;
    let m = u.int_in_range(i64::from(!zero)..=1i64 << 40)?;
    let x = Float::from_integer_exp(&Integer::from(m), &Integer::from(top - 40 - k.min(40)));
    let x = if top > 0 && u.arbitrary()? {
        x
    } else {
        x.mul_2exp(&Integer::from(-u.int_in_range(0i64..=8)?))
    };
    let r = match u.int_in_range(0u8..=3)? {
        0 | 1 => Mag::ZERO,
        2 => Mag::power_of_two(&(x.exponent() - Integer::from(prec))),
        _ => Mag::power_of_two(&(x.exponent() - Integer::from(u.int_in_range(4i64..=40)?))),
    };
    Ok(Ball::new(x, r))
}

fn one() -> Float {
    Float::from(Integer::ONE)
}

/// Whether two balls overlap.
fn overlap(x: &Ball, y: &Ball) -> bool {
    let gap = x.sub(y, 512);
    gap.contains(&Float::ZERO)
}

fn accurate(v: &Ball, prec: u64, what: &str) {
    assert!(v.is_finite(), "{what}: indeterminate");
    if !v.mid().is_zero() {
        let want = i64::try_from(prec).expect("small") - 3;
        assert!(v.rel_accuracy_bits() >= want, "{what}: {v:?}");
    }
}

fn run(u: &mut Unstructured<'_>) -> Result<()> {
    let prec = u.int_in_range(8u64..=256)?;
    if u.arbitrary()? {
        gamma(u, prec)
    } else {
        beta(u, prec)
    }
}

fn gamma(u: &mut Unstructured<'_>, prec: u64) -> Result<()> {
    let a = ball(u, 6, false, prec)?;
    let x = ball(u, 7, true, prec)?;
    let what = format!("({a:?}, {x:?}) at {prec}");
    let p = Ball::gamma_p(&a, &x, prec);
    let q = Ball::gamma_q(&a, &x, prec);
    if a.is_exact() && x.is_exact() {
        accurate(&p, prec, &format!("P{what}"));
        accurate(&q, prec, &format!("Q{what}"));
    }
    if p.is_finite() && q.is_finite() {
        assert!(p.add(&q, 512).contains(&one()), "P + Q{what}: {p:?} {q:?}");
    }
    // Γ(a, x) against MPFR at points of the box.
    let upper = q.mul(&a.gamma(prec + 16), prec + 16);
    if upper.is_finite() {
        for pa in points(u, &a, prec)? {
            for px in points(u, &x, prec)? {
                if pa.is_negative() || pa.is_zero() || px.is_negative() {
                    continue;
                }
                let (Some(am), Some(xm)) = (to_mp(&pa), to_mp(&px)) else {
                    continue;
                };
                let answer = bracket(prec + 64, |w, r| {
                    MpFloat::with_val_round(w, am.gamma_inc_ref(&xm), r).0
                });
                if let Oracle::Bracket(lo, hi) = answer {
                    assert!(
                        meets(&upper, &lo, &hi),
                        "Γ{what} at ({pa:?}, {px:?}): {upper:?} misses [{lo:?}, {hi:?}]"
                    );
                }
            }
        }
    }
    // The recurrence in a.
    let a1 = a.add(&Ball::from_integer(&Integer::ONE), 512);
    let next = Ball::gamma_p(&a1, &x, prec);
    if next.is_finite() && p.is_finite() && !x.mid().is_zero() {
        let w = prec + 64;
        let step = a.mul(&x.log(w), w).sub(&x, w).exp(w).div(&a1.gamma(w), w);
        let rhs = p.sub(&step, w);
        assert!(
            overlap(&next, &rhs),
            "P recurrence{what}: {next:?} vs {rhs:?}"
        );
    }
    Ok(())
}

fn beta(u: &mut Unstructured<'_>, prec: u64) -> Result<()> {
    let a = ball(u, 6, false, prec)?;
    let b = ball(u, 6, false, prec)?;
    let x = ball(u, 0, true, prec)?;
    let what = format!("({a:?}, {b:?}, {x:?}) at {prec}");
    let i = Ball::beta_i(&a, &b, &x, prec);
    if a.is_exact() && b.is_exact() && x.is_exact() {
        accurate(&i, prec, &format!("I{what}"));
    }
    let y = Ball::from_integer(&Integer::ONE).sub(&x, 512);
    let j = Ball::beta_i(&b, &a, &y, prec);
    if i.is_finite() && j.is_finite() {
        assert!(
            i.add(&j, 512).contains(&one()),
            "symmetry{what}: {i:?} {j:?}"
        );
    }
    let a1 = a.add(&Ball::from_integer(&Integer::ONE), 512);
    let next = Ball::beta_i(&a1, &b, &x, prec);
    let inside = !x.mid().is_zero() && x.mid() < &one();
    if next.is_finite() && i.is_finite() && inside {
        let w = prec + 64;
        let whole = a
            .gamma(w)
            .mul(&b.gamma(w), w)
            .div(&a.add(&b, w).gamma(w), w);
        let step = a
            .mul(&x.log(w), w)
            .add(&b.mul(&y.log(w), w), w)
            .exp(w)
            .div(&a.mul(&whole, w), w);
        let rhs = i.sub(&step, w);
        assert!(
            overlap(&next, &rhs),
            "I recurrence{what}: {next:?} vs {rhs:?}"
        );
    }
    Ok(())
}
