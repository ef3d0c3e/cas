//! Prime counting and the n-th prime by a segmented sieve of Eratosthenes
//! over the odd numbers (CAS-271).

#[cfg(test)]
mod tests;

/// The largest `x` [`prime_pi`] counts to: its sieve takes a fraction of a
/// second there.
pub const PRIME_PI_LIMIT: u64 = 100_000_000;

/// The largest `n` [`nth_prime`] finds: the millionth prime is 15 485 863.
pub const NTH_PRIME_LIMIT: u64 = 1_000_000;

/// Odd numbers per segment of the sieve.
const SEGMENT: u64 = 1 << 18;

/// `π(x)`, the number of primes up to `x`; `None` past [`PRIME_PI_LIMIT`].
#[must_use]
pub fn prime_pi(x: u64) -> Option<u64> {
    if x > PRIME_PI_LIMIT {
        return None;
    }
    let mut count = 0;
    sieve(x, |_| {
        count += 1;
        true
    });
    Some(count)
}

/// The `n`-th prime, the first being 2; `None` for 0 and past
/// [`NTH_PRIME_LIMIT`].
#[must_use]
#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "an upper bound for n below 10^6"
)]
pub fn nth_prime(n: u64) -> Option<u64> {
    if n == 0 || n > NTH_PRIME_LIMIT {
        return None;
    }
    // Rosser: p_n < n (ln n + ln ln n) for n ≥ 6.
    let bound = if n < 6 {
        13
    } else {
        let f = n as f64;
        (f * (f.ln() + f.ln().ln())).ceil() as u64
    };
    let (mut k, mut found) = (0, None);
    sieve(bound, |p| {
        k += 1;
        if k == n {
            found = Some(p);
            false
        } else {
            true
        }
    });
    found
}

/// Calls `visit` with each prime up to `limit`, increasing, until it
/// returns `false`.
fn sieve(limit: u64, mut visit: impl FnMut(u64) -> bool) {
    if limit < 2 || !visit(2) {
        return;
    }
    let last = if limit % 2 == 0 { limit - 1 } else { limit };
    if last < 3 {
        return;
    }
    // The odd primes up to √last, by a plain sieve.
    let root = last.isqrt();
    let mut small = vec![true; usize::try_from(root + 1).expect("small")];
    let mut base = Vec::new();
    for p in (3..=root).step_by(2) {
        if small[usize::try_from(p).expect("small")] {
            base.push(p);
            let mut m = p * p;
            while m <= root {
                small[usize::try_from(m).expect("small")] = false;
                m += 2 * p;
            }
        }
    }
    let mut composite = vec![false; usize::try_from(SEGMENT).expect("small")];
    let mut low = 3;
    while low <= last {
        let high = (low + 2 * (SEGMENT - 1)).min(last);
        let len = usize::try_from((high - low) / 2 + 1).expect("a segment");
        composite[..len].fill(false);
        for &p in &base {
            if p * p > high {
                break;
            }
            // The first odd multiple of p from max(p², low).
            let mut m = (p * p).max(low.div_ceil(p) * p);
            if m % 2 == 0 {
                m += p;
            }
            while m <= high {
                composite[usize::try_from((m - low) / 2).expect("in the segment")] = true;
                m += 2 * p;
            }
        }
        for (i, &c) in composite[..len].iter().enumerate() {
            if !c && !visit(low + 2 * u64::try_from(i).expect("small")) {
                return;
            }
        }
        low = high + 2;
    }
}
