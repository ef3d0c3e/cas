# Algorithms of the ball arithmetic

This is the algorithms documentation of `cas_num::ball` and the code under
it, after GMP's "Algorithms" chapter [[GMP]](references.md#gmp) and MPFR's
*Algorithms and Proofs* [[MPFR]](references.md#mpfr). For each operation it
gives:

- what it computes, and the precision and error it promises;
- the base algorithm, and where it comes from;
- the tricks on top of it, and why each is there;
- where the error bound is proved, usually the doc comment of the function
  named, and how it is tested.

The design that chose these algorithms is CAS-DOC-5, and the tickets that
implemented them are listed in CAS-DOC-6 (both in Lific). This
documentation describes the code as it is.

## Scope

It covers the code of the Balls module (CAS-71 … CAS-104):

- the approximate integer kernels for ball midpoints (`kernel/mul/mulhi.rs`
  and `kernel/div/approx.rs`);
- everything under `ball/`;
- the Gaussian rationals of `complex/` (CAS-50);
- the exact functions of [integers.md](integers.md): the factorial
  (CAS-123), and the simplest rational in an interval (CAS-148).

The exact integer arithmetic it builds on (the `kernel` multiplication and
division algorithms, `Natural`, `Integer`) is otherwise outside this
documentation.
When an algorithm here calls into it, the call is named (for example "the
exact integer square root, `Natural::sqrt_rem`") with that algorithm's
reference.

## How the error bounds fit together

The arithmetic is built in three layers. Each layer proves only what the
layer below it cannot:

1. **Floats** (`Float`) are exact, or correctly rounded in a given mode. The
   approximate entries for ball midpoints return a one-sided result within a
   stated number of units. Each such bound is proved next to its code.
2. **Balls** (`Ball`, `ComplexBall`) propagate radii through every operation
   ([balls.md](balls.md)). A computation written as a composition of ball
   operations needs no error analysis of its own: the radius is the
   analysis.
3. **Series** are the only other places with hand-proved bounds:
   - sums in fixed point (`ball/fns/fixed.rs`) carry a run-time error count;
   - binary splitting (`ball/bsplit.rs`) sums exactly, and each series bounds
     its own tail, which then goes into a radius;
   - the machine-ball kernels have a priori error budgets in double-double
     arithmetic.

Every bound is also tested:

- **Oracles.** MPFR, or MPC for complex balls, at `p + 64` bits (or
  `p + 192`) in both directed roundings, so the oracle is a bracket rather
  than a single value.
- **Tail checks.** Each tail bound is compared with the true tail.
- **Few-term runs.** Constants are evaluated with too few terms, where the
  tail dominates the radius and must be there.
- **Stress tunings.** Two build configurations move every algorithm
  threshold to an extreme.
- **Sabotage runs.** Deliberate faults in the code are each checked to fail
  some test.
- **Fuzzing.** libFuzzer targets in `fuzz/` (CAS-113) check the functions of
  real, complex and machine balls against MPFR or MPC at the midpoint, the
  ends and an interior point of fuzzer-chosen balls, and the constants from
  a cleared cache. Their inputs sit next to the points where the functions
  change behaviour (0, ±1, small integers and half-integers, multiples of
  π/2), with tiny to huge exponents and precisions across the table and
  bit-burst thresholds.

## Notation

- `B = 2^64`, the limb base. `p` or `prec` is a precision in bits, and `w`
  (or `wp`) a working precision with guard bits.
- An ulp at `p` bits of a float `x` is `2^(exp(x) - p)`, with
  `x = ±0.m × 2^exp(x)` and `1/2 <= 0.m < 1`.
- `[m ± r]` is the ball of midpoint `m` and radius `r`, the set of reals
  within `r` of `m`. `ρ` is an input's radius.
- A *lower bound with error `e`* (the approximate kernels) is a result `r`
  with `r <= t <= r + e` for the true value `t`, in units of the last limb.
- `u = 2^-53` is the unit roundoff of `f64`. `OP_ERR = 2^-102` is the
  relative error bound used for every double-double operation.
- Citations in brackets, such as [[BZ10]](references.md#bz10), refer to the
  [bibliography](references.md). The code's doc comments use the same keys.

## Contents

| Chapter | Covers |
|---|---|
| [kernels.md](kernels.md) | Short products; approximate reciprocal, quotient, inverse square root and square root |
| [float.md](float.md) | `Float`: representation, rounding, add/sub/mul, division and square root (exact, correctly rounded and approximate), fixed-size paths, exact signs of sums, huge exponents |
| [balls.md](balls.md) | `Mag` and `SmallMag`; ball arithmetic and predicates; decimal and hexadecimal I/O |
| [constants.md](constants.md) | Binary splitting; π, log 2, e, Euler's γ, Catalan's G; the constant cache |
| [functions.md](functions.md) | Real elementary functions: fixed-point series, `exp`, `log`, trigonometric and inverse, `pow`, `root`, hyperbolic, table and bit-burst cores |
| [complex.md](complex.md) | Complex balls: arithmetic, elementary functions, branch cuts; Gaussian rationals |
| [machine.md](machine.md) | `MachineBall`: `f64` balls, error-free transformations, double-double kernels and their budgets |
| [integers.md](integers.md) | Exact integer functions: the factorial |
| [references.md](references.md) | Bibliography |

## Keeping it current

The rule for every ticket that adds or changes an algorithm here:

- The change updates the chapter in the same commit: base algorithm,
  tricks, error bound and its location, and references.
- A new source gets a full bibliography entry, and the code cites it by key.
- A trick of our own is marked as ours, with the ticket that introduced it.
  That lets a reader tell a published method from a local one.
