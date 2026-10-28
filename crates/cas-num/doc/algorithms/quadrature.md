# Gauss–Legendre rules

`gauss_legendre(n, prec)` (`ball/quadrature.rs`, CAS-DOC-26 §2) gives the
`n`-point Gauss–Legendre rule on `[−1, 1]`:

- the nodes `x_{n,k}`, the zeros of the Legendre polynomial `P_n`, in
  descending order;
- the weights `w_{n,k} = 2/((1 − x²) P_n′(x)²)`.

Each is a ball that contains the exact value, of radius about
`2^(−prec)`. `cas-core`'s `NIntegrate` adds to the rule's sum the error
bound of [[Tre08]](references.md#tre08) for an integrand analytic in a
Bernstein ellipse, so the nodes and weights must be enclosures, not
approximations. The method is that of Arb's
`arb_hypgeom_legendre_p_ui_root` [[Joh17]](references.md#joh17).

## The polynomials

`P_n(x)` and `P_n′(x)` come from the three-term recurrence
`(m + 1) P_{m+1} = (2m + 1) x P_m − m P_{m−1}`, which is stable on
`[−1, 1]`, and from `P_n′ = n (x P_n − P_{n−1})/(x² − 1)`. Both are in ball
arithmetic, so at a ball `X` they enclose the values over `X`. The interval
Newton test below uses that enclosure.

The recurrence is stable in floating point, but ball arithmetic does not
know that: it sees no correlation between `P_m` and `P_{m−1}`, and a radius
grows by `1 + √2` at each step, `2^(1.27n)` over the recurrence.

- At a node `x̃`, the rounding errors grow by that factor.
- Over its ball `X`, the radius of `X` grows by it as well.

The interval Newton test needs `P_n′(X)` away from 0, and that costs about
twice the growth. So the work runs at `prec + 32 + 2.6n` bits. With only
`prec + 32`, the zeros of `P_64` failed the test at every precision, since
the ball of `N` and the one allowed for `X` grew together (CAS-322, found
by CI).

## The zeros

Only the zeros in `(0, 1)` are computed. Their negatives are the others,
and 0 is a zero for an odd `n`.

1. **Approximation.** The `k`-th largest zero starts from Tricomi's
   `(1 − (n − 1)/(8n³)) cos(π(4k − 1)/(4n + 2))`, in `f64`. Newton's
   iteration then runs on midpoints, three steps at 64 bits and three at
   each doubling, up to the working precision. The step is the quotient of
   the midpoints of `P_n(x̃)` and `P_n′(x̃)`. At 64 bits the ball of
   `P_64′` is `2³²` wide and holds 0, so dividing by the ball itself gave
   the indeterminate ball and a step to 0 (CAS-322, found by CI).
2. **Validation.** This is the interval Newton test
   [[Moo66]](references.md#moo66). Every zero of `P_n` in `X = x̃ ± r` lies
   in `N = x̃ − P_n(x̃)/P_n′(X)`, so `N ⊆ X` proves a zero in `N`. When
   `0 ∉ P_n′(X)`, it is the only zero in `X`.
   - The radius `r` starts at twice the Newton step's ball,
     `P_n(x̃)/P_n′(x̃)` at the exact node, plus `2^(−w)`, and grows by `2²`
     up to eight times.
   - Past that, the precision doubles, up to six times. A zero still not
     validated is a bug, and panics.
3. **Count.** Each zero must be positive and apart from the previous one.
   Then the `⌊n/2⌋` zeros found are distinct, and they are all of `P_n`'s
   positive zeros.

## The weights

The weights are `2/((1 − x²) P_n′(x)²)` at the node balls, at the same
working precision. Nodes and weights are then rounded to `prec` bits of
accuracy.

## The cache

Rules are cached process-wide per `(n, prec)`, at most 64 of them; past
that the cache starts afresh. `clear_constant_cache` frees them with the
constants.

## Tests

`ball/quadrature/tests.rs` checks:

- the 1-, 2- and 3-point rules against `0`, `±1/√3` and `±√(3/5)`, with
  weights 2, 1, `8/9` and `5/9`;
- that the weighted sums of `x^j` contain `∫_{−1}^{1} x^j` for every
  `j ≤ 2n − 1`, up to `n = 100`;
- that the nodes descend, are symmetric, and make `P_n` vanish;
- the cache.
