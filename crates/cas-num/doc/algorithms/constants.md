# Binary splitting and constants

## Binary splitting

`ball/bsplit.rs` (CAS-78, CAS-95). This is the method of Haible and
Papanikolaou [[HP98]](references.md#hp98) [[BZ10]](references.md#bz10)
(§4.9), for series

```text
S = Σ_{k < n} a(k) / b(k) · p(0) p(1) … p(k) / (q(0) q(1) … q(k))
```

with integer `p`, `q`, `a`, `b` (the `Series` trait). Over a range of terms,
with the products starting at its left end, it keeps `P = Π p`, `Q = Π q`,
`B = Π b` and `T = B Q S`, all integers:

- **Leaf.** `P = p(k)`, `Q = q(k)`, `B = b(k)`, `T = a(k) p(k)`.
- **Merge.** Two halves combine as

  ```text
  P = Pl Pr,  Q = Ql Qr,  B = Bl Br,  T = Br Qr Tl + Bl Pl Tr,
  ```

  since the right half's sum is `Pl / Ql` times its own.
- **Result.** The sum is `T / (B Q)`: one division, after a tree of products
  of growing size. The cost is `O(M(N) log N)` for numbers of `N` bits at the
  top, so its speed is the integer multiplication's.

**Tricks.**
- **Factors of 1.** A factor that is always 1 is `None` and costs nothing.
- **Right spine.** `P` is never formed on the rightmost path of the tree,
  where it isn't needed.
- **Powers of two as shifts (CAS-95, ours).** A power of two in `q(k)`, as in
  the bit-burst series of `x = p / 2^b`, is kept apart as
  `Series::shift(k)`. A range carries `E = Σ shift` instead of multiplying
  it into `Q`, and the merge becomes `T = 2^Er Br Qr Tl + Bl Pl Tr`: a shift
  where there was a product. `Q` and every product with it stay as short as
  the odd part, which made the bit-burst three to four times faster.
- **Final division.** `bsplit::quotient` rounds `T` and `B Q` to the
  working precision before dividing, so the division costs `w`-bit operands
  however long the integers are.

**Tails.** The splitting itself is exact. Each series bounds the tail it
drops, and the bound goes into the radius, so the result is correct for any
term count. The count only sets the accuracy.

**Tests.** The splitting against the series summed term by term in exact
rationals, with each factor present or always 1, and with shifts. Each tail
bound against the true tail, from MPFR at 8000 bits.

## The constants

`ball/consts.rs`, `ball/consts/euler.rs`. Each constant is computed at
`prec + 32` bits (`GUARD`), with its tail bounds in the radius, then
[cached](#the-cache).

### π

Chudnovsky's series [[CC88]](references.md#cc88), as in Arb and
y-cruncher. MPFR's `mpfr_const_pi` uses the Brent–Salamin AGM instead, in
the form of Schönhage, Grotefeld and Vetter [[SGV94]](references.md#sgv94).

The series:

```text
1/π = 12 Σ (-1)^k (6k)! (A + B k) / ((3k)! (k!)^3 C^(3k + 3/2)),
A = 13591409,  B = 545140134,  C = 640320.
```

As a `Series`:
- `p(k) = -(6k-5)(2k-1)(6k-1)`, `q(k) = k^3 C^3 / 24`, `a(k) = A + B k`;
- then `S = Σ a(k) Π p/q` and `π = 426880 √10005 / S`;
- each term gives 47.11 bits, `log2(C^3 / 1728)`.

**Tail.** `(6k)! / ((3k)! (k!)^3) = C(6k, 3k) · (3k)! / (k!)^3 <= 2^(6k) 3^(3k)
= 1728^k`, a binomial and a multinomial coefficient, each at most the sum of
its row. With `r = 1728 / C^3 < 2^-47` and `a(n + j) <= (1 + j) a(n)`, the
tail from `n` on is at most `2 a(n) r^n`.

**Proof.** `chudnovsky_tail`.

### log 2

A Machin-like formula, with three `atanh` series:

```text
log 2 = 18 atanh(1/26) - 2 atanh(1/4801) + 8 atanh(1/8749).
```

As `2 atanh(1/m) = log((m + 1)/(m - 1))`, this is the exact identity
`2 = (27/25)^9 (2401/2400)^-1 (4375/4374)^4`. `atanh(1/m)` is the series
`Σ 1/((2k + 1) m^(2k+1))` (`q(0) = m`, `q(k) = m^2`, `b(k) = 2k + 1`), whose
tail beyond `n` terms is at most `2 m^-(2n+1)`.

**Proof.** `atanh_inv`, `atanh_tail`.

### e

`e = Σ 1/k!`: `q(k) = k`, every other factor 1. With `Q = (n - 1)!` from the
splitting, the tail is at most `2 / n! = 2 / (n Q)`.

**Proof.** `e_uncached`, `e_tail`.

### Euler's constant γ

The Brent-McMillan algorithm B3 [[BM80]](references.md#bm80) with the error
bound of Brent and Johansson [[BJ15]](references.md#bj15), as MPFR's
`mpfr_const_euler` [[MPFR]](references.md#mpfr) (CAS-104).

**Base algorithm.** With `u_k = (n^k / k!)^2`, `I_0(2n) = Σ u_k`,
`S_0(2n) = Σ H_k u_k` and the Bessel function `K_0`,

```text
γ = S_0(2n) / I_0(2n) - K_0(2n) / I_0(2n) - log n
  = S_0 / I_0 - (I_0 K_0) / I_0^2 - log n,
```

and the product `I_0(2n) K_0(2n)` comes from its asymptotic expansion,
truncated at its smallest term:

```text
T = (1/4n) Σ_{k < 2n} ((2k)!)^3 / ((k!)^4 (16n)^(2k)).
```

**The three error terms.** Each goes into a radius:
- **The expansion.** `|I_0(2n) K_0(2n) - T| < 7 e^-4n / (4n)` for every
  integer `n >= 1` [[BJ15]](references.md#bj15) (Lemma 4). As
  `I_0(2n) > e^2n / √(4πn)`, that costs at most `7π e^-8n` in γ, so
  `n ≈ w log(2) / 8` gives `w` bits. The code takes
  `n >= (w + 5) · 0.0866434`, for which the paper's complete bound
  `24 e^-8n` is below `2^-w`.
- **The tails of `I_0` and `S_0`, bounded here.** From `k = N >= 4n` on,
  `u_(k+1) / u_k = n^2 / (k+1)^2 < 1/16` and `H_(k+1) / H_k <= 2`. So the
  tails are below `16/15 u_N` and `8/7 H_N u_N`, taken as `2 u_N` and
  `2 H_N u_N` with `H_N <= 1 + bits(N)` and `u_N = n^(2N) / (N^2 Q)`, `Q`
  from the splitting.

The balls for `S_0` and `I_0` are tested separately against their true
sums. The two tails partly cancel in `S_0 / I_0`, so a missing one would not
show in γ alone.

`N = ⌈4.970626 n⌉ + 1`, from `α (log α - 1) = 3`
[[BJ15]](references.md#bj15) (Corollary 2 and Remark 2), where the terms
have fallen `e^-8n` below `I_0`.

**Splitting with harmonic numbers.** `S_0` and `I_0` share one splitting
that also carries `H_k`, the weighted-partial-sum case of
[[HP98]](references.md#hp98). Over terms `[a, b)`, relative to
`u_(a-1)`, it keeps:
- `P = n^(2(b-a))`, `Q = Π k^2`, `D = Π k`;
- `C / D = Σ 1/k`;
- `T / Q = Σ u_k / u_(a-1)`;
- `V / (Q D) = Σ (H_k - H_(a-1)) u_k / u_(a-1)`.

A leaf is `P = T = V = n^2`, `Q = k^2`, `D = k`, `C = 1`, and two halves
merge as

```text
T = Tl Qr + Pl Tr,   V = Dr (Qr Vl + Cl Pl Tr) + Dl Pl Vr,
C = Cl Dr + Cr Dl,   P = Pl Pr,   Q = Ql Qr,   D = Dl Dr.
```

Then `I_0 ≈ 1 + T / Q` and `S_0 ≈ V / (Q D)`. The asymptotic sum for
`I_0 K_0` is a plain `Series`: `p(k) = (2k-1)^3`, `q(k) = k n^2` with `2^5`
as a shift, and `q(0) = n` with shift 2.

**`log n` (CAS-104, ours).** `n` is free, so the code takes the first
integer of the form `2^a 3^b 5^c 7^d` above the bound, which is at most 5%
above it from 512 on, and 1.1% from `2^19`. Then `log n` is an integer
combination of four `atanh` series. With `x_m = atanh(1/m)`:

```text
log 2 = 144 x_251 + 54 x_449 - 38 x_4801 + 62 x_8749
log 3 = 228 x_251 + 86 x_449 - 60 x_4801 + 98 x_8749
log 5 = 334 x_251 + 126 x_449 - 88 x_4801 + 144 x_8749
log 7 = 404 x_251 + 152 x_449 - 106 x_4801 + 174 x_8749
```

Each `2 x_m` is `log((m + 1)/(m - 1))`, a logarithm of
`126/125`, `225/224`, `2401/2400` or `4375/4374`. The matrix of their
exponents of 2, 3, 5 and 7 has determinant −1, and the coefficients are
twice its inverse. A test checks each identity exactly, as a product of
rational powers. At `10^5` digits this takes 45 ms, against 256 ms for a
general `log`, which would have been a third of γ's time.

**Proof.** `tails`, `bessel_error`, `parameters` and the module comment of
`euler.rs`. **Tests.**
- the splitting against exact sums;
- the tails against the true tails;
- the Bessel bound against `I_0(2n) K_0(2n)` from MPFR for `n <= 40` (the
  true error is 0.6–3.3% of it);
- each ball with too small an `n` or too few terms, where each bound in
  turn dominates the radius;
- MPFR at every precision to 300 bits, up to `10^5` bits, and at `10^6`
  bits in the ignored `million_bits`.

### Catalan's constant G

Lima's series [[Lim12]](references.md#lim12) (Theorem 1). Corollary 1 of
that paper is the same series shifted by one index, which it calls
Guillera's:

```text
G = ½ Σ_{k >= 0} (-8)^k (3k + 2) / ((2k + 1)^3 C(2k, k)^3).
```

The ratio of consecutive `h_k = (-8)^k / ((2k+1)^3 C(2k, k)^3)` is
`-k^3 / (2k+1)^3`, so as a `Series`: `p(k) = -k^3`, `q(k) = (2k+1)^3`,
`a(k) = 3k + 2`. That is three bits per term from one sum, with no π or
logarithm.

MPFR's `mpfr_const_catalan` [[MPFR]](references.md#mpfr) uses Ramanujan's
`G = π/8 log(2 + √3) + 3/8 Σ 1/((2k+1)^2 C(2k, k))` instead, which gives two
bits per term and needs π and a logarithm.

**Tail.** `|h_k| <= 8^-k`, since `C(2k, k) >= 4^k / (2k + 1)`: it is the
largest of the `2k + 1` binomials that sum to `4^k`. From `k = n >= 1` on,
the ratio of `|t_k| = (3k + 2)|h_k|` is at most `(3n + 5)/(3n + 2)/8 <= 1/5`.
So the tail is at most `5/4 (3n + 2) 8^-n`.

**Proof.** `catalan_tail`. **Tests.** The tail against the true tail; G from
1 to 60 terms against MPFR (tail-dominated radii); MPFR at every precision
to 300 bits, to `10^5` bits, and at `10^6` bits.

## The cache

Each constant has one process-wide cache, an `RwLock` holding the ball at
the highest precision computed so far:
- **Lower requests.** A request at or below that precision rounds the
  cached ball, which adds half an ulp to the radius.
- **Higher requests.** A request above it recomputes at
  `max(prec, 2 · cached)`, so growing requests cost a geometric series.
- **Locking.** The computation runs without the lock held.
- **Freeing.** `clear_constant_cache` frees every constant and the function
  tables.

The cache memoises a pure function, so evaluation stays pure; the memory it
holds is bounded by the largest precision requested (CAS-DOC-5 §7.1).
