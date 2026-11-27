# Real roots: Sturm sequences, VAS and LMQ bounds

`src/real_roots.rs` (CAS-538) and `DenseRing<Rationals>::sturm_sequence` in
`src/sturm.rs` (CAS-539), behind cas-core's `RealRootIsolation`,
`PositiveRootIsolation`, `PositiveRootUpperBound`, `PositiveRootLowerBound`
and `SturmSequence` (giac's `VAS`, `VAS_positive`, `posubLMQ`, `poslbdLMQ`,
`sturmseq` and `sturm`). The algorithms are ported from giac's
`csturm.cc` and `poly.h` (GPL-3, as this crate), in this crate's own code,
so that the answers are giac's.

## What is computed

- `sturm_sequence(p)`: giac's Sturm sequence of `p` over ℚ, from
  subresultants, integral when `p` is.
- `posub_lmq(p)`, `poslbd_lmq(p)`: an upper and a lower bound on the
  positive roots, rounded as giac rounds them: an integer, and 0, a power
  of two or `+∞` (`LowerBound`).
- `vas_positive(p)`, `vas(p)`: the positive (all the real) roots of the
  square-free part of `p`, each an `Isolated::Root` when the search meets
  it exactly, otherwise an `Isolated::Interval(a, b)` holding exactly one
  root in `(a, b)`; sorted by left end.

## Sturm sequences [Stu35; BT71, Bro78]

`count_real_roots` (CAS-198, [integrate.md](integrate.md#real-roots)) uses
the classical sequence: `p`, `p′`, then the negated remainders over ℚ. giac
prints another, integral one, and `sturmseq` shows it, so
`sturm_sequence` follows giac's `Tsturm_seq`. With `a` the primitive part of
`p` and `b = p′`, each step takes the pseudo-remainder `r = b₀^(δ+1) a mod
b` (`b₀` the leading coefficient of `b`, `δ = deg a − deg b`) and divides
it by `−g h^δ`, the subresultant PRS's divisor [Bro78], with `g` and `h`
starting at 1, then `g = b₀` and `h = b₀` (`δ = 1`) or `b₀^δ h / h^δ`. The
division keeps the entries integral but not their signs: Sturm's theorem
wants `−(a mod b)`, and dividing by `−g h^δ` and multiplying by `b₀^(δ+1)`
multiplies it by a factor whose sign is that of `h g` for odd `δ` and of
`b₀ g` for even `δ`. giac corrects with two running products, one for the
remainders of even index and one for the odd ones, since an error on a
remainder carries to the one two places on: when `δ` is odd and `h ≠ g`
the slot is multiplied by `h g`, when `δ` is even and `b₀ ≠ g` by `b₀ g`,
and the remainder pushed is `b` times its slot. `x³ − 7x + 7` gives `[x³ −
7x + 7, 3x² − 7, 42x − 63, 49]`; `x⁵ − 3x² + x − 1` ends with
`1317150x − 327150` and `−781633800`.

cas-core's `SturmSequence[p]` is giac's `sturm` wrapper: `p` as `n/d` over
ℤ (`d`'s leading coefficient positive), and for each of `n` and `d` a
constant followed by the sequence of each square-free factor of odd
multiplicity (even ones change no sign). The factors are primitive with a
positive leading coefficient; a negative primitive part puts its sign on
the factor of highest odd multiplicity, or is a factor `−1` of its own when
there is none. As giac takes `x^v` out before its square-free
factorisation (CAS-661), a factor that is `x` alone does not count: it
takes the sign only when `v = 1` and no other factor of odd multiplicity is
left (`sturmseq(-x^3)` keeps `−1` apart, `sturmseq(-x)` does not); the constant is the leading coefficient over the factors'.
With bounds `a ≤ b` the count is `Σ V(a) − V(b)` over the sequences, the
roots and poles in `(a, b]`; with none giac answers −1 when the product of
the first signs at `a` (and the constants' signs) is negative, 0 otherwise.

## VAS: continued fractions [Vin36; CA76; AS05]

Vincent's theorem: repeatedly substituting `x ↦ q + 1/x` with the partial
quotients `q` of a root brings any square-free polynomial to at most one
sign variation, and by Descartes's rule of signs 0 variations means no
positive root, 1 exactly one. VAS keeps the composition of these maps as a
Möbius map `x ↦ (a x + b)/(c x + d)` from the half line onto an interval
between `b/d` and `a/c` (to `+∞` when `c = 0`), with `f` the polynomial in
the new variable; giac's `VAS_positive_roots` numbers the steps of [AS05]:

1. `v = var(p)`: none, no root; one, the interval `[0, posubLMQ(p)]`.
2. Take a map from the work list (last in, first out).
3. `lb = poslbdLMQ(f)`, a lower bound on `f`'s positive roots.
4. If `lb > 16`, `f(x) ← f(lb x)` and `a, c ← lb a, lb c`; `lb ← 1`.
5. If `lb ≥ 1`, `f(x) ← f(x + lb)` (`b, d ← lb a + b, lb c + d`): a root at
   0 now is `b/d`, given exactly and divided out; 0 variations drop the
   map, one gives its interval.
6. `f₁ = f(x + 1)` for `(1, ∞)`, with `b/d` at 0 given exactly; `v₂ = v −
   var(f₁) − [root at 1]` bounds the roots in `(0, 1)`.
7. If `v₂ > 1`, `f₂ = (x + 1)^n f(1/(x + 1))` and `v₂ = var(f₂)`.
8. The branch with more variations is handled first.
9. and 10. Each branch: 0, dropped; 1, its interval, from `b/d` to `a/c`,
   or from `b` to `b + a·posubLMQ(f)` when `c = 0`; more, back on the list.

The entries are sorted by left end, a root before an interval starting at
it (`VAS((x − 1)²(x − 2)) = [1, [1, 4]]`). `vas` runs it on the square-free
part `v` of `p` and on `w = ±v(−x)` with the map `x ↦ −x` for the negative
roots, or, when `w = v` (an even or odd `p`), mirrors the positive ones;
the root 0 is divided out first and put between them. Termination is
Vincent's theorem on a square-free polynomial: `vas_positive` takes the
square-free part too, where giac's `VAS_positive` takes `p` as it is
(a repeated root gives a spurious interval, `[[0,1],1,[1,4]]` for `(x −
1)²(x − 2)`). cas-core leaves polynomials above degree 64.

## LMQ bounds [ASV08]

Local-Max-Quadratic: with the coefficients `cᵢ` of `xⁱ` and the leading
one made positive, each negative `c_m` is paired with every positive `c_n`,
`n > m`, at `(2^t |c_m| / c_n)^(1/(n − m))`, `t` counting the uses of `c_n`
from 1; the least over `n`, then the largest over `m`, bounds the positive
roots. giac's `posubLMQ` returns `⌈65 B / 64⌉`, a small margin, and 0 with
no negative coefficient. `poslbdLMQ` applies the same to the reversed
polynomial (the constant made positive), a bound `B` on `1/x`, and returns
`2^(−⌈log₂ B⌉)`: `+∞` with no negative coefficient, 0 when a negative
coefficient has no positive one below it.

giac works on the logarithms in doubles. Here the candidates `q^(1/k)` are
compared exactly as `q₁^k₂` against `q₂^k₁`; `⌈65 B/64⌉` is the least `n`
with `(64n/65)^k ≥ q`, found from `⌊q^(1/k)⌋` (an integer root) within two
steps; `⌈log₂ B⌉` is the least `e` with `2^(e k) ≥ q`, from the bit lengths
within three. A bound exactly on a rounding edge may thus differ from
giac's by one step where giac's doubles fall on the other side; none of
giac's examples do (`posubLMQ(x − 2) = 5`, `poslbdLMQ(x³ − 7x + 7) = 1/2`).

## Tests

`src/real_roots/tests.rs`: giac's `posubLMQ` and `poslbdLMQ` answers on
twenty-odd polynomials, giac's `VAS` and `VAS_positive` intervals on its
examples (exact roots met by the search, the root 0, even polynomials,
roots scaled past 16), and, for eleven polynomials (Chebyshev's `T₈`,
Wilkinson's to 8, Mignotte's `x¹⁰ − 2(10x − 1)²`, …), that each interval
holds exactly one root in its interior by `count_real_roots` (closed count
less the ends that are roots), each exact entry is a root, and the entries
count the real roots. `src/sturm/tests.rs`: giac's `sturmseq` sequences,
signs corrected, of degree up to 6.

cas-core's `native/polytools/tests.rs` and cas-giac's `tests/real_roots.rs`
check the builtins against giac's answers.
