# References

The bibliography of the [algorithms documentation](README.md). The other
chapters and the doc comments in the code cite these by key. Only works an
algorithm here is taken from, or is compared with, are listed.

## Books

- <a id="coh93"></a>**[Coh93]** Henri Cohen, *A Course in Computational
  Algebraic Number Theory*, Graduate Texts in Mathematics 138, Springer,
  1993. §3.3, "The Sub-Resultant Algorithm": §3.3.1, "Description of the
  Algorithm", the gcd of the contents times the primitive part of the last
  remainder; §3.3.2, "Resultants and Discriminants".
- <a id="gcl92"></a>**[GCL92]** Keith O. Geddes, Stephen R. Czapor and
  George Labahn, *Algorithms for Computer Algebra*, Kluwer, 1992. Chapter 6,
  "Newton's iteration and the Hensel construction": §6.7–6.8, the
  multivariate Hensel construction, with Algorithm 6.2, the multivariate
  diophantine equation, and Algorithm 6.4, multivariate Hensel lifting with
  the leading coefficients given.
- <a id="knu97"></a>**[Knu97]** Donald E. Knuth, *The Art of Computer
  Programming*, volume 2, *Seminumerical Algorithms*, 3rd edition,
  Addison-Wesley, 1997. §4.6.1, "Division of Polynomials": Algorithm D
  (division of polynomials over a field) and Algorithm R (pseudo-division,
  `l(v)^(m−n+1) u(x) = q(x) v(x) + r(x)` with `deg r < n`, without
  division).
- <a id="vzgg13"></a>**[vzGG13]** Joachim von zur Gathen and Jürgen
  Gerhard, *Modern Computer Algebra*, 3rd edition, Cambridge University
  Press, 2013. Chapter 14, "Factoring polynomials over finite fields": §14.2,
  distinct-degree factorisation; §14.3, equal-degree factorisation by Cantor
  and Zassenhaus's method. Chapter 15, "Hensel lifting and factoring
  polynomials": Algorithm 15.10, the Hensel step from `m` to `m²`, and
  §15.5, multifactor Hensel lifting through a factor tree.

## Papers

- <a id="bro71"></a>**[Bro71]** W. S. Brown, "On Euclid's algorithm and the
  computation of polynomial greatest common divisors", *Journal of the ACM*
  18 (4) (1971), 478–504. The modular gcd: images modulo primes and at
  evaluation points, unlucky ones recognised by their degrees.
- <a id="bro78"></a>**[Bro78]** W. S. Brown, "The subresultant PRS
  algorithm", *ACM Transactions on Mathematical Software* 4 (3) (1978),
  237–249.
- <a id="bt71"></a>**[BT71]** W. S. Brown and J. F. Traub, "On Euclid's
  algorithm and the theory of subresultants", *Journal of the ACM* 18 (4)
  (1971), 505–514. The subresultant remainder sequence, and why its
  divisions are exact.
- <a id="cgg89"></a>**[CGG89]** B. W. Char, K. O. Geddes and G. H. Gonnet,
  "GCDHEU: Heuristic polynomial GCD algorithm based on integer GCD
  computation", *Journal of Symbolic Computation* 7 (1989), 31–48; first in
  *EUROSAM 84*, Lecture Notes in Computer Science 174, Springer, 1984. The
  heuristic gcd.
- <a id="cz81"></a>**[CZ81]** David G. Cantor and Hans Zassenhaus, "A new
  algorithm for factoring polynomials over finite fields", *Mathematics of
  Computation* 36 (154) (1981), 587–592. Equal-degree splitting by
  `gcd(a^((qᵈ − 1)/2) − 1, g)` for a random `a`.
- <a id="fat03"></a>**[Fat03]** Richard Fateman, "Comparing the speed of
  programs for sparse polynomial multiplication", *ACM SIGSAM Bulletin* 37
  (1) (2003), 4–15, doi:10.1145/844076.844080. The benchmark `f (f + 1)` for
  `f = (1 + x + y + z + t)^20`.
- <a id="har09"></a>**[Har09]** David Harvey, "Faster polynomial
  multiplication via multipoint Kronecker substitution", *Journal of
  Symbolic Computation* 44 (10) (2009), 1502–1510; arXiv:0712.4046.
  Evaluates at several points (`2^N`, `−2^N`, `2^−N`, …) to replace one large
  integer product by several smaller ones; not used yet.
- <a id="joh74"></a>**[Joh74]** Stephen C. Johnson, "Sparse polynomial
  arithmetic", *ACM SIGSAM Bulletin* 8 (3) (1974), 63–71. Products and
  division of sparse polynomials by merging the rows through a binary heap.
- <a id="ko62"></a>**[KO62]** Anatolii Karatsuba and Yuri Ofman,
  "Multiplication of multidigit numbers on automata", *Doklady Akademii Nauk
  SSSR* 145 (2) (1962), 293–294; English translation in *Soviet Physics
  Doklady* 7 (1963), 595–596. The three-product identity.
- <a id="kro82"></a>**[Kro82]** Leopold Kronecker, "Grundzüge einer
  arithmetischen Theorie der algebraischen Grössen", *Journal für die reine
  und angewandte Mathematik* 92 (1882), 1–122. The substitution of powers of
  one variable for the others, which reduces multivariate problems to
  univariate ones.
- <a id="mig74"></a>**[Mig74]** Maurice Mignotte, "An inequality about
  factors of polynomials", *Mathematics of Computation* 28 (128) (1974),
  1153–1157. For a factor `g` of degree `m` of `f` over ℤ,
  `‖g‖₁ ≤ 2ᵐ |lc(g)/lc(f)| ‖f‖₂`: the bound that sets how far Hensel lifting
  goes.
- <a id="mp07"></a>**[MP07]** Michael Monagan and Roman Pearce, "Polynomial
  division using dynamic arrays, heaps, and packed exponent vectors", in
  *Computer Algebra in Scientific Computing (CASC 2007)*, Lecture Notes in
  Computer Science 4770, Springer, 2007. Terms in dynamic arrays, exponents
  packed into words with the top bit of each field "a guard bit for monomial
  division", Johnson's heap for products and for division by merging `f`
  with `−qᵢ g`, and a divisor heap.
- <a id="mp11"></a>**[MP11]** Michael Monagan and Roman Pearce, "Sparse
  polynomial division using a heap", *Journal of Symbolic Computation* 46
  (7) (2011), 807–822. Heap division as costly as multiplication,
  fraction-free over ℤ.
- <a id="par05"></a>**[Par05]** Bernard Parisse, "A correct proof of the
  heuristic GCD algorithm", note, Institut Fourier, Grenoble, January 2005;
  an earlier version is arXiv:cs/0206032 (2002). Theorem 1: for integers
  `|z| ≥ 2 min(|P|, |Q|) + 2`, `|P|` the largest coefficient, if the
  primitive part of the `z`-adic symmetric reconstruction of `gcd(P(z),
  Q(z))` divides `P` and `Q`, it is their gcd, in any number of variables,
  over ℤ or ℤ\[i\]. It fills the gap of [CGG89]'s proof in several
  variables.
- <a id="sch82"></a>**[Sch82]** Arnold Schönhage, "Asymptotically fast
  algorithms for the numerical multiplication and division of polynomials
  with complex coefficients", in J. Calmet (ed.), *Computer Algebra
  (EUROCAM '82)*, Lecture Notes in Computer Science 144, Springer, 1982,
  3–15, doi:10.1007/3-540-11607-9_1. The reduction of products in ℤ\[x\] to
  products in ℤ, by evaluation at a power of 2 ("binary segmentation"), as
  [Har09] credits it.
- <a id="slf14"></a>**[SLF14]** Guy L. Steele Jr., Doug Lea and Christine
  H. Flood, "Fast splittable pseudorandom number generators", in *Proceedings
  of the 2014 ACM International Conference on Object Oriented Programming
  Systems Languages & Applications (OOPSLA '14)*, ACM, 2014, 453–472. The
  `SplitMix64` generator of the random splitting polynomials.
- <a id="wan78"></a>**[Wan78]** Paul S. Wang, "An improved multivariate
  polynomial factoring algorithm", *Mathematics of Computation* 32 (144)
  (1978), 1215–1231. Evaluation points, the leading coefficients of the
  factors found before lifting, and lifting one variable at a time.
- <a id="yun76"></a>**[Yun76]** David Y. Y. Yun, "On square-free
  decomposition algorithms", in *Proceedings of the Third ACM Symposium on
  Symbolic and Algebraic Computation (SYMSAC '76)*, ACM, 1976, 26–35. The
  decomposition from one gcd with cofactors a round: `gcd(bₖ, cₖ − bₖ′)`.
- <a id="zas69"></a>**[Zas69]** Hans Zassenhaus, "On Hensel factorization,
  I", *Journal of Number Theory* 1 (3) (1969), 291–311. Factoring over ℤ by
  lifting the factors modulo a prime and recombining them.

## Libraries

- <a id="flint"></a>**[FLINT]** The FLINT team, *FLINT: Fast Library for
  Number Theory*, documentation of `fmpq_poly.h`,
  <https://flintlib.org/doc/fmpq_poly.html>. A polynomial over ℚ as an
  integer polynomial over one positive denominator, coprime to it, with zero
  as `0/1`.
- <a id="sympy"></a>**[SymPy]** SymPy 1.14, `sympy/polys/euclidtools.py`:
  `dup_zz_heu_gcd` and `dmp_zz_heu_gcd`, the heuristic gcd, which tries
  `HEU_GCD_MAX = 6` points, the next from `x` as `73794*x*K.sqrt(K.sqrt(x))
  // 27011`. It cites Liao and Fateman's evaluation of the heuristic (ISSAC
  1995). `sympy/polys/factortools.py`: `dmp_zz_wang`, Wang's algorithm, which
  lifts modulo one prime above Mignotte's bound.
