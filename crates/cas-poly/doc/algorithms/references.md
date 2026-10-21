# References

The bibliography of the [algorithms documentation](README.md). The other
chapters and the doc comments in the code cite these by key. Only works an
algorithm here is taken from, or is compared with, are listed.

## Books

- <a id="knu97"></a>**[Knu97]** Donald E. Knuth, *The Art of Computer
  Programming*, volume 2, *Seminumerical Algorithms*, 3rd edition,
  Addison-Wesley, 1997. §4.6.1, "Division of Polynomials": Algorithm D
  (division of polynomials over a field) and Algorithm R (pseudo-division,
  `l(v)^(m−n+1) u(x) = q(x) v(x) + r(x)` with `deg r < n`, without
  division).

## Papers

- <a id="har09"></a>**[Har09]** David Harvey, "Faster polynomial
  multiplication via multipoint Kronecker substitution", *Journal of
  Symbolic Computation* 44 (10) (2009), 1502–1510; arXiv:0712.4046.
  Evaluates at several points (`2^N`, `−2^N`, `2^−N`, …) to replace one large
  integer product by several smaller ones; not used yet.
- <a id="ko62"></a>**[KO62]** Anatolii Karatsuba and Yuri Ofman,
  "Multiplication of multidigit numbers on automata", *Doklady Akademii Nauk
  SSSR* 145 (2) (1962), 293–294; English translation in *Soviet Physics
  Doklady* 7 (1963), 595–596. The three-product identity.
- <a id="kro82"></a>**[Kro82]** Leopold Kronecker, "Grundzüge einer
  arithmetischen Theorie der algebraischen Grössen", *Journal für die reine
  und angewandte Mathematik* 92 (1882), 1–122. The substitution of powers of
  one variable for the others, which reduces multivariate problems to
  univariate ones.
- <a id="sch82"></a>**[Sch82]** Arnold Schönhage, "Asymptotically fast
  algorithms for the numerical multiplication and division of polynomials
  with complex coefficients", in J. Calmet (ed.), *Computer Algebra
  (EUROCAM '82)*, Lecture Notes in Computer Science 144, Springer, 1982,
  3–15, doi:10.1007/3-540-11607-9_1. The reduction of products in ℤ\[x\] to
  products in ℤ, by evaluation at a power of 2 ("binary segmentation"), as
  [Har09] credits it.

## Libraries

- <a id="flint"></a>**[FLINT]** The FLINT team, *FLINT: Fast Library for
  Number Theory*, documentation of `fmpq_poly.h`,
  <https://flintlib.org/doc/fmpq_poly.html>. A polynomial over ℚ as an
  integer polynomial over one positive denominator, coprime to it, with zero
  as `0/1`.
