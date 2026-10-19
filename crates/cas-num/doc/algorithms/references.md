# References

The bibliography of the [algorithms documentation](README.md). The other
chapters and the doc comments in the code cite these by key, as in
"Chudnovsky's series [CC88]". Only works an algorithm here is taken from,
or is checked against, are listed.

## Books

- <a id="bz10"></a>**[BZ10]** Richard P. Brent and Paul Zimmermann, *Modern
  Computer Arithmetic*, Cambridge Monographs on Applied and Computational
  Mathematics 18, Cambridge University Press, 2010. Chapter 1 (integer
  arithmetic), chapter 3 (floating-point arithmetic: rounding, Newton's
  reciprocal and reciprocal square root) and chapter 4 (elementary functions:
  argument reduction, rectangular splitting, Newton's method for inverse
  functions, binary splitting and the bit-burst algorithm in §4.9).
- <a id="cw80"></a>**[CW80]** William J. Cody Jr. and William Waite, *Software
  Manual for the Elementary Functions*, Prentice-Hall, 1980.
- <a id="hig02"></a>**[Hig02]** Nicholas J. Higham, *Accuracy and Stability of
  Numerical Algorithms*, 2nd edition, SIAM, 2002. The `γ_n = n u / (1 - n u)`
  bounds for Horner's rule.
- <a id="hw08"></a>**[HW08]** G. H. Hardy and E. M. Wright, *An Introduction to
  the Theory of Numbers*, 6th edition (revised by D. R. Heath-Brown and J. H.
  Silverman), Oxford University Press, 2008. Chapter XII: the Gaussian integers
  are Euclidean.
- <a id="knu97"></a>**[Knu97]** Donald E. Knuth, *The Art of Computer
  Programming*, volume 2, *Seminumerical Algorithms*, 3rd edition,
  Addison-Wesley, 1997. §4.2.2 (the two-sum), §4.3.1 (Algorithm D).
- <a id="moo66"></a>**[Moo66]** Ramon E. Moore, *Interval Analysis*,
  Prentice-Hall, 1966.
- <a id="mul16"></a>**[Mul16]** Jean-Michel Muller, *Elementary Functions:
  Algorithms and Implementation*, 3rd edition, Birkhäuser, 2016.
- <a id="mul18"></a>**[Mul18]** Jean-Michel Muller, Nicolas Brunie, Florent de
  Dinechin, Claude-Pierre Jeannerod, Mioara Joldes, Vincent Lefèvre, Guillaume
  Melquiond, Nathalie Revol and Serge Torres, *Handbook of Floating-Point
  Arithmetic*, 2nd edition, Birkhäuser, 2018.
- <a id="sgv94"></a>**[SGV94]** Arnold Schönhage, Andreas F. W. Grotefeld and
  Ekkehart Vetter, *Fast Algorithms: A Multitape Turing Machine
  Implementation*, BI-Wissenschaftsverlag, Mannheim, 1994. §8.2.5, the form of
  the Brent–Salamin AGM that MPFR's `mpfr_const_pi` uses.

## Papers and reports

- <a id="bj15"></a>**[BJ15]** Richard P. Brent and Fredrik Johansson, "A bound
  for the error term in the Brent-McMillan algorithm", *Mathematics of
  Computation* 84 (2015), 2351–2359, doi:10.1090/S0025-5718-2015-02931-7;
  arXiv:1312.0039. Lemma 4 bounds the asymptotic series of `I_0 K_0`; Theorem
  1, Corollary 2 and Remark 2 give the complete error bound and the term
  counts.
- <a id="bm80"></a>**[BM80]** Richard P. Brent and Edwin M. McMillan, "Some new
  algorithms for high-precision computation of Euler's constant", *Mathematics
  of Computation* 34 (149) (1980), 305–312. Algorithm B3.
- <a id="bre76"></a>**[Bre76]** Richard P. Brent, "The complexity of
  multiple-precision arithmetic", in R. S. Anderssen and R. P. Brent (eds.),
  *The Complexity of Computational Problem Solving*, University of Queensland
  Press, 1976, 126–165. The bit-burst algorithm.
- <a id="cc88"></a>**[CC88]** David V. Chudnovsky and Gregory V. Chudnovsky,
  "Approximation and complex multiplication according to Ramanujan", in
  *Ramanujan Revisited*, Academic Press, 1988, 375–472. The series for `1/π`.
- <a id="dek71"></a>**[Dek71]** T. J. Dekker, "A floating-point technique for
  extending the available precision", *Numerische Mathematik* 18 (1971),
  224–242. Fast two-sum, the splitting of a double into halves (credited there
  to Veltkamp) and the exact product from it.
- <a id="fou07"></a>**[Fou07]** Laurent Fousse, Guillaume Hanrot, Vincent
  Lefèvre, Patrick Pélissier and Paul Zimmermann, "MPFR: A multiple-precision
  binary floating-point library with correct rounding", *ACM Transactions on
  Mathematical Software* 33 (2) (2007), article 13.
- <a id="hft97"></a>**[HFT97]** T. E. Hull, Thomas F. Fairgrieve and Ping Tak
  Peter Tang, "Implementing the complex arcsine and arccosine functions using
  exception handling", *ACM Transactions on Mathematical Software* 23 (3)
  (1997), 299–335.
- <a id="hoe09"></a>**[Hoe09]** Joris van der Hoeven, "Ball arithmetic",
  technical report, 2009, HAL hal-00432152.
- <a id="hp98"></a>**[HP98]** Bruno Haible and Thomas Papanikolaou, "Fast
  multiprecision evaluation of series of rational numbers", in *Algorithmic
  Number Theory (ANTS-III)*, Lecture Notes in Computer Science 1423, Springer,
  1998, 338–350, doi:10.1007/BFb0054873. Binary splitting, including sums
  weighted by partial sums such as harmonic numbers.
- <a id="hz04"></a>**[HZ04]** Guillaume Hanrot and Paul Zimmermann, "A long
  note on Mulders' short product", *Journal of Symbolic Computation* 37 (3)
  (2004), 391–401.
- <a id="jmp17"></a>**[JMP17]** Mioara Joldes, Jean-Michel Muller and Valentina
  Popescu, "Tight and rigorous error bounds for basic building blocks of
  double-word arithmetic", *ACM Transactions on Mathematical Software* 44 (2)
  (2017), doi:10.1145/3121432. Algorithms 6, 7, 10 and 17 and their bounds.
- <a id="joh15"></a>**[Joh15]** Fredrik Johansson, "Efficient implementation of
  elementary functions in the medium-precision range", in *Proceedings of the
  22nd IEEE Symposium on Computer Arithmetic (ARITH 22)*, 2015, 83–89,
  doi:10.1109/ARITH.2015.16; arXiv:1410.7176. Table-based argument reduction
  and fixed-point series.
- <a id="joh17"></a>**[Joh17]** Fredrik Johansson, "Arb: efficient
  arbitrary-precision midpoint-radius interval arithmetic", *IEEE Transactions
  on Computers* 66 (8) (2017), 1281–1292; arXiv:1611.02831. The `mag_t` radius
  type, radius propagation, and functions evaluated at the midpoint with a
  derivative bound.
- <a id="kah87"></a>**[Kah87]** William Kahan, "Branch cuts for complex
  elementary functions, or much ado about nothing's sign bit", in A. Iserles
  and M. J. D. Powell (eds.), *The State of the Art in Numerical Analysis*,
  Clarendon Press, 1987, 165–211. Principal branches, the side a signed zero
  selects on a cut, and the complex square root used here.
- <a id="lim12"></a>**[Lim12]** F. M. S. Lima, "A rapidly converging
  Ramanujan-type series for Catalan's constant", 2012, arXiv:1207.3139. Theorem
  1; Corollary 1 is the same series shifted by one index, which the paper calls
  Guillera's.
- <a id="mah53"></a>**[Mah53]** Kurt Mahler, "On the approximation of π",
  *Nederl. Akad. Wetensch. Proc. Ser. A* 56 = *Indagationes Mathematicae* 15
  (1953), 30–42. π has a finite irrationality measure, which bounds how close a
  dyadic number comes to a multiple of π/2.
- <a id="mg11"></a>**[MG11]** Niels Möller and Torbjörn Granlund, "Improved
  division by invariant integers", *IEEE Transactions on Computers* 60 (2)
  (2011), 165–175, doi:10.1109/TC.2010.143. The 3-by-2 quotient step of the
  schoolbook division used by the fixed-size float division.
- <a id="mul00"></a>**[Mul00]** Thom Mulders, "On short multiplications and
  divisions", *Applicable Algebra in Engineering, Communication and Computing*
  11 (1) (2000), 69–88, doi:10.1007/s002000000037.
- <a id="ps73"></a>**[PS73]** Michael S. Paterson and Larry J. Stockmeyer, "On
  the number of nonscalar multiplications necessary to evaluate polynomials",
  *SIAM Journal on Computing* 2 (1) (1973), 60–66.
- <a id="rum99"></a>**[Rum99]** Siegfried M. Rump, "Fast and parallel interval
  arithmetic", *BIT Numerical Mathematics* 39 (3) (1999),
  doi:10.1023/A:1022374804152. Midpoint-radius arithmetic in floating point.
- <a id="rzbm09"></a>**[RZBM09]** Siegfried M. Rump, Paul Zimmermann, Sylvie
  Boldo and Guillaume Melquiond, "Computing predecessor and successor in
  rounding to nearest", *BIT Numerical Mathematics* 49 (2) (2009), 419–431,
  doi:10.1007/s10543-009-0218-z.
- <a id="smi89"></a>**[Smi89]** David M. Smith, "Efficient multiple-precision
  evaluation of elementary functions", *Mathematics of Computation* 52 (185)
  (1989), 131–134. Concurrent (rectangular) summation of Taylor series.
- <a id="zim99"></a>**[Zim99]** Paul Zimmermann, "Karatsuba square root", INRIA
  Research Report RR-3805, 1999; Algorithm 1.12 of [BZ10]. The integer square
  root (`Natural::sqrt_rem`) behind the exact float square root and the
  approximate one.
- <a id="ziv91"></a>**[Ziv91]** Abraham Ziv, "Fast evaluation of elementary
  mathematical functions with correctly rounded last bit", *ACM Transactions on
  Mathematical Software* 17 (3) (1991), 410–423. Retrying with more precision
  until the result is decided.

## Libraries, manuals and standards

- <a id="arb"></a>**[Arb]** Fredrik Johansson et al., *Arb, a C library for
  arbitrary-precision ball arithmetic*, https://arblib.org, part of FLINT
  (https://flintlib.org) since version 3.0. Followed here for `arb_printn`
  (decimal printing), `mag_t`, and the bit-burst `arb_exp_arf_bb` and
  `arb_atan_arf_bb`.
- <a id="c99"></a>**[C99]** ISO/IEC 9899:1999, *Programming languages — C*: the
  `%a` hexadecimal float format (§7.19.6.1) and Annex G, the branch cuts and
  signed zeros of the complex functions.
- <a id="fdlibm"></a>**[fdlibm]** Sun Microsystems, *FDLIBM 5.3*, the freely
  distributable C math library (https://www.netlib.org/fdlibm/). The
  cancellation-free forms of `sinh`, `tanh`, `asinh`, `acosh` and `atanh`
  through `expm1` and `log1p` (`e_sinh.c`, `s_tanh.c`, `s_asinh.c`,
  `e_acosh.c`, `e_atanh.c`).
- <a id="gmp"></a>**[GMP]** Torbjörn Granlund and the GMP development team,
  *GNU MP: The GNU Multiple Precision Arithmetic Library*, manual for version
  6.3.0, chapter "Algorithms". The source's `mpn_invertappr` takes its Newton
  residual as a product modulo `B^m - 1`, the wrap-around trick used here.
- <a id="mpc"></a>**[MPC]** Andreas Enge, Philippe Théveny and Paul Zimmermann,
  *GNU MPC: The GNU Multiple Precision Complex Library*, manual for version
  1.4.1. The oracle for complex balls, and the reference for the side of a
  branch cut a point on it takes.
- <a id="mpfr"></a>**[MPFR]** The MPFR team, *The MPFR Library: Algorithms and
  Proofs*, https://www.mpfr.org/algorithms.pdf; and the MPFR 4.2.2 sources
  (`const_euler.c`, `const_catalan.c`). The oracle for real floats and balls.
