# Algorithms of cas-poly

This is the algorithms documentation of the `cas-poly` crate, written as
cas-num's (`crates/cas-num/doc/algorithms`), after GMP's "Algorithms" chapter
and MPFR's *Algorithms and Proofs*. For each operation it gives:

- what it computes;
- the base algorithm, and where it comes from;
- the tricks on top of it, and why each is there (ours are marked with their
  ticket);
- how it is tested.

The design that chose these algorithms is CAS-DOC-8 in Lific, whose §16 lists
the tickets. This documentation describes the code as it is.

## Chapters

- [dense.md](dense.md): dense univariate polynomials (CAS-135). The
  representation; products by schoolbook, Karatsuba's method and Kronecker
  substitution, with their measured thresholds; division over a field,
  pseudo-division and exact division; evaluation, derivative, content and
  primitive part; ℚ and ℚ(i) over a shared denominator.
- [sparse.md](sparse.md): sparse multivariate polynomials (CAS-136). Packed
  exponent vectors with guard bits; sums by merging; Johnson's heap for
  products and Monagan and Pearce's for division, which widens its fields
  when the guards report an overflow; evaluation, substitution, derivative,
  the recursive view; the benchmarks of CAS-DOC-8 §13.
- [gcd.md](gcd.md): greatest common divisors, resultants and discriminants
  (CAS-137). The subresultant reference; the heuristic gcd, with Parisse's
  theorem behind its check; Brown's modular algorithm over ℤ and its dense
  evaluation over prime fields; resultants modulo primes within Hadamard's
  bound; the domains and their choices; the benchmarks against giac.
- [factor.md](factor.md): factoring (CAS-142, planned in CAS-PLAN-1). The
  constant-and-powers form of a factorisation; square-free decomposition by
  Yun's algorithm, in several variables through contents; factors of
  square-free polynomials over prime fields, by distinct degrees and Cantor
  and Zassenhaus's equal-degree splitting; multifactor Hensel lifting;
  univariate factoring over ℤ and ℚ by Zassenhaus's algorithm;
  multivariate factoring over ℤ and ℚ by Wang's algorithm, with its
  Hensel lifting; factoring over ℚ(i) by Trager's norms; van Hoeij's
  recombination by integral LLL.
- [integrate.md](integrate.md): the integration of rational functions over
  ℚ (CAS-196, planned in CAS-PLAN-3). The polynomial part; Hermite's
  reduction in Mack's linear form; the logarithms with rational residues
  from the factors of the residue resultant; Sturm's count of real roots,
  for the absolute values of the logarithms.
- [references.md](references.md): the bibliography. Chapters and doc
  comments cite it by key, as in "Kronecker substitution [Sch82]".
