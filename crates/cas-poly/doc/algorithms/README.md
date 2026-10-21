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
- [references.md](references.md): the bibliography. Chapters and doc
  comments cite it by key, as in "Kronecker substitution [Sch82]".
