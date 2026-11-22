---
name: autonomous-areas
description: Move from one CAS area to the next in my recommended order without asking; ask only when a brand new Lific module is needed
metadata:
  node_type: memory
  type: feedback
  originSessionId: 1b3aa32d-02f2-407f-9868-a62b20c7002c
  modified: 2026-10-28T11:28:42.519Z
---

When an area (CAS-PLAN-n) is done, go straight on to the next one in my recommended order instead of asking which area comes next. Ask the user only when the work needs a brand new module (a Lific module that does not exist yet); the existing ones are Algebra, Balls, Calculus, Evaluator, Expressions, Infrastructure, Integers, Language, Polynomials, and Statistics (approved 2026-10-28).

**Why:** after Laplace transforms the user answered "Go with your recommended order. Don't notify me every time, you can keep working until a brand new module needs work."

**How to apply:** at an area's end, close its tickets and doc as usual, then write the next design doc and plan without an AskUserQuestion; give short progress notes, not area-choice questions. The order recommended on 2026-10-28: inequalities (Algebra), systems of ODEs (Calculus), Gröbner bases (Polynomials, with `eliminate`), all done; then statistics. On 2026-10-28 the user chose "Statistics (Recommended)": a new Statistics module, descriptive statistics and regressions first, then distributions (which need incomplete gamma and beta). See [[work-on-master]] and [[ci-on-github]].
