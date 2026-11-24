//! Tests of the simplex method (CAS-417, CAS-603): a classic program,
//! equations, infeasible and unbounded ones, Beale's, which cycles without
//! Bland's rule, and the work of a call: its iterations, the user's limit,
//! the interrupt and the cells.

use std::sync::atomic::AtomicBool;

use cas_num::{Integer, Rational};

use super::{Outcome, Relation, Row, Work};

/// The outcome of `minimize` with no limit.
fn minimize(c: &[Rational], rows: &[Row]) -> Outcome {
    let stop = AtomicBool::new(false);
    super::minimize(c, rows, &mut Work::new(0, usize::MAX, &stop)).0
}

fn q(n: i64, d: i64) -> Rational {
    Rational::new(Integer::from(n), Integer::from(d))
}

fn ints(v: &[i64]) -> Vec<Rational> {
    v.iter().map(|&n| q(n, 1)).collect()
}

fn row(a: &[i64], rel: Relation, b: i64) -> Row {
    (ints(a), rel, q(b, 1))
}

#[test]
fn programs() {
    // max 3x + 5y on x ≤ 4, 2y ≤ 12, 3x + 2y ≤ 18: (2, 6), 36.
    let rows = [
        row(&[1, 0], Relation::AtMost, 4),
        row(&[0, 2], Relation::AtMost, 12),
        row(&[3, 2], Relation::AtMost, 18),
    ];
    assert_eq!(
        minimize(&ints(&[-3, -5]), &rows),
        Outcome::Optimal(ints(&[2, 6]), q(-36, 1))
    );
    // Equations, and a negative right side.
    let rows = [
        row(&[1, 1], Relation::Equal, 2),
        row(&[1, -1], Relation::Equal, 0),
        row(&[-1, 0], Relation::AtMost, -1),
    ];
    assert_eq!(
        minimize(&ints(&[1, 1]), &rows),
        Outcome::Optimal(ints(&[1, 1]), q(2, 1))
    );
    // x ≥ 2 and x ≤ 1: no point.
    let rows = [
        row(&[1], Relation::AtLeast, 2),
        row(&[1], Relation::AtMost, 1),
    ];
    assert_eq!(minimize(&ints(&[1]), &rows), Outcome::Infeasible);
    // min −x on x ≥ 1.
    let rows = [row(&[1], Relation::AtLeast, 1)];
    assert_eq!(minimize(&ints(&[-1]), &rows), Outcome::Unbounded);
}

/// Beale's program cycles under Dantzig's rule; Bland's ends at
/// `(1, 0, 1, 0)`, `−5/4`.
#[test]
fn beale() {
    let c = vec![q(-3, 4), q(20, 1), q(-1, 2), q(6, 1)];
    let rows = [
        (
            vec![q(1, 4), q(-8, 1), q(-1, 1), q(9, 1)],
            Relation::AtMost,
            q(0, 1),
        ),
        (
            vec![q(1, 2), q(-12, 1), q(-1, 2), q(3, 1)],
            Relation::AtMost,
            q(0, 1),
        ),
        (ints(&[0, 0, 1, 0]), Relation::AtMost, q(1, 1)),
    ];
    assert_eq!(
        minimize(&c, &rows),
        Outcome::Optimal(ints(&[1, 0, 1, 0]), q(-5, 4))
    );
}

/// giac's `icount`: each search for an entering column counts, the last
/// one too; past the user's limit the method stops, on an interrupt too,
/// and a tableau past the cells gives up before any work.
#[test]
fn work() {
    let rows = [
        row(&[1, 1], Relation::Equal, 2),
        row(&[1, -1], Relation::Equal, 0),
        row(&[-1, 0], Relation::AtMost, -1),
    ];
    let c = ints(&[1, 1]);
    let stop = AtomicBool::new(false);
    let mut w = Work::new(0, usize::MAX, &stop);
    assert!(matches!(
        super::minimize(&c, &rows, &mut w).0,
        Outcome::Optimal(..)
    ));
    assert_eq!(w.iterations, 4);
    let mut w = Work::new(4, usize::MAX, &stop);
    assert!(matches!(
        super::minimize(&c, &rows, &mut w).0,
        Outcome::Optimal(..)
    ));
    let mut w = Work::new(3, usize::MAX, &stop);
    assert_eq!(super::minimize(&c, &rows, &mut w).0, Outcome::Limit);
    assert_eq!(w.iterations, 4);
    let mut w = Work::new(0, 8, &stop);
    assert_eq!(super::minimize(&c, &rows, &mut w).0, Outcome::GaveUp);
    assert_eq!(w.iterations, 0);
    let stop = AtomicBool::new(true);
    let mut w = Work::new(0, usize::MAX, &stop);
    assert_eq!(super::minimize(&c, &rows, &mut w).0, Outcome::Interrupted);
}
