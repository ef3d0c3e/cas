//! Whole-interner checks. One test per process: nothing else may intern
//! concurrently, so the global counts are exact.

use cas_core::{Expr, Symbol, interner_stats};

#[test]
fn interner_empties_and_shrinks() {
    assert_eq!(interner_stats().entries, 0);
    let (f, g, x) = (
        Expr::from(Symbol::new("f")),
        Expr::from(Symbol::new("g")),
        Expr::from(Symbol::new("x")),
    );
    let terms: Vec<Expr> = (0..200_000)
        .map(|i| {
            let inner = Expr::apply(&g, &[Expr::from(i), x.clone()]);
            Expr::apply(&f, &[inner, Expr::string("s")])
        })
        .collect();
    // 200 000 outer and inner nodes each, plus the one shared string.
    assert_eq!(interner_stats().entries, 400_001);
    let sum = Expr::apply(&f, &terms);
    drop(terms);
    assert_eq!(interner_stats().entries, 400_002);
    drop(sum);
    let stats = interner_stats();
    assert_eq!(stats.entries, 0);
    // Without shrinking this would stay above 400 000.
    assert!(
        stats.capacity <= 64 * 64,
        "capacity {} after emptying",
        stats.capacity
    );
}
