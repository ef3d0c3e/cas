//! The input form on what lowering makes of arbitrary text: every expression
//! that prints reads back as itself, up to the arithmetic the text writes
//! out, and printing what it reads back is stable.
#![no_main]

use cas_lang::{Names, fold_arithmetic, lower, parse, print};
use libfuzzer_sys::fuzz_target;

/// The expression `text` lowers to: one statement, without errors or notes.
fn read(text: &str, names: &Names) -> cas_lang::Statement {
    let parse = parse(text);
    assert!(parse.errors().is_empty(), "{text:?}: {:?}", parse.errors());
    let lowered = lower(&parse, names);
    assert!(lowered.diagnostics.is_empty(), "{text:?}: {:?}", lowered.diagnostics);
    let mut statements = lowered.statements;
    assert_eq!(statements.len(), 1, "{text:?}");
    statements.pop().expect("one statement")
}

fuzz_target!(|data: &[u8]| {
    let text = String::from_utf8_lossy(data);
    let names = Names::english();
    for s in lower(&parse(&text), names).statements {
        // Statements, patterns and rules have no input form.
        let Ok(printed) = print(&s.expr, names) else {
            continue;
        };
        let back = read(&printed, names).expr;
        assert_eq!(fold_arithmetic(&back), fold_arithmetic(&s.expr), "{printed:?}");
        let again = print(&back, names).expect("what the input form reads back as prints");
        let back_again = read(&again, names).expr;
        assert_eq!(print(&back_again, names).as_ref(), Ok(&again), "{printed:?}");
    }
});
