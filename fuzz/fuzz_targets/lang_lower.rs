//! Lowering on arbitrary text: total, with every range it reports within the
//! text, and every statement printable.
#![no_main]

use cas_lang::{Names, lower, parse};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let text = String::from_utf8_lossy(data);
    let lowered = lower(&parse(&text), Names::english());
    for d in &lowered.diagnostics {
        assert!(d.range.start <= d.range.end && d.range.end <= text.len());
        if let Some(fix) = &d.fix {
            assert!(fix.range.start <= fix.range.end && fix.range.end <= text.len());
        }
        let _ = d.to_string();
    }
    for s in &lowered.statements {
        let range = s.spans.range();
        assert!(range.start <= range.end && range.end <= text.len());
        let _ = format!("{:?}", s.expr);
    }
});
