//! The full-form reader of cas-core on arbitrary text: it never panics, and
//! whatever it reads prints in full form as text that reads back as the same
//! expression.
#![no_main]

use cas_core::Expr;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let text = String::from_utf8_lossy(data);
    match Expr::parse_full_form(&text) {
        Ok(e) => {
            let printed = format!("{e:?}");
            let back = Expr::parse_full_form(&printed).expect("printed full form reads");
            assert!(back == e, "{text:?} printed as {printed:?}");
        }
        Err(e) => {
            assert!(e.position <= text.len() && text.is_char_boundary(e.position));
            let _ = e.to_string();
        }
    }
});
