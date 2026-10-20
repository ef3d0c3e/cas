//! The session on arbitrary text (CAS-121): running it never panics, every
//! output and message lies within the text on character boundaries, and
//! every result it shows in input form parses and lowers without errors.
//!
//! Each input runs with small limits, and a watchdog interrupts it after a
//! second and then every 10 ms: a program can branch without end, which no
//! limit stops (CAS-DOC-4 R8), and the session goes on with the next
//! statement after an interrupted one.
#![no_main]

use std::ops::Range;
use std::sync::mpsc;
use std::time::Duration;

use cas_core::Limits;
use cas_lang::{Session, Shown, lower, parse};
use libfuzzer_sys::fuzz_target;

fn within(text: &str, range: &Range<usize>) -> bool {
    range.start <= range.end
        && range.end <= text.len()
        && text.is_char_boundary(range.start)
        && text.is_char_boundary(range.end)
}

fuzz_target!(|data: &[u8]| {
    let text = String::from_utf8_lossy(data);
    let mut session = Session::new();
    let mut limits = Limits::default();
    limits.recursion = 200;
    limits.iterations = 256;
    limits.match_steps = 20_000;
    limits.number_bits = 1 << 12;
    session.set_limits(limits);
    let interrupt = session.interrupt();
    let (done, finished) = mpsc::channel::<()>();
    let watchdog = std::thread::spawn(move || {
        let mut wait = Duration::from_secs(1);
        while finished.recv_timeout(wait) == Err(mpsc::RecvTimeoutError::Timeout) {
            interrupt.interrupt();
            wait = Duration::from_millis(10);
        }
    });
    let outputs = session.run(&text);
    for o in &outputs {
        assert!(within(&text, &o.range), "{:?}", o.range);
        for m in &o.messages {
            assert!(within(&text, &m.range), "{m:?}");
            if let Some(fix) = &m.fix {
                assert!(within(&text, &fix.range), "{m:?}");
            }
        }
        if let Some(Shown::Input(result)) = &o.result {
            let parsed = parse(result);
            assert!(
                parsed.errors().is_empty(),
                "{result:?}: {:?}",
                parsed.errors()
            );
            let lowered = lower(&parsed, session.names());
            assert!(
                !lowered.has_errors(),
                "{result:?}: {:?}",
                lowered.diagnostics
            );
        }
    }
    drop(done);
    watchdog.join().unwrap();
});
