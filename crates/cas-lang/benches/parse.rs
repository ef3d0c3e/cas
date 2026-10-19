//! Lexing and parsing time against the size of the text (CAS-118): both
//! should be linear, so the time per byte stays flat from 10 kB to 1 MB.
//!
//! The text repeats a program that uses every construct: definitions with
//! guards, lambdas, blocks and loops, lists, implicit products, powers and
//! comparisons. `lex/bytes` times the lexer alone, `parse/bytes` the lexer
//! and the parser. Pass a filter to run part of it:
//! `cargo bench -p cas-lang --bench parse -- '^parse/'`.
#![allow(missing_docs, reason = "benchmark binary")]

use std::hint::black_box;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

const PROGRAM: &str = r#"// A sample of the language.
f(n) := n*f(n - 1) when n > 0
f(0) := 1
g := (x, y) -> x^2 + 2x y + y²
roots := [(-b + √(b^2 - 4a c))/(2a), (-b - √(b^2 - 4a c))/(2a)]
collatz(n) := begin
  m := n
  steps := 0
  while m != 1 do
    m := if even(m) then m/2 else 3m + 1 end
    steps := steps + 1
  end
  steps
end
for k in 0..10 do s := s + k^2 end
p := 3x^2(x - 1) + sin(x)cos(x) - 1/2x
ok := 0 < x <= 1 and not (y == 2)
msg := "sum = \u{3a3}"
"#;

/// `PROGRAM` repeated to at least `size` bytes.
fn text(size: usize) -> String {
    PROGRAM.repeat(size.div_ceil(PROGRAM.len()))
}

fn bench(c: &mut Criterion) {
    // The benchmark measures parsing, not error recovery.
    let errors = cas_lang::parse(PROGRAM).errors().to_vec();
    assert!(errors.is_empty(), "{errors:?}");
    for (group, run) in [
        (
            "lex",
            (|t: &str| cas_lang::lex(t).len()) as fn(&str) -> usize,
        ),
        ("parse", |t: &str| cas_lang::parse(t).errors().len()),
    ] {
        let mut group = c.benchmark_group(group);
        for size in [10_000, 100_000, 1_000_000] {
            let text = text(size);
            group.throughput(Throughput::Bytes(text.len() as u64));
            group.bench_with_input(BenchmarkId::new("bytes", size), &text, |b, text| {
                b.iter(|| run(black_box(text)));
            });
        }
        group.finish();
    }
}

criterion_group!(benches, bench);
criterion_main!(benches);
