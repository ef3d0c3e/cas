//! Permutations, bits and transforms against giac (CAS-DOC-45,
//! CAS-PLAN-39): giac's permutation commands on its examples and on
//! pseudo-random permutations, its bitwise operations, `hamdist`,
//! `convert` and `perm`, and `fft`, `ifft` and `convolution`, read back:
//! its decimals equal to this CAS's balls or exact values where both are
//! integers or Gaussian integers, and its modular transforms exactly.
//! `is_permu`'s 1 and 0 are this CAS's `True` and `False`.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the checks skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: skipped");
    }
    giac
}

/// giac's answers to `inputs`, read back, are this CAS's.
fn check(inputs: &[String]) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| giac.agrees(input).err().map(|e| e.to_string()))
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

fn owned(inputs: &[&str]) -> Vec<String> {
    inputs.iter().map(|s| (*s).to_owned()).collect()
}

/// Pseudo-random permutations of sizes 1 to 8, by Fisher and Yates on a
/// linear congruential sequence.
fn permutations() -> Vec<Vec<usize>> {
    let mut state: u64 = 0x2545_f491_4f6c_dd1d;
    let mut next = |bound: usize| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        usize::try_from((state >> 33) % u64::try_from(bound).expect("small")).expect("small")
    };
    (1..=8)
        .flat_map(|n| std::iter::repeat_n(n, 3))
        .map(|n| {
            let mut p: Vec<usize> = (0..n).collect();
            for i in (1..n).rev() {
                let j = next(i + 1);
                p.swap(i, j);
            }
            p
        })
        .collect()
}

fn written(p: &[usize]) -> String {
    let entries: Vec<String> = p.iter().map(ToString::to_string).collect();
    format!("[{}]", entries.join(","))
}

/// giac's examples of the permutation commands (CAS-463).
#[test]
fn permutation_examples() {
    check(&owned(&[
        "permu2cycles([1,3,4,5,2,0])",
        "permu2cycles([0,2,1,3])",
        "permu2cycles([1,0,2])",
        "permu2cycles([4,3,2,1,0])",
        "cycles2permu([[1,3,5],[2,4]])",
        "cycles2permu([[2],[0,1]])",
        "cycle2perm([1,3,5])",
        "cycle2perm([3,1])",
        "perminv([1,3,4,5,2,0])",
        "cycleinv([1,3,5])",
        "p1op2([1,2,0],[2,1,0])",
        "p1oc2([1,2,0],[0,2])",
        "c1op2([0,2],[1,2,0])",
        "c1oc2([0,1],[1,2])",
        "permuorder([1,3,4,5,2,0])",
        "signature([1,3,4,5,2,0])",
        "permu2mat([1,2,0])",
    ]));
    let Some(giac) = giac() else {
        return;
    };
    for (input, expected) in [("is_permu([2,0,1])", true), ("is_permu([1,1,0])", false)] {
        let theirs = giac.value(input).unwrap();
        let ours = ours(input).unwrap();
        assert_eq!(theirs, Expr::from(i64::from(expected)), "{input}: giac's");
        let truth = if expected {
            Symbol::TRUE
        } else {
            Symbol::FALSE
        };
        assert_eq!(ours, Expr::from(truth), "{input}: ours");
    }
}

/// Pseudo-random permutations, and products of pairs (CAS-463).
#[test]
fn random_permutations() {
    let all = permutations();
    let mut inputs = Vec::new();
    for (k, p) in all.iter().enumerate() {
        let p = written(p);
        for command in ["permu2cycles", "perminv", "permuorder", "signature"] {
            inputs.push(format!("{command}({p})"));
        }
        if let Some(q) = all.get(k + 1).filter(|q| q.len() == all[k].len()) {
            inputs.push(format!("p1op2({p},{})", written(q)));
        }
    }
    check(&inputs);
}

/// Bits, digits and falling factorials (CAS-464).
#[test]
fn bits_and_digits() {
    check(&owned(&[
        "bitand(12,10)",
        "bitor(12,10)",
        "bitxor(12,10)",
        "bitand(12,10,6)",
        "bitand(-12,10)",
        "bitor(-12,10)",
        "bitxor(-12,10)",
        "bitxor(255,-256)",
        "hamdist(5,3)",
        "hamdist(15,0)",
        "convert(255,base,16)",
        "convert(-255,base,16)",
        "convert(0,base,2)",
        "convert(1000,base,7)",
        "convert([6,2,6,2],base,7)",
        "convert([15,15],base,16)",
        "perm(10,3)",
        "perm(7,7)",
        "perm(3,5)",
    ]));
}

/// Transforms and convolutions (CAS-465): Gaussian integers, so giac's
/// decimals are exact.
#[test]
fn transforms() {
    check(&owned(&[
        "fft([1,2,3,4])",
        "fft([1,0,-1,2])",
        "fft([3,-1])",
        "ifft([4,2,0,2])",
        "ifft([1,2])",
        "fft([1,2,3,4],2,5)",
        "ifft([1,2,3,4],2,5)",
        "fft([1,2,3,4,5,6,7,8],9,17)",
        "ifft([1,2,3,4,5,6,7,8],9,17)",
        "fft([5,1,4],7,19)",
    ]));
    // giac's convolutions are decimals, by its fft; ours exact.
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    for input in ["convolution([1,2,3],[1,-1])", "convolution([1,2],[3,4,5])"] {
        let (theirs, ours) = (giac.value(input).unwrap(), ours(input).unwrap());
        let gap = Expr::call(
            Symbol::MAX,
            [Expr::call(
                Symbol::ABS,
                [Expr::call(Symbol::N, [ours.clone() - theirs.clone()])],
            )],
        );
        let small = Expr::call(Symbol::LESS, [gap, Expr::from(10).pow(-9)]);
        let holds = ctx.eval(&small).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE);
        assert!(holds, "{input}: ours {ours:?}, giac's {theirs:?}");
    }
}
