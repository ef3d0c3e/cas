//! Xcas's integer functions against giac (CAS-548 to CAS-553): `chrem`,
//! `lgcd`, `fracmod` and `iratrecon`, `icomp` and `bernoulli_mod`, giac's
//! answers read back and compared with ours (`Giac::agrees`: the same
//! expression, or `==`); and the departures, each with giac's answer as it
//! prints it: its errors, where ours leaves the call with a warning, its
//! Horner form, its sign of a polynomial gcd, and the inputs that crash
//! libgiac, checked on our side only.

use cas_giac::{Giac, ours};

/// The shared library, or `None` with a note: the comparisons skip.
fn giac() -> Option<&'static Giac> {
    let giac = Giac::shared();
    if giac.is_none() {
        eprintln!("libgiac is missing: comparisons skipped");
    }
    giac
}

fn check(inputs: &[&str]) {
    let Some(giac) = giac() else {
        return;
    };
    let wrong: Vec<String> = inputs
        .iter()
        .filter_map(|input| giac.agrees(input).err().map(|e| format!("{e}")))
        .collect();
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// `chrem` of integers, and of polynomials coefficient by coefficient
/// modulo the product of the moduli (CAS-548).
#[test]
fn chinese_remainders() {
    check(&[
        "chrem([2,3],[5,7])",
        "chrem([2,3],[7,5])",
        "chrem([2,4,6],[3,5,7])",
        "chrem([2,4,6,7],[3,5,7,11])",
        "chrem([2,3,1],[5,7,11])",
        "chrem([1,3],[2,4])",
        "chrem([-2,3],[5,7])",
        "chrem([-1000,3],[5,7])",
        "chrem([3,-1000],[5,7])",
        "chrem([2,3],[2^40,7])",
        "chrem([[1,2],[3,4]])",
        "chrem([2],[5])",
        "chrem([-12],[5])",
        "chrem([x],[5])",
        "chrem([2*x+1,4*x+2,6*x-1,x+1],[3,5,7,11])",
        "chrem([x,1],[3,5])",
        "chrem([3*x+1,1],[3,5])",
        "chrem([7*x,8],[3,5])",
        "chrem([x,0],[3,5])",
        "chrem([x,x],[3,5])",
        "chrem([x*y,1],[3,5])",
        "chrem([x+1,2,3],[3,5,7])",
        "chrem([x+1,x,3],[3,5,7])",
        "chrem([x+1,x+4],[3,6])",
        "chrem([x+1,x+3],[4,6])",
        "chrem([x+5,x+11],[6,15])",
        "chrem([x+7,x+1],[10,4])",
        "chrem([x-1,x+5],[4,10])",
        "chrem([x+1,x+2],[2^40,3])",
        "chrem([x+2^39,x],[2^40,3])",
        "chrem([x-2^39,x],[2^40,3])",
        "chrem([x+2^35,x+2],[2^40,6])",
        "chrem([2147483647*x,1],[3,5])",
        "chrem([x,1],[2147483647,5])",
        "chrem([x+1,2],[2147483647,2147483629])",
    ]);
}

/// `lgcd` of lists and sequences (CAS-549).
#[test]
fn list_gcds() {
    check(&[
        "lgcd([])",
        "lgcd()",
        "lgcd([4])",
        "lgcd([-4])",
        "lgcd([0,0])",
        "lgcd([-4,6])",
        "lgcd(-4,-6)",
        "lgcd(18,12,30)",
        "lgcd([18,12,30])",
        "lgcd([45,75,20,15])",
        "lgcd([1/2,1/3])",
        "lgcd([2*x+2,4*x+4])",
        "lgcd([-2*x-2,4*x+4])",
        "lgcd([6,x])",
        "lgcd([x])",
        "lgcd(x,y)",
        "lgcd([x^2-1,x^2-2*x+1])",
        "lgcd([x^2-2*x+1,x^3-1,x-1])",
        "lgcd(x^2-2*x+1,x^3-1,x-1)",
        "lgcd(6)",
    ]);
}

/// `fracmod` and its alias `iratrecon` (CAS-550).
#[test]
fn fractions_modulo() {
    check(&[
        "fracmod(13,121)",
        "fracmod(41,121)",
        "fracmod(3,7)",
        "fracmod(10,7)",
        "fracmod(-3,7)",
        "fracmod(0,7)",
        "fracmod(5,1)",
        "fracmod([13,41],121)",
        "fracmod([13,[41]],121)",
        "fracmod(123456789,10^20+39)",
        "iratrecon(41,121)",
        "iratrecon(13,121)",
        "fracmod(13)",
        "fracmod(13,121,3)",
    ]);
}

/// `icomp` in NEXCOM's order, with `zeros` (CAS-553).
#[test]
fn compositions() {
    check(&[
        "icomp(4,1)",
        "icomp(3,3)",
        "icomp(4,2)",
        "icomp(4,2,zeros=true)",
        "icomp(5,2,zeros=false)",
        "icomp(3,3,zeros=false)",
        "icomp(6,3,zeros=false)",
        "icomp(3,2,zeros=1)",
    ]);
}

/// `bernoulli_mod`, giac's signed values (CAS-552).
#[test]
fn bernoulli_numbers_modulo() {
    check(&[
        "bernoulli_mod(10,13)",
        "bernoulli_mod(14,13)",
        "bernoulli_mod(1234,1009)",
        "bernoulli_mod(20,10007)",
        "bernoulli_mod(2,7)",
        "bernoulli_mod(8,17)",
        "bernoulli_mod(2,5)",
        "bernoulli_mod(1234,10^6+3)",
        "bernoulli_mod(50,1000003)",
    ]);
}

/// Where this CAS departs from giac: the input, giac's answer as it prints
/// it, ours as decided.
const DEPARTURES: [(&str, &str, &str); 24] = [
    // giac's errors; ours leaves the call, with the error as a warning.
    (
        "chrem([2,3],[4,6])",
        "GIAC_ERROR: No Integer Solution Error: Bad Argument Value",
        "chrem([2,3],[4,6])",
    ),
    (
        "chrem([x^2+1,2*x],[3,6])",
        "GIAC_ERROR: No Integer Solution Error: Bad Argument Value",
        "chrem([x^2+1,2*x],[3,6])",
    ),
    (
        "chrem([2,3,4],[5,7])",
        "GIAC_ERROR: chrem([2,3,4],[5,7]) \n Error: Bad Argument Value",
        "chrem([2,3,4],[5,7])",
    ),
    (
        "chrem([],[])",
        "GIAC_ERROR: chrem([],[]) \n Error: Bad Argument Value",
        "chrem([],[])",
    ),
    (
        "chrem(2,5)",
        "GIAC_ERROR: chrem(2,5) \n Error: Bad Argument Value",
        "chrem(2,5)",
    ),
    (
        "chrem([1/2,1],[3,5])",
        "GIAC_ERROR: Non integer argument Error: Bad Argument Type",
        "chrem([1/2,1],[3,5])",
    ),
    (
        "chrem([1,2],[1,5])",
        "GIAC_ERROR: Bad mod value Error: Bad Argument Type",
        "chrem([1,2],[1,5])",
    ),
    (
        "fracmod(2,4)",
        "GIAC_ERROR: Reconstructed denominator is not prime with modulo Error: Bad Argument Value",
        "fracmod(2,4)",
    ),
    (
        "fracmod(5,8)",
        "GIAC_ERROR: Reconstructed denominator is not prime with modulo Error: Bad Argument Value",
        "fracmod(5,8)",
    ),
    (
        "fracmod(13,x)",
        "GIAC_ERROR: Reconstructed denominator is not prime with modulo Error: Bad Argument Value",
        "fracmod(13,x)",
    ),
    (
        "fracmod(1/2,7)",
        "GIAC_ERROR: Reconstructed denominator is not prime with modulo Error: Bad Argument Value",
        "fracmod(1/2,7)",
    ),
    (
        "fracmod(1,0)",
        "GIAC_ERROR: Reconstructed denominator is not prime with modulo Error: Bad Argument Value",
        "fracmod(1,0)",
    ),
    (
        "icomp(3,4)",
        "GIAC_ERROR: icomp(3,4) \n Error: Bad Argument Value",
        "icomp(3,4)",
    ),
    (
        "icomp(0,1)",
        "GIAC_ERROR: icomp(0,1) \n Error: Bad Argument Value",
        "icomp(0,1)",
    ),
    (
        "icomp(4,0)",
        "GIAC_ERROR: icomp(4,0) \n Error: Bad Argument Value",
        "icomp(4,0)",
    ),
    (
        "icomp(3,2,zeros=x)",
        "GIAC_ERROR: icomp(3,2,'zeros'=x) \n Error: Bad Argument Value",
        "icomp(3,2,zeros=x)",
    ),
    (
        "icomp(5)",
        "GIAC_ERROR: icomp(5) \n Error: Bad Argument Type",
        "icomp(5)",
    ),
    (
        "bernoulli_mod(4,5)",
        "GIAC_ERROR: Not invertible Error: Bad Argument Value",
        "bernoulli_mod(4,5)",
    ),
    (
        "bernoulli_mod(2,3)",
        "GIAC_ERROR: Not invertible Error: Bad Argument Value",
        "bernoulli_mod(2,3)",
    ),
    (
        "bernoulli_mod(3,13)",
        "GIAC_ERROR: bernoulli_mod(3,13) \n Error: Bad Argument Type",
        "bernoulli_mod(3,13)",
    ),
    (
        "bernoulli_mod(10,12)",
        "GIAC_ERROR: bernoulli_mod(10,12) \n Error: Bad Argument Type",
        "bernoulli_mod(10,12)",
    ),
    (
        "bernoulli_mod(0,7)",
        "GIAC_ERROR: bernoulli_mod(0,7) \n Error: Bad Argument Type",
        "bernoulli_mod(0,7)",
    ),
    (
        "bernoulli_mod(4,9)",
        "GIAC_ERROR: bernoulli_mod(4,9) \n Error: Bad Argument Type",
        "bernoulli_mod(4,9)",
    ),
    // GCD makes the leading coefficient positive; giac keeps the sign
    // of a single polynomial.
    ("lgcd([-2*x-2])", "-2*x-2", "2*x+2"),
];

/// The departures of [`DEPARTURES`], then giac's Horner form and the
/// inputs that crash libgiac.
#[test]
fn departures() {
    let Some(giac) = giac() else {
        return;
    };
    let mut wrong = Vec::new();
    for (input, answer, decided) in DEPARTURES {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        if ours(input) != ours(decided) {
            wrong.push(format!("{input}: ours is not {decided}"));
        }
        if giac.agrees(input).is_ok() {
            wrong.push(format!("{input} agrees now"));
        }
    }
    // giac prints the polynomial in Horner form; expanded, it is ours.
    let (input, answer) = ("chrem([x^2+1,2*x],[3,5])", "[(-5*x+12)*x-5,15]");
    let printed = giac.eval(input).unwrap();
    if printed != answer {
        wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
    }
    if ours(input) != ours("[-5*x^2+12*x-5,15]")
        || ours(input) != ours(&format!("expand({answer})"))
    {
        wrong.push(format!("{input}: ours is not giac's expanded"));
    }
    // libgiac crashes on polynomial moduli and on p = 2, so giac is not
    // asked: ours folds through chinrem, and leaves bernoulli_mod(2,2).
    for (input, decided) in [
        ("chrem([x+1,x],[x^2+1,x+2])", "chinrem([x+1,x^2+1],[x,x+2])"),
        ("chrem([x+1,x],[x^2+1,x-1])", "chinrem([x+1,x^2+1],[x,x-1])"),
        ("bernoulli_mod(2,2)", "bernoulli_mod(2,2)"),
    ] {
        if ours(input) != ours(decided) {
            wrong.push(format!("{input}: ours is not {decided}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
