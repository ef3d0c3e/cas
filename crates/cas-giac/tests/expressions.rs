//! Xcas's expression commands against giac (CAS-562 to CAS-567): `pow`,
//! `comDenom`, `getNum`, `getDenom`, `simp2`, `regroup`, `sq`, `logb`,
//! `log2`, `log10`, `alog10`, `normalize`, `unitV` and `algsubs`, giac's answers read
//! back and compared with ours (`Giac::agrees`); and the departures, each
//! with giac's answer as it prints it: `powexpand`, whose answer ours
//! holds, as evaluation would merge the powers again, and `radians` and
//! `degrees`, exact where giac's are decimals.

use cas_core::{Context, Expr, Symbol};
use cas_giac::{Giac, ours};
use cas_lang::{Names, print};

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

/// `pow` (CAS-562); `comDenom`, `getNum`, `getDenom` and `simp2`
/// (CAS-563).
#[test]
fn powers_and_fractions() {
    check(&[
        "pow(x,3)",
        "pow(2,3)",
        "pow(x+1,452)",
        "powexpand(exp(x+y))",
        "comDenom(1/x+1/y^2+1)",
        "comDenom(1/x+1/y^2+1,x)",
        "comDenom((x^2-1)/(x-1)+1/x)",
        "comDenom(1/(x+1)+1/(x-1))",
        "getNum((x^3-1)/(x^2-1))",
        "getDenom((x^3-1)/(x^2-1))",
        "getNum(1+(x^3-1)/x^2)",
        "getDenom(1+(x^3-1)/x^2)",
        "getNum(25/15)",
        "getDenom(25/15)",
        "getNum(2*x/(3*y))",
        "getDenom(2*x/(3*y))",
        "getNum(-x/y)",
        "getDenom(-x/y)",
        "getNum([1/x,2/3])",
        "getNum(x^(-1/2)*y)",
        "getDenom(x^(-1/2)*y)",
        "simp2(12,18)",
        "simp2(4,6)",
        "simp2(2/3,4/9)",
        "simp2(x^3-1,x^2-1)",
        "simp2(x^2-1,x+1)",
    ]);
}

/// `regroup` (CAS-566), `sq` (CAS-567), the logarithms (CAS-565) and
/// `normalize` (CAS-564).
#[test]
fn regroup_squares_logarithms_and_norms() {
    check(&[
        "regroup(x+3*x+5*4/x)",
        "regroup((x^2-1)/(x-1))",
        "regroup(sin(x)^2+cos(x)^2)",
        "sq(x)",
        "sq(5)",
        "sq(1,2,3)",
        "sq([1,2,3])",
        "sq(x,y)",
        "logb(x,b)",
        "logb(x)",
        "logb(24,2)",
        "logb(5,2)",
        "logb(7,10)",
        "log2(8)",
        "log2(x)",
        "log10(10)",
        "log10(1000)",
        "log10(2000)",
        "log10(x)",
        "alog10(3)",
        "alog10(x)",
        "alog10([1,2])",
        "normalize(3+4*i)",
        "normalize([3,4])",
        "normalize(x)",
        "unitV([3,4])",
        "unitV(3+4*i)",
    ]);
}

/// `algsubs` (CAS-566): giac's elimination of an algebraic relation.
#[test]
fn algebraic_substitution() {
    check(&[
        "algsubs(x^2=u,1+x^2+x^4)",
        "algsubs(a*b/c=d,2*a*b^2/c)",
        "algsubs(2*a=p^2-q^2,algsubs(2*c=p^2+q^2,c^2-a^2))",
        "algsubs(x*y=u,x^2*y^2+x)",
        "algsubs(x+y=u,(x+y)^2+x+y)",
    ]);
}

/// `powexpand` (CAS-562): ours is giac's answer held, so it does not read
/// back equal. giac's answer as it prints it; ours, once released, its
/// value; and where the canonical order is giac's, ours prints as
/// `hold(` giac's text `)`.
#[test]
fn power_expand_departures() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let cases = [
        ("powexpand(2^(x+y))", "2^x*2^y", true),
        ("powexpand(3^(2*x))", "(3^x)^2", true),
        ("powexpand(2^(2*x))", "(2^x)^2", true),
        ("powexpand(a^(3*x*y))", "(a^(x*y))^3", true),
        // The canonical order puts `x^1` first, and sums and products
        // print otherwise.
        ("powexpand(x^(2*y+1))", "(x^y)^2*x^1", false),
        ("powexpand(2^(-x-y))", "1/2^x/2^y", false),
        ("powexpand(2^(x+y)+1)", "2^x*2^y+1", false),
        ("powexpand(2^(x+y)=3)", "2^x*2^y=3", false),
    ];
    let mut wrong = Vec::new();
    for (input, answer, same_text) in cases {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        let held = ours(input).unwrap();
        if held.head() != Some(&Expr::from(Symbol::HOLD)) {
            wrong.push(format!("{input}: ours is not held"));
            continue;
        }
        let released = ctx.eval(&held.args()[0]).unwrap();
        if released != giac.value(input).unwrap() {
            wrong.push(format!("{input}: ours released is not giac's value"));
        }
        let text = print(&held, Names::english()).unwrap();
        if same_text && text != format!("hold({answer})") {
            wrong.push(format!("{input}: ours prints {text}"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

/// Whether `a` and `b` at `x = 1` are within `10⁻⁹ (1 + |b|)`, giac's
/// decimals having 12 digits (C0006).
fn close(ctx: &Context, a: &Expr, b: &Expr) -> bool {
    let point = [(Expr::symbol("x"), Expr::from(1))];
    let (a, b) = (a.substitute(&point), b.substitute(&point));
    let gap = Expr::call(Symbol::ABS, [a - b.clone()]);
    let bound = (Expr::from(1) + Expr::call(Symbol::ABS, [b])) / Expr::from(1_000_000_000_i64);
    let less = Expr::call(Symbol::LESS, [Expr::call(Symbol::N, [gap]), bound]);
    ctx.eval(&less).ok().and_then(|e| e.as_symbol()) == Some(Symbol::TRUE)
}

/// `radians` and `degrees` (CAS-564): exact where giac's are decimals.
/// giac's answer as it prints it, ours as decided, and the two within
/// giac's digits.
#[test]
fn angle_departures() {
    let Some(giac) = giac() else {
        return;
    };
    let ctx = Context::new();
    let cases = [
        ("radians(180)", "3.14159265359", "pi"),
        ("radians(pi)", "0.0548311355616", "pi^2/180"),
        ("radians(x)", "0.0174532925199*x", "pi*x/180"),
        ("degrees(pi)", "180.0", "180"),
        ("degrees(1)", "57.2957795131", "180/pi"),
        ("degrees(x)", "57.2957795131*x", "180*x/pi"),
    ];
    let mut wrong = Vec::new();
    for (input, answer, decided) in cases {
        let printed = giac.eval(input).unwrap();
        if printed != answer {
            wrong.push(format!("{input}: giac prints {printed}, recorded {answer}"));
        }
        let value = ours(input).unwrap();
        if value != ours(decided).unwrap() {
            wrong.push(format!("{input}: ours is not {decided}"));
        }
        if !close(&ctx, &value, &giac.value(input).unwrap()) {
            wrong.push(format!("{input}: ours is not within giac's digits"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}
