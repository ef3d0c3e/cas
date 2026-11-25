"""Emit crates/cas-core/src/native/wavelet/filters.rs from GSL's wavelet sources (argv[1]) to stdout. Origin: phase 1.7, T9.c1 scratch."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gslfilters import arrays, BSPLINE

d = sys.argv[1]
a = {}
for f in ("daubechies.c", "haar.c", "bspline.c"):
    a.update(arrays(d + "/" + f))

used = ["ch_2", "cg_2"] + [f"{p}_{k}" for k in range(4, 21, 2) for p in ("h", "g")]
for k, (nc, h1, (g1, _), (h2, _), g2) in BSPLINE.items():
    for n in (h1, g1, h2, g2):
        if n not in used:
            used.append(n)

out = []
w = out.append
w("//! GSL 2.7.1's wavelet filters (wavelet/haar.c, daubechies.c, bspline.c;")
w("//! GPL-3), their decimal literals as written, `M_SQRT1_2` as gsl_math.h's")
w("//! `0.70710678118654752440084436210`. Transcribed by a script from the C")
w("//! arrays once (crates/cas-core/proto/wavelet/proto/gen_filters.py, T9.c1); each `*_init` function's choice of arrays and offsets is")
w("//! [`filters`].")
w("")
w("use super::Kind;")
w("")
w("/// A filter: decimal literals.")
w("type Taps = &'static [&'static str];")
w("")
w("/// One wavelet's filters, `nc` decimal literals each: `h1` and `g1` for the")
w("/// forward transform, `h2` and `g2` for the inverse (GSL's `gsl_wavelet`).")
w("pub(super) struct Filters {")
w("    pub(super) h1: Taps,")
w("    pub(super) g1: Taps,")
w("    pub(super) h2: Taps,")
w("    pub(super) g2: Taps,")
w("}")
w("")
for n in used:
    vals = a[n]
    w(f"const {n.upper()}: [&str; {len(vals)}] = [")
    for v in vals:
        w(f'    "{v}",')
    w("];")
    w("")
w("/// The filters of `kind` with member `k` (GSL's `haar_init`,")
w("/// `daubechies_init` and `bspline_init`, the centred variants alike), or")
w("/// `None` for a member GSL does not implement.")
w("pub(super) fn filters(kind: Kind, k: i64) -> Option<Filters> {")
w("    let f = |h1: Taps, g1: Taps, h2: Taps, g2: Taps| Filters { h1, g1, h2, g2 };")
w("    Some(match (kind, k) {")
w("        (Kind::Haar, 2) => f(&CH_2, &CG_2, &CH_2, &CG_2),")
for k in range(4, 21, 2):
    w(f"        (Kind::Daubechies, {k}) => f(&H_{k}, &G_{k}, &H_{k}, &G_{k}),")
for k, (nc, h1, (g1, o1), (h2, o2), g2) in BSPLINE.items():
    w(f"        (Kind::Bspline, {k}) => f(&{h1.upper()}, &{g1.upper()}[{o1}..{o1 + nc}], &{h2.upper()}[{o2}..{o2 + nc}], &{g2.upper()}),")
w("        _ => return None,")
w("    })")
w("}")
print("\n".join(out))
