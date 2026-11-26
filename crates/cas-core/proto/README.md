# Statistics and signal prototypes (Python ports of giac, phase 1.7)

Our answers are seen only on CI (R74), so the expected rows of phase 1.7's natives come from these
ports, each checked row for row against giac 2.0.0.21 through libgiac (`scripts/caseval.py`, C0018).
Not compiled, not shipped. Moved from phase scratch in 1.7/T12; the copies differ from the originals
only in their paths and a line naming the task they came from.

- `giac/` (the planner's probes): `probe1.cas`..`probe15.cas`, giac's answers read by the prototypes
  below and cited by `crates/cas-giac/tests/*.rs`: probe1, probe10 markov; probe2, probe14 randvar;
  probe3 kernel_density, fitdistr; probe4, probe8 cluster, kmeans; probe5, probe9, probe10, probe11
  the sound clips.
- `randvar/` (T2, CAS-624): `proto.py` ports randvar's parameter fits and discrete forms
  (graphtheory.cc:612-866) in exact arithmetic, prints `rows off: 0`; `probe.cas` the fits, messages
  and multinomial forms; `c2/rows.cas`, `c2/extra.cas` the T2.c2 oracle rows. Cited by
  `crates/cas-giac/tests/random_variables.rs`.
- `markov/` (T4.c1, CAS-628): `markov.py` ports markov (misc.cc:7082-7400) with exact solves and
  prints each row against the `.out` files given; `probe_a.cas`..`probe_f.cas`. Cited by
  `src/native/statistics/markov/tests.rs`, `crates/cas-giac/tests/markov.rs`.
- `density/` (T5.c1, CAS-629): `proto.py` ports kernel_density and fitdistr
  (optimization.cc:7160-7546, misc.cc:3108-3171), reproduces `giac/probe3.out` and `extra.out`,
  writes `expected.txt`; `extra.cas`. Cited by `src/native/statistics/{density,fitting}/tests.rs`,
  `crates/cas-giac/tests/{fitdistr,kernel_density}.rs`, C0084.
- `clustering/` (T6.c1, CAS-630): `proto.py` ports kmeans (AS 136) and cluster (hclust) from
  optimization.cc against `giac/probe4.out`, `probe8.out`; 1.8/T5 `t5` (`probe_clustering2`) and
  `t5c3` (`probe_clustering3`: kmeans' `cluster_crit` indices, every row over 200 seeds). Cited by
  `src/native/statistics/clustering{.rs,/tests.rs}`, `crates/cas-giac/tests/clustering.rs`, R1.7-012.
- `sound/clip/` (T7.c1, CAS-631): `proto.py` the audio clip with giac's 47-bit gen doubles,
  `proto_f64.py` the same in plain f64; rows from `giac/probe5,9,10,11.out` and `dur.out`; each writes
  `table.txt` / `table_f64.txt`. Cited by `src/native/sound/tests.rs`,
  `crates/cas-giac/tests/sound.rs`, C0088.
- `sound/edit/` (T8, CAS-632, CAS-633): `proto.py` set_channel_data, stereo2mono, mixdown, splice
  (imports `../clip/`), writes `table.txt`; `resample.py` the resample filters, writes
  `resample.log`; `probe1.cas`, `probe2.cas`, `probe_rs.cas`. Cited by `src/native/sound/{resample,tests}.rs`,
  `crates/cas-giac/tests/sound.rs`, C0090.
- `wavelet/` (T9, CAS-634, CAS-635): `proto_dwt.py` ports GSL wavelet/dwt.c against `giac/p1.out`,
  `p2.out` (filters through `proto/gslfilters.py`); `proto/gen_filters.py` emitted
  `src/native/wavelet/filters.rs` (its header since reworded); `proto_emd.py` giac's emd
  with GSL's spline; `proto_analytic.py` dht, instfreq, instphase against `giac/p3.out`; `giac/p4.cas`
  decimal read-back probes. Cited by `src/native/wavelet/tests.rs`,
  `crates/cas-giac/tests/transforms.rs`.
- `leftovers/` (1.8/T4, CAS-640, CAS-641; moved from scratch in 1.8/T12): `lin_period.py` ports lin
  (lin.cc) and period (signalprocessing.cc); `proto_linabs.py` linabs' reduce_Heaviside tail;
  `proto_env.py` the envelope curves (giac answers a plot, so these are the expected rows);
  `probe1.cas`, `probe2.cas` giac's shift_phase, linabs, lin, period answers. Cited by
  `crates/cas-giac/tests/gap_leftovers_two.rs` (README ## Tests).

The wavelet scripts read GSL 2.7.1's sources (GPL-3, C0091): fetch
`https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz` (sha256
`dcb0fbd43048832b757ff9942691a8dd70026d5da0ff85601e52687f6deeb34b`) to `.run/<phase>/dl/`,
extract `wavelet/`, and pass `<dir>/gsl-2.7.1/wavelet` as the argument.

Rerun (repo root). Answers: `/usr/bin/python3 scripts/caseval.py <x>.cas > <x>.out` for every `.cas`
(`markov/probe_b.cas` with a 15 s limit: `caseval.py markov/probe_b.cas 15`, so its hang reads
`<timeout 15s>` as markov.py expects). Prototypes, by path, any working directory:
`/usr/bin/python3 crates/cas-core/proto/<dir>/proto.py` for randvar, density, clustering,
sound/clip (and `proto_f64.py`), sound/edit (and `resample.py`); `markov/markov.py
crates/cas-core/proto/giac/probe1.out ... crates/cas-core/proto/markov/probe_f.out`;
`/usr/bin/python3 -I wavelet/proto_dwt.py <gsl wavelet dir>`, `wavelet/proto_emd.py`,
`wavelet/proto_analytic.py`, `wavelet/proto/gen_filters.py <gsl wavelet dir>` (to stdout);
`leftovers/lin_period.py`, `leftovers/proto_linabs.py`, `leftovers/proto_env.py` (sympy, to stdout).
Generated, not committed (`.gitignore`): every `.out`, `table.txt`, `table_f64.txt`, `expected.txt`,
`resample.log`.
