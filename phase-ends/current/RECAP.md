## Recap
Phase 1.7 Statistics III and signal processing II: CAS-DOC-59 (still a draft in Lific) covers 37 giac entries; tickets CAS-624..635 are done, and CAS-636 (diff(e,x,n) doesn't take the nth derivative) is filed and stays todo. Built:
- Random variables: randvar/random_variable. A named law is its density call left unevaluated (normald(5,3/2)); the binomial is randvar(binomial,n,p); discrete laws are given as lists. sample, randvector and randmultinomial use the seeded generator.
- Probabilities and fits: UTPN/T/C/F (exact where possible, UTPC(2,6)=exp(-3)), cdf, icdf, mgf, logistic, logistic_regression. cdf of a data list gives count(≤x)/n, a departure from giac. kernel_density/kde and fitdistr: closed forms are exact; density lists compare within giac's FFT noise; the gamma fit is a departure.
- markov and randmarkov, exact on exact input. kmeans (Hartigan-Wong) and cluster (hierarchical).
- Sound clips: SoundClip expressions holding exact integer codes; createwav, the property commands, about, set_channel_data, stereo2mono, mixdown, splice and resample. Resample qualities 2-4 are our own windowed-sinc filters, stopband 102-147 dB.
- Signal transforms: dwt/idwt (GSL filters), hilbert of a list, instfreq, instphase, emd (bit for bit with giac).
Fixes outside the plan: the reader's f(x)→f*x rule skips function names; a simplify walk that could take 2^80 steps; the fuzz harness reseeds before each compared evaluation; nested LowPass no longer grows exponentially. Final gate green on 0e22fdf (CI 45637752766, Fuzz 45637752779). Eval peak memory 128 MB of 2048, so the symbol_array memory watch is closed. The cited prototypes and probes moved to crates/cas-core/proto/. Cookbook C0077-C0093.
Open, carried to 1.8 in the inbox:
- giac's Weibull fit drops a square, and ours copies it. I recommend the correct fit plus a departure row.
- kmeans/cluster output=part, index= and string clustering, and hilbert(f,t), go to the gap check.
- Coder briefs restate R74 (broken twice this phase, in T6 and T7; nothing from those runs was used).
Also open: publish CAS-DOC-57/58/59; commit_task.sh still sweeps logs/ (R73/R96); docs-only commits stay unpushed until the next phase's first push.
