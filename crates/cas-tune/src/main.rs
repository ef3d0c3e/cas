//! Measures cas-num's algorithm thresholds on this machine and writes a tuning
//! configuration (see `crates/cas-num/build.rs`).
//!
//! ```text
//! cargo run --release -p cas-tune -- [--quick] [--only key,key] [--base FILE] --out FILE
//! ```
//!
//! Every candidate value is measured as compiled: for each one the driver
//! writes a configuration, rebuilds `cas-tune-probe` against it (in
//! `target/cas-tune`, so other builds are untouched) and times the operation
//! the threshold controls over a window of sizes around the plausible
//! crossovers. Building with threshold `T` puts the new algorithm on top from
//! `T` up, over whatever the earlier thresholds chose below, so comparing the
//! candidates (including "never") compares the algorithms as they are really
//! combined. A candidate's cost is the mean, over the window, of its time
//! divided by the best candidate's time at that size: no single size
//! dominates, and a flat optimum shows up as near-equal costs. The best grid
//! point is then refined with the midpoints to its neighbours.
//!
//! Timings need an idle machine: the driver warns when the load average
//! shows other work, and reports the noise it measured (how much the two
//! passes over the same binary disagree), which bounds how far the costs
//! can be trusted.
//!
//! Thresholds are tuned in dependency order, each with the values already
//! chosen, starting from `--base` (default `tuning/default.toml`).
#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "sizes and counts far below 2^52, rounded on purpose"
)]

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;
use std::{env, fs};

/// A candidate value that disables the algorithm in every window.
const NEVER: u64 = 1_000_000_000;

/// What to measure for one threshold.
struct Param {
    key: &'static str,
    /// The probe operation.
    op: &'static str,
    /// Candidate values, increasing.
    grid: Vec<u64>,
    /// Sizes to time, from the values already chosen; by default the grid's
    /// range up to 2.5 times its largest finite value, so every candidate
    /// is timed both where it takes over and inside the recursion above.
    window: Option<fn(&Config) -> Vec<u64>>,
}

impl Param {
    /// The window, plus the powers of two inside it: transform costs step
    /// there, and they are common sizes, which a geometric sweep rarely hits.
    fn sizes(&self, config: &Config) -> Vec<u64> {
        let mut sizes = self.window.map_or_else(
            || {
                let top = self.grid.iter().filter(|&&v| v != NEVER).max();
                geometric(self.grid[0], top.expect("a finite candidate") * 5 / 2, 1.25)
            },
            |window| window(config),
        );
        if let (Some(&lo), Some(&hi)) = (sizes.first(), sizes.last()) {
            sizes.extend((0..64).map(|k| 1u64 << k).filter(|&p| lo < p && p < hi));
            sizes.sort_unstable();
            sizes.dedup();
        }
        sizes
    }
}

/// `start`, `start * ratio`, ... up to `end`, rounded and deduplicated.
fn geometric(start: u64, end: u64, ratio: f64) -> Vec<u64> {
    let mut v = Vec::new();
    let mut x = start as f64;
    while x <= end as f64 * 1.000_001 {
        let r = x.round() as u64;
        if v.last() != Some(&r) {
            v.push(r);
        }
        x *= ratio;
    }
    v
}

fn with_never(mut grid: Vec<u64>) -> Vec<u64> {
    grid.push(NEVER);
    grid
}

#[expect(
    clippy::too_many_lines,
    reason = "a table: one entry per key, in tuning order"
)]
fn params() -> Vec<Param> {
    vec![
        Param {
            key: "mul_karatsuba_threshold",
            op: "mul",
            grid: geometric(4, 64, 1.19),
            window: None,
        },
        Param {
            key: "mul_toom3_threshold",
            op: "mul",
            grid: with_never(geometric(48, 640, 1.19)),
            window: None,
        },
        Param {
            key: "mul_toom4_threshold",
            op: "mul",
            grid: with_never(geometric(96, 1600, 1.19)),
            window: None,
        },
        Param {
            key: "mul_toom32_threshold",
            op: "mul32",
            grid: with_never(geometric(16, 320, 1.19)),
            window: None,
        },
        Param {
            key: "mul_toom42_threshold",
            op: "mul42",
            grid: with_never(geometric(16, 320, 1.19)),
            window: None,
        },
        Param {
            key: "ntt_efficiency_exponent",
            op: "mul",
            grid: vec![0, 10, 15, 20, 25, 30, 40, 50, 60, 80],
            // Products from the threshold up, where the rule decides.
            window: Some(|c| {
                let from = c["mul_ntt_threshold"].clamp(512, 16384);
                geometric(from, from * 6, 1.1)
            }),
        },
        Param {
            key: "mul_ntt_threshold",
            op: "mul",
            grid: with_never(geometric(512, 16384, 1.19)),
            window: None,
        },
        Param {
            key: "ntt_unbalanced_exponent",
            op: "mul_wide",
            grid: vec![6, 8, 10, 12, 14, 16, 18, 20, 25],
            // Shorter operands around the balanced crossover, where the
            // rule decides between one transform and Toom pieces.
            window: Some(|c| {
                let t = c["mul_ntt_threshold"].clamp(512, 16384);
                geometric(t / 4, t * 2, 1.25)
            }),
        },
        Param {
            key: "sqr_karatsuba_threshold",
            op: "sqr",
            grid: geometric(4, 128, 1.19),
            window: None,
        },
        Param {
            key: "sqr_toom3_threshold",
            op: "sqr",
            grid: with_never(geometric(64, 768, 1.19)),
            window: None,
        },
        Param {
            key: "sqr_toom4_threshold",
            op: "sqr",
            grid: with_never(geometric(96, 1600, 1.19)),
            window: None,
        },
        Param {
            key: "sqr_ntt_threshold",
            op: "sqr",
            grid: with_never(geometric(512, 16384, 1.19)),
            window: None,
        },
        Param {
            key: "ntt_block",
            op: "mul",
            grid: (5..=15).map(|k| 1 << k).collect(),
            window: Some(|c| {
                let from = c["mul_ntt_threshold"].clamp(4096, 16384);
                geometric(from, from * 32, 1.5)
            }),
        },
        Param {
            key: "ntt_wrap_fraction",
            op: "mul",
            grid: vec![1, 2, 4, 8, NEVER],
            // Balanced products just above powers of two, where wrapping
            // applies: len = 2n - 1 runs from 1x to 1.5x a power of two.
            window: Some(|c| {
                let from = c["mul_ntt_threshold"].clamp(4096, 16384);
                let mut v = Vec::new();
                // n >= p / 2 >= the NTT threshold, so every size is an NTT product.
                let mut p = 2 * from.next_power_of_two();
                while v.len() < 12 {
                    v.extend([p / 2 + p / 32, p / 2 + p / 16, p / 2 + p / 8, p / 2 + p / 5]);
                    p *= 2;
                }
                v
            }),
        },
        Param {
            key: "dc_div_threshold",
            op: "div",
            grid: geometric(8, 320, 1.19),
            window: None,
        },
        Param {
            key: "mul_wrap_ntt_threshold",
            op: "div",
            grid: with_never((6..=14).map(|k| 1 << k).collect()),
            // Divisions long enough to take Barrett's estimates, whose
            // products modulo B^m - 1 have every power-of-two m up to n.
            window: Some(|_| geometric(2048, 32768, 1.5)),
        },
        Param {
            key: "inv_newton_threshold",
            op: "div",
            grid: geometric(8, 2400, 1.3),
            window: Some(|_| geometric(2048, 32768, 1.5)),
        },
        Param {
            key: "mu_div_threshold",
            op: "div",
            grid: with_never(geometric(256, 16384, 1.19)),
            window: None,
        },
        // Hensel division with remainder, timed where exact division with
        // its check uses it (all sizes when the base disables
        // divexact_div_threshold).
        Param {
            key: "bdiv_qr_dc_threshold",
            op: "divexact",
            grid: geometric(8, 640, 1.25),
            window: None,
        },
        Param {
            key: "bdiv_q_dc_threshold",
            op: "bdiv_q",
            grid: with_never(geometric(16, 1600, 1.19)),
            window: None,
        },
        Param {
            key: "bdiv_mu_threshold",
            op: "bdiv_q",
            grid: with_never(geometric(64, 8192, 1.19)),
            window: None,
        },
        Param {
            key: "binv_newton_threshold",
            op: "bdiv_q",
            grid: geometric(8, 2000, 1.25),
            // Exact divisions long enough to work in blocks with the inverse.
            window: Some(|c| {
                let from = c["bdiv_mu_threshold"].clamp(64, 4096);
                geometric(from, from * 16, 1.4)
            }),
        },
        Param {
            key: "divexact_div_threshold",
            op: "divexact",
            grid: with_never(geometric(16, 6400, 1.19)),
            window: None,
        },
        Param {
            key: "mullo_dc_threshold",
            op: "mullo",
            grid: geometric(8, 320, 1.19),
            window: None,
        },
        Param {
            key: "mullo_mul_threshold",
            op: "mullo",
            grid: with_never(geometric(512, 16384, 1.19)),
            window: None,
        },
        Param {
            key: "redc_n_threshold",
            op: "redc",
            grid: with_never(geometric(16, 640, 1.19)),
            window: None,
        },
        Param {
            key: "matrix22_strassen_threshold",
            op: "gcd",
            grid: with_never(geometric(4, 600, 1.25)),
            // Sizes whose half-GCD composes matrices of all entry lengths.
            window: Some(|c| {
                let from = c["gcd_dc_threshold"].clamp(64, 1000);
                geometric(from, from.max(1600) * 5 / 2, 1.25)
            }),
        },
        Param {
            key: "hgcd_threshold",
            op: "gcd",
            grid: geometric(16, 400, 1.19),
            // Sizes that take half-GCD steps, where the recursion runs.
            window: Some(|c| {
                let from = c["gcd_dc_threshold"].clamp(64, 1000);
                geometric(from, from.max(1600) * 5 / 2, 1.25)
            }),
        },
        Param {
            key: "gcd_dc_threshold",
            op: "gcd",
            grid: with_never(geometric(32, 1000, 1.19)),
            window: None,
        },
        Param {
            key: "to_string_dc_threshold",
            op: "to_string",
            grid: geometric(2, 160, 1.25),
            window: None,
        },
        Param {
            key: "from_string_dc_threshold",
            op: "from_str",
            grid: geometric(100, 8000, 1.25),
            window: None,
        },
    ]
}

/// A configuration: the base file's text, and its values by key.
#[derive(Clone)]
struct Config {
    text: String,
    values: BTreeMap<String, u64>,
}

impl std::ops::Index<&str> for Config {
    type Output = u64;

    fn index(&self, key: &str) -> &u64 {
        &self.values[key]
    }
}

impl Config {
    fn parse(text: String) -> Config {
        let values = text
            .lines()
            .filter_map(|line| {
                let line = line.split('#').next()?.trim();
                let (k, v) = line.split_once('=')?;
                Some((k.trim().to_owned(), v.trim().replace('_', "").parse().ok()?))
            })
            .collect();
        Config { text, values }
    }

    /// The text with `key`'s value replaced.
    fn set(&mut self, key: &str, value: u64) {
        let mut out = String::new();
        for line in self.text.lines() {
            let code = line.split('#').next().unwrap_or("");
            match code.split_once('=') {
                Some((k, _)) if k.trim() == key => writeln!(out, "{key} = {value}"),
                _ => writeln!(out, "{line}"),
            }
            .expect("writing to a String");
        }
        self.text = out;
        self.values.insert(key.to_owned(), value);
    }
}

struct Tuner {
    root: PathBuf,
    target: PathBuf,
    min_us: u64,
    builds: usize,
    /// The median disagreement between the two passes, per ranking.
    noise: Vec<f64>,
}

impl Tuner {
    /// Rebuilds the probe against `config`; returns a copy of the binary, so
    /// it can be run again after later builds.
    fn build(&mut self, config: &Config) -> PathBuf {
        let file = self.target.join("candidate.toml");
        fs::write(&file, &config.text).expect("writing the candidate configuration");
        let cargo = env::var_os("CARGO").unwrap_or_else(|| "cargo".into());
        let status = Command::new(cargo)
            .current_dir(&self.root)
            .env("CAS_NUM_TUNING", &file)
            .args([
                "build",
                "--release",
                "-q",
                "-p",
                "cas-tune",
                "--bin",
                "cas-tune-probe",
            ])
            .arg("--target-dir")
            .arg(&self.target)
            .status()
            .expect("running cargo");
        assert!(status.success(), "building the probe failed");
        self.builds += 1;
        let copy = self.target.join(format!("probe-{}", self.builds));
        fs::copy(self.target.join("release").join("cas-tune-probe"), &copy)
            .expect("copying the probe");
        copy
    }

    /// Times `op` at `sizes` with a built probe.
    fn run(&self, probe: &Path, op: &str, sizes: &[u64]) -> Vec<f64> {
        let output = Command::new(probe)
            .arg(op)
            .arg(self.min_us.to_string())
            .args(sizes.iter().map(u64::to_string))
            .output()
            .expect("running the probe");
        assert!(output.status.success(), "the probe failed");
        String::from_utf8(output.stdout)
            .expect("probe output")
            .lines()
            .map(|l| l.split_once(' ').expect("`size ns`").1.parse().expect("ns"))
            .collect()
    }

    /// Times every candidate and returns `(value, cost)` pairs, lowest cost
    /// first, with the costs relative to the best time at each size. Each
    /// candidate is timed twice, the passes in opposite orders, keeping the
    /// faster time per size, so drift during the run favours none of them.
    fn rank(
        &mut self,
        config: &Config,
        param: &Param,
        candidates: &[u64],
        sizes: &[u64],
    ) -> Vec<(u64, f64)> {
        let probes: Vec<PathBuf> = candidates
            .iter()
            .map(|&v| {
                let mut c = config.clone();
                c.set(param.key, v);
                self.build(&c)
            })
            .collect();
        let mut passes = [
            vec![Vec::new(); candidates.len()],
            vec![Vec::new(); candidates.len()],
        ];
        for (pass, runs) in passes.iter_mut().enumerate() {
            for k in 0..candidates.len() {
                let k = if pass == 0 {
                    k
                } else {
                    candidates.len() - 1 - k
                };
                runs[k] = self.run(&probes[k], param.op, sizes);
            }
        }
        // The two passes of the same binary disagree only by noise.
        let mut spread: Vec<f64> = passes[0]
            .iter()
            .flatten()
            .zip(passes[1].iter().flatten())
            .map(|(a, b)| (a - b).abs() / a.min(*b))
            .collect();
        spread.sort_by(f64::total_cmp);
        self.noise.push(spread[spread.len() / 2]);
        let times: Vec<Vec<f64>> = passes[0]
            .iter()
            .zip(&passes[1])
            .map(|(a, b)| a.iter().zip(b).map(|(a, b)| a.min(*b)).collect())
            .collect();
        for probe in &probes {
            fs::remove_file(probe).expect("removing a probe copy");
        }
        let best: Vec<f64> = (0..sizes.len())
            .map(|i| times.iter().map(|t| t[i]).fold(f64::INFINITY, f64::min))
            .collect();
        let mut ranked: Vec<(u64, f64)> = candidates
            .iter()
            .zip(&times)
            .map(|(&v, t)| {
                let cost = t.iter().zip(&best).map(|(t, b)| t / b).sum::<f64>();
                (v, cost / sizes.len() as f64)
            })
            .collect();
        ranked.sort_by(|a, b| a.1.total_cmp(&b.1));
        ranked
    }
}

/// The one-minute load average, where the system reports one (Linux).
fn load_average() -> Option<f64> {
    fs::read_to_string("/proc/loadavg")
        .ok()?
        .split_whitespace()
        .next()?
        .parse()
        .ok()
}

fn show(v: u64) -> String {
    if v == NEVER {
        "never".to_owned()
    } else {
        v.to_string()
    }
}

/// Tunes one threshold in `config`; returns the report line.
fn tune(tuner: &mut Tuner, config: &mut Config, param: &Param, quick: bool) -> String {
    let t0 = Instant::now();
    let sizes = param.sizes(config);
    let mut grid: Vec<u64> = if quick {
        let every_other = param.grid.iter().copied().step_by(2);
        every_other.chain(param.grid.last().copied()).collect()
    } else {
        param.grid.clone()
    };
    grid.dedup();
    let coarse = tuner.rank(config, param, &grid, &sizes);
    let best = coarse[0].0;
    let mut fine = Vec::new();
    if !quick && best != NEVER {
        // Refine between the best grid point and its neighbours.
        let at = grid.iter().position(|&v| v == best).expect("on the grid");
        let mut refined: Vec<u64> = [at.checked_sub(1), Some(at + 1)]
            .into_iter()
            .flatten()
            .filter_map(|i| grid.get(i))
            .filter(|&&v| v != NEVER)
            .map(|&v| v.midpoint(best))
            .filter(|v| !grid.contains(v))
            .collect();
        if !refined.is_empty() {
            refined.push(best);
            fine = tuner.rank(config, param, &refined, &sizes);
        }
    }
    let chosen = fine.first().unwrap_or(&coarse[0]).0;
    let previous = config[param.key];
    config.set(param.key, chosen);
    let costs = |ranked: &[(u64, f64)]| {
        ranked
            .iter()
            .map(|(v, c)| format!("{}={c:.3}", show(*v)))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let noise = tuner.noise.last().copied().unwrap_or(0.0);
    let mut line = format!(
        "{}: {} (was {}), {:.1} s, sizes {}..={}, noise {:.1}%\n  grid:   {}",
        param.key,
        show(chosen),
        show(previous),
        t0.elapsed().as_secs_f64(),
        sizes[0],
        sizes[sizes.len() - 1],
        100.0 * noise,
        costs(&coarse)
    );
    if !fine.is_empty() {
        write!(line, "\n  refine: {}", costs(&fine)).expect("writing to a String");
    }
    line
}

fn main() {
    let mut args = env::args().skip(1);
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut base = root.join("crates/cas-num/tuning/default.toml");
    let (mut out, mut quick, mut only) = (None, false, None::<Vec<String>>);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--quick" => quick = true,
            "--base" => base = args.next().expect("--base FILE").into(),
            "--out" => out = Some(PathBuf::from(args.next().expect("--out FILE"))),
            "--only" => {
                only = Some(
                    args.next()
                        .expect("--only key,key")
                        .split(',')
                        .map(str::to_owned)
                        .collect(),
                );
            }
            _ => panic!(
                "usage: cas-tune [--quick] [--only key,key] [--base FILE] --out FILE \
                 (unknown argument `{arg}`)"
            ),
        }
    }
    let out = out.expect("--out FILE is required");
    let mut config = Config::parse(fs::read_to_string(&base).expect("reading the base"));
    let target = root.join("target/cas-tune");
    fs::create_dir_all(&target).expect("creating target/cas-tune");
    let mut tuner = Tuner {
        root,
        target,
        min_us: if quick { 1000 } else { 4000 },
        builds: 0,
        noise: Vec::new(),
    };

    let started = Instant::now();
    let mut report = String::new();
    if let Some(load) = load_average().filter(|&l| l > 1.5) {
        let warning = format!(
            "warning: load average {load:.1}: other work is running, so timings will be noisy \
             and the chosen values unreliable; tune on an idle machine\n"
        );
        print!("{warning}");
        report.push_str(&warning);
    }
    for param in params() {
        if only
            .as_ref()
            .is_some_and(|o| !o.iter().any(|k| k == param.key))
        {
            continue;
        }
        let line = tune(&mut tuner, &mut config, &param, quick);
        println!("{line}");
        report.push_str(&line);
        report.push('\n');
    }
    let mut noise = tuner.noise.clone();
    noise.sort_by(f64::total_cmp);
    let noise = noise.get(noise.len() / 2).copied().unwrap_or(0.0);
    let header = format!(
        "# Written by cas-tune in {:.0} s ({} builds), median timing noise {:.1}%.\n",
        started.elapsed().as_secs_f64(),
        tuner.builds,
        100.0 * noise
    );
    if noise > 0.03 {
        let warning = format!(
            "warning: identical runs differed by {:.1}% (median), as much as most \
             candidates differ; rerun on an idle machine\n",
            100.0 * noise
        );
        print!("{warning}");
        report.push_str(&warning);
    }
    let mut text = header.clone();
    text.push_str(&config.text);
    fs::write(&out, text).expect("writing the configuration");
    let report_path = out.with_extension("report.txt");
    fs::write(&report_path, header + &report).expect("writing the report");
    println!(
        "wrote {} and {} in {:.0} s, {} builds",
        out.display(),
        report_path.display(),
        started.elapsed().as_secs_f64(),
        tuner.builds
    );
}
