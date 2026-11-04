//! Certified eigenpairs of numeric matrices (CAS-DOC-51 §4, CAS-498), after
//! Rump's `verifyeig` \[Rum01\].
//!
//! For a matrix of balls `A` and an approximate eigenvalue `λ̃`:
//!
//! 1. An approximate eigenvector `x̃` by inverse iteration on the midpoints,
//!    then Newton's method on `F(x, λ) = A x − λ x` with `x`'s largest
//!    component fixed at 1, `k` its index, all in midpoint arithmetic.
//! 2. The unknowns `y`: `x_j` for `j ≠ k`, and `λ` in place `k`. Exactly,
//!    `F(ỹ + y) = F(ỹ) + J̃(y) y`, with `J̃(y)` the matrix `A − λ̃ I` whose
//!    column `k` is replaced by `−(x̃ + y_x)`. With `R` an approximate
//!    inverse of `J̃(0)`'s midpoint, Krawczyk's operator `K(Y) = −R F(ỹ) +
//!    (I − R J̃(Y)) Y` maps a box `Y` around 0 into its interior only if
//!    `ỹ + Y` holds one zero of `F` for each matrix in `A`'s balls, which
//!    is then in `ỹ + K(Y)`: an eigenpair, the eigenvalue simple.
//!
//! `Y` starts from `|R F(ỹ)|` inflated (ε-inflation) and grows a few times
//! before the check gives up: a multiple eigenvalue, or a cluster closer
//! than the working precision tells, has no such box. A real `A` with a
//! real `λ̃` keeps every box real: the eigenpair is then real.

use cas_num::Integer;
use cas_num::ball::{Ball, ComplexBall, Float, Mag};

/// The tries at a box before giving up.
const INFLATIONS: usize = 6;

/// A ball's midpoint, exact.
fn mid(z: &ComplexBall) -> ComplexBall {
    ComplexBall::new(
        Ball::exact(z.re().mid().clone()),
        Ball::exact(z.im().mid().clone()),
    )
}

fn zero() -> ComplexBall {
    ComplexBall::from_real(Ball::from_integer(&Integer::ZERO))
}

fn one() -> ComplexBall {
    ComplexBall::from_real(Ball::from_integer(&Integer::ONE))
}

/// The larger of the midpoint's parts in magnitude, for pivoting.
fn size(z: &ComplexBall) -> Float {
    z.re().mid().abs().max(z.im().mid().abs())
}

/// An upper bound of `|x|` over the real ball.
fn upper(b: &Ball) -> Mag {
    Mag::from_float(b.mid()).add(b.rad())
}

/// Arithmetic at `prec` bits: on midpoints for the approximations, on balls
/// for the check.
struct Arith {
    prec: u64,
}

impl Arith {
    /// `m⁻¹ b` on midpoints, by Gaussian elimination with partial pivoting;
    /// `None` when a pivot is 0.
    fn solve(&self, m: &[Vec<ComplexBall>], b: &[ComplexBall]) -> Option<Vec<ComplexBall>> {
        let n = b.len();
        let mut a: Vec<Vec<ComplexBall>> = m
            .iter()
            .zip(b)
            .map(|(row, y)| row.iter().chain([y]).map(mid).collect())
            .collect();
        for k in 0..n {
            let best = (k..n).max_by(|&i, &j| size(&a[i][k]).cmp(&size(&a[j][k])))?;
            if size(&a[best][k]).is_zero() {
                return None;
            }
            a.swap(k, best);
            for i in k + 1..n {
                let c = mid(&a[i][k].div(&a[k][k], self.prec));
                for j in k..=n {
                    let t = mid(&c.mul(&a[k][j], self.prec));
                    a[i][j] = mid(&a[i][j].sub(&t, self.prec));
                }
            }
        }
        let mut x = vec![zero(); n];
        for k in (0..n).rev() {
            let mut s = a[k][n].clone();
            for j in k + 1..n {
                s = mid(&s.sub(&a[k][j].mul(&x[j], self.prec), self.prec));
            }
            x[k] = mid(&s.div(&a[k][k], self.prec));
        }
        Some(x)
    }

    /// The approximate inverse of `m`'s midpoint, column by column.
    fn inverse(&self, m: &[Vec<ComplexBall>]) -> Option<Vec<Vec<ComplexBall>>> {
        let n = m.len();
        let mut columns = Vec::with_capacity(n);
        for j in 0..n {
            let e: Vec<ComplexBall> = (0..n)
                .map(|i| if i == j { one() } else { zero() })
                .collect();
            columns.push(self.solve(m, &e)?);
        }
        Some(
            (0..n)
                .map(|i| columns.iter().map(|c| c[i].clone()).collect())
                .collect(),
        )
    }

    /// `m v` in balls.
    fn apply(&self, m: &[Vec<ComplexBall>], v: &[ComplexBall]) -> Vec<ComplexBall> {
        m.iter()
            .map(|row| {
                row.iter()
                    .zip(v)
                    .fold(zero(), |s, (a, x)| s.add(&a.mul(x, self.prec), self.prec))
            })
            .collect()
    }

    /// `m p` in balls.
    fn times(&self, m: &[Vec<ComplexBall>], p: &[Vec<ComplexBall>]) -> Vec<Vec<ComplexBall>> {
        let n = p.first().map_or(0, Vec::len);
        m.iter()
            .map(|row| {
                (0..n)
                    .map(|j| {
                        row.iter().zip(p).fold(zero(), |s, (a, r)| {
                            s.add(&a.mul(&r[j], self.prec), self.prec)
                        })
                    })
                    .collect()
            })
            .collect()
    }

    /// `A − λ I` with column `k` replaced by `−x`.
    fn jacobian(
        &self,
        a: &[Vec<ComplexBall>],
        lambda: &ComplexBall,
        x: &[ComplexBall],
        k: usize,
    ) -> Vec<Vec<ComplexBall>> {
        a.iter()
            .enumerate()
            .map(|(i, row)| {
                row.iter()
                    .enumerate()
                    .map(|(j, aij)| {
                        if j == k {
                            x[i].neg()
                        } else if i == j {
                            aij.sub(lambda, self.prec)
                        } else {
                            aij.clone()
                        }
                    })
                    .collect()
            })
            .collect()
    }

    /// `A x − λ x` in balls.
    fn residual(
        &self,
        a: &[Vec<ComplexBall>],
        lambda: &ComplexBall,
        x: &[ComplexBall],
    ) -> Vec<ComplexBall> {
        self.apply(a, x)
            .into_iter()
            .zip(x)
            .map(|(ax, xi)| ax.sub(&lambda.mul(xi, self.prec), self.prec))
            .collect()
    }
}

/// The index of `x`'s largest component, and `x` scaled to 1 there.
fn normalized(arith: &Arith, x: &[ComplexBall]) -> Option<(usize, Vec<ComplexBall>)> {
    let k = (0..x.len()).max_by(|&i, &j| size(&x[i]).cmp(&size(&x[j])))?;
    if size(&x[k]).is_zero() {
        return None;
    }
    let pivot = x[k].clone();
    let mut y: Vec<ComplexBall> = x.iter().map(|z| mid(&z.div(&pivot, arith.prec))).collect();
    y[k] = one();
    Some((k, y))
}

/// The box of half-widths `rho` around 0, real when `real`.
fn centred(rho: &[Mag], real: bool) -> Vec<ComplexBall> {
    rho.iter()
        .map(|r| {
            let part = Ball::new(Float::ZERO, r.clone());
            let im = if real {
                Ball::from_integer(&Integer::ZERO)
            } else {
                part.clone()
            };
            ComplexBall::new(part, im)
        })
        .collect()
}

/// Whether the ball `z` lies within the half-width `r`, strictly.
fn inside(z: &ComplexBall, r: &Mag) -> bool {
    upper(z.re()) < *r && upper(z.im()) < *r
}

/// A certified eigenpair of `a` near `approximate`: the eigenvalue and an
/// eigenvector with 1 at its largest component, as balls; `None` when the
/// check fails.
pub(super) fn eigenpair(
    a: &[Vec<ComplexBall>],
    approximate: &ComplexBall,
    prec: u64,
) -> Option<(ComplexBall, Vec<ComplexBall>)> {
    let arith = Arith { prec };
    let n = a.len();
    let real = approximate.im().mid().is_zero()
        && a.iter()
            .flatten()
            .all(|z| z.im().mid().is_zero() && z.im().is_exact());
    let mut lambda = mid(approximate);
    if real {
        lambda = ComplexBall::from_real(lambda.re().clone());
    }
    let mids: Vec<Vec<ComplexBall>> = a.iter().map(|r| r.iter().map(mid).collect()).collect();
    // Inverse iteration, the shift nudged off an exact eigenvalue.
    let nudge = ComplexBall::from_real(Ball::exact(
        Mag::power_of_two(&Integer::from(-i64::try_from(prec / 2).ok()?)).to_float()?,
    ));
    // Column `n` replaced: none, `A − (λ̃ + ν) I`.
    let shifted = arith.jacobian(&mids, &lambda.add(&nudge, prec), &[], n);
    let mut x = vec![one(); n];
    for _ in 0..3 {
        x = arith.solve(&shifted, &x)?;
        x = normalized(&arith, &x)?.1;
    }
    // Newton's method on (x, λ), x_k = 1.
    let (mut k, mut xt) = normalized(&arith, &x)?;
    for _ in 0..4 {
        let j = arith.jacobian(&mids, &lambda, &xt, k);
        let f: Vec<ComplexBall> = arith
            .residual(&mids, &lambda, &xt)
            .iter()
            .map(mid)
            .collect();
        let step = arith.solve(&j, &f)?;
        for (i, s) in step.iter().enumerate() {
            if i == k {
                lambda = mid(&lambda.sub(s, prec));
            } else {
                xt[i] = mid(&xt[i].sub(s, prec));
            }
        }
        (k, xt) = normalized(&arith, &xt)?;
    }
    // Krawczyk's check.
    let r = arith.inverse(&arith.jacobian(&mids, &lambda, &xt, k))?;
    let z: Vec<ComplexBall> = arith
        .apply(&r, &arith.residual(a, &lambda, &xt))
        .iter()
        .map(ComplexBall::neg)
        .collect();
    let floor = Mag::power_of_two(&Integer::from(-i64::try_from(prec).ok()?));
    let mut rho: Vec<Mag> = z
        .iter()
        .map(|w| {
            let m = upper(w.re()).max(upper(w.im()));
            m.add(&m).add(&floor)
        })
        .collect();
    for _ in 0..INFLATIONS {
        let y = centred(&rho, real);
        // x̃ + y_x, its component k fixed.
        let moved: Vec<ComplexBall> = (0..n)
            .map(|i| {
                if i == k {
                    one()
                } else {
                    xt[i].add(&y[i], prec)
                }
            })
            .collect();
        let jy = arith.jacobian(a, &lambda, &moved, k);
        let rj = arith.times(&r, &jy);
        let ry: Vec<ComplexBall> = (0..n)
            .map(|i| {
                let mut s = y[i].clone();
                for (j, yj) in y.iter().enumerate() {
                    s = s.sub(&rj[i][j].mul(yj, prec), prec);
                }
                s
            })
            .collect();
        let kv: Vec<ComplexBall> = z.iter().zip(&ry).map(|(a, b)| a.add(b, prec)).collect();
        if kv.iter().zip(&rho).all(|(w, r)| inside(w, r)) {
            let value = lambda.add(&kv[k], prec);
            let vector = (0..n)
                .map(|i| {
                    if i == k {
                        one()
                    } else {
                        xt[i].add(&kv[i], prec)
                    }
                })
                .collect();
            return Some((value, vector));
        }
        rho = kv
            .iter()
            .zip(&rho)
            .map(|(w, r)| {
                let m = upper(w.re()).max(upper(w.im())).max(r.clone());
                m.add(&m).add(&floor)
            })
            .collect();
    }
    None
}

#[cfg(test)]
mod tests;
