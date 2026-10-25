//! LLL reduction of integer lattices, in exact integer arithmetic: Cohen's
//! integral LLL \[Coh93 Algorithm 2.6.7\], after de Weger, with `δ = 3/4`.
//! van Hoeij's recombination (CAS-188) reduces its knapsack lattices with
//! it.

#![allow(
    clippy::many_single_char_names,
    reason = "the names of Cohen's algorithm: b, d, λ, k, l, q"
)]

use cas_num::Integer;

use crate::cancel::{Cancel, Cancelled};

/// A basis being reduced: its vectors, and the integral Gram-Schmidt data
/// of the first `k_max` of them. `d[i]` is the Gram determinant of the
/// first `i` vectors, so that `‖b*ᵢ‖² = d[i + 1] / d[i]` (0-based), and
/// `lambda[i][j] = d[j + 1] μᵢⱼ`.
pub(crate) struct Lattice {
    basis: Vec<Vec<Integer>>,
    d: Vec<Integer>,
    lambda: Vec<Vec<Integer>>,
}

fn dot(a: &[Integer], b: &[Integer]) -> Integer {
    a.iter()
        .zip(b)
        .fold(Integer::ZERO, |acc, (x, y)| &acc + &(x * y))
}

/// The integer nearest to `a / b`, `b > 0`, halves rounded up.
fn nearest(a: &Integer, b: &Integer) -> Integer {
    let twice = a + a;
    let (q, _) = (&twice + b).div_rem_floor(&(b + b));
    q
}

impl Lattice {
    /// The LLL-reduced basis of the lattice the rows of `basis` generate,
    /// which must be linearly independent, polling `cancel` at each step.
    ///
    /// # Panics
    ///
    /// If the rows are linearly dependent.
    pub(crate) fn reduce(
        basis: Vec<Vec<Integer>>,
        cancel: &dyn Cancel,
    ) -> Result<Lattice, Cancelled> {
        let n = basis.len();
        let mut lattice = Lattice {
            basis,
            d: vec![Integer::ONE; n + 1],
            lambda: vec![vec![Integer::ZERO; n]; n],
        };
        if n == 0 {
            return Ok(lattice);
        }
        lattice.d[1] = dot(&lattice.basis[0], &lattice.basis[0]);
        let (mut k, mut k_max) = (1, 0);
        while k < n {
            if cancel.cancelled() {
                return Err(Cancelled);
            }
            if k > k_max {
                k_max = k;
                lattice.gram_schmidt(k);
            }
            lattice.reduce_one(k, k - 1);
            if lattice.lovasz_fails(k) {
                lattice.swap(k, k_max);
                k = k.saturating_sub(1).max(1);
            } else {
                for l in (0..k - 1).rev() {
                    lattice.reduce_one(k, l);
                }
                k += 1;
            }
        }
        Ok(lattice)
    }

    /// The incremental Gram-Schmidt step for vector `k`.
    fn gram_schmidt(&mut self, k: usize) {
        for j in 0..=k {
            let mut u = dot(&self.basis[k], &self.basis[j]);
            for i in 0..j {
                u = (&(&self.d[i + 1] * &u) - &(&self.lambda[k][i] * &self.lambda[j][i]))
                    .div_exact(&self.d[i]);
            }
            if j < k {
                self.lambda[k][j] = u;
            } else {
                assert!(!u.is_zero(), "linearly independent rows");
                self.d[k + 1] = u;
            }
        }
    }

    /// Size reduction of vector `k` against vector `l < k`.
    fn reduce_one(&mut self, k: usize, l: usize) {
        let d = &self.d[l + 1];
        let twice = &self.lambda[k][l] + &self.lambda[k][l];
        if twice.unsigned_abs() <= d.unsigned_abs() {
            return;
        }
        let q = nearest(&self.lambda[k][l], d);
        let row = self.basis[l].clone();
        for (x, y) in self.basis[k].iter_mut().zip(&row) {
            *x = &*x - &(&q * y);
        }
        self.lambda[k][l] = &self.lambda[k][l] - &(&q * d);
        for i in 0..l {
            let t = &q * &self.lambda[l][i];
            self.lambda[k][i] = &self.lambda[k][i] - &t;
        }
    }

    /// Whether the Lovász condition with `δ = 3/4` fails at `k`:
    /// `4 d_k d_(k−2) < 3 d_(k−1)² − 4 λ²` in Cohen's 1-based terms.
    fn lovasz_fails(&self, k: usize) -> bool {
        let lambda = &self.lambda[k][k - 1];
        let left = &Integer::from(4) * &(&self.d[k + 1] * &self.d[k - 1]);
        let right = &(&Integer::from(3) * &(&self.d[k] * &self.d[k]))
            - &(&Integer::from(4) * &(lambda * lambda));
        left < right
    }

    /// Exchanges vectors `k − 1` and `k`, updating the data up to `k_max`.
    fn swap(&mut self, k: usize, k_max: usize) {
        self.basis.swap(k, k - 1);
        for j in 0..k - 1 {
            let t = std::mem::take(&mut self.lambda[k][j]);
            self.lambda[k][j] = std::mem::replace(&mut self.lambda[k - 1][j], t);
        }
        let lambda = self.lambda[k][k - 1].clone();
        let b = (&(&self.d[k - 1] * &self.d[k + 1]) + &(&lambda * &lambda)).div_exact(&self.d[k]);
        for i in k + 1..=k_max {
            let t = self.lambda[i][k].clone();
            self.lambda[i][k] = (&(&self.d[k + 1] * &self.lambda[i][k - 1]) - &(&lambda * &t))
                .div_exact(&self.d[k]);
            self.lambda[i][k - 1] =
                (&(&b * &t) + &(&lambda * &self.lambda[i][k])).div_exact(&self.d[k + 1]);
        }
        self.d[k] = b;
    }

    /// The reduced basis.
    pub(crate) fn basis(&self) -> &[Vec<Integer>] {
        &self.basis
    }

    pub(crate) fn into_basis(self) -> Vec<Vec<Integer>> {
        self.basis
    }

    /// Whether `‖b*ᵢ‖² > bound` for the `i`-th vector, 0-based: `d_(i+1) >
    /// bound d_i`.
    pub(crate) fn gram_schmidt_exceeds(&self, i: usize, bound: &Integer) -> bool {
        self.d[i + 1] > &self.d[i] * bound
    }
}
