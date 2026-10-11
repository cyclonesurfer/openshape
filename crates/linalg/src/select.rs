//! Choosing a maximal set of linearly independent columns in a given
//! priority order: Householder QR in which the next column in order is
//! accepted when its distance from the span of the columns already accepted
//! exceeds a threshold, and deferred otherwise (Golub & Van Loan, *Matrix
//! Computations*, 4th ed., §5.4.1 with a restricted pivot choice; the
//! matroid greedy algorithm gives the lexicographically first basis).

use crate::householder::{apply_reflector, make_reflector};
use crate::matrix::Matrix;
use crate::vector::norm2;

/// Result of [`select_columns_in_order`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ColumnSelection {
    /// Accepted columns, in increasing index order.
    pub selected: Vec<usize>,
    /// The other columns, in increasing index order.
    pub rejected: Vec<usize>,
}

/// Selects linearly independent columns of `a`, earlier columns first.
///
/// Columns are visited in index order; column `j` is accepted when its
/// distance from the span of the columns accepted so far is greater than
/// `rel_tol · max_k ‖a_k‖₂`, and while fewer than `max_count` (if given) are
/// accepted. If `max_count` is given and the first pass accepted fewer
/// columns, the rejected columns farthest from the current span are added
/// (largest distance first, ties to the lowest index) until `max_count` is
/// reached or every remaining column lies exactly in the span. This lets a
/// caller make the selection agree with a rank decided elsewhere (for
/// example by an SVD).
pub fn select_columns_in_order(
    a: &Matrix,
    rel_tol: f64,
    max_count: Option<usize>,
) -> ColumnSelection {
    let (m, n) = (a.rows(), a.cols());
    let mut w = a.clone();
    let scale = (0..n).map(|j| norm2(a.col(j))).fold(0.0f64, f64::max);
    let cap = max_count.unwrap_or(usize::MAX).min(m).min(n);
    let mut accepted = vec![false; n];
    let mut k = 0usize;
    if scale > 0.0 && scale.is_finite() {
        let thr = rel_tol * scale;
        for j in 0..n {
            if k >= cap {
                break;
            }
            let d = norm2(&w.col(j)[k..]);
            if d > thr {
                accept(&mut w, &mut accepted, j, k);
                k += 1;
            }
        }
        while max_count.is_some() && k < cap {
            let mut best: Option<(usize, f64)> = None;
            for (j, &taken) in accepted.iter().enumerate() {
                if taken {
                    continue;
                }
                let d = norm2(&w.col(j)[k..]);
                if d > 0.0 && best.is_none_or(|(_, bd)| d > bd) {
                    best = Some((j, d));
                }
            }
            match best {
                Some((j, _)) => {
                    accept(&mut w, &mut accepted, j, k);
                    k += 1;
                }
                None => break,
            }
        }
    }
    let selected = (0..n).filter(|&j| accepted[j]).collect();
    let rejected = (0..n).filter(|&j| !accepted[j]).collect();
    ColumnSelection { selected, rejected }
}

/// Accepts column `j` as pivot `k`: makes its reflector and applies it to
/// every column not yet accepted.
fn accept(w: &mut Matrix, accepted: &mut [bool], j: usize, k: usize) {
    let tau = make_reflector(&mut w.col_mut(j)[k..]);
    accepted[j] = true;
    for (c, &taken) in accepted.iter().enumerate() {
        if taken {
            continue;
        }
        let (cj, cc) = w.two_cols_mut(j, c);
        apply_reflector(tau, &cj[k + 1..], &mut cc[k..]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prefers_earlier_columns() {
        // Columns: e0, e0 + e1, e1, 2 e0. Rank 2; the first two are chosen.
        let a = Matrix::from_row_major(2, 4, &[1.0, 1.0, 0.0, 2.0, 0.0, 1.0, 1.0, 0.0]);
        let s = select_columns_in_order(&a, 1e-10, None);
        assert_eq!(s.selected, vec![0, 1]);
        assert_eq!(s.rejected, vec![2, 3]);
    }

    #[test]
    fn cap_and_top_up() {
        let a = Matrix::from_row_major(2, 3, &[1.0, 0.0, 1.0, 0.0, 1e-12, 0.0]);
        // The second column is below the threshold: rejected in pass 1.
        let s = select_columns_in_order(&a, 1e-10, None);
        assert_eq!(s.selected, vec![0]);
        // Forcing two columns adds it back (the third is in the span).
        let s = select_columns_in_order(&a, 1e-10, Some(2));
        assert_eq!(s.selected, vec![0, 1]);
        // Capping at one keeps only the first.
        let s = select_columns_in_order(&a, 0.0, Some(1));
        assert_eq!(s.selected, vec![0]);
    }

    #[test]
    fn zero_matrix_selects_nothing() {
        let s = select_columns_in_order(&Matrix::zeros(3, 2), 1e-10, Some(2));
        assert!(s.selected.is_empty());
        assert_eq!(s.rejected, vec![0, 1]);
    }
}
