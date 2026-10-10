//! Dense column-major matrices.

use std::fmt;
use std::ops::{Index, IndexMut};

use crate::vector::{axpy, dot, norm2};

/// A dense `rows × cols` matrix of `f64`, stored column by column.
#[derive(Clone, PartialEq)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

impl Matrix {
    /// The `rows × cols` zero matrix.
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Matrix {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }

    /// The `n × n` identity.
    pub fn identity(n: usize) -> Self {
        let mut m = Matrix::zeros(n, n);
        for i in 0..n {
            m[(i, i)] = 1.0;
        }
        m
    }

    /// A matrix whose entry `(i, j)` is `f(i, j)`.
    pub fn from_fn(rows: usize, cols: usize, mut f: impl FnMut(usize, usize) -> f64) -> Self {
        let mut data = Vec::with_capacity(rows * cols);
        for j in 0..cols {
            for i in 0..rows {
                data.push(f(i, j));
            }
        }
        Matrix { rows, cols, data }
    }

    /// A matrix from entries listed row by row.
    ///
    /// # Panics
    /// If `values.len() != rows * cols`.
    pub fn from_row_major(rows: usize, cols: usize, values: &[f64]) -> Self {
        assert_eq!(values.len(), rows * cols, "from_row_major: wrong length");
        Matrix::from_fn(rows, cols, |i, j| values[i * cols + j])
    }

    /// A matrix from entries listed column by column.
    ///
    /// # Panics
    /// If `values.len() != rows * cols`.
    pub fn from_col_major(rows: usize, cols: usize, values: Vec<f64>) -> Self {
        assert_eq!(values.len(), rows * cols, "from_col_major: wrong length");
        Matrix {
            rows,
            cols,
            data: values,
        }
    }

    /// Number of rows.
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns.
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// The entries, column by column.
    pub fn as_col_major(&self) -> &[f64] {
        &self.data
    }

    /// Column `j`.
    pub fn col(&self, j: usize) -> &[f64] {
        &self.data[j * self.rows..(j + 1) * self.rows]
    }

    /// Column `j`, mutably.
    pub fn col_mut(&mut self, j: usize) -> &mut [f64] {
        &mut self.data[j * self.rows..(j + 1) * self.rows]
    }

    /// Two distinct columns, mutably (`i != j`).
    ///
    /// # Panics
    /// If `i == j`.
    pub fn two_cols_mut(&mut self, i: usize, j: usize) -> (&mut [f64], &mut [f64]) {
        assert_ne!(i, j, "two_cols_mut: same column");
        let r = self.rows;
        if i < j {
            let (a, b) = self.data.split_at_mut(j * r);
            (&mut a[i * r..(i + 1) * r], &mut b[..r])
        } else {
            let (a, b) = self.data.split_at_mut(i * r);
            let (bj, bi) = (&mut a[j * r..(j + 1) * r], &mut b[..r]);
            (bi, bj)
        }
    }

    /// Swaps columns `i` and `j`.
    pub fn swap_cols(&mut self, i: usize, j: usize) {
        if i != j {
            let (a, b) = self.two_cols_mut(i, j);
            a.swap_with_slice(b);
        }
    }

    /// Row `i` as a new vector.
    pub fn row(&self, i: usize) -> Vec<f64> {
        (0..self.cols).map(|j| self[(i, j)]).collect()
    }

    /// The transpose.
    pub fn transpose(&self) -> Matrix {
        Matrix::from_fn(self.cols, self.rows, |i, j| self[(j, i)])
    }

    /// The product `self · other`.
    ///
    /// # Panics
    /// If the inner dimensions differ.
    pub fn matmul(&self, other: &Matrix) -> Matrix {
        assert_eq!(self.cols, other.rows, "matmul: dimension mismatch");
        let mut out = Matrix::zeros(self.rows, other.cols);
        for j in 0..other.cols {
            let bj = other.col(j);
            let oj = &mut out.data[j * self.rows..(j + 1) * self.rows];
            for (k, &b) in bj.iter().enumerate() {
                if b != 0.0 {
                    axpy(b, &self.data[k * self.rows..(k + 1) * self.rows], oj);
                }
            }
        }
        out
    }

    /// The product `selfᵀ · other`.
    ///
    /// # Panics
    /// If the row counts differ.
    pub fn tr_matmul(&self, other: &Matrix) -> Matrix {
        assert_eq!(self.rows, other.rows, "tr_matmul: dimension mismatch");
        Matrix::from_fn(self.cols, other.cols, |i, j| dot(self.col(i), other.col(j)))
    }

    /// The product `self · x`.
    ///
    /// # Panics
    /// If `x.len() != self.cols()`.
    pub fn mul_vec(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.cols, "mul_vec: dimension mismatch");
        let mut y = vec![0.0; self.rows];
        for (j, &xj) in x.iter().enumerate() {
            if xj != 0.0 {
                axpy(xj, self.col(j), &mut y);
            }
        }
        y
    }

    /// The product `selfᵀ · y`.
    ///
    /// # Panics
    /// If `y.len() != self.rows()`.
    pub fn tr_mul_vec(&self, y: &[f64]) -> Vec<f64> {
        assert_eq!(y.len(), self.rows, "tr_mul_vec: dimension mismatch");
        (0..self.cols).map(|j| dot(self.col(j), y)).collect()
    }

    /// Multiplies column `j` by `s[j]`.
    ///
    /// # Panics
    /// If `s.len() != self.cols()`.
    pub fn scale_cols(&mut self, s: &[f64]) {
        assert_eq!(s.len(), self.cols, "scale_cols: dimension mismatch");
        for (j, &sj) in s.iter().enumerate() {
            for v in self.col_mut(j) {
                *v *= sj;
            }
        }
    }

    /// The matrix formed by the listed columns, in that order.
    pub fn select_cols(&self, cols: &[usize]) -> Matrix {
        let mut data = Vec::with_capacity(self.rows * cols.len());
        for &j in cols {
            data.extend_from_slice(self.col(j));
        }
        Matrix {
            rows: self.rows,
            cols: cols.len(),
            data,
        }
    }

    /// The matrix formed by the listed rows, in that order.
    pub fn select_rows(&self, rows: &[usize]) -> Matrix {
        Matrix::from_fn(rows.len(), self.cols, |i, j| self[(rows[i], j)])
    }

    /// Frobenius norm.
    pub fn frobenius_norm(&self) -> f64 {
        norm2(&self.data)
    }

    /// Largest entry magnitude.
    pub fn max_abs(&self) -> f64 {
        crate::vector::norm_inf(&self.data)
    }

    /// True if every entry is finite.
    pub fn is_finite(&self) -> bool {
        self.data.iter().all(|v| v.is_finite())
    }
}

impl Index<(usize, usize)> for Matrix {
    type Output = f64;
    fn index(&self, (i, j): (usize, usize)) -> &f64 {
        debug_assert!(i < self.rows && j < self.cols, "index out of bounds");
        &self.data[j * self.rows + i]
    }
}

impl IndexMut<(usize, usize)> for Matrix {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut f64 {
        debug_assert!(i < self.rows && j < self.cols, "index out of bounds");
        &mut self.data[j * self.rows + i]
    }
}

impl fmt::Debug for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Matrix {}x{} [", self.rows, self.cols)?;
        for i in 0..self.rows {
            write!(f, "  ")?;
            for j in 0..self.cols {
                write!(f, "{:>12.5e} ", self[(i, j)])?;
            }
            writeln!(f)?;
        }
        write!(f, "]")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn construction_and_products() {
        let a = Matrix::from_row_major(2, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(a[(1, 0)], 4.0);
        assert_eq!(a.col(2), &[3.0, 6.0]);
        assert_eq!(a.row(0), vec![1.0, 2.0, 3.0]);
        let at = a.transpose();
        assert_eq!(at[(2, 1)], 6.0);
        let p = a.matmul(&at);
        assert_eq!(p, Matrix::from_row_major(2, 2, &[14.0, 32.0, 32.0, 77.0]));
        assert_eq!(a.tr_matmul(&a), at.matmul(&a));
        assert_eq!(a.mul_vec(&[1.0, 0.0, -1.0]), vec![-2.0, -2.0]);
        assert_eq!(a.tr_mul_vec(&[1.0, 1.0]), vec![5.0, 7.0, 9.0]);
        assert_eq!(Matrix::identity(2).matmul(&a), a);
    }

    #[test]
    fn column_helpers() {
        let mut a = Matrix::from_row_major(2, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        a.swap_cols(0, 2);
        assert_eq!(a.row(0), vec![3.0, 2.0, 1.0]);
        let (x, y) = a.two_cols_mut(2, 1);
        x[0] = 10.0;
        y[0] = 20.0;
        assert_eq!(a.row(0), vec![3.0, 20.0, 10.0]);
        a.scale_cols(&[1.0, 0.5, 2.0]);
        assert_eq!(a.row(1), vec![6.0, 2.5, 8.0]);
        assert_eq!(a.select_cols(&[2, 0]).row(0), vec![20.0, 3.0]);
        assert_eq!(a.select_rows(&[1]).row(0), vec![6.0, 2.5, 8.0]);
        assert_eq!(a.max_abs(), 20.0);
        assert!(a.is_finite());
    }
}
