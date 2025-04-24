use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::{
    gp::{NGP, NXYZ},
    nature_errors::NErrors,
};

// Trait to define the behavior of a 3x3 matrix
pub trait Mat {
    fn new() -> Self;
    fn new_with_elements(
        a11: f64,
        a12: f64,
        a13: f64,
        a21: f64,
        a22: f64,
        a23: f64,
        a31: f64,
        a32: f64,
        a33: f64,
    ) -> Self;
    fn new_with_columns(col1: &NXYZ, col2: &NXYZ, col3: &NXYZ) -> Self;
    fn set_col(&mut self, col: i32, value: &NXYZ) -> Result<(), NErrors>;
    fn set_cols(&mut self, col1: &NXYZ, col2: &NXYZ, col3: &NXYZ);
    fn set_cross(&mut self, ref_vec: &NXYZ);
    fn set_diagonal(&mut self, x1: f64, x2: f64, x3: f64);
    fn set_dot(&mut self, ref_vec: &NXYZ);
    fn set_identity(&mut self);
    fn set_rotation(&mut self, axis: &NXYZ, ang: f64) -> Result<(), NErrors>;
    fn set_row(&mut self, row: i32, value: &NXYZ) -> Result<(), NErrors>;
    fn set_rows(&mut self, row1: &NXYZ, row2: &NXYZ, row3: &NXYZ);
    fn set_scale(&mut self, s: f64);
    fn set_value(&mut self, row: i32, col: i32, value: f64) -> Result<(), NErrors>;
    fn column(&self, col: i32) -> Result<NXYZ, NErrors>;
    fn determinant(&self) -> f64;
    fn diagonal(&self) -> NXYZ;
    fn row(&self, row: i32) -> Result<NXYZ, NErrors>;
    fn value(&self, row: i32, col: i32) -> Result<f64, NErrors>;
    fn change_value(&mut self, row: i32, col: i32) -> Result<&mut f64, NErrors>;
    fn is_singular(&self) -> bool;
    fn add(&mut self, other: &Self);
    fn added(&self, other: &Self) -> Self
    where
        Self: Sized;
    fn divide(&mut self, scalar: f64) -> Result<(), NErrors>;
    fn divided(&self, scalar: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn invert(&mut self) -> Result<(), NErrors>;
    fn inverted(&self) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn multiply_mat(&mut self, other: &Self);
    fn multiplied_mat(&self, other: &Self) -> Self
    where
        Self: Sized;
    fn pre_multiply(&mut self, other: &Self);
    fn multiply_scalar(&mut self, scalar: f64);
    fn multiplied_scalar(&self, scalar: f64) -> Self
    where
        Self: Sized;
    fn power(&mut self, n: i32) -> Result<(), NErrors>;
    fn powered(&self, n: i32) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn subtract(&mut self, other: &Self);
    fn subtracted(&self, other: &Self) -> Self
    where
        Self: Sized;
    fn transpose(&mut self);
    fn transposed(&self) -> Self
    where
        Self: Sized;
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
}

// Struct representing a 3x3 matrix
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NMat {
    mat: [[f64; 3]; 3],
}

impl Mat for NMat {
    /// Creates a matrix with null coefficients.
    fn new() -> Self {
        NMat { mat: [[0.0; 3]; 3] }
    }

    /// Creates a matrix with specified elements.
    fn new_with_elements(
        a11: f64,
        a12: f64,
        a13: f64,
        a21: f64,
        a22: f64,
        a23: f64,
        a31: f64,
        a32: f64,
        a33: f64,
    ) -> Self {
        NMat {
            mat: [[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]],
        }
    }

    /// Creates a matrix with columns specified by XYZ vectors.
    fn new_with_columns(col1: &NXYZ, col2: &NXYZ, col3: &NXYZ) -> Self {
        NMat {
            mat: [
                [col1.x(), col2.x(), col3.x()],
                [col1.y(), col2.y(), col3.y()],
                [col1.z(), col2.z(), col3.z()],
            ],
        }
    }

    /// Assigns the coordinates of value to the column of index col.
    fn set_col(&mut self, col: i32, value: &NXYZ) -> Result<(), NErrors> {
        if col < 1 || col > 3 {
            return Err(NErrors::OutOfRange);
        }
        let col_idx = (col - 1) as usize;
        self.mat[0][col_idx] = value.x();
        self.mat[1][col_idx] = value.y();
        self.mat[2][col_idx] = value.z();
        Ok(())
    }

    /// Assigns the three columns to the matrix.
    fn set_cols(&mut self, col1: &NXYZ, col2: &NXYZ, col3: &NXYZ) {
        self.mat[0][0] = col1.x();
        self.mat[1][0] = col1.y();
        self.mat[2][0] = col1.z();
        self.mat[0][1] = col2.x();
        self.mat[1][1] = col2.y();
        self.mat[2][1] = col2.z();
        self.mat[0][2] = col3.x();
        self.mat[1][2] = col3.y();
        self.mat[2][2] = col3.z();
    }

    /// Sets the matrix to represent the cross product with ref_vec.
    fn set_cross(&mut self, ref_vec: &NXYZ) {
        let x = ref_vec.x();
        let y = ref_vec.y();
        let z = ref_vec.z();
        self.mat = [[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]];
    }

    /// Modifies the main diagonal of the matrix.
    fn set_diagonal(&mut self, x1: f64, x2: f64, x3: f64) {
        self.mat[0][0] = x1;
        self.mat[1][1] = x2;
        self.mat[2][2] = x3;
    }

    /// Sets the matrix to represent the dot product with ref_vec.
    fn set_dot(&mut self, ref_vec: &NXYZ) {
        let x = ref_vec.x();
        let y = ref_vec.y();
        let z = ref_vec.z();
        self.mat = [
            [x * x, x * y, x * z],
            [x * y, y * y, y * z],
            [x * z, y * z, z * z],
        ];
    }

    /// Sets the matrix to the identity matrix.
    fn set_identity(&mut self) {
        self.mat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }

    /// Sets the matrix to represent a rotation around axis by ang radians.
    fn set_rotation(&mut self, axis: &NXYZ, ang: f64) -> Result<(), NErrors> {
        let v = axis.normalized()?;
        self.set_cross(&v);
        self.multiply_scalar(ang.sin());
        let mut temp = NMat::new();
        temp.set_scale(1.0);
        self.add(&temp);
        let a = v.x();
        let b = v.y();
        let c = v.z();
        temp.set_row(1, &NXYZ::new(-c * c - b * b, a * b, a * c))?;
        temp.set_row(2, &NXYZ::new(a * b, -a * a - c * c, b * c))?;
        temp.set_row(3, &NXYZ::new(a * c, b * c, -a * a - b * b))?;
        temp.multiply_scalar(1.0 - ang.cos());
        self.add(&temp);
        Ok(())
    }

    /// Assigns the coordinates of value to the row of index row.
    fn set_row(&mut self, row: i32, value: &NXYZ) -> Result<(), NErrors> {
        if row < 1 || row > 3 {
            return Err(NErrors::OutOfRange);
        }
        let row_idx = (row - 1) as usize;
        self.mat[row_idx] = [value.x(), value.y(), value.z()];
        Ok(())
    }

    /// Assigns the three rows to the matrix.
    fn set_rows(&mut self, row1: &NXYZ, row2: &NXYZ, row3: &NXYZ) {
        self.mat = [
            [row1.x(), row1.y(), row1.z()],
            [row2.x(), row2.y(), row2.z()],
            [row3.x(), row3.y(), row3.z()],
        ];
    }

    /// Sets the matrix to represent a scaling transformation.
    fn set_scale(&mut self, s: f64) {
        self.mat = [[s, 0.0, 0.0], [0.0, s, 0.0], [0.0, 0.0, s]];
    }

    /// Assigns value to the coefficient at (row, col).
    fn set_value(&mut self, row: i32, col: i32, value: f64) -> Result<(), NErrors> {
        if row < 1 || row > 3 || col < 1 || col > 3 {
            return Err(NErrors::OutOfRange);
        }
        self.mat[(row - 1) as usize][(col - 1) as usize] = value;
        Ok(())
    }

    /// Returns the column of index col.
    fn column(&self, col: i32) -> Result<NXYZ, NErrors> {
        if col < 1 || col > 3 {
            return Err(NErrors::OutOfRange);
        }
        let col_idx = (col - 1) as usize;
        Ok(NXYZ::new(
            self.mat[0][col_idx],
            self.mat[1][col_idx],
            self.mat[2][col_idx],
        ))
    }

    /// Computes the determinant of the matrix.
    fn determinant(&self) -> f64 {
        self.mat[0][0] * (self.mat[1][1] * self.mat[2][2] - self.mat[2][1] * self.mat[1][2])
            - self.mat[0][1] * (self.mat[1][0] * self.mat[2][2] - self.mat[2][0] * self.mat[1][2])
            + self.mat[0][2] * (self.mat[1][0] * self.mat[2][1] - self.mat[2][0] * self.mat[1][1])
    }

    /// Returns the main diagonal of the matrix.
    fn diagonal(&self) -> NXYZ {
        NXYZ::new(self.mat[0][0], self.mat[1][1], self.mat[2][2])
    }

    /// Returns the row of index row.
    fn row(&self, row: i32) -> Result<NXYZ, NErrors> {
        if row < 1 || row > 3 {
            return Err(NErrors::OutOfRange);
        }
        let row_idx = (row - 1) as usize;
        Ok(NXYZ::new(
            self.mat[row_idx][0],
            self.mat[row_idx][1],
            self.mat[row_idx][2],
        ))
    }

    /// Returns the coefficient at (row, col).
    fn value(&self, row: i32, col: i32) -> Result<f64, NErrors> {
        if row < 1 || row > 3 || col < 1 || col > 3 {
            return Err(NErrors::OutOfRange);
        }
        Ok(self.mat[(row - 1) as usize][(col - 1) as usize])
    }

    /// Returns a mutable reference to the coefficient at (row, col).
    fn change_value(&mut self, row: i32, col: i32) -> Result<&mut f64, NErrors> {
        if row < 1 || row > 3 || col < 1 || col > 3 {
            return Err(NErrors::OutOfRange);
        }
        Ok(&mut self.mat[(row - 1) as usize][(col - 1) as usize])
    }

    /// Checks if the matrix is singular.
    fn is_singular(&self) -> bool {
        self.determinant().abs() <= NGP::resolution()
    }

    /// Adds another matrix to this one.
    fn add(&mut self, other: &Self) {
        for i in 0..3 {
            for j in 0..3 {
                self.mat[i][j] += other.mat[i][j];
            }
        }
    }

    /// Returns the sum of this matrix and another.
    fn added(&self, other: &Self) -> Self {
        let mut new_mat = self.clone();
        new_mat.add(other);
        new_mat
    }

    /// Divides all coefficients by a scalar.
    fn divide(&mut self, scalar: f64) -> Result<(), NErrors> {
        if scalar.abs() <= NGP::resolution() {
            return Err(NErrors::ConstructionError);
        }
        self.multiply_scalar(1.0 / scalar);
        Ok(())
    }

    /// Returns a matrix with all coefficients divided by a scalar.
    fn divided(&self, scalar: f64) -> Result<Self, NErrors> {
        if scalar.abs() <= NGP::resolution() {
            return Err(NErrors::ConstructionError);
        }
        Ok(self.multiplied_scalar(1.0 / scalar))
    }

    /// Inverts the matrix.
    fn invert(&mut self) -> Result<(), NErrors> {
        let mut new_mat = [[0.0; 3]; 3];
        new_mat[0][0] = self.mat[1][1] * self.mat[2][2] - self.mat[1][2] * self.mat[2][1];
        new_mat[1][0] = -(self.mat[1][0] * self.mat[2][2] - self.mat[2][0] * self.mat[1][2]);
        new_mat[2][0] = self.mat[1][0] * self.mat[2][1] - self.mat[2][0] * self.mat[1][1];
        new_mat[0][1] = -(self.mat[0][1] * self.mat[2][2] - self.mat[2][1] * self.mat[0][2]);
        new_mat[1][1] = self.mat[0][0] * self.mat[2][2] - self.mat[2][0] * self.mat[0][2];
        new_mat[2][1] = -(self.mat[0][0] * self.mat[2][1] - self.mat[2][0] * self.mat[0][1]);
        new_mat[0][2] = self.mat[0][1] * self.mat[1][2] - self.mat[1][1] * self.mat[0][2];
        new_mat[1][2] = -(self.mat[0][0] * self.mat[1][2] - self.mat[1][0] * self.mat[0][2]);
        new_mat[2][2] = self.mat[0][0] * self.mat[1][1] - self.mat[0][1] * self.mat[1][0];
        let det = self.mat[0][0] * new_mat[0][0]
            + self.mat[0][1] * new_mat[1][0]
            + self.mat[0][2] * new_mat[2][0];
        if det.abs() <= NGP::resolution() {
            return Err(NErrors::ConstructionError);
        }
        self.mat = new_mat;
        self.multiply_scalar(1.0 / det);
        Ok(())
    }

    /// Returns the inverted matrix.
    fn inverted(&self) -> Result<Self, NErrors> {
        let mut new_mat = self.clone();
        new_mat.invert()?;
        Ok(new_mat)
    }

    /// Multiplies this matrix by another (self = self * other).
    fn multiply_mat(&mut self, other: &Self) {
        let t = [
            [
                self.mat[0][0] * other.mat[0][0]
                    + self.mat[0][1] * other.mat[1][0]
                    + self.mat[0][2] * other.mat[2][0],
                self.mat[0][0] * other.mat[0][1]
                    + self.mat[0][1] * other.mat[1][1]
                    + self.mat[0][2] * other.mat[2][1],
                self.mat[0][0] * other.mat[0][2]
                    + self.mat[0][1] * other.mat[1][2]
                    + self.mat[0][2] * other.mat[2][2],
            ],
            [
                self.mat[1][0] * other.mat[0][0]
                    + self.mat[1][1] * other.mat[1][0]
                    + self.mat[1][2] * other.mat[2][0],
                self.mat[1][0] * other.mat[0][1]
                    + self.mat[1][1] * other.mat[1][1]
                    + self.mat[1][2] * other.mat[2][1],
                self.mat[1][0] * other.mat[0][2]
                    + self.mat[1][1] * other.mat[1][2]
                    + self.mat[1][2] * other.mat[2][2],
            ],
            [
                self.mat[2][0] * other.mat[0][0]
                    + self.mat[2][1] * other.mat[1][0]
                    + self.mat[2][2] * other.mat[2][0],
                self.mat[2][0] * other.mat[0][1]
                    + self.mat[2][1] * other.mat[1][1]
                    + self.mat[2][2] * other.mat[2][1],
                self.mat[2][0] * other.mat[0][2]
                    + self.mat[2][1] * other.mat[1][2]
                    + self.mat[2][2] * other.mat[2][2],
            ],
        ];
        self.mat = t;
    }

    /// Returns the product of this matrix and another.
    fn multiplied_mat(&self, other: &Self) -> Self {
        let mut new_mat = self.clone();
        new_mat.multiply_mat(other);
        new_mat
    }

    /// Multiplies this matrix by another (self = other * self).
    fn pre_multiply(&mut self, other: &Self) {
        let t = [
            [
                other.mat[0][0] * self.mat[0][0]
                    + other.mat[0][1] * self.mat[1][0]
                    + other.mat[0][2] * self.mat[2][0],
                other.mat[0][0] * self.mat[0][1]
                    + other.mat[0][1] * self.mat[1][1]
                    + other.mat[0][2] * self.mat[2][1],
                other.mat[0][0] * self.mat[0][2]
                    + other.mat[0][1] * self.mat[1][2]
                    + other.mat[0][2] * self.mat[2][2],
            ],
            [
                other.mat[1][0] * self.mat[0][0]
                    + other.mat[1][1] * self.mat[1][0]
                    + other.mat[1][2] * self.mat[2][0],
                other.mat[1][0] * self.mat[0][1]
                    + other.mat[1][1] * self.mat[1][1]
                    + other.mat[1][2] * self.mat[2][1],
                other.mat[1][0] * self.mat[0][2]
                    + other.mat[1][1] * self.mat[1][2]
                    + other.mat[1][2] * self.mat[2][2],
            ],
            [
                other.mat[2][0] * self.mat[0][0]
                    + other.mat[2][1] * self.mat[1][0]
                    + other.mat[2][2] * self.mat[2][0],
                other.mat[2][0] * self.mat[0][1]
                    + other.mat[2][1] * self.mat[1][1]
                    + other.mat[2][2] * self.mat[2][1],
                other.mat[2][0] * self.mat[0][2]
                    + other.mat[2][1] * self.mat[1][2]
                    + other.mat[2][2] * self.mat[2][2],
            ],
        ];
        self.mat = t;
    }

    /// Multiplies all coefficients by a scalar.
    fn multiply_scalar(&mut self, scalar: f64) {
        for i in 0..3 {
            for j in 0..3 {
                self.mat[i][j] *= scalar;
            }
        }
    }

    /// Returns a matrix with all coefficients multiplied by a scalar.
    fn multiplied_scalar(&self, scalar: f64) -> Self {
        let mut new_mat = self.clone();
        new_mat.multiply_scalar(scalar);
        new_mat
    }

    /// Computes the matrix raised to the power n.
    fn power(&mut self, n: i32) -> Result<(), NErrors> {
        if n == 1 {
            // No-op
        } else if n == 0 {
            self.set_identity();
        } else if n == -1 {
            self.invert()?;
        } else {
            let mut n_power = n;
            if n_power < 0 {
                self.invert()?;
                n_power = -n_power;
            }
            n_power -= 1;
            let temp = self.clone();
            loop {
                if n_power % 2 == 1 {
                    self.multiply_mat(&temp);
                }
                if n_power == 1 {
                    break;
                }
                let mut temp2 = temp.clone();
                temp2.multiply_mat(&temp);
                n_power /= 2;
            }
        }
        Ok(())
    }

    /// Returns the matrix raised to the power n.
    fn powered(&self, n: i32) -> Result<Self, NErrors> {
        let mut new_mat = self.clone();
        new_mat.power(n)?;
        Ok(new_mat)
    }

    /// Subtracts another matrix from this one.
    fn subtract(&mut self, other: &Self) {
        for i in 0..3 {
            for j in 0..3 {
                self.mat[i][j] -= other.mat[i][j];
            }
        }
    }

    /// Returns the difference of this matrix and another.
    fn subtracted(&self, other: &Self) -> Self {
        let mut new_mat = self.clone();
        new_mat.subtract(other);
        new_mat
    }

    /// Transposes the matrix.
    fn transpose(&mut self) {
        let t = [
            [self.mat[0][0], self.mat[1][0], self.mat[2][0]],
            [self.mat[0][1], self.mat[1][1], self.mat[2][1]],
            [self.mat[0][2], self.mat[1][2], self.mat[2][2]],
        ];
        self.mat = t;
    }

    /// Returns the transposed matrix.
    fn transposed(&self) -> Self {
        let mut new_mat = self.clone();
        new_mat.transpose();
        new_mat
    }

    /// Dumps the matrix as JSON.
    fn dump_json(&self, out: &mut dyn Write, _depth: i32) {
        writeln!(out, "{{\n  \"type\": \"NMat\",").unwrap();
        writeln!(
            out,
            "  \"elements\": [{}, {}, {}, {}, {}, {}, {}, {}, {}]",
            self.mat[0][0],
            self.mat[0][1],
            self.mat[0][2],
            self.mat[1][0],
            self.mat[1][1],
            self.mat[1][2],
            self.mat[2][0],
            self.mat[2][1],
            self.mat[2][2]
        )
        .unwrap();
        writeln!(out, "}}").unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn identity() -> NMat {
        let mut m = NMat::new();
        m.set_identity();
        m
    }

    #[test]
    fn test_new() {
        let m = NMat::new();
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(m.mat[i][j], 0.0);
            }
        }
    }

    #[test]
    fn test_new_with_elements() {
        let m = NMat::new_with_elements(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
        assert_eq!(m.mat, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
    }

    #[test]
    fn test_new_with_columns() {
        let m = NMat::new_with_columns(
            &NXYZ::new(1.0, 4.0, 7.0),
            &NXYZ::new(2.0, 5.0, 8.0),
            &NXYZ::new(3.0, 6.0, 9.0),
        );
        assert_eq!(m.mat, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
    }

    #[test]
    fn test_set_col() {
        let mut m = NMat::new();
        m.set_col(1, &NXYZ::new(1.0, 2.0, 3.0)).unwrap();
        assert_eq!(m.column(1).unwrap(), NXYZ::new(1.0, 2.0, 3.0));
        assert!(matches!(
            m.set_col(4, &NXYZ::new(0.0, 0.0, 0.0)),
            Err(NErrors::OutOfRange)
        ));
    }

    #[test]
    fn test_set_cross() {
        let mut m = NMat::new();
        m.set_cross(&NXYZ::new(1.0, 0.0, 0.0));
        assert_eq!(m.mat, [[0.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]]);
    }

    #[test]
    fn test_set_dot() {
        let mut m = NMat::new();
        m.set_dot(&NXYZ::new(1.0, 2.0, 3.0));
        assert_eq!(m.mat, [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [3.0, 6.0, 9.0]]);
    }

    #[test]
    fn test_set_rotation() {
        let mut m = NMat::new();
        m.set_rotation(&NXYZ::new(0.0, 0.0, 1.0), std::f64::consts::PI / 2.0)
            .unwrap();
        let expected = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        for i in 0..3 {
            for j in 0..3 {
                assert!((m.mat[i][j] - expected[i][j]).abs() < 1e-9);
            }
        }
        assert!(matches!(
            m.set_rotation(&NXYZ::new(0.0, 0.0, 0.0), 0.0),
            Err(NErrors::ConstructionError)
        ));
    }

    #[test]
    fn test_invert() {
        let mut m = identity();
        m.invert().unwrap();
        assert_eq!(m, identity());
        let mut m = NMat::new_with_elements(1.0, 2.0, 3.0, 0.0, 1.0, 4.0, 5.0, 6.0, 0.0);
        m.invert().unwrap();
        let expected = [[-24.0, 18.0, 5.0], [20.0, -15.0, -4.0], [-5.0, 4.0, 1.0]];
        for i in 0..3 {
            for j in 0..3 {
                assert!((m.mat[i][j] - expected[i][j]).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn test_multiply_mat() {
        let mut m1 = identity();
        let m2 = NMat::new_with_elements(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
        m1.multiply_mat(&m2);
        assert_eq!(m1, m2);
    }

    #[test]
    fn test_power() {
        let mut m = identity();
        m.power(0).unwrap();
        assert_eq!(m, identity());
        let mut m = NMat::new_with_elements(1.0, 2.0, 3.0, 0.0, 1.0, 4.0, 5.0, 6.0, 0.0);
        m.power(2).unwrap();
        let expected = [[16.0, 20.0, 11.0], [20.0, 17.0, 4.0], [5.0, 12.0, 29.0]];
        for i in 0..3 {
            for j in 0..3 {
                assert!((m.mat[i][j] - expected[i][j]).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn test_transpose() {
        let mut m = NMat::new_with_elements(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
        m.transpose();
        assert_eq!(m.mat, [[1.0, 4.0, 7.0], [2.0, 5.0, 8.0], [3.0, 6.0, 9.0]]);
    }

    #[test]
    fn test_dump_json() {
        let m = identity();
        let mut output = Vec::new();
        m.dump_json(&mut output, -1);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NMat\""));
        assert!(json.contains("1, 0, 0, 0, 1, 0, 0, 0, 1"));
    }
}
