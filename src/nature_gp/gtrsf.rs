use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::{
    gp::{NAx1, NAx2, NGP, NMat, NTrsf, NTrsfForm, NXYZ},
    nature_errors::NErrors,
};

// Placeholder for NMat4 (4x4 matrix), assuming it's not provided
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NMat4 {
    data: [[f64; 4]; 4],
}

impl NMat4 {
    pub fn new() -> Self {
        NMat4 {
            data: [[0.0; 4]; 4],
        }
    }

    pub fn init_identity(&mut self) {
        self.data = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
    }

    pub fn set_value(&mut self, row: usize, col: usize, value: f64) {
        self.data[row][col] = value;
    }

    pub fn get_value(&self, row: usize, col: usize) -> f64 {
        self.data[row][col]
    }
}

// Trait to define the behavior of a general transformation in 3D space
pub trait GTrsf {
    fn new() -> Self;
    fn from_trsf(t: &NTrsf) -> Self;
    fn from_matrix_vector(m: &NMat, v: &NXYZ) -> Self;
    fn set_affinity_ax1(&mut self, a1: &NAx1, ratio: f64);
    fn set_affinity_ax2(&mut self, a2: &NAx2, ratio: f64);
    fn set_value(&mut self, row: usize, col: usize, value: f64) -> Result<(), NErrors>;
    fn set_vectorial_part(&mut self, matrix: &NMat);
    fn set_translation_part(&mut self, coord: &NXYZ);
    fn set_trsf(&mut self, t: &NTrsf);
    fn is_negative(&self) -> bool;
    fn is_singular(&self) -> bool;
    fn form(&self) -> NTrsfForm;
    fn set_form(&mut self);
    fn translation_part(&self) -> &NXYZ;
    fn vectorial_part(&self) -> &NMat;
    fn value(&self, row: usize, col: usize) -> Result<f64, NErrors>;
    fn invert(&mut self) -> Result<(), NErrors>;
    fn inverted(&self) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn multiply(&mut self, t: &Self);
    fn multiplied(&self, t: &Self) -> Self
    where
        Self: Sized;
    fn pre_multiply(&mut self, t: &Self);
    fn power(&mut self, n: i32) -> Result<(), NErrors>;
    fn powered(&self, n: i32) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn transforms_xyz(&self, coord: &mut NXYZ);
    fn transforms_coords(&self, x: &mut f64, y: &mut f64, z: &mut f64);
    fn to_trsf(&self) -> Result<NTrsf, NErrors>;
    fn get_mat4(&self) -> NMat4;
    fn set_mat4(&mut self, mat: &NMat4);
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
}

// Struct representing a general transformation in 3D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NGTrsf {
    matrix: NMat,
    loc: NXYZ,
    shape: NTrsfForm,
    scale: f64,
}

impl GTrsf for NGTrsf {
    /// Returns the identity transformation.
    fn new() -> Self {
        let mut matrix = NMat::new();
        matrix.set_scale(1.0);
        NGTrsf {
            matrix,
            loc: NXYZ::new(0.0, 0.0, 0.0),
            shape: NTrsfForm::Identity,
            scale: 1.0,
        }
    }

    /// Converts an `NTrsf` into a general transformation.
    fn from_trsf(t: &NTrsf) -> Self {
        NGTrsf {
            shape: t.form(),
            matrix: t.matrix(),
            loc: t.translation_part(),
            scale: t.scale_factor(),
        }
    }

    /// Creates a transformation from a matrix and vector.
    fn from_matrix_vector(m: &NMat, v: &NXYZ) -> Self {
        NGTrsf {
            matrix: m.clone(),
            loc: v.clone(),
            shape: NTrsfForm::Other,
            scale: 0.0,
        }
    }

    /// Sets this transformation as an affinity with respect to an axis.
    fn set_affinity_ax1(&mut self, a1: &NAx1, ratio: f64) {
        self.shape = NTrsfForm::Other;
        self.scale = 0.0;
        self.matrix.set_dot(&a1.direction().xyz());
        self.matrix.multiply(1.0 - ratio);
        self.matrix.set_diagonal(
            self.matrix.value(1, 1) + ratio,
            self.matrix.value(2, 2) + ratio,
            self.matrix.value(3, 3) + ratio,
        );
        self.loc = a1.location().xyz();
        self.loc.reverse();
        let loc_copy = self.loc.clone();
        self.loc.multiply_matrix(&self.matrix);
        self.loc.add(&loc_copy);
    }

    /// Sets this transformation as an affinity with respect to a plane.
    fn set_affinity_ax2(&mut self, a2: &NAx2, ratio: f64) {
        self.shape = NTrsfForm::Other;
        self.scale = 0.0;
        self.matrix.set_dot(&a2.direction().xyz());
        self.matrix.multiply(ratio - 1.0);
        self.loc = a2.location().xyz();
        self.loc.reverse();
        let loc_copy = self.loc.clone();
        self.loc.multiply_matrix(&self.matrix);
        self.matrix.set_diagonal(
            self.matrix.value(1, 1) + 1.0,
            self.matrix.value(2, 2) + 1.0,
            self.matrix.value(3, 3) + 1.0,
        );
        self.loc.add(&loc_copy);
    }

    /// Sets a specific coefficient in the transformation matrix.
    fn set_value(&mut self, row: usize, col: usize, value: f64) -> Result<(), NErrors> {
        if row < 1 || row > 3 || col < 1 || col > 4 {
            return Err(NErrors::OutOfRange);
        }
        if col == 4 {
            self.loc.set_coord(row, value);
            if self.shape == NTrsfForm::Identity {
                self.shape = NTrsfForm::Translation;
            }
        } else {
            if self.shape != NTrsfForm::Other && self.scale != 1.0 {
                self.matrix.multiply(self.scale);
            }
            self.matrix.set_value(row, col, value);
            self.shape = NTrsfForm::Other;
            self.scale = 0.0;
        }
        Ok(())
    }

    /// Sets the vectorial part of the transformation.
    fn set_vectorial_part(&mut self, matrix: &NMat) {
        self.matrix = matrix.clone();
        self.shape = NTrsfForm::Other;
        self.scale = 0.0;
    }

    /// Sets the translation part of the transformation.
    fn set_translation_part(&mut self, coord: &NXYZ) {
        self.loc = coord.clone();
        match self.form() {
            NTrsfForm::CompoundTrsf | NTrsfForm::Other | NTrsfForm::Translation => {}
            NTrsfForm::Identity => self.shape = NTrsfForm::Translation,
            _ => self.shape = NTrsfForm::CompoundTrsf,
        }
    }

    /// Assigns the vectorial and translation parts from an `NTrsf`.
    fn set_trsf(&mut self, t: &NTrsf) {
        self.shape = t.form();
        self.matrix = t.matrix();
        self.loc = t.translation_part();
        self.scale = t.scale_factor();
    }

    /// Returns true if the determinant of the vectorial part is negative.
    fn is_negative(&self) -> bool {
        self.matrix.determinant() < 0.0
    }

    /// Returns true if the transformation is singular (non-invertible).
    fn is_singular(&self) -> bool {
        self.matrix.is_singular()
    }

    /// Returns the nature of the transformation.
    fn form(&self) -> NTrsfForm {
        self.shape
    }

    /// Verifies and sets the shape of the transformation.
    fn set_form(&mut self) {
        let tol = 1e-12;
        let mut m = self.matrix.clone();
        let s = m.determinant();

        if s.abs() < NGP::resolution() {
            panic!("Null determinant in set_form");
        }

        let s = if s > 0.0 {
            s.powf(1.0 / 3.0)
        } else {
            (-s).powf(1.0 / 3.0)
        };
        m.divide(s);

        let mut tm = m.clone();
        tm.transpose();
        tm.multiply(&m);
        let mut identity = NMat::new();
        identity.set_identity();
        tm.subtract(&identity);

        if self.shape == NTrsfForm::Other {
            self.shape = NTrsfForm::CompoundTrsf;
        }

        for i in 1..=3 {
            for j in 1..=3 {
                if tm.value(i, j).abs() > tol {
                    self.shape = NTrsfForm::Other;
                    return;
                }
            }
        }
    }

    /// Returns the translation part of the transformation.
    fn translation_part(&self) -> &NXYZ {
        &self.loc
    }

    /// Returns the vectorial part of the transformation.
    fn vectorial_part(&self) -> &NMat {
        &self.matrix
    }

    /// Returns a coefficient of the transformation matrix.
    fn value(&self, row: usize, col: usize) -> Result<f64, NErrors> {
        if row < 1 || row > 3 || col < 1 || col > 4 {
            return Err(NErrors::OutOfRange);
        }
        if col == 4 {
            Ok(self.loc.coord(row))
        } else if self.shape == NTrsfForm::Other {
            Ok(self.matrix.value(row, col))
        } else {
            Ok(self.scale * self.matrix.value(row, col))
        }
    }

    /// Inverts the transformation.
    fn invert(&mut self) -> Result<(), NErrors> {
        if self.shape == NTrsfForm::Other {
            self.matrix.invert()?;
            let loc_copy = self.loc.clone();
            self.loc.multiply_matrix(&self.matrix);
            self.loc.reverse();
            Ok(())
        } else {
            let mut t = self.to_trsf()?;
            t.invert()?;
            self.set_trsf(&t);
            Ok(())
        }
    }

    /// Returns the inverted transformation.
    fn inverted(&self) -> Result<Self, NErrors> {
        let mut result = self.clone();
        result.invert()?;
        Ok(result)
    }

    /// Composes this transformation with another (self = self * t).
    fn multiply(&mut self, t: &Self) {
        if self.form() == NTrsfForm::Other || t.form() == NTrsfForm::Other {
            self.shape = NTrsfForm::Other;
            let mut t_loc = t.loc.clone();
            t_loc.multiply_matrix(&self.matrix);
            self.loc.add(&t_loc);
            self.matrix.multiply(&t.matrix);
        } else {
            let mut t1 = self.to_trsf().unwrap();
            let t2 = t.to_trsf().unwrap();
            t1.multiply(&t2);
            self.matrix = t1.matrix();
            self.loc = t1.translation_part();
            self.scale = t1.scale_factor();
            self.shape = t1.form();
        }
    }

    /// Returns the composed transformation (self * t).
    fn multiplied(&self, t: &Self) -> Self {
        let mut result = self.clone();
        result.multiply(t);
        result
    }

    /// Composes another transformation with this (self = t * self).
    fn pre_multiply(&mut self, t: &Self) {
        if self.form() == NTrsfForm::Other || t.form() == NTrsfForm::Other {
            self.shape = NTrsfForm::Other;
            self.loc.multiply_matrix(&t.matrix);
            self.loc.add(&t.loc);
            self.matrix.pre_multiply(&t.matrix);
        } else {
            let mut t1 = self.to_trsf().unwrap();
            let t2 = t.to_trsf().unwrap();
            t1.pre_multiply(&t2);
            self.matrix = t1.matrix();
            self.loc = t1.translation_part();
            self.scale = t1.scale_factor();
            self.shape = t1.form();
        }
    }

    /// Computes the transformation raised to a power.
    fn power(&mut self, n: i32) -> Result<(), NErrors> {
        if n == 0 {
            self.scale = 1.0;
            self.shape = NTrsfForm::Identity;
            self.matrix.set_identity();
            self.loc = NXYZ::new(0.0, 0.0, 0.0);
        } else if n == 1 {
            // No change
        } else if n == -1 {
            self.invert()?;
        } else if self.shape == NTrsfForm::Other {
            let mut n_power = n.abs();
            if n < 0 {
                self.invert()?;
            }
            n_power -= 1;
            let temp_loc = self.loc.clone();
            let temp_matrix = self.matrix.clone();
            loop {
                if n_power % 2 == 1 {
                    let mut loc_copy = temp_loc.clone();
                    loc_copy.multiply_matrix(&self.matrix);
                    self.loc.add(&loc_copy);
                    self.matrix.multiply(&temp_matrix);
                }
                if n_power == 1 {
                    break;
                }
                let mut loc_copy = temp_loc.clone();
                loc_copy.multiply_matrix(&temp_matrix);
                let mut temp_loc = temp_loc.clone();
                temp_loc.add(&loc_copy);
                let mut temp_matrix = temp_matrix.clone();
                temp_matrix.multiply(&temp_matrix);
                n_power /= 2;
            }
        } else {
            let mut t = self.to_trsf()?;
            t.power(n)?;
            self.set_trsf(&t);
        }
        Ok(())
    }

    /// Returns the transformation raised to a power.
    fn powered(&self, n: i32) -> Result<Self, NErrors> {
        let mut result = self.clone();
        result.power(n)?;
        Ok(result)
    }

    /// Transforms a coordinate triplet.
    fn transforms_xyz(&self, coord: &mut NXYZ) {
        coord.multiply_matrix(&self.matrix);
        if self.shape != NTrsfForm::Other && self.scale != 1.0 {
            coord.multiply(self.scale);
        }
        coord.add(&self.loc);
    }

    /// Transforms individual coordinates.
    fn transforms_coords(&self, x: &mut f64, y: &mut f64, z: &mut f64) {
        let mut triplet = NXYZ::new(*x, *y, *z);
        triplet.multiply_matrix(&self.matrix);
        if self.shape != NTrsfForm::Other && self.scale != 1.0 {
            triplet.multiply(self.scale);
        }
        triplet.add(&self.loc);
        *x = triplet.x();
        *y = triplet.y();
        *z = triplet.z();
    }

    /// Converts to an `NTrsf` if possible.
    fn to_trsf(&self) -> Result<NTrsf, NErrors> {
        if self.form() == NTrsfForm::Other {
            return Err(NErrors::ConstructionError);
        }
        let mut t = NTrsf::new();
        t.set_form(self.shape);
        t.set_scale(self.scale);
        t.set_matrix(&self.matrix);
        t.set_translation_part(&self.loc);
        Ok(t)
    }

    /// Converts the transformation to a 4x4 matrix.
    fn get_mat4(&self) -> NMat4 {
        let mut mat = NMat4::new();
        if self.shape == NTrsfForm::Identity {
            mat.init_identity();
            return mat;
        }
        mat.set_value(0, 0, self.value(1, 1).unwrap());
        mat.set_value(0, 1, self.value(1, 2).unwrap());
        mat.set_value(0, 2, self.value(1, 3).unwrap());
        mat.set_value(0, 3, self.value(1, 4).unwrap());
        mat.set_value(1, 0, self.value(2, 1).unwrap());
        mat.set_value(1, 1, self.value(2, 2).unwrap());
        mat.set_value(1, 2, self.value(2, 3).unwrap());
        mat.set_value(1, 3, self.value(2, 4).unwrap());
        mat.set_value(2, 0, self.value(3, 1).unwrap());
        mat.set_value(2, 1, self.value(3, 2).unwrap());
        mat.set_value(2, 2, self.value(3, 3).unwrap());
        mat.set_value(2, 3, self.value(3, 4).unwrap());
        mat.set_value(3, 0, 0.0);
        mat.set_value(3, 1, 0.0);
        mat.set_value(3, 2, 0.0);
        mat.set_value(3, 3, 1.0);
        mat
    }

    /// Sets the transformation from a 4x4 matrix.
    fn set_mat4(&mut self, mat: &NMat4) {
        self.shape = NTrsfForm::Other;
        self.scale = 0.0;
        self.matrix.set_value(1, 1, mat.get_value(0, 0));
        self.matrix.set_value(1, 2, mat.get_value(0, 1));
        self.matrix.set_value(1, 3, mat.get_value(0, 2));
        self.matrix.set_value(2, 1, mat.get_value(1, 0));
        self.matrix.set_value(2, 2, mat.get_value(1, 1));
        self.matrix.set_value(2, 3, mat.get_value(1, 2));
        self.matrix.set_value(3, 1, mat.get_value(2, 0));
        self.matrix.set_value(3, 2, mat.get_value(2, 1));
        self.matrix.set_value(3, 3, mat.get_value(2, 2));
        self.loc.set_coord(
            mat.get_value(0, 3),
            mat.get_value(1, 3),
            mat.get_value(2, 3),
        );
    }

    /// Dumps the transformation as JSON.
    fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        writeln!(out, "{{\n  \"type\": \"NGTrsf\",").unwrap();
        writeln!(out, "  \"matrix\": {{").unwrap();
        for i in 1..=3 {
            writeln!(
                out,
                "    \"row{}\": [{}, {}, {}],",
                i,
                self.matrix.value(i, 1),
                self.matrix.value(i, 2),
                self.matrix.value(i, 3)
            )
            .unwrap();
        }
        writeln!(out, "  }},").unwrap();
        writeln!(
            out,
            "  \"loc\": [{}, {}, {}],",
            self.loc.x(),
            self.loc.y(),
            self.loc.z()
        )
        .unwrap();
        writeln!(out, "  \"shape\": \"{:?}\",", self.shape).unwrap();
        writeln!(out, "  \"scale\": {}", self.scale).unwrap();
        writeln!(out, "}}").unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gtrsf() -> NGTrsf {
        NGTrsf::new()
    }

    fn ax1() -> NAx1 {
        NAx1::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
        )
    }

    fn ax2() -> NAx2 {
        NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn test_new() {
        let g = gtrsf();
        assert_eq!(g.form(), NTrsfForm::Identity);
        assert_eq!(g.scale, 1.0);
        assert_eq!(g.translation_part(), &NXYZ::new(0.0, 0.0, 0.0));
        let m = g.vectorial_part();
        assert_eq!(m.value(1, 1), 1.0);
    }

    #[test]
    fn test_from_trsf() {
        let t = NTrsf::new(); // Placeholder, assumes NTrsf::new() exists
        let g = NGTrsf::from_trsf(&t);
        assert_eq!(g.form(), t.form());
        assert_eq!(g.scale, t.scale_factor());
    }

    #[test]
    fn test_from_matrix_vector() {
        let m = NMat::new();
        let v = NXYZ::new(1.0, 2.0, 3.0);
        let g = NGTrsf::from_matrix_vector(&m, &v);
        assert_eq!(g.form(), NTrsfForm::Other);
        assert_eq!(g.scale, 0.0);
        assert_eq!(g.translation_part(), &v);
    }

    #[test]
    fn test_set_affinity() {
        let mut g = gtrsf();
        g.set_affinity_ax1(&ax1(), 2.0);
        assert_eq!(g.form(), NTrsfForm::Other);
        assert_eq!(g.scale, 0.0);

        let mut g = gtrsf();
        g.set_affinity_ax2(&ax2(), 2.0);
        assert_eq!(g.form(), NTrsfForm::Other);
        assert_eq!(g.scale, 0.0);
    }

    #[test]
    fn test_set_value() {
        let mut g = gtrsf();
        g.set_value(1, 1, 2.0).unwrap();
        assert_eq!(g.value(1, 1).unwrap(), 2.0);
        assert_eq!(g.form(), NTrsfForm::Other);
        g.set_value(1, 4, 3.0).unwrap();
        assert_eq!(g.value(1, 4).unwrap(), 3.0);
        assert!(matches!(g.set_value(4, 1, 1.0), Err(NErrors::OutOfRange)));
    }

    #[test]
    fn test_set_translation_part() {
        let mut g = gtrsf();
        let coord = NXYZ::new(1.0, 2.0, 3.0);
        g.set_translation_part(&coord);
        assert_eq!(g.translation_part(), &coord);
        assert_eq!(g.form(), NTrsfForm::Translation);
    }

    #[test]
    fn test_is_negative() {
        let mut g = gtrsf();
        assert!(!g.is_negative());
        let mut m = NMat::new();
        m.set_value(1, 1, -1.0);
        g.set_vectorial_part(&m);
        assert!(g.is_negative());
    }

    #[test]
    fn test_form_and_set_form() {
        let mut g = gtrsf();
        assert_eq!(g.form(), NTrsfForm::Identity);
        g.set_form();
        assert_eq!(g.form(), NTrsfForm::CompoundTrsf);
    }

    #[test]
    fn test_transforms() {
        let g = gtrsf();
        let mut xyz = NXYZ::new(1.0, 2.0, 3.0);
        g.transforms_xyz(&mut xyz);
        assert_eq!(xyz, NXYZ::new(1.0, 2.0, 3.0));

        let mut x = 1.0;
        let mut y = 2.0;
        let mut z = 3.0;
        g.transforms_coords(&mut x, &mut y, &mut z);
        assert_eq!(x, 1.0);
        assert_eq!(y, 2.0);
        assert_eq!(z, 3.0);
    }

    #[test]
    fn test_mat4() {
        let g = gtrsf();
        let mat = g.get_mat4();
        assert_eq!(mat.get_value(0, 0), 1.0);
        assert_eq!(mat.get_value(3, 3), 1.0);

        let mut g2 = gtrsf();
        g2.set_mat4(&mat);
        assert_eq!(g2.form(), NTrsfForm::Other);
        assert_eq!(g2.scale, 0.0);
    }

    #[test]
    fn test_dump_json() {
        let g = gtrsf();
        let mut output = Vec::new();
        g.dump_json(&mut output, -1);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NGTrsf\""));
        assert!(json.contains("\"shape\": \"Identity\""));
    }
}
