use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::{
    gp::{NAx2d, NMat2d, NTrsf2d, NTrsfForm, NXY},
    nature_errors::NErrors,
};

// Assume NPrecision provides angular precision
pub mod NPrecision {
    pub fn angular() -> f64 {
        1e-12 // Typical value for Precision::Angular()
    }
}

// Trait to define the behavior of a general transformation in 2D space
pub trait GTrsf2d {
    fn new() -> Self;
    fn from_trsf2d(t: &NTrsf2d) -> Self;
    fn from_matrix_vector(m: &NMat2d, v: &NXY) -> Self;
    fn set_affinity(&mut self, a: &NAx2d, ratio: f64);
    fn set_value(&mut self, row: usize, col: usize, value: f64) -> Result<(), NErrors>;
    fn set_vectorial_part(&mut self, matrix: &NMat2d);
    fn set_translation_part(&mut self, coord: &NXY);
    fn set_trsf2d(&mut self, t: &NTrsf2d);
    fn is_negative(&self) -> bool;
    fn is_singular(&self) -> bool;
    fn form(&self) -> NTrsfForm;
    fn translation_part(&self) -> &NXY;
    fn vectorial_part(&self) -> &NMat2d;
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
    fn transforms_xy(&self, coord: &mut NXY);
    fn transforms_coords(&self, x: &mut f64, y: &mut f64);
    fn to_trsf2d(&self) -> Result<NTrsf2d, NErrors>;
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
}

// Struct representing a general transformation in 2D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NGTrsf2d {
    matrix: NMat2d,
    loc: NXY,
    shape: NTrsfForm,
    scale: f64,
}

impl GTrsf2d for NGTrsf2d {
    /// Returns the identity transformation.
    fn new() -> Self {
        let mut matrix = NMat2d::new();
        matrix.set_scale(1.0);
        NGTrsf2d {
            matrix,
            loc: NXY::new(0.0, 0.0),
            shape: NTrsfForm::Identity,
            scale: 1.0,
        }
    }

    /// Converts an `NTrsf2d` into a general transformation.
    fn from_trsf2d(t: &NTrsf2d) -> Self {
        NGTrsf2d {
            shape: t.form(),
            matrix: t.matrix(),
            loc: t.translation_part(),
            scale: t.scale_factor(),
        }
    }

    /// Creates a transformation from a matrix and vector.
    fn from_matrix_vector(m: &NMat2d, v: &NXY) -> Self {
        NGTrsf2d {
            matrix: m.clone(),
            loc: v.clone(),
            shape: NTrsfForm::Other,
            scale: 0.0,
        }
    }

    /// Sets this transformation as an affinity with respect to an axis.
    fn set_affinity(&mut self, a: &NAx2d, ratio: f64) {
        self.shape = NTrsfForm::Other;
        self.scale = 0.0;
        let (a_x, a_y) = a.direction().coords();
        self.matrix
            .set_value(1, 1, (1.0 - ratio) * a_x * a_x + ratio);
        self.matrix
            .set_value(2, 2, (1.0 - ratio) * a_y * a_y + ratio);
        self.matrix.set_value(1, 2, (1.0 - ratio) * a_x * a_y);
        self.matrix.set_value(2, 1, self.matrix.value(1, 2));
        self.loc = a.location().xy();
        self.loc.reverse();
        let loc_copy = self.loc.clone();
        self.loc.multiply_matrix(&self.matrix);
        self.loc.add(&loc_copy);
    }

    /// Sets a specific coefficient in the transformation matrix.
    fn set_value(&mut self, row: usize, col: usize, value: f64) -> Result<(), NErrors> {
        if row < 1 || row > 2 || col < 1 || col > 3 {
            return Err(NErrors::OutOfRange);
        }
        if col == 3 {
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
    fn set_vectorial_part(&mut self, matrix: &NMat2d) {
        self.matrix = matrix.clone();
        self.shape = NTrsfForm::Other;
        self.scale = 0.0;
    }

    /// Sets the translation part of the transformation.
    fn set_translation_part(&mut self, coord: &NXY) {
        self.loc = coord.clone();
        match self.form() {
            NTrsfForm::CompoundTrsf | NTrsfForm::Other | NTrsfForm::Translation => {}
            NTrsfForm::Identity => self.shape = NTrsfForm::Translation,
            _ => self.shape = NTrsfForm::CompoundTrsf,
        }
    }

    /// Assigns the vectorial and translation parts from an `NTrsf2d`.
    fn set_trsf2d(&mut self, t: &NTrsf2d) {
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

    /// Returns the translation part of the transformation.
    fn translation_part(&self) -> &NXY {
        &self.loc
    }

    /// Returns the vectorial part of the transformation.
    fn vectorial_part(&self) -> &NMat2d {
        &self.matrix
    }

    /// Returns a coefficient of the transformation matrix.
    fn value(&self, row: usize, col: usize) -> Result<f64, NErrors> {
        if row < 1 || row > 2 || col < 1 || col > 3 {
            return Err(NErrors::OutOfRange);
        }
        if col == 3 {
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
            let mut t = self.to_trsf2d()?;
            t.invert()?;
            self.set_trsf2d(&t);
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
            let mut t1 = self.to_trsf2d().unwrap();
            let t2 = t.to_trsf2d().unwrap();
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
            let mut t1 = self.to_trsf2d().unwrap();
            let t2 = t.to_trsf2d().unwrap();
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
            self.loc = NXY::new(0.0, 0.0);
        } else if n == 1 {
            // No change
        } else if n == -1 {
            self.invert()?;
        } else {
            if n < 0 {
                self.invert()?;
            }
            let mut n_power = n.abs();
            if self.shape == NTrsfForm::Other {
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
                let mut t = self.to_trsf2d()?;
                t.power(n)?;
                self.set_trsf2d(&t);
            }
        }
        Ok(())
    }

    /// Returns the transformation raised to a power.
    fn powered(&self, n: i32) -> Result<Self, NErrors> {
        let mut result = self.clone();
        result.power(n)?;
        Ok(result)
    }

    /// Transforms a coordinate pair.
    fn transforms_xy(&self, coord: &mut NXY) {
        coord.multiply_matrix(&self.matrix);
        if self.shape != NTrsfForm::Other && self.scale != 1.0 {
            coord.multiply(self.scale);
        }
        coord.add(&self.loc);
    }

    /// Transforms individual coordinates.
    fn transforms_coords(&self, x: &mut f64, y: &mut f64) {
        let mut pair = NXY::new(*x, *y);
        pair.multiply_matrix(&self.matrix);
        if self.shape != NTrsfForm::Other && self.scale != 1.0 {
            pair.multiply(self.scale);
        }
        pair.add(&self.loc);
        *x = pair.x();
        *y = pair.y();
    }

    /// Converts to an `NTrsf2d` if possible.
    fn to_trsf2d(&self) -> Result<NTrsf2d, NErrors> {
        if self.form() == NTrsfForm::Other {
            return Err(NErrors::ConstructionError);
        }

        let tolerance = NPrecision::angular();
        let tolerance2 = 2.0 * tolerance;

        let value = self.matrix.value(1, 1) * self.matrix.value(1, 1)
            + self.matrix.value(2, 1) * self.matrix.value(2, 1);
        if (value - 1.0).abs() > tolerance2 {
            return Err(NErrors::ConstructionError);
        }

        let value = self.matrix.value(1, 2) * self.matrix.value(1, 2)
            + self.matrix.value(2, 2) * self.matrix.value(2, 2);
        if (value - 1.0).abs() > tolerance2 {
            return Err(NErrors::ConstructionError);
        }

        let value = self.matrix.value(1, 1) * self.matrix.value(1, 2)
            + self.matrix.value(2, 1) * self.matrix.value(2, 2);
        if value.abs() > tolerance {
            return Err(NErrors::ConstructionError);
        }

        let mut t = NTrsf2d::new();
        t.set_matrix(&self.matrix);
        t.set_form(self.shape);
        t.set_scale(self.scale);
        t.set_translation_part(&self.loc);
        Ok(t)
    }

    /// Dumps the transformation as JSON.
    fn dump_json(&self, out: &mut dyn Write, _depth: i32) {
        writeln!(out, "{{\n  \"type\": \"NGTrsf2d\",").unwrap();
        writeln!(out, "  \"matrix\": {{").unwrap();
        for i in 1..=2 {
            writeln!(
                out,
                "    \"row{}\": [{}, {}],",
                i,
                self.matrix.value(i, 1),
                self.matrix.value(i, 2)
            )
            .unwrap();
        }
        writeln!(out, "  }},").unwrap();
        writeln!(out, "  \"loc\": [{}, {}],", self.loc.x(), self.loc.y()).unwrap();
        writeln!(out, "  \"shape\": \"{:?}\",", self.shape).unwrap();
        writeln!(out, "  \"scale\": {}", self.scale).unwrap();
        writeln!(out, "}}").unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gtrsf2d() -> NGTrsf2d {
        NGTrsf2d::new()
    }

    fn ax2d() -> NAx2d {
        NAx2d::new(&NPoint2d::new(0.0, 0.0), &NDir2d::new(1.0, 0.0).unwrap())
    }

    #[test]
    fn test_new() {
        let g = gtrsf2d();
        assert_eq!(g.form(), NTrsfForm::Identity);
        assert_eq!(g.scale, 1.0);
        assert_eq!(g.translation_part(), &NXY::new(0.0, 0.0));
        let m = g.vectorial_part();
        assert_eq!(m.value(1, 1), 1.0);
    }

    #[test]
    fn test_from_trsf2d() {
        let t = NTrsf2d::new();
        let g = NGTrsf2d::from_trsf2d(&t);
        assert_eq!(g.form(), t.form());
        assert_eq!(g.scale, t.scale_factor());
    }

    #[test]
    fn test_from_matrix_vector() {
        let m = NMat2d::new();
        let v = NXY::new(1.0, 2.0);
        let g = NGTrsf2d::from_matrix_vector(&m, &v);
        assert_eq!(g.form(), NTrsfForm::Other);
        assert_eq!(g.scale, 0.0);
        assert_eq!(g.translation_part(), &v);
    }

    #[test]
    fn test_set_affinity() {
        let mut g = gtrsf2d();
        g.set_affinity(&ax2d(), 2.0);
        assert_eq!(g.form(), NTrsfForm::Other);
        assert_eq!(g.scale, 0.0);
    }

    #[test]
    fn test_set_value() {
        let mut g = gtrsf2d();
        g.set_value(1, 1, 2.0).unwrap();
        assert_eq!(g.value(1, 1).unwrap(), 2.0);
        assert_eq!(g.form(), NTrsfForm::Other);
        g.set_value(1, 3, 3.0).unwrap();
        assert_eq!(g.value(1, 3).unwrap(), 3.0);
        assert!(matches!(g.set_value(3, 1, 1.0), Err(NErrors::OutOfRange)));
    }

    #[test]
    fn test_set_translation_part() {
        let mut g = gtrsf2d();
        let coord = NXY::new(1.0, 2.0);
        g.set_translation_part(&coord);
        assert_eq!(g.translation_part(), &coord);
        assert_eq!(g.form(), NTrsfForm::Translation);
    }

    #[test]
    fn test_is_negative() {
        let mut g = gtrsf2d();
        assert!(!g.is_negative());
        let mut m = NMat2d::new();
        m.set_value(1, 1, -1.0);
        g.set_vectorial_part(&m);
        assert!(g.is_negative());
    }

    #[test]
    fn test_transforms() {
        let g = gtrsf2d();
        let mut xy = NXY::new(1.0, 2.0);
        g.transforms_xy(&mut xy);
        assert_eq!(xy, NXY::new(1.0, 2.0));

        let mut x = 1.0;
        let mut y = 2.0;
        g.transforms_coords(&mut x, &mut y);
        assert_eq!(x, 1.0);
        assert_eq!(y, 2.0);
    }

    #[test]
    fn test_to_trsf2d() {
        let g = gtrsf2d();
        let t = g.to_trsf2d().unwrap();
        assert_eq!(t.form(), NTrsfForm::Identity);
    }

    #[test]
    fn test_dump_json() {
        let g = gtrsf2d();
        let mut output = Vec::new();
        g.dump_json(&mut output, -1);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NGTrsf2d\""));
        assert!(json.contains("\"shape\": \"Identity\""));
    }
}
