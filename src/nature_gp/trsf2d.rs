use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::gp::{NAx2d, NMat2d, NPoint2d, NTrsf, NTrsfForm, NVec2d, NXY, NErrors};

mod gp {
    pub fn resolution() -> f64 {
        1e-12 // Placeholder; replace with actual value if defined
    }
}

// Trait to define the behavior of a 2D transformation
pub trait Trsf2d {
    fn new() -> Self;
    fn new_from_trsf(t: &NTrsf) -> Self;
    fn set_mirror_pnt(&mut self, p: &NPoint2d);
    fn set_mirror_ax2d(&mut self, a: &NAx2d);
    fn set_rotation(&mut self, p: &NPoint2d, ang: f64);
    fn set_scale(&mut self, p: &NPoint2d, s: f64);
    fn set_transformation_ax2d(&mut self, from_a1: &NAx2d, to_a2: &NAx2d);
    fn set_transformation_ax2d_single(&mut self, to_a: &NAx2d);
    fn set_translation_vec(&mut self, v: &NVec2d);
    fn set_translation_pnts(&mut self, p1: &NPoint2d, p2: &NPoint2d);
    fn set_translation_part(&mut self, v: &NVec2d);
    fn set_scale_factor(&mut self, s: f64);
    fn set_values(&mut self, a11: f64, a12: f64, a13: f64, a21: f64, a22: f64, a23: f64) -> Result<(), NErrors>;
    fn is_negative(&self) -> bool;
    fn form(&self) -> NTrsfForm;
    fn scale_factor(&self) -> f64;
    fn translation_part(&self) -> NXY;
    fn vectorial_part(&self) -> NMat2d;
    fn h_vectorial_part(&self) -> NMat2d;
    fn rotation_part(&self) -> f64;
    fn value(&self, row: i32, col: i32) -> Result<f64, NErrors>;
    fn invert(&mut self) -> Result<(), NErrors>;
    fn inverted(&self) -> Result<Self, NErrors> where Self: Sized;
    fn multiply(&mut self, t: &Self);
    fn multiplied(&self, t: &Self) -> Self where Self: Sized;
    fn pre_multiply(&mut self, t: &Self);
    fn power(&mut self, n: i32) -> Result<(), NErrors>;
    fn powered(&self, n: i32) -> Result<Self, NErrors> where Self: Sized;
    fn transforms_xy(&self, coord: &mut NXY);
    fn transforms_coords(&self, x: &mut f64, y: &mut f64);
    fn orthogonalize(&mut self);
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
    fn init_from_json(&mut self, json: &str, pos: &mut usize) -> bool;
}

// Struct representing a 2D transformation
#[derive(Clone, Default,Debug, PartialEq, Serialize, Deserialize)]
pub struct NTrsf2d {
    scale: f64,
    shape: NTrsfForm,
    matrix: NMat2d,
    loc: NXY,
}

impl Trsf2d for NTrsf2d {
    /// Returns the identity transformation.
    fn new() -> Self {
        NTrsf2d {
            scale: 1.0,
            shape: NTrsfForm::Identity,
            matrix: NMat2d::new_identity(),
            loc: NXY::new(0.0, 0.0),
        }
    }

    /// Creates a 2D transformation from a 3D transformation.
    fn new_from_trsf(t: &NTrsf) -> Self {
        let mut result = NTrsf2d::new();
        result.scale = t.scale_factor();
        result.shape = t.form();
        let m = t.h_vectorial_part();
        result.matrix.set_value(1, 1, m.value(1, 1).unwrap()).unwrap();
        result.matrix.set_value(1, 2, m.value(1, 2).unwrap()).unwrap();
        result.matrix.set_value(2, 1, m.value(2, 1).unwrap()).unwrap();
        result.matrix.set_value(2, 2, m.value(2, 2).unwrap()).unwrap();
        result.loc = NXY::new(t.translation_part().x(), t.translation_part().y());
        result
    }

    /// Sets the transformation to a point mirror.
    fn set_mirror_pnt(&mut self, p: &NPoint2d) {
        self.shape = NTrsfForm::PntMirror;
        self.scale = -1.0;
        self.matrix = NMat2d::new_identity();
        self.loc = p.xy();
        self.loc.multiply_scalar(2.0);
    }

    /// Sets the transformation to an axial mirror.
    fn set_mirror_ax2d(&mut self, a: &NAx2d) {
        self.shape = NTrsfForm::Ax1Mirror;
        self.scale = -1.0;
        let v = a.direction();
        let p = a.location();
        let vx = v.x();
        let vy = v.y();
        let x0 = p.x();
        let y0 = p.y();
        self.matrix.set_col(1, &NXY::new(1.0 - 2.0 * vx * vx, -2.0 * vx * vy));
        self.matrix.set_col(2, &NXY::new(-2.0 * vx * vy, 1.0 - 2.0 * vy * vy));
        self.loc.set_coord(
            -2.0 * ((vx * vx - 1.0) * x0 + (vx * vy * y0)),
            -2.0 * ((vx * vy * x0) + (vy * vy - 1.0) * y0),
        );
    }

    /// Sets the transformation to a rotation.
    fn set_rotation(&mut self, p: &NPoint2d, ang: f64) {
        self.shape = NTrsfForm::Rotation;
        self.scale = 1.0;
        self.loc = p.xy();
        let mut loc_copy = self.loc.clone();
        loc_copy.reverse();
        self.matrix.set_rotation(ang);
        loc_copy.multiply(&self.matrix);
        loc_copy.add(&p.xy());
        self.loc = loc_copy;
    }

    /// Sets the transformation to a scale.
    fn set_scale(&mut self, p: &NPoint2d, s: f64) {
        self.shape = NTrsfForm::Scale;
        self.scale = s;
        self.matrix = NMat2d::new_identity();
        self.loc = p.xy();
        self.loc.multiply_scalar(1.0 - s);
    }

    /// Sets the transformation for passage between two coordinate systems.
    fn set_transformation_ax2d(&mut self, from_a1: &NAx2d, to_a2: &NAx2d) {
        self.shape = NTrsfForm::CompoundTrsf;
        self.scale = 1.0;
        let v1 = to_a2.direction().xy();
        let v2 = NXY::new(-v1.y(), v1.x());
        self.matrix.set_col(1, &v1);
        self.matrix.set_col(2, &v2);
        self.loc = to_a2.location().xy();
        let mut matrix_copy = self.matrix.clone();
        matrix_copy.transpose();
        self.loc.multiply(&matrix_copy);
        self.loc.reverse();
        let v3 = from_a1.direction().xy();
        let v4 = NXY::new(-v3.y(), v3.x());
        let ma1 = NMat2d::new_from_cols(&v3, &v4);
        let mut ma1_loc = from_a1.location().xy();
        ma1_loc.multiply(&self.matrix);
        self.loc.add(&ma1_loc);
        self.matrix.multiply(&ma1);
    }

    /// Sets the transformation to a coordinate system from the default system.
    fn set_transformation_ax2d_single(&mut self, to_a: &NAx2d) {
        self.shape = NTrsfForm::CompoundTrsf;
        self.scale = 1.0;
        let v1 = to_a.direction().xy();
        let v2 = NXY::new(-v1.y(), v1.x());
        self.matrix.set_col(1, &v1);
        self.matrix.set_col(2, &v2);
        self.loc = to_a.location().xy();
        let mut matrix_copy = self.matrix.clone();
        matrix_copy.transpose();
        self.loc.multiply(&matrix_copy);
        self.loc.reverse();
    }

    /// Sets the transformation to a translation by a vector.
    fn set_translation_vec(&mut self, v: &NVec2d) {
        self.shape = NTrsfForm::Translation;
        self.scale = 1.0;
        self.matrix = NMat2d::new_identity();
        self.loc = v.xy();
    }

    /// Sets the transformation to a translation between two points.
    fn set_translation_pnts(&mut self, p1: &NPoint2d, p2: &NPoint2d) {
        self.shape = NTrsfForm::Translation;
        self.scale = 1.0;
        self.matrix = NMat2d::new_identity();
        self.loc = p2.xy().subtracted(&p1.xy());
    }

    /// Replaces the translation part with the given vector.
    fn set_translation_part(&mut self, v: &NVec2d) {
        self.loc = v.xy();
        let x = self.loc.x().abs();
        let y = self.loc.y().abs();
        let loc_null = x <= gp::resolution() && y <= gp::resolution();
        if loc_null {
            if matches!(self.shape, NTrsfForm::Identity | NTrsfForm::PntMirror | NTrsfForm::Scale | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror) {
            } else if self.shape == NTrsfForm::Translation {
                self.shape = NTrsfForm::Identity;
            } else {
                self.shape = NTrsfForm::CompoundTrsf;
            }
        } else {
            if matches!(self.shape, NTrsfForm::Translation | NTrsfForm::Scale | NTrsfForm::PntMirror) {
            } else if self.shape == NTrsfForm::Identity {
                self.shape = NTrsfForm::Translation;
            } else {
                self.shape = NTrsfForm::CompoundTrsf;
            }
        }
    }

    /// Modifies the scale factor.
    fn set_scale_factor(&mut self, s: f64) {
        if s == 1.0 {
            let x = self.loc.x().abs();
            let y = self.loc.y().abs();
            let loc_null = x <= gp::resolution() && y <= gp::resolution();
            if loc_null {
                if matches!(self.shape, NTrsfForm::Identity | NTrsfForm::Rotation) {
                } else if self.shape == NTrsfForm::Scale {
                    self.shape = NTrsfForm::Identity;
                } else if self.shape == NTrsfForm::PntMirror {
                    self.shape = NTrsfForm::Translation;
                } else {
                    self.shape = NTrsfForm::CompoundTrsf;
                }
            } else {
                if matches!(self.shape, NTrsfForm::Identity | NTrsfForm::Rotation | NTrsfForm::Scale) {
                } else if self.shape == NTrsfForm::PntMirror {
                    self.shape = NTrsfForm::Translation;
                } else {
                    self.shape = NTrsfForm::CompoundTrsf;
                }
            }
        } else if s == -1.0 {
            if matches!(self.shape, NTrsfForm::PntMirror | NTrsfForm::Ax1Mirror) {
            } else if matches!(self.shape, NTrsfForm::Identity | NTrsfForm::Scale) {
                self.shape = NTrsfForm::PntMirror;
            } else {
                self.shape = NTrsfForm::CompoundTrsf;
            }
        } else {
            if self.shape == NTrsfForm::Scale {
            } else if matches!(self.shape, NTrsfForm::Identity | NTrsfForm::Translation | NTrsfForm::PntMirror) {
                self.shape = NTrsfForm::Scale;
            } else {
                self.shape = NTrsfForm::CompoundTrsf;
            }
        }
        self.scale = s;
    }

    /// Sets the transformation coefficients.
    fn set_values(&mut self, a11: f64, a12: f64, a13: f64, a21: f64, a22: f64, a23: f64) -> Result<(), NErrors> {
        let col1 = NXY::new(a11, a21);
        let col2 = NXY::new(a12, a22);
        let col3 = NXY::new(a13, a23);
        let mut m = NMat2d::new_from_cols(&col1, &col2);
        let s = m.determinant();
        let as = s.abs();
        if as < gp::resolution() {
            return Err(NErrors::InvalidConstructionParameters);
        }
        self.scale = if s > 0.0 { s.sqrt() } else { (-s).sqrt() };
        m.divide(self.scale);
        self.shape = NTrsfForm::CompoundTrsf;
        self.matrix = m;
        self.orthogonalize();
        self.loc = col3;
        Ok(())
    }

    /// Returns true if the determinant of the vectorial part is negative.
    fn is_negative(&self) -> bool {
        self.matrix.determinant() < 0.0
    }

    /// Returns the transformation form.
    fn form(&self) -> NTrsfForm {
        self.shape
    }

    /// Returns the scale factor.
    fn scale_factor(&self) -> f64 {
        self.scale
    }

    /// Returns the translation part.
    fn translation_part(&self) -> NXY {
        self.loc.clone()
    }

    /// Returns the vectorial part of the transformation.
    fn vectorial_part(&self) -> NMat2d {
        if self.scale == 1.0 {
            self.matrix.clone()
        } else if self.shape == NTrsfForm::Scale || self.shape == NTrsfForm::PntMirror {
            let mut m = self.matrix.clone();
            m.set_diagonal(
                self.matrix.value(1, 1).unwrap() * self.scale,
                self.matrix.value(2, 2).unwrap() * self.scale,
            );
            m
        } else {
            let mut m = self.matrix.clone();
            m.multiply_scalar(self.scale);
            m
        }
    }

    /// Returns the homogeneous vectorial part.
    fn h_vectorial_part(&self) -> NMat2d {
        self.matrix.clone()
    }

    /// Returns the rotation angle.
    fn rotation_part(&self) -> f64 {
        self.matrix.value(2, 1).unwrap().atan2(self.matrix.value(1, 1).unwrap())
    }

    /// Returns the coefficient at the specified row and column.
    fn value(&self, row: i32, col: i32) -> Result<f64, NErrors> {
        if row < 1 || row > 2 || col < 1 || col > 3 {
            return Err(NErrors::OutOfRange);
        }
        if col < 3 {
            Ok(self.scale * self.matrix.value(row as usize, col as usize).unwrap())
        } else {
            Ok(self.loc.coord(row as usize))
        }
    }

    /// Inverts the transformation.
    fn invert(&mut self) -> Result<(), NErrors> {
        if self.shape == NTrsfForm::Identity {
            Ok(())
        } else if self.shape == NTrsfForm::Translation || self.shape == NTrsfForm::PntMirror {
            self.loc.reverse();
            Ok(())
        } else if self.shape == NTrsfForm::Scale {
            if self.scale.abs() <= gp::resolution() {
                return Err(NErrors::InvalidConstructionParameters);
            }
            self.scale = 1.0 / self.scale;
            self.loc.multiply_scalar(-self.scale);
            Ok(())
        } else {
            if self.scale.abs() <= gp::resolution() {
                return Err(NErrors::InvalidConstructionParameters);
            }
            self.scale = 1.0 / self.scale;
            self.matrix.transpose();
            let mut loc_copy = self.loc.clone();
            loc_copy.multiply(&self.matrix);
            self.loc = loc_copy;
            self.loc.multiply_scalar(-self.scale);
            Ok(())
        }
    }

    /// Returns the inverted transformation.
    fn inverted(&self) -> Result<Self, NErrors> {
        let mut result = self.clone();
        result.invert()?;
        Ok(result)
    }

    /// Composes the transformation with another.
    fn multiply(&mut self, t: &Self) {
        if t.shape == NTrsfForm::Identity {
            return;
        }
        if self.shape == NTrsfForm::Identity {
            self.shape = t.shape;
            self.scale = t.scale;
            self.loc = t.loc.clone();
            self.matrix = t.matrix.clone();
            return;
        }
        if self.shape == NTrsfForm::Rotation && t.shape == NTrsfForm::Rotation {
            if self.loc.x() != 0.0 || self.loc.y() != 0.0 {
                let mut t_loc = t.loc.clone();
                t_loc.multiply(&self.matrix);
                self.loc.add(&t_loc);
            }
            self.matrix.multiply(&t.matrix);
        } else if self.shape == NTrsfForm::Translation && t.shape == NTrsfForm::Translation {
            self.loc.add(&t.loc);
        } else if self.shape == NTrsfForm::Scale && t.shape == NTrsfForm::Scale {
            let mut t_loc = t.loc.clone();
            t_loc.multiply_scalar(self.scale);
            self.loc.add(&t_loc);
            self.scale *= t.scale;
        } else if self.shape == NTrsfForm::PntMirror && t.shape == NTrsfForm::PntMirror {
            self.scale = 1.0;
            self.shape = NTrsfForm::Translation;
            let mut t_loc = t.loc.clone();
            t_loc.reverse();
            self.loc.add(&t_loc);
        } else if self.shape == NTrsfForm::Ax1Mirror && t.shape == NTrsfForm::Ax1Mirror {
            self.shape = NTrsfForm::Rotation;
            let mut t_loc = t.loc.clone();
            t_loc.multiply(&self.matrix);
            t_loc.multiply_scalar(self.scale);
            self.loc.add(&t_loc);
            self.scale *= t.scale;
            self.matrix.multiply(&t.matrix);
        } else if matches!(self.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror) && t.shape == NTrsfForm::Translation {
            let mut t_loc = t.loc.clone();
            t_loc.multiply(&self.matrix);
            if self.scale != 1.0 {
                t_loc.multiply_scalar(self.scale);
            }
            self.loc.add(&t_loc);
        } else if (self.shape == NTrsfForm::Scale || self.shape == NTrsfForm::PntMirror) && t.shape == NTrsfForm::Translation {
            let mut t_loc = t.loc.clone();
            t_loc.multiply_scalar(self.scale);
            self.loc.add(&t_loc);
        } else if self.shape == NTrsfForm::Translation && matches!(t.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror) {
            self.shape = NTrsfForm::CompoundTrsf;
            self.scale = t.scale;
            self.loc.add(&t.loc);
            self.matrix = t.matrix.clone();
        } else if self.shape == NTrsfForm::Translation && (t.shape == NTrsfForm::Scale || t.shape == NTrsfForm::PntMirror) {
            self.shape = t.shape;
            self.loc.add(&t.loc);
            self.scale = t.scale;
        } else if (self.shape == NTrsfForm::PntMirror || self.shape == NTrsfForm::Scale) && (t.shape == NTrsfForm::PntMirror || t.shape == NTrsfForm::Scale) {
            self.shape = NTrsfForm::CompoundTrsf;
            let mut t_loc = t.loc.clone();
            t_loc.multiply_scalar(self.scale);
            self.loc.add(&t_loc);
            self.scale *= t.scale;
        } else if matches!(self.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror) && (t.shape == NTrsfForm::Scale || t.shape == NTrsfForm::PntMirror) {
            self.shape = NTrsfForm::CompoundTrsf;
            let mut t_loc = t.loc.clone();
            t_loc.multiply(&self.matrix);
            if self.scale == 1.0 {
                self.scale = t.scale;
            } else {
                t_loc.multiply_scalar(self.scale);
                self.scale *= t.scale;
            }
            self.loc.add(&t_loc);
        } else if matches!(t.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror) && (self.shape == NTrsfForm::Scale || self.shape == NTrsfForm::PntMirror) {
            self.shape = NTrsfForm::CompoundTrsf;
            let mut t_loc = t.loc.clone();
            t_loc.multiply_scalar(self.scale);
            self.scale *= t.scale;
            self.loc.add(&t_loc);
            self.matrix = t.matrix.clone();
        } else {
            self.shape = NTrsfForm::CompoundTrsf;
            let mut t_loc = t.loc.clone();
            t_loc.multiply(&self.matrix);
            if self.scale != 1.0 {
                t_loc.multiply_scalar(self.scale);
                self.scale *= t.scale;
            } else {
                self.scale = t.scale;
            }
            self.loc.add(&t_loc);
            self.matrix.multiply(&t.matrix);
        }
    }

    /// Returns the transformation composed with another.
    fn multiplied(&self, t: &Self) -> Self {
        let mut result = self.clone();
        result.multiply(t);
        result
    }

    /// Composes the transformation with another (pre-multiply).
    fn pre_multiply(&mut self, t: &Self) {
        if t.shape == NTrsfForm::Identity {
            return;
        }
        if self.shape == NTrsfForm::Identity {
            self.shape = t.shape;
            self.scale = t.scale;
            self.loc = t.loc.clone();
            self.matrix = t.matrix.clone();
            return;
        }
        if self.shape == NTrsfForm::Rotation && t.shape == NTrsfForm::Rotation {
            let mut loc_copy = self.loc.clone();
            loc_copy.multiply(&t.matrix);
            self.loc = loc_copy;
            self.loc.add(&t.loc);
            self.matrix.pre_multiply(&t.matrix);
        } else if self.shape == NTrsfForm::Translation && t.shape == NTrsfForm::Translation {
            self.loc.add(&t.loc);
        } else if self.shape == NTrsfForm::Scale && t.shape == NTrsfForm::Scale {
            let mut loc_copy = self.loc.clone();
            loc_copy.multiply_scalar(t.scale);
            self.loc = loc_copy;
            self.loc.add(&t.loc);
            self.scale *= t.scale;
        } else if self.shape == NTrsfForm::PntMirror && t.shape == NTrsfForm::PntMirror {
            self.scale = 1.0;
            self.shape = NTrsfForm::Translation;
            let mut loc_copy = self.loc.clone();
            loc_copy.reverse();
            self.loc = loc_copy;
            self.loc.add(&t.loc);
        } else if self.shape == NTrsfForm::Ax1Mirror && t.shape == NTrsfForm::Ax1Mirror {
            self.shape = NTrsfForm::Rotation;
            let mut loc_copy = self.loc.clone();
            loc_copy.multiply(&t.matrix);
            loc_copy.multiply_scalar(t.scale);
            self.loc = loc_copy;
            self.scale *= t.scale;
            self.loc.add(&t.loc);
            self.matrix.pre_multiply(&t.matrix);
        } else if matches!(self.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror) && t.shape == NTrsfForm::Translation {
            self.loc.add(&t.loc);
        } else if (self.shape == NTrsfForm::Scale || self.shape == NTrsfForm::PntMirror) && t.shape == NTrsfForm::Translation {
            self.loc.add(&t.loc);
        } else if self.shape == NTrsfForm::Translation && matches!(t.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror) {
            self.shape = NTrsfForm::CompoundTrsf;
            self.matrix = t.matrix.clone();
            let mut loc_copy = self.loc.clone();
            if t.scale == 1.0 {
                loc_copy.multiply(&t.matrix);
            } else {
                self.scale = t.scale;
                loc_copy.multiply(&self.matrix);
                loc_copy.multiply_scalar(self.scale);
            }
            self.loc = loc_copy;
            self.loc.add(&t.loc);
        } else if (t.shape == NTrsfForm::Scale || t.shape == NTrsfForm::PntMirror) && self.shape == NTrsfForm::Translation {
            let mut loc_copy = self.loc.clone();
            loc_copy.multiply_scalar(t.scale);
            self.loc = loc_copy;
            self.loc.add(&t.loc);
            self.scale = t.scale;
            self.shape = t.shape;
        } else if (self.shape == NTrsfForm::PntMirror || self.shape == NTrsfForm::Scale) && (t.shape == NTrsfForm::PntMirror || t.shape == NTrsfForm::Scale) {
            self.shape = NTrsfForm::CompoundTrsf;
            let mut loc_copy = self.loc.clone();
            loc_copy.multiply_scalar(t.scale);
            self.loc = loc_copy;
            self.loc.add(&t.loc);
            self.scale *= t.scale;
        } else if matches!(self.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror) && (t.shape == NTrsfForm::Scale || t.shape == NTrsfForm::PntMirror) {
            self.shape = NTrsfForm::CompoundTrsf;
            let mut loc_copy = self.loc.clone();
            loc_copy.multiply_scalar(t.scale);
            self.loc = loc_copy;
            self.loc.add(&t.loc);
            self.scale *= t.scale;
        } else if matches!(t.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror) && (self.shape == NTrsfForm::Scale || self.shape == NTrsfForm::PntMirror) {
            self.shape = NTrsfForm::CompoundTrsf;
            self.matrix = t.matrix.clone();
            let mut loc_copy = self.loc.clone();
            if t.scale == 1.0 {
                loc_copy.multiply(&t.matrix);
            } else {
                loc_copy.multiply(&self.matrix);
                loc_copy.multiply_scalar(t.scale);
                self.scale *= t.scale;
            }
            self.loc = loc_copy;
            self.loc.add(&t.loc);
        } else {
            self.shape = NTrsfForm::CompoundTrsf;
            let mut loc_copy = self.loc.clone();
            loc_copy.multiply(&t.matrix);
            if t.scale != 1.0 {
                loc_copy.multiply_scalar(t.scale);
                self.scale *= t.scale;
            }
            self.loc = loc_copy;
            self.loc.add(&t.loc);
            self.matrix.pre_multiply(&t.matrix);
        }
    }

    /// Computes the transformation raised to a power.
    fn power(&mut self, n: i32) -> Result<(), NErrors> {
        if self.shape == NTrsfForm::Identity {
            return Ok(());
        }
        if n == 0 {
            self.scale = 1.0;
            self.shape = NTrsfForm::Identity;
            self.matrix = NMat2d::new_identity();
            self.loc = NXY::new(0.0, 0.0);
            return Ok(());
        }
        if n == 1 {
            return Ok(());
        }
        if n == -1 {
            return self.invert();
        }
        let mut n_power = if n < 0 {
            self.invert()?;
            -n
        } else {
            n
        };
        n_power -= 1;
        if self.shape == NTrsfForm::Translation {
            let temp_loc = self.loc.clone();
            for _ in 0..n_power {
                if n_power % 2 == 1 {
                    self.loc.add(&temp_loc);
                }
                let mut temp = temp_loc.clone();
                temp.add(&temp_loc);
                n_power /= 2;
                if n_power <= 1 {
                    break;
                }
            }
        } else if self.shape == NTrsfForm::Scale {
            let mut temp_loc = self.loc.clone();
            let mut temp_scale = self.scale;
            for _ in 0..n_power {
                if n_power % 2 == 1 {
                    let mut scaled_loc = temp_loc.clone();
                    scaled_loc.multiply_scalar(self.scale);
                    self.loc.add(&scaled_loc);
                    self.scale *= temp_scale;
                }
                let mut temp = temp_loc.clone();
                temp.multiply_scalar(temp_scale);
                temp_loc.add(&temp);
                temp_scale *= temp_scale;
                n_power /= 2;
                if n_power <= 1 {
                    break;
                }
            }
        } else if self.shape == NTrsfForm::Rotation {
            let mut temp_matrix = self.matrix.clone();
            if self.loc.x() == 0.0 && self.loc.y() == 0.0 {
                for _ in 0..n_power {
                    if n_power % 2 == 1 {
                        self.matrix.multiply(&temp_matrix);
                    }
                    let mut temp = temp_matrix.clone();
                    temp.multiply(&temp_matrix);
                    temp_matrix = temp;
                    n_power /= 2;
                    if n_power <= 1 {
                        break;
                    }
                }
            } else {
                let mut temp_loc = self.loc.clone();
                for _ in 0..n_power {
                    if n_power % 2 == 1 {
                        let mut scaled_loc = temp_loc.clone();
                        scaled_loc.multiply(&self.matrix);
                        self.loc.add(&scaled_loc);
                        self.matrix.multiply(&temp_matrix);
                    }
                    let mut temp = temp_loc.clone();
                    temp.multiply(&temp_matrix);
                    temp_loc.add(&temp);
                    let mut temp_mat = temp_matrix.clone();
                    temp_mat.multiply(&temp_matrix);
                    temp_matrix = temp_mat;
                    n_power /= 2;
                    if n_power <= 1 {
                        break;
                    }
                }
            }
        } else if matches!(self.shape, NTrsfForm::PntMirror | NTrsfForm::Ax1Mirror) {
            if n % 2 == 0 {
                self.shape = NTrsfForm::Identity;
                self.scale = 1.0;
                self.matrix = NMat2d::new_identity();
                self.loc = NXY::new(0.0, 0.0);
            }
        } else {
            self.shape = NTrsfForm::CompoundTrsf;
            self.matrix.set_diagonal(
                self.scale * self.matrix.value(1, 1).unwrap(),
                self.scale * self.matrix.value(2, 2).unwrap(),
            );
            let mut temp_loc = self.loc.clone();
            let mut temp_scale = self.scale;
            let mut temp_matrix = self.matrix.clone();
            for _ in 0..n_power {
                if n_power % 2 == 1 {
                    let mut scaled_loc = temp_loc.clone();
                    scaled_loc.multiply(&self.matrix);
                    scaled_loc.multiply_scalar(self.scale);
                    self.loc.add(&scaled_loc);
                    self.scale *= temp_scale;
                    self.matrix.multiply(&temp_matrix);
                }
                let mut temp = temp_loc.clone();
                temp.multiply(&temp_matrix);
                temp.multiply_scalar(temp_scale);
                temp_loc.add(&temp);
                temp_scale *= temp_scale;
                let mut temp_mat = temp_matrix.clone();
                temp_mat.multiply(&temp_matrix);
                temp_matrix = temp_mat;
                n_power /= 2;
                if n_power <= 1 {
                    break;
                }
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

    /// Transforms a doublet of coordinates.
    fn transforms_xy(&self, coord: &mut NXY) {
        coord.multiply(&self.matrix);
        if self.scale != 1.0 {
            coord.multiply_scalar(self.scale);
        }
        coord.add(&self.loc);
    }

    /// Transforms individual coordinates.
    fn transforms_coords(&self, x: &mut f64, y: &mut f64) {
        let mut doublet = NXY::new(*x, *y);
        doublet.multiply(&self.matrix);
        if self.scale != 1.0 {
            doublet.multiply_scalar(self.scale);
        }
        doublet.add(&self.loc);
        *x = doublet.x();
        *y = doublet.y();
    }

    /// Orthogonalizes the matrix.
    fn orthogonalize(&mut self) {
        let mut tm = self.matrix.clone();
        let mut v1 = tm.column(1);
        let mut v2 = tm.column(2);

        v1.normalize();
        v2 = v2.subtracted(&v1.scaled(v2.dot(&v1)));
        v2.normalize();
        tm.set_cols(&v1, &v2);

        v1 = tm.row(1);
        v2 = tm.row(2);

        v1.normalize();
        v2 = v2.subtracted(&v1.scaled(v2.dot(&v1)));
        v2.normalize();
        tm.set_rows(&v1, &v2);

        self.matrix = tm;
    }

    /// Dumps the transformation as JSON.
    fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(out, "{} {{", indent).unwrap();
        writeln!(out, "{}   \"type\": \"NTrsf2d\",", indent).unwrap();
        writeln!(
            out,
            "{}   \"location\": [{}, {}],",
            indent,
            self.loc.x(),
            self.loc.y()
        ).unwrap();
        writeln!(
            out,
            "{}   \"matrix\": [{}, {}, {}, {}],",
            indent,
            self.matrix.value(1, 1).unwrap(),
            self.matrix.value(1, 2).unwrap(),
            self.matrix.value(2, 1).unwrap(),
            self.matrix.value(2, 2).unwrap()
        ).unwrap();
        writeln!(out, "{}   \"shape\": {},", indent, self.shape as i32).unwrap();
        writeln!(out, "{}   \"scale\": {}", indent, self.scale).unwrap();
        writeln!(out, "{} }}", indent).unwrap();
    }

    /// Initializes the transformation from JSON.
    fn init_from_json(&mut self, json: &str, pos: &mut usize) -> bool {
        let mut loc = NXY::new(0.0, 0.0);
        if !Self::init_vector_class(json, pos, "location", 2, &mut [&mut loc.x(), &mut loc.y()]) {
            return false;
        }
        self.set_translation_vec(&NVec2d::new(loc.x(), loc.y()));

        let mut matrix_vals = [0.0; 4];
        if !Self::init_vector_class(json, pos, "matrix", 4, &mut matrix_vals.iter_mut().collect::<Vec<_>>()) {
            return false;
        }
        for i in 0..2 {
            for j in 0..2 {
                self.matrix.set_value(i + 1, j + 1, matrix_vals[i * 2 + j]).unwrap();
            }
        }

        let mut shape_val: i32 = 0;
        if !Self::init_field_value(json, pos, &mut shape_val) {
            return false;
        }
        self.shape = match shape_val {
            0 => NTrsfForm::Identity,
            1 => NTrsfForm::Rotation,
            2 => NTrsfForm::Translation,
            3 => NTrsfForm::PntMirror,
            4 => NTrsfForm::Ax1Mirror,
            5 => NTrsfForm::Ax2Mirror,
            6 => NTrsfForm::Scale,
            7 => NTrsfForm::CompoundTrsf,
            _ => NTrsfForm::Other,
        };

        let mut scale: f64 = 0.0;
        if !Self::init_field_value(json, pos, &mut scale) {
            return false;
        }
        self.scale = scale;

        true
    }
}

impl NTrsf2d {
    /// Helper method to parse vector class from JSON.
    fn init_vector_class(json: &str, pos: &mut usize, name: &str, count: usize, vals: &mut [&mut f64]) -> bool {
        let search = format!("\"{}\": [", name);
        if let Some(start) = json[*pos..].find(&search) {
            *pos += start + search.len();
            let mut num_str = String::new();
            let mut val_idx = 0;
            while val_idx < count {
                let c = json.chars().nth(*pos).unwrap();
                if c.is_digit(10) || c == '.' || c == '-' {
                    num_str.push(c);
                } else if c == ',' || c == ']' {
                    if !num_str.is_empty() {
                        *vals[val_idx] = num_str.parse().unwrap_or(0.0);
                        val_idx += 1;
                        num_str.clear();
                    }
                }
                *pos += 1;
                if c == ']' {
                    break;
                }
            }
            val_idx == count
        } else {
            false
        }
    }

    /// Helper method to parse field value from JSON.
    fn init_field_value<T: std::str::FromStr>(json: &str, pos: &mut usize, val: &mut T) -> bool {
        let mut num_str = String::new();
        let mut found_value = false;
        while *pos < json.len() {
            let c = json.chars().nth(*pos).unwrap();
            if c.is_digit(10) || c == '.' || c == '-' {
                num_str.push(c);
                found_value = true;
            } else if found_value && (c == ',' || c == '}' || c.is_whitespace()) {
                if let Ok(parsed) = num_str.parse() {
                    *val = parsed;
                    return true;
                }
                return false;
            }
            *pos += 1;
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let t = NTrsf2d::new();
        assert_eq!(t.scale_factor(), 1.0);
        assert_eq!(t.form(), NTrsfForm::Identity);
        assert_eq!(t.translation_part(), NXY::new(0.0, 0.0));
    }

    #[test]
    fn test_set_mirror_pnt() {
        let mut t = NTrsf2d::new();
        let p = NPoint2d::new(1.0, 2.0);
        t.set_mirror_pnt(&p);
        assert_eq!(t.form(), NTrsfForm::PntMirror);
        assert_eq!(t.scale_factor(), -1.0);
        assert_eq!(t.translation_part(), NXY::new(2.0, 4.0));
    }

    #[test]
    fn test_set_scale() {
        let mut t = NTrsf2d::new();
        let p = NPoint2d::new(1.0, 0.0);
        t.set_scale(&p, 2.0);
        assert_eq!(t.form(), NTrsfForm::Scale);
        assert_eq!(t.scale_factor(), 2.0);
        assert_eq!(t.translation_part(), NXY::new(-1.0, 0.0));
    }

    #[test]
    fn test_invert() {
        let mut t = NTrsf2d::new();
        t.set_translation_vec(&NVec2d::new(1.0, 2.0));
        t.invert().unwrap();
        assert_eq!(t.translation_part(), NXY::new(-1.0, -2.0));
    }

    #[test]
    fn test_multiply() {
        let mut t1 = NTrsf2d::new();
        let mut t2 = NTrsf2d::new();
        t1.set_translation_vec(&NVec2d::new(1.0, 0.0));
        t2.set_translation_vec(&NVec2d::new(0.0, 1.0));
        t1.multiply(&t2);
        assert_eq!(t1.form(), NTrsfForm::Translation);
        assert_eq!(t1.translation_part(), NXY::new(1.0, 1.0));
    }

    #[test]
    fn test_power() {
        let mut t = NTrsf2d::new();
        t.set_translation_vec(&NVec2d::new(1.0, 0.0));
        t.power(3).unwrap();
        assert_eq!(t.translation_part(), NXY::new(3.0, 0.0));
    }

    #[test]
    fn test_transforms() {
        let mut t = NTrsf2d::new();
        t.set_translation_vec(&NVec2d::new(1.0, 2.0));
        let mut coord = NXY::new(1.0, 1.0);
        t.transforms_xy(&mut coord);
        assert_eq!(coord, NXY::new(2.0, 3.0));
    }

    #[test]
    fn test_dump_json() {
        let t = NTrsf2d::new();
        let mut output = Vec::new();
        t.dump_json(&mut output, 0);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NTrsf2d\""));
        assert!(json.contains("\"scale\": 1"));
    }
}