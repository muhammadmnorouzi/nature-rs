use std::io::Write;

use serde::{Deserialize, Serialize};
use super::prelude::*;

mod gp {
    pub fn resolution() -> f64 {
        1e-12 // Placeholder; replace with actual value if defined
    }
}

// Trait to define the behavior of a 3D transformation
pub trait Trsf {
    fn new() -> Self;
    fn new_from_trsf2d(t: &NTrsf2d) -> Self;
    fn set_mirror_pnt(&mut self, p: &NPoint3d);
    fn set_mirror_ax1(&mut self, a1: &NAx1);
    fn set_mirror_ax2(&mut self, a2: &NAx2);
    fn set_rotation_ax1(&mut self, a1: &NAx1, ang: f64);
    fn set_rotation_quat(&mut self, r: &NQuaternion);
    fn set_rotation_part(&mut self, r: &NQuaternion);
    fn set_scale(&mut self, p: &NPoint3d, s: f64) -> Result<(), NErrors>;
    fn set_displacement(&mut self, from_a1: &NAx3, to_a2: &NAx3);
    fn set_transformation_ax3(&mut self, from_a1: &NAx3, to_a2: &NAx3);
    fn set_transformation_ax3_single(&mut self, to_a2: &NAx3);
    fn set_transformation_quat_vec(&mut self, r: &NQuaternion, t: &NVec);
    fn set_translation_vec(&mut self, v: &NVec);
    fn set_translation_pnts(&mut self, p1: &NPoint3d, p2: &NPoint3d);
    fn set_translation_part(&mut self, v: &NVec);
    fn set_scale_factor(&mut self, s: f64) -> Result<(), NErrors>;
    fn set_form(&mut self, p: NTrsfForm);
    fn set_values(&mut self, a11: f64, a12: f64, a13: f64, a14: f64, a21: f64, a22: f64, a23: f64, a24: f64, a31: f64, a32: f64, a33: f64, a34: f64) -> Result<(), NErrors>;
    fn is_negative(&self) -> bool;
    fn form(&self) -> NTrsfForm;
    fn scale_factor(&self) -> f64;
    fn translation_part(&self) -> NXYZ;
    fn get_rotation(&self) -> NQuaternion;
    fn get_rotation_with_axis(&self) -> (NXYZ, f64);
    fn vectorial_part(&self) -> NMat;
    fn h_vectorial_part(&self) -> NMat;
    fn value(&self, row: i32, col: i32) -> Result<f64, NErrors>;
    fn invert(&mut self) -> Result<(), NErrors>;
    fn inverted(&self) -> Result<Self, NErrors> where Self: Sized;
    fn multiply(&mut self, t: &Self);
    fn multiplied(&self, t: &Self) -> Self where Self: Sized;
    fn pre_multiply(&mut self, t: &Self);
    fn power(&mut self, n: i32) -> Result<(), NErrors>;
    fn powered(&self, n: i32) -> Result<Self, NErrors> where Self: Sized;
    fn transforms_xyz(&self, coord: &mut NXYZ);
    fn transforms_coords(&self, x: &mut f64, y: &mut f64, z: &mut f64);
    fn orthogonalize(&mut self);
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
    fn init_from_json(&mut self, json: &str, pos: &mut usize) -> bool;
}

// Struct representing a 3D transformation
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize, Default)]
pub struct NTrsf {
    scale: f64,
    shape: NTrsfForm,
    matrix: NMat,
    loc: NXYZ,
}

impl Trsf for NTrsf {
    /// Returns the identity transformation.
    fn new() -> Self {
        NTrsf {
            scale: 1.0,
            shape: NTrsfForm::Identity,
            matrix: NMat::new_identity(),
            loc: NXYZ::new(0.0, 0.0, 0.0),
        }
    }

    /// Creates a 3D transformation from a 2D transformation.
    fn new_from_trsf2d(t: &NTrsf2d) -> Self {
        let mut result = NTrsf::new();
        result.scale = t.scale_factor();
        result.shape = t.form();
        let m = t.h_vectorial_part();
        result.matrix.set_value(1, 1, m.value(1, 1)).unwrap();
        result.matrix.set_value(1, 2, m.value(1, 2)).unwrap();
        result.matrix.set_value(2, 1, m.value(2, 1)).unwrap();
        result.matrix.set_value(2, 2, m.value(2, 2)).unwrap();
        result.matrix.set_value(3, 3, 1.0).unwrap();
        result.loc = NXYZ::new(t.translation_part().x(), t.translation_part().y(), 0.0);
        if result.shape == NTrsfForm::Ax1Mirror {
            result.scale = 1.0;
            result.matrix.multiply_scalar(-1.0);
        }
        result
    }

    /// Sets the transformation to a point mirror.
    fn set_mirror_pnt(&mut self, p: &NPoint3d) {
        self.shape = NTrsfForm::PntMirror;
        self.scale = -1.0;
        self.loc = p.xyz();
        self.matrix = NMat::new_identity();
        self.loc.multiply_scalar(2.0);
    }

    /// Sets the transformation to an axial mirror.
    fn set_mirror_ax1(&mut self, a1: &NAx1) {
        self.shape = NTrsfForm::Ax1Mirror;
        self.scale = 1.0;
        self.loc = a1.location().xyz();
        let dir = a1.direction().xyz();
        self.matrix.set_dot(&dir);
        self.matrix.multiply_scalar(-2.0);
        self.matrix.set_diagonal(
            self.matrix.value(1, 1).unwrap() + 1.0,
            self.matrix.value(2, 2).unwrap() + 1.0,
            self.matrix.value(3, 3).unwrap() + 1.0,
        );
        let loc_copy = self.loc.clone();
        self.loc.multiply(&self.matrix);
        self.loc.add(&a1.location().xyz());
        self.matrix.multiply_scalar(-1.0);
    }

    /// Sets the transformation to a planar mirror.
    fn set_mirror_ax2(&mut self, a2: &NAx2) {
        self.shape = NTrsfForm::Ax2Mirror;
        self.scale = -1.0;
        self.loc = a2.location().xyz();
        let dir = a2.direction().xyz();
        self.matrix.set_dot(&dir);
        self.matrix.multiply_scalar(2.0);
        self.matrix.set_diagonal(
            self.matrix.value(1, 1).unwrap() - 1.0,
            self.matrix.value(2, 2).unwrap() - 1.0,
            self.matrix.value(3, 3).unwrap() - 1.0,
        );
        let loc_copy = self.loc.clone();
        self.loc.multiply(&self.matrix);
        self.loc.add(&a2.location().xyz());
    }

    /// Sets the transformation to a rotation around an axis.
    fn set_rotation_ax1(&mut self, a1: &NAx1, ang: f64) {
        self.shape = NTrsfForm::Rotation;
        self.scale = 1.0;
        self.loc = a1.location().xyz();
        self.matrix.set_rotation(&a1.direction().xyz(), ang);
        let mut loc_copy = self.loc.clone();
        loc_copy.reverse();
        loc_copy.multiply(&self.matrix);
        loc_copy.add(&a1.location().xyz());
        self.loc = loc_copy;
    }

    /// Sets the transformation to a rotation defined by a quaternion.
    fn set_rotation_quat(&mut self, r: &NQuaternion) {
        self.shape = NTrsfForm::Rotation;
        self.scale = 1.0;
        self.loc = NXYZ::new(0.0, 0.0, 0.0);
        self.matrix = r.get_matrix();
    }

    /// Replaces the rotation part with the specified quaternion.
    fn set_rotation_part(&mut self, r: &NQuaternion) {
        let has_rotation = !r.is_equal(&NQuaternion::new());
        self.matrix = if has_rotation {
            r.get_matrix()
        } else {
            NMat::new_identity()
        };

        match self.shape {
            NTrsfForm::Identity => {
                if has_rotation {
                    self.shape = NTrsfForm::Rotation;
                }
            }
            NTrsfForm::Rotation => {
                if !has_rotation {
                    self.shape = NTrsfForm::Identity;
                }
            }
            NTrsfForm::Translation | NTrsfForm::PntMirror | NTrsfForm::Ax1Mirror | NTrsfForm::Ax2Mirror | NTrsfForm::Scale | NTrsfForm::CompoundTrsf | NTrsfForm::Other => {
                if has_rotation {
                    self.shape = NTrsfForm::CompoundTrsf;
                }
            }
        }
    }

    /// Sets the transformation to a scale.
    fn set_scale(&mut self, p: &NPoint3d, s: f64) -> Result<(), NErrors> {
        let _as = s.abs();
        if _as <= gp::resolution() {
            return Err(NErrors::InvalidConstructionParameters);
        }
        self.shape = NTrsfForm::Scale;
        self.scale = s;
        self.loc = p.xyz();
        self.matrix = NMat::new_identity();
        self.loc.multiply_scalar(1.0 - s);
        Ok(())
    }

    /// Sets the transformation for displacement between two coordinate systems.
    fn set_displacement(&mut self, from_a1: &NAx3, to_a2: &NAx3) {
        self.shape = NTrsfForm::CompoundTrsf;
        self.scale = 1.0;
        self.matrix.set_cols(
            &to_a2.x_direction().xyz(),
            &to_a2.y_direction().xyz(),
            &to_a2.direction().xyz(),
        );
        self.loc = to_a2.location().xyz();
        let mut ma1 = NMat::new_from_cols(
            &from_a1.x_direction().xyz(),
            &from_a1.y_direction().xyz(),
            &from_a1.direction().xyz(),
        );
        ma1.transpose();
        let mut ma1_loc = from_a1.location().xyz();
        ma1_loc.multiply(&ma1);
        ma1_loc.reverse();
        let mut loc_copy = ma1_loc.clone();
        loc_copy.multiply(&self.matrix);
        self.loc.add(&loc_copy);
        self.matrix.multiply(&ma1);
    }

    /// Sets the transformation between two coordinate systems.
    fn set_transformation_ax3(&mut self, from_a1: &NAx3, to_a2: &NAx3) {
        self.shape = NTrsfForm::CompoundTrsf;
        self.scale = 1.0;
        self.matrix.set_rows(
            &to_a2.x_direction().xyz(),
            &to_a2.y_direction().xyz(),
            &to_a2.direction().xyz(),
        );
        self.loc = to_a2.location().xyz();
        let mut loc_copy = self.loc.clone();
        loc_copy.multiply(&self.matrix);
        loc_copy.reverse();
        let ma1 = NMat::new_from_cols(
            &from_a1.x_direction().xyz(),
            &from_a1.y_direction().xyz(),
            &from_a1.direction().xyz(),
        );
        let mut ma1_loc = from_a1.location().xyz();
        ma1_loc.multiply(&self.matrix);
        self.loc.add(&ma1_loc);
        self.matrix.multiply(&ma1);
    }

    /// Sets the transformation to a coordinate system from the default system.
    fn set_transformation_ax3_single(&mut self, to_a2: &NAx3) {
        self.shape = NTrsfForm::CompoundTrsf;
        self.scale = 1.0;
        self.matrix.set_rows(
            &to_a2.x_direction().xyz(),
            &to_a2.y_direction().xyz(),
            &to_a2.direction().xyz(),
        );
        self.loc = to_a2.location().xyz();
        let mut loc_copy = self.loc.clone();
        loc_copy.multiply(&self.matrix);
        loc_copy.reverse();
        self.loc = loc_copy;
    }

    /// Sets the transformation with a quaternion and translation.
    fn set_transformation_quat_vec(&mut self, r: &NQuaternion, t: &NVec) {
        self.shape = NTrsfForm::CompoundTrsf;
        self.scale = 1.0;
        self.loc = t.xyz();
        self.matrix = r.get_matrix();
    }

    /// Sets the transformation to a translation by a vector.
    fn set_translation_vec(&mut self, v: &NVec) {
        self.shape = NTrsfForm::Translation;
        self.scale = 1.0;
        self.matrix = NMat::new_identity();
        self.loc = v.xyz();
    }

    /// Sets the transformation to a translation between two points.
    fn set_translation_pnts(&mut self, p1: &NPoint3d, p2: &NPoint3d) {
        self.shape = NTrsfForm::Translation;
        self.scale = 1.0;
        self.matrix = NMat::new_identity();
        self.loc = p2.xyz().subtracted(&p1.xyz());
    }

    /// Replaces the translation part with the given vector.
    fn set_translation_part(&mut self, v: &NVec) {
        self.loc = v.xyz();
        let loc_null = self.loc.square_modulus() < gp::resolution();
        match self.shape {
            NTrsfForm::Identity => {
                if !loc_null {
                    self.shape = NTrsfForm::Translation;
                }
            }
            NTrsfForm::Translation => {
                if loc_null {
                    self.shape = NTrsfForm::Identity;
                }
            }
            NTrsfForm::Rotation | NTrsfForm::PntMirror | NTrsfForm::Ax1Mirror | NTrsfForm::Ax2Mirror | NTrsfForm::Scale | NTrsfForm::CompoundTrsf | NTrsfForm::Other => {
                if !loc_null {
                    self.shape = NTrsfForm::CompoundTrsf;
                }
            }
        }
    }

    /// Modifies the scale factor.
    fn set_scale_factor(&mut self, s: f64) -> Result<(), NErrors> {
        let as = s.abs();
        if as <= gp::resolution() {
            return Err(NErrors::InvalidConstructionParameters);
        }
        self.scale = s;
        let unit = (s - 1.0).abs() <= gp::resolution();
        let munit = (s + 1.0).abs() <= gp::resolution();
        match self.shape {
            NTrsfForm::Identity | NTrsfForm::Translation => {
                if !unit {
                    self.shape = NTrsfForm::Scale;
                }
                if munit {
                    self.shape = NTrsfForm::PntMirror;
                }
            }
            NTrsfForm::Rotation => {
                if !unit {
                    self.shape = NTrsfForm::CompoundTrsf;
                }
            }
            NTrsfForm::PntMirror | NTrsfForm::Ax1Mirror | NTrsfForm::Ax2Mirror => {
                if !munit {
                    self.shape = NTrsfForm::Scale;
                }
                if unit {
                    self.shape = NTrsfForm::Identity;
                }
            }
            NTrsfForm::Scale => {
                if unit {
                    self.shape = NTrsfForm::Identity;
                }
                if munit {
                    self.shape = NTrsfForm::PntMirror;
                }
            }
            NTrsfForm::CompoundTrsf | NTrsfForm::Other => {}
        }
        Ok(())
    }

    /// Sets the transformation form.
    fn set_form(&mut self, p: NTrsfForm) {
        self.shape = p;
    }

    /// Sets the transformation coefficients.
    fn set_values(&mut self, a11: f64, a12: f64, a13: f64, a14: f64, a21: f64, a22: f64, a23: f64, a24: f64, a31: f64, a32: f64, a33: f64, a34: f64) -> Result<(), NErrors> {
        let col1 = NXYZ::new(a11, a21, a31);
        let col2 = NXYZ::new(a12, a22, a32);
        let col3 = NXYZ::new(a13, a23, a33);
        let col4 = NXYZ::new(a14, a24, a34);
        let mut m = NMat::new_from_cols(&col1, &col2, &col3);
        let s = m.determinant();
        let as = s.abs();
        if as < gp::resolution() {
            return Err(NErrors::InvalidConstructionParameters);
        }
        self.scale = if s > 0.0 { s.powf(1.0 / 3.0) } else { -((-s).powf(1.0 / 3.0)) };
        m.divide(self.scale);
        self.shape = NTrsfForm::CompoundTrsf;
        self.matrix = m;
        self.orthogonalize();
        self.loc = col4;
        Ok(())
    }

    /// Returns true if the determinant of the vectorial part is negative.
    fn is_negative(&self) -> bool {
        self.scale < 0.0
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
    fn translation_part(&self) -> NXYZ {
        self.loc.clone()
    }

    /// Returns the quaternion representing the rotational part.
    fn get_rotation(&self) -> NQuaternion {
        NQuaternion::new_from_matrix(&self.matrix)
    }

    /// Returns the rotation axis and angle.
    fn get_rotation_with_axis(&self) -> (NXYZ, f64) {
        let q = self.get_rotation();
        let (vec, angle) = q.get_vector_and_angle();
        (vec.xyz(), angle)
    }

    /// Returns the vectorial part of the transformation.
    fn vectorial_part(&self) -> NMat {
        if self.scale == 1.0 {
            self.matrix.clone()
        } else if self.shape == NTrsfForm::Scale || self.shape == NTrsfForm::PntMirror {
            let mut m = self.matrix.clone();
            m.set_diagonal(
                self.scale * m.value(1, 1).unwrap(),
                self.scale * m.value(2, 2).unwrap(),
                self.scale * m.value(3, 3).unwrap(),
            );
            m
        } else {
            let mut m = self.matrix.clone();
            m.multiply_scalar(self.scale);
            m
        }
    }

    /// Returns the homogeneous vectorial part.
    fn h_vectorial_part(&self) -> NMat {
        self.matrix.clone()
    }

    /// Returns the coefficient at the specified row and column.
    fn value(&self, row: i32, col: i32) -> Result<f64, NErrors> {
        if row < 1 || row > 3 || col < 1 || col > 4 {
            return Err(NErrors::OutOfRange);
        }
        if col < 4 {
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
            if t.loc.x() != 0.0 || t.loc.y() != 0.0 || t.loc.z() != 0.0 {
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
            self.loc.add(&t_loc);
            self.matrix.multiply(&t.matrix);
        } else if matches!(self.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror | NTrsfForm::Ax2Mirror) && t.shape == NTrsfForm::Translation {
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
        } else if self.shape == NTrsfForm::Translation && matches!(t.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror | NTrsfForm::Ax2Mirror) {
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
        } else if matches!(self.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror | NTrsfForm::Ax2Mirror) && (t.shape == NTrsfForm::Scale || t.shape == NTrsfForm::PntMirror) {
            self.shape = NTrsfForm::CompoundTrsf;
            let mut t_loc = t.loc.clone();
            if self.scale == 1.0 {
                self.scale = t.scale;
                t_loc.multiply(&self.matrix);
            } else {
                t_loc.multiply(&self.matrix);
                t_loc.multiply_scalar(self.scale);
                self.scale *= t.scale;
            }
            self.loc.add(&t_loc);
        } else if matches!(t.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror | NTrsfForm::Ax2Mirror) && (self.shape == NTrsfForm::Scale || self.shape == NTrsfForm::PntMirror) {
            self.shape = NTrsfForm::CompoundTrsf;
            let mut t_loc = t.loc.clone();
            t_loc.multiply_scalar(self.scale);
            self.loc.add(&t_loc);
            self.scale *= t.scale;
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
            self.loc = loc_copy;
            self.loc.add(&t.loc);
            self.matrix.pre_multiply(&t.matrix);
        } else if matches!(self.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror | NTrsfForm::Ax2Mirror) && t.shape == NTrsfForm::Translation {
            self.loc.add(&t.loc);
        } else if (self.shape == NTrsfForm::Scale || self.shape == NTrsfForm::PntMirror) && t.shape == NTrsfForm::Translation {
            self.loc.add(&t.loc);
        } else if self.shape == NTrsfForm::Translation && matches!(t.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror | NTrsfForm::Ax2Mirror) {
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
        } else if matches!(self.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror | NTrsfForm::Ax2Mirror) && (t.shape == NTrsfForm::Scale || t.shape == NTrsfForm::PntMirror) {
            self.shape = NTrsfForm::CompoundTrsf;
            let mut loc_copy = self.loc.clone();
            loc_copy.multiply_scalar(t.scale);
            self.loc = loc_copy;
            self.loc.add(&t.loc);
            self.scale *= t.scale;
        } else if matches!(t.shape, NTrsfForm::CompoundTrsf | NTrsfForm::Rotation | NTrsfForm::Ax1Mirror | NTrsfForm::Ax2Mirror) && (self.shape == NTrsfForm::Scale || self.shape == NTrsfForm::PntMirror) {
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
            self.matrix = NMat::new_identity();
            self.loc = NXYZ::new(0.0, 0.0, 0.0);
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
            if self.loc.x() == 0.0 && self.loc.y() == 0.0 && self.loc.z() == 0.0 {
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
        } else if matches!(self.shape, NTrsfForm::PntMirror | NTrsfForm::Ax1Mirror | NTrsfForm::Ax2Mirror) {
            if n % 2 == 0 {
                self.shape = NTrsfForm::Identity;
                self.scale = 1.0;
                self.matrix = NMat::new_identity();
                self.loc = NXYZ::new(0.0, 0.0, 0.0);
            }
        } else {
            self.shape = NTrsfForm::CompoundTrsf;
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

    /// Transforms a triplet of coordinates.
    fn transforms_xyz(&self, coord: &mut NXYZ) {
        coord.multiply(&self.matrix);
        if self.scale != 1.0 {
            coord.multiply_scalar(self.scale);
        }
        coord.add(&self.loc);
    }

    /// Transforms individual coordinates.
    fn transforms_coords(&self, x: &mut f64, y: &mut f64, z: &mut f64) {
        let mut triplet = NXYZ::new(*x, *y, *z);
        triplet.multiply(&self.matrix);
        if self.scale != 1.0 {
            triplet.multiply_scalar(self.scale);
        }
        triplet.add(&self.loc);
        *x = triplet.x();
        *y = triplet.y();
        *z = triplet.z();
    }

    /// Orthogonalizes the matrix.
    fn orthogonalize(&mut self) {
        let mut tm = self.matrix.clone();
        let mut v1 = tm.column(1);
        let mut v2 = tm.column(2);
        let mut v3 = tm.column(3);

        v1.normalize();
        v2 = v2.subtracted(&v1.scaled(v2.dot(&v1)));
        v2.normalize();
        v3 = v3.subtracted(&v1.scaled(v3.dot(&v1))).subtracted(&v2.scaled(v3.dot(&v2)));
        v3.normalize();
        tm.set_cols(&v1, &v2, &v3);

        v1 = tm.row(1);
        v2 = tm.row(2);
        v3 = tm.row(3);

        v1.normalize();
        v2 = v2.subtracted(&v1.scaled(v2.dot(&v1)));
        v2.normalize();
        v3 = v3.subtracted(&v1.scaled(v3.dot(&v1))).subtracted(&v2.scaled(v3.dot(&v2)));
        v3.normalize();
        tm.set_rows(&v1, &v2, &v3);

        self.matrix = tm;
    }

    /// Dumps the transformation as JSON.
    fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(out, "{} {{", indent).unwrap();
        writeln!(out, "{}   \"type\": \"NTrsf\",", indent).unwrap();
        writeln!(
            out,
            "{}   \"location\": [{}, {}, {}],",
            indent,
            self.loc.x(),
            self.loc.y(),
            self.loc.z()
        ).unwrap();
        writeln!(
            out,
            "{}   \"matrix\": [{}, {}, {}, {}, {}, {}, {}, {}, {}],",
            indent,
            self.matrix.value(1, 1).unwrap(),
            self.matrix.value(1, 2).unwrap(),
            self.matrix.value(1, 3).unwrap(),
            self.matrix.value(2, 1).unwrap(),
            self.matrix.value(2, 2).unwrap(),
            self.matrix.value(2, 3).unwrap(),
            self.matrix.value(3, 1).unwrap(),
            self.matrix.value(3, 2).unwrap(),
            self.matrix.value(3, 3).unwrap()
        ).unwrap();
        writeln!(out, "{}   \"shape\": {},", indent, self.shape as i32).unwrap();
        writeln!(out, "{}   \"scale\": {}", indent, self.scale).unwrap();
        writeln!(out, "{} }}", indent).unwrap();
    }

    /// Initializes the transformation from JSON.
    fn init_from_json(&mut self, json: &str, pos: &mut usize) -> bool {
        let mut loc = NXYZ::new(0.0, 0.0, 0.0);
        if !Self::init_vector_class(json, pos, "location", 3, &mut [&mut loc.x(), &mut loc.y(), &mut loc.z()]) {
            return false;
        }
        self.set_translation_vec(&NVec::new(loc.x(), loc.y(), loc.z()));

        let mut matrix_vals = [0.0; 9];
        if !Self::init_vector_class(json, pos, "matrix", 9, &mut matrix_vals.iter_mut().collect::<Vec<_>>()) {
            return false;
        }
        for i in 0..3 {
            for j in 0..3 {
                self.matrix.set_value(i + 1, j + 1, matrix_vals[i * 3 + j]).unwrap();
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

impl NTrsf {
    /// Helper method to parse vector class from JSON.
    fn init_vector_class(json: &str, pos: &mut usize, name: &str, count: usize, vals: &mut [&mut f64]) -> bool {
        // Simplified JSON parsing (assumes format like "name": [v1, v2, v3])
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
        // Simplified JSON parsing for single value
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
        let t = NTrsf::new();
        assert_eq!(t.scale_factor(), 1.0);
        assert_eq!(t.form(), NTrsfForm::Identity);
        assert_eq!(t.translation_part(), NXYZ::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn test_set_mirror_pnt() {
        let mut t = NTrsf::new();
        let p = NPoint3d::new(1.0, 2.0, 3.0);
        t.set_mirror_pnt(&p);
        assert_eq!(t.form(), NTrsfForm::PntMirror);
        assert_eq!(t.scale_factor(), -1.0);
        assert_eq!(t.translation_part(), NXYZ::new(2.0, 4.0, 6.0));
    }

    #[test]
    fn test_set_scale() {
        let mut t = NTrsf::new();
        let p = NPoint3d::new(1.0, 0.0, 0.0);
        t.set_scale(&p, 2.0).unwrap();
        assert_eq!(t.form(), NTrsfForm::Scale);
        assert_eq!(t.scale_factor(), 2.0);
        assert_eq!(t.translation_part(), NXYZ::new(-1.0, 0.0, 0.0));
        assert!(t.set_scale(&p, 0.0).is_err());
    }

    #[test]
    fn test_invert() {
        let mut t = NTrsf::new();
        t.set_translation_vec(&NVec::new(1.0, 2.0, 3.0));
        t.invert().unwrap();
        assert_eq!(t.translation_part(), NXYZ::new(-1.0, -2.0, -3.0));
    }

    #[test]
    fn test_multiply() {
        let mut t1 = NTrsf::new();
        let mut t2 = NTrsf::new();
        t1.set_translation_vec(&NVec::new(1.0, 0.0, 0.0));
        t2.set_translation_vec(&NVec::new(0.0, 1.0, 0.0));
        t1.multiply(&t2);
        assert_eq!(t1.form(), NTrsfForm::Translation);
        assert_eq!(t1.translation_part(), NXYZ::new(1.0, 1.0, 0.0));
    }

    #[test]
    fn test_power() {
        let mut t = NTrsf::new();
        t.set_translation_vec(&NVec::new(1.0, 0.0, 0.0));
        t.power(3).unwrap();
        assert_eq!(t.translation_part(), NXYZ::new(3.0, 0.0, 0.0));
    }

    #[test]
    fn test_transforms() {
        let mut t = NTrsf::new();
        t.set_translation_vec(&NVec::new(1.0, 2.0, 3.0));
        let mut coord = NXYZ::new(1.0, 1.0, 1.0);
        t.transforms_xyz(&mut coord);
        assert_eq!(coord, NXYZ::new(2.0, 3.0, 4.0));
    }

    #[test]
    fn test_dump_json() {
        let t = NTrsf::new();
        let mut output = Vec::new();
        t.dump_json(&mut output, 0);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NTrsf\""));
        assert!(json.contains("\"scale\": 1"));
    }
}