use std::f64::consts::{PI, SQRT_2};

use serde::{Deserialize, Serialize};

use crate::{
    gp::{NAx2d, NTrsf2d, NVec2d, NXY},
    nature_errors::NErrors,
};

// Trait to define the behavior of a unit vector (direction) in 2D space
pub trait Dir2d {
    fn new(x: f64, y: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn from_vec2d(vec: &NVec2d) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn from_xy(xy: &NXY) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn set_coord(&mut self, index: usize, value: f64) -> Result<(), NErrors>;
    fn set_coords(&mut self, x: f64, y: f64) -> Result<(), NErrors>;
    fn set_x(&mut self, x: f64) -> Result<(), NErrors>;
    fn set_y(&mut self, y: f64) -> Result<(), NErrors>;
    fn set_xy(&mut self, xy: &NXY) -> Result<(), NErrors>;
    fn coord(&self, index: usize) -> Result<f64, NErrors>;
    fn coords(&self) -> (f64, f64);
    fn x(&self) -> f64;
    fn y(&self) -> f64;
    fn xy(&self) -> &NXY;
    fn is_equal(&self, other: &Self, angular_tolerance: f64) -> bool;
    fn is_normal(&self, other: &Self, angular_tolerance: f64) -> bool;
    fn is_opposite(&self, other: &Self, angular_tolerance: f64) -> bool;
    fn is_parallel(&self, other: &Self, angular_tolerance: f64) -> bool;
    fn angle(&self, other: &Self) -> f64;
    fn dot(&self, other: &Self) -> f64;
    fn crossed(&self, right: &Self) -> f64;
    fn reverse(&mut self);
    fn reversed(&self) -> Self;
    fn mirror_dir(&mut self, v: &Self);
    fn mirrored_dir(&self, v: &Self) -> Self;
    fn mirror_ax2d(&mut self, a: &NAx2d);
    fn mirrored_ax2d(&self, a: &NAx2d) -> Self;
    fn rotate(&mut self, angle: f64);
    fn rotated(&self, angle: f64) -> Self;
    fn transform(&mut self, trsf: &NTrsf2d);
    fn transformed(&self, trsf: &NTrsf2d) -> Self;
}

// Struct representing a unit vector (direction) in 2D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NDir2d {
    coord: NXY,
}

impl Dir2d for NDir2d {
    fn new(x: f64, y: f64) -> Result<Self, NErrors> {
        let norm = (x * x + y * y).sqrt();
        if norm <= crate::gp::NGP::resolution() {
            return Err(NErrors::ZeroNorm);
        }
        Ok(NDir2d {
            coord: NXY::new(x / norm, y / norm),
        })
    }

    fn from_vec2d(vec: &NVec2d) -> Result<Self, NErrors> {
        let xy = vec.xy();
        Self::from_xy(&xy)
    }

    fn from_xy(xy: &NXY) -> Result<Self, NErrors> {
        let norm = xy.modulus();
        if norm <= crate::gp::NGP::resolution() {
            return Err(NErrors::ZeroNorm);
        }
        Ok(NDir2d {
            coord: NXY::new(xy.x() / norm, xy.y() / norm),
        })
    }

    fn set_coord(&mut self, index: usize, value: f64) -> Result<(), NErrors> {
        let mut coords = [self.coord.x(), self.coord.y()];
        if index < 1 || index > 2 {
            return Err(NErrors::OutOfRange);
        }
        coords[index - 1] = value;
        let norm = (coords[0] * coords[0] + coords[1] * coords[1]).sqrt();
        if norm <= crate::gp::NGP::resolution() {
            return Err(NErrors::ZeroNorm);
        }
        self.coord = NXY::new(coords[0] / norm, coords[1] / norm);
        Ok(())
    }

    fn set_coords(&mut self, x: f64, y: f64) -> Result<(), NErrors> {
        let norm = (x * x + y * y).sqrt();
        if norm <= crate::gp::NGP::resolution() {
            return Err(NErrors::ZeroNorm);
        }
        self.coord = NXY::new(x / norm, y / norm);
        Ok(())
    }

    fn set_x(&mut self, x: f64) -> Result<(), NErrors> {
        self.set_coord(1, x)
    }

    fn set_y(&mut self, y: f64) -> Result<(), NErrors> {
        self.set_coord(2, y)
    }

    fn set_xy(&mut self, xy: &NXY) -> Result<(), NErrors> {
        self.set_coords(xy.x(), xy.y())
    }

    fn coord(&self, index: usize) -> Result<f64, NErrors> {
        match index {
            1 => Ok(self.coord.x()),
            2 => Ok(self.coord.y()),
            _ => Err(NErrors::OutOfRange),
        }
    }

    fn coords(&self) -> (f64, f64) {
        (self.coord.x(), self.coord.y())
    }

    fn x(&self) -> f64 {
        self.coord.x()
    }

    fn y(&self) -> f64 {
        self.coord.y()
    }

    fn xy(&self) -> &NXY {
        &self.coord
    }

    fn is_equal(&self, other: &Self, angular_tolerance: f64) -> bool {
        self.angle(other).abs() <= angular_tolerance
    }

    fn is_normal(&self, other: &Self, angular_tolerance: f64) -> bool {
        let ang = self.angle(other).abs();
        (PI / 2.0 - ang).abs() <= angular_tolerance
    }

    fn is_opposite(&self, other: &Self, angular_tolerance: f64) -> bool {
        let ang = self.angle(other).abs();
        (PI - ang) <= angular_tolerance
    }

    fn is_parallel(&self, other: &Self, angular_tolerance: f64) -> bool {
        let ang = self.angle(other).abs();
        ang <= angular_tolerance || (PI - ang) <= angular_tolerance
    }

    fn angle(&self, other: &Self) -> f64 {
        let cosinus = self.coord.dot(&other.coord);
        let sinus = self.coord.crossed(&other.coord);
        if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
            if sinus > 0.0 {
                cosinus.acos()
            } else {
                -cosinus.acos()
            }
        } else if cosinus > 0.0 {
            sinus.asin()
        } else if sinus > 0.0 {
            PI - sinus.asin()
        } else {
            -PI - sinus.asin()
        }
    }

    fn dot(&self, other: &Self) -> f64 {
        self.coord.dot(&other.coord)
    }

    fn crossed(&self, right: &Self) -> f64 {
        self.coord.crossed(&right.coord)
    }

    fn reverse(&mut self) {
        self.coord.reverse();
    }

    fn reversed(&self) -> Self {
        let mut result = self.clone();
        result.reverse();
        result
    }

    fn mirror_dir(&mut self, v: &Self) {
        let a = v.coord.x();
        let b = v.coord.y();
        let x = self.coord.x();
        let y = self.coord.y();
        let m1 = 2.0 * a * b;
        let xx = ((2.0 * a * a) - 1.0) * x + m1 * y;
        let yy = m1 * x + ((2.0 * b * b) - 1.0) * y;
        self.coord = NXY::new(xx, yy);
        let norm = self.coord.modulus();
        if norm <= crate::gp::NGP::resolution() {
            self.coord = NXY::new(1.0, 0.0); // Fallback to avoid zero norm
        } else {
            self.coord.divide(norm);
        }
    }

    fn mirrored_dir(&self, v: &Self) -> Self {
        let mut result = self.clone();
        result.mirror_dir(v);
        result
    }

    fn mirror_ax2d(&mut self, a: &NAx2d) {
        self.mirror_dir(&a.direction());
    }

    fn mirrored_ax2d(&self, a: &NAx2d) -> Self {
        let mut result = self.clone();
        result.mirror_ax2d(a);
        result
    }

    fn rotate(&mut self, angle: f64) {
        let trsf =
            NTrsf2d::new_rotation(&NPoint2d::new(0.0, 0.0), angle).expect("Invalid rotation");
        self.transform(&trsf);
    }

    fn rotated(&self, angle: f64) -> Self {
        let mut result = self.clone();
        result.rotate(angle);
        result
    }

    fn transform(&mut self, trsf: &NTrsf2d) {
        match trsf.form() {
            NTrsf2dForm::Identity | NTrsf2dForm::Translation => {}
            NTrsf2dForm::PntMirror => self.reverse(),
            NTrsf2dForm::Scale => {
                if trsf.scale_factor() < 0.0 {
                    self.reverse();
                }
            }
            _ => {
                self.coord.multiply(&trsf.hvectorial_part());
                let norm = self.coord.modulus();
                if norm <= crate::gp::NGP::resolution() {
                    self.coord = NXY::new(1.0, 0.0); // Fallback to avoid zero norm
                } else {
                    self.coord.divide(norm);
                }
                if trsf.scale_factor() < 0.0 {
                    self.reverse();
                }
            }
        }
    }

    fn transformed(&self, trsf: &NTrsf2d) -> Self {
        let mut result = self.clone();
        result.transform(trsf);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn dir2d(x: f64, y: f64) -> NDir2d {
        NDir2d::new(x, y).expect("Invalid direction")
    }

    #[test]
    fn test_new() {
        let dir = NDir2d::new(1.0, 0.0).unwrap();
        assert_eq!(dir.coords(), (1.0, 0.0));

        assert!(matches!(NDir2d::new(0.0, 0.0), Err(NErrors::ZeroNorm)));
    }

    #[test]
    fn test_from_vec2d() {
        let vec = NVec2d::new(2.0, 0.0);
        let dir = NDir2d::from_vec2d(&vec).unwrap();
        assert_eq!(dir.coords(), (1.0, 0.0));

        let vec_zero = NVec2d::new(0.0, 0.0);
        assert!(matches!(
            NDir2d::from_vec2d(&vec_zero),
            Err(NErrors::ZeroNorm)
        ));
    }

    #[test]
    fn test_from_xy() {
        let xy = NXY::new(3.0, 4.0);
        let dir = NDir2d::from_xy(&xy).unwrap();
        assert!((dir.x() - 0.6).abs() < 1e-10);
        assert!((dir.y() - 0.8).abs() < 1e-10);

        let xy_zero = NXY::new(0.0, 0.0);
        assert!(matches!(NDir2d::from_xy(&xy_zero), Err(NErrors::ZeroNorm)));
    }

    #[test]
    fn test_set_coord() {
        let mut dir = dir2d(1.0, 0.0);
        dir.set_coord(2, 1.0).unwrap();
        assert!((dir.x() - 1.0 / SQRT_2).abs() < 1e-10);
        assert!((dir.y() - 1.0 / SQRT_2).abs() < 1e-10);

        assert!(matches!(dir.set_coord(3, 1.0), Err(NErrors::OutOfRange)));
        assert!(matches!(dir.set_coord(1, 0.0), Err(NErrors::ZeroNorm)));
    }

    #[test]
    fn test_set_coords() {
        let mut dir = dir2d(1.0, 0.0);
        dir.set_coords(0.0, 1.0).unwrap();
        assert_eq!(dir.coords(), (0.0, 1.0));

        assert!(matches!(dir.set_coords(0.0, 0.0), Err(NErrors::ZeroNorm)));
    }

    #[test]
    fn test_set_xy() {
        let mut dir = dir2d(1.0, 0.0);
        let xy = NXY::new(3.0, 4.0);
        dir.set_xy(&xy).unwrap();
        assert!((dir.x() - 0.6).abs() < 1e-10);
        assert!((dir.y() - 0.8).abs() < 1e-10);
    }

    #[test]
    fn test_getters() {
        let dir = dir2d(3.0, 4.0);
        assert!((dir.x() - 0.6).abs() < 1e-10);
        assert!((dir.y() - 0.8).abs() < 1e-10);
        assert_eq!(dir.coord(1).unwrap(), dir.x());
        assert_eq!(dir.coords(), (dir.x(), dir.y()));
        assert_eq!(dir.xy(), &NXY::new(dir.x(), dir.y()));
    }

    #[test]
    fn test_relations() {
        let d1 = dir2d(1.0, 0.0);
        let d2 = dir2d(0.0, 1.0);
        let d3 = dir2d(-1.0, 0.0);
        let tol = 1e-6;

        assert!(d1.is_normal(&d2, tol));
        assert!(d1.is_opposite(&d3, tol));
        assert!(d1.is_parallel(&d3, tol));
        assert!(!d1.is_equal(&d2, tol));
        assert!(d1.is_equal(&d1, tol));
    }

    #[test]
    fn test_angle() {
        let d1 = dir2d(1.0, 0.0);
        let d2 = dir2d(0.0, 1.0);
        let d3 = dir2d(-1.0, 0.0);
        let d4 = dir2d(0.0, -1.0);
        assert!((d1.angle(&d2) - PI / 2.0).abs() < 1e-10);
        assert!((d1.angle(&d3) - PI).abs() < 1e-10);
        assert!((d1.angle(&d4) - (-PI / 2.0)).abs() < 1e-10);
    }

    #[test]
    fn test_dot_crossed() {
        let d1 = dir2d(1.0, 0.0);
        let d2 = dir2d(0.0, 1.0);
        let d3 = dir2d(3.0, 4.0);
        assert_eq!(d1.dot(&d2), 0.0);
        assert!((d1.dot(&d3) - 0.6).abs() < 1e-10);
        assert_eq!(d1.crossed(&d2), 1.0);
        assert_eq!(d2.crossed(&d1), -1.0);
    }

    #[test]
    fn test_reverse() {
        let mut d1 = dir2d(1.0, 0.0);
        d1.reverse();
        assert_eq!(d1.coords(), (-1.0, 0.0));

        let d2 = d1.reversed();
        assert_eq!(d2.coords(), (1.0, 0.0));
    }

    #[test]
    fn test_mirror_dir() {
        let mut d1 = dir2d(1.0, 0.0);
        let v = dir2d(0.0, 1.0);
        d1.mirror_dir(&v);
        assert_eq!(d1.coords(), (1.0, 0.0));

        let d2 = dir2d(0.0, 1.0);
        let mirrored = d2.mirrored_dir(&v);
        assert_eq!(mirrored.coords(), (0.0, -1.0));
    }

    #[test]
    fn test_mirror_ax2d() {
        let mut d1 = dir2d(0.0, 1.0);
        let ax2d = NAx2d::new(NPoint2d::new(0.0, 0.0), dir2d(0.0, 1.0));
        d1.mirror_ax2d(&ax2d);
        assert_eq!(d1.coords(), (0.0, -1.0));

        let mirrored = d1.mirrored_ax2d(&ax2d);
        assert_eq!(mirrored.coords(), (0.0, 1.0));
    }

    #[test]
    fn test_rotate() {
        let mut d1 = dir2d(1.0, 0.0);
        d1.rotate(PI / 2.0);
        assert!((d1.x() - 0.0).abs() < 1e-10);
        assert!((d1.y() - 1.0).abs() < 1e-10);

        let rotated = d1.rotated(-PI / 2.0);
        assert_eq!(rotated.coords(), (1.0, 0.0));
    }

    #[test]
    fn test_transform() {
        let mut d1 = dir2d(1.0, 0.0);
        let trsf = NTrsf2d::new_scale(&NPoint2d::new(0.0, 0.0), -2.0).unwrap();
        d1.transform(&trsf);
        assert_eq!(d1.coords(), (-1.0, 0.0));

        let trsf = NTrsf2d::new_rotation(&NPoint2d::new(0.0, 0.0), PI / 2.0).unwrap();
        let transformed = d1.transformed(&trsf);
        assert!((transformed.x() - 0.0).abs() < 1e-10);
        assert!((transformed.y() - 1.0).abs() < 1e-10);
    }
}
