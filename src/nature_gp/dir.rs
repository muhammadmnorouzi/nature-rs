use std::f64::consts::{PI, SQRT_2};

use serde::{Deserialize, Serialize};

use crate::{
    gp::{NAx1, NAx2, NTrsf, NVec, NXYZ},
    nature_errors::NErrors,
};

use super::
    xyz::{NXYZ, XYZ}
;

// Trait to define the behavior of a unit vector (direction) in 3D space
pub trait Dir {
    fn new(x: f64, y: f64, z: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn from_vec(vec: &NVec) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn from_xyz(xyz: &NXYZ) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn set_coord(&mut self, index: usize, value: f64) -> Result<(), NErrors>;
    fn set_coords(&mut self, x: f64, y: f64, z: f64) -> Result<(), NErrors>;
    fn set_x(&mut self, x: f64) -> Result<(), NErrors>;
    fn set_y(&mut self, y: f64) -> Result<(), NErrors>;
    fn set_z(&mut self, z: f64) -> Result<(), NErrors>;
    fn set_xyz(&mut self, xyz: &NXYZ) -> Result<(), NErrors>;
    fn coord(&self, index: usize) -> Result<f64, NErrors>;
    fn coords(&self) -> (f64, f64, f64);
    fn x(&self) -> f64;
    fn y(&self) -> f64;
    fn z(&self) -> f64;
    fn xyz(&self) -> &NXYZ;
    fn is_equal(&self, other: &Self, angular_tolerance: f64) -> bool;
    fn is_normal(&self, other: &Self, angular_tolerance: f64) -> bool;
    fn is_opposite(&self, other: &Self, angular_tolerance: f64) -> bool;
    fn is_parallel(&self, other: &Self, angular_tolerance: f64) -> bool;
    fn angle(&self, other: &Self) -> f64;
    fn angle_with_ref(&self, other: &Self, vref: &Self) -> Result<f64, NErrors>;
    fn cross(&mut self, right: &Self) -> Result<(), NErrors>;
    fn crossed(&self, right: &Self) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn cross_cross(&mut self, v1: &Self, v2: &Self) -> Result<(), NErrors>;
    fn cross_crossed(&self, v1: &Self, v2: &Self) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn dot(&self, other: &Self) -> f64;
    fn dot_cross(&self, v1: &Self, v2: &Self) -> f64;
    fn reverse(&mut self);
    fn reversed(&self) -> Self;
    fn mirror_dir(&mut self, v: &Self);
    fn mirrored_dir(&self, v: &Self) -> Self;
    fn mirror_ax1(&mut self, a1: &NAx1);
    fn mirrored_ax1(&self, a1: &NAx1) -> Self;
    fn mirror_ax2(&mut self, a2: &NAx2);
    fn mirrored_ax2(&self, a2: &NAx2) -> Self;
    fn rotate(&mut self, axis: &NAx1, angle: f64);
    fn rotated(&self, axis: &NAx1, angle: f64) -> Self;
    fn transform(&mut self, trsf: &NTrsf);
    fn transformed(&self, trsf: &NTrsf) -> Self;
}

// Struct representing a unit vector (direction) in 3D space
#[derive(Clone, PartialEq, Default, Debug, Serialize, Deserialize)]
pub struct NDir {
    coord: NXYZ,
}

impl Dir for NDir {
    fn new(x: f64, y: f64, z: f64) -> Result<Self, NErrors> {
        let mut xyz = NXYZ::new(x, y, z);

        if let Err(error) = xyz.normalize() {
            Err(error)
        } else {
            Ok(NDir { coord: xyz })
        }
    }

    fn from_vec(vec: &NVec) -> Result<Self, NErrors> {
        let xyz = vec.xyz();
        Self::from_xyz(&xyz)
    }

    fn from_xyz(xyz: &NXYZ) -> Result<Self, NErrors> {
        let norm = xyz.modulus();
        if norm <= crate::gp::NGP::resolution() {
            return Err(NErrors::ZeroNorm);
        }
        Ok(NDir {
            coord: NXYZ::new(xyz.x() / norm, xyz.y() / norm, xyz.z() / norm),
        })
    }

    fn set_coord(&mut self, index: usize, value: f64) -> Result<(), NErrors> {
        let mut coords = [self.coord.x(), self.coord.y(), self.coord.z()];
        if index < 1 || index > 3 {
            return Err(NErrors::OutOfRange);
        }
        coords[index - 1] = value;
        let norm = (coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2]).sqrt();
        if norm <= crate::gp::NGP::resolution() {
            return Err(NErrors::ZeroNorm);
        }
        self.coord = NXYZ::new(coords[0] / norm, coords[1] / norm, coords[2] / norm);
        Ok(())
    }

    fn set_coords(&mut self, x: f64, y: f64, z: f64) -> Result<(), NErrors> {
        let norm = (x * x + y * y + z * z).sqrt();
        if norm <= crate::gp::NGP::resolution() {
            return Err(NErrors::ZeroNorm);
        }
        self.coord = NXYZ::new(x / norm, y / norm, z / norm);
        Ok(())
    }

    fn set_x(&mut self, x: f64) -> Result<(), NErrors> {
        self.set_coord(1, x)
    }

    fn set_y(&mut self, y: f64) -> Result<(), NErrors> {
        self.set_coord(2, y)
    }

    fn set_z(&mut self, z: f64) -> Result<(), NErrors> {
        self.set_coord(3, z)
    }

    fn set_xyz(&mut self, xyz: &NXYZ) -> Result<(), NErrors> {
        self.set_coords(xyz.x(), xyz.y(), xyz.z())
    }

    fn coord(&self, index: usize) -> Result<f64, NErrors> {
        match index {
            1 => Ok(self.coord.x()),
            2 => Ok(self.coord.y()),
            3 => Ok(self.coord.z()),
            _ => Err(NErrors::OutOfRange),
        }
    }

    fn coords(&self) -> (f64, f64, f64) {
        (self.coord.x(), self.coord.y(), self.coord.z())
    }

    fn x(&self) -> f64 {
        self.coord.x()
    }

    fn y(&self) -> f64 {
        self.coord.y()
    }

    fn z(&self) -> f64 {
        self.coord.z()
    }

    fn xyz(&self) -> &NXYZ {
        &self.coord
    }

    fn is_equal(&self, other: &Self, angular_tolerance: f64) -> bool {
        self.angle(other) <= angular_tolerance
    }

    fn is_normal(&self, other: &Self, angular_tolerance: f64) -> bool {
        (PI / 2.0 - self.angle(other)).abs() <= angular_tolerance
    }

    fn is_opposite(&self, other: &Self, angular_tolerance: f64) -> bool {
        (PI - self.angle(other)) <= angular_tolerance
    }

    fn is_parallel(&self, other: &Self, angular_tolerance: f64) -> bool {
        let ang = self.angle(other);
        ang <= angular_tolerance || (PI - ang) <= angular_tolerance
    }

    fn angle(&self, other: &Self) -> f64 {
        let cosinus = self.coord.dot(&other.coord);
        if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
            cosinus.acos()
        } else {
            let sinus = self.coord.crossed(&other.coord).modulus();
            if cosinus < 0.0 {
                PI - sinus.asin()
            } else {
                sinus.asin()
            }
        }
    }

    fn angle_with_ref(&self, other: &Self, vref: &Self) -> Result<f64, NErrors> {
        let xyz = self.coord.crossed(&other.coord);
        let cosinus = self.coord.dot(&other.coord);
        let sinus = xyz.modulus();
        let ang = if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
            cosinus.acos()
        } else if cosinus < 0.0 {
            PI - sinus.asin()
        } else {
            sinus.asin()
        };
        if sinus <= crate::gp::NGP::resolution() {
            return Err(NErrors::ParallelVectors);
        }
        if xyz.dot(&vref.coord) >= 0.0 {
            Ok(ang)
        } else {
            Ok(-ang)
        }
    }

    fn cross(&mut self, right: &Self) -> Result<(), NErrors> {
        self.coord = self.coord.crossed(&right.coord);
        let norm = self.coord.modulus();
        if norm <= crate::gp::NGP::resolution() {
            return Err(NErrors::ZeroNorm);
        }
        self.coord.divide(norm);
        Ok(())
    }

    fn crossed(&self, right: &Self) -> Result<Self, NErrors> {
        let mut result = self.clone();
        result.cross(right)?;
        Ok(result)
    }

    fn cross_cross(&mut self, v1: &Self, v2: &Self) -> Result<(), NErrors> {
        self.coord = self.coord.cross_crossed(&v1.coord, &v2.coord);
        let norm = self.coord.modulus();
        if norm <= crate::gp::NGP::resolution() {
            return Err(NErrors::ZeroNorm);
        }
        self.coord.divide(norm);
        Ok(())
    }

    fn cross_crossed(&self, v1: &Self, v2: &Self) -> Result<Self, NErrors> {
        let mut result = self.clone();
        result.cross_cross(v1, v2)?;
        Ok(result)
    }

    fn dot(&self, other: &Self) -> f64 {
        self.coord.dot(&other.coord)
    }

    fn dot_cross(&self, v1: &Self, v2: &Self) -> f64 {
        self.coord.dot(&v1.coord.crossed(&v2.coord))
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
        let c = v.coord.z();
        let x = self.coord.x();
        let y = self.coord.y();
        let z = self.coord.z();
        let m1 = 2.0 * a * b;
        let m2 = 2.0 * a * c;
        let m3 = 2.0 * b * c;
        let xx = ((2.0 * a * a) - 1.0) * x + m1 * y + m2 * z;
        let yy = m1 * x + ((2.0 * b * b) - 1.0) * y + m3 * z;
        let zz = m2 * x + m3 * y + ((2.0 * c * c) - 1.0) * z;
        self.coord = NXYZ::new(xx, yy, zz);
        let norm = self.coord.modulus();
        if norm <= crate::gp::NGP::resolution() {
            self.coord = NXYZ::new(1.0, 0.0, 0.0); // Fallback to avoid zero norm
        } else {
            self.coord.divide(norm);
        }
    }

    fn mirrored_dir(&self, v: &Self) -> Self {
        let mut result = self.clone();
        result.mirror_dir(v);
        result
    }

    fn mirror_ax1(&mut self, a1: &NAx1) {
        self.mirror_dir(&a1.direction());
    }

    fn mirrored_ax1(&self, a1: &NAx1) -> Self {
        let mut result = self.clone();
        result.mirror_ax1(a1);
        result
    }

    fn mirror_ax2(&mut self, a2: &NAx2) {
        self.mirror_dir(&a2.direction());
        self.reverse();
    }

    fn mirrored_ax2(&self, a2: &NAx2) -> Self {
        let mut result = self.clone();
        result.mirror_ax2(a2);
        result
    }

    fn rotate(&mut self, axis: &NAx1, angle: f64) {
        let trsf = NTrsf::new_rotation(axis, angle).expect("Invalid rotation");
        self.transform(&trsf);
    }

    fn rotated(&self, axis: &NAx1, angle: f64) -> Self {
        let mut result = self.clone();
        result.rotate(axis, angle);
        result
    }

    fn transform(&mut self, trsf: &NTrsf) {
        match trsf.form() {
            NTrsfForm::Identity | NTrsfForm::Translation => {}
            NTrsfForm::PntMirror => self.reverse(),
            NTrsfForm::Scale => {
                if trsf.scale_factor() < 0.0 {
                    self.reverse();
                }
            }
            _ => {
                self.coord.multiply(&trsf.hvectorial_part());
                let norm = self.coord.modulus();
                if norm <= crate::gp::NGP::resolution() {
                    self.coord = NXYZ::new(1.0, 0.0, 0.0); // Fallback to avoid zero norm
                } else {
                    self.coord.divide(norm);
                }
                if trsf.scale_factor() < 0.0 {
                    self.reverse();
                }
            }
        }
    }

    fn transformed(&self, trsf: &NTrsf) -> Self {
        let mut result = self.clone();
        result.transform(trsf);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn dir(x: f64, y: f64, z: f64) -> NDir {
        NDir::new(x, y, z).expect("Invalid direction")
    }

    #[test]
    fn test_new() {
        let dir = NDir::new(1.0, 0.0, 0.0).unwrap();
        assert_eq!(dir.coords(), (1.0, 0.0, 0.0));

        assert!(matches!(NDir::new(0.0, 0.0, 0.0), Err(NErrors::ZeroNorm)));
    }

    fn test_from_vec() {
        let vec = NVec::new(2.0, 0.0, 0.0);
        let dir = NDir::from_vec(&vec).unwrap();
        assert_eq!(dir.coords(), (1.0, 0.0, 0.0));

        let vec_zero = NVec::new(0.0, 0.0, 0.0);
        assert!(matches!(NDir::from_vec(&vec_zero), Err(NErrors::ZeroNorm)));
    }

    #[test]
    fn test_from_xyz() {
        let xyz = NXYZ::new(3.0, 4.0, 0.0);
        let dir = NDir::from_xyz(&xyz).unwrap();
        assert!((dir.x() - 0.6).abs() < 1e-10);
        assert!((dir.y() - 0.8).abs() < 1e-10);
        assert_eq!(dir.z(), 0.0);

        let xyz_zero = NXYZ::new(0.0, 0.0, 0.0);
        assert!(matches!(NDir::from_xyz(&xyz_zero), Err(NErrors::ZeroNorm)));
    }

    #[test]
    fn test_set_coord() {
        let mut dir = dir(1.0, 0.0, 0.0);
        dir.set_coord(2, 1.0).unwrap();
        assert!((dir.x() - 1.0 / SQRT_2).abs() < 1e-10);
        assert!((dir.y() - 1.0 / SQRT_2).abs() < 1e-10);
        assert_eq!(dir.z(), 0.0);

        assert!(matches!(dir.set_coord(4, 1.0), Err(NErrors::OutOfRange)));
        assert!(matches!(dir.set_coord(1, 0.0), Err(NErrors::ZeroNorm)));
    }

    #[test]
    fn test_set_coords() {
        let mut dir = dir(1.0, 0.0, 0.0);
        dir.set_coords(0.0, 1.0, 0.0).unwrap();
        assert_eq!(dir.coords(), (0.0, 1.0, 0.0));

        assert!(matches!(
            dir.set_coords(0.0, 0.0, 0.0),
            Err(NErrors::ZeroNorm)
        ));
    }

    #[test]
    fn test_set_xyz() {
        let mut dir = dir(1.0, 0.0, 0.0);
        let xyz = NXYZ::new(0.0, 3.0, 4.0);
        dir.set_xyz(&xyz).unwrap();
        assert!((dir.x() - 0.0).abs() < 1e-10);
        assert!((dir.y() - 0.6).abs() < 1e-10);
        assert!((dir.z() - 0.8).abs() < 1e-10);
    }

    #[test]
    fn test_getters() {
        let dir = dir(3.0, 4.0, 0.0);
        assert!((dir.x() - 0.6).abs() < 1e-10);
        assert!((dir.y() - 0.8).abs() < 1e-10);
        assert_eq!(dir.z(), 0.0);
        assert_eq!(dir.coord(1).unwrap(), dir.x());
        assert_eq!(dir.coords(), (dir.x(), dir.y(), dir.z()));
        assert_eq!(dir.xyz(), &NXYZ::new(dir.x(), dir.y(), dir.z()));
    }

    #[test]
    fn test_relations() {
        let d1 = dir(1.0, 0.0, 0.0);
        let d2 = dir(0.0, 1.0, 0.0);
        let d3 = dir(-1.0, 0.0, 0.0);
        let tol = 1e-6;

        assert!(d1.is_normal(&d2, tol));
        assert!(d1.is_opposite(&d3, tol));
        assert!(d1.is_parallel(&d3, tol));
        assert!(!d1.is_equal(&d2, tol));
        assert!(d1.is_equal(&d1, tol));
    }

    #[test]
    fn test_angle() {
        let d1 = dir(1.0, 0.0, 0.0);
        let d2 = dir(0.0, 1.0, 0.0);
        let d3 = dir(-1.0, 0.0, 0.0);
        assert!((d1.angle(&d2) - PI / 2.0).abs() < 1e-10);
        assert!((d1.angle(&d3) - PI).abs() < 1e-10);

        let vref = dir(0.0, 0.0, 1.0);
        assert!((d1.angle_with_ref(&d2, &vref).unwrap() - PI / 2.0).abs() < 1e-10);
        assert!((d2.angle_with_ref(&d1, &vref).unwrap() - (-PI / 2.0)).abs() < 1e-10);
        assert!(matches!(
            d1.angle_with_ref(&d1, &vref),
            Err(NErrors::ParallelVectors)
        ));
    }

    #[test]
    fn test_cross() {
        let mut d1 = dir(1.0, 0.0, 0.0);
        let d2 = dir(0.0, 1.0, 0.0);
        d1.cross(&d2).unwrap();
        assert_eq!(d1.coords(), (0.0, 0.0, 1.0));

        let d3 = dir(1.0, 0.0, 0.0);
        assert!(matches!(d1.cross(&d3), Err(NErrors::ZeroNorm)));

        let crossed = d3.crossed(&d2).unwrap();
        assert_eq!(crossed.coords(), (0.0, 0.0, 1.0));
    }

    #[test]
    fn test_cross_cross() {
        let mut d1 = dir(1.0, 0.0, 0.0);
        let v1 = dir(0.0, 1.0, 0.0);
        let v2 = dir(0.0, 0.0, 1.0);
        d1.cross_cross(&v1, &v2).unwrap();
        assert!((d1.x() - 1.0).abs() < 1e-10);
        assert_eq!(d1.y(), 0.0);
        assert_eq!(d1.z(), 0.0);

        let crossed = d1.cross_crossed(&v1, &v2).unwrap();
        assert_eq!(crossed.coords(), d1.coords());
    }

    #[test]
    fn test_dot() {
        let d1 = dir(1.0, 0.0, 0.0);
        let d2 = dir(0.0, 1.0, 0.0);
        let d3 = dir(3.0, 4.0, 0.0);
        assert_eq!(d1.dot(&d2), 0.0);
        assert!((d1.dot(&d3) - 0.6).abs() < 1e-10);
        assert_eq!(d1.dot_cross(&d2, &dir(0.0, 0.0, 1.0)), 1.0);
    }

    #[test]
    fn test_reverse() {
        let mut d1 = dir(1.0, 0.0, 0.0);
        d1.reverse();
        assert_eq!(d1.coords(), (-1.0, 0.0, 0.0));

        let d2 = d1.reversed();
        assert_eq!(d2.coords(), (1.0, 0.0, 0.0));
    }

    #[test]
    fn test_mirror_dir() {
        let mut d1 = dir(1.0, 0.0, 0.0);
        let v = dir(0.0, 1.0, 0.0);
        d1.mirror_dir(&v);
        assert_eq!(d1.coords(), (1.0, 0.0, 0.0));

        let d2 = dir(0.0, 1.0, 0.0);
        let mirrored = d2.mirrored_dir(&v);
        assert_eq!(mirrored.coords(), (0.0, -1.0, 0.0));
    }

    #[test]
    fn test_mirror_ax1() {
        let mut d1 = dir(0.0, 1.0, 0.0);
        let ax1 = NAx1::new(NPnt::new(0.0, 0.0, 0.0), dir(0.0, 1.0, 0.0));
        d1.mirror_ax1(&ax1);
        assert_eq!(d1.coords(), (0.0, -1.0, 0.0));

        let mirrored = d1.mirrored_ax1(&ax1);
        assert_eq!(mirrored.coords(), (0.0, 1.0, 0.0));
    }

    #[test]
    fn test_mirror_ax2() {
        let mut d1 = dir(0.0, 1.0, 0.0);
        let ax2 = NAx2::new(
            NPnt::new(0.0, 0.0, 0.0),
            dir(0.0, 0.0, 1.0),
            dir(1.0, 0.0, 0.0),
        )
        .unwrap();
        d1.mirror_ax2(&ax2);
        assert_eq!(d1.coords(), (0.0, -1.0, 0.0));

        let mirrored = d1.mirrored_ax2(&ax2);
        assert_eq!(mirrored.coords(), (0.0, 1.0, 0.0));
    }

    #[test]
    fn test_rotate() {
        let mut d1 = dir(1.0, 0.0, 0.0);
        let axis = NAx1::new(NPnt::new(0.0, 0.0, 0.0), dir(0.0, 0.0, 1.0));
        d1.rotate(&axis, PI / 2.0);
        assert!((d1.x() - 0.0).abs() < 1e-10);
        assert!((d1.y() - 1.0).abs() < 1e-10);
        assert_eq!(d1.z(), 0.0);

        let rotated = d1.rotated(&axis, -PI / 2.0);
        assert_eq!(rotated.coords(), (1.0, 0.0, 0.0));
    }

    #[test]
    fn test_transform() {
        let mut d1 = dir(1.0, 0.0, 0.0);
        let trsf = NTrsf::new_scale(&NPnt::new(0.0, 0.0, 0.0), -2.0).unwrap();
        d1.transform(&trsf);
        assert_eq!(d1.coords(), (-1.0, 0.0, 0.0));

        let trsf = NTrsf::new_rotation(
            &NAx1::new(NPnt::new(0.0, 0.0, 0.0), dir(0.0, 0.0, 1.0)),
            PI / 2.0,
        )
        .unwrap();
        let transformed = d1.transformed(&trsf);
        assert!((transformed.x() - 0.0).abs() < 1e-10);
        assert!((transformed.y() - 1.0).abs() < 1e-10);
        assert_eq!(transformed.z(), 0.0);
    }
}
