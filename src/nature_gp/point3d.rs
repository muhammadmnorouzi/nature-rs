use super::prelude::*;
use crate::nature_common::prelude::*;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

// Trait to define the behavior of a 3D Cartesian point
pub trait Point3d
where
    Self: Default + Clone + PartialEq + Eq + From<f64> + Debug,
{
    fn new(x: f64, y: f64, z: f64) -> Self;
    fn change_coord_xyz(&mut self) -> &mut NXYZ;
    fn change_coord(&mut self, index: usize) -> Result<&mut f64, NErrors>;
    fn set_coord_by_index(&mut self, index: usize, scalar: f64) -> Result<(), NErrors>;
    fn set_coords(&mut self, x: f64, y: f64, z: f64);
    fn set_x(&mut self, x: f64);
    fn set_y(&mut self, y: f64);
    fn set_z(&mut self, z: f64);
    fn set_xyz(&mut self, coord: NXYZ);
    fn coord_by_index(&self, index: usize) -> Result<f64, NErrors>;
    fn coords(&self) -> (f64, f64, f64);
    fn x(&self) -> f64;
    fn y(&self) -> f64;
    fn z(&self) -> f64;
    fn xyz(&self) -> NXYZ;
    fn bary_center(&mut self, alpha: f64, other: &Self, beta: f64);
    fn is_equal(&self, other: &Self, linear_tolerance: f64) -> bool;
    fn distance(&self, other: &Self) -> f64;
    fn square_distance(&self, other: &Self) -> f64;
    fn mirror_pnt(&mut self, p: &Self);
    fn mirrored_pnt(&self, p: &Self) -> Self
    where
        Self: Sized;
    fn mirror_ax1(&mut self, a1: &NAx1);
    fn mirrored_ax1(&self, a1: &NAx1) -> Self
    where
        Self: Sized;
    fn mirror_ax2(&mut self, a2: &NAx2);
    fn mirrored_ax2(&self, a2: &NAx2) -> Self
    where
        Self: Sized;
    fn rotate(&mut self, a1: &NAx1, ang: f64);
    fn rotated(&self, a1: &NAx1, ang: f64) -> Self
    where
        Self: Sized;
    fn scale(&mut self, p: &Self, s: f64);
    fn scaled(&self, p: &Self, s: f64) -> Self
    where
        Self: Sized;
    fn transform(&mut self, t: &NTrsf);
    fn transformed(&self, t: &NTrsf) -> Self
    where
        Self: Sized;
    fn translate_vec(&mut self, v: &NVec);
    fn translated_vec(&self, v: &NVec) -> Self
    where
        Self: Sized;
    fn translate_point3d(&mut self, p1: &Self, p2: &Self);
    fn translated_point3d(&self, p1: &Self, p2: &Self) -> Self
    where
        Self: Sized;
}

#[derive(Clone, Default, PartialEq, Debug, Serialize, Deserialize)]
pub struct NPoint3d {
    coord: NXYZ,
}

impl Eq for NPoint3d {}

impl From<f64> for NPoint3d {
    fn from(value: f64) -> Self {
        Self::new(value, value, value)
    }
}

impl From<(f64, f64, f64)> for NPoint3d {
    fn from(value: (f64, f64, f64)) -> Self {
        Self::new(value.0, value.1, value.2)
    }
}

impl From<NXYZ> for NPoint3d {
    fn from(value: NXYZ) -> Self {
        Self { coord: value }
    }
}

impl Point3d for NPoint3d {
    fn set_coord_by_index(&mut self, index: usize, value: f64) -> Result<(), NErrors> {
        self.coord.set_coord_by_index(index, value)
    }

    /// Assigns the given values to the coordinates.
    fn set_coords(&mut self, x: f64, y: f64, z: f64) {
        self.coord.set_coords(x, y, z);
    }

    /// Assigns the given value to the X coordinate.
    fn set_x(&mut self, x: f64) {
        self.coord.set_x(x);
    }

    /// Assigns the given value to the Y coordinate.
    fn set_y(&mut self, y: f64) {
        self.coord.set_y(y);
    }

    /// Assigns the given value to the Z coordinate.
    fn set_z(&mut self, z: f64) {
        self.coord.set_z(z);
    }

    /// Assigns the coordinates of the given XYZ object.
    fn set_xyz(&mut self, coord: NXYZ) {
        self.coord = coord;
    }

    /// Returns the coordinate at the given index (1=X, 2=Y, 3=Z).
    fn coord_by_index(&self, index: usize) -> Result<f64, NErrors> {
        self.coord.coord_by_index(index)
    }

    /// Returns the three coordinates as a tuple.
    fn coords(&self) -> (f64, f64, f64) {
        (self.x(), self.y(), self.z())
    }

    /// Returns the X coordinate.
    fn x(&self) -> f64 {
        self.coord.x()
    }

    /// Returns the Y coordinate.
    fn y(&self) -> f64 {
        self.coord.y()
    }

    /// Returns the Z coordinate.
    fn z(&self) -> f64 {
        self.coord.z()
    }

    /// Returns the coordinates as an XYZ object.
    fn xyz(&self) -> NXYZ {
        self.coord.clone()
    }

    /// Assigns the result of (alpha*this + beta*other)/(alpha + beta) to this point.
    fn bary_center(&mut self, alpha: f64, other: &Self, beta: f64) {
        let mut new_coord = NXYZ::default();
        new_coord.set_linear_form22(alpha, &self.coord, beta, &other.coord);
        new_coord.divide(alpha + beta);
        self.coord = new_coord;
    }

    /// Returns true if the distance to the other point is within the linear tolerance.
    fn is_equal(&self, other: &Self, linear_tolerance: f64) -> bool {
        self.distance(other) <= linear_tolerance
    }

    /// Computes the Euclidean distance to another point.
    fn distance(&self, other: &Self) -> f64 {
        f64::sqrt(self.square_distance(other))
    }

    /// Computes the square of the Euclidean distance to another point.
    fn square_distance(&self, other: &Self) -> f64 {
        let dx = self.coord.x() - other.coord.x();
        let dy = self.coord.y() - other.coord.y();
        let dz = self.coord.z() - other.coord.z();
        (dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt()
    }

    /// Mirrors the point with respect to another point.
    fn mirror_pnt(&mut self, p: &Self) {
        let mut xyz = p.coord.clone();
        xyz.multiply(2.0);
        self.coord.reverse();
        self.coord.add(&xyz);
    }

    /// Returns the point mirrored with respect to another point.
    fn mirrored_pnt(&self, p: &Self) -> Self {
        let mut res = self.clone();
        res.mirror_pnt(p);
        res
    }

    /// Mirrors the point with respect to an axis.
    fn mirror_ax1(&mut self, a1: &NAx1) {
        let mut trsf = NTrsf::default();
        trsf.set_mirror_ax1(a1);
        trsf.transforms_xyz(&mut self.coord);
    }

    /// Returns the point mirrored with respect to an axis.
    fn mirrored_ax1(&self, a1: &NAx1) -> Self {
        let mut res = self.clone();
        res.mirror_ax1(a1);
        res
    }

    /// Mirrors the point with respect to a plane defined by Ax2.
    fn mirror_ax2(&mut self, a2: &NAx2) {
        let mut trsf = NTrsf::new();
        trsf.set_mirror_ax2(a2);
        trsf.transforms_xyz(&mut self.coord);
    }

    /// Returns the point mirrored with respect to a plane defined by Ax2.
    fn mirrored_ax2(&self, a2: &NAx2) -> Self {
        let mut res = self.clone();
        res.mirror_ax2(a2);
        res
    }

    /// Rotates the point around an axis by an angle.
    fn rotate(&mut self, a1: &NAx1, ang: f64) {
        let mut t = NTrsf::default();
        t.set_rotation_ax1(a1, ang);
        t.transforms_xyz(&mut self.coord);
    }

    /// Returns the point rotated around an axis by an angle.
    fn rotated(&self, a1: &NAx1, ang: f64) -> Self {
        let mut res = self.clone();
        res.rotate(a1, ang);
        res
    }

    /// Scales the point with respect to another point.
    fn scale(&mut self, p: &Self, scalar: f64) {
        let mut p_coord = p.coord.clone();
        p_coord.multiply(1.0 - scalar);
        self.coord.multiply(scalar);
        self.coord.add(&p_coord);
    }

    /// Returns the point scaled with respect to another point.
    fn scaled(&self, p: &Self, s: f64) -> Self {
        let mut res = self.clone();
        res.scale(p, s);
        res
    }

    /// Transforms the point with a transformation.
    fn transform(&mut self, t: &NTrsf) {
        match t.form() {
            NTrsfForm::Identity => {}
            NTrsfForm::Translation => {
                self.coord.add(&t.translation_part());
            }
            NTrsfForm::Scale => {
                self.coord.multiply(t.scale_factor());
                self.coord.add(&t.translation_part());
            }
            NTrsfForm::PointMirror => {
                self.coord.reverse();
                self.coord.add(&t.translation_part());
            }
            _ => {
                t.transforms_xyz(&mut self.coord);
            }
        }
    }

    /// Returns the point transformed with a transformation.
    fn transformed(&self, t: &NTrsf) -> Self {
        let mut res = self.clone();
        res.transform(t);
        res
    }

    /// Translates the point by a vector.
    fn translate_vec(&mut self, v: &NVec) {
        self.coord.add(&v.xyz());
    }

    /// Returns the point translated by a vector.
    fn translated_vec(&self, v: &NVec) -> Self {
        let mut res = self.clone();
        res.translate_vec(v);
        res
    }

    /// Translates the point from one point to another.
    fn translate_point3d(&mut self, p1: &Self, p2: &Self) {
        self.coord.add(&p2.coord);
        self.coord.subtract(&p1.coord);
    }

    /// Returns the point translated from one point to another.
    fn translated_point3d(&self, p1: &Self, p2: &Self) -> Self {
        let mut res = self.clone();
        res.translate_point3d(p1, p2);
        res
    }

    fn new(x: f64, y: f64, z: f64) -> Self {
        Self {
            coord: NXYZ::new(x, y, z),
        }
    }

    fn change_coord_xyz(&mut self) -> &mut NXYZ {
        &mut self.coord
    }

    fn change_coord(&mut self, index: usize) -> Result<&mut f64, NErrors> {
        self.coord.change_coord(index)
    }
}

// Implement std::hash::Hash for NPoint3d
impl std::hash::Hash for NPoint3d {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let x_bits = self.x().to_bits();
        let y_bits = self.y().to_bits();
        let z_bits = self.z().to_bits();
        // Simplified hashing based on C++ logic
        (x_bits / 23 + y_bits / 19 + z_bits / 17).hash(state);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_point() -> NPoint3d {
        NPoint3d::new(1.0, 2.0, 3.0)
    }

    #[test]
    fn test_new() {
        let p = NPoint3d::default();
        assert_eq!(p.coords(), (0.0, 0.0, 0.0));
    }

    #[test]
    fn test_new_with_coords() {
        let p = NPoint3d::new(1.0, 2.0, 3.0);
        assert_eq!(p.x(), 1.0);
        assert_eq!(p.y(), 2.0);
        assert_eq!(p.z(), 3.0);
    }

    #[test]
    fn test_set_coord() {
        let mut p = create_test_point();
        p.set_coord_by_index(1, 4.0).unwrap();
        assert_eq!(p.x(), 4.0);
        assert!(matches!(
            p.set_coord_by_index(4, 5.0),
            Err(NErrors::IndexOutOfRange)
        ));
    }

    #[test]
    fn test_coords() {
        let p = create_test_point();
        assert_eq!(p.coords(), (1.0, 2.0, 3.0));
        assert_eq!(p.coord_by_index(1).unwrap(), 1.0);
        assert!(matches!(p.coord_by_index(4), Err(NErrors::IndexOutOfRange)));
    }

    #[test]
    fn test_distance() {
        let p1 = NPoint3d::new(0.0, 0.0, 0.0);
        let p2 = NPoint3d::new(3.0, 4.0, 0.0);
        assert!((p1.distance(&p2) - 5.0).abs() < 1e-9);
        assert!((p1.square_distance(&p2) - 25.0).abs() < 1e-9);
    }

    #[test]
    fn test_is_equal() {
        let p1 = NPoint3d::new(1.0, 2.0, 3.0);
        let p2 = NPoint3d::new(1.0 + 1e-6, 2.0, 3.0);
        assert!(p1.is_equal(&p2, 1e-5));
        assert!(!p1.is_equal(&p2, 1e-7));
    }

    #[test]
    fn test_bary_center() {
        let mut p1 = NPoint3d::new(1.0, 0.0, 0.0);
        let p2 = NPoint3d::new(0.0, 1.0, 0.0);
        p1.bary_center(1.0, &p2, 1.0);
        let (x, y, z) = p1.coords();
        assert!((x - 0.5).abs() < 1e-9);
        assert!((y - 0.5).abs() < 1e-9);
        assert!((z - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_mirror_pnt() {
        let mut p = NPoint3d::new(1.0, 1.0, 1.0);
        let center = NPoint3d::new(0.0, 0.0, 0.0);
        p.mirror_pnt(&center);
        assert_eq!(p.coords(), (-1.0, -1.0, -1.0));
    }

    #[test]
    fn test_scale() {
        let mut p = NPoint3d::new(2.0, 2.0, 2.0);
        let center = NPoint3d::new(0.0, 0.0, 0.0);
        p.scale(&center, 2.0);
        assert_eq!(p.coords(), (4.0, 4.0, 4.0));
    }

    #[test]
    fn test_translate_vec() {
        let mut p = NPoint3d::new(1.0, 1.0, 1.0);
        let v = NVec::new_from_coords(1.0, 2.0, 3.0);
        p.translate_vec(&v);
        assert_eq!(p.coords(), (2.0, 3.0, 4.0));
    }
}
