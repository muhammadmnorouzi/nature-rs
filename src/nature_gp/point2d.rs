use serde::{Deserialize, Serialize};
use std::fmt::Debug;

use super::prelude::*;
use crate::nature_errors::NErrors;

/// Trait ponts, providing coordinate manipulation, geometric operations, and transformations.
pub trait Point2d
where
    Self: Default + Clone + PartialEq + Eq + From<f64> + Debug + Sized,
{
    fn new(x: f64, y: f64) -> Self;
    fn change_coord_xy(&mut self) -> &mut NXY;
    fn change_coord(&mut self, index: usize) -> Result<&mut f64, NErrors>;
    fn set_coord_by_index(&mut self, index: usize, scalar: f64) -> Result<(), NErrors>;
    fn set_coords(&mut self, x: f64, y: f64);
    fn set_x(&mut self, x: f64);
    fn set_y(&mut self, y: f64);
    fn set_xy(&mut self, coord: NXY);
    fn coord_by_index(&self, index: usize) -> Result<f64, NErrors>;
    fn coords(&self) -> (f64, f64);
    fn x(&self) -> f64;
    fn y(&self) -> f64;
    fn xy(&self) -> NXY;
    fn is_equal(&self, other: &Self, linear_tolerance: f64) -> bool;
    fn distance(&self, other: &Self) -> f64;
    fn square_distance(&self, other: &Self) -> f64;
    fn mirror_pnt(&mut self, p: &Self);
    fn mirrored_pnt(&self, p: &Self) -> Self;
    fn mirror_ax2d(&mut self, a: &NAx2d);
    fn mirrored_ax2d(&self, a: &NAx2d) -> Self;
    fn rotate(&mut self, p: &Self, ang: f64);
    fn rotated(&self, p: &Self, ang: f64) -> Self;
    fn scale(&mut self, p: &Self, s: f64);
    fn scaled(&self, p: &Self, s: f64) -> Self;
    fn transform(&mut self, t: &NTrsf2d);
    fn transformed(&self, t: &NTrsf2d) -> Self;
    fn translate_vec(&mut self, v: &NVec2d);
    fn translated_vec(&self, v: &NVec2d) -> Self;
    fn translate_point2d(&mut self, p1: &Self, p2: &Self);
    fn translated_point2d(&self, p1: &Self, p2: &Self) -> Self;
}

#[derive(Clone, Default, PartialEq, Debug, Serialize, Deserialize)]
pub struct NPoint2d {
    coord: NXY,
}

impl Eq for NPoint2d {}

impl From<f64> for NPoint2d {
    fn from(value: f64) -> Self {
        Self::new(value, value)
    }
}

impl From<(f64, f64)> for NPoint2d {
    fn from(value: (f64, f64)) -> Self {
        Self::new(value.0, value.1)
    }
}

impl From<NXY> for NPoint2d {
    fn from(value: NXY) -> Self {
        Self { coord: value }
    }
}

impl Point2d for NPoint2d {
    fn set_coord_by_index(&mut self, index: usize, value: f64) -> Result<(), NErrors> {
        self.coord.set_coord_by_index(index, value)
    }

    fn set_coords(&mut self, x: f64, y: f64) {
        self.coord.set_coords(x, y);
    }

    fn set_x(&mut self, x: f64) {
        self.coord.set_x(x);
    }

    fn set_y(&mut self, y: f64) {
        self.coord.set_y(y);
    }

    fn set_xy(&mut self, coord: NXY) {
        self.coord = coord;
    }

    fn coord_by_index(&self, index: usize) -> Result<f64, NErrors> {
        self.coord.coord_by_index(index)
    }

    fn coords(&self) -> (f64, f64) {
        (self.x(), self.y())
    }

    fn x(&self) -> f64 {
        self.coord.x()
    }

    fn y(&self) -> f64 {
        self.coord.y()
    }

    fn xy(&self) -> NXY {
        self.coord.clone()
    }

    fn is_equal(&self, other: &Self, linear_tolerance: f64) -> bool {
        self.distance(other) <= linear_tolerance
    }

    fn distance(&self, other: &Self) -> f64 {
        f64::sqrt(self.square_distance(other))
    }

    fn square_distance(&self, other: &Self) -> f64 {
        let dx = self.coord.x() - other.coord.x();
        let dy = self.coord.y() - other.coord.y();

        dx.powi(2) + dy.powi(2)
    }

    fn mirror_pnt(&mut self, p: &Self) {
        let mut xyz = p.coord.clone();
        xyz.multiply(2.0);
        self.coord.reverse();
        self.coord.add(&xyz);
    }

    fn mirrored_pnt(&self, p: &Self) -> Self {
        let mut res = self.clone();
        res.mirror_pnt(p);
        res
    }

    fn mirror_ax2d(&mut self, a: &NAx2d) {
        let mut trsf = NTrsf2d::default();
        trsf.set_mirror_ax2d(a);
        trsf.transforms_xy(&mut self.coord);
    }

    fn mirrored_ax2d(&self, a: &NAx2d) -> Self {
        let mut res = self.clone();
        res.mirror_ax2d(a);
        res
    }

    fn rotate(&mut self, p: &Self, ang: f64) {
        let mut t = NTrsf2d::default();
        t.set_rotation(p, ang);
        t.transforms_xy(&mut self.coord);
    }

    fn rotated(&self, p: &Self, ang: f64) -> Self {
        let mut res = self.clone();
        res.rotate(p, ang);
        res
    }

    fn scale(&mut self, p: &Self, scalar: f64) {
        let mut p_coord = p.coord.clone();
        p_coord.multiply(1.0 - scalar);
        self.coord.multiply(scalar);
        self.coord.add(&p_coord);
    }

    fn scaled(&self, p: &Self, s: f64) -> Self {
        let mut res = self.clone();
        res.scale(p, s);
        res
    }

    fn transform(&mut self, t: &NTrsf2d) {
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
                t.transforms_xy(&mut self.coord);
            }
        }
    }

    fn transformed(&self, t: &NTrsf2d) -> Self {
        let mut res = self.clone();
        res.transform(t);
        res
    }

    fn translate_vec(&mut self, v: &NVec2d) {
        self.coord.add(&v.xy());
    }

    fn translated_vec(&self, v: &NVec2d) -> Self {
        let mut res = self.clone();
        res.translate_vec(v);
        res
    }

    fn translate_point2d(&mut self, p1: &Self, p2: &Self) {
        self.coord.add(&p2.coord);
        self.coord.subtract(&p1.coord);
    }

    fn translated_point2d(&self, p1: &Self, p2: &Self) -> Self {
        let mut res = self.clone();
        res.translate_point2d(p1, p2);
        res
    }

    fn new(x: f64, y: f64) -> Self {
        Self {
            coord: NXY::new(x, y),
        }
    }

    fn change_coord_xy(&mut self) -> &mut NXY {
        &mut self.coord
    }

    fn change_coord(&mut self, index: usize) -> Result<&mut f64, NErrors> {
        self.coord.change_coord(index)
    }
}

// Implement std::hash::Hash for NPoint3d
impl std::hash::Hash for NPoint2d {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let x_bits = self.x().to_bits();
        let y_bits = self.y().to_bits();
        // Simplified hashing based on C++ logic
        (x_bits / 23 + y_bits / 19).hash(state);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_point() -> NPoint2d {
        NPoint2d::new(1.0, 2.0)
    }

    #[test]
    fn test_new() {
        let p = NPoint2d::default();
        assert_eq!(p.coords(), (0.0, 0.0));
    }

    #[test]
    fn test_new_with_coords() {
        let p = NPoint2d::new(1.0, 2.0);
        assert_eq!(p.x(), 1.0);
        assert_eq!(p.y(), 2.0);
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
        assert_eq!(p.coords(), (1.0, 2.0));
        assert_eq!(p.coord_by_index(1).unwrap(), 1.0);
        assert!(matches!(p.coord_by_index(4), Err(NErrors::IndexOutOfRange)));
    }

    #[test]
    fn test_distance() {
        let p1 = NPoint2d::new(0.0, 0.0);
        let p2 = NPoint2d::new(3.0, 4.0);
        assert!((p1.distance(&p2) - 5.0).abs() < 1e-9);
        assert!((p1.square_distance(&p2) - 25.0).abs() < 1e-9);
    }

    #[test]
    fn test_is_equal() {
        let p1 = NPoint2d::new(1.0, 2.0);
        let p2 = NPoint2d::new(1.0 + 1e-6, 2.0);
        assert!(p1.is_equal(&p2, 1e-5));
        assert!(!p1.is_equal(&p2, 1e-7));
    }

    #[test]
    fn test_mirror_pnt() {
        let mut p = NPoint2d::new(1.0, 1.0);
        let center = NPoint2d::new(0.0, 0.0);
        p.mirror_pnt(&center);
        assert_eq!(p.coords(), (-1.0, -1.0));
    }

    #[test]
    fn test_scale() {
        let mut p = NPoint2d::new(2.0, 2.0);
        let center = NPoint2d::new(0.0, 0.0);
        p.scale(&center, 2.0);
        assert_eq!(p.coords(), (4.0, 4.0));
    }

    #[test]
    fn test_translate_vec() {
        let mut p = NPoint2d::new(1.0, 1.0);
        let v = NVec2d::new_from_coords(1.0, 2.0);
        p.translate_vec(&v);
        assert_eq!(p.coords(), (2.0, 3.0));
    }
}
