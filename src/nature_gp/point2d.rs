use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::{
    gp::{NAx2d, NTrsf2d, NVec2d, NXY},
    nature_errors::NErrors,
};

// Trait to define the behavior of a 2D Cartesian point
pub trait Pnt2d {
    fn new() -> Self;
    fn new_with_xy(coord: &NXY) -> Self;
    fn new_with_coords(x: f64, y: f64) -> Self;
    fn set_coord(&mut self, index: i32, xi: f64) -> Result<(), NErrors>;
    fn set_coords(&mut self, x: f64, y: f64);
    fn set_x(&mut self, x: f64);
    fn set_y(&mut self, y: f64);
    fn set_xy(&mut self, coord: &NXY);
    fn coord(&self, index: i32) -> Result<f64, NErrors>;
    fn coords(&self) -> (f64, f64);
    fn x(&self) -> f64;
    fn y(&self) -> f64;
    fn xy(&self) -> NXY;
    fn is_equal(&self, other: &Self, linear_tolerance: f64) -> bool;
    fn distance(&self, other: &Self) -> f64;
    fn square_distance(&self, other: &Self) -> f64;
    fn mirror_pnt(&mut self, p: &Self);
    fn mirrored_pnt(&self, p: &Self) -> Self
    where
        Self: Sized;
    fn mirror_ax2d(&mut self, a: &NAx2d);
    fn mirrored_ax2d(&self, a: &NAx2d) -> Self
    where
        Self: Sized;
    fn rotate(&mut self, p: &Self, ang: f64);
    fn rotated(&self, p: &Self, ang: f64) -> Self
    where
        Self: Sized;
    fn scale(&mut self, p: &Self, s: f64);
    fn scaled(&self, p: &Self, s: f64) -> Self
    where
        Self: Sized;
    fn transform(&mut self, t: &NTrsf2d);
    fn transformed(&self, t: &NTrsf2d) -> Self
    where
        Self: Sized;
    fn translate_vec(&mut self, v: &NVec2d);
    fn translated_vec(&self, v: &NVec2d) -> Self
    where
        Self: Sized;
    fn translate_pnts(&mut self, p1: &Self, p2: &Self);
    fn translated_pnts(&self, p1: &Self, p2: &Self) -> Self
    where
        Self: Sized;
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
}

// Struct representing a 2D Cartesian point
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NPnt2d {
    coord: NXY,
}

impl Pnt2d for NPnt2d {
    /// Creates a point with zero coordinates.
    fn new() -> Self {
        NPnt2d {
            coord: NXY::new(0.0, 0.0),
        }
    }

    /// Creates a point from an XY object.
    fn new_with_xy(coord: &NXY) -> Self {
        NPnt2d {
            coord: coord.clone(),
        }
    }

    /// Creates a point with its 2 Cartesian coordinates.
    fn new_with_coords(x: f64, y: f64) -> Self {
        NPnt2d {
            coord: NXY::new(x, y),
        }
    }

    /// Changes the coordinate at the given index (1=X, 2=Y).
    fn set_coord(&mut self, index: i32, xi: f64) -> Result<(), NErrors> {
        match index {
            1 => self.coord.set_x(xi),
            2 => self.coord.set_y(xi),
            _ => return Err(NErrors::OutOfRange),
        }
        Ok(())
    }

    /// Assigns the given values to the coordinates.
    fn set_coords(&mut self, x: f64, y: f64) {
        self.coord.set_coords(x, y);
    }

    /// Assigns the given value to the X coordinate.
    fn set_x(&mut self, x: f64) {
        self.coord.set_x(x);
    }

    /// Assigns the given value to the Y coordinate.
    fn set_y(&mut self, y: f64) {
        self.coord.set_y(y);
    }

    /// Assigns the coordinates of the given XY object.
    fn set_xy(&mut self, coord: &NXY) {
        self.coord = coord.clone();
    }

    /// Returns the coordinate at the given index (1=X, 2=Y).
    fn coord(&self, index: i32) -> Result<f64, NErrors> {
        match index {
            1 => Ok(self.coord.x()),
            2 => Ok(self.coord.y()),
            _ => Err(NErrors::OutOfRange),
        }
    }

    /// Returns the two coordinates as a tuple.
    fn coords(&self) -> (f64, f64) {
        (self.coord.x(), self.coord.y())
    }

    /// Returns the X coordinate.
    fn x(&self) -> f64 {
        self.coord.x()
    }

    /// Returns the Y coordinate.
    fn y(&self) -> f64 {
        self.coord.y()
    }

    /// Returns the coordinates as an XY object.
    fn xy(&self) -> NXY {
        self.coord.clone()
    }

    /// Returns true if the distance to the other point is within the linear tolerance.
    fn is_equal(&self, other: &Self, linear_tolerance: f64) -> bool {
        self.distance(other) <= linear_tolerance
    }

    /// Computes the Euclidean distance to another point.
    fn distance(&self, other: &Self) -> f64 {
        let dx = self.coord.x() - other.coord.x();
        let dy = self.coord.y() - other.coord.y();
        (dx * dx + dy * dy).sqrt()
    }

    /// Computes the square of the Euclidean distance to another point.
    fn square_distance(&self, other: &Self) -> f64 {
        let dx = self.coord.x() - other.coord.x();
        let dy = self.coord.y() - other.coord.y();
        dx * dx + dy * dy
    }

    /// Mirrors the point with respect to another point.
    fn mirror_pnt(&mut self, p: &Self) {
        let mut p_coord = p.coord.clone();
        p_coord.multiply(2.0);
        self.coord.reverse();
        self.coord.add(&p_coord);
    }

    /// Returns the point mirrored with respect to another point.
    fn mirrored_pnt(&self, p: &Self) -> Self {
        let mut res = self.clone();
        res.mirror_pnt(p);
        res
    }

    /// Mirrors the point with respect to an axis.
    fn mirror_ax2d(&mut self, a: &NAx2d) {
        let mut t = NTrsf2d::new();
        t.set_mirror(a);
        t.transforms(&mut self.coord);
    }

    /// Returns the point mirrored with respect to an axis.
    fn mirrored_ax2d(&self, a: &NAx2d) -> Self {
        let mut res = self.clone();
        res.mirror_ax2d(a);
        res
    }

    /// Rotates the point around a center point by an angle.
    fn rotate(&mut self, p: &Self, ang: f64) {
        let mut t = NTrsf2d::new();
        t.set_rotation(p, ang);
        t.transforms(&mut self.coord);
    }

    /// Returns the point rotated around a center point by an angle.
    fn rotated(&self, p: &Self, ang: f64) -> Self {
        let mut res = self.clone();
        res.rotate(p, ang);
        res
    }

    /// Scales the point with respect to another point.
    fn scale(&mut self, p: &Self, s: f64) {
        let mut p_coord = p.coord.clone();
        p_coord.multiply(1.0 - s);
        self.coord.multiply(s);
        self.coord.add(&p_coord);
    }

    /// Returns the point scaled with respect to another point.
    fn scaled(&self, p: &Self, s: f64) -> Self {
        let mut res = self.clone();
        res.scale(p, s);
        res
    }

    /// Transforms the point with a 2D transformation.
    fn transform(&mut self, t: &NTrsf2d) {
        match t.form() {
            NTrsfForm2d::Identity => {}
            NTrsfForm2d::Translation => {
                self.coord.add(&t.translation_part());
            }
            NTrsfForm2d::Scale => {
                self.coord.multiply(t.scale_factor());
                self.coord.add(&t.translation_part());
            }
            NTrsfForm2d::PntMirror => {
                self.coord.reverse();
                self.coord.add(&t.translation_part());
            }
            _ => {
                t.transforms(&mut self.coord);
            }
        }
    }

    /// Returns the point transformed with a 2D transformation.
    fn transformed(&self, t: &NTrsf2d) -> Self {
        let mut res = self.clone();
        res.transform(t);
        res
    }

    /// Translates the point by a 2D vector.
    fn translate_vec(&mut self, v: &NVec2d) {
        self.coord.add(&v.xy());
    }

    /// Returns the point translated by a 2D vector.
    fn translated_vec(&self, v: &NVec2d) -> Self {
        let mut res = self.clone();
        res.translate_vec(v);
        res
    }

    /// Translates the point from one point to another.
    fn translate_pnts(&mut self, p1: &Self, p2: &Self) {
        self.coord.add(&p2.coord);
        self.coord.subtract(&p1.coord);
    }

    /// Returns the point translated from one point to another.
    fn translated_pnts(&self, p1: &Self, p2: &Self) -> Self {
        let mut res = self.clone();
        res.translate_pnts(p1, p2);
        res
    }

    /// Dumps the point as JSON.
    fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(out, "{} {{", indent).unwrap();
        writeln!(out, "{}   \"type\": \"NPnt2d\",", indent).unwrap();
        writeln!(
            out,
            "{}   \"coordinates\": [{}, {}]",
            indent,
            self.coord.x(),
            self.coord.y()
        )
        .unwrap();
        writeln!(out, "{} }}", indent).unwrap();
    }
}

// Implement std::hash::Hash for NPnt2d
impl std::hash::Hash for NPnt2d {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let x_bits = self.x().to_bits();
        let y_bits = self.y().to_bits();
        // Simplified hashing based on 2D coordinates
        (x_bits / 23 + y_bits / 19).hash(state);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_point() -> NPnt2d {
        NPnt2d::new_with_coords(1.0, 2.0)
    }

    #[test]
    fn test_new() {
        let p = NPnt2d::new();
        assert_eq!(p.coords(), (0.0, 0.0));
    }

    #[test]
    fn test_new_with_coords() {
        let p = NPnt2d::new_with_coords(1.0, 2.0);
        assert_eq!(p.x(), 1.0);
        assert_eq!(p.y(), 2.0);
    }

    #[test]
    fn test_set_coord() {
        let mut p = create_test_point();
        p.set_coord(1, 4.0).unwrap();
        assert_eq!(p.x(), 4.0);
        assert!(matches!(p.set_coord(3, 5.0), Err(NErrors::OutOfRange)));
    }

    #[test]
    fn test_coords() {
        let p = create_test_point();
        assert_eq!(p.coords(), (1.0, 2.0));
        assert_eq!(p.coord(1).unwrap(), 1.0);
        assert!(matches!(p.coord(3), Err(NErrors::OutOfRange)));
    }

    #[test]
    fn test_distance() {
        let p1 = NPnt2d::new_with_coords(0.0, 0.0);
        let p2 = NPnt2d::new_with_coords(3.0, 4.0);
        assert!((p1.distance(&p2) - 5.0).abs() < 1e-9);
        assert!((p1.square_distance(&p2) - 25.0).abs() < 1e-9);
    }

    #[test]
    fn test_is_equal() {
        let p1 = NPnt2d::new_with_coords(1.0, 2.0);
        let p2 = NPnt2d::new_with_coords(1.0 + 1e-6, 2.0);
        assert!(p1.is_equal(&p2, 1e-5));
        assert!(!p1.is_equal(&p2, 1e-7));
    }

    #[test]
    fn test_mirror_pnt() {
        let mut p = NPnt2d::new_with_coords(1.0, 1.0);
        let center = NPnt2d::new_with_coords(0.0, 0.0);
        p.mirror_pnt(&center);
        assert_eq!(p.coords(), (-1.0, -1.0));
    }

    #[test]
    fn test_scale() {
        let mut p = NPnt2d::new_with_coords(2.0, 2.0);
        let center = NPnt2d::new_with_coords(0.0, 0.0);
        p.scale(&center, 2.0);
        assert_eq!(p.coords(), (4.0, 4.0));
    }

    #[test]
    fn test_translate_vec() {
        let mut p = NPnt2d::new_with_coords(1.0, 1.0);
        let v = NVec2d::new_with_xy(&NXY::new(1.0, 2.0));
        p.translate_vec(&v);
        assert_eq!(p.coords(), (2.0, 3.0));
    }

    #[test]
    fn test_dump_json() {
        let p = create_test_point();
        let mut output = Vec::new();
        p.dump_json(&mut output, 0);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NPnt2d\""));
        assert!(json.contains("[1, 2]"));
    }
}
