use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::{
    gp::{NAx1, NAx2, NTrsf, NVec, NXYZ},
    nature_errors::NErrors,
};

// Trait to define the behavior of a 3D Cartesian point
pub trait Pnt {
    fn new() -> Self;
    fn new_with_xyz(coord: &NXYZ) -> Self;
    fn new_with_coords(x: f64, y: f64, z: f64) -> Self;
    fn set_coord(&mut self, index: i32, xi: f64) -> Result<(), NErrors>;
    fn set_coords(&mut self, x: f64, y: f64, z: f64);
    fn set_x(&mut self, x: f64);
    fn set_y(&mut self, y: f64);
    fn set_z(&mut self, z: f64);
    fn set_xyz(&mut self, coord: &NXYZ);
    fn coord(&self, index: i32) -> Result<f64, NErrors>;
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
    fn translate_pnts(&mut self, p1: &Self, p2: &Self);
    fn translated_pnts(&self, p1: &Self, p2: &Self) -> Self
    where
        Self: Sized;
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
    // Note: InitFromJson is not included as it requires stream parsing, which is less idiomatic in Rust
}

// Struct representing a 3D Cartesian point
#[derive(Clone, Default, PartialEq, Debug, Serialize, Deserialize)]
pub struct NPnt {
    coord: NXYZ,
}

impl Pnt for NPnt {
    /// Creates a point with zero coordinates.
    fn new() -> Self {
        NPnt {
            coord: NXYZ::new(0.0, 0.0, 0.0),
        }
    }

    /// Creates a point from an XYZ object.
    fn new_with_xyz(coord: &NXYZ) -> Self {
        NPnt {
            coord: coord.clone(),
        }
    }

    /// Creates a point with its 3 Cartesian coordinates.
    fn new_with_coords(x: f64, y: f64, z: f64) -> Self {
        NPnt {
            coord: NXYZ::new(x, y, z),
        }
    }

    /// Changes the coordinate at the given index (1=X, 2=Y, 3=Z).
    fn set_coord(&mut self, index: i32, xi: f64) -> Result<(), NErrors> {
        match index {
            1 => self.coord.set_x(xi),
            2 => self.coord.set_y(xi),
            3 => self.coord.set_z(xi),
            _ => return Err(NErrors::OutOfRange),
        }
        Ok(())
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
    fn set_xyz(&mut self, coord: &NXYZ) {
        self.coord = coord.clone();
    }

    /// Returns the coordinate at the given index (1=X, 2=Y, 3=Z).
    fn coord(&self, index: i32) -> Result<f64, NErrors> {
        match index {
            1 => Ok(self.coord.x()),
            2 => Ok(self.coord.y()),
            3 => Ok(self.coord.z()),
            _ => Err(NErrors::OutOfRange),
        }
    }

    /// Returns the three coordinates as a tuple.
    fn coords(&self) -> (f64, f64, f64) {
        (self.coord.x(), self.coord.y(), self.coord.z())
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
        let mut new_coord = NXYZ::new(0.0, 0.0, 0.0);
        new_coord.set_linear_form(alpha, &self.coord, beta, &other.coord);
        new_coord.divide(alpha + beta);
        self.coord = new_coord;
    }

    /// Returns true if the distance to the other point is within the linear tolerance.
    fn is_equal(&self, other: &Self, linear_tolerance: f64) -> bool {
        self.distance(other) <= linear_tolerance
    }

    /// Computes the Euclidean distance to another point.
    fn distance(&self, other: &Self) -> f64 {
        let dx = self.coord.x() - other.coord.x();
        let dy = self.coord.y() - other.coord.y();
        let dz = self.coord.z() - other.coord.z();
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Computes the square of the Euclidean distance to another point.
    fn square_distance(&self, other: &Self) -> f64 {
        let dx = self.coord.x() - other.coord.x();
        let dy = self.coord.y() - other.coord.y();
        let dz = self.coord.z() - other.coord.z();
        dx * dx + dy * dy + dz * dz
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
    fn mirror_ax1(&mut self, a1: &NAx1) {
        let mut t = NTrsf::new();
        t.set_mirror_ax1(a1);
        t.transforms(&mut self.coord);
    }

    /// Returns the point mirrored with respect to an axis.
    fn mirrored_ax1(&self, a1: &NAx1) -> Self {
        let mut res = self.clone();
        res.mirror_ax1(a1);
        res
    }

    /// Mirrors the point with respect to a plane defined by Ax2.
    fn mirror_ax2(&mut self, a2: &NAx2) {
        let mut t = NTrsf::new();
        t.set_mirror_ax2(a2);
        t.transforms(&mut self.coord);
    }

    /// Returns the point mirrored with respect to a plane defined by Ax2.
    fn mirrored_ax2(&self, a2: &NAx2) -> Self {
        let mut res = self.clone();
        res.mirror_ax2(a2);
        res
    }

    /// Rotates the point around an axis by an angle.
    fn rotate(&mut self, a1: &NAx1, ang: f64) {
        let mut t = NTrsf::new();
        t.set_rotation(a1, ang);
        t.transforms(&mut self.coord);
    }

    /// Returns the point rotated around an axis by an angle.
    fn rotated(&self, a1: &NAx1, ang: f64) -> Self {
        let mut res = self.clone();
        res.rotate(a1, ang);
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
            NTrsfForm::PntMirror => {
                self.coord.reverse();
                self.coord.add(&t.translation_part());
            }
            _ => {
                t.transforms(&mut self.coord);
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
        writeln!(out, "{}   \"type\": \"NPnt\",", indent).unwrap();
        writeln!(
            out,
            "{}   \"coordinates\": [{}, {}, {}]",
            indent,
            self.coord.x(),
            self.coord.y(),
            self.coord.z()
        )
        .unwrap();
        writeln!(out, "{} }}", indent).unwrap();
    }
}

// Implement std::hash::Hash for NPnt
impl std::hash::Hash for NPnt {
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

    fn create_test_point() -> NPnt {
        NPnt::new_with_coords(1.0, 2.0, 3.0)
    }

    #[test]
    fn test_new() {
        let p = NPnt::new();
        assert_eq!(p.coords(), (0.0, 0.0, 0.0));
    }

    #[test]
    fn test_new_with_coords() {
        let p = NPnt::new_with_coords(1.0, 2.0, 3.0);
        assert_eq!(p.x(), 1.0);
        assert_eq!(p.y(), 2.0);
        assert_eq!(p.z(), 3.0);
    }

    #[test]
    fn test_set_coord() {
        let mut p = create_test_point();
        p.set_coord(1, 4.0).unwrap();
        assert_eq!(p.x(), 4.0);
        assert!(matches!(p.set_coord(4, 5.0), Err(NErrors::OutOfRange)));
    }

    #[test]
    fn test_coords() {
        let p = create_test_point();
        assert_eq!(p.coords(), (1.0, 2.0, 3.0));
        assert_eq!(p.coord(1).unwrap(), 1.0);
        assert!(matches!(p.coord(4), Err(NErrors::OutOfRange)));
    }

    #[test]
    fn test_distance() {
        let p1 = NPnt::new_with_coords(0.0, 0.0, 0.0);
        let p2 = NPnt::new_with_coords(3.0, 4.0, 0.0);
        assert!((p1.distance(&p2) - 5.0).abs() < 1e-9);
        assert!((p1.square_distance(&p2) - 25.0).abs() < 1e-9);
    }

    #[test]
    fn test_is_equal() {
        let p1 = NPnt::new_with_coords(1.0, 2.0, 3.0);
        let p2 = NPnt::new_with_coords(1.0 + 1e-6, 2.0, 3.0);
        assert!(p1.is_equal(&p2, 1e-5));
        assert!(!p1.is_equal(&p2, 1e-7));
    }

    #[test]
    fn test_bary_center() {
        let mut p1 = NPnt::new_with_coords(1.0, 0.0, 0.0);
        let p2 = NPnt::new_with_coords(0.0, 1.0, 0.0);
        p1.bary_center(1.0, &p2, 1.0);
        let (x, y, z) = p1.coords();
        assert!((x - 0.5).abs() < 1e-9);
        assert!((y - 0.5).abs() < 1e-9);
        assert!((z - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_mirror_pnt() {
        let mut p = NPnt::new_with_coords(1.0, 1.0, 1.0);
        let center = NPnt::new_with_coords(0.0, 0.0, 0.0);
        p.mirror_pnt(&center);
        assert_eq!(p.coords(), (-1.0, -1.0, -1.0));
    }

    #[test]
    fn test_scale() {
        let mut p = NPnt::new_with_coords(2.0, 2.0, 2.0);
        let center = NPnt::new_with_coords(0.0, 0.0, 0.0);
        p.scale(&center, 2.0);
        assert_eq!(p.coords(), (4.0, 4.0, 4.0));
    }

    #[test]
    fn test_translate_vec() {
        let mut p = NPnt::new_with_coords(1.0, 1.0, 1.0);
        let v = NVec::new_with_xyz(&NXYZ::new(1.0, 2.0, 3.0));
        p.translate_vec(&v);
        assert_eq!(p.coords(), (2.0, 3.0, 4.0));
    }

    #[test]
    fn test_dump_json() {
        let p = create_test_point();
        let mut output = Vec::new();
        p.dump_json(&mut output, 0);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NPnt\""));
        assert!(json.contains("[1, 2, 3]"));
    }
}
