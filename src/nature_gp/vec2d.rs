use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::gp::{NAx2d, NDir2d, NPnt2d, NTrsf2d, NXY, NatureError};

mod gp {
    pub fn resolution() -> f64 {
        1e-12 // Consistent with trsf.rs and vec.rs
    }
}

/// Represents a non-persistent 2D vector.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NVec2d {
    coord: NXY,
}

impl NVec2d {
    /// Creates a zero vector.
    pub fn new() -> Self {
        NVec2d {
            coord: NXY::new(0.0, 0.0),
        }
    }

    /// Creates a unitary vector from a direction.
    pub fn new_from_dir(dir: &NDir2d) -> Self {
        NVec2d {
            coord: dir.xy(),
        }
    }

    /// Creates a vector with the given coordinates.
    pub fn new_from_coords(x: f64, y: f64) -> Self {
        NVec2d {
            coord: NXY::new(x, y),
        }
    }

    /// Creates a vector from two points (P2 - P1).
    pub fn new_from_points(p1: &NPnt2d, p2: &NPnt2d) -> Self {
        NVec2d {
            coord: p2.xy().subtracted(&p1.xy()),
        }
    }

    /// Sets the coordinate at the given index (1=X, 2=Y).
    pub fn set_coord(&mut self, index: i32, value: f64) -> Result<(), NatureError> {
        match index {
            1 => self.coord.set_x(value),
            2 => self.coord.set_y(value),
            _ => return Err(NatureError::OutOfRange),
        }
        Ok(())
    }

    /// Sets all coordinates.
    pub fn set_coords(&mut self, x: f64, y: f64) {
        self.coord.set_coord(x, y);
    }

    /// Sets the X coordinate.
    pub fn set_x(&mut self, x: f64) {
        self.coord.set_x(x);
    }

    /// Sets the Y coordinate.
    pub fn set_y(&mut self, y: f64) {
        self.coord.set_y(y);
    }

    /// Sets the coordinates from an NXY.
    pub fn set_xy(&mut self, xy: &NXY) {
        self.coord = xy.clone();
    }

    /// Gets the coordinate at the given index (1=X, 2=Y).
    pub fn coord(&self, index: i32) -> Result<f64, NatureError> {
        match index {
            1 => Ok(self.coord.x()),
            2 => Ok(self.coord.y()),
            _ => Err(NatureError::OutOfRange),
        }
    }

    /// Gets all coordinates.
    pub fn coords(&self) -> (f64, f64) {
        (self.coord.x(), self.coord.y())
    }

    /// Gets the X coordinate.
    pub fn x(&self) -> f64 {
        self.coord.x()
    }

    /// Gets the Y coordinate.
    pub fn y(&self) -> f64 {
        self.coord.y()
    }

    /// Gets the coordinates as an NXY.
    pub fn xy(&self) -> NXY {
        self.coord.clone()
    }

    /// Checks if two vectors are equal within tolerances.
    pub fn is_equal(&self, other: &NVec2d, linear_tolerance: f64, angular_tolerance: f64) -> bool {
        let norm = self.magnitude();
        let other_norm = other.magnitude();
        let val = (norm - other_norm).abs();
        let is_equal_length = val <= linear_tolerance;
        if norm > linear_tolerance && other_norm > linear_tolerance {
            let mut ang = self.angle(other).unwrap_or(std::f64::INFINITY).abs();
            if ang > std::f64::consts::PI {
                ang = 2.0 * std::f64::consts::PI - ang;
            }
            is_equal_length && ang <= angular_tolerance
        } else {
            is_equal_length
        }
    }

    /// Checks if the vector is normal to another within angular tolerance.
    pub fn is_normal(&self, other: &NVec2d, angular_tolerance: f64) -> Result<bool, NatureError> {
        let angle = (std::f64::consts::PI / 2.0 - self.angle(other)?.abs()).abs();
        Ok(angle <= angular_tolerance)
    }

    /// Checks if the vector is opposite to another within angular tolerance.
    pub fn is_opposite(&self, other: &NVec2d, angular_tolerance: f64) -> Result<bool, NatureError> {
        let mut angle = self.angle(other)?.abs();
        if angle > std::f64::consts::PI {
            angle = 2.0 * std::f64::consts::PI - angle;
        }
        Ok((std::f64::consts::PI - angle) <= angular_tolerance)
    }

    /// Checks if the vector is parallel to another within angular tolerance.
    pub fn is_parallel(&self, other: &NVec2d, angular_tolerance: f64) -> Result<bool, NatureError> {
        let mut angle = self.angle(other)?.abs();
        if angle > std::f64::consts::PI {
            angle = 2.0 * std::f64::consts::PI - angle;
        }
        Ok(angle <= angular_tolerance || (std::f64::consts::PI - angle) <= angular_tolerance)
    }

    /// Computes the angle between two vectors (-PI to PI radians).
    pub fn angle(&self, other: &NVec2d) -> Result<f64, NatureError> {
        let norm = self.magnitude();
        let other_norm = other.magnitude();
        if norm <= gp::resolution() || other_norm <= gp::resolution() {
            return Err(NatureError::VectorWithNullMagnitude);
        }
        let d = norm * other_norm;
        let cosinus = self.coord.dot(&other.coord) / d;
        let sinus = self.coord.crossed(&other.coord) / d;
        if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
            Ok(if sinus > 0.0 { cosinus.acos() } else { -cosinus.acos() })
        } else if cosinus > 0.0 {
            Ok(sinus.asin())
        } else if sinus > 0.0 {
            Ok(std::f64::consts::PI - sinus.asin())
        } else {
            Ok(-std::f64::consts::PI - sinus.asin())
        }
    }

    /// Computes the magnitude of the vector.
    pub fn magnitude(&self) -> f64 {
        self.coord.modulus()
    }

    /// Computes the square magnitude of the vector.
    pub fn square_magnitude(&self) -> f64 {
        self.coord.square_modulus()
    }

    /// Adds another vector to this one.
    pub fn add(&mut self, other: &NVec2d) {
        self.coord.add(&other.coord);
    }

    /// Returns the sum of two vectors.
    pub fn added(&self, other: &NVec2d) -> NVec2d {
        NVec2d {
            coord: self.coord.added(&other.coord),
        }
    }

    /// Subtracts another vector from this one.
    pub fn subtract(&mut self, other: &NVec2d) {
        self.coord.subtract(&other.coord);
    }

    /// Returns the difference of two vectors.
    pub fn subtracted(&self, other: &NVec2d) -> NVec2d {
        NVec2d {
            coord: self.coord.subtracted(&other.coord),
        }
    }

    /// Multiplies the vector by a scalar.
    pub fn multiply(&mut self, scalar: f64) {
        self.coord.multiply(scalar);
    }

    /// Returns the vector multiplied by a scalar.
    pub fn multiplied(&self, scalar: f64) -> NVec2d {
        NVec2d {
            coord: self.coord.multiplied(scalar),
        }
    }

    /// Divides the vector by a scalar.
    pub fn divide(&mut self, scalar: f64) -> Result<(), NatureError> {
        if scalar.abs() <= gp::resolution() {
            return Err(NatureError::InvalidConstructionParameters);
        }
        self.coord.divide(scalar);
        Ok(())
    }

    /// Returns the vector divided by a scalar.
    pub fn divided(&self, scalar: f64) -> Result<NVec2d, NatureError> {
        if scalar.abs() <= gp::resolution() {
            return Err(NatureError::InvalidConstructionParameters);
        }
        Ok(NVec2d {
            coord: self.coord.divided(scalar),
        })
    }

    /// Computes the cross product with another vector (returns a scalar).
    pub fn crossed(&self, other: &NVec2d) -> f64 {
        self.coord.crossed(&other.coord)
    }

    /// Computes the magnitude of the cross product with another vector.
    pub fn cross_magnitude(&self, other: &NVec2d) -> f64 {
        self.coord.cross_magnitude(&other.coord)
    }

    /// Computes the square magnitude of the cross product with another vector.
    pub fn cross_square_magnitude(&self, other: &NVec2d) -> f64 {
        self.coord.cross_square_magnitude(&other.coord)
    }

    /// Computes the dot product with another vector.
    pub fn dot(&self, other: &NVec2d) -> f64 {
        self.coord.dot(&other.coord)
    }

    /// Returns a vector normal to this one (rotated 90 degrees counterclockwise).
    pub fn get_normal(&self) -> NVec2d {
        NVec2d {
            coord: NXY::new(self.y(), -self.x()),
        }
    }

    /// Normalizes the vector.
    pub fn normalize(&mut self) -> Result<(), NatureError> {
        let mag = self.magnitude();
        if mag <= gp::resolution() {
            return Err(NatureError::VectorWithNullMagnitude);
        }
        self.coord.divide(mag);
        Ok(())
    }

    /// Returns a normalized copy of the vector.
    pub fn normalized(&self) -> Result<NVec2d, NatureError> {
        let mag = self.magnitude();
        if mag <= gp::resolution() {
            return Err(NatureError::VectorWithNullMagnitude);
        }
        Ok(NVec2d {
            coord: self.coord.divided(mag),
        })
    }

    /// Reverses the direction of the vector.
    pub fn reverse(&mut self) {
        self.coord.reverse();
    }

    /// Returns a reversed copy of the vector.
    pub fn reversed(&self) -> NVec2d {
        NVec2d {
            coord: self.coord.reversed(),
        }
    }

    /// Sets the vector to a linear combination: a1*v1 + a2*v2 + v3.
    pub fn set_linear_form(&mut self, a1: f64, v1: &NVec2d, a2: f64, v2: &NVec2d, v3: &NVec2d) {
        self.coord.set_linear_form(a1, &v1.coord, a2, &v2.coord, &v3.coord);
    }

    /// Sets the vector to a linear combination: a1*v1 + a2*v2.
    pub fn set_linear_form_2(&mut self, a1: f64, v1: &NVec2d, a2: f64, v2: &NVec2d) {
        self.coord.set_linear_form(a1, &v1.coord, a2, &v2.coord);
    }

    /// Sets the vector to a linear combination: a1*v1 + v2.
    pub fn set_linear_form_1_plus(&mut self, a1: f64, v1: &NVec2d, v2: &NVec2d) {
        self.coord.set_linear_form(a1, &v1.coord, &v2.coord);
    }

    /// Sets the vector to a linear combination: v1 + v2.
    pub fn set_linear_form_sum(&mut self, v1: &NVec2d, v2: &NVec2d) {
        self.coord.set_linear_form(&v1.coord, &v2.coord);
    }

    /// Mirrors the vector with respect to another vector.
    pub fn mirror_vec(&mut self, v: &NVec2d) -> Result<(), NatureError> {
        let d = v.magnitude();
        if d <= gp::resolution() {
            return Err(NatureError::VectorWithNullMagnitude);
        }
        let xy = &v.coord;
        let a = xy.x() / d;
        let b = xy.y() / d;
        let m1 = 2.0 * a * b;
        let x = self.coord.x();
        let y = self.coord.y();
        self.coord.set_coord(
            ((2.0 * a * a) - 1.0) * x + m1 * y,
            m1 * x + ((2.0 * b * b) - 1.0) * y,
        );
        Ok(())
    }

    /// Returns a mirrored copy with respect to another vector.
    pub fn mirrored_vec(&self, v: &NVec2d) -> Result<NVec2d, NatureError> {
        let mut result = self.clone();
        result.mirror_vec(v)?;
        Ok(result)
    }

    /// Mirrors the vector with respect to an axis.
    pub fn mirror_ax2d(&mut self, a1: &NAx2d) {
        let xy = a1.direction().xy();
        let a = xy.x();
        let b = xy.y();
        let m1 = 2.0 * a * b;
        let x = self.coord.x();
        let y = self.coord.y();
        self.coord.set_coord(
            ((2.0 * a * a) - 1.0) * x + m1 * y,
            m1 * x + ((2.0 * b * b) - 1.0) * y,
        );
    }

    /// Returns a mirrored copy with respect to an axis.
    pub fn mirrored_ax2d(&self, a1: &NAx2d) -> NVec2d {
        let mut result = self.clone();
        result.mirror_ax2d(a1);
        result
    }

    /// Rotates the vector by an angle.
    pub fn rotate(&mut self, ang: f64) -> Result<(), NatureError> {
        let mut t = NTrsf2d::new();
        t.set_rotation(&NPnt2d::new(0.0, 0.0), ang)?;
        self.coord.multiply(&t.vectorial_part());
        Ok(())
    }

    /// Returns a rotated copy of the vector.
    pub fn rotated(&self, ang: f64) -> Result<NVec2d, NatureError> {
        let mut result = self.clone();
        result.rotate(ang)?;
        Ok(result)
    }

    /// Scales the vector by a factor.
    pub fn scale(&mut self, s: f64) {
        self.coord.multiply(s);
    }

    /// Returns a scaled copy of the vector.
    pub fn scaled(&self, s: f64) -> NVec2d {
        NVec2d {
            coord: self.coord.multiplied(s),
        }
    }

    /// Transforms the vector with a transformation.
    pub fn transform(&mut self, t: &NTrsf2d) {
        match t.form() {
            NTrsfForm::Identity | NTrsfForm::Translation => {}
            NTrsfForm::PntMirror => self.coord.reverse(),
            NTrsfForm::Scale => self.coord.multiply(t.scale_factor()),
            _ => self.coord.multiply(&t.vectorial_part()),
        }
    }

    /// Returns a transformed copy of the vector.
    pub fn transformed(&self, t: &NTrsf2d) -> NVec2d {
        let mut result = self.clone();
        result.transform(t);
        result
    }

    /// Dumps the vector as JSON.
    pub fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(
            out,
            "{} {{ \"type\": \"NVec2d\", \"coordinates\": [{}, {}] }}",
            indent,
            self.coord.x(),
            self.coord.y()
        ).unwrap();
    }

    /// Initializes the vector from JSON.
    pub fn init_from_json(&mut self, json: &str, pos: &mut usize) -> bool {
        let search = "\"coordinates\": [";
        if let Some(start) = json[*pos..].find(search) {
            *pos += start + search.len();
            let mut coords = [0.0; 2];
            let mut num_str = String::new();
            let mut idx = 0;
            while idx < 2 {
                let c = json.chars().nth(*pos).unwrap();
                if c.is_digit(10) || c == '.' || c == '-' {
                    num_str.push(c);
                } else if c == ',' || c == ']' {
                    if !num_str.is_empty() {
                        coords[idx] = num_str.parse().unwrap_or(0.0);
                        idx += 1;
                        num_str.clear();
                    }
                }
                *pos += 1;
                if c == ']' {
                    break;
                }
            }
            if idx == 2 {
                self.set_coords(coords[0], coords[1]);
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::NTrsfForm;

    #[test]
    fn test_new() {
        let v = NVec2d::new();
        assert_eq!(v.coords(), (0.0, 0.0));
    }

    #[test]
    fn test_new_from_coords() {
        let v = NVec2d::new_from_coords(1.0, 2.0);
        assert_eq!(v.coords(), (1.0, 2.0));
    }

    #[test]
    fn test_is_equal() {
        let v1 = NVec2d::new_from_coords(1.0, 0.0);
        let v2 = NVec2d::new_from_coords(1.0, 0.0);
        assert!(v1.is_equal(&v2, 1e-6, 1e-6));
        let v3 = NVec2d::new_from_coords(0.0, 0.0);
        assert!(v3.is_equal(&v3, 1e-6, 1e-6));
    }

    #[test]
    fn test_angle() {
        let v1 = NVec2d::new_from_coords(1.0, 0.0);
        let v2 = NVec2d::new_from_coords(0.0, 1.0);
        let angle = v1.angle(&v2).unwrap();
        assert!((angle - std::f64::consts::PI / 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_get_normal() {
        let v = NVec2d::new_from_coords(1.0, 0.0);
        let normal = v.get_normal();
        assert_eq!(normal.coords(), (0.0, -1.0));
    }

    #[test]
    fn test_mirror_vec() {
        let mut v = NVec2d::new_from_coords(1.0, 1.0);
        let mirror = NVec2d::new_from_coords(1.0, 0.0);
        v.mirror_vec(&mirror).unwrap();
        let (x, y) = v.coords();
        assert!((x - 1.0).abs() < 1e-6 && (y + 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_transform() {
        let mut v = NVec2d::new_from_coords(1.0, 0.0);
        let mut t = NTrsf2d::new();
        t.set_scale(&NPnt2d::new(0.0, 0.0), 2.0).unwrap();
        v.transform(&t);
        assert_eq!(v.coords(), (2.0, 0.0));
    }

    #[test]
    fn test_json() {
        let v = NVec2d::new_from_coords(1.0, 2.0);
        let mut output = Vec::new();
        v.dump_json(&mut output, 0);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"coordinates\": [1, 2]"));

        let mut v2 = NVec2d::new();
        let mut pos = 0;
        assert!(v2.init_from_json(&json, &mut pos));
        assert_eq!(v2.coords(), (1.0, 2.0));
    }
}