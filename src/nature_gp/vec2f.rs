use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::gp::NErrors;

mod gp {
    pub fn resolution_f32() -> f32 {
        1e-6 // Resolution for f32, less precise than f64
    }
}

/// Represents a 2D coordinate pair with single-precision floats.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NXYf {
    x: f32,
    y: f32,
}

impl NXYf {
    pub fn new(x: f32, y: f32) -> Self {
        NXYf { x, y }
    }

    pub fn x(&self) -> f32 {
        self.x
    }

    pub fn y(&self) -> f32 {
        self.y
    }

    pub fn set_x(&mut self, x: f32) {
        self.x = x;
    }

    pub fn set_y(&mut self, y: f32) {
        self.y = y;
    }

    pub fn set_coord(&mut self, x: f32, y: f32) {
        self.x = x;
        self.y = y;
    }

    pub fn add(&mut self, other: &NXYf) {
        self.x += other.x;
        self.y += other.y;
    }

    pub fn added(&self, other: &NXYf) -> NXYf {
        NXYf {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }

    pub fn subtract(&mut self, other: &NXYf) {
        self.x -= other.x;
        self.y -= other.y;
    }

    pub fn subtracted(&self, other: &NXYf) -> NXYf {
        NXYf {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }

    pub fn multiply(&mut self, scalar: f32) {
        self.x *= scalar;
        self.y *= scalar;
    }

    pub fn multiplied(&self, scalar: f32) -> NXYf {
        NXYf {
            x: self.x * scalar,
            y: self.y * scalar,
        }
    }

    pub fn divide(&mut self, scalar: f32) {
        self.x /= scalar;
        self.y /= scalar;
    }

    pub fn divided(&self, scalar: f32) -> NXYf {
        NXYf {
            x: self.x / scalar,
            y: self.y / scalar,
        }
    }

    pub fn dot(&self, other: &NXYf) -> f32 {
        self.x * other.x + self.y * other.y
    }

    pub fn crossed(&self, other: &NXYf) -> f32 {
        self.x * other.y - self.y * other.x
    }

    pub fn cross_magnitude(&self, other: &NXYf) -> f32 {
        self.crossed(other).abs()
    }

    pub fn cross_square_magnitude(&self, other: &NXYf) -> f32 {
        let cross = self.crossed(other);
        cross * cross
    }

    pub fn modulus(&self) -> f32 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    pub fn square_modulus(&self) -> f32 {
        self.x * self.x + self.y * self.y
    }

    pub fn reverse(&mut self) {
        self.x = -self.x;
        self.y = -self.y;
    }

    pub fn reversed(&self) -> NXYf {
        NXYf {
            x: -self.x,
            y: -self.y,
        }
    }
}

/// Represents a non-persistent 2D vector with single-precision floats.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NVec2f {
    coord: NXYf,
}

impl NVec2f {
    /// Creates a zero vector.
    pub fn new() -> Self {
        NVec2f {
            coord: NXYf::new(0.0, 0.0),
        }
    }

    /// Creates a vector with the given coordinates.
    pub fn new_from_coords(x: f32, y: f32) -> Self {
        NVec2f {
            coord: NXYf::new(x, y),
        }
    }

    /// Sets the coordinate at the given index (1=X, 2=Y).
    pub fn set_coord(&mut self, index: i32, value: f32) -> Result<(), NErrors> {
        match index {
            1 => self.coord.set_x(value),
            2 => self.coord.set_y(value),
            _ => return Err(NErrors::OutOfRange),
        }
        Ok(())
    }

    /// Sets all coordinates.
    pub fn set_coords(&mut self, x: f32, y: f32) {
        self.coord.set_coord(x, y);
    }

    /// Sets the X coordinate.
    pub fn set_x(&mut self, x: f32) {
        self.coord.set_x(x);
    }

    /// Sets the Y coordinate.
    pub fn set_y(&mut self, y: f32) {
        self.coord.set_y(y);
    }

    /// Gets the coordinate at the given index (1=X, 2=Y).
    pub fn coord(&self, index: i32) -> Result<f32, NErrors> {
        match index {
            1 => Ok(self.coord.x()),
            2 => Ok(self.coord.y()),
            _ => Err(NErrors::OutOfRange),
        }
    }

    /// Gets all coordinates.
    pub fn coords(&self) -> (f32, f32) {
        (self.coord.x(), self.coord.y())
    }

    /// Gets the X coordinate.
    pub fn x(&self) -> f32 {
        self.coord.x()
    }

    /// Gets the Y coordinate.
    pub fn y(&self) -> f32 {
        self.coord.y()
    }

    /// Checks if two vectors are equal within a tolerance.
    pub fn is_equal(&self, other: &NVec2f, tolerance: f32) -> bool {
        let dx = self.x() - other.x();
        let dy = self.y() - other.y();
        (dx * dx + dy * dy).sqrt() <= tolerance
    }

    /// Computes the angle between two vectors (-PI to PI radians).
    pub fn angle(&self, other: &NVec2f) -> Result<f32, NErrors> {
        let norm = self.magnitude();
        let other_norm = other.magnitude();
        if norm <= gp::resolution_f32() || other_norm <= gp::resolution_f32() {
            return Err(NErrors::VectorWithNullMagnitude);
        }
        let d = norm * other_norm;
        let cosinus = self.coord.dot(&other.coord) / d;
        let sinus = self.coord.crossed(&other.coord) / d;
        if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
            Ok(if sinus > 0.0 { cosinus.acos() } else { -cosinus.acos() })
        } else if cosinus > 0.0 {
            Ok(sinus.asin())
        } else if sinus > 0.0 {
            Ok(std::f32::consts::PI - sinus.asin())
        } else {
            Ok(-std::f32::consts::PI - sinus.asin())
        }
    }

    /// Computes the magnitude of the vector.
    pub fn magnitude(&self) -> f32 {
        self.coord.modulus()
    }

    /// Computes the square magnitude of the vector.
    pub fn square_magnitude(&self) -> f32 {
        self.coord.square_modulus()
    }

    /// Adds another vector to this one.
    pub fn add(&mut self, other: &NVec2f) {
        self.coord.add(&other.coord);
    }

    /// Returns the sum of two vectors.
    pub fn added(&self, other: &NVec2f) -> NVec2f {
        NVec2f {
            coord: self.coord.added(&other.coord),
        }
    }

    /// Subtracts another vector from this one.
    pub fn subtract(&mut self, other: &NVec2f) {
        self.coord.subtract(&other.coord);
    }

    /// Returns the difference of two vectors.
    pub fn subtracted(&self, other: &NVec2f) -> NVec2f {
        NVec2f {
            coord: self.coord.subtracted(&other.coord),
        }
    }

    /// Multiplies the vector by a scalar.
    pub fn multiply(&mut self, scalar: f32) {
        self.coord.multiply(scalar);
    }

    /// Returns the vector multiplied by a scalar.
    pub fn multiplied(&self, scalar: f32) -> NVec2f {
        NVec2f {
            coord: self.coord.multiplied(scalar),
        }
    }

    /// Divides the vector by a scalar.
    pub fn divide(&mut self, scalar: f32) -> Result<(), NErrors> {
        if scalar.abs() <= gp::resolution_f32() {
            return Err(NErrors::InvalidConstructionParameters);
        }
        self.coord.divide(scalar);
        Ok(())
    }

    /// Returns the vector divided by a scalar.
    pub fn divided(&self, scalar: f32) -> Result<NVec2f, NErrors> {
        if scalar.abs() <= gp::resolution_f32() {
            return Err(NErrors::InvalidConstructionParameters);
        }
        Ok(NVec2f {
            coord: self.coord.divided(scalar),
        })
    }

    /// Computes the cross product with another vector (returns a scalar).
    pub fn crossed(&self, other: &NVec2f) -> f32 {
        self.coord.crossed(&other.coord)
    }

    /// Computes the magnitude of the cross product with another vector.
    pub fn cross_magnitude(&self, other: &NVec2f) -> f32 {
        self.coord.cross_magnitude(&other.coord)
    }

    /// Computes the square magnitude of the cross product with another vector.
    pub fn cross_square_magnitude(&self, other: &NVec2f) -> f32 {
        self.coord.cross_square_magnitude(&other.coord)
    }

    /// Computes the dot product with another vector.
    pub fn dot(&self, other: &NVec2f) -> f32 {
        self.coord.dot(&other.coord)
    }

    /// Normalizes the vector.
    pub fn normalize(&mut self) -> Result<(), NErrors> {
        let mag = self.magnitude();
        if mag <= gp::resolution_f32() {
            return Err(NErrors::VectorWithNullMagnitude);
        }
        self.coord.divide(mag);
        Ok(())
    }

    /// Returns a normalized copy of the vector.
    pub fn normalized(&self) -> Result<NVec2f, NErrors> {
        let mag = self.magnitude();
        if mag <= gp::resolution_f32() {
            return Err(NErrors::VectorWithNullMagnitude);
        }
        Ok(NVec2f {
            coord: self.coord.divided(mag),
        })
    }

    /// Reverses the direction of the vector.
    pub fn reverse(&mut self) {
        self.coord.reverse();
    }

    /// Returns a reversed copy of the vector.
    pub fn reversed(&self) -> NVec2f {
        NVec2f {
            coord: self.coord.reversed(),
        }
    }

    /// Dumps the vector as JSON.
    pub fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(
            out,
            "{} {{ \"type\": \"NVec2f\", \"coordinates\": [{}, {}] }}",
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

    #[test]
    fn test_new() {
        let v = NVec2f::new();
        assert_eq!(v.coords(), (0.0, 0.0));
    }

    #[test]
    fn test_new_from_coords() {
        let v = NVec2f::new_from_coords(1.0, 2.0);
        assert_eq!(v.coords(), (1.0, 2.0));
    }

    #[test]
    fn test_is_equal() {
        let v1 = NVec2f::new_from_coords(1.0, 0.0);
        let v2 = NVec2f::new_from_coords(1.0, 0.0);
        assert!(v1.is_equal(&v2, 1e-6));
        let v3 = NVec2f::new_from_coords(0.0, 0.0);
        assert!(v3.is_equal(&v3, 1e-6));
    }

    #[test]
    fn test_angle() {
        let v1 = NVec2f::new_from_coords(1.0, 0.0);
        let v2 = NVec2f::new_from_coords(0.0, 1.0);
        let angle = v1.angle(&v2).unwrap();
        assert!((angle - std::f32::consts::PI / 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_operations() {
        let v1 = NVec2f::new_from_coords(1.0, 2.0);
        let v2 = NVec2f::new_from_coords(3.0, 4.0);
        assert_eq!(v1.added(&v2).coords(), (4.0, 6.0));
        assert_eq!(v1.subtracted(&v2).coords(), (-2.0, -2.0));
        assert_eq!(v1.multiplied(2.0).coords(), (2.0, 4.0));
        assert_eq!(v1.dot(&v2), 11.0);
        assert_eq!(v1.crossed(&v2), -2.0);
    }

    #[test]
    fn test_normalize() {
        let mut v = NVec2f::new_from_coords(3.0, 4.0);
        v.normalize().unwrap();
        let mag = v.magnitude();
        assert!((mag - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_json() {
        let v = NVec2f::new_from_coords(1.0, 2.0);
        let mut output = Vec::new();
        v.dump_json(&mut output, 0);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"coordinates\": [1, 2]"));

        let mut v2 = NVec2f::new();
        let mut pos = 0;
        assert!(v2.init_from_json(&json, &mut pos));
        assert_eq!(v2.coords(), (1.0, 2.0));
    }
}