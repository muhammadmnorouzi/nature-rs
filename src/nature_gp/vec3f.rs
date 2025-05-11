use std::io::Write;

use serde::{Deserialize, Serialize};

use super::prelude::*;

mod gp {
    pub fn resolution_f32() -> f32 {
        1e-6 // Resolution for f32, consistent with vec2f.rs
    }
}

/// Represents a 3D coordinate triplet with single-precision floats.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NXYZf {
    x: f32,
    y: f32,
    z: f32,
}

impl NXYZf {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        NXYZf { x, y, z }
    }

    pub fn x(&self) -> f32 {
        self.x
    }

    pub fn y(&self) -> f32 {
        self.y
    }

    pub fn z(&self) -> f32 {
        self.z
    }

    pub fn set_x(&mut self, x: f32) {
        self.x = x;
    }

    pub fn set_y(&mut self, y: f32) {
        self.y = y;
    }

    pub fn set_z(&mut self, z: f32) {
        self.z = z;
    }

    pub fn set_coord(&mut self, x: f32, y: f32, z: f32) {
        self.x = x;
        self.y = y;
        self.z = z;
    }

    pub fn add(&mut self, other: &NXYZf) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }

    pub fn added(&self, other: &NXYZf) -> NXYZf {
        NXYZf {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }

    pub fn subtract(&mut self, other: &NXYZf) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }

    pub fn subtracted(&self, other: &NXYZf) -> NXYZf {
        NXYZf {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    pub fn multiply(&mut self, scalar: f32) {
        self.x *= scalar;
        self.y *= scalar;
        self.z *= scalar;
    }

    pub fn multiplied(&self, scalar: f32) -> NXYZf {
        NXYZf {
            x: self.x * scalar,
            y: self.y * scalar,
            z: self.z * scalar,
        }
    }

    pub fn divide(&mut self, scalar: f32) {
        self.x /= scalar;
        self.y /= scalar;
        self.z /= scalar;
    }

    pub fn divided(&self, scalar: f32) -> NXYZf {
        NXYZf {
            x: self.x / scalar,
            y: self.y / scalar,
            z: self.z / scalar,
        }
    }

    pub fn dot(&self, other: &NXYZf) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(&mut self, other: &NXYZf) {
        let x = self.y * other.z - self.z * other.y;
        let y = self.z * other.x - self.x * other.z;
        let z = self.x * other.y - self.y * other.x;
        self.x = x;
        self.y = y;
        self.z = z;
    }

    pub fn crossed(&self, other: &NXYZf) -> NXYZf {
        NXYZf {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    pub fn cross_magnitude(&self, other: &NXYZf) -> f32 {
        self.crossed(other).modulus()
    }

    pub fn cross_square_magnitude(&self, other: &NXYZf) -> f32 {
        self.crossed(other).square_modulus()
    }

    pub fn modulus(&self) -> f32 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn square_modulus(&self) -> f32 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub fn reverse(&mut self) {
        self.x = -self.x;
        self.y = -self.y;
        self.z = -self.z;
    }

    pub fn reversed(&self) -> NXYZf {
        NXYZf {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

/// Represents a non-persistent 3D vector with single-precision floats.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NVec3f {
    coord: NXYZf,
}

impl NVec3f {
    /// Creates a zero vector.
    pub fn new() -> Self {
        NVec3f {
            coord: NXYZf::new(0.0, 0.0, 0.0),
        }
    }

    /// Creates a vector with the given coordinates.
    pub fn new_from_coords(x: f32, y: f32, z: f32) -> Self {
        NVec3f {
            coord: NXYZf::new(x, y, z),
        }
    }

    /// Sets the coordinate at the given index (1=X, 2=Y, 3=Z).
    pub fn set_coord(&mut self, index: i32, value: f32) -> Result<(), NErrors> {
        match index {
            1 => self.coord.set_x(value),
            2 => self.coord.set_y(value),
            3 => self.coord.set_z(value),
            _ => return Err(NErrors::OutOfRange),
        }
        Ok(())
    }

    /// Sets all coordinates.
    pub fn set_coords(&mut self, x: f32, y: f32, z: f32) {
        self.coord.set_coord(x, y, z);
    }

    /// Sets the X coordinate.
    pub fn set_x(&mut self, x: f32) {
        self.coord.set_x(x);
    }

    /// Sets the Y coordinate.
    pub fn set_y(&mut self, y: f32) {
        self.coord.set_y(y);
    }

    /// Sets the Z coordinate.
    pub fn set_z(&mut self, z: f32) {
        self.coord.set_z(z);
    }

    /// Gets the coordinate at the given index (1=X, 2=Y, 3=Z).
    pub fn coord(&self, index: i32) -> Result<f32, NErrors> {
        match index {
            1 => Ok(self.coord.x()),
            2 => Ok(self.coord.y()),
            3 => Ok(self.coord.z()),
            _ => Err(NErrors::OutOfRange),
        }
    }

    /// Gets all coordinates.
    pub fn coords(&self) -> (f32, f32, f32) {
        (self.coord.x(), self.coord.y(), self.coord.z())
    }

    /// Gets the X coordinate.
    pub fn x(&self) -> f32 {
        self.coord.x()
    }

    /// Gets the Y coordinate.
    pub fn y(&self) -> f32 {
        self.coord.y()
    }

    /// Gets the Z coordinate.
    pub fn z(&self) -> f32 {
        self.coord.z()
    }

    /// Checks if two vectors are equal within a tolerance.
    pub fn is_equal(&self, other: &NVec3f, tolerance: f32) -> bool {
        let dx = self.x() - other.x();
        let dy = self.y() - other.y();
        let dz = self.z() - other.z();
        (dx * dx + dy * dy + dz * dz).sqrt() <= tolerance
    }

    /// Computes the angle between two vectors (0 to PI radians).
    pub fn angle(&self, other: &NVec3f) -> Result<f32, NErrors> {
        let norm = self.magnitude();
        let other_norm = other.magnitude();
        if norm <= gp::resolution_f32() || other_norm <= gp::resolution_f32() {
            return Err(NErrors::VectorWithNullMagnitude);
        }
        let dot = self.dot(other) / (norm * other_norm);
        // Clamp to [-1, 1] to handle floating-point errors
        let dot = dot.max(-1.0).min(1.0);
        Ok(dot.acos())
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
    pub fn add(&mut self, other: &NVec3f) {
        self.coord.add(&other.coord);
    }

    /// Returns the sum of two vectors.
    pub fn added(&self, other: &NVec3f) -> NVec3f {
        NVec3f {
            coord: self.coord.added(&other.coord),
        }
    }

    /// Subtracts another vector from this one.
    pub fn subtract(&mut self, other: &NVec3f) {
        self.coord.subtract(&other.coord);
    }

    /// Returns the difference of two vectors.
    pub fn subtracted(&self, other: &NVec3f) -> NVec3f {
        NVec3f {
            coord: self.coord.subtracted(&other.coord),
        }
    }

    /// Multiplies the vector by a scalar.
    pub fn multiply(&mut self, scalar: f32) {
        self.coord.multiply(scalar);
    }

    /// Returns the vector multiplied by a scalar.
    pub fn multiplied(&self, scalar: f32) -> NVec3f {
        NVec3f {
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
    pub fn divided(&self, scalar: f32) -> Result<NVec3f, NErrors> {
        if scalar.abs() <= gp::resolution_f32() {
            return Err(NErrors::InvalidConstructionParameters);
        }
        Ok(NVec3f {
            coord: self.coord.divided(scalar),
        })
    }

    /// Computes the cross product with another vector.
    pub fn cross(&mut self, other: &NVec3f) {
        self.coord.cross(&other.coord);
    }

    /// Returns the cross product with another vector.
    pub fn crossed(&self, other: &NVec3f) -> NVec3f {
        NVec3f {
            coord: self.coord.crossed(&other.coord),
        }
    }

    /// Computes the magnitude of the cross product with another vector.
    pub fn cross_magnitude(&self, other: &NVec3f) -> f32 {
        self.coord.cross_magnitude(&other.coord)
    }

    /// Computes the square magnitude of the cross product with another vector.
    pub fn cross_square_magnitude(&self, other: &NVec3f) -> f32 {
        self.coord.cross_square_magnitude(&other.coord)
    }

    /// Computes the dot product with another vector.
    pub fn dot(&self, other: &NVec3f) -> f32 {
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
    pub fn normalized(&self) -> Result<NVec3f, NErrors> {
        let mag = self.magnitude();
        if mag <= gp::resolution_f32() {
            return Err(NErrors::VectorWithNullMagnitude);
        }
        Ok(NVec3f {
            coord: self.coord.divided(mag),
        })
    }

    /// Reverses the direction of the vector.
    pub fn reverse(&mut self) {
        self.coord.reverse();
    }

    /// Returns a reversed copy of the vector.
    pub fn reversed(&self) -> NVec3f {
        NVec3f {
            coord: self.coord.reversed(),
        }
    }

    /// Dumps the vector as JSON.
    pub fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(
            out,
            "{} {{ \"type\": \"NVec3f\", \"coordinates\": [{}, {}, {}] }}",
            indent,
            self.coord.x(),
            self.coord.y(),
            self.coord.z()
        ).unwrap();
    }

    /// Initializes the vector from JSON.
    pub fn init_from_json(&mut self, json: &str, pos: &mut usize) -> bool {
        let search = "\"coordinates\": [";
        if let Some(start) = json[*pos..].find(search) {
            *pos += start + search.len();
            let mut coords = [0.0; 3];
            let mut num_str = String::new();
            let mut idx = 0;
            while idx < 3 {
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
            if idx == 3 {
                self.set_coords(coords[0], coords[1], coords[2]);
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
        let v = NVec3f::new();
        assert_eq!(v.coords(), (0.0, 0.0, 0.0));
    }

    #[test]
    fn test_new_from_coords() {
        let v = NVec3f::new_from_coords(1.0, 2.0, 3.0);
        assert_eq!(v.coords(), (1.0, 2.0, 3.0));
    }

    #[test]
    fn test_is_equal() {
        let v1 = NVec3f::new_from_coords(1.0, 0.0, 0.0);
        let v2 = NVec3f::new_from_coords(1.0, 0.0, 0.0);
        assert!(v1.is_equal(&v2, 1e-6));
        let v3 = NVec3f::new_from_coords(0.0, 0.0, 0.0);
        assert!(v3.is_equal(&v3, 1e-6));
    }

    #[test]
    fn test_angle() {
        let v1 = NVec3f::new_from_coords(1.0, 0.0, 0.0);
        let v2 = NVec3f::new_from_coords(0.0, 1.0, 0.0);
        let angle = v1.angle(&v2).unwrap();
        assert!((angle - std::f32::consts::PI / 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_operations() {
        let v1 = NVec3f::new_from_coords(1.0, 2.0, 3.0);
        let v2 = NVec3f::new_from_coords(4.0, 5.0, 6.0);
        assert_eq!(v1.added(&v2).coords(), (5.0, 7.0, 9.0));
        assert_eq!(v1.subtracted(&v2).coords(), (-3.0, -3.0, -3.0));
        assert_eq!(v1.multiplied(2.0).coords(), (2.0, 4.0, 6.0));
        assert_eq!(v1.dot(&v2), 32.0);
        assert_eq!(
            v1.crossed(&v2).coords(),
            (-3.0, 6.0, -3.0)
        );
    }

    #[test]
    fn test_normalize() {
        let mut v = NVec3f::new_from_coords(3.0, 4.0, 0.0);
        v.normalize().unwrap();
        let mag = v.magnitude();
        assert!((mag - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_json() {
        let v = NVec3f::new_from_coords(1.0, 2.0, 3.0);
        let mut output = Vec::new();
        v.dump_json(&mut output, 0);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"coordinates\": [1, 2, 3]"));

        let mut v2 = NVec3f::new();
        let mut pos = 0;
        assert!(v2.init_from_json(&json, &mut pos));
        assert_eq!(v2.coords(), (1.0, 2.0, 3.0));
    }
}