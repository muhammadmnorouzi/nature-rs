use std::io::Write;

use crate::nature_common::prelude::*;
use serde::{Deserialize, Serialize};

use super::prelude::*;

/// Represents a non-persistent 3D vector.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NVec {
    coord: NXYZ,
}

pub trait Vec {
    /// Creates a zero vector.
    fn new() -> Self;

    /// Creates a unitary vector from a direction.
    fn new_from_dir(dir: &NDir) -> Self;

    /// Creates a vector with the given coordinates.
    fn new_from_coords(x: f64, y: f64, z: f64) -> Self;

    /// Creates a vector from two points (P2 - P1).
    fn new_from_points(p1: &NPoint3d, p2: &NPoint3d) -> Self;

    /// Sets the coordinate at the given index (1=X, 2=Y, 3=Z).
    fn set_coord(&mut self, index: i32, value: f64) -> Result<(), NErrors>;

    /// Sets all coordinates.
    fn set_coords(&mut self, x: f64, y: f64, z: f64);

    /// Sets the X coordinate.
    fn set_x(&mut self, x: f64);

    /// Sets the Y coordinate.
    fn set_y(&mut self, y: f64);

    /// Sets the Z coordinate.
    fn set_z(&mut self, z: f64);

    /// Sets the coordinates from an NXYZ.
    fn set_xyz(&mut self, xyz: &NXYZ);

    /// Gets the coordinate at the given index (1=X, 2=Y, 3=Z).
    fn coord(&self, index: i32) -> Result<f64, NErrors>;

    /// Gets all coordinates.
    fn coords(&self) -> (f64, f64, f64);

    /// Gets the X coordinate.
    fn x(&self) -> f64;

    /// Gets the Y coordinate.
    fn y(&self) -> f64;

    /// Gets the Z coordinate.
    fn z(&self) -> f64;

    /// Gets the coordinates as an NXYZ.
    fn xyz(&self) -> NXYZ;

    /// Checks if two vectors are equal within tolerances.
    fn is_equal(&self, other: &NVec, linear_tolerance: f64, angular_tolerance: f64) -> bool;

    /// Checks if the vector is normal to another within angular tolerance.
    fn is_normal(&self, other: &NVec, angular_tolerance: f64) -> Result<bool, NErrors>;

    /// Checks if the vector is opposite to another within angular tolerance.
    fn is_opposite(&self, other: &NVec, angular_tolerance: f64) -> Result<bool, NErrors>;

    /// Checks if the vector is parallel to another within angular tolerance.
    fn is_parallel(&self, other: &NVec, angular_tolerance: f64) -> Result<bool, NErrors>;

    /// Computes the angle between two vectors (0 to PI radians).
    fn angle(&self, other: &NVec) -> Result<f64, NErrors>;

    /// Computes the signed angle between two vectors with respect to a reference vector (-PI to PI radians).
    fn angle_with_ref(&self, other: &NVec, v_ref: &NVec) -> Result<f64, NErrors>;

    /// Computes the magnitude of the vector.
    fn magnitude(&self) -> f64;

    /// Computes the square magnitude of the vector.
    fn square_magnitude(&self) -> f64;

    /// Adds another vector to this one.
    fn add(&mut self, other: &NVec);

    /// Returns the sum of two vectors.
    fn added(&self, other: &NVec) -> NVec;

    /// Subtracts another vector from this one.
    fn subtract(&mut self, other: &NVec);

    /// Returns the difference of two vectors.
    fn subtracted(&self, other: &NVec) -> NVec;

    /// Multiplies the vector by a scalar.
    fn multiply(&mut self, scalar: f64);

    /// Returns the vector multiplied by a scalar.
    fn multiplied(&self, scalar: f64) -> NVec;

    /// Divides the vector by a scalar.
    fn divide(&mut self, scalar: f64) -> Result<(), NErrors>;

    /// Returns the vector divided by a scalar.
    fn divided(&self, scalar: f64) -> Result<NVec, NErrors>;

    /// Computes the cross product with another vector.
    fn cross(&mut self, other: &NVec);

    /// Returns the cross product with another vector.
    fn crossed(&self, other: &NVec) -> NVec;

    /// Computes the magnitude of the cross product with another vector.
    fn cross_magnitude(&self, other: &NVec) -> f64;

    /// Computes the square magnitude of the cross product with another vector.
    fn cross_square_magnitude(&self, other: &NVec) -> f64;

    /// Computes the triple vector product (self = self ^ (v1 ^ v2)).
    fn cross_cross(&mut self, v1: &NVec, v2: &NVec);

    /// Returns the triple vector product (self ^ (v1 ^ v2)).
    fn cross_crossed(&self, v1: &NVec, v2: &NVec) -> NVec;

    /// Computes the dot product with another vector.
    fn dot(&self, other: &NVec) -> f64;

    /// Computes the triple scalar product (self * (v1 ^ v2)).
    fn dot_cross(&self, v1: &NVec, v2: &NVec) -> f64;

    /// Normalizes the vector.
    fn normalize(&mut self) -> Result<(), NErrors>;

    /// Returns a normalized copy of the vector.
    fn normalized(&self) -> Result<NVec, NErrors>;

    /// Reverses the direction of the vector.
    fn reverse(&mut self);

    /// Returns a reversed copy of the vector.
    fn reversed(&self) -> NVec;

    /// Sets the vector to a linear combination: a1*v1 + a2*v2 + a3*v3 + v4.
    fn set_linear_form(
        &mut self,
        a1: f64,
        v1: &NVec,
        a2: f64,
        v2: &NVec,
        a3: f64,
        v3: &NVec,
        v4: &NVec,
    );

    /// Sets the vector to a linear combination: a1*v1 + a2*v2 + a3*v3.
    fn set_linear_form_3(&mut self, a1: f64, v1: &NVec, a2: f64, v2: &NVec, a3: f64, v3: &NVec);

    /// Sets the vector to a linear combination: a1*v1 + a2*v2 + v3.
    fn set_linear_form_2_plus(&mut self, a1: f64, v1: &NVec, a2: f64, v2: &NVec, v3: &NVec);

    /// Sets the vector to a linear combination: a1*v1 + a2*v2.
    fn set_linear_form_2(&mut self, a1: f64, v1: &NVec, a2: f64, v2: &NVec);

    /// Sets the vector to a linear combination: a1*v1 + v2.
    fn set_linear_form_1_plus(&mut self, a1: f64, v1: &NVec, v2: &NVec);

    /// Sets the vector to a linear combination: v1 + v2.
    fn set_linear_form_sum(&mut self, v1: &NVec, v2: &NVec);

    /// Mirrors the vector with respect to another vector.
    fn mirror_vec(&mut self, v: &NVec) -> Result<(), NErrors>;

    /// Returns a mirrored copy with respect to another vector.
    fn mirrored_vec(&self, v: &NVec) -> Result<NVec, NErrors>;

    /// Mirrors the vector with respect to an axis.
    fn mirror_ax1(&mut self, a1: &NAx1);

    /// Returns a mirrored copy with respect to an axis.
    fn mirrored_ax1(&self, a1: &NAx1) -> NVec;

    /// Mirrors the vector with respect to a plane.
    fn mirror_ax2(&mut self, a2: &NAx2) -> Result<(), NErrors>;

    /// Returns a mirrored copy with respect to a plane.
    fn mirrored_ax2(&self, a2: &NAx2) -> Result<NVec, NErrors>;

    /// Rotates the vector around an axis by an angle.
    fn rotate(&mut self, a1: &NAx1, ang: f64) -> Result<(), NErrors>;

    /// Returns a rotated copy of the vector.
    fn rotated(&self, a1: &NAx1, ang: f64) -> Result<NVec, NErrors>;

    /// Scales the vector by a factor.
    fn scale(&mut self, s: f64);

    /// Returns a scaled copy of the vector.
    fn scaled(&self, s: f64) -> NVec;

    /// Transforms the vector with a transformation.
    fn transform(&mut self, t: &NTrsf);

    /// Returns a transformed copy of the vector.
    fn transformed(&self, t: &NTrsf) -> NVec;

    /// Dumps the vector as JSON.
    fn dump_json(&self, out: &mut dyn Write, depth: i32);

    /// Initializes the vector from JSON.
    fn init_from_json(&mut self, json: &str, pos: &mut usize) -> bool;
}

pub impl Vec for NVec {
    /// Creates a zero vector.
    fn new() -> Self {
        NVec {
            coord: NXYZ::new(0.0, 0.0, 0.0),
        }
    }

    /// Creates a unitary vector from a direction.
    fn new_from_dir(dir: &NDir) -> Self {
        NVec { coord: dir.xyz() }
    }

    /// Creates a vector with the given coordinates.
    fn new_from_coords(x: f64, y: f64, z: f64) -> Self {
        NVec {
            coord: NXYZ::new(x, y, z),
        }
    }

    /// Creates a vector from two points (P2 - P1).
    fn new_from_points(p1: &NPoint3d, p2: &NPoint3d) -> Self {
        NVec {
            coord: p2.xyz().subtracted(&p1.xyz()),
        }
    }

    /// Sets the coordinate at the given index (1=X, 2=Y, 3=Z).
    fn set_coord(&mut self, index: i32, value: f64) -> Result<(), NErrors> {
        match index {
            1 => self.coord.set_x(value),
            2 => self.coord.set_y(value),
            3 => self.coord.set_z(value),
            _ => return Err(NErrors::OutOfRange),
        }
        Ok(())
    }

    /// Sets all coordinates.
    fn set_coords(&mut self, x: f64, y: f64, z: f64) {
        self.coord.set_coord(x, y, z);
    }

    /// Sets the X coordinate.
    fn set_x(&mut self, x: f64) {
        self.coord.set_x(x);
    }

    /// Sets the Y coordinate.
    fn set_y(&mut self, y: f64) {
        self.coord.set_y(y);
    }

    /// Sets the Z coordinate.
    fn set_z(&mut self, z: f64) {
        self.coord.set_z(z);
    }

    /// Sets the coordinates from an NXYZ.
    fn set_xyz(&mut self, xyz: &NXYZ) {
        self.coord = xyz.clone();
    }

    /// Gets the coordinate at the given index (1=X, 2=Y, 3=Z).
    fn coord(&self, index: i32) -> Result<f64, NErrors> {
        match index {
            1 => Ok(self.coord.x()),
            2 => Ok(self.coord.y()),
            3 => Ok(self.coord.z()),
            _ => Err(NErrors::OutOfRange),
        }
    }

    /// Gets all coordinates.
    fn coords(&self) -> (f64, f64, f64) {
        (self.coord.x(), self.coord.y(), self.coord.z())
    }

    /// Gets the X coordinate.
    fn x(&self) -> f64 {
        self.coord.x()
    }

    /// Gets the Y coordinate.
    fn y(&self) -> f64 {
        self.coord.y()
    }

    /// Gets the Z coordinate.
    fn z(&self) -> f64 {
        self.coord.z()
    }

    /// Gets the coordinates as an NXYZ.
    fn xyz(&self) -> NXYZ {
        self.coord.clone()
    }

    /// Checks if two vectors are equal within tolerances.
    fn is_equal(&self, other: &NVec, linear_tolerance: f64, angular_tolerance: f64) -> bool {
        if self.magnitude() <= linear_tolerance || other.magnitude() <= linear_tolerance {
            (self.magnitude() - other.magnitude()).abs() <= linear_tolerance
        } else {
            (self.magnitude() - other.magnitude()).abs() <= linear_tolerance
                && self.angle(other).unwrap_or(std::f64::INFINITY) <= angular_tolerance
        }
    }

    /// Checks if the vector is normal to another within angular tolerance.
    fn is_normal(&self, other: &NVec, angular_tolerance: f64) -> Result<bool, NErrors> {
        let angle = (std::f64::consts::PI / 2.0 - self.angle(other)?).abs();
        Ok(angle <= angular_tolerance)
    }

    /// Checks if the vector is opposite to another within angular tolerance.
    fn is_opposite(&self, other: &NVec, angular_tolerance: f64) -> Result<bool, NErrors> {
        let angle = std::f64::consts::PI - self.angle(other)?;
        Ok(angle <= angular_tolerance)
    }

    /// Checks if the vector is parallel to another within angular tolerance.
    fn is_parallel(&self, other: &NVec, angular_tolerance: f64) -> Result<bool, NErrors> {
        let angle = self.angle(other)?;
        Ok(angle <= angular_tolerance || (std::f64::consts::PI - angle) <= angular_tolerance)
    }

    /// Computes the angle between two vectors (0 to PI radians).
    fn angle(&self, other: &NVec) -> Result<f64, NErrors> {
        if self.magnitude() <= gp::resolution() || other.magnitude() <= gp::resolution() {
            return Err(NErrors::VectorWithNullMagnitude);
        }
        let dir = NDir::new_from_xyz(&self.coord)?;
        dir.angle(&NDir::new_from_xyz(&other.coord))
    }

    /// Computes the signed angle between two vectors with respect to a reference vector (-PI to PI radians).
    fn angle_with_ref(&self, other: &NVec, v_ref: &NVec) -> Result<f64, NErrors> {
        if self.magnitude() <= gp::resolution()
            || other.magnitude() <= gp::resolution()
            || v_ref.magnitude() <= gp::resolution()
        {
            return Err(NErrors::VectorWithNullMagnitude);
        }
        let dir = NDir::new_from_xyz(&self.coord)?;
        let other_dir = NDir::new_from_xyz(&other.coord)?;
        dir.angle_with_ref(&other_dir, v_ref)
    }

    /// Computes the magnitude of the vector.
    fn magnitude(&self) -> f64 {
        self.coord.modulus()
    }

    /// Computes the square magnitude of the vector.
    fn square_magnitude(&self) -> f64 {
        self.coord.square_modulus()
    }

    /// Adds another vector to this one.
    fn add(&mut self, other: &NVec) {
        self.coord.add(&other.coord);
    }

    /// Returns the sum of two vectors.
    fn added(&self, other: &NVec) -> NVec {
        NVec {
            coord: self.coord.added(&other.coord),
        }
    }

    /// Subtracts another vector from this one.
    fn subtract(&mut self, other: &NVec) {
        self.coord.subtract(&other.coord);
    }

    /// Returns the difference of two vectors.
    fn subtracted(&self, other: &NVec) -> NVec {
        NVec {
            coord: self.coord.subtracted(&other.coord),
        }
    }

    /// Multiplies the vector by a scalar.
    fn multiply(&mut self, scalar: f64) {
        self.coord.multiply(scalar);
    }

    /// Returns the vector multiplied by a scalar.
    fn multiplied(&self, scalar: f64) -> NVec {
        NVec {
            coord: self.coord.multiplied(scalar),
        }
    }

    /// Divides the vector by a scalar.
    fn divide(&mut self, scalar: f64) -> Result<(), NErrors> {
        if scalar.abs() <= gp::resolution() {
            return Err(NErrors::InvalidConstructionParameters);
        }
        self.coord.divide(scalar);
        Ok(())
    }

    /// Returns the vector divided by a scalar.
    fn divided(&self, scalar: f64) -> Result<NVec, NErrors> {
        if scalar.abs() <= gp::resolution() {
            return Err(NErrors::InvalidConstructionParameters);
        }
        Ok(NVec {
            coord: self.coord.divided(scalar),
        })
    }

    /// Computes the cross product with another vector.
    fn cross(&mut self, other: &NVec) {
        self.coord.cross(&other.coord);
    }

    /// Returns the cross product with another vector.
    fn crossed(&self, other: &NVec) -> NVec {
        NVec {
            coord: self.coord.crossed(&other.coord),
        }
    }

    /// Computes the magnitude of the cross product with another vector.
    fn cross_magnitude(&self, other: &NVec) -> f64 {
        self.coord.cross_magnitude(&other.coord)
    }

    /// Computes the square magnitude of the cross product with another vector.
    fn cross_square_magnitude(&self, other: &NVec) -> f64 {
        self.coord.cross_square_magnitude(&other.coord)
    }

    /// Computes the triple vector product (self = self ^ (v1 ^ v2)).
    fn cross_cross(&mut self, v1: &NVec, v2: &NVec) {
        self.coord.cross_cross(&v1.coord, &v2.coord);
    }

    /// Returns the triple vector product (self ^ (v1 ^ v2)).
    fn cross_crossed(&self, v1: &NVec, v2: &NVec) -> NVec {
        NVec {
            coord: self.coord.cross_crossed(&v1.coord, &v2.coord),
        }
    }

    /// Computes the dot product with another vector.
    fn dot(&self, other: &NVec) -> f64 {
        self.coord.dot(&other.coord)
    }

    /// Computes the triple scalar product (self * (v1 ^ v2)).
    fn dot_cross(&self, v1: &NVec, v2: &NVec) -> f64 {
        self.coord.dot_cross(&v1.coord, &v2.coord)
    }

    /// Normalizes the vector.
    fn normalize(&mut self) -> Result<(), NErrors> {
        let mag = self.magnitude();
        if mag <= gp::resolution() {
            return Err(NErrors::VectorWithNullMagnitude);
        }
        self.coord.divide(mag);
        Ok(())
    }

    /// Returns a normalized copy of the vector.
    fn normalized(&self) -> Result<NVec, NErrors> {
        let mag = self.magnitude();
        if mag <= gp::resolution() {
            return Err(NErrors::VectorWithNullMagnitude);
        }
        Ok(NVec {
            coord: self.coord.divided(mag),
        })
    }

    /// Reverses the direction of the vector.
    fn reverse(&mut self) {
        self.coord.reverse();
    }

    /// Returns a reversed copy of the vector.
    fn reversed(&self) -> NVec {
        NVec {
            coord: self.coord.reversed(),
        }
    }

    /// Sets the vector to a linear combination: a1*v1 + a2*v2 + a3*v3 + v4.
    fn set_linear_form(
        &mut self,
        a1: f64,
        v1: &NVec,
        a2: f64,
        v2: &NVec,
        a3: f64,
        v3: &NVec,
        v4: &NVec,
    ) {
        self.coord
            .set_linear_form(a1, &v1.coord, a2, &v2.coord, a3, &v3.coord, &v4.coord);
    }

    /// Sets the vector to a linear combination: a1*v1 + a2*v2 + a3*v3.
    fn set_linear_form_3(&mut self, a1: f64, v1: &NVec, a2: f64, v2: &NVec, a3: f64, v3: &NVec) {
        self.coord
            .set_linear_form(a1, &v1.coord, a2, &v2.coord, a3, &v3.coord);
    }

    /// Sets the vector to a linear combination: a1*v1 + a2*v2 + v3.
    fn set_linear_form_2_plus(&mut self, a1: f64, v1: &NVec, a2: f64, v2: &NVec, v3: &NVec) {
        self.coord
            .set_linear_form(a1, &v1.coord, a2, &v2.coord, &v3.coord);
    }

    /// Sets the vector to a linear combination: a1*v1 + a2*v2.
    fn set_linear_form_2(&mut self, a1: f64, v1: &NVec, a2: f64, v2: &NVec) {
        self.coord.set_linear_form(a1, &v1.coord, a2, &v2.coord);
    }

    /// Sets the vector to a linear combination: a1*v1 + v2.
    fn set_linear_form_1_plus(&mut self, a1: f64, v1: &NVec, v2: &NVec) {
        self.coord.set_linear_form(a1, &v1.coord, &v2.coord);
    }

    /// Sets the vector to a linear combination: v1 + v2.
    fn set_linear_form_sum(&mut self, v1: &NVec, v2: &NVec) {
        self.coord.set_linear_form(&v1.coord, &v2.coord);
    }

    /// Mirrors the vector with respect to another vector.
    fn mirror_vec(&mut self, v: &NVec) -> Result<(), NErrors> {
        let d = v.magnitude();
        if d <= gp::resolution() {
            return Err(NErrors::VectorWithNullMagnitude);
        }
        let xyz = &v.coord;
        let a = xyz.x() / d;
        let b = xyz.y() / d;
        let c = xyz.z() / d;
        let m1 = 2.0 * a * b;
        let m2 = 2.0 * a * c;
        let m3 = 2.0 * b * c;
        let x = self.coord.x();
        let y = self.coord.y();
        let z = self.coord.z();
        self.coord.set_coord(
            ((2.0 * a * a) - 1.0) * x + m1 * y + m2 * z,
            m1 * x + ((2.0 * b * b) - 1.0) * y + m3 * z,
            m2 * x + m3 * y + ((2.0 * c * c) - 1.0) * z,
        );
        Ok(())
    }

    /// Returns a mirrored copy with respect to another vector.
    fn mirrored_vec(&self, v: &NVec) -> Result<NVec, NErrors> {
        let mut result = self.clone();
        result.mirror_vec(v)?;
        Ok(result)
    }

    /// Mirrors the vector with respect to an axis.
    fn mirror_ax1(&mut self, a1: &NAx1) {
        let v = a1.direction().xyz();
        let a = v.x();
        let b = v.y();
        let c = v.z();
        let m1 = 2.0 * a * b;
        let m2 = 2.0 * a * c;
        let m3 = 2.0 * b * c;
        let x = self.coord.x();
        let y = self.coord.y();
        let z = self.coord.z();
        self.coord.set_coord(
            ((2.0 * a * a) - 1.0) * x + m1 * y + m2 * z,
            m1 * x + ((2.0 * b * b) - 1.0) * y + m3 * z,
            m2 * x + m3 * y + ((2.0 * c * c) - 1.0) * z,
        );
    }

    /// Returns a mirrored copy with respect to an axis.
    fn mirrored_ax1(&self, a1: &NAx1) -> NVec {
        let mut result = self.clone();
        result.mirror_ax1(a1);
        result
    }

    /// Mirrors the vector with respect to a plane.
    fn mirror_ax2(&mut self, a2: &NAx2) -> Result<(), NErrors> {
        let z = a2.direction().xyz();
        let mut mir_xyz = z.crossed(&self.coord);
        if mir_xyz.modulus() <= gp::resolution() {
            self.coord.reverse();
        } else {
            let z_cross = z.crossed(&mir_xyz);
            self.mirror_vec(&NVec { coord: z_cross })?;
        }
        Ok(())
    }

    /// Returns a mirrored copy with respect to a plane.
    fn mirrored_ax2(&self, a2: &NAx2) -> Result<NVec, NErrors> {
        let mut result = self.clone();
        result.mirror_ax2(a2)?;
        Ok(result)
    }

    /// Rotates the vector around an axis by an angle.
    fn rotate(&mut self, a1: &NAx1, ang: f64) -> Result<(), NErrors> {
        let mut t = NTrsf::new();
        t.set_rotation(a1, ang)?;
        self.coord.multiply(&t.vectorial_part());
        Ok(())
    }

    /// Returns a rotated copy of the vector.
    fn rotated(&self, a1: &NAx1, ang: f64) -> Result<NVec, NErrors> {
        let mut result = self.clone();
        result.rotate(a1, ang)?;
        Ok(result)
    }

    /// Scales the vector by a factor.
    fn scale(&mut self, s: f64) {
        self.coord.multiply(s);
    }

    /// Returns a scaled copy of the vector.
    fn scaled(&self, s: f64) -> NVec {
        NVec {
            coord: self.coord.multiplied(s),
        }
    }

    /// Transforms the vector with a transformation.
    fn transform(&mut self, t: &NTrsf) {
        match t.form() {
            NTrsfForm::Identity | NTrsfForm::Translation => {}
            NTrsfForm::PntMirror => self.coord.reverse(),
            NTrsfForm::Scale => self.coord.multiply(t.scale_factor()),
            _ => self.coord.multiply(&t.vectorial_part()),
        }
    }

    /// Returns a transformed copy of the vector.
    fn transformed(&self, t: &NTrsf) -> NVec {
        let mut result = self.clone();
        result.transform(t);
        result
    }

    /// Dumps the vector as JSON.
    fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(
            out,
            "{} {{ \"type\": \"NVec\", \"coordinates\": [{}, {}, {}] }}",
            indent,
            self.coord.x(),
            self.coord.y(),
            self.coord.z()
        )
        .unwrap();
    }

    /// Initializes the vector from JSON.
    fn init_from_json(&mut self, json: &str, pos: &mut usize) -> bool {
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
    use crate::gp::NTrsfForm;

    #[test]
    fn test_new() {
        let v = NVec::new();
        assert_eq!(v.coords(), (0.0, 0.0, 0.0));
    }

    #[test]
    fn test_new_from_coords() {
        let v = NVec::new_from_coords(1.0, 2.0, 3.0);
        assert_eq!(v.coords(), (1.0, 2.0, 3.0));
    }

    #[test]
    fn test_is_equal() {
        let v1 = NVec::new_from_coords(1.0, 0.0, 0.0);
        let v2 = NVec::new_from_coords(1.0, 0.0, 0.0);
        assert!(v1.is_equal(&v2, 1e-6, 1e-6));
        let v3 = NVec::new_from_coords(0.0, 0.0, 0.0);
        assert!(v3.is_equal(&v3, 1e-6, 1e-6));
    }

    #[test]
    fn test_angle() {
        let v1 = NVec::new_from_coords(1.0, 0.0, 0.0);
        let v2 = NVec::new_from_coords(0.0, 1.0, 0.0);
        let angle = v1.angle(&v2).unwrap();
        assert!((angle - std::f64::consts::PI / 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_mirror_vec() {
        let mut v = NVec::new_from_coords(1.0, 1.0, 0.0);
        let mirror = NVec::new_from_coords(1.0, 0.0, 0.0);
        v.mirror_vec(&mirror).unwrap();
        let (x, y, z) = v.coords();
        assert!((x - 1.0).abs() < 1e-6 && (y + 1.0).abs() < 1e-6 && z.abs() < 1e-6);
    }

    #[test]
    fn test_transform() {
        let mut v = NVec::new_from_coords(1.0, 0.0, 0.0);
        let mut t = NTrsf::new();
        t.set_scale(&NPoint3d::new(0.0, 0.0, 0.0), 2.0).unwrap();
        v.transform(&t);
        assert_eq!(v.coords(), (2.0, 0.0, 0.0));
    }

    #[test]
    fn test_json() {
        let v = NVec::new_from_coords(1.0, 2.0, 3.0);
        let mut output = Vec::new();
        v.dump_json(&mut output, 0);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"coordinates\": [1, 2, 3]"));

        let mut v2 = NVec::new();
        let mut pos = 0;
        assert!(v2.init_from_json(&json, &mut pos));
        assert_eq!(v2.coords(), (1.0, 2.0, 3.0));
    }
}
