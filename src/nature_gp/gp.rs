use lazy_static::lazy_static;
use serde::{Deserialize, Serialize};
use std::io::Write;

use crate::gp::NatureError;

use super::xyz::NXYZ;

pub mod gp {
    /// Tolerance for floating-point comparisons.
    pub fn resolution() -> f64 {
        1e-12 // Consistent with vec.rs, vec2d.rs, etc.
    }

    pub fn resolution_f32() -> f32 {
        1e-6 // Consistent with vec2f.rs, vec3f.rs
    }
}

/// Represents a 3D point with double-precision coordinates.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NPnt {
    coord: NXYZ,
}

impl NPnt {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        NPnt {
            coord: NXYZ::new(x, y, z),
        }
    }

    pub fn x(&self) -> f64 {
        self.coord.x()
    }

    pub fn y(&self) -> f64 {
        self.coord.y()
    }

    pub fn z(&self) -> f64 {
        self.coord.z()
    }

    pub fn coords(&self) -> (f64, f64, f64) {
        (self.x(), self.y(), self.z())
    }

    pub fn xyz(&self) -> NXYZ {
        self.coord.clone()
    }

    pub fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(
            out,
            "{} {{ \"type\": \"NPnt\", \"coordinates\": [{}, {}, {}] }}",
            indent,
            self.x(),
            self.y(),
            self.z()
        ).unwrap();
    }

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
                self.coord = NXYZ::new(coords[0], coords[1], coords[2]);
                return true;
            }
        }
        false
    }
}

/// Represents a 3D direction (unit vector) with double-precision coordinates.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NDir {
    coord: NXYZ,
}

impl NDir {
    pub fn new(x: f64, y: f64, z: f64) -> Result<Self, NatureError> {
        let mag = (x * x + y * y + z * z).sqrt();
        raise_vector_with_null_magnitude_if!(
            mag <= gp::resolution(),
            "Direction cannot have zero magnitude"
        );
        Ok(NDir {
            coord: NXYZ::new(x / mag, y / mag, z / mag),
        })
    }

    pub fn x(&self) -> f64 {
        self.coord.x()
    }

    pub fn y(&self) -> f64 {
        self.coord.y()
    }

    pub fn z(&self) -> f64 {
        self.coord.z()
    }

    pub fn xyz(&self) -> NXYZ {
        self.coord.clone()
    }

    pub fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(
            out,
            "{} {{ \"type\": \"NDir\", \"coordinates\": [{}, {}, {}] }}",
            indent,
            self.x(),
            self.y(),
            self.z()
        ).unwrap();
    }

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
                match NDir::new(coords[0], coords[1], coords[2]) {
                    Ok(dir) => {
                        self.coord = dir.coord;
                        return true;
                    }
                    Err(_) => return false,
                }
            }
        }
        false
    }
}

/// Represents a 3D axis with a point and direction.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NAx1 {
    location: NPnt,
    direction: NDir,
}

impl NAx1 {
    pub fn new(location: NPnt, direction: NDir) -> Self {
        NAx1 { location, direction }
    }

    pub fn location(&self) -> &NPnt {
        &self.location
    }

    pub fn direction(&self) -> &NDir {
        &self.direction
    }

    pub fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        write!(out, "{} {{ \"type\": \"NAx1\", \"location\": ", indent).unwrap();
        self.location.dump_json(out, depth + 1);
        write!(out, ", \"direction\": ").unwrap();
        self.direction.dump_json(out, depth + 1);
        writeln!(out, "{}}}", indent).unwrap();
    }

    pub fn init_from_json(&mut self, json: &str, pos: &mut usize) -> bool {
        let search = "\"type\": \"NAx1\"";
        if let Some(start) = json[*pos..].find(search) {
            *pos += start + search.len();
            let mut location = NPnt::new(0.0, 0.0, 0.0);
            let mut direction = NDir::new(1.0, 0.0, 0.0).unwrap();
            if location.init_from_json(json, pos) && direction.init_from_json(json, pos) {
                self.location = location;
                self.direction = direction;
                return true;
            }
        }
        false
    }
}

/// Represents a 3D coordinate system with a point, main direction, and X direction.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NAx2 {
    location: NPnt,
    direction: NDir,
    x_direction: NDir,
}

impl NAx2 {
    pub fn new(location: NPnt, direction: NDir, x_direction: NDir) -> Result<Self, NatureError> {
        // Ensure direction and x_direction are not parallel
        let dot = direction.xyz().dot(&x_direction.xyz());
        if dot.abs() > 1.0 - gp::resolution() {
            return Err(NatureError::InvalidConstructionParameters);
        }
        Ok(NAx2 {
            location,
            direction,
            x_direction,
        })
    }

    pub fn location(&self) -> &NPnt {
        &self.location
    }

    pub fn direction(&self) -> &NDir {
        &self.direction
    }

    pub fn x_direction(&self) -> &NDir {
        &self.x_direction
    }

    pub fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        write!(out, "{} {{ \"type\": \"NAx2\", \"location\": ", indent).unwrap();
        self.location.dump_json(out, depth + 1);
        write!(out, ", \"direction\": ").unwrap();
        self.direction.dump_json(out, depth + 1);
        write!(out, ", \"x_direction\": ").unwrap();
        self.x_direction.dump_json(out, depth + 1);
        writeln!(out, "{}}}", indent).unwrap();
    }

    pub fn init_from_json(&mut self, json: &str, pos: &mut usize) -> bool {
        let search = "\"type\": \"NAx2\"";
        if let Some(start) = json[*pos..].find(search) {
            *pos += start + search.len();
            let mut location = NPnt::new(0.0, 0.0, 0.0);
            let mut direction = NDir::new(1.0, 0.0, 0.0).unwrap();
            let mut x_direction = NDir::new(0.0, 1.0, 0.0).unwrap();
            if location.init_from_json(json, pos)
                && direction.init_from_json(json, pos)
                && x_direction.init_from_json(json, pos)
            {
                match NAx2::new(location, direction, x_direction) {
                    Ok(ax2) => {
                        self.location = ax2.location;
                        self.direction = ax2.direction;
                        self.x_direction = ax2.x_direction;
                        return true;
                    }
                    Err(_) => return false,
                }
            }
        }
        false
    }
}

/// Represents a 2D point with double-precision coordinates.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NPnt2d {
    coord: NXY,
}

impl NPnt2d {
    pub fn new(x: f64, y: f64) -> Self {
        NPnt2d {
            coord: NXY::new(x, y),
        }
    }

    pub fn x(&self) -> f64 {
        self.coord.x()
    }

    pub fn y(&self) -> f64 {
        self.coord.y()
    }

    pub fn coords(&self) -> (f64, f64) {
        (self.x(), self.y())
    }

    pub fn xy(&self) -> NXY {
        self.coord.clone()
    }

    pub fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(
            out,
            "{} {{ \"type\": \"NPnt2d\", \"coordinates\": [{}, {}] }}",
            indent,
            self.x(),
            self.y()
        ).unwrap();
    }

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
                self.coord = NXY::new(coords[0], coords[1]);
                return true;
            }
        }
        false
    }
}

/// Represents a 2D direction (unit vector) with double-precision coordinates.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NDir2d {
    coord: NXY,
}

impl NDir2d {
    pub fn new(x: f64, y: f64) -> Result<Self, NatureError> {
        let mag = (x * x + y * y).sqrt();
        raise_vector_with_null_magnitude_if!(
            mag <= gp::resolution(),
            "Direction cannot have zero magnitude"
        );
        Ok(NDir2d {
            coord: NXY::new(x / mag, y / mag),
        })
    }

    pub fn x(&self) -> f64 {
        self.coord.x()
    }

    pub fn y(&self) -> f64 {
        self.coord.y()
    }

    pub fn xy(&self) -> NXY {
        self.coord.clone()
    }

    pub fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(
            out,
            "{} {{ \"type\": \"NDir2d\", \"coordinates\": [{}, {}] }}",
            indent,
            self.x(),
            self.y()
        ).unwrap();
    }

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
                match NDir2d::new(coords[0], coords[1]) {
                    Ok(dir) => {
                        self.coord = dir.coord;
                        return true;
                    }
                    Err(_) => return false,
                }
            }
        }
        false
    }
}

/// Represents a 2D axis with a point and direction.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NAx2d {
    location: NPnt2d,
    direction: NDir2d,
}

impl NAx2d {
    pub fn new(location: NPnt2d, direction: NDir2d) -> Self {
        NAx2d { location, direction }
    }

    pub fn location(&self) -> &NPnt2d {
        &self.location
    }

    pub fn direction(&self) -> &NDir2d {
        &self.direction
    }

    pub fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        write!(out, "{} {{ \"type\": \"NAx2d\", \"location\": ", indent).unwrap();
        self.location.dump_json(out, depth + 1);
        write!(out, ", \"direction\": ").unwrap();
        self.direction.dump_json(out, depth + 1);
        writeln!(out, "{}}}", indent).unwrap();
    }

    pub fn init_from_json(&mut self, json: &str, pos: &mut usize) -> bool {
        let search = "\"type\": \"NAx2d\"";
        if let Some(start) = json[*pos..].find(search) {
            *pos += start + search.len();
            let mut location = NPnt2d::new(0.0, 0.0);
            let mut direction = NDir2d::new(1.0, 0.0).unwrap();
            if location.init_from_json(json, pos) && direction.init_from_json(json, pos) {
                self.location = location;
                self.direction = direction;
                return true;
            }
        }
        false
    }
}

lazy_static! {
    /// Cartesian point at (0, 0, 0).
    pub static ref ORIGIN: NPnt = NPnt::new(0.0, 0.0, 0.0);

    /// Unit vector (1, 0, 0).
    pub static ref DX: NDir = NDir::new(1.0, 0.0, 0.0).unwrap();

    /// Unit vector (0, 1, 0).
    pub static ref DY: NDir = NDir::new(0.0, 1.0, 0.0).unwrap();

    /// Unit vector (0, 0, 1).
    pub static ref DZ: NDir = NDir::new(0.0, 0.0, 1.0).unwrap();

    /// Axis at origin with direction (1, 0, 0).
    pub static ref OX: NAx1 = NAx1::new(ORIGIN.clone(), DX.clone());

    /// Axis at origin with direction (0, 1, 0).
    pub static ref OY: NAx1 = NAx1::new(ORIGIN.clone(), DY.clone());

    /// Axis at origin with direction (0, 0, 1).
    pub static ref OZ: NAx1 = NAx1::new(ORIGIN.clone(), DZ.clone());

    /// Coordinate system at origin with main direction (0, 0, 1) and X direction (1, 0, 0).
    pub static ref XOY: NAx2 = NAx2::new(ORIGIN.clone(), DZ.clone(), DX.clone()).unwrap();

    /// Coordinate system at origin with main direction (0, 1, 0) and X direction (0, 0, 1).
    pub static ref ZOX: NAx2 = NAx2::new(ORIGIN.clone(), DY.clone(), DZ.clone()).unwrap();

    /// Coordinate system at origin with main direction (1, 0, 0) and X direction (0, 1, 0).
    pub static ref YOZ: NAx2 = NAx2::new(ORIGIN.clone(), DX.clone(), DY.clone()).unwrap();

    /// Cartesian point at (0, 0) in 2D.
    pub static ref ORIGIN2D: NPnt2d = NPnt2d::new(0.0, 0.0);

    /// Unit vector (1, 0) in 2D.
    pub static ref DX2D: NDir2d = NDir2d::new(1.0, 0.0).unwrap();

    /// Unit vector (0, 1) in 2D.
    pub static ref DY2D: NDir2d = NDir2d::new(0.0, 1.0).unwrap();

    /// Axis at origin with direction (1, 0) in 2D.
    pub static ref OX2D: NAx2d = NAx2d::new(ORIGIN2D.clone(), DX2D.clone());

    /// Axis at origin with direction (0, 1) in 2D.
    pub static ref OY2D: NAx2d = NAx2d::new(ORIGIN2D.clone(), DY2D.clone());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_origin() {
        assert_eq!(ORIGIN.coords(), (0.0, 0.0, 0.0));
    }

    #[test]
    fn test_dx() {
        assert_eq!(DX.xyz().coords(), (1.0, 0.0, 0.0));
    }

    #[test]
    fn test_ox() {
        assert_eq!(OX.location().coords(), (0.0, 0.0, 0.0));
        assert_eq!(OX.direction().xyz().coords(), (1.0, 0.0, 0.0));
    }

    #[test]
    fn test_xoy() {
        assert_eq!(XOY.location().coords(), (0.0, 0.0, 0.0));
        assert_eq!(XOY.direction().xyz().coords(), (0.0, 0.0, 1.0));
        assert_eq!(XOY.x_direction().xyz().coords(), (1.0, 0.0, 0.0));
    }

    #[test]
    fn test_origin2d() {
        assert_eq!(ORIGIN2D.coords(), (0.0, 0.0));
    }

    #[test]
    fn test_dx2d() {
        assert_eq!(DX2D.xy().coords(), (1.0, 0.0));
    }

    #[test]
    fn test_ox2d() {
        assert_eq!(OX2D.location().coords(), (0.0, 0.0));
        assert_eq!(OX2D.direction().xy().coords(), (1.0, 0.0));
    }

    #[test]
    fn test_json_serialization() {
        let mut output = Vec::new();
        ORIGIN.dump_json(&mut output, 0);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"coordinates\": [0, 0, 0]"));

        let mut p = NPnt::new(0.0, 0.0, 0.0);
        let mut pos = 0;
        assert!(p.init_from_json(&json, &mut pos));
        assert_eq!(p.coords(), (0.0, 0.0, 0.0));
    }

    #[test]
    fn test_ndir_invalid() {
        assert!(matches!(
            NDir::new(0.0, 0.0, 0.0),
            Err(NatureError::VectorWithNullMagnitude(_))
        ));
    }

    #[test]
    fn test_nax2_invalid() {
        let p = NPnt::new(0.0, 0.0, 0.0);
        let d = NDir::new(1.0, 0.0, 0.0).unwrap();
        assert!(matches!(
            NAx2::new(p.clone(), d.clone(), d),
            Err(NatureError::InvalidConstructionParameters)
        ));
    }
}