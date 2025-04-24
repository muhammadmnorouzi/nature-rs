use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::{
    gp::{NAx2d, NDir2d, NGP, NPnt2d, NTrsf2d, NVec2d, NXY},
    nature_errors::NErrors,
};

// Trait to define the behavior of a line in 2D space
pub trait Lin2d {
    fn new() -> Self;
    fn new_with_axis(a: &NAx2d) -> Self;
    fn new_with_point_dir(p: &NPnt2d, v: &NDir2d) -> Self;
    fn new_with_coefficients(a: f64, b: f64, c: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn reverse(&mut self);
    fn reversed(&self) -> Self
    where
        Self: Sized;
    fn set_direction(&mut self, v: &NDir2d);
    fn set_location(&mut self, p: &NPnt2d);
    fn set_position(&mut self, a: &NAx2d);
    fn coefficients(&self) -> (f64, f64, f64);
    fn direction(&self) -> &NDir2d;
    fn location(&self) -> &NPnt2d;
    fn position(&self) -> &NAx2d;
    fn angle(&self, other: &Self) -> f64;
    fn contains(&self, p: &NPnt2d, linear_tolerance: f64) -> bool;
    fn distance_to_point(&self, p: &NPnt2d) -> f64;
    fn distance_to_line(&self, other: &Self) -> f64;
    fn square_distance_to_point(&self, p: &NPnt2d) -> f64;
    fn square_distance_to_line(&self, other: &Self) -> f64;
    fn normal(&self, p: &NPnt2d) -> Self
    where
        Self: Sized;
    fn mirror_pnt(&mut self, p: &NPnt2d);
    fn mirrored_pnt(&self, p: &NPnt2d) -> Self
    where
        Self: Sized;
    fn mirror_ax2d(&mut self, a: &NAx2d);
    fn mirrored_ax2d(&self, a: &NAx2d) -> Self
    where
        Self: Sized;
    fn rotate(&mut self, p: &NPnt2d, ang: f64);
    fn rotated(&self, p: &NPnt2d, ang: f64) -> Self
    where
        Self: Sized;
    fn scale(&mut self, p: &NPnt2d, s: f64);
    fn scaled(&self, p: &NPnt2d, s: f64) -> Self
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
    fn translate_pnts(&mut self, p1: &NPnt2d, p2: &NPnt2d);
    fn translated_pnts(&self, p1: &NPnt2d, p2: &NPnt2d) -> Self
    where
        Self: Sized;
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
}

// Struct representing a line in 2D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NLin2d {
    pos: NAx2d,
}

impl Lin2d for NLin2d {
    /// Creates a line corresponding to the X-axis of the reference coordinate system.
    fn new() -> Self {
        NLin2d {
            pos: NAx2d::new(&NPnt2d::new(0.0, 0.0), &NDir2d::new(1.0, 0.0).unwrap()),
        }
    }

    /// Creates a line defined by an axis.
    fn new_with_axis(a: &NAx2d) -> Self {
        NLin2d { pos: a.clone() }
    }

    /// Creates a line passing through a point with a given direction.
    fn new_with_point_dir(p: &NPnt2d, v: &NDir2d) -> Self {
        NLin2d {
            pos: NAx2d::new(p, v),
        }
    }

    /// Creates a line from the equation a*X + b*Y + c = 0.
    fn new_with_coefficients(a: f64, b: f64, c: f64) -> Result<Self, NErrors> {
        let norm2 = a * a + b * b;
        if norm2 <= NGP::resolution() {
            return Err(NErrors::ConstructionError);
        }
        let p = NPnt2d::new(-a * c / norm2, -b * c / norm2);
        let v = NDir2d::new(-b, a)?;
        Ok(NLin2d {
            pos: NAx2d::new(&p, &v),
        })
    }

    /// Reverses the direction of the line.
    fn reverse(&mut self) {
        self.pos.reverse();
    }

    /// Returns a line with reversed direction.
    fn reversed(&self) -> Self {
        let mut l = self.clone();
        l.reverse();
        l
    }

    /// Sets the direction of the line.
    fn set_direction(&mut self, v: &NDir2d) {
        self.pos.set_direction(v);
    }

    /// Sets the location point (origin) of the line.
    fn set_location(&mut self, p: &NPnt2d) {
        self.pos.set_location(p);
    }

    /// Sets the position (axis) of the line.
    fn set_position(&mut self, a: &NAx2d) {
        self.pos = a.clone();
    }

    /// Returns the normalized coefficients of the line: a*X + b*Y + c = 0.
    fn coefficients(&self) -> (f64, f64, f64) {
        let a = self.pos.direction().y();
        let b = -self.pos.direction().x();
        let c = -(a * self.pos.location().x() + b * self.pos.location().y());
        (a, b, c)
    }

    /// Returns the direction of the line.
    fn direction(&self) -> &NDir2d {
        self.pos.direction()
    }

    /// Returns the location point (origin) of the line.
    fn location(&self) -> &NPnt2d {
        self.pos.location()
    }

    /// Returns the axis defining the line.
    fn position(&self) -> &NAx2d {
        &self.pos
    }

    /// Computes the angle between two lines in radians.
    fn angle(&self, other: &Self) -> f64 {
        self.pos.direction().angle(other.pos.direction())
    }

    /// Checks if the line contains a point within a linear tolerance.
    fn contains(&self, p: &NPnt2d, linear_tolerance: f64) -> bool {
        self.distance_to_point(p) <= linear_tolerance
    }

    /// Computes the distance from the line to a point.
    fn distance_to_point(&self, p: &NPnt2d) -> f64 {
        let mut coord = p.xy();
        coord.subtract(&self.pos.location().xy());
        let val = coord.crossed(&self.pos.direction().xy());
        val.abs()
    }

    /// Computes the distance between two lines.
    fn distance_to_line(&self, other: &Self) -> f64 {
        if self.pos.is_parallel(&other.pos, NGP::resolution()) {
            other.distance_to_point(&self.pos.location())
        } else {
            0.0
        }
    }

    /// Computes the square distance from the line to a point.
    fn square_distance_to_point(&self, p: &NPnt2d) -> f64 {
        let mut coord = p.xy();
        coord.subtract(&self.pos.location().xy());
        let d = coord.crossed(&self.pos.direction().xy());
        d * d
    }

    /// Computes the square distance between two lines.
    fn square_distance_to_line(&self, other: &Self) -> f64 {
        if self.pos.is_parallel(&other.pos, NGP::resolution()) {
            other.square_distance_to_point(&self.pos.location())
        } else {
            0.0
        }
    }

    /// Computes a line normal to this line passing through a point.
    fn normal(&self, p: &NPnt2d) -> Self {
        let dir = NDir2d::new(-self.pos.direction().y(), self.pos.direction().x()).unwrap();
        NLin2d::new_with_point_dir(p, &dir)
    }

    /// Mirrors the line about a point.
    fn mirror_pnt(&mut self, p: &NPnt2d) {
        self.pos.mirror_pnt(p);
    }

    /// Returns a line mirrored about a point.
    fn mirrored_pnt(&self, p: &NPnt2d) -> Self {
        let mut l = self.clone();
        l.mirror_pnt(p);
        l
    }

    /// Mirrors the line about an axis.
    fn mirror_ax2d(&mut self, a: &NAx2d) {
        self.pos.mirror_ax2d(a);
    }

    /// Returns a line mirrored about an axis.
    fn mirrored_ax2d(&self, a: &NAx2d) -> Self {
        let mut l = self.clone();
        l.mirror_ax2d(a);
        l
    }

    /// Rotates the line about a point.
    fn rotate(&mut self, p: &NPnt2d, ang: f64) {
        self.pos.rotate(p, ang);
    }

    /// Returns a rotated line.
    fn rotated(&self, p: &NPnt2d, ang: f64) -> Self {
        let mut l = self.clone();
        l.rotate(p, ang);
        l
    }

    /// Scales the line about a point.
    fn scale(&mut self, p: &NPnt2d, s: f64) {
        self.pos.scale(p, s);
    }

    /// Returns a scaled line.
    fn scaled(&self, p: &NPnt2d, s: f64) -> Self {
        let mut l = self.clone();
        l.scale(p, s);
        l
    }

    /// Transforms the line with a transformation.
    fn transform(&mut self, t: &NTrsf2d) {
        self.pos.transform(t);
    }

    /// Returns a transformed line.
    fn transformed(&self, t: &NTrsf2d) -> Self {
        let mut l = self.clone();
        l.transform(t);
        l
    }

    /// Translates the line by a vector.
    fn translate_vec(&mut self, v: &NVec2d) {
        self.pos.translate_vec(v);
    }

    /// Returns a translated line by a vector.
    fn translated_vec(&self, v: &NVec2d) -> Self {
        let mut l = self.clone();
        l.translate_vec(v);
        l
    }

    /// Translates the line from one point to another.
    fn translate_pnts(&mut self, p1: &NPnt2d, p2: &NPnt2d) {
        self.pos.translate_pnts(p1, p2);
    }

    /// Returns a translated line from one point to another.
    fn translated_pnts(&self, p1: &NPnt2d, p2: &NPnt2d) -> Self {
        let mut l = self.clone();
        l.translate_pnts(p1, p2);
        l
    }

    /// Dumps the line as JSON.
    fn dump_json(&self, out: &mut dyn Write, _depth: i32) {
        writeln!(out, "{{\n  \"type\": \"NLin2d\",").unwrap();
        writeln!(
            out,
            "  \"pos\": {{ \"location\": [{}, {}], \"direction\": [{}, {}] }},",
            self.pos.location().x(),
            self.pos.location().y(),
            self.pos.direction().x(),
            self.pos.direction().y()
        )
        .unwrap();
        writeln!(out, "}}").unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lin2d() -> NLin2d {
        NLin2d::new_with_point_dir(&NPnt2d::new(0.0, 0.0), &NDir2d::new(1.0, 0.0).unwrap())
    }

    #[test]
    fn test_new() {
        let l = NLin2d::new();
        assert_eq!(l.direction(), &NDir2d::new(1.0, 0.0).unwrap());
        assert_eq!(l.location(), &NPnt2d::new(0.0, 0.0));
    }

    #[test]
    fn test_new_with_point_dir() {
        let l = lin2d();
        assert_eq!(l.direction(), &NDir2d::new(1.0, 0.0).unwrap());
        assert_eq!(l.location(), &NPnt2d::new(0.0, 0.0));
    }

    #[test]
    fn test_new_with_coefficients() {
        let l = NLin2d::new_with_coefficients(1.0, 0.0, 0.0).unwrap();
        assert!(
            l.direction()
                .is_parallel(&NDir2d::new(0.0, 1.0).unwrap(), 1e-9)
        );
        assert_eq!(l.location(), &NPnt2d::new(0.0, 0.0));
        assert!(matches!(
            NLin2d::new_with_coefficients(0.0, 0.0, 1.0),
            Err(NErrors::ConstructionError)
        ));
    }

    #[test]
    fn test_reverse() {
        let mut l = lin2d();
        l.reverse();
        assert_eq!(l.direction(), &NDir2d::new(-1.0, 0.0).unwrap());
    }

    #[test]
    fn test_setters() {
        let mut l = lin2d();
        let p = NPnt2d::new(1.0, 2.0);
        l.set_location(&p);
        assert_eq!(l.location(), &p);
        let v = NDir2d::new(0.0, 1.0).unwrap();
        l.set_direction(&v);
        assert_eq!(l.direction(), &v);
    }

    #[test]
    fn test_coefficients() {
        let l = lin2d();
        let (a, b, c) = l.coefficients();
        assert!((a - 0.0).abs() < 1e-9);
        assert!((b - -1.0).abs() < 1e-9);
        assert!((c - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_distances() {
        let l = lin2d();
        let p = NPnt2d::new(0.0, 1.0);
        assert!((l.distance_to_point(&p) - 1.0).abs() < 1e-9);
        assert!((l.square_distance_to_point(&p) - 1.0).abs() < 1e-9);

        let l2 =
            NLin2d::new_with_point_dir(&NPnt2d::new(0.0, 1.0), &NDir2d::new(1.0, 0.0).unwrap());
        assert!((l.distance_to_line(&l2) - 1.0).abs() < 1e-9);
        assert!((l.square_distance_to_line(&l2) - 1.0).abs() < 1e-9);

        let l3 =
            NLin2d::new_with_point_dir(&NPnt2d::new(0.0, 1.0), &NDir2d::new(0.0, 1.0).unwrap());
        assert!((l.distance_to_line(&l3) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_contains() {
        let l = lin2d();
        assert!(l.contains(&NPnt2d::new(1.0, 0.0), 1e-9));
        assert!(!l.contains(&NPnt2d::new(1.0, 1.0), 1e-9));
    }

    #[test]
    fn test_angle() {
        let l1 = lin2d();
        let l2 =
            NLin2d::new_with_point_dir(&NPnt2d::new(0.0, 0.0), &NDir2d::new(0.0, 1.0).unwrap());
        assert!((l1.angle(&l2) - std::f64::consts::PI / 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_normal() {
        let l = lin2d();
        let p = NPnt2d::new(0.0, 1.0);
        let n = l.normal(&p);
        assert!(
            n.direction()
                .is_parallel(&NDir2d::new(0.0, 1.0).unwrap(), 1e-9)
        );
        assert_eq!(n.location(), &p);
    }

    #[test]
    fn test_transformations() {
        let l = lin2d();
        let mut l_scaled = l.scaled(&NPnt2d::new(0.0, 0.0), 2.0);
        assert_eq!(l_scaled.location(), &NPnt2d::new(0.0, 0.0));

        let mut l_mirrored = l.mirrored_pnt(&NPnt2d::new(0.0, 1.0));
        assert_eq!(l_mirrored.location(), &NPnt2d::new(0.0, 2.0));
    }

    #[test]
    fn test_dump_json() {
        let l = lin2d();
        let mut output = Vec::new();
        l.dump_json(&mut output, -1);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NLin2d\""));
        assert!(json.contains("\"direction\": [1, 0]"));
    }
}
