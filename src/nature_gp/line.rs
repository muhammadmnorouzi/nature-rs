use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::nature_gp::{NAx1, NAx2, NDir, NGP, NPoint3d, NTrsf, NVec, NXYZ};

// Trait to define the behavior of a line in 3D space
pub trait Lin {
    fn new() -> Self;
    fn new_with_axis(a1: &NAx1) -> Self;
    fn new_with_point_dir(p: &NPoint3d, v: &NDir) -> Self;
    fn reverse(&mut self);
    fn reversed(&self) -> Self
    where
        Self: Sized;
    fn set_direction(&mut self, v: &NDir);
    fn set_location(&mut self, p: &NPoint3d);
    fn set_position(&mut self, a1: &NAx1);
    fn direction(&self) -> &NDir;
    fn location(&self) -> &NPoint3d;
    fn position(&self) -> &NAx1;
    fn angle(&self, other: &Self) -> f64;
    fn contains(&self, p: &NPoint3d, linear_tolerance: f64) -> bool;
    fn distance_to_point(&self, p: &NPoint3d) -> f64;
    fn distance_to_line(&self, other: &Self) -> f64;
    fn square_distance_to_point(&self, p: &NPoint3d) -> f64;
    fn square_distance_to_line(&self, other: &Self) -> f64;
    fn normal(&self, p: &NPoint3d) -> Self
    where
        Self: Sized;
    fn mirror_point3d(&mut self, p: &NPoint3d);
    fn mirrored_point3d(&self, p: &NPoint3d) -> Self
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
    fn scale(&mut self, p: &NPoint3d, s: f64);
    fn scaled(&self, p: &NPoint3d, s: f64) -> Self
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
    fn translate_point3d(&mut self, p1: &NPoint3d, p2: &NPoint3d);
    fn translated_point3d(&self, p1: &NPoint3d, p2: &NPoint3d) -> Self
    where
        Self: Sized;
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
}

// Struct representing a line in 3D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NLin {
    pos: NAx1,
}

impl Lin for NLin {
    /// Creates a line corresponding to the Z-axis of the reference coordinate system.
    fn new() -> Self {
        NLin {
            pos: NAx1::new(
                &NPoint3d::new(0.0, 0.0, 0.0),
                &NDir::new(0.0, 0.0, 1.0).unwrap(),
            ),
        }
    }

    /// Creates a line defined by an axis.
    fn new_with_axis(a1: &NAx1) -> Self {
        NLin { pos: a1.clone() }
    }

    /// Creates a line passing through a point with a given direction.
    fn new_with_point_dir(p: &NPoint3d, v: &NDir) -> Self {
        NLin {
            pos: NAx1::new(p, v),
        }
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
    fn set_direction(&mut self, v: &NDir) {
        self.pos.set_direction(v);
    }

    /// Sets the location point (origin) of the line.
    fn set_location(&mut self, p: &NPoint3d) {
        self.pos.set_location(p);
    }

    /// Sets the position (axis) of the line.
    fn set_position(&mut self, a1: &NAx1) {
        self.pos = a1.clone();
    }

    /// Returns the direction of the line.
    fn direction(&self) -> &NDir {
        self.pos.direction()
    }

    /// Returns the location point (origin) of the line.
    fn location(&self) -> &NPoint3d {
        self.pos.location()
    }

    /// Returns the axis defining the line.
    fn position(&self) -> &NAx1 {
        &self.pos
    }

    /// Computes the angle between two lines in radians.
    fn angle(&self, other: &Self) -> f64 {
        self.pos.direction().angle(other.pos.direction())
    }

    /// Checks if the line contains a point within a linear tolerance.
    fn contains(&self, p: &NPoint3d, linear_tolerance: f64) -> bool {
        self.distance_to_point(p) <= linear_tolerance
    }

    /// Computes the distance from the line to a point.
    fn distance_to_point(&self, p: &NPoint3d) -> f64 {
        let mut coord = p.xyz();
        coord.subtract(&self.pos.location().xyz());
        coord.cross(&self.pos.direction().xyz());
        coord.modulus()
    }

    /// Computes the distance between two lines.
    fn distance_to_line(&self, other: &Self) -> f64 {
        if self.pos.is_parallel(&other.pos, NGP::resolution()) {
            other.distance_to_point(&self.pos.location())
        } else {
            let dir = self.pos.direction().crossed(&other.pos.direction());
            let mut d = NVec::new_from_points(&self.pos.location(), &other.pos.location())
                .dot(&NVec::from_dir(&dir));
            if d < 0.0 {
                d = -d;
            }
            d
        }
    }

    /// Computes the square distance from the line to a point.
    fn square_distance_to_point(&self, p: &NPoint3d) -> f64 {
        let loc = self.pos.location();
        let mut v = NVec::new(p.x() - loc.x(), p.y() - loc.y(), p.z() - loc.z());
        v.cross(&self.pos.direction());
        v.square_magnitude()
    }

    /// Computes the square distance between two lines.
    fn square_distance_to_line(&self, other: &Self) -> f64 {
        let d = self.distance_to_line(other);
        d * d
    }

    /// Computes a line normal to this line passing through a point.
    fn normal(&self, p: &NPoint3d) -> Self {
        let loc = self.pos.location();
        let v = NDir::new(p.x() - loc.x(), p.y() - loc.y(), p.z() - loc.z()).unwrap();
        let dir = self
            .pos
            .direction()
            .cross_crossed(&v, &self.pos.direction());
        NLin::new_with_point_dir(p, &dir)
    }

    /// Mirrors the line about a point.
    fn mirror_point3d(&mut self, p: &NPoint3d) {
        self.pos.mirror_point3d(p);
    }

    /// Returns a line mirrored about a point.
    fn mirrored_point3d(&self, p: &NPoint3d) -> Self {
        let mut l = self.clone();
        l.mirror_point3d(p);
        l
    }

    /// Mirrors the line about an axis.
    fn mirror_ax1(&mut self, a1: &NAx1) {
        self.pos.mirror_ax1(a1);
    }

    /// Returns a line mirrored about an axis.
    fn mirrored_ax1(&self, a1: &NAx1) -> Self {
        let mut l = self.clone();
        l.mirror_ax1(a1);
        l
    }

    /// Mirrors the line about a plane.
    fn mirror_ax2(&mut self, a2: &NAx2) {
        self.pos.mirror_ax2(a2);
    }

    /// Returns a line mirrored about a plane.
    fn mirrored_ax2(&self, a2: &NAx2) -> Self {
        let mut l = self.clone();
        l.mirror_ax2(a2);
        l
    }

    /// Rotates the line about an axis.
    fn rotate(&mut self, a1: &NAx1, ang: f64) {
        self.pos.rotate(a1, ang);
    }

    /// Returns a rotated line.
    fn rotated(&self, a1: &NAx1, ang: f64) -> Self {
        let mut l = self.clone();
        l.rotate(a1, ang);
        l
    }

    /// Scales the line about a point.
    fn scale(&mut self, p: &NPoint3d, s: f64) {
        self.pos.scale(p, s);
    }

    /// Returns a scaled line.
    fn scaled(&self, p: &NPoint3d, s: f64) -> Self {
        let mut l = self.clone();
        l.scale(p, s);
        l
    }

    /// Transforms the line with a transformation.
    fn transform(&mut self, t: &NTrsf) {
        self.pos.transform(t);
    }

    /// Returns a transformed line.
    fn transformed(&self, t: &NTrsf) -> Self {
        let mut l = self.clone();
        l.transform(t);
        l
    }

    /// Translates the line by a vector.
    fn translate_vec(&mut self, v: &NVec) {
        self.pos.translate_vec(v);
    }

    /// Returns a translated line by a vector.
    fn translated_vec(&self, v: &NVec) -> Self {
        let mut l = self.clone();
        l.translate_vec(v);
        l
    }

    /// Translates the line from one point to another.
    fn translate_point3d(&mut self, p1: &NPoint3d, p2: &NPoint3d) {
        self.pos.translate_point3d(p1, p2);
    }

    /// Returns a translated line from one point to another.
    fn translated_point3d(&self, p1: &NPoint3d, p2: &NPoint3d) -> Self {
        let mut l = self.clone();
        l.translate_point3d(p1, p2);
        l
    }

    /// Dumps the line as JSON.
    fn dump_json(&self, out: &mut dyn Write, _depth: i32) {
        writeln!(out, "{{\n  \"type\": \"NLin\",").unwrap();
        writeln!(
            out,
            "  \"pos\": {{ \"location\": [{}, {}, {}], \"direction\": [{}, {}, {}] }},",
            self.pos.location().x(),
            self.pos.location().y(),
            self.pos.location().z(),
            self.pos.direction().x(),
            self.pos.direction().y(),
            self.pos.direction().z()
        )
        .unwrap();
        writeln!(out, "}}").unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lin() -> NLin {
        NLin::new_with_point_dir(
            &NPoint3d::new(0.0, 0.0, 0.0),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
    }

    #[test]
    fn test_new() {
        let l = NLin::new();
        assert_eq!(l.direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
        assert_eq!(l.location(), &NPoint3d::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn test_new_with_point_dir() {
        let l = lin();
        assert_eq!(l.direction(), &NDir::new(1.0, 0.0, 0.0).unwrap());
        assert_eq!(l.location(), &NPoint3d::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn test_reverse() {
        let mut l = lin();
        l.reverse();
        assert_eq!(l.direction(), &NDir::new(-1.0, 0.0, 0.0).unwrap());
    }

    #[test]
    fn test_setters() {
        let mut l = lin();
        let p = NPoint3d::new(1.0, 2.0, 3.0);
        l.set_location(&p);
        assert_eq!(l.location(), &p);
        let v = NDir::new(0.0, 1.0, 0.0).unwrap();
        l.set_direction(&v);
        assert_eq!(l.direction(), &v);
    }

    #[test]
    fn test_distances() {
        let l = lin();
        let p = NPoint3d::new(0.0, 1.0, 0.0);
        assert!((l.distance_to_point(&p) - 1.0).abs() < 1e-9);
        assert!((l.square_distance_to_point(&p) - 1.0).abs() < 1e-9);

        let l2 = NLin::new_with_point_dir(
            &NPoint3d::new(0.0, 1.0, 0.0),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        );
        assert!((l.distance_to_line(&l2) - 1.0).abs() < 1e-9);
        assert!((l.square_distance_to_line(&l2) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_contains() {
        let l = lin();
        assert!(l.contains(&NPoint3d::new(1.0, 0.0, 0.0), 1e-9));
        assert!(!l.contains(&NPoint3d::new(1.0, 1.0, 0.0), 1e-9));
    }

    #[test]
    fn test_angle() {
        let l1 = lin();
        let l2 = NLin::new_with_point_dir(
            &NPoint3d::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 1.0, 0.0).unwrap(),
        );
        assert!((l1.angle(&l2) - std::f64::consts::PI / 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_normal() {
        let l = lin();
        let p = NPoint3d::new(0.0, 1.0, 0.0);
        let n = l.normal(&p);
        assert!(
            n.direction()
                .is_parallel(&NDir::new(0.0, 1.0, 0.0).unwrap(), 1e-9)
        );
        assert_eq!(n.location(), &p);
    }

    #[test]
    fn test_transformations() {
        let l = lin();
        let mut l_scaled = l.scaled(&NPoint3d::new(0.0, 0.0, 0.0), 2.0);
        assert_eq!(l_scaled.location(), &NPoint3d::new(0.0, 0.0, 0.0));

        let mut l_mirrored = l.mirrored_point3d(&NPoint3d::new(0.0, 1.0, 0.0));
        assert_eq!(l_mirrored.location(), &NPoint3d::new(0.0, 2.0, 0.0));
    }

    #[test]
    fn test_dump_json() {
        let l = lin();
        let mut output = Vec::new();
        l.dump_json(&mut output, -1);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NLin\""));
        assert!(json.contains("\"direction\": [1, 0, 0]"));
    }
}
