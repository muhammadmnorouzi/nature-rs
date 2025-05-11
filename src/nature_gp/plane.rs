use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::{
    gp::{NAx1, NAx2, NAx3, NDir, NGP, NLin, NPoint3d, NTrsf, NVec},
    nature_errors::NErrors,
};

// Trait to define the behavior of a plane in 3D space
pub trait Pln {
    fn new() -> Self;
    fn new_with_ax3(a3: &NAx3) -> Self;
    fn new_with_point_dir(p: &NPoint3d, v: &NDir) -> Result<Self, NErrors>;
    fn new_with_coefficients(a: f64, b: f64, c: f64, d: f64) -> Result<Self, NErrors>;
    fn coefficients(&self) -> (f64, f64, f64, f64);
    fn set_axis(&mut self, a1: &NAx1) -> Result<(), NErrors>;
    fn set_location(&mut self, loc: &NPoint3d);
    fn set_position(&mut self, a3: &NAx3);
    fn u_reverse(&mut self);
    fn v_reverse(&mut self);
    fn direct(&self) -> bool;
    fn axis(&self) -> NAx1;
    fn location(&self) -> NPoint3d;
    fn position(&self) -> NAx3;
    fn distance_pnt(&self, p: &NPoint3d) -> f64;
    fn distance_lin(&self, l: &NLin) -> f64;
    fn distance_pln(&self, other: &Self) -> f64;
    fn square_distance_pnt(&self, p: &NPoint3d) -> f64;
    fn square_distance_lin(&self, l: &NLin) -> f64;
    fn square_distance_pln(&self, other: &Self) -> f64;
    fn x_axis(&self) -> NAx1;
    fn y_axis(&self) -> NAx1;
    fn contains_pnt(&self, p: &NPoint3d, linear_tolerance: f64) -> bool;
    fn contains_lin(&self, l: &NLin, linear_tolerance: f64, angular_tolerance: f64) -> bool;
    fn mirror_pnt(&mut self, p: &NPoint3d);
    fn mirrored_pnt(&self, p: &NPoint3d) -> Self
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

// Struct representing a plane in 3D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NPln {
    pos: NAx3,
}

impl Pln for NPln {
    /// Creates a plane coincident with the OXY plane of the reference coordinate system.
    fn new() -> Self {
        NPln {
            pos: NAx3::new(
                &NPoint3d::new(0.0, 0.0, 0.0),
                &NDir::new(0.0, 0.0, 1.0),
                &NDir::new(1.0, 0.0, 0.0),
            )
            .unwrap(),
        }
    }

    /// Creates a plane with the specified coordinate system.
    fn new_with_ax3(a3: &NAx3) -> Self {
        NPln { pos: a3.clone() }
    }

    /// Creates a plane with a location point and normal direction.
    fn new_with_point_dir(p: &NPoint3d, v: &NDir) -> Result<Self, NErrors> {
        let a = v.x();
        let b = v.y();
        let c = v.z();
        let a_abs = a.abs();
        let b_abs = b.abs();
        let c_abs = c.abs();

        let pos = if b_abs <= a_abs && b_abs <= c_abs {
            if a_abs > c_abs {
                NAx3::new(p, v, &NDir::new(-c, 0.0, a)).map_err(|_| NErrors::ConstructionError)?
            } else {
                NAx3::new(p, v, &NDir::new(c, 0.0, -a)).map_err(|_| NErrors::ConstructionError)?
            }
        } else if a_abs <= b_abs && a_abs <= c_abs {
            if b_abs > c_abs {
                NAx3::new(p, v, &NDir::new(0.0, -c, b)).map_err(|_| NErrors::ConstructionError)?
            } else {
                NAx3::new(p, v, &NDir::new(0.0, c, -b)).map_err(|_| NErrors::ConstructionError)?
            }
        } else {
            if a_abs > b_abs {
                NAx3::new(p, v, &NDir::new(-b, a, 0.0)).map_err(|_| NErrors::ConstructionError)?
            } else {
                NAx3::new(p, v, &NDir::new(b, -a, 0.0)).map_err(|_| NErrors::ConstructionError)?
            }
        };
        Ok(NPln { pos })
    }

    /// Creates a plane from its cartesian equation: a*X + b*Y + c*Z + d = 0.
    fn new_with_coefficients(a: f64, b: f64, c: f64, d: f64) -> Result<Self, NErrors> {
        let norm = (a * a + b * b + c * c).sqrt();
        if norm <= NGP::resolution() {
            return Err(NErrors::ConstructionError);
        }
        let a_abs = a.abs();
        let b_abs = b.abs();
        let c_abs = c.abs();
        let pos = if b_abs <= a_abs && b_abs <= c_abs {
            if a_abs > c_abs {
                NAx3::new(
                    &NPoint3d::new(-d / a, 0.0, 0.0),
                    &NDir::new(a, b, c),
                    &NDir::new(-c, 0.0, a),
                )
            } else {
                NAx3::new(
                    &NPoint3d::new(0.0, 0.0, -d / c),
                    &NDir::new(a, b, c),
                    &NDir::new(c, 0.0, -a),
                )
            }
        } else if a_abs <= b_abs && a_abs <= c_abs {
            if b_abs > c_abs {
                NAx3::new(
                    &NPoint3d::new(0.0, -d / b, 0.0),
                    &NDir::new(a, b, c),
                    &NDir::new(0.0, -c, b),
                )
            } else {
                NAx3::new(
                    &NPoint3d::new(0.0, 0.0, -d / c),
                    &NDir::new(a, b, c),
                    &NDir::new(0.0, c, -b),
                )
            }
        } else {
            if a_abs > b_abs {
                NAx3::new(
                    &NPoint3d::new(-d / a, 0.0, 0.0),
                    &NDir::new(a, b, c),
                    &NDir::new(-b, a, 0.0),
                )
            } else {
                NAx3::new(
                    &NPoint3d::new(0.0, -d / b, 0.0),
                    &NDir::new(a, b, c),
                    &NDir::new(b, -a, 0.0),
                )
            }
        }
        .map_err(|_| NErrors::ConstructionError)?;
        Ok(NPln { pos })
    }

    /// Returns the coefficients of the plane's cartesian equation: a*X + b*Y + c*Z + d = 0.
    fn coefficients(&self) -> (f64, f64, f64, f64) {
        let dir = self.pos.direction();
        let (a, b, c) = if self.pos.direct() {
            (dir.x(), dir.y(), dir.z())
        } else {
            (-dir.x(), -dir.y(), -dir.z())
        };
        let p = self.pos.location();
        let d = -(a * p.x() + b * p.y() + c * p.z());
        (a, b, c, d)
    }

    /// Sets the main axis of the plane.
    fn set_axis(&mut self, a1: &NAx1) -> Result<(), NErrors> {
        self.pos.set_axis(a1)
    }

    /// Changes the origin of the plane.
    fn set_location(&mut self, loc: &NPoint3d) {
        self.pos.set_location(loc);
    }

    /// Changes the local coordinate system of the plane.
    fn set_position(&mut self, a3: &NAx3) {
        self.pos = a3.clone();
    }

    /// Reverses the U parametrization (X-axis).
    fn u_reverse(&mut self) {
        self.pos.x_reverse();
    }

    /// Reverses the V parametrization (Y-axis).
    fn v_reverse(&mut self) {
        self.pos.y_reverse();
    }

    /// Returns true if the coordinate system is right-handed.
    fn direct(&self) -> bool {
        self.pos.direct()
    }

    /// Returns the plane's normal axis.
    fn axis(&self) -> NAx1 {
        self.pos.axis()
    }

    /// Returns the plane's location (origin).
    fn location(&self) -> NPoint3d {
        self.pos.location()
    }

    /// Returns the local coordinate system of the plane.
    fn position(&self) -> NAx3 {
        self.pos.clone()
    }

    /// Computes the distance between the plane and a point.
    fn distance_pnt(&self, p: &NPoint3d) -> f64 {
        let loc = self.pos.location();
        let dir = self.pos.direction();
        let d =
            dir.x() * (p.x() - loc.x()) + dir.y() * (p.y() - loc.y()) + dir.z() * (p.z() - loc.z());
        d.abs()
    }

    /// Computes the distance between the plane and a line.
    fn distance_lin(&self, l: &NLin) -> f64 {
        if self
            .pos
            .direction()
            .is_normal(&l.direction(), NGP::resolution())
        {
            let p = l.location();
            let loc = self.pos.location();
            let dir = self.pos.direction();
            let d = dir.x() * (p.x() - loc.x())
                + dir.y() * (p.y() - loc.y())
                + dir.z() * (p.z() - loc.z());
            d.abs()
        } else {
            0.0
        }
    }

    /// Computes the distance between two planes.
    fn distance_pln(&self, other: &Self) -> f64 {
        if self
            .pos
            .direction()
            .is_parallel(&other.pos.direction(), NGP::resolution())
        {
            let p = other.pos.location();
            let loc = self.pos.location();
            let dir = self.pos.direction();
            let d = dir.x() * (p.x() - loc.x())
                + dir.y() * (p.y() - loc.y())
                + dir.z() * (p.z() - loc.z());
            d.abs()
        } else {
            0.0
        }
    }

    /// Computes the square distance between the plane and a point.
    fn square_distance_pnt(&self, p: &NPoint3d) -> f64 {
        let d = self.distance_pnt(p);
        d * d
    }

    /// Computes the square distance between the plane and a line.
    fn square_distance_lin(&self, l: &NLin) -> f64 {
        let d = self.distance_lin(l);
        d * d
    }

    /// Computes the square distance between two planes.
    fn square_distance_pln(&self, other: &Self) -> f64 {
        let d = self.distance_pln(other);
        d * d
    }

    /// Returns the X-axis of the plane.
    fn x_axis(&self) -> NAx1 {
        NAx1::new(&self.pos.location(), &self.pos.x_direction())
    }

    /// Returns the Y-axis of the plane.
    fn y_axis(&self) -> NAx1 {
        NAx1::new(&self.pos.location(), &self.pos.y_direction())
    }

    /// Returns true if the plane contains the point within the linear tolerance.
    fn contains_pnt(&self, p: &NPoint3d, linear_tolerance: f64) -> bool {
        self.distance_pnt(p) <= linear_tolerance
    }

    /// Returns true if the plane contains the line within the specified tolerances.
    fn contains_lin(&self, l: &NLin, linear_tolerance: f64, angular_tolerance: f64) -> bool {
        self.contains_pnt(&l.location(), linear_tolerance)
            && self
                .pos
                .direction()
                .is_normal(&l.direction(), angular_tolerance)
    }

    /// Mirrors the plane with respect to a point.
    fn mirror_pnt(&mut self, p: &NPoint3d) {
        self.pos.mirror_pnt(p);
    }

    /// Returns the plane mirrored with respect to a point.
    fn mirrored_pnt(&self, p: &NPoint3d) -> Self {
        let mut pl = self.clone();
        pl.mirror_pnt(p);
        pl
    }

    /// Mirrors the plane with respect to an axis.
    fn mirror_ax1(&mut self, a1: &NAx1) {
        self.pos.mirror_ax1(a1);
    }

    /// Returns the plane mirrored with respect to an axis.
    fn mirrored_ax1(&self, a1: &NAx1) -> Self {
        let mut pl = self.clone();
        pl.mirror_ax1(a1);
        pl
    }

    /// Mirrors the plane with respect to a plane defined by Ax2.
    fn mirror_ax2(&mut self, a2: &NAx2) {
        self.pos.mirror_ax2(a2);
    }

    /// Returns the plane mirrored with respect to a plane defined by Ax2.
    fn mirrored_ax2(&self, a2: &NAx2) -> Self {
        let mut pl = self.clone();
        pl.mirror_ax2(a2);
        pl
    }

    /// Rotates the plane around an axis by an angle.
    fn rotate(&mut self, a1: &NAx1, ang: f64) {
        self.pos.rotate(a1, ang);
    }

    /// Returns the plane rotated around an axis by an angle.
    fn rotated(&self, a1: &NAx1, ang: f64) -> Self {
        let mut pl = self.clone();
        pl.rotate(a1, ang);
        pl
    }

    /// Scales the plane with respect to a point.
    fn scale(&mut self, p: &NPoint3d, s: f64) {
        self.pos.scale(p, s);
    }

    /// Returns the plane scaled with respect to a point.
    fn scaled(&self, p: &NPoint3d, s: f64) -> Self {
        let mut pl = self.clone();
        pl.scale(p, s);
        pl
    }

    /// Transforms the plane with a transformation.
    fn transform(&mut self, t: &NTrsf) {
        self.pos.transform(t);
    }

    /// Returns the plane transformed with a transformation.
    fn transformed(&self, t: &NTrsf) -> Self {
        let mut pl = self.clone();
        pl.transform(t);
        pl
    }

    /// Translates the plane by a vector.
    fn translate_vec(&mut self, v: &NVec) {
        self.pos.translate_vec(v);
    }

    /// Returns the plane translated by a vector.
    fn translated_vec(&self, v: &NVec) -> Self {
        let mut pl = self.clone();
        pl.translate_vec(v);
        pl
    }

    /// Translates the plane from one point to another.
    fn translate_point3d(&mut self, p1: &NPoint3d, p2: &NPoint3d) {
        self.pos.translate_point3d(p1, p2);
    }

    /// Returns the plane translated from one point to another.
    fn translated_point3d(&self, p1: &NPoint3d, p2: &NPoint3d) -> Self {
        let mut pl = self.clone();
        pl.translate_point3d(p1, p2);
        pl
    }

    /// Dumps the plane as JSON.
    fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(out, "{} {{", indent).unwrap();
        writeln!(out, "{}   \"type\": \"NPln\",", indent).unwrap();
        writeln!(out, "{}   \"position\":", indent).unwrap();
        self.pos.dump_json(out, depth + 2);
        writeln!(out, "{} }}", indent).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_plane() -> NPln {
        let a3 = NAx3::new(
            &NPoint3d::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0),
            &NDir::new(1.0, 0.0, 0.0),
        )
        .unwrap();
        NPln::new_with_ax3(&a3)
    }

    #[test]
    fn test_new() {
        let pln = NPln::new();
        assert_eq!(pln.location(), NPoint3d::new(0.0, 0.0, 0.0));
        assert_eq!(pln.pos.direction(), NDir::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn test_new_with_point_dir() {
        let p = NPoint3d::new(1.0, 2.0, 3.0);
        let v = NDir::new(0.0, 0.0, 1.0);
        let pln = NPln::new_with_point_dir(&p, &v).unwrap();
        assert_eq!(pln.location(), p);
        assert_eq!(pln.pos.direction(), v);
    }

    #[test]
    fn test_new_with_coefficients() {
        let pln = NPln::new_with_coefficients(0.0, 0.0, 1.0, -5.0).unwrap();
        let (a, b, c, d) = pln.coefficients();
        assert!((a - 0.0).abs() < 1e-9);
        assert!((b - 0.0).abs() < 1e-9);
        assert!((c - 1.0).abs() < 1e-9);
        assert!((d - -5.0).abs() < 1e-9);
        assert!(matches!(
            NPln::new_with_coefficients(0.0, 0.0, 0.0, 1.0),
            Err(NErrors::ConstructionError)
        ));
    }

    #[test]
    fn test_coefficients() {
        let pln = create_test_plane();
        let (a, b, c, d) = pln.coefficients();
        assert!((a - 0.0).abs() < 1e-9);
        assert!((b - 0.0).abs() < 1e-9);
        assert!((c - 1.0).abs() < 1e-9);
        assert!((d - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_distance_pnt() {
        let pln = create_test_plane();
        let p = NPoint3d::new(0.0, 0.0, 5.0);
        assert!((pln.distance_pnt(&p) - 5.0).abs() < 1e-9);
        assert!((pln.square_distance_pnt(&p) - 25.0).abs() < 1e-9);
    }

    #[test]
    fn test_distance_lin() {
        let pln = create_test_plane();
        let l = NLin::new_with_ax1(&NAx1::new(
            &NPoint3d::new(0.0, 0.0, 5.0),
            &NDir::new(1.0, 0.0, 0.0),
        ));
        assert!((pln.distance_lin(&l) - 5.0).abs() < 1e-9);
        let l = NLin::new_with_ax1(&NAx1::new(
            &NPoint3d::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0),
        ));
        assert!((pln.distance_lin(&l) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_contains() {
        let pln = create_test_plane();
        let p = NPoint3d::new(1.0, 1.0, 0.0);
        assert!(pln.contains_pnt(&p, 1e-6));
        let l = NLin::new_with_ax1(&NAx1::new(
            &NPoint3d::new(0.0, 0.0, 0.0),
            &NDir::new(1.0, 0.0, 0.0),
        ));
        assert!(pln.contains_lin(&l, 1e-6, 1e-6));
    }

    #[test]
    fn test_direct() {
        let pln = create_test_plane();
        assert!(pln.direct());
        let mut pln = pln;
        pln.u_reverse();
        assert!(!pln.direct());
    }

    #[test]
    fn test_dump_json() {
        let pln = create_test_plane();
        let mut output = Vec::new();
        pln.dump_json(&mut output, 0);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NPln\""));
    }
}
