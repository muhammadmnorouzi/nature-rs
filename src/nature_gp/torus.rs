use std::f64::consts::PI;
use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::gp::{NAx1, NAx2, NAx3, NPoint3d, NTrsf, NVec, NErrors};

// Assuming gp::resolution and constants are defined elsewhere
mod gp {
    pub fn resolution() -> f64 {
        1e-12 // Placeholder; replace with actual value if defined
    }
    pub fn real_last() -> f64 {
        f64::MAX
    }
    pub fn real_small() -> f64 {
        f64::MIN_POSITIVE
    }
}

// Trait to define the behavior of a torus
pub trait Torus {
    fn new() -> Self;
    fn new_with_params(a3: NAx3, major_radius: f64, minor_radius: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn set_axis(&mut self, a1: NAx1) -> Result<(), NErrors>;
    fn set_location(&mut self, loc: NPoint3d);
    fn set_major_radius(&mut self, major_radius: f64) -> Result<(), NErrors>;
    fn set_minor_radius(&mut self, minor_radius: f64) -> Result<(), NErrors>;
    fn set_position(&mut self, a3: NAx3);
    fn area(&self) -> f64;
    fn u_reverse(&mut self);
    fn v_reverse(&mut self);
    fn direct(&self) -> bool;
    fn axis(&self) -> NAx1;
    fn coefficients(&self, coef: &mut [f64]) -> Result<(), NErrors>;
    fn location(&self) -> NPoint3d;
    fn position(&self) -> NAx3;
    fn major_radius(&self) -> f64;
    fn minor_radius(&self) -> f64;
    fn volume(&self) -> f64;
    fn x_axis(&self) -> NAx1;
    fn y_axis(&self) -> NAx1;
    fn mirror_pnt(&mut self, p: NPoint3d);
    fn mirrored_pnt(&self, p: NPoint3d) -> Self;
    fn mirror_ax1(&mut self, a1: NAx1);
    fn mirrored_ax1(&self, a1: NAx1) -> Self;
    fn mirror_ax2(&mut self, a2: NAx2);
    fn mirrored_ax2(&self, a2: NAx2) -> Self;
    fn rotate(&mut self, a1: NAx1, ang: f64);
    fn rotated(&self, a1: NAx1, ang: f64) -> Self;
    fn scale(&mut self, p: NPoint3d, s: f64);
    fn scaled(&self, p: NPoint3d, s: f64) -> Self;
    fn transform(&mut self, t: NTrsf);
    fn transformed(&self, t: NTrsf) -> Self;
    fn translate_vec(&mut self, v: NVec);
    fn translated_vec(&self, v: NVec) -> Self;
    fn translate_pnts(&mut self, p1: NPoint3d, p2: NPoint3d);
    fn translated_pnts(&self, p1: NPoint3d, p2: NPoint3d) -> Self;
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
}

// Struct representing a torus in 3D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NTorus {
    pos: NAx3,
    major_radius: f64,
    minor_radius: f64,
}

impl Torus for NTorus {
    /// Creates an indefinite torus.
    fn new() -> Self {
        NTorus {
            pos: NAx3::new(NPoint3d::new(0.0, 0.0, 0.0), NVec::new(0.0, 0.0, 1.0), NVec::new(1.0, 0.0, 0.0)),
            major_radius: gp::real_last(),
            minor_radius: gp::real_small(),
        }
    }

    /// Creates a torus with the given position and radii.
    fn new_with_params(a3: NAx3, major_radius: f64, minor_radius: f64) -> Result<Self, NErrors> {
        if minor_radius < 0.0 || major_radius < 0.0 {
            return Err(NErrors::InvalidConstructionParameters);
        }
        Ok(NTorus {
            pos: a3,
            major_radius,
            minor_radius,
        })
    }

    /// Sets the symmetry axis of the torus.
    fn set_axis(&mut self, a1: NAx1) -> Result<(), NErrors> {
        // Note: C++ raises if a1's direction is parallel to XDirection.
        // Assuming NAx3::set_axis handles this validation.
        self.pos.set_axis(a1);
        Ok(())
    }

    /// Sets the location of the torus.
    fn set_location(&mut self, loc: NPoint3d) {
        self.pos.set_location(loc);
    }

    /// Sets the major radius.
    fn set_major_radius(&mut self, major_radius: f64) -> Result<(), NErrors> {
        if major_radius - self.minor_radius <= gp::resolution() {
            return Err(NErrors::InvalidConstructionParameters);
        }
        self.major_radius = major_radius;
        Ok(())
    }

    /// Sets the minor radius.
    fn set_minor_radius(&mut self, minor_radius: f64) -> Result<(), NErrors> {
        if minor_radius < 0.0 || self.major_radius - minor_radius <= gp::resolution() {
            return Err(NErrors::InvalidConstructionParameters);
        }
        self.minor_radius = minor_radius;
        Ok(())
    }

    /// Sets the local coordinate system.
    fn set_position(&mut self, a3: NAx3) {
        self.pos = a3;
    }

    /// Computes the surface area of the torus.
    fn area(&self) -> f64 {
        4.0 * PI * PI * self.minor_radius * self.major_radius
    }

    /// Reverses the U parametrization (Y-axis).
    fn u_reverse(&mut self) {
        self.pos.y_reverse();
    }

    /// Reverses the V parametrization (Z-axis).
    fn v_reverse(&mut self) {
        self.pos.z_reverse();
    }

    /// Returns true if the coordinate system is right-handed.
    fn direct(&self) -> bool {
        self.pos.direct()
    }

    /// Returns the symmetry axis.
    fn axis(&self) -> NAx1 {
        self.pos.axis()
    }

    /// Computes the coefficients of the implicit equation.
    fn coefficients(&self, coef: &mut [f64]) -> Result<(), NErrors> {
        if coef.len() < 35 {
            return Err(NErrors::DimensionMismatch);
        }

        let mut tr = NTrsf::new();
        tr.set_transformation(self.pos.clone());
        let t11 = tr.value(1, 1);
        let t12 = tr.value(1, 2);
        let t13 = tr.value(1, 3);
        let t14 = tr.value(1, 4);
        let t21 = tr.value(2, 1);
        let t22 = tr.value(2, 2);
        let t23 = tr.value(2, 3);
        let t24 = tr.value(2, 4);
        let t31 = tr.value(3, 1);
        let t32 = tr.value(3, 2);
        let t33 = tr.value(3, 3);
        let t34 = tr.value(3, 4);

        let tcol1sq = t11 * t11 + t21 * t21 + t31 * t31;
        let tcol2sq = t12 * t12 + t22 * t22 + t32 * t32;
        let tcol3sq = t13 * t13 + t23 * t23 + t33 * t33;
        let tcol4sq = t14 * t14 + t24 * t24 + t34 * t34;
        let tcol1_tcol2 = t11 * t12 + t21 * t22 + t31 * t32;
        let tcol1_tcol3 = t11 * t13 + t21 * t23 + t31 * t33;
        let tcol2_tcol3 = t12 * t13 + t22 * t23 + t32 * t33;
        let tcol1_tcol4 = t11 * t14 + t21 * t24 + t31 * t34;
        let tcol2_tcol4 = t12 * t14 + t22 * t24 + t32 * t34;
        let tcol3_tcol4 = t13 * t14 + t23 * t24 + t33 * t34;

        let sum_radius = self.major_radius * self.major_radius + self.minor_radius * self.minor_radius;
        let sub_radius = self.major_radius * self.major_radius - self.minor_radius * self.minor_radius;

        coef[0] = tcol1sq * tcol1sq; // X^4
        coef[1] = tcol2sq * tcol2sq; // Y^4
        coef[2] = tcol3sq * tcol3sq; // Z^4
        coef[3] = 4.0 * tcol1sq * tcol1_tcol2; // X^3*Y
        coef[4] = 4.0 * tcol1sq * tcol1_tcol3; // X^3*Z
        coef[5] = 4.0 * tcol2sq * tcol1_tcol2; // X*Y^3
        coef[6] = 4.0 * tcol2sq * tcol2_tcol3; // Y^3*Z
        coef[7] = 4.0 * tcol3sq * tcol1_tcol3; // X*Z^3
        coef[8] = 4.0 * tcol3sq * tcol2_tcol3; // Y*Z^3
        coef[9] = 2.0 * (tcol1sq * tcol2sq + 2.0 * tcol1_tcol2 * tcol1_tcol2); // X^2*Y^2
        coef[10] = 2.0 * (tcol1sq * tcol3sq + 2.0 * tcol1_tcol3 * tcol1_tcol3); // X^2*Z^2
        coef[11] = 2.0 * (tcol2sq * tcol3sq + 2.0 * tcol2_tcol3 * tcol2_tcol3); // Y^2*Z^2
        coef[12] = 4.0 * (tcol1sq * tcol2_tcol3 + 2.0 * tcol1_tcol2 * tcol1_tcol3); // X^2*Y*Z
        coef[13] = 4.0 * (tcol2sq * tcol1_tcol3 + 2.0 * tcol1_tcol2 * tcol2_tcol3); // X*Y^2*Z
        coef[14] = 4.0 * (tcol3sq * tcol1_tcol2 + 2.0 * tcol1_tcol3 * tcol2_tcol3); // X*Y*Z^2
        coef[15] = 4.0 * tcol1sq * tcol1_tcol4; // X^3
        coef[16] = 4.0 * tcol2sq * tcol2_tcol4; // Y^3
        coef[17] = 4.0 * tcol3sq * tcol3_tcol4; // Z^3
        coef[18] = 4.0 * (tcol1sq * tcol2_tcol4 + 2.0 * tcol1_tcol4 * tcol1_tcol2); // X^2*Y
        coef[19] = 4.0 * (tcol1sq * tcol3_tcol4 + 2.0 * tcol1_tcol4 * tcol1_tcol3); // X^2*Z
        coef[20] = 4.0 * (tcol2sq * tcol1_tcol4 + 2.0 * tcol2_tcol4 * tcol1_tcol2); // X*Y^2
        coef[21] = 4.0 * (tcol2sq * tcol3_tcol4 + 2.0 * tcol2_tcol4 * tcol2_tcol3); // Y^2*Z
        coef[22] = 4.0 * (tcol3sq * tcol1_tcol4 + 2.0 * tcol3_tcol4 * tcol1_tcol3); // X*Z^2
        coef[23] = 4.0 * (tcol3sq * tcol2_tcol4 + 2.0 * tcol3_tcol4 * tcol2_tcol3); // Y*Z^2
        coef[24] = 8.0 * (tcol1_tcol2 * tcol3_tcol4 + tcol2_tcol3 * tcol1_tcol4 + tcol2_tcol4 * tcol1_tcol3); // X*Y*Z
        coef[25] = 2.0 * (sub_radius * t31 * t31 - sum_radius * (t11 * t11 + t21 * t21)
                         + tcol4sq * tcol1sq + 2.0 * tcol1_tcol4 * tcol1_tcol4); // X^2
        coef[26] = 2.0 * (sub_radius * t32 * t32 - sum_radius * (t12 * t12 + t22 * t22)
                         + tcol4sq * tcol2sq + 2.0 * tcol2_tcol4 * tcol2_tcol4); // Y^2
        coef[27] = 2.0 * (sub_radius * t33 * t33 - sum_radius * (t13 * t13 + t23 * t23)
                         + tcol4sq * tcol3sq + 2.0 * tcol3_tcol4 * tcol3_tcol4); // Z^2
        coef[28] = 4.0 * (sub_radius * t31 * t32 - sum_radius * (t11 * t12 + t21 * t22)
                         + tcol4sq * tcol1_tcol2 + 2.0 * tcol1_tcol4 * tcol2_tcol4); // X*Y
        coef[29] = 4.0 * (sub_radius * t31 * t33 - sum_radius * (t11 * t13 + t21 * t23)
                         + tcol4sq * tcol1_tcol3 + 2.0 * tcol1_tcol4 * tcol3_tcol4); // X*Z
        coef[30] = 4.0 * (sub_radius * t32 * t33 - sum_radius * (t12 * t13 + t22 * t23)
                         + tcol4sq * tcol2_tcol3 + 2.0 * tcol2_tcol4 * tcol3_tcol4); // Y*Z
        coef[31] = 4.0 * (tcol4sq * tcol1_tcol4 + sub_radius * t31 * t34
                         - sum_radius * (t11 * t14 + t21 * t24)); // X
        coef[32] = 4.0 * (tcol4sq * tcol2_tcol4 + sub_radius * t32 * t34
                         - sum_radius * (t12 * t14 + t22 * t24)); // Y
        coef[33] = 4.0 * (tcol4sq * tcol3_tcol4 + sub_radius * t33 * t34
                         - sum_radius * (t13 * t14 + t23 * t24)); // Z
        coef[34] = 2.0 * sub_radius * t34 * t34
                  - 2.0 * sum_radius * (t14 * t14 + t24 * t24)
                  + tcol4sq * tcol4sq + sub_radius * sub_radius;

        Ok(())
    }

    /// Returns the torus's location.
    fn location(&self) -> NPoint3d {
        self.pos.location()
    }

    /// Returns the local coordinate system.
    fn position(&self) -> NAx3 {
        self.pos.clone()
    }

    /// Returns the major radius.
    fn major_radius(&self) -> f64 {
        self.major_radius
    }

    /// Returns the minor radius.
    fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// Computes the volume of the torus.
    fn volume(&self) -> f64 {
        (PI * self.minor_radius * self.minor_radius) * (2.0 * PI * self.major_radius)
    }

    /// Returns the X-axis of the torus.
    fn x_axis(&self) -> NAx1 {
        NAx1::new(self.pos.location(), self.pos.x_direction())
    }

    /// Returns the Y-axis of the torus.
    fn y_axis(&self) -> NAx1 {
        NAx1::new(self.pos.location(), self.pos.y_direction())
    }

    /// Mirrors the torus with respect to a point.
    fn mirror_pnt(&mut self, p: NPoint3d) {
        self.pos.mirror(p);
    }

    /// Returns a mirrored torus with respect to a point.
    fn mirrored_pnt(&self, p: NPoint3d) -> Self {
        let mut c = self.clone();
        c.pos.mirror(p);
        c
    }

    /// Mirrors the torus with respect to an axis.
    fn mirror_ax1(&mut self, a1: NAx1) {
        self.pos.mirror(a1);
    }

    /// Returns a mirrored torus with respect to an axis.
    fn mirrored_ax1(&self, a1: NAx1) -> Self {
        let mut c = self.clone();
        c.pos.mirror(a1);
        c
    }

    /// Mirrors the torus with respect to a plane.
    fn mirror_ax2(&mut self, a2: NAx2) {
        self.pos.mirror(a2);
    }

    /// Returns a mirrored torus with respect to a plane.
    fn mirrored_ax2(&self, a2: NAx2) -> Self {
        let mut c = self.clone();
        c.pos.mirror(a2);
        c
    }

    /// Rotates the torus around an axis.
    fn rotate(&mut self, a1: NAx1, ang: f64) {
        self.pos.rotate(a1, ang);
    }

    /// Returns a rotated torus.
    fn rotated(&self, a1: NAx1, ang: f64) -> Self {
        let mut c = self.clone();
        c.pos.rotate(a1, ang);
        c
    }

    /// Scales the torus.
    fn scale(&mut self, p: NPoint3d, s: f64) {
        self.pos.scale(p, s);
        let scale = s.abs();
        self.major_radius *= scale;
        self.minor_radius *= scale;
    }

    /// Returns a scaled torus.
    fn scaled(&self, p: NPoint3d, s: f64) -> Self {
        let mut c = self.clone();
        c.pos.scale(p, s);
        c.major_radius *= s.abs();
        c.minor_radius *= s.abs();
        c
    }

    /// Transforms the torus.
    fn transform(&mut self, t: NTrsf) {
        self.pos.transform(t);
        let scale = t.scale_factor().abs();
        self.major_radius *= scale;
        self.minor_radius *= scale;
    }

    /// Returns a transformed torus.
    fn transformed(&self, t: NTrsf) -> Self {
        let mut c = self.clone();
        c.pos.transform(t);
        c.major_radius *= t.scale_factor().abs();
        c.minor_radius *= t.scale_factor().abs();
        c
    }

    /// Translates the torus by a vector.
    fn translate_vec(&mut self, v: NVec) {
        self.pos.translate(v);
    }

    /// Returns a translated torus by a vector.
    fn translated_vec(&self, v: NVec) -> Self {
        let mut c = self.clone();
        c.pos.translate(v);
        c
    }

    /// Translates the torus from one point to another.
    fn translate_pnts(&mut self, p1: NPoint3d, p2: NPoint3d) {
        self.pos.translate(p1, p2);
    }

    /// Returns a translated torus from one point to another.
    fn translated_pnts(&self, p1: NPoint3d, p2: NPoint3d) -> Self {
        let mut c = self.clone();
        c.pos.translate(p1, p2);
        c
    }

    /// Dumps the torus as JSON.
    fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(out, "{} {{", indent).unwrap();
        writeln!(out, "{}   \"type\": \"NTorus\",", indent).unwrap();
        writeln!(out, "{}   \"position\":", indent).unwrap();
        self.pos.dump_json(out, depth + 2);
        writeln!(out, "{}   \"major_radius\": {},", indent, self.major_radius).unwrap();
        writeln!(out, "{}   \"minor_radius\": {}", indent, self.minor_radius).unwrap();
        writeln!(out, "{} }}", indent).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_torus() -> NTorus {
        NTorus::new_with_params(
            NAx3::new(NPoint3d::new(0.0, 0.0, 0.0), NVec::new(0.0, 0.0, 1.0), NVec::new(1.0, 0.0, 0.0)),
            10.0,
            2.0,
        ).unwrap()
    }

    #[test]
    fn test_new() {
        let t = NTorus::new();
        assert_eq!(t.major_radius(), gp::real_last());
        assert_eq!(t.minor_radius(), gp::real_small());
    }

    #[test]
    fn test_new_with_params() {
        let t = create_test_torus();
        assert_eq!(t.major_radius(), 10.0);
        assert_eq!(t.minor_radius(), 2.0);
    }

    #[test]
    fn test_set_major_radius() {
        let mut t = create_test_torus();
        t.set_major_radius(15.0).unwrap();
        assert_eq!(t.major_radius(), 15.0);
        assert!(t.set_major_radius(2.0).is_err());
    }

    #[test]
    fn test_set_minor_radius() {
        let mut t = create_test_torus();
        t.set_minor_radius(3.0).unwrap();
        assert_eq!(t.minor_radius(), 3.0);
        assert!(t.set_minor_radius(-1.0).is_err());
        assert!(t.set_minor_radius(10.0).is_err());
    }

    #[test]
    fn test_area() {
        let t = create_test_torus();
        assert!((t.area() - 4.0 * PI * PI * 2.0 * 10.0).abs() < 1e-9);
    }

    #[test]
    fn test_volume() {
        let t = create_test_torus();
        assert!((t.volume() - (PI * 2.0 * 2.0) * (2.0 * PI * 10.0)).abs() < 1e-9);
    }

    #[test]
    fn test_coefficients() {
        let t = create_test_torus();
        let mut coef = vec![0.0; 35];
        t.coefficients(&mut coef).unwrap();
        assert!((coef[0] - 1.0).abs() < 1e-9); // X^4 (identity transformation)
        assert!(t.coefficients(&mut vec![0.0; 34]).is_err());
    }

    #[test]
    fn test_mirror_pnt() {
        let mut t = create_test_torus();
        let p = NPoint3d::new(1.0, 2.0, 3.0);
        t.mirror_pnt(p);
        let mirrored = t.mirrored_pnt(p);
        assert_eq!(t.position(), mirrored.position());
    }

    #[test]
    fn test_scale() {
        let mut t = create_test_torus();
        t.scale(NPoint3d::new(0.0, 0.0, 0.0), 2.0);
        assert_eq!(t.major_radius(), 20.0);
        assert_eq!(t.minor_radius(), 4.0);
    }

    #[test]
    fn test_dump_json() {
        let t = create_test_torus();
        let mut output = Vec::new();
        t.dump_json(&mut output, 0);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NTorus\""));
        assert!(json.contains("\"major_radius\": 10"));
    }
}