use std::f64::consts::PI;
use std::io::Write;

use serde::{Deserialize, Serialize};

use nature_errors::NErrors;

use crate::gp::{NPoint3d, NVec, NAx1, NAx2, NAx3, NTrsf};

// Trait to define the behavior of a sphere
pub trait Sphere {
    fn new() -> Self;
    fn new_with_position_and_radius(pos: &NAx3, radius: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn set_location(&mut self, loc: &NPoint3d);
    fn set_position(&mut self, pos: &NAx3);
    fn set_radius(&mut self, radius: f64) -> Result<(), NErrors>;
    fn area(&self) -> f64;
    fn coefficients(
        &self,
        a1: &mut f64,
        a2: &mut f64,
        a3: &mut f64,
        b1: &mut f64,
        b2: &mut f64,
        b3: &mut f64,
        c1: &mut f64,
        c2: &mut f64,
        c3: &mut f64,
        d: &mut f64,
    );
    fn u_reverse(&mut self);
    fn v_reverse(&mut self);
    fn direct(&self) -> bool;
    fn location(&self) -> NPoint3d;
    fn position(&self) -> NAx3;
    fn radius(&self) -> f64;
    fn volume(&self) -> f64;
    fn x_axis(&self) -> NAx1;
    fn y_axis(&self) -> NAx1;
    fn mirror_point(&mut self, point: &NPoint3d);
    fn mirrored_point(&self, point: &NPoint3d) -> Self
    where
        Self: Sized;
    fn mirror_axis(&mut self, axis: &NAx1);
    fn mirrored_axis(&self, axis: &NAx1) -> Self
    where
        Self: Sized;
    fn mirror_plane(&mut self, plane: &NAx2);
    fn mirrored_plane(&self, plane: &NAx2) -> Self
    where
        Self: Sized;
    fn rotate(&mut self, axis: &NAx1, angle: f64);
    fn rotated(&self, axis: &NAx1, angle: f64) -> Self
    where
        Self: Sized;
    fn scale(&mut self, point: &NPoint3d, scale: f64);
    fn scaled(&self, point: &NPoint3d, scale: f64) -> Self
    where
        Self: Sized;
    fn transform(&mut self, trsf: &NTrsf);
    fn transformed(&self, trsf: &NTrsf) -> Self
    where
        Self: Sized;
    fn translate_vec(&mut self, vec: &NVec);
    fn translated_vec(&self, vec: &NVec) -> Self
    where
        Self: Sized;
    fn translate_points(&mut self, p1: &NPoint3d, p2: &NPoint3d);
    fn translated_points(&self, p1: &NPoint3d, p2: &NPoint3d) -> Self
    where
        Self: Sized;
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
}

// Struct representing a sphere in 3D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NSphere {
    pos: NAx3,
    radius: f64,
}

impl Sphere for NSphere {
    /// Creates an indefinite sphere (infinite radius).
    fn new() -> Self {
        NSphere {
            pos: NAx3::new(),
            radius: f64::INFINITY,
        }
    }

    /// Constructs a sphere with a position and radius.
    fn new_with_position_and_radius(pos: &NAx3, radius: f64) -> Result<Self, NErrors> {
        if radius < 0.0 {
            return Err(NErrors::InvalidGeometry(
                "Sphere radius must be non-negative".to_string(),
            ));
        }
        Ok(NSphere {
            pos: pos.clone(),
            radius,
        })
    }

    /// Changes the center of the sphere.
    fn set_location(&mut self, loc: &NPoint3d) {
        self.pos.set_location(loc);
    }

    /// Changes the local coordinate system of the sphere.
    fn set_position(&mut self, pos: &NAx3) {
        self.pos = pos.clone();
    }

    /// Sets the radius of the sphere.
    fn set_radius(&mut self, radius: f64) -> Result<(), NErrors> {
        if radius < 0.0 {
            return Err(NErrors::InvalidGeometry(
                "Sphere radius must be non-negative".to_string(),
            ));
        }
        self.radius = radius;
        Ok(())
    }

    /// Computes the surface area of the sphere.
    fn area(&self) -> f64 {
        4.0 * PI * self.radius * self.radius
    }

    /// Computes the coefficients of the implicit equation:
    /// a1*X^2 + a2*Y^2 + a3*Z^2 + 2*(b1*X*Y + b2*X*Z + b3*Y*Z) + 2*(c1*X + c2*Y + c3*Z) + d = 0
    fn coefficients(
        &self,
        a1: &mut f64,
        a2: &mut f64,
        a3: &mut f64,
        b1: &mut f64,
        b2: &mut f64,
        b3: &mut f64,
        c1: &mut f64,
        c2: &mut f64,
        c3: &mut f64,
        d: &mut f64,
    ) {
        let mut trsf = NTrsf::new();
        trsf.set_transformation(&self.pos);
        let t11 = trsf.value(1, 1);
        let t12 = trsf.value(1, 2);
        let t13 = trsf.value(1, 3);
        let t14 = trsf.value(1, 4);
        let t21 = trsf.value(2, 1);
        let t22 = trsf.value(2, 2);
        let t23 = trsf.value(2, 3);
        let t24 = trsf.value(2, 4);
        let t31 = trsf.value(3, 1);
        let t32 = trsf.value(3, 2);
        let t33 = trsf.value(3, 3);
        let t34 = trsf.value(3, 4);
        *a1 = t11 * t11 + t21 * t21 + t31 * t31;
        *a2 = t12 * t12 + t22 * t22 + t32 * t32;
        *a3 = t13 * t13 + t23 * t23 + t33 * t33;
        *b1 = t11 * t12 + t21 * t22 + t31 * t32;
        *b2 = t11 * t13 + t21 * t23 + t31 * t33;
        *b3 = t12 * t13 + t22 * t23 + t32 * t33;
        *c1 = t11 * t14 + t21 * t24 + t31 * t34;
        *c2 = t12 * t14 + t22 * t24 + t32 * t34;
        *c3 = t13 * t14 + t23 * t24 + t33 * t34;
        *d = t14 * t14 + t24 * t24 + t34 * t34 - self.radius * self.radius;
    }

    /// Reverses the U parametrization (reverses Y axis).
    fn u_reverse(&mut self) {
        self.pos.y_reverse();
    }

    /// Reverses the V parametrization (reverses Z axis).
    fn v_reverse(&mut self) {
        self.pos.z_reverse();
    }

    /// Returns true if the coordinate system is right-handed.
    fn direct(&self) -> bool {
        self.pos.direct()
    }

    /// Returns the center of the sphere.
    fn location(&self) -> NPoint3d {
        self.pos.location()
    }

    /// Returns the local coordinate system.
    fn position(&self) -> NAx3 {
        self.pos.clone()
    }

    /// Returns the radius.
    fn radius(&self) -> f64 {
        self.radius
    }

    /// Computes the volume of the sphere.
    fn volume(&self) -> f64 {
        (4.0 * PI * self.radius * self.radius * self.radius) / 3.0
    }

    /// Returns the X axis of the sphere.
    fn x_axis(&self) -> NAx1 {
        NAx1::new(&self.pos.location(), &self.pos.x_direction())
    }

    /// Returns the Y axis of the sphere.
    fn y_axis(&self) -> NAx1 {
        NAx1::new(&self.pos.location(), &self.pos.y_direction())
    }

    /// Mirrors the sphere with respect to a point.
    fn mirror_point(&mut self, point: &NPoint3d) {
        self.pos.mirror_point(point);
    }

    /// Returns a mirrored sphere with respect to a point.
    fn mirrored_point(&self, point: &NPoint3d) -> Self {
        let mut sphere = self.clone();
        sphere.pos.mirror_point(point);
        sphere
    }

    /// Mirrors the sphere with respect to an axis.
    fn mirror_axis(&mut self, axis: &NAx1) {
        self.pos.mirror_axis(axis);
    }

    /// Returns a mirrored sphere with respect to an axis.
    fn mirrored_axis(&self, axis: &NAx1) -> Self {
        let mut sphere = self.clone();
        sphere.pos.mirror_axis(axis);
        sphere
    }

    /// Mirrors the sphere with respect to a plane.
    fn mirror_plane(&mut self, plane: &NAx2) {
        self.pos.mirror_plane(plane);
    }

    /// Returns a mirrored sphere with respect to a plane.
    fn mirrored_plane(&self, plane: &NAx2) -> Self {
        let mut sphere = self.clone();
        sphere.pos.mirror_plane(plane);
        sphere
    }

    /// Rotates the sphere around an axis by an angle.
    fn rotate(&mut self, axis: &NAx1, angle: f64) {
        self.pos.rotate(axis, angle);
    }

    /// Returns a rotated sphere.
    fn rotated(&self, axis: &NAx1, angle: f64) -> Self {
        let mut sphere = self.clone();
        sphere.pos.rotate(axis, angle);
        sphere
    }

    /// Scales the sphere with respect to a point.
    fn scale(&mut self, point: &NPoint3d, scale: f64) {
        self.pos.scale(point, scale);
        self.radius *= scale.abs();
    }

    /// Returns a scaled sphere.
    fn scaled(&self, point: &NPoint3d, scale: f64) -> Self {
        let mut sphere = self.clone();
        sphere.pos.scale(point, scale);
        sphere.radius *= scale.abs();
        sphere
    }

    /// Transforms the sphere with a transformation.
    fn transform(&mut self, trsf: &NTrsf) {
        self.pos.transform(trsf);
        self.radius *= trsf.scale_factor().abs();
    }

    /// Returns a transformed sphere.
    fn transformed(&self, trsf: &NTrsf) -> Self {
        let mut sphere = self.clone();
        sphere.pos.transform(trsf);
        sphere.radius *= trsf.scale_factor().abs();
        sphere
    }

    /// Translates the sphere by a vector.
    fn translate_vec(&mut self, vec: &NVec) {
        self.pos.translate_vec(vec);
    }

    /// Returns a translated sphere by a vector.
    fn translated_vec(&self, vec: &NVec) -> Self {
        let mut sphere = self.clone();
        sphere.pos.translate_vec(vec);
        sphere
    }

    /// Translates the sphere from one point to another.
    fn translate_points(&mut self, p1: &NPoint3d, p2: &NPoint3d) {
        self.pos.translate_points(p1, p2);
    }

    /// Returns a translated sphere from one point to another.
    fn translated_points(&self, p1: &NPoint3d, p2: &NPoint3d) -> Self {
        let mut sphere = self.clone();
        sphere.pos.translate_points(p1, p2);
        sphere
    }

    /// Dumps the sphere as JSON.
    fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(out, "{} {{", indent).unwrap();
        writeln!(out, "{}   \"type\": \"NSphere\",", indent).unwrap();
        writeln!(out, "{}   \"position\":", indent).unwrap();
        self.pos.dump_json(out, depth + 1);
        writeln!(out, "{}   \"radius\": {}", indent, self.radius).unwrap();
        writeln!(out, "{} }}", indent).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_sphere() -> NSphere {
        let pos = NAx3::new();
        NSphere::new_with_position_and_radius(&pos, 1.0).unwrap()
    }

    #[test]
    fn test_new() {
        let sphere = NSphere::new();
        assert_eq!(sphere.radius(), f64::INFINITY);
    }

    #[test]
    fn test_new_with_position_and_radius() {
        let pos = NAx3::new();
        let sphere = NSphere::new_with_position_and_radius(&pos, 2.0).unwrap();
        assert_eq!(sphere.radius(), 2.0);
        assert!(sphere.position().is_equal(&pos));
    }

    #[test]
    fn test_new_negative_radius() {
        let pos = NAx3::new();
        let result = NSphere::new_with_position_and_radius(&pos, -1.0);
        assert!(matches!(result, Err(NErrors::InvalidGeometry(_))));
    }

    #[test]
    fn test_area() {
        let sphere = create_test_sphere();
        assert!((sphere.area() - 4.0 * PI).abs() < 1e-9);
    }

    #[test]
    fn test_volume() {
        let sphere = create_test_sphere();
        assert!((sphere.volume() - (4.0 * PI / 3.0)).abs() < 1e-9);
    }

    #[test]
    fn test_coefficients() {
        let sphere = create_test_sphere();
        let mut a1 = 0.0;
        let mut a2 = 0.0;
        let mut a3 = 0.0;
        let mut b1 = 0.0;
        let mut b2 = 0.0;
        let mut b3 = 0.0;
        let mut c1 = 0.0;
        let mut c2 = 0.0;
        let mut c3 = 0.0;
        let mut d = 0.0;
        sphere.coefficients(&mut a1, &mut a2, &mut a3, &mut b1, &mut b2, &mut b3, &mut c1, &mut c2, &mut c3, &mut d);
        assert!((a1 - 1.0).abs() < 1e-9);
        assert!((a2 - 1.0).abs() < 1e-9);
        assert!((a3 - 1.0).abs() < 1e-9);
        assert!(b1.abs() < 1e-9);
        assert!(d - (-1.0).abs() < 1e-9);
    }

    #[test]
    fn test_scale() {
        let mut sphere = create_test_sphere();
        let point = NPoint3d::new_with_coords(0.0, 0.0, 0.0);
        sphere.scale(&point, 2.0);
        assert_eq!(sphere.radius(), 2.0);
        sphere.scale(&point, -2.0);
        assert_eq!(sphere.radius(), 4.0);
    }

    #[test]
    fn test_transform() {
        let mut sphere = create_test_sphere();
        let mut trsf = NTrsf::new();
        trsf.set_scale_factor(2.0);
        sphere.transform(&trsf);
        assert_eq!(sphere.radius(), 2.0);
    }

    #[test]
    fn test_dump_json() {
        let sphere = create_test_sphere();
        let mut output = Vec::new();
        sphere.dump_json(&mut output, 0);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NSphere\""));
        assert!(json.contains("\"radius\": 1"));
    }
}