use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::{
    nature_errors::NErrors,
    nature_gp::{NAx2d, NAx22d, NDir2d, NGP, NPoint2d, NTrsf2d, NVec2d, NXY},
};

// Trait to define the behavior of a hyperbola branch in 2D space
pub trait Hypr2d {
    fn new() -> Self;
    fn new_with_major_axis(
        major_axis: &NAx2d,
        major_radius: f64,
        minor_radius: f64,
        is_sense: bool,
    ) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn new_with_axis(a: &NAx22d, major_radius: f64, minor_radius: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn set_location(&mut self, p: &NPoint2d);
    fn set_major_radius(&mut self, major_radius: f64) -> Result<(), NErrors>;
    fn set_minor_radius(&mut self, minor_radius: f64) -> Result<(), NErrors>;
    fn set_axis(&mut self, a: &NAx22d);
    fn set_x_axis(&mut self, a: &NAx2d);
    fn set_y_axis(&mut self, a: &NAx2d);
    fn asymptote1(&self) -> Result<NAx2d, NErrors>;
    fn asymptote2(&self) -> Result<NAx2d, NErrors>;
    fn coefficients(&self) -> (f64, f64, f64, f64, f64, f64);
    fn conjugate_branch1(&self) -> Self
    where
        Self: Sized;
    fn conjugate_branch2(&self) -> Self
    where
        Self: Sized;
    fn directrix1(&self) -> NAx2d;
    fn directrix2(&self) -> NAx2d;
    fn eccentricity(&self) -> Result<f64, NErrors>;
    fn focal(&self) -> f64;
    fn focus1(&self) -> NPoint2d;
    fn focus2(&self) -> NPoint2d;
    fn location(&self) -> &NPoint2d;
    fn major_radius(&self) -> f64;
    fn minor_radius(&self) -> f64;
    fn other_branch(&self) -> Self
    where
        Self: Sized;
    fn parameter(&self) -> Result<f64, NErrors>;
    fn axis(&self) -> &NAx22d;
    fn x_axis(&self) -> NAx2d;
    fn y_axis(&self) -> NAx2d;
    fn reverse(&mut self);
    fn reversed(&self) -> Self
    where
        Self: Sized;
    fn is_direct(&self) -> bool;
    fn mirror_point3d(&mut self, p: &NPoint2d);
    fn mirrored_point3d(&self, p: &NPoint2d) -> Self
    where
        Self: Sized;
    fn mirror_ax2d(&mut self, a: &NAx2d);
    fn mirrored_ax2d(&self, a: &NAx2d) -> Self
    where
        Self: Sized;
    fn rotate(&mut self, p: &NPoint2d, ang: f64);
    fn rotated(&self, p: &NPoint2d, ang: f64) -> Self
    where
        Self: Sized;
    fn scale(&mut self, p: &NPoint2d, s: f64);
    fn scaled(&self, p: &NPoint2d, s: f64) -> Self
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
    fn translate_point3d(&mut self, p1: &NPoint2d, p2: &NPoint2d);
    fn translated_point3d(&self, p1: &NPoint2d, p2: &NPoint2d) -> Self
    where
        Self: Sized;
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
}

// Struct representing a branch of a hyperbola in 2D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NHypr2d {
    pos: NAx22d,
    major_radius: f64,
    minor_radius: f64,
}

impl Hypr2d for NHypr2d {
    /// Creates an indefinite hyperbola.
    fn new() -> Self {
        NHypr2d {
            pos: NAx22d::new(
                &NPoint2d::new(0.0, 0.0),
                &NDir2d::new(1.0, 0.0).unwrap(),
                &NDir2d::new(0.0, 1.0).unwrap(),
            )
            .unwrap(),
            major_radius: f64::MAX,
            minor_radius: f64::MAX,
        }
    }

    /// Creates a hyperbola with major axis, radii, and sense.
    fn new_with_major_axis(
        major_axis: &NAx2d,
        major_radius: f64,
        minor_radius: f64,
        is_sense: bool,
    ) -> Result<Self, NErrors> {
        if major_radius < 0.0 || minor_radius < 0.0 {
            return Err(NErrors::ConstructionError);
        }
        Ok(NHypr2d {
            pos: NAx22d::new_from_axis(major_axis, is_sense),
            major_radius,
            minor_radius,
        })
    }

    /// Creates a hyperbola with coordinate system and radii.
    fn new_with_axis(a: &NAx22d, major_radius: f64, minor_radius: f64) -> Result<Self, NErrors> {
        if major_radius < 0.0 || minor_radius < 0.0 {
            return Err(NErrors::ConstructionError);
        }
        Ok(NHypr2d {
            pos: a.clone(),
            major_radius,
            minor_radius,
        })
    }

    /// Sets the location (center) of the hyperbola.
    fn set_location(&mut self, p: &NPoint2d) {
        self.pos.set_location(p);
    }

    /// Sets the major radius of the hyperbola.
    fn set_major_radius(&mut self, major_radius: f64) -> Result<(), NErrors> {
        if major_radius < 0.0 {
            return Err(NErrors::ConstructionError);
        }
        self.major_radius = major_radius;
        Ok(())
    }

    /// Sets the minor radius of the hyperbola.
    fn set_minor_radius(&mut self, minor_radius: f64) -> Result<(), NErrors> {
        if minor_radius < 0.0 {
            return Err(NErrors::ConstructionError);
        }
        self.minor_radius = minor_radius;
        Ok(())
    }

    /// Sets the coordinate system of the hyperbola.
    fn set_axis(&mut self, a: &NAx22d) {
        self.pos.set_axis(a);
    }

    /// Sets the major axis of the hyperbola.
    fn set_x_axis(&mut self, a: &NAx2d) {
        self.pos.set_x_axis(a);
    }

    /// Sets the minor axis of the hyperbola.
    fn set_y_axis(&mut self, a: &NAx2d) {
        self.pos.set_y_axis(a);
    }

    /// Returns the first asymptote (Y = (minor/major)*X).
    fn asymptote1(&self) -> Result<NAx2d, NErrors> {
        if self.major_radius <= NGP::resolution() {
            return Err(NErrors::ConstructionError);
        }
        let mut v_dir = self.pos.x_direction();
        let mut coord1 = self.pos.y_direction().xy();
        let coord2 = coord1.multiplied(self.minor_radius / self.major_radius);
        coord1.add(&coord2);
        v_dir.set_xy(&coord1);
        Ok(NAx2d::new(&self.pos.location(), &v_dir))
    }

    /// Returns the second asymptote (Y = -(minor/major)*X).
    fn asymptote2(&self) -> Result<NAx2d, NErrors> {
        if self.major_radius <= NGP::resolution() {
            return Err(NErrors::ConstructionError);
        }
        let mut v_dir = self.pos.x_direction();
        let mut coord1 = self.pos.y_direction().xy();
        let coord2 = coord1.multiplied(-self.minor_radius / self.major_radius);
        coord1.add(&coord2);
        v_dir.set_xy(&coord1);
        Ok(NAx2d::new(&self.pos.location(), &v_dir))
    }

    /// Computes the coefficients of the implicit equation A*X^2 + B*Y^2 + 2*C*X*Y + 2*D*X + 2*E*Y + F = 0.
    fn coefficients(&self) -> (f64, f64, f64, f64, f64, f64) {
        let d_min = self.minor_radius * self.minor_radius;
        let d_maj = self.major_radius * self.major_radius;
        if d_min <= NGP::resolution() && d_maj <= NGP::resolution() {
            return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
        let mut t = NTrsf2d::new();
        t.set_transformation(&self.pos.x_axis());
        let t11 = t.value(1, 1).unwrap();
        let t12 = t.value(1, 2).unwrap();
        let t13 = t.value(1, 3).unwrap();
        if d_min <= NGP::resolution() {
            let a = t11 * t11;
            let b = t12 * t12;
            let c = t11 * t12;
            let d = t11 * t13;
            let e = t12 * t13;
            let f = t13 * t13 - d_maj;
            return (a, b, c, d, e, f);
        }
        let t21 = t.value(2, 1).unwrap();
        let t22 = t.value(2, 2).unwrap();
        let t23 = t.value(2, 3).unwrap();
        let a = (t11 * t11 / d_maj) - (t21 * t21 / d_min);
        let b = (t12 * t12 / d_maj) - (t22 * t22 / d_min);
        let c = (t11 * t12 / d_maj) - (t21 * t22 / d_min);
        let d = (t11 * t13 / d_maj) - (t21 * t23 / d_min);
        let e = (t12 * t13 / d_maj) - (t22 * t23 / d_min);
        let f = (t13 * t13 / d_maj) - (t23 * t23 / d_min) - 1.0;
        (a, b, c, d, e, f)
    }

    /// Returns the conjugate branch on the positive Y-axis side.
    fn conjugate_branch1(&self) -> Self {
        let is_sign = self.is_direct();
        NHypr2d {
            pos: NAx22d::new(
                &self.pos.location(),
                &self.pos.y_direction(),
                &self.pos.x_direction(),
            )
            .unwrap(),
            major_radius: self.minor_radius,
            minor_radius: self.major_radius,
        }
    }

    /// Returns the conjugate branch on the negative Y-axis side.
    fn conjugate_branch2(&self) -> Self {
        let mut v = self.pos.y_direction();
        v.reverse();
        let is_sign = self.is_direct();
        NHypr2d {
            pos: NAx22d::new(&self.pos.location(), &v, &self.pos.x_direction()).unwrap(),
            major_radius: self.minor_radius,
            minor_radius: self.major_radius,
        }
    }

    /// Returns the first directrix (positive side of X-axis).
    fn directrix1(&self) -> NAx2d {
        let e = self.eccentricity().unwrap();
        let mut orig = self.pos.x_direction().xy();
        orig.multiply(self.major_radius / e);
        orig.add(&self.pos.location().xy());
        NAx2d::new(&NPoint2d::from_xy(&orig), &self.pos.y_direction())
    }

    /// Returns the second directrix (negative side of X-axis).
    fn directrix2(&self) -> NAx2d {
        let e = self.eccentricity().unwrap();
        let mut orig = self.pos.x_direction().xy();
        orig.multiply(self.parameter().unwrap() / e);
        orig.add(&self.focus1().xy());
        NAx2d::new(&NPoint2d::from_xy(&orig), &self.pos.y_direction())
    }

    /// Returns the eccentricity of the hyperbola (e > 1).
    fn eccentricity(&self) -> Result<f64, NErrors> {
        if self.major_radius <= NGP::resolution() {
            return Err(NErrors::DomainError);
        }
        Ok(
            (self.major_radius * self.major_radius + self.minor_radius * self.minor_radius).sqrt()
                / self.major_radius,
        )
    }

    /// Returns the focal distance (distance between foci).
    fn focal(&self) -> f64 {
        2.0 * (self.major_radius * self.major_radius + self.minor_radius * self.minor_radius).sqrt()
    }

    /// Returns the first focus (positive X-axis side).
    fn focus1(&self) -> NPoint2d {
        let c =
            (self.major_radius * self.major_radius + self.minor_radius * self.minor_radius).sqrt();
        NPoint2d::new(
            self.pos.location().x() + c * self.pos.x_direction().x(),
            self.pos.location().y() + c * self.pos.x_direction().y(),
        )
    }

    /// Returns the second focus (negative X-axis side).
    fn focus2(&self) -> NPoint2d {
        let c =
            (self.major_radius * self.major_radius + self.minor_radius * self.minor_radius).sqrt();
        NPoint2d::new(
            self.pos.location().x() - c * self.pos.x_direction().x(),
            self.pos.location().y() - c * self.pos.x_direction().y(),
        )
    }

    /// Returns the location (center) of the hyperbola.
    fn location(&self) -> &NPoint2d {
        self.pos.location()
    }

    /// Returns the major radius of the hyperbola.
    fn major_radius(&self) -> f64 {
        self.major_radius
    }

    /// Returns the minor radius of the hyperbola.
    fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// Returns the other branch (symmetric about Y-axis).
    fn other_branch(&self) -> Self {
        let is_sign = self.is_direct();
        NHypr2d {
            pos: NAx22d::new(
                &self.pos.location(),
                &self.pos.x_direction().reversed(),
                &self.pos.y_direction(),
            )
            .unwrap(),
            major_radius: self.major_radius,
            minor_radius: self.minor_radius,
        }
    }

    /// Returns the parameter p = (e^2 - 1) * major_radius.
    fn parameter(&self) -> Result<f64, NErrors> {
        if self.major_radius <= NGP::resolution() {
            return Err(NErrors::DomainError);
        }
        Ok((self.minor_radius * self.minor_radius) / self.major_radius)
    }

    /// Returns the coordinate system of the hyperbola.
    fn axis(&self) -> &NAx22d {
        &self.pos
    }

    /// Returns the major axis (X-axis) of the hyperbola.
    fn x_axis(&self) -> NAx2d {
        self.pos.x_axis()
    }

    /// Returns the minor axis (Y-axis) of the hyperbola.
    fn y_axis(&self) -> NAx2d {
        self.pos.y_axis()
    }

    /// Reverses the orientation of the hyperbola.
    fn reverse(&mut self) {
        let mut temp = self.pos.y_direction();
        temp.reverse();
        self.pos
            .set_axis(&NAx22d::new(&self.pos.location(), &self.pos.x_direction(), &temp).unwrap());
    }

    /// Returns a reversed hyperbola.
    fn reversed(&self) -> Self {
        let mut h = self.clone();
        h.reverse();
        h
    }

    /// Returns true if the coordinate system is direct.
    fn is_direct(&self) -> bool {
        self.pos.x_direction().crossed(&self.pos.y_direction()) >= 0.0
    }

    /// Mirrors the hyperbola about a point.
    fn mirror_point3d(&mut self, p: &NPoint2d) {
        self.pos.mirror_point3d(p);
    }

    /// Returns a hyperbola mirrored about a point.
    fn mirrored_point3d(&self, p: &NPoint2d) -> Self {
        let mut h = self.clone();
        h.mirror_point3d(p);
        h
    }

    /// Mirrors the hyperbola about an axis.
    fn mirror_ax2d(&mut self, a: &NAx2d) {
        self.pos.mirror_ax2d(a);
    }

    /// Returns a hyperbola mirrored about an axis.
    fn mirrored_ax2d(&self, a: &NAx2d) -> Self {
        let mut h = self.clone();
        h.mirror_ax2d(a);
        h
    }

    /// Rotates the hyperbola about a point.
    fn rotate(&mut self, p: &NPoint2d, ang: f64) {
        self.pos.rotate(p, ang);
    }

    /// Returns a rotated hyperbola.
    fn rotated(&self, p: &NPoint2d, ang: f64) -> Self {
        let mut h = self.clone();
        h.rotate(p, ang);
        h
    }

    /// Scales the hyperbola about a point.
    fn scale(&mut self, p: &NPoint2d, s: f64) {
        self.major_radius *= s;
        if self.major_radius < 0.0 {
            self.major_radius = -self.major_radius;
        }
        self.minor_radius *= s;
        if self.minor_radius < 0.0 {
            self.minor_radius = -self.minor_radius;
        }
        self.pos.scale(p, s);
    }

    /// Returns a scaled hyperbola.
    fn scaled(&self, p: &NPoint2d, s: f64) -> Self {
        let mut h = self.clone();
        h.scale(p, s);
        h
    }

    /// Transforms the hyperbola with a transformation.
    fn transform(&mut self, t: &NTrsf2d) {
        let scale = t.scale_factor();
        self.major_radius *= scale;
        if self.major_radius < 0.0 {
            self.major_radius = -self.major_radius;
        }
        self.minor_radius *= scale;
        if self.minor_radius < 0.0 {
            self.minor_radius = -self.minor_radius;
        }
        self.pos.transform(t);
    }

    /// Returns a transformed hyperbola.
    fn transformed(&self, t: &NTrsf2d) -> Self {
        let mut h = self.clone();
        h.transform(t);
        h
    }

    /// Translates the hyperbola by a vector.
    fn translate_vec(&mut self, v: &NVec2d) {
        self.pos.translate_vec(v);
    }

    /// Returns a translated hyperbola by a vector.
    fn translated_vec(&self, v: &NVec2d) -> Self {
        let mut h = self.clone();
        h.translate_vec(v);
        h
    }

    /// Translates the hyperbola from one point to another.
    fn translate_point3d(&mut self, p1: &NPoint2d, p2: &NPoint2d) {
        self.pos.translate_point3d(p1, p2);
    }

    /// Returns a translated hyperbola from one point to another.
    fn translated_point3d(&self, p1: &NPoint2d, p2: &NPoint2d) -> Self {
        let mut h = self.clone();
        h.translate_point3d(p1, p2);
        h
    }

    /// Dumps the hyperbola as JSON.
    fn dump_json(&self, out: &mut dyn Write, _depth: i32) {
        writeln!(out, "{{\n  \"type\": \"NHypr2d\",").unwrap();
        writeln!(
            out,
            "  \"pos\": {{ \"location\": [{}, {}], \"x_direction\": [{}, {}], \"y_direction\": [{}, {}] }},",
            self.pos.location().x(), self.pos.location().y(),
            self.pos.x_direction().x(), self.pos.x_direction().y(),
            self.pos.y_direction().x(), self.pos.y_direction().y()
        ).unwrap();
        writeln!(out, "  \"major_radius\": {},", self.major_radius).unwrap();
        writeln!(out, "  \"minor_radius\": {}", self.minor_radius).unwrap();
        writeln!(out, "}}").unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hypr2d() -> NHypr2d {
        NHypr2d::new_with_major_axis(
            &NAx2d::new(&NPoint2d::new(0.0, 0.0), &NDir2d::new(1.0, 0.0).unwrap()),
            2.0,
            1.0,
            true,
        )
        .unwrap()
    }

    #[test]
    fn test_new() {
        let h = NHypr2d::new();
        assert_eq!(h.major_radius(), f64::MAX);
        assert_eq!(h.minor_radius(), f64::MAX);
    }

    #[test]
    fn test_new_with_major_axis() {
        let h = hypr2d();
        assert_eq!(h.major_radius(), 2.0);
        assert_eq!(h.minor_radius(), 1.0);
        assert!(matches!(
            NHypr2d::new_with_major_axis(
                &NAx2d::new(&NPoint2d::new(0.0, 0.0), &NDir2d::new(1.0, 0.0).unwrap()),
                -1.0,
                1.0,
                true
            ),
            Err(NErrors::ConstructionError)
        ));
    }

    #[test]
    fn test_setters() {
        let mut h = hypr2d();
        h.set_major_radius(3.0).unwrap();
        assert_eq!(h.major_radius(), 3.0);
        h.set_minor_radius(2.0).unwrap();
        assert_eq!(h.minor_radius(), 2.0);
        let p = NPoint2d::new(1.0, 2.0);
        h.set_location(&p);
        assert_eq!(h.location(), &p);
    }

    #[test]
    fn test_asymptotes() {
        let h = hypr2d();
        let a1 = h.asymptote1().unwrap();
        let a2 = h.asymptote2().unwrap();
        assert!(
            a1.direction()
                .is_parallel(&NDir2d::new(2.0, 1.0).unwrap(), 1e-9)
        );
        assert!(
            a2.direction()
                .is_parallel(&NDir2d::new(2.0, -1.0).unwrap(), 1e-9)
        );
    }

    #[test]
    fn test_coefficients() {
        let h = hypr2d();
        let (a, b, c, d, e, f) = h.coefficients();
        assert!((a - 0.25).abs() < 1e-9);
        assert!((b + 1.0).abs() < 1e-9);
        assert!(c.abs() < 1e-9);
        assert!(d.abs() < 1e-9);
        assert!(e.abs() < 1e-9);
        assert!((f + 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_foci() {
        let h = hypr2d();
        let f1 = h.focus1();
        let f2 = h.focus2();
        assert!(f1.distance(&NPoint2d::new(5.0f64.sqrt(), 0.0)) < 1e-9);
        assert!(f2.distance(&NPoint2d::new(-5.0f64.sqrt(), 0.0)) < 1e-9);
    }

    #[test]
    fn test_eccentricity() {
        let h = hypr2d();
        assert!((h.eccentricity().unwrap() - 5.0f64.sqrt() / 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_directrices() {
        let h = hypr2d();
        let d1 = h.directrix1();
        let d2 = h.directrix2();
        assert!(
            d1.location()
                .distance(&NPoint2d::new(4.0 / 5.0f64.sqrt(), 0.0))
                < 1e-9
        );
        assert!(
            d2.location()
                .distance(&NPoint2d::new(5.0f64.sqrt() - 0.5, 0.0))
                < 1e-9
        );
    }

    #[test]
    fn test_transformations() {
        let h = hypr2d();
        let mut h_scaled = h.scaled(&NPoint2d::new(0.0, 0.0), 2.0);
        assert_eq!(h_scaled.major_radius(), 4.0);
        assert_eq!(h_scaled.minor_radius(), 2.0);

        let mut h_mirrored = h.mirrored_point3d(&NPoint2d::new(1.0, 0.0));
        assert!(h_mirrored.location().distance(&NPoint2d::new(2.0, 0.0)) < 1e-9);
    }

    #[test]
    fn test_reverse() {
        let mut h = hypr2d();
        let is_direct = h.is_direct();
        h.reverse();
        assert_eq!(h.is_direct(), !is_direct);
    }

    #[test]
    fn test_dump_json() {
        let h = hypr2d();
        let mut output = Vec::new();
        h.dump_json(&mut output, -1);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NHypr2d\""));
        assert!(json.contains("\"major_radius\": 2"));
    }
}
