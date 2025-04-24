use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::{
    gp::{NAx1, NAx2, NDir, NGP, NPnt, NTrsf, NVec, NXYZ},
    nature_errors::NErrors,
};

// Trait to define the behavior of a hyperbola branch in 3D space
pub trait Hypr {
    fn new() -> Self;
    fn new_with_params(a2: &NAx2, major_radius: f64, minor_radius: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn set_axis(&mut self, a1: &NAx1) -> Result<(), NErrors>;
    fn set_location(&mut self, p: &NPnt);
    fn set_major_radius(&mut self, major_radius: f64) -> Result<(), NErrors>;
    fn set_minor_radius(&mut self, minor_radius: f64) -> Result<(), NErrors>;
    fn set_position(&mut self, a2: &NAx2);
    fn asymptote1(&self) -> Result<NAx1, NErrors>;
    fn asymptote2(&self) -> Result<NAx1, NErrors>;
    fn axis(&self) -> &NAx1;
    fn conjugate_branch1(&self) -> Self
    where
        Self: Sized;
    fn conjugate_branch2(&self) -> Self
    where
        Self: Sized;
    fn directrix1(&self) -> NAx1;
    fn directrix2(&self) -> NAx1;
    fn eccentricity(&self) -> Result<f64, NErrors>;
    fn focal(&self) -> f64;
    fn focus1(&self) -> NPnt;
    fn focus2(&self) -> NPnt;
    fn location(&self) -> &NPnt;
    fn major_radius(&self) -> f64;
    fn minor_radius(&self) -> f64;
    fn other_branch(&self) -> Self
    where
        Self: Sized;
    fn parameter(&self) -> Result<f64, NErrors>;
    fn position(&self) -> &NAx2;
    fn x_axis(&self) -> NAx1;
    fn y_axis(&self) -> NAx1;
    fn mirror_pnt(&mut self, p: &NPnt);
    fn mirrored_pnt(&self, p: &NPnt) -> Self
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
    fn scale(&mut self, p: &NPnt, s: f64);
    fn scaled(&self, p: &NPnt, s: f64) -> Self
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
    fn translate_pnts(&mut self, p1: &NPnt, p2: &NPnt);
    fn translated_pnts(&self, p1: &NPnt, p2: &NPnt) -> Self
    where
        Self: Sized;
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
}

// Struct representing a branch of a hyperbola in 3D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NHypr {
    pos: NAx2,
    major_radius: f64,
    minor_radius: f64,
}

impl Hypr for NHypr {
    /// Creates an indefinite hyperbola.
    fn new() -> Self {
        NHypr {
            pos: NAx2::new(
                &NPnt::new(0.0, 0.0, 0.0),
                &NDir::new(0.0, 0.0, 1.0).unwrap(),
                &NDir::new(1.0, 0.0, 0.0).unwrap(),
            )
            .unwrap(),
            major_radius: f64::MAX,
            minor_radius: f64::MIN,
        }
    }

    /// Creates a hyperbola with given position and radii.
    fn new_with_params(a2: &NAx2, major_radius: f64, minor_radius: f64) -> Result<Self, NErrors> {
        if major_radius < 0.0 || minor_radius < 0.0 {
            return Err(NErrors::ConstructionError);
        }
        Ok(NHypr {
            pos: a2.clone(),
            major_radius,
            minor_radius,
        })
    }

    /// Sets the main axis of the hyperbola.
    fn set_axis(&mut self, a1: &NAx1) -> Result<(), NErrors> {
        self.pos.set_axis(a1)?;
        Ok(())
    }

    /// Sets the location (center) of the hyperbola.
    fn set_location(&mut self, p: &NPnt) {
        self.pos = NAx2::new(p, &self.pos.direction(), &self.pos.x_direction()).unwrap();
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

    /// Sets the position (local coordinate system) of the hyperbola.
    fn set_position(&mut self, a2: &NAx2) {
        self.pos = a2.clone();
    }

    /// Returns the first asymptote (Y = (minor/major)*X).
    fn asymptote1(&self) -> Result<NAx1, NErrors> {
        if self.major_radius <= NGP::resolution() {
            return Err(NErrors::ConstructionError);
        }
        let mut v1 = NVec::from_dir(&self.pos.y_direction());
        v1.multiply(self.minor_radius / self.major_radius);
        let mut v = NVec::from_dir(&self.pos.x_direction());
        v.add(&v1);
        Ok(NAx1::new(&self.pos.location(), &NDir::new_from_vec(&v)?))
    }

    /// Returns the second asymptote (Y = -(minor/major)*X).
    fn asymptote2(&self) -> Result<NAx1, NErrors> {
        if self.major_radius <= NGP::resolution() {
            return Err(NErrors::ConstructionError);
        }
        let mut v1 = NVec::from_dir(&self.pos.y_direction());
        v1.multiply(-self.minor_radius / self.major_radius);
        let mut v = NVec::from_dir(&self.pos.x_direction());
        v.add(&v1);
        Ok(NAx1::new(&self.pos.location(), &NDir::new_from_vec(&v)?))
    }

    /// Returns the main axis of the hyperbola (normal to its plane).
    fn axis(&self) -> &NAx1 {
        self.pos.axis()
    }

    /// Returns the conjugate branch on the positive Y-axis side.
    fn conjugate_branch1(&self) -> Self {
        NHypr {
            pos: NAx2::new(
                &self.pos.location(),
                &self.pos.direction(),
                &self.pos.y_direction(),
            )
            .unwrap(),
            major_radius: self.minor_radius,
            minor_radius: self.major_radius,
        }
    }

    /// Returns the conjugate branch on the negative Y-axis side.
    fn conjugate_branch2(&self) -> Self {
        let mut d = self.pos.y_direction();
        d.reverse();
        NHypr {
            pos: NAx2::new(&self.pos.location(), &self.pos.direction(), &d).unwrap(),
            major_radius: self.minor_radius,
            minor_radius: self.major_radius,
        }
    }

    /// Returns the first directrix (positive side of X-axis).
    fn directrix1(&self) -> NAx1 {
        let e = self.eccentricity().unwrap();
        let mut orig = self.pos.x_direction().xyz();
        orig.multiply(self.major_radius / e);
        orig.add(&self.pos.location().xyz());
        NAx1::new(&NPnt::from_xyz(&orig), &self.pos.y_direction())
    }

    /// Returns the second directrix (negative side of X-axis).
    fn directrix2(&self) -> NAx1 {
        let e = self.eccentricity().unwrap();
        let mut orig = self.pos.x_direction().xyz();
        orig.multiply(-self.major_radius / e);
        orig.add(&self.pos.location().xyz());
        NAx1::new(&NPnt::from_xyz(&orig), &self.pos.y_direction())
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
    fn focus1(&self) -> NPnt {
        let c =
            (self.major_radius * self.major_radius + self.minor_radius * self.minor_radius).sqrt();
        let pp = self.pos.location();
        let dd = self.pos.x_direction();
        NPnt::new(
            pp.x() + c * dd.x(),
            pp.y() + c * dd.y(),
            pp.z() + c * dd.z(),
        )
    }

    /// Returns the second focus (negative X-axis side).
    fn focus2(&self) -> NPnt {
        let c =
            (self.major_radius * self.major_radius + self.minor_radius * self.minor_radius).sqrt();
        let pp = self.pos.location();
        let dd = self.pos.x_direction();
        NPnt::new(
            pp.x() - c * dd.x(),
            pp.y() - c * dd.y(),
            pp.z() - c * dd.z(),
        )
    }

    /// Returns the location (center) of the hyperbola.
    fn location(&self) -> &NPnt {
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
        let mut d = self.pos.x_direction();
        d.reverse();
        NHypr {
            pos: NAx2::new(&self.pos.location(), &self.pos.direction(), &d).unwrap(),
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
    fn position(&self) -> &NAx2 {
        &self.pos
    }

    /// Returns the major axis (X-axis) of the hyperbola.
    fn x_axis(&self) -> NAx1 {
        NAx1::new(&self.pos.location(), &self.pos.x_direction())
    }

    /// Returns the minor axis (Y-axis) of the hyperbola.
    fn y_axis(&self) -> NAx1 {
        NAx1::new(&self.pos.location(), &self.pos.y_direction())
    }

    /// Mirrors the hyperbola about a point.
    fn mirror_pnt(&mut self, p: &NPnt) {
        self.pos.mirror_pnt(p);
    }

    /// Returns a hyperbola mirrored about a point.
    fn mirrored_pnt(&self, p: &NPnt) -> Self {
        let mut h = self.clone();
        h.mirror_pnt(p);
        h
    }

    /// Mirrors the hyperbola about an axis.
    fn mirror_ax1(&mut self, a1: &NAx1) {
        self.pos.mirror_ax1(a1);
    }

    /// Returns a hyperbola mirrored about an axis.
    fn mirrored_ax1(&self, a1: &NAx1) -> Self {
        let mut h = self.clone();
        h.mirror_ax1(a1);
        h
    }

    /// Mirrors the hyperbola about a plane.
    fn mirror_ax2(&mut self, a2: &NAx2) {
        self.pos.mirror_ax2(a2);
    }

    /// Returns a hyperbola mirrored about a plane.
    fn mirrored_ax2(&self, a2: &NAx2) -> Self {
        let mut h = self.clone();
        h.mirror_ax2(a2);
        h
    }

    /// Rotates the hyperbola about an axis.
    fn rotate(&mut self, a1: &NAx1, ang: f64) {
        self.pos.rotate(a1, ang);
    }

    /// Returns a rotated hyperbola.
    fn rotated(&self, a1: &NAx1, ang: f64) -> Self {
        let mut h = self.clone();
        h.rotate(a1, ang);
        h
    }

    /// Scales the hyperbola about a point.
    fn scale(&mut self, p: &NPnt, s: f64) {
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
    fn scaled(&self, p: &NPnt, s: f64) -> Self {
        let mut h = self.clone();
        h.scale(p, s);
        h
    }

    /// Transforms the hyperbola with a transformation.
    fn transform(&mut self, t: &NTrsf) {
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
    fn transformed(&self, t: &NTrsf) -> Self {
        let mut h = self.clone();
        h.transform(t);
        h
    }

    /// Translates the hyperbola by a vector.
    fn translate_vec(&mut self, v: &NVec) {
        self.pos.translate_vec(v);
    }

    /// Returns a translated hyperbola by a vector.
    fn translated_vec(&self, v: &NVec) -> Self {
        let mut h = self.clone();
        h.translate_vec(v);
        h
    }

    /// Translates the hyperbola from one point to another.
    fn translate_pnts(&mut self, p1: &NPnt, p2: &NPnt) {
        self.pos.translate_pnts(p1, p2);
    }

    /// Returns a translated hyperbola from one point to another.
    fn translated_pnts(&self, p1: &NPnt, p2: &NPnt) -> Self {
        let mut h = self.clone();
        h.translate_pnts(p1, p2);
        h
    }

    /// Dumps the hyperbola as JSON.
    fn dump_json(&self, out: &mut dyn Write, _depth: i32) {
        writeln!(out, "{{\n  \"type\": \"NHypr\",").unwrap();
        writeln!(
            out,
            "  \"pos\": {{ \"location\": [{}, {}, {}], \"x_direction\": [{}, {}, {}], \"y_direction\": [{}, {}, {}] }},",
            self.pos.location().x(), self.pos.location().y(), self.pos.location().z(),
            self.pos.x_direction().x(), self.pos.x_direction().y(), self.pos.x_direction().z(),
            self.pos.y_direction().x(), self.pos.y_direction().y(), self.pos.y_direction().z()
        ).unwrap();
        writeln!(out, "  \"major_radius\": {},", self.major_radius).unwrap();
        writeln!(out, "  \"minor_radius\": {}", self.minor_radius).unwrap();
        writeln!(out, "}}").unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hypr() -> NHypr {
        NHypr::new_with_params(
            &NAx2::new(
                &NPnt::new(0.0, 0.0, 0.0),
                &NDir::new(0.0, 0.0, 1.0).unwrap(),
                &NDir::new(1.0, 0.0, 0.0).unwrap(),
            )
            .unwrap(),
            2.0,
            1.0,
        )
        .unwrap()
    }

    #[test]
    fn test_new() {
        let h = NHypr::new();
        assert_eq!(h.major_radius(), f64::MAX);
        assert_eq!(h.minor_radius(), f64::MIN);
    }

    #[test]
    fn test_new_with_params() {
        let h = hypr();
        assert_eq!(h.major_radius(), 2.0);
        assert_eq!(h.minor_radius(), 1.0);
        assert!(matches!(
            NHypr::new_with_params(
                &NAx2::new(
                    &NPnt::new(0.0, 0.0, 0.0),
                    &NDir::new(0.0, 0.0, 1.0).unwrap(),
                    &NDir::new(1.0, 0.0, 0.0).unwrap()
                )
                .unwrap(),
                -1.0,
                1.0
            ),
            Err(NErrors::ConstructionError)
        ));
    }

    #[test]
    fn test_setters() {
        let mut h = hypr();
        h.set_major_radius(3.0).unwrap();
        assert_eq!(h.major_radius(), 3.0);
        h.set_minor_radius(2.0).unwrap();
        assert_eq!(h.minor_radius(), 2.0);
        let p = NPnt::new(1.0, 2.0, 3.0);
        h.set_location(&p);
        assert_eq!(h.location(), &p);
    }

    #[test]
    fn test_asymptotes() {
        let h = hypr();
        let a1 = h.asymptote1().unwrap();
        let a2 = h.asymptote2().unwrap();
        assert!(
            a1.direction()
                .is_parallel(&NDir::new(2.0, 1.0, 0.0).unwrap(), 1e-9)
        );
        assert!(
            a2.direction()
                .is_parallel(&NDir::new(2.0, -1.0, 0.0).unwrap(), 1e-9)
        );
    }

    #[test]
    fn test_foci() {
        let h = hypr();
        let f1 = h.focus1();
        let f2 = h.focus2();
        assert!(f1.distance(&NPnt::new(5.0f64.sqrt(), 0.0, 0.0)) < 1e-9);
        assert!(f2.distance(&NPnt::new(-5.0f64.sqrt(), 0.0, 0.0)) < 1e-9);
    }

    #[test]
    fn test_eccentricity() {
        let h = hypr();
        assert!((h.eccentricity().unwrap() - 5.0f64.sqrt() / 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_directrices() {
        let h = hypr();
        let d1 = h.directrix1();
        let d2 = h.directrix2();
        assert!(
            d1.location()
                .distance(&NPnt::new(4.0 / 5.0f64.sqrt(), 0.0, 0.0))
                < 1e-9
        );
        assert!(
            d2.location()
                .distance(&NPnt::new(-4.0 / 5.0f64.sqrt(), 0.0, 0.0))
                < 1e-9
        );
    }

    #[test]
    fn test_transformations() {
        let h = hypr();
        let mut h_scaled = h.scaled(&NPnt::new(0.0, 0.0, 0.0), 2.0);
        assert_eq!(h_scaled.major_radius(), 4.0);
        assert_eq!(h_scaled.minor_radius(), 2.0);

        let mut h_mirrored = h.mirrored_pnt(&NPnt::new(1.0, 0.0, 0.0));
        assert!(h_mirrored.location().distance(&NPnt::new(2.0, 0.0, 0.0)) < 1e-9);
    }

    #[test]
    fn test_dump_json() {
        let h = hypr();
        let mut output = Vec::new();
        h.dump_json(&mut output, -1);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NHypr\""));
        assert!(json.contains("\"major_radius\": 2"));
    }
}
