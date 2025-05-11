use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::nature_common::prelude::*;

use crate::nature_gp::{NAx2d, NAx22d, NDir2d, NPoint2d, NTrsf2d, NVec2d};

// Trait to define the behavior of an ellipse in 2D space
pub trait Elips2d {
    fn new(
        major_axis: &NAx2d,
        major_radius: f64,
        minor_radius: f64,
        is_sense: bool,
    ) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn new_with_ax22d(a: &NAx22d, major_radius: f64, minor_radius: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn set_location(&mut self, p: &NPoint2d);
    fn set_major_radius(&mut self, major_radius: f64) -> Result<(), NErrors>;
    fn set_minor_radius(&mut self, minor_radius: f64) -> Result<(), NErrors>;
    fn set_axis(&mut self, a: &NAx22d);
    fn set_x_axis(&mut self, a: &NAx2d);
    fn set_y_axis(&mut self, a: &NAx2d);
    fn area(&self) -> f64;
    fn coefficients(&self) -> (f64, f64, f64, f64, f64, f64);
    fn directrix1(&self) -> Result<NAx2d, NErrors>;
    fn directrix2(&self) -> Result<NAx2d, NErrors>;
    fn eccentricity(&self) -> f64;
    fn focal(&self) -> f64;
    fn focus1(&self) -> NPoint2d;
    fn focus2(&self) -> NPoint2d;
    fn location(&self) -> &NPoint2d;
    fn major_radius(&self) -> f64;
    fn minor_radius(&self) -> f64;
    fn parameter(&self) -> f64;
    fn axis(&self) -> &NAx22d;
    fn x_axis(&self) -> NAx2d;
    fn y_axis(&self) -> NAx2d;
    fn reverse(&mut self);
    fn reversed(&self) -> Self;
    fn is_direct(&self) -> bool;
    fn mirror_point3d(&mut self, p: &NPoint2d);
    fn mirrored_point3d(&self, p: &NPoint2d) -> Self;
    fn mirror_ax2d(&mut self, a: &NAx2d);
    fn mirrored_ax2d(&self, a: &NAx2d) -> Self;
    fn rotate(&mut self, p: &NPoint2d, angle: f64);
    fn rotated(&self, p: &NPoint2d, angle: f64) -> Self;
    fn scale(&mut self, p: &NPoint2d, s: f64);
    fn scaled(&self, p: &NPoint2d, s: f64) -> Self;
    fn transform(&mut self, t: &NTrsf2d);
    fn transformed(&self, t: &NTrsf2d) -> Self;
    fn translate_vec(&mut self, v: &NVec2d);
    fn translated_vec(&self, v: &NVec2d) -> Self;
    fn translate_point3d(&mut self, p1: &NPoint2d, p2: &NPoint2d);
    fn translated_point3d(&self, p1: &NPoint2d, p2: &NPoint2d) -> Self;
}

// Struct representing an ellipse in 2D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NElips2d {
    pos: NAx22d,
    major_radius: f64,
    minor_radius: f64,
}

impl Elips2d for NElips2d {
    fn new(
        major_axis: &NAx2d,
        major_radius: f64,
        minor_radius: f64,
        is_sense: bool,
    ) -> Result<Self, NErrors> {
        if minor_radius < 0.0 || major_radius < minor_radius {
            return Err(NErrors::ConstructionError);
        }
        Ok(NElips2d {
            pos: NAx22d::new_from_ax2d(major_axis, is_sense)?,
            major_radius,
            minor_radius,
        })
    }

    fn new_with_ax22d(a: &NAx22d, major_radius: f64, minor_radius: f64) -> Result<Self, NErrors> {
        if minor_radius < 0.0 || major_radius < minor_radius {
            return Err(NErrors::ConstructionError);
        }
        Ok(NElips2d {
            pos: a.clone(),
            major_radius,
            minor_radius,
        })
    }

    fn set_location(&mut self, p: &NPoint2d) {
        self.pos.set_location(p);
    }

    fn set_major_radius(&mut self, major_radius: f64) -> Result<(), NErrors> {
        if major_radius < self.minor_radius {
            return Err(NErrors::ConstructionError);
        }
        self.major_radius = major_radius;
        Ok(())
    }

    fn set_minor_radius(&mut self, minor_radius: f64) -> Result<(), NErrors> {
        if minor_radius < 0.0 || self.major_radius < minor_radius {
            return Err(NErrors::ConstructionError);
        }
        self.minor_radius = minor_radius;
        Ok(())
    }

    fn set_axis(&mut self, a: &NAx22d) {
        self.pos.set_axis(a);
    }

    fn set_x_axis(&mut self, a: &NAx2d) {
        self.pos.set_x_axis(a);
    }

    fn set_y_axis(&mut self, a: &NAx2d) {
        self.pos.set_y_axis(a);
    }

    fn area(&self) -> f64 {
        PI * self.major_radius * self.minor_radius
    }

    fn coefficients(&self) -> (f64, f64, f64, f64, f64, f64) {
        let d_min = self.minor_radius * self.minor_radius;
        let d_maj = self.major_radius * self.major_radius;
        if d_min <= crate::gp::NGP::resolution() && d_maj <= crate::gp::NGP::resolution() {
            return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        }

        let t = NTrsf2d::new_transformation(&self.pos.x_axis()).expect("Invalid transformation");
        let t11 = t.value(1, 1);
        let t12 = t.value(1, 2);
        let t13 = t.value(1, 3);

        if d_min <= crate::gp::NGP::resolution() {
            let a = t11 * t11;
            let b = t12 * t12;
            let c = t11 * t12;
            let d = t11 * t13;
            let e = t12 * t13;
            let f = t13 * t13 - d_maj;
            (a, b, c, d, e, f)
        } else {
            let t21 = t.value(2, 1);
            let t22 = t.value(2, 2);
            let t23 = t.value(2, 3);
            let a = (t11 * t11 / d_maj) + (t21 * t21 / d_min);
            let b = (t12 * t12 / d_maj) + (t22 * t22 / d_min);
            let c = (t11 * t12 / d_maj) + (t21 * t22 / d_min);
            let d = (t11 * t13 / d_maj) + (t21 * t23 / d_min);
            let e = (t12 * t13 / d_maj) + (t22 * t23 / d_min);
            let f = (t13 * t13 / d_maj) + (t23 * t23 / d_min) - 1.0;
            (a, b, c, d, e, f)
        }
    }

    fn directrix1(&self) -> Result<NAx2d, NErrors> {
        let e = self.eccentricity();
        if e <= crate::gp::NGP::resolution() {
            return Err(NErrors::ConstructionError);
        }
        let orig = self
            .pos
            .x_direction()
            .xy()
            .multiplied(self.major_radius / e)
            .added(&self.pos.location().xy());
        Ok(NAx2d::new(
            &NPoint2d::new_from_xy(&orig),
            &self.pos.y_direction(),
        ))
    }

    fn directrix2(&self) -> Result<NAx2d, NErrors> {
        let e = self.eccentricity();
        if e <= crate::gp::NGP::resolution() {
            return Err(NErrors::ConstructionError);
        }
        let orig = self
            .pos
            .x_direction()
            .xy()
            .multiplied(-self.major_radius / e)
            .added(&self.pos.location().xy());
        Ok(NAx2d::new(
            &NPoint2d::new_from_xy(&orig),
            &self.pos.y_direction(),
        ))
    }

    fn eccentricity(&self) -> f64 {
        if self.major_radius == 0.0 {
            0.0
        } else {
            ((self.major_radius * self.major_radius - self.minor_radius * self.minor_radius).sqrt())
                / self.major_radius
        }
    }

    fn focal(&self) -> f64 {
        2.0 * (self.major_radius * self.major_radius - self.minor_radius * self.minor_radius).sqrt()
    }

    fn focus1(&self) -> NPoint2d {
        let c =
            (self.major_radius * self.major_radius - self.minor_radius * self.minor_radius).sqrt();
        let p = self.pos.location();
        let d = self.pos.x_direction();
        NPoint2d::new(p.x() + c * d.x(), p.y() + c * d.y())
    }

    fn focus2(&self) -> NPoint2d {
        let c =
            (self.major_radius * self.major_radius - self.minor_radius * self.minor_radius).sqrt();
        let p = self.pos.location();
        let d = self.pos.x_direction();
        NPoint2d::new(p.x() - c * d.x(), p.y() - c * d.y())
    }

    fn location(&self) -> &NPoint2d {
        self.pos.location()
    }

    fn major_radius(&self) -> f64 {
        self.major_radius
    }

    fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    fn parameter(&self) -> f64 {
        if self.major_radius == 0.0 {
            0.0
        } else {
            (self.minor_radius * self.minor_radius) / self.major_radius
        }
    }

    fn axis(&self) -> &NAx22d {
        &self.pos
    }

    fn x_axis(&self) -> NAx2d {
        self.pos.x_axis()
    }

    fn y_axis(&self) -> NAx2d {
        self.pos.y_axis()
    }

    fn reverse(&mut self) {
        let mut y_dir = self.pos.y_direction();
        y_dir.reverse();
        self.pos.set_axis(
            &NAx22d::new(self.pos.location(), &self.pos.x_direction(), &y_dir)
                .expect("Invalid axis"),
        );
    }

    fn reversed(&self) -> Self {
        let mut result = self.clone();
        result.reverse();
        result
    }

    fn is_direct(&self) -> bool {
        self.pos.x_direction().crossed(&self.pos.y_direction()) >= 0.0
    }

    fn mirror_point3d(&mut self, p: &NPoint2d) {
        self.pos.mirror_point3d(p);
    }

    fn mirrored_point3d(&self, p: &NPoint2d) -> Self {
        let mut result = self.clone();
        result.mirror_point3d(p);
        result
    }

    fn mirror_ax2d(&mut self, a: &NAx2d) {
        self.pos.mirror_ax2d(a);
    }

    fn mirrored_ax2d(&self, a: &NAx2d) -> Self {
        let mut result = self.clone();
        result.mirror_ax2d(a);
        result
    }

    fn rotate(&mut self, p: &NPoint2d, angle: f64) {
        self.pos.rotate(p, angle);
    }

    fn rotated(&self, p: &NPoint2d, angle: f64) -> Self {
        let mut result = self.clone();
        result.rotate(p, angle);
        result
    }

    fn scale(&mut self, p: &NPoint2d, s: f64) {
        self.major_radius *= s.abs();
        self.minor_radius *= s.abs();
        self.pos.scale(p, s);
    }

    fn scaled(&self, p: &NPoint2d, s: f64) -> Self {
        let mut result = self.clone();
        result.scale(p, s);
        result
    }

    fn transform(&mut self, t: &NTrsf2d) {
        let scale = t.scale_factor().abs();
        self.major_radius *= scale;
        self.minor_radius *= scale;
        self.pos.transform(t);
    }

    fn transformed(&self, t: &NTrsf2d) -> Self {
        let mut result = self.clone();
        result.transform(t);
        result
    }

    fn translate_vec(&mut self, v: &NVec2d) {
        self.pos.translate_vec(v);
    }

    fn translated_vec(&self, v: &NVec2d) -> Self {
        let mut result = self.clone();
        result.translate_vec(v);
        result
    }

    fn translate_point3d(&mut self, p1: &NPoint2d, p2: &NPoint2d) {
        self.pos.translate_point3d(p1, p2);
    }

    fn translated_point3d(&self, p1: &NPoint2d, p2: &NPoint2d) -> Self {
        let mut result = self.clone();
        result.translate_point3d(p1, p2);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn elips2d(major_axis: &NAx2d, major: f64, minor: f64, is_sense: bool) -> NElips2d {
        NElips2d::new(major_axis, major, minor, is_sense).expect("Invalid ellipse")
    }

    fn ax2d(p: &NPoint2d, d: &NDir2d) -> NAx2d {
        NAx2d::new(p, d)
    }

    #[test]
    fn test_new() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let elips = NElips2d::new(&axis, 5.0, 3.0, true).unwrap();
        assert_eq!(elips.major_radius(), 5.0);
        assert_eq!(elips.minor_radius(), 3.0);
        assert_eq!(elips.location(), &p);

        assert!(matches!(
            NElips2d::new(&axis, 2.0, 3.0, true),
            Err(NErrors::ConstructionError)
        ));
        assert!(matches!(
            NElips2d::new(&axis, 2.0, -1.0, true),
            Err(NErrors::ConstructionError)
        ));
    }

    #[test]
    fn test_new_with_ax22d() {
        let p = NPoint2d::new(0.0, 0.0);
        let x_dir = NDir2d::new(1.0, 0.0).unwrap();
        let y_dir = NDir2d::new(0.0, 1.0).unwrap();
        let ax22d = NAx22d::new(&p, &x_dir, &y_dir).unwrap();
        let elips = NElips2d::new_with_ax22d(&ax22d, 5.0, 3.0).unwrap();
        assert_eq!(elips.major_radius(), 5.0);
        assert_eq!(elips.minor_radius(), 3.0);
        assert_eq!(elips.axis(), &ax22d);
    }

    #[test]
    fn test_setters() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let mut elips = elips2d(&axis, 5.0, 3.0, true);

        elips.set_major_radius(6.0).unwrap();
        assert_eq!(elips.major_radius(), 6.0);
        assert!(matches!(
            elips.set_major_radius(2.0),
            Err(NErrors::ConstructionError)
        ));

        elips.set_minor_radius(2.0).unwrap();
        assert_eq!(elips.minor_radius(), 2.0);
        assert!(matches!(
            elips.set_minor_radius(7.0),
            Err(NErrors::ConstructionError)
        ));
        assert!(matches!(
            elips.set_minor_radius(-1.0),
            Err(NErrors::ConstructionError)
        ));

        let new_p = NPoint2d::new(1.0, 1.0);
        elips.set_location(&new_p);
        assert_eq!(elips.location(), &new_p);

        let new_axis = ax2d(&NPoint2d::new(0.0, 0.0), &NDir2d::new(0.0, 1.0).unwrap());
        elips.set_x_axis(&new_axis);
        assert_eq!(elips.x_axis(), new_axis);
    }

    #[test]
    fn test_area() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let elips = elips2d(&axis, 5.0, 3.0, true);
        assert!((elips.area() - PI * 5.0 * 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_coefficients() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let elips = elips2d(&axis, 5.0, 3.0, true);
        let (a, b, c, d, e, f) = elips.coefficients();
        assert!((a - 1.0 / (5.0 * 5.0)).abs() < 1e-10);
        assert!((b - 1.0 / (3.0 * 3.0)).abs() < 1e-10);
        assert_eq!(c, 0.0);
        assert_eq!(d, 0.0);
        assert_eq!(e, 0.0);
        assert!((f - -1.0).abs() < 1e-10);

        let degenerate = elips2d(&axis, 0.0, 0.0, true);
        let (a, b, c, d, e, f) = degenerate.coefficients();
        assert_eq!((a, b, c, d, e, f), (0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    }

    #[test]
    fn test_eccentricity() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let elips = elips2d(&axis, 5.0, 3.0, true);
        let e = ((5.0 * 5.0 - 3.0 * 3.0).sqrt()) / 5.0;
        assert!((elips.eccentricity() - e).abs() < 1e-10);

        let circle = elips2d(&axis, 5.0, 5.0, true);
        assert_eq!(circle.eccentricity(), 0.0);
    }

    #[test]
    fn test_focal() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let elips = elips2d(&axis, 5.0, 3.0, true);
        let f = 2.0 * (5.0 * 5.0 - 3.0 * 3.0).sqrt();
        assert!((elips.focal() - f).abs() < 1e-10);
    }

    #[test]
    fn test_focus() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let elips = elips2d(&axis, 5.0, 3.0, true);
        let c = (5.0 * 5.0 - 3.0 * 3.0).sqrt();
        let f1 = elips.focus1();
        let f2 = elips.focus2();
        assert!((f1.x() - c).abs() < 1e-10);
        assert_eq!(f1.y(), 0.0);
        assert!((f2.x() + c).abs() < 1e-10);
        assert_eq!(f2.y(), 0.0);
    }

    #[test]
    fn test_directrix() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let elips = elips2d(&axis, 5.0, 3.0, true);
        let e = elips.eccentricity();
        let d1 = elips.directrix1().unwrap();
        let d2 = elips.directrix2().unwrap();
        let d1_loc = d1.location();
        let d2_loc = d2.location();
        assert!((d1_loc.x() - 5.0 / e).abs() < 1e-10);
        assert!((d2_loc.x() + 5.0 / e).abs() < 1e-10);
        assert_eq!(d1_loc.y(), 0.0);
        assert_eq!(d2_loc.y(), 0.0);

        let circle = elips2d(&axis, 5.0, 5.0, true);
        assert!(matches!(
            circle.directrix1(),
            Err(NErrors::ConstructionError)
        ));
        assert!(matches!(
            circle.directrix2(),
            Err(NErrors::ConstructionError)
        ));
    }

    #[test]
    fn test_parameter() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let elips = elips2d(&axis, 5.0, 3.0, true);
        assert!((elips.parameter() - (3.0 * 3.0) / 5.0).abs() < 1e-10);

        let zero = elips2d(&axis, 0.0, 0.0, true);
        assert_eq!(zero.parameter(), 0.0);
    }

    #[test]
    fn test_axes() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let elips = elips2d(&axis, 5.0, 3.0, true);
        let x_axis = elips.x_axis();
        let y_axis = elips.y_axis();
        assert_eq!(x_axis.location(), &NPoint2d::new(0.0, 0.0));
        assert_eq!(y_axis.location(), &NPoint2d::new(0.0, 0.0));
        assert_eq!(x_axis.direction().coords(), (1.0, 0.0));
        assert_eq!(y_axis.direction().coords(), (0.0, 1.0));
    }

    #[test]
    fn test_reverse() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let mut elips = elips2d(&axis, 5.0, 3.0, true);
        elips.reverse();
        assert_eq!(elips.y_axis().direction().coords(), (0.0, -1.0));

        let reversed = elips.reversed();
        assert_eq!(reversed.y_axis().direction().coords(), (0.0, 1.0));
        assert!(reversed.is_direct());
    }

    #[test]
    fn test_is_direct() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let elips = elips2d(&axis, 5.0, 3.0, true);
        assert!(elips.is_direct());

        let elips_inv = elips2d(&axis, 5.0, 3.0, false);
        assert!(!elips_inv.is_direct());
    }

    #[test]
    fn test_mirror() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let mut elips = elips2d(&axis, 5.0, 3.0, true);
        let p_mirror = NPoint2d::new(1.0, 1.0);
        elips.mirror_point3d(&p_mirror);
        assert_eq!(elips.location(), &NPoint2d::new(2.0, 2.0));

        let mirrored = elips.mirrored_point3d(&p_mirror);
        assert_eq!(mirrored.location(), &NPoint2d::new(0.0, 0.0));

        let ax_mirror = ax2d(&NPoint2d::new(0.0, 0.0), &NDir2d::new(0.0, 1.0).unwrap());
        elips.mirror_ax2d(&ax_mirror);
        assert_eq!(elips.location(), &NPoint2d::new(2.0, -2.0));
    }

    #[test]
    fn test_rotate() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let mut elips = elips2d(&axis, 5.0, 3.0, true);
        let p_rot = NPoint2d::new(0.0, 0.0);
        elips.rotate(&p_rot, PI / 2.0);
        let x_dir = elips.x_axis().direction();
        assert!((x_dir.x() - 0.0).abs() < 1e-10);
        assert!((x_dir.y() - 1.0).abs() < 1e-10);

        let rotated = elips.rotated(&p_rot, -PI / 2.0);
        assert_eq!(rotated.x_axis().direction().coords(), (1.0, 0.0));
    }

    #[test]
    fn test_scale() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let mut elips = elips2d(&axis, 5.0, 3.0, true);
        let p_scale = NPoint2d::new(1.0, 1.0);
        elips.scale(&p_scale, -2.0);
        assert_eq!(elips.major_radius(), 10.0);
        assert_eq!(elips.minor_radius(), 6.0);
        assert_eq!(elips.location(), &NPoint2d::new(3.0, 3.0));

        let scaled = elips.scaled(&p_scale, 0.5);
        assert_eq!(scaled.major_radius(), 5.0);
        assert_eq!(scaled.minor_radius(), 3.0);
        assert_eq!(scaled.location(), &NPoint2d::new(2.0, 2.0));
    }

    #[test]
    fn test_transform() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let mut elips = elips2d(&axis, 5.0, 3.0, true);
        let trsf = NTrsf2d::new_scale(&NPoint2d::new(0.0, 0.0), -2.0).unwrap();
        elips.transform(&trsf);
        assert_eq!(elips.major_radius(), 10.0);
        assert_eq!(elips.minor_radius(), 6.0);
        assert_eq!(elips.location(), &NPoint2d::new(0.0, 0.0));

        let trsf = NTrsf2d::new_rotation(&NPoint2d::new(0.0, 0.0), PI / 2.0).unwrap();
        let transformed = elips.transformed(&trsf);
        let x_dir = transformed.x_axis().direction();
        assert!((x_dir.x() - 0.0).abs() < 1e-10);
        assert!((x_dir.y() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_translate() {
        let p = NPoint2d::new(0.0, 0.0);
        let d = NDir2d::new(1.0, 0.0).unwrap();
        let axis = ax2d(&p, &d);
        let mut elips = elips2d(&axis, 5.0, 3.0, true);
        let v = NVec2d::new(1.0, 2.0);
        elips.translate_vec(&v);
        assert_eq!(elips.location(), &NPoint2d::new(1.0, 2.0));

        let p1 = NPoint2d::new(1.0, 1.0);
        let p2 = NPoint2d::new(2.0, 3.0);
        let translated = elips.translated_point3d(&p1, &p2);
        assert_eq!(translated.location(), &NPoint2d::new(2.0, 3.0));
    }
}
