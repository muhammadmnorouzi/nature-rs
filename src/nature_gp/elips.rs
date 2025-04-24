use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::{
    gp::{NAx1, NAx2, NPnt, NTrsf, NVec},
    nature_errors::NErrors,
};

// Trait to define the behavior of an ellipse in 3D space
pub trait Elips {
    fn new(a2: &NAx2, major_radius: f64, minor_radius: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn set_axis(&mut self, a1: &NAx1) -> Result<(), NErrors>;
    fn set_location(&mut self, p: &NPnt);
    fn set_major_radius(&mut self, major_radius: f64) -> Result<(), NErrors>;
    fn set_minor_radius(&mut self, minor_radius: f64) -> Result<(), NErrors>;
    fn set_position(&mut self, a2: &NAx2);
    fn area(&self) -> f64;
    fn axis(&self) -> NAx1;
    fn directrix1(&self) -> Result<NAx1, NErrors>;
    fn directrix2(&self) -> Result<NAx1, NErrors>;
    fn eccentricity(&self) -> f64;
    fn focal(&self) -> f64;
    fn focus1(&self) -> NPnt;
    fn focus2(&self) -> NPnt;
    fn location(&self) -> &NPnt;
    fn major_radius(&self) -> f64;
    fn minor_radius(&self) -> f64;
    fn parameter(&self) -> f64;
    fn position(&self) -> &NAx2;
    fn x_axis(&self) -> NAx1;
    fn y_axis(&self) -> NAx1;
    fn mirror_pnt(&mut self, p: &NPnt);
    fn mirrored_pnt(&self, p: &NPnt) -> Self;
    fn mirror_ax1(&mut self, a1: &NAx1);
    fn mirrored_ax1(&self, a1: &NAx1) -> Self;
    fn mirror_ax2(&mut self, a2: &NAx2);
    fn mirrored_ax2(&self, a2: &NAx2) -> Self;
    fn rotate(&mut self, a1: &NAx1, angle: f64);
    fn rotated(&self, a1: &NAx1, angle: f64) -> Self;
    fn scale(&mut self, p: &NPnt, s: f64);
    fn scaled(&self, p: &NPnt, s: f64) -> Self;
    fn transform(&mut self, t: &NTrsf);
    fn transformed(&self, t: &NTrsf) -> Self;
    fn translate_vec(&mut self, v: &NVec);
    fn translated_vec(&self, v: &NVec) -> Self;
    fn translate_pnts(&mut self, p1: &NPnt, p2: &NPnt);
    fn translated_pnts(&self, p1: &NPnt, p2: &NPnt) -> Self;
}

// Struct representing an ellipse in 3D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NElips {
    pos: NAx2,
    major_radius: f64,
    minor_radius: f64,
}

impl Elips for NElips {
    fn new(a2: &NAx2, major_radius: f64, minor_radius: f64) -> Result<Self, NErrors> {
        if minor_radius < 0.0 || major_radius < minor_radius {
            return Err(NErrors::ConstructionError);
        }
        Ok(NElips {
            pos: a2.clone(),
            major_radius,
            minor_radius,
        })
    }

    fn set_axis(&mut self, a1: &NAx1) -> Result<(), NErrors> {
        self.pos.set_axis(a1)?;
        Ok(())
    }

    fn set_location(&mut self, p: &NPnt) {
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

    fn set_position(&mut self, a2: &NAx2) {
        self.pos = a2.clone();
    }

    fn area(&self) -> f64 {
        PI * self.major_radius * self.minor_radius
    }

    fn axis(&self) -> NAx1 {
        self.pos.axis()
    }

    fn directrix1(&self) -> Result<NAx1, NErrors> {
        let e = self.eccentricity();
        if e <= crate::gp::NGP::resolution() {
            return Err(NErrors::ConstructionError);
        }
        let orig = self
            .pos
            .x_direction()
            .xyz()
            .multiplied(self.major_radius / e)
            .added(&self.pos.location().xyz());
        Ok(NAx1::new(
            &NPnt::new_from_xyz(&orig),
            &self.pos.y_direction(),
        ))
    }

    fn directrix2(&self) -> Result<NAx1, NErrors> {
        let e = self.eccentricity();
        if e <= crate::gp::NGP::resolution() {
            return Err(NErrors::ConstructionError);
        }
        let orig = self
            .pos
            .x_direction()
            .xyz()
            .multiplied(-self.major_radius / e)
            .added(&self.pos.location().xyz());
        Ok(NAx1::new(
            &NPnt::new_from_xyz(&orig),
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

    fn focus1(&self) -> NPnt {
        let c =
            (self.major_radius * self.major_radius - self.minor_radius * self.minor_radius).sqrt();
        let p = self.pos.location();
        let d = self.pos.x_direction();
        NPnt::new(p.x() + c * d.x(), p.y() + c * d.y(), p.z() + c * d.z())
    }

    fn focus2(&self) -> NPnt {
        let c =
            (self.major_radius * self.major_radius - self.minor_radius * self.minor_radius).sqrt();
        let p = self.pos.location();
        let d = self.pos.x_direction();
        NPnt::new(p.x() - c * d.x(), p.y() - c * d.y(), p.z() - c * d.z())
    }

    fn location(&self) -> &NPnt {
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

    fn position(&self) -> &NAx2 {
        &self.pos
    }

    fn x_axis(&self) -> NAx1 {
        NAx1::new(self.pos.location(), &self.pos.x_direction())
    }

    fn y_axis(&self) -> NAx1 {
        NAx1::new(self.pos.location(), &self.pos.y_direction())
    }

    fn mirror_pnt(&mut self, p: &NPnt) {
        self.pos.mirror_pnt(p);
    }

    fn mirrored_pnt(&self, p: &NPnt) -> Self {
        let mut result = self.clone();
        result.mirror_pnt(p);
        result
    }

    fn mirror_ax1(&mut self, a1: &NAx1) {
        self.pos.mirror_ax1(a1);
    }

    fn mirrored_ax1(&self, a1: &NAx1) -> Self {
        let mut result = self.clone();
        result.mirror_ax1(a1);
        result
    }

    fn mirror_ax2(&mut self, a2: &NAx2) {
        self.pos.mirror_ax2(a2);
    }

    fn mirrored_ax2(&self, a2: &NAx2) -> Self {
        let mut result = self.clone();
        result.mirror_ax2(a2);
        result
    }

    fn rotate(&mut self, a1: &NAx1, angle: f64) {
        self.pos.rotate(a1, angle);
    }

    fn rotated(&self, a1: &NAx1, angle: f64) -> Self {
        let mut result = self.clone();
        result.rotate(a1, angle);
        result
    }

    fn scale(&mut self, p: &NPnt, s: f64) {
        self.major_radius *= s.abs();
        self.minor_radius *= s.abs();
        self.pos.scale(p, s);
    }

    fn scaled(&self, p: &NPnt, s: f64) -> Self {
        let mut result = self.clone();
        result.scale(p, s);
        result
    }

    fn transform(&mut self, t: &NTrsf) {
        self.major_radius *= t.scale_factor().abs();
        self.minor_radius *= t.scale_factor().abs();
        self.pos.transform(t);
    }

    fn transformed(&self, t: &NTrsf) -> Self {
        let mut result = self.clone();
        result.transform(t);
        result
    }

    fn translate_vec(&mut self, v: &NVec) {
        self.pos.translate_vec(v);
    }

    fn translated_vec(&self, v: &NVec) -> Self {
        let mut result = self.clone();
        result.translate_vec(v);
        result
    }

    fn translate_pnts(&mut self, p1: &NPnt, p2: &NPnt) {
        self.pos.translate_pnts(p1, p2);
    }

    fn translated_pnts(&self, p1: &NPnt, p2: &NPnt) -> Self {
        let mut result = self.clone();
        result.translate_pnts(p1, p2);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn elips(a2: &NAx2, major: f64, minor: f64) -> NElips {
        NElips::new(a2, major, minor).expect("Invalid ellipse")
    }

    #[test]
    fn test_new() {
        let pos = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let elips = NElips::new(&pos, 5.0, 3.0).unwrap();
        assert_eq!(elips.major_radius(), 5.0);
        assert_eq!(elips.minor_radius(), 3.0);
        assert_eq!(elips.position(), &pos);

        assert!(matches!(
            NElips::new(&pos, 2.0, 3.0),
            Err(NErrors::ConstructionError)
        ));
        assert!(matches!(
            NElips::new(&pos, 2.0, -1.0),
            Err(NErrors::ConstructionError)
        ));
    }

    #[test]
    fn test_setters() {
        let pos = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let mut elips = elips(&pos, 5.0, 3.0);

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

        let new_pos = NAx2::new(
            &NPnt::new(1.0, 1.0, 1.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        elips.set_position(&new_pos);
        assert_eq!(elips.position(), &new_pos);

        let p = NPnt::new(2.0, 2.0, 2.0);
        elips.set_location(&p);
        assert_eq!(elips.location(), &p);
    }

    #[test]
    fn test_area() {
        let pos = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let elips = elips(&pos, 5.0, 3.0);
        assert!((elips.area() - PI * 5.0 * 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_eccentricity() {
        let pos = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let elips = elips(&pos, 5.0, 3.0);
        let e = ((5.0 * 5.0 - 3.0 * 3.0).sqrt()) / 5.0;
        assert!((elips.eccentricity() - e).abs() < 1e-10);

        let circle = elips(&pos, 5.0, 5.0);
        assert_eq!(circle.eccentricity(), 0.0);
    }

    #[test]
    fn test_focal() {
        let pos = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let elips = elips(&pos, 5.0, 3.0);
        let f = 2.0 * (5.0 * 5.0 - 3.0 * 3.0).sqrt();
        assert!((elips.focal() - f).abs() < 1e-10);
    }

    #[test]
    fn test_focus() {
        let pos = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let elips = elips(&pos, 5.0, 3.0);
        let c = (5.0 * 5.0 - 3.0 * 3.0).sqrt();
        let f1 = elips.focus1();
        let f2 = elips.focus2();
        assert!((f1.x() - c).abs() < 1e-10);
        assert_eq!(f1.y(), 0.0);
        assert_eq!(f1.z(), 0.0);
        assert!((f2.x() + c).abs() < 1e-10);
        assert_eq!(f2.y(), 0.0);
        assert_eq!(f2.z(), 0.0);
    }

    #[test]
    fn test_directrix() {
        let pos = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let elips = elips(&pos, 5.0, 3.0);
        let e = elips.eccentricity();
        let d1 = elips.directrix1().unwrap();
        let d2 = elips.directrix2().unwrap();
        let d1_loc = d1.location();
        let d2_loc = d2.location();
        assert!((d1_loc.x() - 5.0 / e).abs() < 1e-10);
        assert!((d2_loc.x() + 5.0 / e).abs() < 1e-10);
        assert_eq!(d1_loc.y(), 0.0);
        assert_eq!(d2_loc.y(), 0.0);
        assert_eq!(d1_loc.z(), 0.0);
        assert_eq!(d2_loc.z(), 0.0);

        let circle = elips(&pos, 5.0, 5.0);
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
        let pos = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let elips = elips(&pos, 5.0, 3.0);
        assert!((elips.parameter() - (3.0 * 3.0) / 5.0).abs() < 1e-10);

        let zero = elips(&pos, 0.0, 0.0);
        assert_eq!(zero.parameter(), 0.0);
    }

    #[test]
    fn test_axes() {
        let pos = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let elips = elips(&pos, 5.0, 3.0);
        let x_axis = elips.x_axis();
        let y_axis = elips.y_axis();
        assert_eq!(x_axis.location(), &NPnt::new(0.0, 0.0, 0.0));
        assert_eq!(y_axis.location(), &NPnt::new(0.0, 0.0, 0.0));
        assert_eq!(x_axis.direction().coords(), (1.0, 0.0, 0.0));
        assert_eq!(y_axis.direction().coords(), (0.0, 1.0, 0.0));
    }

    #[test]
    fn test_mirror() {
        let pos = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let mut elips = elips(&pos, 5.0, 3.0);
        let p = NPnt::new(1.0, 1.0, 1.0);
        elips.mirror_pnt(&p);
        let mirrored = elips.mirrored_pnt(&p);
        assert_eq!(elips.location(), &NPnt::new(2.0, 2.0, 2.0));
        assert_eq!(mirrored.location(), &NPnt::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn test_rotate() {
        let pos = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let mut elips = elips(&pos, 5.0, 3.0);
        let axis = NAx1::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
        );
        elips.rotate(&axis, PI / 2.0);
        let x_dir = elips.x_axis().direction();
        assert!((x_dir.x() - 0.0).abs() < 1e-10);
        assert!((x_dir.y() - 1.0).abs() < 1e-10);
        assert_eq!(x_dir.z(), 0.0);

        let rotated = elips.rotated(&axis, -PI / 2.0);
        assert_eq!(rotated.x_axis().direction().coords(), (1.0, 0.0, 0.0));
    }

    #[test]
    fn test_scale() {
        let pos = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let mut elips = elips(&pos, 5.0, 3.0);
        let p = NPnt::new(1.0, 1.0, 1.0);
        elips.scale(&p, -2.0);
        assert_eq!(elips.major_radius(), 10.0);
        assert_eq!(elips.minor_radius(), 6.0);
        assert_eq!(elips.location(), &NPnt::new(3.0, 3.0, 3.0));

        let scaled = elips.scaled(&p, 0.5);
        assert_eq!(scaled.major_radius(), 5.0);
        assert_eq!(scaled.minor_radius(), 3.0);
        assert_eq!(scaled.location(), &NPnt::new(2.0, 2.0, 2.0));
    }

    #[test]
    fn test_transform() {
        let pos = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let mut elips = elips(&pos, 5.0, 3.0);
        let trsf = NTrsf::new_scale(&NPnt::new(0.0, 0.0, 0.0), -2.0).unwrap();
        elips.transform(&trsf);
        assert_eq!(elips.major_radius(), 10.0);
        assert_eq!(elips.minor_radius(), 6.0);
        assert_eq!(elips.location(), &NPnt::new(0.0, 0.0, 0.0));

        let trsf = NTrsf::new_rotation(
            &NAx1::new(
                &NPnt::new(0.0, 0.0, 0.0),
                &NDir::new(0.0, 0.0, 1.0).unwrap(),
            ),
            PI / 2.0,
        )
        .unwrap();
        let transformed = elips.transformed(&trsf);
        let x_dir = transformed.x_axis().direction();
        assert!((x_dir.x() - 0.0).abs() < 1e-10);
        assert!((x_dir.y() - 1.0).abs() < 1e-10);
        assert_eq!(x_dir.z(), 0.0);
    }

    #[test]
    fn test_translate() {
        let pos = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0).unwrap(),
            &NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let mut elips = elips(&pos, 5.0, 3.0);
        let v = NVec::new(1.0, 2.0, 3.0);
        elips.translate_vec(&v);
        assert_eq!(elips.location(), &NPnt::new(1.0, 2.0, 3.0));

        let p1 = NPnt::new(1.0, 1.0, 1.0);
        let p2 = NPnt::new(2.0, 3.0, 4.0);
        let translated = elips.translated_pnts(&p1, &p2);
        assert_eq!(translated.location(), &NPnt::new(2.0, 3.0, 4.0));
    }
}
