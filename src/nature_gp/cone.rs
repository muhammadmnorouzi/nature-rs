use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::{
    nature_errors::NErrors,
    nature_gp::{NAx1, NAx2, NAx3, NPnt, NTrsf, NVec},
};

// Trait to define the behavior of a conical surface in 3D space
pub trait Cone {
    fn new(position: NAx3, semi_angle: f64, radius: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn set_axis(&mut self, axis: &NAx1) -> Result<(), NErrors>;
    fn set_location(&mut self, location: NPnt);
    fn set_position(&mut self, position: NAx3);
    fn set_radius(&mut self, radius: f64) -> Result<(), NErrors>;
    fn set_semi_angle(&mut self, semi_angle: f64) -> Result<(), NErrors>;
    fn apex(&self) -> NPnt;
    fn u_reverse(&mut self);
    fn v_reverse(&mut self);
    fn direct(&self) -> bool;
    fn axis(&self) -> &NAx1;
    fn location(&self) -> &NPnt;
    fn position(&self) -> &NAx3;
    fn ref_radius(&self) -> f64;
    fn semi_angle(&self) -> f64;
    fn x_axis(&self) -> NAx1;
    fn y_axis(&self) -> NAx1;
    fn coefficients(&self) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64);
    fn mirror_pnt(&mut self, point: &NPnt);
    fn mirrored_pnt(&self, point: &NPnt) -> Self;
    fn mirror_ax1(&mut self, axis: &NAx1);
    fn mirrored_ax1(&self, axis: &NAx1) -> Self;
    fn mirror_ax2(&mut self, plane: &NAx2);
    fn mirrored_ax2(&self, plane: &NAx2) -> Self;
    fn rotate(&mut self, axis: &NAx1, angle: f64);
    fn rotated(&self, axis: &NAx1, angle: f64) -> Self;
    fn scale(&mut self, point: &NPnt, factor: f64);
    fn scaled(&self, point: &NPnt, factor: f64) -> Self;
    fn transform(&mut self, transformation: &NTrsf);
    fn transformed(&self, transformation: &NTrsf) -> Self;
    fn translate_vec(&mut self, vector: &NVec);
    fn translated_vec(&self, vector: &NVec) -> Self;
    fn translate_pnts(&mut self, from: &NPnt, to: &NPnt);
    fn translated_pnts(&self, from: &NPnt, to: &NPnt) -> Self;
}

// Struct representing an infinite conical surface in 3D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NCone {
    pos: NAx3,
    radius: f64,
    semi_angle: f64,
}

impl Cone for NCone {
    fn new(position: NAx3, semi_angle: f64, radius: f64) -> Result<Self, NErrors> {
        if radius < 0.0 {
            return Err(NErrors::NegativeRadius);
        }
        let abs_angle = semi_angle.abs();
        if abs_angle <= crate::gp::NGP::resolution()
            || PI * 0.5 - abs_angle <= crate::gp::NGP::resolution()
        {
            return Err(NErrors::InvalidAngle);
        }
        Ok(NCone {
            pos: position,
            radius,
            semi_angle,
        })
    }

    fn set_axis(&mut self, axis: &NAx1) -> Result<(), NErrors> {
        self.pos.set_axis(axis)
    }

    fn set_location(&mut self, location: NPnt) {
        self.pos.set_location(location);
    }

    fn set_position(&mut self, position: NAx3) {
        self.pos = position;
    }

    fn set_radius(&mut self, radius: f64) -> Result<(), NErrors> {
        if radius < 0.0 {
            return Err(NErrors::NegativeRadius);
        }
        self.radius = radius;
        Ok(())
    }

    fn set_semi_angle(&mut self, semi_angle: f64) -> Result<(), NErrors> {
        let abs_angle = semi_angle.abs();
        if abs_angle <= crate::gp::NGP::resolution()
            || PI * 0.5 - abs_angle <= crate::gp::NGP::resolution()
        {
            return Err(NErrors::InvalidAngle);
        }
        self.semi_angle = semi_angle;
        Ok(())
    }

    fn apex(&self) -> NPnt {
        let mut coord = self.pos.direction().xyz();
        coord.multiply(-self.radius / self.semi_angle.tan());
        coord.add(&self.pos.location().xyz());
        NPnt::new(coord.x(), coord.y(), coord.z())
    }

    fn u_reverse(&mut self) {
        self.pos.y_reverse();
    }

    fn v_reverse(&mut self) {
        self.pos.z_reverse();
        self.semi_angle = -self.semi_angle;
    }

    fn direct(&self) -> bool {
        self.pos.direct()
    }

    fn axis(&self) -> &NAx1 {
        self.pos.axis()
    }

    fn location(&self) -> &NPnt {
        self.pos.location()
    }

    fn position(&self) -> &NAx3 {
        &self.pos
    }

    fn ref_radius(&self) -> f64 {
        self.radius
    }

    fn semi_angle(&self) -> f64 {
        self.semi_angle
    }

    fn x_axis(&self) -> NAx1 {
        NAx1::new(self.pos.location().clone(), self.pos.x_direction().clone())
    }

    fn y_axis(&self) -> NAx1 {
        NAx1::new(self.pos.location().clone(), self.pos.y_direction().clone())
    }

    fn coefficients(&self) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
        let trsf = NTrsf::new_transformation(&self.pos).expect("Invalid transformation");
        let k_ang = self.semi_angle.tan();
        let t11 = trsf.value(1, 1);
        let t12 = trsf.value(1, 2);
        let t13 = trsf.value(1, 3);
        let t14 = trsf.value(1, 4);
        let t21 = trsf.value(2, 1);
        let t22 = trsf.value(2, 2);
        let t23 = trsf.value(2, 3);
        let t24 = trsf.value(2, 4);
        let t31 = trsf.value(3, 1) * k_ang;
        let t32 = trsf.value(3, 2) * k_ang;
        let t33 = trsf.value(3, 3) * k_ang;
        let t34 = trsf.value(3, 4) * k_ang;
        let a1 = t11 * t11 + t21 * t21 - t31 * t31;
        let a2 = t12 * t12 + t22 * t22 - t32 * t32;
        let a3 = t13 * t13 + t23 * t23 - t33 * t33;
        let b1 = t11 * t12 + t21 * t22 - t31 * t32;
        let b2 = t11 * t13 + t21 * t23 - t31 * t33;
        let b3 = t12 * t13 + t22 * t23 - t32 * t33;
        let c1 = t11 * t14 + t21 * t24 - t31 * (self.radius + t34);
        let c2 = t12 * t14 + t22 * t24 - t32 * (self.radius + t34);
        let c3 = t13 * t14 + t23 * t24 - t33 * (self.radius + t34);
        let d =
            t14 * t14 + t24 * t24 - self.radius * self.radius - t34 * t34 - 2.0 * self.radius * t34;
        (a1, a2, a3, b1, b2, b3, c1, c2, c3, d)
    }

    fn mirror_pnt(&mut self, point: &NPnt) {
        self.pos.mirror_pnt(point);
    }

    fn mirrored_pnt(&self, point: &NPnt) -> Self {
        let mut result = self.clone();
        result.mirror_pnt(point);
        result
    }

    fn mirror_ax1(&mut self, axis: &NAx1) {
        self.pos.mirror_ax1(axis);
    }

    fn mirrored_ax1(&self, axis: &NAx1) -> Self {
        let mut result = self.clone();
        result.mirror_ax1(axis);
        result
    }

    fn mirror_ax2(&mut self, plane: &NAx2) {
        self.pos.mirror_ax2(plane);
    }

    fn mirrored_ax2(&self, plane: &NAx2) -> Self {
        let mut result = self.clone();
        result.mirror_ax2(plane);
        result
    }

    fn rotate(&mut self, axis: &NAx1, angle: f64) {
        self.pos.rotate(axis, angle);
    }

    fn rotated(&self, axis: &NAx1, angle: f64) -> Self {
        let mut result = self.clone();
        result.rotate(axis, angle);
        result
    }

    fn scale(&mut self, point: &NPnt, factor: f64) {
        self.pos.scale(point, factor);
        self.radius *= factor.abs();
    }

    fn scaled(&self, point: &NPnt, factor: f64) -> Self {
        let mut result = self.clone();
        result.scale(point, factor);
        result
    }

    fn transform(&mut self, transformation: &NTrsf) {
        self.pos.transform(transformation);
        self.radius *= transformation.scale_factor().abs();
    }

    fn transformed(&self, transformation: &NTrsf) -> Self {
        let mut result = self.clone();
        result.transform(transformation);
        result
    }

    fn translate_vec(&mut self, vector: &NVec) {
        self.pos.translate_vec(vector);
    }

    fn translated_vec(&self, vector: &NVec) -> Self {
        let mut result = self.clone();
        result.translate_vec(vector);
        result
    }

    fn translate_pnts(&mut self, from: &NPnt, to: &NPnt) {
        self.pos.translate_pnts(from, to);
    }

    fn translated_pnts(&self, from: &NPnt, to: &NPnt) -> Self {
        let mut result = self.clone();
        result.translate_pnts(from, to);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn cone(
        pos: (f64, f64, f64),
        dir: (f64, f64, f64),
        x_dir: (f64, f64, f64),
        semi_angle: f64,
        radius: f64,
    ) -> NCone {
        NCone::new(
            NAx3::new(
                NPnt::new(pos.0, pos.1, pos.2),
                NDir::new(dir.0, dir.1, dir.2).expect("Invalid direction"),
                NDir::new(x_dir.0, x_dir.1, x_dir.2).expect("Invalid X direction"),
            )
            .expect("Invalid NAx3"),
            semi_angle,
            radius,
        )
        .expect("Invalid NCone")
    }

    #[test]
    fn test_new() {
        let pos = NAx3::new(
            NPnt::new(1.0, 2.0, 3.0),
            NDir::new(0.0, 0.0, 1.0).unwrap(),
            NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let cone = NCone::new(pos.clone(), PI / 4.0, 5.0).unwrap();
        assert_eq!(cone.position(), &pos);
        assert_eq!(cone.semi_angle(), PI / 4.0);
        assert_eq!(cone.ref_radius(), 5.0);

        assert!(matches!(
            NCone::new(pos.clone(), -1.0, 5.0),
            Err(NErrors::InvalidAngle)
        ));
        assert!(matches!(
            NCone::new(pos.clone(), PI / 2.0, 5.0),
            Err(NErrors::InvalidAngle)
        ));
        assert!(matches!(
            NCone::new(pos, 0.0, -1.0),
            Err(NErrors::NegativeRadius)
        ));
    }

    #[test]
    fn test_setters() {
        let mut cone = cone(
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0),
            PI / 4.0,
            5.0,
        );
        let new_axis = NAx1::new(NPnt::new(0.0, 0.0, 0.0), NDir::new(0.0, 1.0, 0.0).unwrap());
        cone.set_axis(&new_axis).unwrap();
        assert_eq!(cone.axis().direction(), &NDir::new(0.0, 1.0, 0.0).unwrap());

        cone.set_location(NPnt::new(1.0, 2.0, 3.0));
        assert_eq!(cone.location(), &NPnt::new(1.0, 2.0, 3.0));

        let new_pos = NAx3::new(
            NPnt::new(4.0, 5.0, 6.0),
            NDir::new(0.0, 0.0, 1.0).unwrap(),
            NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        cone.set_position(new_pos.clone());
        assert_eq!(cone.position(), &new_pos);

        cone.set_radius(10.0).unwrap();
        assert_eq!(cone.ref_radius(), 10.0);
        assert!(matches!(
            cone.set_radius(-1.0),
            Err(NErrors::NegativeRadius)
        ));

        cone.set_semi_angle(-PI / 4.0).unwrap();
        assert_eq!(cone.semi_angle(), -PI / 4.0);
        assert!(matches!(
            cone.set_semi_angle(PI / 2.0),
            Err(NErrors::InvalidAngle)
        ));
    }

    #[test]
    fn test_apex() {
        let cone = cone(
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0),
            PI / 4.0,
            5.0,
        );
        let apex = cone.apex();
        assert!((apex.z() - (-5.0)).abs() < 1e-10);
        assert_eq!(apex.x(), 0.0);
        assert_eq!(apex.y(), 0.0);

        let cone_neg = cone(
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0),
            -PI / 4.0,
            5.0,
        );
        let apex_neg = cone_neg.apex();
        assert!((apex_neg.z() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_reverse() {
        let mut cone = cone(
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0),
            PI / 4.0,
            5.0,
        );
        assert!(cone.direct());
        cone.u_reverse();
        assert_eq!(
            cone.y_axis().direction(),
            &NDir::new(0.0, -1.0, 0.0).unwrap()
        );

        cone.v_reverse();
        assert_eq!(cone.axis().direction(), &NDir::new(0.0, 0.0, -1.0).unwrap());
        assert_eq!(cone.semi_angle(), -PI / 4.0);
    }

    #[test]
    fn test_axes() {
        let cone = cone(
            (1.0, 2.0, 3.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0),
            PI / 4.0,
            5.0,
        );
        assert_eq!(cone.axis().location(), &NPnt::new(1.0, 2.0, 3.0));
        assert_eq!(cone.axis().direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
        assert_eq!(
            cone.x_axis().direction(),
            &NDir::new(1.0, 0.0, 0.0).unwrap()
        );
        assert_eq!(
            cone.y_axis().direction(),
            &NDir::new(0.0, 1.0, 0.0).unwrap()
        );
    }

    #[test]
    fn test_coefficients() {
        let cone = cone(
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0),
            PI / 4.0,
            5.0,
        );
        let (a1, a2, a3, b1, b2, b3, c1, c2, c3, d) = cone.coefficients();
        assert!((a1 - 1.0).abs() < 1e-10);
        assert!((a2 - 1.0).abs() < 1e-10);
        assert!((a3 - (-1.0)).abs() < 1e-10);
        assert_eq!(b1, 0.0);
        assert_eq!(b2, 0.0);
        assert_eq!(b3, 0.0);
        assert_eq!(c1, 0.0);
        assert_eq!(c2, 0.0);
        assert_eq!(c3, 0.0);
        assert!((d - (-25.0)).abs() < 1e-10);
    }

    #[test]
    fn test_mirror_pnt() {
        let mut cone = cone(
            (1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0),
            PI / 4.0,
            5.0,
        );
        let point = NPnt::new(0.0, 0.0, 0.0);
        cone.mirror_pnt(&point);
        assert_eq!(cone.location(), &NPnt::new(-1.0, 0.0, 0.0));
        assert_eq!(cone.semi_angle(), PI / 4.0);
        assert_eq!(cone.ref_radius(), 5.0);
    }

    #[test]
    fn test_mirror_ax1() {
        let mut cone = cone(
            (1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0),
            PI / 4.0,
            5.0,
        );
        let axis = NAx1::new(NPnt::new(0.0, 0.0, 0.0), NDir::new(0.0, 1.0, 0.0).unwrap());
        cone.mirror_ax1(&axis);
        assert_eq!(cone.location(), &NPnt::new(1.0, 0.0, 0.0));
        assert_eq!(cone.axis().direction(), &NDir::new(0.0, 0.0, -1.0).unwrap());
        assert_eq!(cone.ref_radius(), 5.0);
    }

    #[test]
    fn test_mirror_ax2() {
        let mut cone = cone(
            (1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0),
            PI / 4.0,
            5.0,
        );
        let plane = NAx2::new(
            NPnt::new(0.0, 0.0, 0.0),
            NDir::new(0.0, 1.0, 0.0).unwrap(),
            NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        cone.mirror_ax2(&plane);
        assert_eq!(cone.location(), &NPnt::new(1.0, 0.0, 0.0));
        assert_eq!(cone.axis().direction(), &NDir::new(0.0, 0.0, -1.0).unwrap());
        assert_eq!(cone.ref_radius(), 5.0);
    }

    #[test]
    fn test_rotate() {
        let mut cone = cone(
            (1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0),
            PI / 4.0,
            5.0,
        );
        let axis = NAx1::new(NPnt::new(0.0, 0.0, 0.0), NDir::new(0.0, 0.0, 1.0).unwrap());
        cone.rotate(&axis, PI / 2.0);
        assert!((cone.location().x() + 1.0).abs() < 1e-5);
        assert!((cone.location().y() - 0.0).abs() < 1e-5);
        assert_eq!(cone.axis().direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
        assert_eq!(cone.ref_radius(), 5.0);
    }

    #[test]
    fn test_scale() {
        let mut cone = cone(
            (1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0),
            PI / 4.0,
            5.0,
        );
        let point = NPnt::new(0.0, 0.0, 0.0);
        cone.scale(&point, 2.0);
        assert_eq!(cone.location(), &NPnt::new(2.0, 0.0, 0.0));
        assert_eq!(cone.ref_radius(), 10.0);
        assert_eq!(cone.semi_angle(), PI / 4.0);

        cone.scale(&point, -2.0);
        assert_eq!(cone.location(), &NPnt::new(-4.0, 0.0, 0.0));
        assert_eq!(cone.ref_radius(), 20.0);
        assert_eq!(
            cone.x_axis().direction(),
            &NDir::new(-1.0, 0.0, 0.0).unwrap()
        );
    }

    #[test]
    fn test_transform() {
        let mut cone = cone(
            (1.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0),
            PI / 4.0,
            5.0,
        );
        let trsf = NTrsf::new_scale(&NPnt::new(0.0, 0.0, 0.0), 2.0).unwrap();
        cone.transform(&trsf);
        assert_eq!(cone.location(), &NPnt::new(2.0, 0.0, 0.0));
        assert_eq!(cone.ref_radius(), 10.0);
        assert_eq!(cone.semi_angle(), PI / 4.0);
    }

    #[test]
    fn test_translate_vec() {
        let mut cone = cone(
            (1.0, 2.0, 3.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0),
            PI / 4.0,
            5.0,
        );
        let vec = NVec::new(1.0, 1.0, 1.0);
        cone.translate_vec(&vec);
        assert_eq!(cone.location(), &NPnt::new(2.0, 3.0, 4.0));
        assert_eq!(cone.ref_radius(), 5.0);
        assert_eq!(cone.semi_angle(), PI / 4.0);
    }

    #[test]
    fn test_translate_pnts() {
        let mut cone = cone(
            (1.0, 2.0, 3.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 0.0),
            PI / 4.0,
            5.0,
        );
        let p1 = NPnt::new(0.0, 0.0, 0.0);
        let p2 = NPnt::new(1.0, 1.0, 1.0);
        cone.translate_pnts(&p1, &p2);
        assert_eq!(cone.location(), &NPnt::new(2.0, 3.0, 4.0));
        assert_eq!(cone.ref_radius(), 5.0);
        assert_eq!(cone.semi_angle(), PI / 4.0);
    }
}
