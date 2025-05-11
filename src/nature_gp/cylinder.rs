use serde::{Deserialize, Serialize};

use crate::{
    gp::{NAx1, NAx2, NAx3, NPoint3d, NTrsf, NVec},
    nature_errors::NErrors,
};

// Trait to define the behavior of a cylindrical surface in 3D space
pub trait Cylinder {
    fn new(position: NAx3, radius: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn set_axis(&mut self, axis: &NAx1) -> Result<(), NErrors>;
    fn set_location(&mut self, location: NPoint3d);
    fn set_position(&mut self, position: NAx3);
    fn set_radius(&mut self, radius: f64) -> Result<(), NErrors>;
    fn u_reverse(&mut self);
    fn v_reverse(&mut self);
    fn direct(&self) -> bool;
    fn axis(&self) -> &NAx1;
    fn location(&self) -> &NPoint3d;
    fn position(&self) -> &NAx3;
    fn radius(&self) -> f64;
    fn x_axis(&self) -> NAx1;
    fn y_axis(&self) -> NAx1;
    fn coefficients(&self) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64);
    fn mirror_point3d(&mut self, point: &NPoint3d);
    fn mirrored_point3d(&self, point: &NPoint3d) -> Self;
    fn mirror_ax1(&mut self, axis: &NAx1);
    fn mirrored_ax1(&self, axis: &NAx1) -> Self;
    fn mirror_ax2(&mut self, plane: &NAx2);
    fn mirrored_ax2(&self, plane: &NAx2) -> Self;
    fn rotate(&mut self, axis: &NAx1, angle: f64);
    fn rotated(&self, axis: &NAx1, angle: f64) -> Self;
    fn scale(&mut self, point: &NPoint3d, factor: f64);
    fn scaled(&self, point: &NPoint3d, factor: f64) -> Self;
    fn transform(&mut self, transformation: &NTrsf);
    fn transformed(&self, transformation: &NTrsf) -> Self;
    fn translate_vec(&mut self, vector: &NVec);
    fn translated_vec(&self, vector: &NVec) -> Self;
    fn translate_point3d(&mut self, from: &NPoint3d, to: &NPoint3d);
    fn translated_point3d(&self, from: &NPoint3d, to: &NPoint3d) -> Self;
}

// Struct representing an infinite cylindrical surface in 3D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NCylinder {
    pos: NAx3,
    radius: f64,
}

impl Cylinder for NCylinder {
    fn new(position: NAx3, radius: f64) -> Result<Self, NErrors> {
        if radius < 0.0 {
            return Err(NErrors::NegativeRadius);
        }
        Ok(NCylinder {
            pos: position,
            radius,
        })
    }

    fn set_axis(&mut self, axis: &NAx1) -> Result<(), NErrors> {
        self.pos.set_axis(axis)
    }

    fn set_location(&mut self, location: NPoint3d) {
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

    fn u_reverse(&mut self) {
        self.pos.y_reverse();
    }

    fn v_reverse(&mut self) {
        self.pos.z_reverse();
    }

    fn direct(&self) -> bool {
        self.pos.direct()
    }

    fn axis(&self) -> &NAx1 {
        self.pos.axis()
    }

    fn location(&self) -> &NPoint3d {
        self.pos.location()
    }

    fn position(&self) -> &NAx3 {
        &self.pos
    }

    fn radius(&self) -> f64 {
        self.radius
    }

    fn x_axis(&self) -> NAx1 {
        NAx1::new(self.pos.location().clone(), self.pos.x_direction().clone())
    }

    fn y_axis(&self) -> NAx1 {
        NAx1::new(self.pos.location().clone(), self.pos.y_direction().clone())
    }

    fn coefficients(&self) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
        let trsf = NTrsf::new_transformation(&self.pos).expect("Invalid transformation");
        let t11 = trsf.value(1, 1);
        let t12 = trsf.value(1, 2);
        let t13 = trsf.value(1, 3);
        let t14 = trsf.value(1, 4);
        let t21 = trsf.value(2, 1);
        let t22 = trsf.value(2, 2);
        let t23 = trsf.value(2, 3);
        let t24 = trsf.value(2, 4);
        let a1 = t11 * t11 + t21 * t21;
        let a2 = t12 * t12 + t22 * t22;
        let a3 = t13 * t13 + t23 * t23;
        let b1 = t11 * t12 + t21 * t22;
        let b2 = t11 * t13 + t21 * t23;
        let b3 = t12 * t13 + t22 * t23;
        let c1 = t11 * t14 + t21 * t24;
        let c2 = t12 * t14 + t22 * t24;
        let c3 = t13 * t14 + t23 * t24;
        let d = t14 * t14 + t24 * t24 - self.radius * self.radius;
        (a1, a2, a3, b1, b2, b3, c1, c2, c3, d)
    }

    fn mirror_point3d(&mut self, point: &NPoint3d) {
        self.pos.mirror_point3d(point);
    }

    fn mirrored_point3d(&self, point: &NPoint3d) -> Self {
        let mut result = self.clone();
        result.mirror_point3d(point);
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

    fn scale(&mut self, point: &NPoint3d, factor: f64) {
        self.pos.scale(point, factor);
        self.radius *= factor.abs();
    }

    fn scaled(&self, point: &NPoint3d, factor: f64) -> Self {
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

    fn translate_point3d(&mut self, from: &NPoint3d, to: &NPoint3d) {
        self.pos.translate_point3d(from, to);
    }

    fn translated_point3d(&self, from: &NPoint3d, to: &NPoint3d) -> Self {
        let mut result = self.clone();
        result.translate_point3d(from, to);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn cylinder(
        pos: (f64, f64, f64),
        dir: (f64, f64, f64),
        x_dir: (f64, f64, f64),
        radius: f64,
    ) -> NCylinder {
        NCylinder::new(
            NAx3::new(
                NPoint3d::new(pos.0, pos.1, pos.2),
                NDir::new(dir.0, dir.1, dir.2).expect("Invalid direction"),
                NDir::new(x_dir.0, x_dir.1, x_dir.2).expect("Invalid X direction"),
            )
            .expect("Invalid NAx3"),
            radius,
        )
        .expect("Invalid NCylinder")
    }

    #[test]
    fn test_new() {
        let pos = NAx3::new(
            NPoint3d::new(1.0, 2.0, 3.0),
            NDir::new(0.0, 0.0, 1.0).unwrap(),
            NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let cyl = NCylinder::new(pos.clone(), 5.0).unwrap();
        assert_eq!(cyl.position(), &pos);
        assert_eq!(cyl.radius(), 5.0);

        assert!(matches!(
            NCylinder::new(pos, -1.0),
            Err(NErrors::NegativeRadius)
        ));
    }

    #[test]
    fn test_setters() {
        let mut cyl = cylinder((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let new_axis = NAx1::new(
            NPoint3d::new(0.0, 0.0, 0.0),
            NDir::new(0.0, 1.0, 0.0).unwrap(),
        );
        cyl.set_axis(&new_axis).unwrap();
        assert_eq!(cyl.axis().direction(), &NDir::new(0.0, 1.0, 0.0).unwrap());

        cyl.set_location(NPoint3d::new(1.0, 2.0, 3.0));
        assert_eq!(cyl.location(), &NPoint3d::new(1.0, 2.0, 3.0));

        let new_pos = NAx3::new(
            NPoint3d::new(4.0, 5.0, 6.0),
            NDir::new(0.0, 0.0, 1.0).unwrap(),
            NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        cyl.set_position(new_pos.clone());
        assert_eq!(cyl.position(), &new_pos);

        cyl.set_radius(10.0).unwrap();
        assert_eq!(cyl.radius(), 10.0);
        assert!(matches!(cyl.set_radius(-1.0), Err(NErrors::NegativeRadius)));
    }

    #[test]
    fn test_reverse() {
        let mut cyl = cylinder((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        assert!(cyl.direct());
        cyl.u_reverse();
        assert_eq!(
            cyl.y_axis().direction(),
            &NDir::new(0.0, -1.0, 0.0).unwrap()
        );

        cyl.v_reverse();
        assert_eq!(cyl.axis().direction(), &NDir::new(0.0, 0.0, -1.0).unwrap());
        assert_eq!(cyl.radius(), 5.0);
    }

    #[test]
    fn test_axes() {
        let cyl = cylinder((1.0, 2.0, 3.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        assert_eq!(cyl.axis().location(), &NPoint3d::new(1.0, 2.0, 3.0));
        assert_eq!(cyl.axis().direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
        assert_eq!(cyl.x_axis().direction(), &NDir::new(1.0, 0.0, 0.0).unwrap());
        assert_eq!(cyl.y_axis().direction(), &NDir::new(0.0, 1.0, 0.0).unwrap());
    }

    #[test]
    fn test_coefficients() {
        let cyl = cylinder((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let (a1, a2, a3, b1, b2, b3, c1, c2, c3, d) = cyl.coefficients();
        assert!((a1 - 1.0).abs() < 1e-10);
        assert!((a2 - 1.0).abs() < 1e-10);
        assert_eq!(a3, 0.0);
        assert_eq!(b1, 0.0);
        assert_eq!(b2, 0.0);
        assert_eq!(b3, 0.0);
        assert_eq!(c1, 0.0);
        assert_eq!(c2, 0.0);
        assert_eq!(c3, 0.0);
        assert!((d - (-25.0)).abs() < 1e-10);
    }

    #[test]
    fn test_mirror_point3d() {
        let mut cyl = cylinder((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let point = NPoint3d::new(0.0, 0.0, 0.0);
        cyl.mirror_point3d(&point);
        assert_eq!(cyl.location(), &NPoint3d::new(-1.0, 0.0, 0.0));
        assert_eq!(cyl.radius(), 5.0);
    }

    #[test]
    fn test_mirror_ax1() {
        let mut cyl = cylinder((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let axis = NAx1::new(
            NPoint3d::new(0.0, 0.0, 0.0),
            NDir::new(0.0, 1.0, 0.0).unwrap(),
        );
        cyl.mirror_ax1(&axis);
        assert_eq!(cyl.location(), &NPoint3d::new(1.0, 0.0, 0.0));
        assert_eq!(cyl.axis().direction(), &NDir::new(0.0, 0.0, -1.0).unwrap());
        assert_eq!(cyl.radius(), 5.0);
    }

    #[test]
    fn test_mirror_ax2() {
        let mut cyl = cylinder((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let plane = NAx2::new(
            NPoint3d::new(0.0, 0.0, 0.0),
            NDir::new(0.0, 1.0, 0.0).unwrap(),
            NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        cyl.mirror_ax2(&plane);
        assert_eq!(cyl.location(), &NPoint3d::new(1.0, 0.0, 0.0));
        assert_eq!(cyl.axis().direction(), &NDir::new(0.0, 0.0, -1.0).unwrap());
        assert_eq!(cyl.radius(), 5.0);
    }

    #[test]
    fn test_rotate() {
        let mut cyl = cylinder((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let axis = NAx1::new(
            NPoint3d::new(0.0, 0.0, 0.0),
            NDir::new(0.0, 0.0, 1.0).unwrap(),
        );
        cyl.rotate(&axis, PI / 2.0);
        assert!((cyl.location().x() + 1.0).abs() < 1e-5);
        assert!((cyl.location().y() - 0.0).abs() < 1e-5);
        assert_eq!(cyl.axis().direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
        assert_eq!(cyl.radius(), 5.0);
    }

    #[test]
    fn test_scale() {
        let mut cyl = cylinder((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let point = NPoint3d::new(0.0, 0.0, 0.0);
        cyl.scale(&point, 2.0);
        assert_eq!(cyl.location(), &NPoint3d::new(2.0, 0.0, 0.0));
        assert_eq!(cyl.radius(), 10.0);

        cyl.scale(&point, -2.0);
        assert_eq!(cyl.location(), &NPoint3d::new(-4.0, 0.0, 0.0));
        assert_eq!(cyl.radius(), 20.0);
        assert_eq!(
            cyl.x_axis().direction(),
            &NDir::new(-1.0, 0.0, 0.0).unwrap()
        );
    }

    #[test]
    fn test_transform() {
        let mut cyl = cylinder((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let trsf = NTrsf::new_scale(&NPoint3d::new(0.0, 0.0, 0.0), 2.0).unwrap();
        cyl.transform(&trsf);
        assert_eq!(cyl.location(), &NPoint3d::new(2.0, 0.0, 0.0));
        assert_eq!(cyl.radius(), 10.0);
    }

    #[test]
    fn test_translate_vec() {
        let mut cyl = cylinder((1.0, 2.0, 3.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let vec = NVec::new(1.0, 1.0, 1.0);
        cyl.translate_vec(&vec);
        assert_eq!(cyl.location(), &NPoint3d::new(2.0, 3.0, 4.0));
        assert_eq!(cyl.radius(), 5.0);
    }

    #[test]
    fn test_translate_point3d() {
        let mut cyl = cylinder((1.0, 2.0, 3.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let p1 = NPoint3d::new(0.0, 0.0, 0.0);
        let p2 = NPoint3d::new(1.0, 1.0, 1.0);
        cyl.translate_point3d(&p1, &p2);
        assert_eq!(cyl.location(), &NPoint3d::new(2.0, 3.0, 4.0));
        assert_eq!(cyl.radius(), 5.0);
    }
}
