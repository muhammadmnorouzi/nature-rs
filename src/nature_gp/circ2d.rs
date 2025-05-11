use std::f64::consts::PI;

use crate::{
    gp::{NAx2d, NAx22d, NPoint2d, NTrsf2d, NVec2d},
    nature_errors::NErrors,
};

// Trait to define the behavior of a circle in 2D space
pub trait Circ2d {
    fn new(axis: NAx22d, radius: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn new_with_sense(x_axis: NAx2d, radius: f64, is_right_handed: bool) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn set_location(&mut self, location: NPoint2d);
    fn set_x_axis(&mut self, axis: NAx2d);
    fn set_y_axis(&mut self, axis: NAx2d);
    fn set_axis(&mut self, axis: NAx22d);
    fn set_radius(&mut self, radius: f64) -> Result<(), NErrors>;
    fn area(&self) -> f64;
    fn length(&self) -> f64;
    fn coefficients(&self) -> (f64, f64, f64, f64, f64, f64);
    fn contains(&self, point: &NPoint2d, linear_tolerance: f64) -> bool;
    fn distance(&self, point: &NPoint2d) -> f64;
    fn square_distance(&self, point: &NPoint2d) -> f64;
    fn location(&self) -> &NPoint2d;
    fn radius(&self) -> f64;
    fn axis(&self) -> &NAx22d;
    fn position(&self) -> &NAx22d;
    fn x_axis(&self) -> NAx2d;
    fn y_axis(&self) -> NAx2d;
    fn reverse(&mut self);
    fn reversed(&self) -> Self;
    fn is_direct(&self) -> bool;
    fn mirror_pnt(&mut self, point: &NPoint2d);
    fn mirrored_pnt(&self, point: &NPoint2d) -> Self;
    fn mirror_ax2d(&mut self, axis: &NAx2d);
    fn mirrored_ax2d(&self, axis: &NAx2d) -> Self;
    fn rotate(&mut self, point: &NPoint2d, angle: f64);
    fn rotated(&self, point: &NPoint2d, angle: f64) -> Self;
    fn scale(&mut self, point: &NPoint2d, factor: f64);
    fn scaled(&self, point: &NPoint2d, factor: f64) -> Self;
    fn transform(&mut self, transformation: &NTrsf2d);
    fn transformed(&self, transformation: &NTrsf2d) -> Self;
    fn translate_vec(&mut self, vector: &NVec2d);
    fn translated_vec(&self, vector: &NVec2d) -> Self;
    fn translate_point3d(&mut self, from: &NPoint2d, to: &NPoint2d);
    fn translated_point3d(&self, from: &NPoint2d, to: &NPoint2d) -> Self;
}

// Struct representing a circle in 2D space
#[derive(Clone, PartialEq, Debug)]
pub struct NCirc2d {
    pos: NAx22d,
    radius: f64,
}

impl Circ2d for NCirc2d {
    fn new(axis: NAx22d, radius: f64) -> Result<Self, NErrors> {
        if radius < 0.0 {
            return Err(NErrors::NegativeRadius);
        }
        Ok(NCirc2d { pos: axis, radius })
    }

    fn new_with_sense(x_axis: NAx2d, radius: f64, is_right_handed: bool) -> Result<Self, NErrors> {
        if radius < 0.0 {
            return Err(NErrors::NegativeRadius);
        }
        let pos = NAx22d::from_ax2d_with_sense(&x_axis, is_right_handed);
        Ok(NCirc2d { pos, radius })
    }

    fn set_location(&mut self, location: NPoint2d) {
        self.pos.set_location(location);
    }

    fn set_x_axis(&mut self, axis: NAx2d) {
        self.pos.set_x_axis(&axis);
    }

    fn set_y_axis(&mut self, axis: NAx2d) {
        self.pos.set_y_axis(&axis);
    }

    fn set_axis(&mut self, axis: NAx22d) {
        self.pos = axis;
    }

    fn set_radius(&mut self, radius: f64) -> Result<(), NErrors> {
        if radius < 0.0 {
            return Err(NErrors::NegativeRadius);
        }
        self.radius = radius;
        Ok(())
    }

    fn area(&self) -> f64 {
        PI * self.radius * self.radius
    }

    fn length(&self) -> f64 {
        2.0 * PI * self.radius
    }

    fn coefficients(&self) -> (f64, f64, f64, f64, f64, f64) {
        let xc = self.pos.location().x();
        let yc = self.pos.location().y();
        (
            1.0,                                           // A
            1.0,                                           // B
            0.0,                                           // C
            -xc,                                           // D
            -yc,                                           // E
            xc * xc + yc * yc - self.radius * self.radius, // F
        )
    }

    fn contains(&self, point: &NPoint2d, linear_tolerance: f64) -> bool {
        self.distance(point) <= linear_tolerance
    }

    fn distance(&self, point: &NPoint2d) -> f64 {
        let v = NVec2d::new_from_points(self.location(), point);
        let d = self.radius - v.magnitude();
        d.abs()
    }

    fn square_distance(&self, point: &NPoint2d) -> f64 {
        let v = NVec2d::new_from_points(self.location(), point);
        let d = self.radius - v.magnitude();
        d * d
    }

    fn location(&self) -> &NPoint2d {
        self.pos.location()
    }

    fn radius(&self) -> f64 {
        self.radius
    }

    fn axis(&self) -> &NAx22d {
        &self.pos
    }

    fn position(&self) -> &NAx22d {
        &self.pos
    }

    fn x_axis(&self) -> NAx2d {
        self.pos.x_axis()
    }

    fn y_axis(&self) -> NAx2d {
        self.pos.y_axis()
    }

    fn reverse(&mut self) {
        let mut y_dir = self.pos.y_direction().clone();
        y_dir.reverse();
        self.pos.set_axis(
            &NAx22d::new(
                self.pos.location().clone(),
                self.pos.x_direction().clone(),
                y_dir,
            )
            .expect("Invalid NAx22d"),
        );
    }

    fn reversed(&self) -> Self {
        let mut result = self.clone();
        result.reverse();
        result
    }

    fn is_direct(&self) -> bool {
        self.pos.x_direction().crossed(self.pos.y_direction()) >= 0.0
    }

    fn mirror_pnt(&mut self, point: &NPoint2d) {
        self.pos.mirror_pnt(point);
    }

    fn mirrored_pnt(&self, point: &NPoint2d) -> Self {
        let mut result = self.clone();
        result.mirror_pnt(point);
        result
    }

    fn mirror_ax2d(&mut self, axis: &NAx2d) {
        self.pos.mirror_ax2d(axis);
    }

    fn mirrored_ax2d(&self, axis: &NAx2d) -> Self {
        let mut result = self.clone();
        result.mirror_ax2d(axis);
        result
    }

    fn rotate(&mut self, point: &NPoint2d, angle: f64) {
        self.pos.rotate(point, angle);
    }

    fn rotated(&self, point: &NPoint2d, angle: f64) -> Self {
        let mut result = self.clone();
        result.rotate(point, angle);
        result
    }

    fn scale(&mut self, point: &NPoint2d, factor: f64) {
        self.radius *= factor.abs();
        self.pos.scale(point, factor);
    }

    fn scaled(&self, point: &NPoint2d, factor: f64) -> Self {
        let mut result = self.clone();
        result.scale(point, factor);
        result
    }

    fn transform(&mut self, transformation: &NTrsf2d) {
        self.radius *= transformation.scale_factor().abs();
        self.pos.transform(transformation);
    }

    fn transformed(&self, transformation: &NTrsf2d) -> Self {
        let mut result = self.clone();
        result.transform(transformation);
        result
    }

    fn translate_vec(&mut self, vector: &NVec2d) {
        self.pos.translate_vec(vector);
    }

    fn translated_vec(&self, vector: &NVec2d) -> Self {
        let mut result = self.clone();
        result.translate_vec(vector);
        result
    }

    fn translate_point3d(&mut self, from: &NPoint2d, to: &NPoint2d) {
        self.pos.translate_point3d(from, to);
    }

    fn translated_point3d(&self, from: &NPoint2d, to: &NPoint2d) -> Self {
        let mut result = self.clone();
        result.translate_point3d(from, to);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn circ2d(pos: (f64, f64), x_dir: (f64, f64), y_dir: (f64, f64), radius: f64) -> NCirc2d {
        NCirc2d::new(
            NAx22d::new(
                NPoint2d::new(pos.0, pos.1),
                NDir2d::new(x_dir.0, x_dir.1).expect("Invalid X direction"),
                NDir2d::new(y_dir.0, y_dir.1).expect("Invalid Y direction"),
            )
            .expect("Invalid NAx22d"),
            radius,
        )
        .expect("Invalid NCirc2d")
    }

    #[test]
    fn test_new() {
        let pos = NAx22d::new(
            NPoint2d::new(1.0, 2.0),
            NDir2d::new(1.0, 0.0).unwrap(),
            NDir2d::new(0.0, 1.0).unwrap(),
        )
        .unwrap();
        let circ = NCirc2d::new(pos.clone(), 5.0).unwrap();
        assert_eq!(circ.position(), &pos);
        assert_eq!(circ.radius(), 5.0);

        let result = NCirc2d::new(pos, -1.0);
        assert!(matches!(result, Err(NErrors::NegativeRadius)));
    }

    #[test]
    fn test_new_with_sense() {
        let x_axis = NAx2d::new(NPoint2d::new(1.0, 2.0), NDir2d::new(1.0, 0.0).unwrap());
        let circ = NCirc2d::new_with_sense(x_axis.clone(), 5.0, true).unwrap();
        assert_eq!(circ.location(), &NPoint2d::new(1.0, 2.0));
        assert_eq!(circ.x_axis().direction(), &NDir2d::new(1.0, 0.0).unwrap());
        assert_eq!(circ.y_axis().direction(), &NDir2d::new(0.0, 1.0).unwrap());
        assert_eq!(circ.radius(), 5.0);

        let circ_left = NCirc2d::new_with_sense(x_axis, 5.0, false).unwrap();
        assert_eq!(
            circ_left.y_axis().direction(),
            &NDir2d::new(0.0, -1.0).unwrap()
        );
    }

    #[test]
    fn test_setters() {
        let mut circ = circ2d((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), 5.0);
        circ.set_location(NPoint2d::new(1.0, 2.0));
        assert_eq!(circ.location(), &NPoint2d::new(1.0, 2.0));

        let new_x_axis = NAx2d::new(NPoint2d::new(3.0, 4.0), NDir2d::new(0.0, 1.0).unwrap());
        circ.set_x_axis(new_x_axis.clone());
        assert_eq!(circ.x_axis(), new_x_axis);

        let new_y_axis = NAx2d::new(NPoint2d::new(3.0, 4.0), NDir2d::new(1.0, 0.0).unwrap());
        circ.set_y_axis(new_y_axis.clone());
        assert_eq!(circ.y_axis(), new_y_axis);

        let new_axis = NAx22d::new(
            NPoint2d::new(5.0, 6.0),
            NDir2d::new(1.0, 0.0).unwrap(),
            NDir2d::new(0.0, 1.0).unwrap(),
        )
        .unwrap();
        circ.set_axis(new_axis.clone());
        assert_eq!(circ.axis(), &new_axis);

        circ.set_radius(10.0).unwrap();
        assert_eq!(circ.radius(), 10.0);
        let result = circ.set_radius(-1.0);
        assert!(matches!(result, Err(NErrors::NegativeRadius)));
    }

    #[test]
    fn test_area_length() {
        let circ = circ2d((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), 5.0);
        assert!((circ.area() - 25.0 * PI).abs() < 1e-10);
        assert!((circ.length() - 10.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn test_coefficients() {
        let circ = circ2d((1.0, 2.0), (1.0, 0.0), (0.0, 1.0), 5.0);
        let (a, b, c, d, e, f) = circ.coefficients();
        assert_eq!(a, 1.0);
        assert_eq!(b, 1.0);
        assert_eq!(c, 0.0);
        assert_eq!(d, -1.0);
        assert_eq!(e, -2.0);
        assert!((f - (1.0 + 4.0 - 25.0)).abs() < 1e-10);
    }

    #[test]
    fn test_distance() {
        let circ = circ2d((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), 5.0);
        let p = NPoint2d::new(3.0, 4.0); // Point on circle (5^2 = 3^2 + 4^2)
        assert!(circ.distance(&p).abs() < 1e-10);
        assert!(circ.square_distance(&p).abs() < 1e-10);
        assert!(circ.contains(&p, 1e-5));

        let p2 = NPoint2d::new(10.0, 0.0); // Point outside circle
        assert!((circ.distance(&p2) - 5.0).abs() < 1e-10);
        assert!((circ.square_distance(&p2) - 25.0).abs() < 1e-10);
        assert!(!circ.contains(&p2, 1e-5));
    }

    #[test]
    fn test_reverse() {
        let mut circ = circ2d((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), 5.0);
        assert!(circ.is_direct());
        circ.reverse();
        assert!(!circ.is_direct());
        assert_eq!(circ.y_axis().direction(), &NDir2d::new(0.0, -1.0).unwrap());
        assert_eq!(circ.radius(), 5.0);
    }

    #[test]
    fn test_mirror_pnt() {
        let mut circ = circ2d((1.0, 0.0), (1.0, 0.0), (0.0, 1.0), 5.0);
        let point = NPoint2d::new(0.0, 0.0);
        circ.mirror_pnt(&point);
        assert_eq!(circ.location(), &NPoint2d::new(-1.0, 0.0));
        assert_eq!(circ.radius(), 5.0);
    }

    #[test]
    fn test_mirror_ax2d() {
        let mut circ = circ2d((1.0, 0.0), (1.0, 0.0), (0.0, 1.0), 5.0);
        let axis = NAx2d::new(NPoint2d::new(0.0, 0.0), NDir2d::new(0.0, 1.0).unwrap());
        circ.mirror_ax2d(&axis);
        assert_eq!(circ.location(), &NPoint2d::new(1.0, 0.0));
        assert_eq!(circ.x_axis().direction(), &NDir2d::new(-1.0, 0.0).unwrap());
        assert_eq!(circ.radius(), 5.0);
    }

    #[test]
    fn test_rotate() {
        let mut circ = circ2d((1.0, 0.0), (1.0, 0.0), (0.0, 1.0), 5.0);
        let point = NPoint2d::new(0.0, 0.0);
        circ.rotate(&point, PI / 2.0);
        assert!((circ.location().x() - 0.0).abs() < 1e-5);
        assert!((circ.location().y() - 1.0).abs() < 1e-5);
        assert_eq!(circ.radius(), 5.0);
    }

    #[test]
    fn test_scale() {
        let mut circ = circ2d((1.0, 0.0), (1.0, 0.0), (0.0, 1.0), 5.0);
        let point = NPoint2d::new(0.0, 0.0);
        circ.scale(&point, 2.0);
        assert_eq!(circ.location(), &NPoint2d::new(2.0, 0.0));
        assert_eq!(circ.radius(), 10.0);

        circ.scale(&point, -2.0);
        assert_eq!(circ.location(), &NPoint2d::new(-4.0, 0.0));
        assert_eq!(circ.radius(), 20.0);
        assert_eq!(circ.x_axis().direction(), &NDir2d::new(-1.0, 0.0).unwrap());
    }

    #[test]
    fn test_transform() {
        let mut circ = circ2d((1.0, 0.0), (1.0, 0.0), (0.0, 1.0), 5.0);
        let trsf = NTrsf2d::new_scale(&NPoint2d::new(0.0, 0.0), 2.0).unwrap();
        circ.transform(&trsf);
        assert_eq!(circ.location(), &NPoint2d::new(2.0, 0.0));
        assert_eq!(circ.radius(), 10.0);
    }

    #[test]
    fn test_translate_vec() {
        let mut circ = circ2d((1.0, 2.0), (1.0, 0.0), (0.0, 1.0), 5.0);
        let vec = NVec2d::new(1.0, 1.0);
        circ.translate_vec(&vec);
        assert_eq!(circ.location(), &NPoint2d::new(2.0, 3.0));
        assert_eq!(circ.radius(), 5.0);
    }

    #[test]
    fn test_translate_point3d() {
        let mut circ = circ2d((1.0, 2.0), (1.0, 0.0), (0.0, 1.0), 5.0);
        let p1 = NPoint2d::new(0.0, 0.0);
        let p2 = NPoint2d::new(1.0, 1.0);
        circ.translate_point3d(&p1, &p2);
        assert_eq!(circ.location(), &NPoint2d::new(2.0, 3.0));
        assert_eq!(circ.radius(), 5.0);
    }
}
