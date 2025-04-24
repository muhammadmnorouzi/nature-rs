use std::f64::consts::PI;

use crate::{
    gp::{NAx1, NAx2, NPnt, NTrsf, NVec},
    nature_errors::NErrors,
};

// Trait to define the behavior of a circle in 3D space
pub trait Circ {
    fn new(position: NAx2, radius: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn set_axis(&mut self, axis: &NAx1) -> Result<(), NErrors>;
    fn set_location(&mut self, location: NPnt);
    fn set_position(&mut self, position: NAx2);
    fn set_radius(&mut self, radius: f64) -> Result<(), NErrors>;
    fn area(&self) -> f64;
    fn length(&self) -> f64;
    fn axis(&self) -> &NAx1;
    fn location(&self) -> &NPnt;
    fn position(&self) -> &NAx2;
    fn radius(&self) -> f64;
    fn x_axis(&self) -> NAx1;
    fn y_axis(&self) -> NAx1;
    fn distance(&self, point: &NPnt) -> f64;
    fn square_distance(&self, point: &NPnt) -> f64;
    fn contains(&self, point: &NPnt, linear_tolerance: f64) -> bool;
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

// Struct representing a circle in 3D space
#[derive(Clone, PartialEq, Debug)]
pub struct NCirc {
    pos: NAx2,
    radius: f64,
}

impl Circ for NCirc {
    fn new(position: NAx2, radius: f64) -> Result<Self, NErrors> {
        if radius < 0.0 {
            return Err(NErrors::NegativeRadius);
        }
        Ok(NCirc {
            pos: position,
            radius,
        })
    }

    fn set_axis(&mut self, axis: &NAx1) -> Result<(), NErrors> {
        self.pos.set_axis(axis)
    }

    fn set_location(&mut self, location: NPnt) {
        self.pos.set_location(location);
    }

    fn set_position(&mut self, position: NAx2) {
        self.pos = position;
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

    fn axis(&self) -> &NAx1 {
        self.pos.axis()
    }

    fn location(&self) -> &NPnt {
        self.pos.location()
    }

    fn position(&self) -> &NAx2 {
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

    fn distance(&self, point: &NPnt) -> f64 {
        self.square_distance(point).sqrt()
    }

    fn square_distance(&self, point: &NPnt) -> f64 {
        let v = NVec::new_from_points(self.location(), point);
        let x = v.dot(self.pos.x_direction());
        let y = v.dot(self.pos.y_direction());
        let z = v.dot(self.pos.direction());
        let t = (x * x + y * y).sqrt() - self.radius;
        t * t + z * z
    }

    fn contains(&self, point: &NPnt, linear_tolerance: f64) -> bool {
        self.distance(point) <= linear_tolerance
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
        self.radius *= factor.abs();
        self.pos.scale(point, factor);
    }

    fn scaled(&self, point: &NPnt, factor: f64) -> Self {
        let mut result = self.clone();
        result.scale(point, factor);
        result
    }

    fn transform(&mut self, transformation: &NTrsf) {
        self.radius *= transformation.scale_factor().abs();
        self.pos.transform(transformation);
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

    fn circ(
        pos: (f64, f64, f64),
        dir: (f64, f64, f64),
        x_dir: (f64, f64, f64),
        radius: f64,
    ) -> NCirc {
        NCirc::new(
            NAx2::new(
                NPnt::new(pos.0, pos.1, pos.2),
                NDir::new(dir.0, dir.1, dir.2).expect("Invalid direction"),
                NDir::new(x_dir.0, x_dir.1, x_dir.2).expect("Invalid X direction"),
            )
            .expect("Invalid NAx2"),
            radius,
        )
        .expect("Invalid NCirc")
    }

    #[test]
    fn test_new() {
        let pos = NAx2::new(
            NPnt::new(1.0, 2.0, 3.0),
            NDir::new(0.0, 0.0, 1.0).unwrap(),
            NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let circ = NCirc::new(pos.clone(), 5.0).unwrap();
        assert_eq!(circ.position(), &pos);
        assert_eq!(circ.radius(), 5.0);

        let result = NCirc::new(pos, -1.0);
        assert!(matches!(result, Err(NErrors::NegativeRadius)));
    }

    #[test]
    fn test_setters() {
        let mut circ = circ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let new_axis = NAx1::new(NPnt::new(0.0, 0.0, 0.0), NDir::new(0.0, 1.0, 0.0).unwrap());
        circ.set_axis(&new_axis).unwrap();
        assert_eq!(circ.axis().direction(), &NDir::new(0.0, 1.0, 0.0).unwrap());

        circ.set_location(NPnt::new(1.0, 2.0, 3.0));
        assert_eq!(circ.location(), &NPnt::new(1.0, 2.0, 3.0));

        let new_pos = NAx2::new(
            NPnt::new(4.0, 5.0, 6.0),
            NDir::new(0.0, 0.0, 1.0).unwrap(),
            NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        circ.set_position(new_pos.clone());
        assert_eq!(circ.position(), &new_pos);

        circ.set_radius(10.0).unwrap();
        assert_eq!(circ.radius(), 10.0);
        let result = circ.set_radius(-1.0);
        assert!(matches!(result, Err(NErrors::NegativeRadius)));
    }

    #[test]
    fn test_area_length() {
        let circ = circ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        assert!((circ.area() - 25.0 * PI).abs() < 1e-10);
        assert!((circ.length() - 10.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn test_axes() {
        let circ = circ((1.0, 2.0, 3.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        assert_eq!(circ.axis().location(), &NPnt::new(1.0, 2.0, 3.0));
        assert_eq!(circ.axis().direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
        assert_eq!(
            circ.x_axis().direction(),
            &NDir::new(1.0, 0.0, 0.0).unwrap()
        );
        assert_eq!(
            circ.y_axis().direction(),
            &NDir::new(0.0, 1.0, 0.0).unwrap()
        );
    }

    #[test]
    fn test_distance() {
        let circ = circ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let p = NPnt::new(3.0, 4.0, 0.0); // Point on circle (5^2 = 3^2 + 4^2)
        assert!(circ.distance(&p).abs() < 1e-10);
        assert!(circ.contains(&p, 1e-5));

        let p2 = NPnt::new(0.0, 0.0, 2.0); // Point above center
        assert!((circ.distance(&p2) - 2.0).abs() < 1e-10);
        assert!((circ.square_distance(&p2) - 4.0).abs() < 1e-10);
        assert!(!circ.contains(&p2, 1e-5));
    }

    #[test]
    fn test_mirror_pnt() {
        let mut circ = circ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let point = NPnt::new(0.0, 0.0, 0.0);
        circ.mirror_pnt(&point);
        assert_eq!(circ.location(), &NPnt::new(-1.0, 0.0, 0.0));
        assert_eq!(circ.radius(), 5.0);
    }

    #[test]
    fn test_mirror_ax1() {
        let mut circ = circ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let axis = NAx1::new(NPnt::new(0.0, 0.0, 0.0), NDir::new(0.0, 1.0, 0.0).unwrap());
        circ.mirror_ax1(&axis);
        assert_eq!(circ.location(), &NPnt::new(1.0, 0.0, 0.0));
        assert_eq!(circ.axis().direction(), &NDir::new(0.0, 0.0, -1.0).unwrap());
        assert_eq!(circ.radius(), 5.0);
    }

    #[test]
    fn test_mirror_ax2() {
        let mut circ = circ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let plane = NAx2::new(
            NPnt::new(0.0, 0.0, 0.0),
            NDir::new(0.0, 1.0, 0.0).unwrap(),
            NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        circ.mirror_ax2(&plane);
        assert_eq!(circ.location(), &NPnt::new(1.0, 0.0, 0.0));
        assert_eq!(circ.axis().direction(), &NDir::new(0.0, 0.0, -1.0).unwrap());
        assert_eq!(circ.radius(), 5.0);
    }

    #[test]
    fn test_rotate() {
        let mut circ = circ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let axis = NAx1::new(NPnt::new(0.0, 0.0, 0.0), NDir::new(0.0, 0.0, 1.0).unwrap());
        circ.rotate(&axis, PI / 2.0);
        assert!((circ.location().x() + 1.0).abs() < 1e-5);
        assert!((circ.location().y() - 0.0).abs() < 1e-5);
        assert_eq!(circ.axis().direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
        assert_eq!(circ.radius(), 5.0);
    }

    #[test]
    fn test_scale() {
        let mut circ = circ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let point = NPnt::new(0.0, 0.0, 0.0);
        circ.scale(&point, 2.0);
        assert_eq!(circ.location(), &NPnt::new(2.0, 0.0, 0.0));
        assert_eq!(circ.radius(), 10.0);

        circ.scale(&point, -2.0);
        assert_eq!(circ.location(), &NPnt::new(-4.0, 0.0, 0.0));
        assert_eq!(circ.radius(), 20.0);
        assert_eq!(
            circ.x_axis().direction(),
            &NDir::new(-1.0, 0.0, 0.0).unwrap()
        );
    }

    #[test]
    fn test_transform() {
        let mut circ = circ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let trsf = NTrsf::new_scale(&NPnt::new(0.0, 0.0, 0.0), 2.0).unwrap();
        circ.transform(&trsf);
        assert_eq!(circ.location(), &NPnt::new(2.0, 0.0, 0.0));
        assert_eq!(circ.radius(), 10.0);
    }

    #[test]
    fn test_translate_vec() {
        let mut circ = circ((1.0, 2.0, 3.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let vec = NVec::new(1.0, 1.0, 1.0);
        circ.translate_vec(&vec);
        assert_eq!(circ.location(), &NPnt::new(2.0, 3.0, 4.0));
        assert_eq!(circ.radius(), 5.0);
    }

    #[test]
    fn test_translate_pnts() {
        let mut circ = circ((1.0, 2.0, 3.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 5.0);
        let p1 = NPnt::new(0.0, 0.0, 0.0);
        let p2 = NPnt::new(1.0, 1.0, 1.0);
        circ.translate_pnts(&p1, &p2);
        assert_eq!(circ.location(), &NPnt::new(2.0, 3.0, 4.0));
        assert_eq!(circ.radius(), 5.0);
    }
}
