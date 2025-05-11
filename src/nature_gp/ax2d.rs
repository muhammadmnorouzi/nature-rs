use serde::{Deserialize, Serialize};

use crate::{
    nature_errors::NErrors,
    nature_gp::{NDir2d, NGP, NPoint2d, NTrsf2d, NVec2d},
};

// Trait to define the behavior of an axis in 2D space
pub trait Ax2d {
    fn new(location: NPoint2d, direction: NDir2d) -> Self;
    fn x_axis() -> Self;
    fn set_location(&mut self, location: NPoint2d);
    fn set_direction(&mut self, direction: NDir2d);
    fn location(&self) -> &NPoint2d;
    fn direction(&self) -> &NDir2d;
    fn is_coaxial(&self, other: &Self, angular_tolerance: f64, linear_tolerance: f64) -> bool;
    fn is_normal(&self, other: &Self, angular_tolerance: f64) -> bool;
    fn is_opposite(&self, other: &Self, angular_tolerance: f64) -> bool;
    fn is_parallel(&self, other: &Self, angular_tolerance: f64) -> bool;
    fn angle(&self, other: &Self) -> f64;
    fn reverse(&mut self);
    fn reversed(&self) -> Self;
    fn mirror_pnt(&mut self, point: &NPoint2d);
    fn mirrored_pnt(&self, point: &NPoint2d) -> Self;
    fn mirror_ax2d(&mut self, axis: &Self);
    fn mirrored_ax2d(&self, axis: &Self) -> Self;
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

// Struct representing an axis in 2D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NAx2d {
    loc: NPoint2d,
    vdir: NDir2d,
}

impl Default for NAx2d {
    fn default() -> Self {
        Self::x_axis()
    }
}

impl Ax2d for NAx2d {
    fn new(location: NPoint2d, direction: NDir2d) -> Self {
        NAx2d {
            loc: location,
            vdir: direction,
        }
    }

    fn x_axis() -> Self {
        NAx2d {
            loc: NPoint2d::new(0.0, 0.0),
            vdir: NDir2d::new(1.0, 0.0).expect("Invalid direction"),
        }
    }

    fn set_location(&mut self, location: NPoint2d) {
        self.loc = location;
    }

    fn set_direction(&mut self, direction: NDir2d) {
        self.vdir = direction;
    }

    fn location(&self) -> &NPoint2d {
        &self.loc
    }

    fn direction(&self) -> &NDir2d {
        &self.vdir
    }

    fn is_coaxial(&self, other: &Self, angular_tolerance: f64, linear_tolerance: f64) -> bool {
        let mut xy1 = self.loc.xy();
        xy1.subtract(&other.loc.xy());
        let d1 = xy1.crossed(&other.vdir.xy()).abs();

        let mut xy2 = other.loc.xy();
        xy2.subtract(&self.loc.xy());
        let d2 = xy2.crossed(&self.vdir.xy()).abs();

        self.vdir.is_parallel(&other.vdir, angular_tolerance)
            && d1 <= linear_tolerance
            && d2 <= linear_tolerance
    }

    fn is_normal(&self, other: &Self, angular_tolerance: f64) -> bool {
        self.vdir.is_normal(&other.vdir, angular_tolerance)
    }

    fn is_opposite(&self, other: &Self, angular_tolerance: f64) -> bool {
        self.vdir.is_opposite(&other.vdir, angular_tolerance)
    }

    fn is_parallel(&self, other: &Self, angular_tolerance: f64) -> bool {
        self.vdir.is_parallel(&other.vdir, angular_tolerance)
    }

    fn angle(&self, other: &Self) -> f64 {
        self.vdir.angle(&other.vdir)
    }

    fn reverse(&mut self) {
        self.vdir.reverse();
    }

    fn reversed(&self) -> Self {
        NAx2d {
            loc: self.loc.clone(),
            vdir: self.vdir.reversed(),
        }
    }

    fn mirror_pnt(&mut self, point: &NPoint2d) {
        self.loc.mirror_pnt(point);
        self.vdir.reverse();
    }

    fn mirrored_pnt(&self, point: &NPoint2d) -> Self {
        let mut result = self.clone();
        result.mirror_pnt(point);
        result
    }

    fn mirror_ax2d(&mut self, axis: &Self) {
        self.loc.mirror_ax2d(axis);
        self.vdir.mirror(&axis.vdir);
    }

    fn mirrored_ax2d(&self, axis: &Self) -> Self {
        let mut result = self.clone();
        result.mirror_ax2d(axis);
        result
    }

    fn rotate(&mut self, point: &NPoint2d, angle: f64) {
        self.loc.rotate(point, angle);
        self.vdir.rotate(angle);
    }

    fn rotated(&self, point: &NPoint2d, angle: f64) -> Self {
        let mut result = self.clone();
        result.rotate(point, angle);
        result
    }

    fn scale(&mut self, point: &NPoint2d, factor: f64) {
        self.loc.scale(point, factor);
        if factor < 0.0 {
            self.vdir.reverse();
        }
    }

    fn scaled(&self, point: &NPoint2d, factor: f64) -> Self {
        let mut result = self.clone();
        result.scale(point, factor);
        result
    }

    fn transform(&mut self, transformation: &NTrsf2d) {
        self.loc.transform(transformation);
        self.vdir.transform(transformation);
    }

    fn transformed(&self, transformation: &NTrsf2d) -> Self {
        let mut result = self.clone();
        result.transform(transformation);
        result
    }

    fn translate_vec(&mut self, vector: &NVec2d) {
        self.loc.translate_vec(vector);
    }

    fn translated_vec(&self, vector: &NVec2d) -> Self {
        let mut result = self.clone();
        result.translate_vec(vector);
        result
    }

    fn translate_point3d(&mut self, from: &NPoint2d, to: &NPoint2d) {
        self.loc.translate_point3d(from, to);
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

    fn ax2d(x: f64, y: f64, dx: f64, dy: f64) -> NAx2d {
        NAx2d::new(
            NPoint2d::new(x, y),
            NDir2d::new(dx, dy).expect("Invalid direction"),
        )
    }

    #[test]
    fn test_x_axis() {
        let x_axis = NAx2d::x_axis();
        assert_eq!(x_axis.location(), &NPoint2d::new(0.0, 0.0));
        assert_eq!(x_axis.direction(), &NDir2d::new(1.0, 0.0).unwrap());
    }

    #[test]
    fn test_setters() {
        let mut ax = NAx2d::x_axis();
        ax.set_location(NPoint2d::new(1.0, 2.0));
        ax.set_direction(NDir2d::new(0.0, 1.0).unwrap());
        assert_eq!(ax.location(), &NPoint2d::new(1.0, 2.0));
        assert_eq!(ax.direction(), &NDir2d::new(0.0, 1.0).unwrap());
    }

    #[test]
    fn test_is_coaxial() {
        let ax1 = ax2d(0.0, 0.0, 1.0, 0.0);
        let ax2 = ax2d(1.0, 0.0, 1.0, 0.0);
        assert!(ax1.is_coaxial(&ax2, 1e-5, 1e-5));
        let ax3 = ax2d(0.0, 0.0, 0.0, 1.0);
        assert!(!ax1.is_coaxial(&ax3, 1e-5, 1e-5));
    }

    #[test]
    fn test_is_normal_opposite_parallel() {
        let ax1 = ax2d(0.0, 0.0, 1.0, 0.0);
        let ax2 = ax2d(0.0, 0.0, 0.0, 1.0);
        let ax3 = ax2d(0.0, 0.0, -1.0, 0.0);
        assert!(ax1.is_normal(&ax2, 1e-5));
        assert!(ax1.is_opposite(&ax3, 1e-5));
        assert!(ax1.is_parallel(&ax3, 1e-5));
    }

    #[test]
    fn test_angle() {
        let ax1 = ax2d(0.0, 0.0, 1.0, 0.0);
        let ax2 = ax2d(0.0, 0.0, 0.0, 1.0);
        assert!((ax1.angle(&ax2) - PI / 2.0).abs() < 1e-5);
    }

    #[test]
    fn test_reverse_reversed() {
        let mut ax = ax2d(1.0, 2.0, 1.0, 0.0);
        let reversed = ax.reversed();
        ax.reverse();
        assert_eq!(ax.direction(), &NDir2d::new(-1.0, 0.0).unwrap());
        assert_eq!(reversed.direction(), &NDir2d::new(-1.0, 0.0).unwrap());
        assert_eq!(ax.location(), reversed.location());
    }

    #[test]
    fn test_mirror_pnt() {
        let mut ax = ax2d(1.0, 0.0, 1.0, 0.0);
        let point = NPoint2d::new(0.0, 0.0);
        ax.mirror_pnt(&point);
        assert_eq!(ax.location(), &NPoint2d::new(-1.0, 0.0));
        assert_eq!(ax.direction(), &NDir2d::new(-1.0, 0.0).unwrap());
    }

    #[test]
    fn test_mirrored_pnt() {
        let ax = ax2d(1.0, 0.0, 1.0, 0.0);
        let point = NPoint2d::new(0.0, 0.0);
        let mirrored = ax.mirrored_pnt(&point);
        assert_eq!(mirrored.location(), &NPoint2d::new(-1.0, 0.0));
        assert_eq!(mirrored.direction(), &NDir2d::new(-1.0, 0.0).unwrap());
    }

    #[test]
    fn test_mirror_ax2d() {
        let mut ax = ax2d(1.0, 0.0, 1.0, 0.0);
        let axis = ax2d(0.0, 0.0, 0.0, 1.0);
        ax.mirror_ax2d(&axis);
        assert_eq!(ax.location(), &NPoint2d::new(1.0, 0.0));
        assert_eq!(ax.direction(), &NDir2d::new(-1.0, 0.0).unwrap());
    }

    #[test]
    fn test_rotate() {
        let mut ax = ax2d(1.0, 0.0, 1.0, 0.0);
        let point = NPoint2d::new(0.0, 0.0);
        ax.rotate(&point, PI / 2.0);
        assert!((ax.location().x() - 0.0).abs() < 1e-5);
        assert!((ax.location().y() - 1.0).abs() < 1e-5);
        assert!((ax.direction().x() - 0.0).abs() < 1e-5);
        assert!((ax.direction().y() - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_scale() {
        let mut ax = ax2d(1.0, 0.0, 1.0, 0.0);
        let point = NPoint2d::new(0.0, 0.0);
        ax.scale(&point, 2.0);
        assert_eq!(ax.location(), &NPoint2d::new(2.0, 0.0));
        assert_eq!(ax.direction(), &NDir2d::new(1.0, 0.0).unwrap());

        ax.scale(&point, -2.0);
        assert_eq!(ax.location(), &NPoint2d::new(-4.0, 0.0));
        assert_eq!(ax.direction(), &NDir2d::new(-1.0, 0.0).unwrap());
    }

    #[test]
    fn test_translate_vec() {
        let mut ax = ax2d(1.0, 2.0, 1.0, 0.0);
        let vec = NVec2d::new(1.0, 1.0);
        ax.translate_vec(&vec);
        assert_eq!(ax.location(), &NPoint2d::new(2.0, 3.0));
        assert_eq!(ax.direction(), &NDir2d::new(1.0, 0.0).unwrap());
    }

    #[test]
    fn test_translate_point3d() {
        let mut ax = ax2d(1.0, 2.0, 1.0, 0.0);
        let p1 = NPoint2d::new(0.0, 0.0);
        let p2 = NPoint2d::new(1.0, 1.0);
        ax.translate_point3d(&p1, &p2);
        assert_eq!(ax.location(), &NPoint2d::new(2.0, 3.0));
        assert_eq!(ax.direction(), &NDir2d::new(1.0, 0.0).unwrap());
    }
}
