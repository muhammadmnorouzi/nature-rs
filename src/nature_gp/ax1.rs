use super::prelude::*;
use serde::{Deserialize, Serialize};

// Trait to define the behavior of an axis in 3D space
pub trait Ax1 {
    fn new(location: NPoint3d, direction: NDir) -> Self;
    fn z_axis() -> Self;
    fn set_direction(&mut self, direction: NDir);
    fn set_location(&mut self, location: NPoint3d);
    fn direction(&self) -> &NDir;
    fn location(&self) -> &NPoint3d;
    fn is_coaxial(&self, other: &Self, angular_tolerance: f64, linear_tolerance: f64) -> bool;
    fn is_normal(&self, other: &Self, angular_tolerance: f64) -> bool;
    fn is_opposite(&self, other: &Self, angular_tolerance: f64) -> bool;
    fn is_parallel(&self, other: &Self, angular_tolerance: f64) -> bool;
    fn angle(&self, other: &Self) -> f64;
    fn reverse(&mut self);
    fn reversed(&self) -> Self;
    fn mirror_point3d(&mut self, point: &NPoint3d);
    fn mirrored_point3d(&self, point: &NPoint3d) -> Self;
    fn mirror_ax1(&mut self, axis: &Self);
    fn mirrored_ax1(&self, axis: &Self) -> Self;
    fn mirror_ax2(&mut self, plane: &NAx2);
    fn mirrored_ax2(&self, plane: &NAx2) -> Self;
    fn rotate(&mut self, axis: &Self, angle: f64);
    fn rotated(&self, axis: &Self, angle: f64) -> Self;
    fn scale(&mut self, point: &NPoint3d, factor: f64);
    fn scaled(&self, point: &NPoint3d, factor: f64) -> Self;
    fn transform(&mut self, transformation: &NTrsf);
    fn transformed(&self, transformation: &NTrsf) -> Self;
    fn translate_vec(&mut self, vector: &NVec);
    fn translated_vec(&self, vector: &NVec) -> Self;
    fn translate_point3d(&mut self, from: &NPoint3d, to: &NPoint3d);
    fn translated_point3d(&self, from: &NPoint3d, to: &NPoint3d) -> Self;
}

// Struct representing an axis in 3D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NAx1 {
    location: NPoint3d,
    vertical_direction: NDir,
}

impl Default for NAx1 {
    fn default() -> Self {
        Self::z_axis()
    }
}

impl Ax1 for NAx1 {
    fn new(location: NPoint3d, direction: NDir) -> Self {
        NAx1 {
            location,
            vertical_direction: direction,
        }
    }

    fn z_axis() -> Self {
        NAx1 {
            location: NPoint3d::new(0.0, 0.0, 0.0),
            vertical_direction: NDir::new(0.0, 0.0, 1.0).unwrap(),
        }
    }

    fn set_direction(&mut self, direction: NDir) {
        self.vertical_direction = direction;
    }

    fn set_location(&mut self, location: NPoint3d) {
        self.location = location;
    }

    fn direction(&self) -> &NDir {
        &self.vertical_direction
    }

    fn location(&self) -> &NPoint3d {
        &self.location
    }

    fn is_coaxial(&self, other: &Self, angular_tolerance: f64, linear_tolerance: f64) -> bool {
        let mut xyz1 = self.location.xyz();
        xyz1.subtract(&other.location.xyz());
        xyz1.cross(&other.vertical_direction.xyz());
        let d1 = xyz1.modulus();

        let mut xyz2 = other.location.xyz();
        xyz2.subtract(&self.location.xyz());
        xyz2.cross(&self.vertical_direction.xyz());
        let d2 = xyz2.modulus();

        self.vertical_direction
            .is_equal(&other.vertical_direction, angular_tolerance)
            && d1 <= linear_tolerance
            && d2 <= linear_tolerance
    }

    fn is_normal(&self, other: &Self, angular_tolerance: f64) -> bool {
        self.vertical_direction
            .is_normal(&other.vertical_direction, angular_tolerance)
    }

    fn is_opposite(&self, other: &Self, angular_tolerance: f64) -> bool {
        self.vertical_direction
            .is_opposite(&other.vertical_direction, angular_tolerance)
    }

    fn is_parallel(&self, other: &Self, angular_tolerance: f64) -> bool {
        self.vertical_direction
            .is_parallel(&other.vertical_direction, angular_tolerance)
    }

    fn angle(&self, other: &Self) -> f64 {
        self.vertical_direction.angle(&other.vertical_direction)
    }

    fn reverse(&mut self) {
        self.vertical_direction.reverse();
    }

    fn reversed(&self) -> Self {
        NAx1 {
            location: self.location.clone(),
            vertical_direction: self.vertical_direction.reversed(),
        }
    }

    fn mirror_point3d(&mut self, point: &NPoint3d) {
        self.location.mirror_point3d(point);
        self.vertical_direction.reverse();
    }

    fn mirrored_point3d(&self, point: &NPoint3d) -> Self {
        let mut result = self.clone();
        result.mirror_point3d(point);
        result
    }

    fn mirror_ax1(&mut self, axis: &Self) {
        self.location.mirror_ax1(axis);
        self.vertical_direction.mirror_dir(&axis.vertical_direction);
    }

    fn mirrored_ax1(&self, axis: &Self) -> Self {
        let mut result = self.clone();
        result.mirror_ax1(axis);
        result
    }

    fn mirror_ax2(&mut self, plane: &NAx2) {
        self.location.mirror_ax2(plane);
        self.vertical_direction.mirror_ax2(plane);
    }

    fn mirrored_ax2(&self, plane: &NAx2) -> Self {
        let mut result = self.clone();
        result.mirror_ax2(plane);
        result
    }

    fn rotate(&mut self, axis: &Self, angle: f64) {
        self.location.rotate(axis, angle);
        self.vertical_direction.rotate(axis, angle);
    }

    fn rotated(&self, axis: &Self, angle: f64) -> Self {
        let mut result = self.clone();
        result.rotate(axis, angle);
        result
    }

    fn scale(&mut self, point: &NPoint3d, factor: f64) {
        self.location.scale(point, factor);
        if factor < 0.0 {
            self.vertical_direction.reverse();
        }
    }

    fn scaled(&self, point: &NPoint3d, factor: f64) -> Self {
        let mut result = self.clone();
        result.scale(point, factor);
        result
    }

    fn transform(&mut self, transformation: &NTrsf) {
        self.location.transform(transformation);
        self.vertical_direction.transform(transformation);
    }

    fn transformed(&self, transformation: &NTrsf) -> Self {
        let mut result = self.clone();
        result.transform(transformation);
        result
    }

    fn translate_vec(&mut self, vector: &NVec) {
        self.location.translate_vec(vector);
    }

    fn translated_vec(&self, vector: &NVec) -> Self {
        let mut result = self.clone();
        result.translate_vec(vector);
        result
    }

    fn translate_point3d(&mut self, from: &NPoint3d, to: &NPoint3d) {
        self.location.translate_point3d(from, to);
    }

    fn translated_point3d(&self, from: &NPoint3d, to: &NPoint3d) -> Self {
        let mut result = self.clone();
        result.translate_point3d(from, to);
        result
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;

    fn new_ax1(x: f64, y: f64, z: f64, dx: f64, dy: f64, dz: f64) -> NAx1 {
        NAx1::new(
            NPoint3d::new(x, y, z),
            NDir::new(dx, dy, dz).expect("Invalid direction"),
        )
    }

    #[test]
    fn test_new_z_axis() {
        let z_axis = NAx1::z_axis();
        assert_eq!(z_axis.location(), &NPoint3d::new(0.0, 0.0, 0.0));
        assert_eq!(z_axis.direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
    }

    #[test]
    fn test_setters() {
        let mut ax = NAx1::z_axis();
        ax.set_location(NPoint3d::new(1.0, 2.0, 3.0));
        ax.set_direction(NDir::new(1.0, 0.0, 0.0).unwrap());
        assert_eq!(ax.location(), &NPoint3d::new(1.0, 2.0, 3.0));
        assert_eq!(ax.direction(), &NDir::new(1.0, 0.0, 0.0).unwrap());
    }

    #[test]
    fn test_is_coaxial() {
        let ax1 = new_ax1(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let ax2 = new_ax1(1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        assert!(ax1.is_coaxial(&ax2, 1e-5, 1e-5));
        let ax3 = new_ax1(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        assert!(!ax1.is_coaxial(&ax3, 1e-5, 1e-5));
    }

    #[test]
    fn test_is_normal_opposite_parallel() {
        let ax1 = new_ax1(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let ax2 = new_ax1(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        let ax3 = new_ax1(0.0, 0.0, 0.0, 0.0, 0.0, -1.0);
        assert!(ax1.is_normal(&ax2, 1e-5));
        assert!(ax1.is_opposite(&ax3, 1e-5));
        assert!(ax1.is_parallel(&ax3, 1e-5));
    }

    #[test]
    fn test_angle() {
        let ax1 = new_ax1(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let ax2 = new_ax1(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        assert!((ax1.angle(&ax2) - PI / 2.0).abs() < 1e-5);
    }

    #[test]
    fn test_reverse_reversed() {
        let mut ax = new_ax1(1.0, 2.0, 3.0, 0.0, 0.0, 1.0);
        let reversed = ax.reversed();
        ax.reverse();
        assert_eq!(ax.direction(), &NDir::new(0.0, 0.0, -1.0).unwrap());
        assert_eq!(reversed.direction(), &NDir::new(0.0, 0.0, -1.0).unwrap());
        assert_eq!(ax.location(), reversed.location());
    }

    #[test]
    fn test_mirror_point3d() {
        let mut ax = new_ax1(1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let point = NPoint3d::new(0.0, 0.0, 0.0);
        ax.mirror_point3d(&point);
        assert_eq!(ax.location(), &NPoint3d::new(-1.0, 0.0, 0.0));
        assert_eq!(ax.direction(), &NDir::new(0.0, 0.0, -1.0).unwrap());
    }

    #[test]
    fn test_mirrored_point3d() {
        let ax = new_ax1(1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let point = NPoint3d::new(0.0, 0.0, 0.0);
        let mirrored = ax.mirrored_point3d(&point);
        assert_eq!(mirrored.location(), &NPoint3d::new(-1.0, 0.0, 0.0));
        assert_eq!(mirrored.direction(), &NDir::new(0.0, 0.0, -1.0).unwrap());
    }

    #[test]
    fn test_scale() {
        let mut ax = new_ax1(1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let point = NPoint3d::new(0.0, 0.0, 0.0);
        ax.scale(&point, 2.0);
        assert_eq!(ax.location(), &NPoint3d::new(2.0, 0.0, 0.0));
        assert_eq!(ax.direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());

        ax.scale(&point, -2.0);
        assert_eq!(ax.location(), &NPoint3d::new(-4.0, 0.0, 0.0));
        assert_eq!(ax.direction(), &NDir::new(0.0, 0.0, -1.0).unwrap());
    }

    #[test]
    fn test_translate_vec() {
        let mut ax = new_ax1(1.0, 2.0, 3.0, 0.0, 0.0, 1.0);
        let vec = NVec::new_from_coords(1.0, 1.0, 1.0);
        ax.translate_vec(&vec);
        assert_eq!(ax.location(), &NPoint3d::new(2.0, 3.0, 4.0));
        assert_eq!(ax.direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
    }

    #[test]
    fn test_translate_point3d() {
        let mut ax = new_ax1(1.0, 2.0, 3.0, 0.0, 0.0, 1.0);
        let p1 = NPoint3d::new(0.0, 0.0, 0.0);
        let p2 = NPoint3d::new(1.0, 1.0, 1.0);
        ax.translate_point3d(&p1, &p2);
        assert_eq!(ax.location(), &NPoint3d::new(2.0, 3.0, 4.0));
        assert_eq!(ax.direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
    }
}
