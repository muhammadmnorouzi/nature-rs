use serde::{Deserialize, Serialize};

use crate::{
    nature_errors::NErrors,
    nature_gp::{NAx2d, NDir2d, NPoint2d, NTrsf2d, NVec2d},
};

// Trait to define the behavior of a 2D coordinate system
pub trait Ax22d {
    fn new(location: NPoint2d, x_direction: NDir2d, y_direction: NDir2d) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn new_with_sense(location: NPoint2d, direction: NDir2d, is_right_handed: bool) -> Self;
    fn from_ax2d_with_sense(axis: &NAx2d, is_right_handed: bool) -> Self;
    fn origin() -> Self;
    fn set_axis(&mut self, axis: &Self);
    fn set_x_axis(&mut self, axis: &NAx2d);
    fn set_y_axis(&mut self, axis: &NAx2d);
    fn set_location(&mut self, location: NPoint2d);
    fn set_x_direction(&mut self, x_direction: NDir2d);
    fn set_y_direction(&mut self, y_direction: NDir2d);
    fn x_axis(&self) -> NAx2d;
    fn y_axis(&self) -> NAx2d;
    fn location(&self) -> &NPoint2d;
    fn x_direction(&self) -> &NDir2d;
    fn y_direction(&self) -> &NDir2d;
    fn mirror_point3d(&mut self, point: &NPoint2d);
    fn mirrored_point3d(&self, point: &NPoint2d) -> Self;
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

// Struct representing a 2D coordinate system (right- or left-handed)
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NAx22d {
    point: NPoint2d,
    vxdir: NDir2d,
    vydir: NDir2d,
}

impl Default for NAx22d {
    fn default() -> Self {
        Self::origin()
    }
}

impl Ax22d for NAx22d {
    fn new(location: NPoint2d, x_direction: NDir2d, y_direction: NDir2d) -> Result<Self, NErrors> {
        // Check if x_direction and y_direction are parallel
        if x_direction.is_parallel(&y_direction, crate::gp::NGP::angular()) {
            return Err(NErrors::ParallelVectors);
        }
        let is_right_handed = x_direction.crossed(&y_direction) >= 0.0;
        let mut result = NAx22d {
            point: location,
            vxdir: x_direction,
            vydir: y_direction,
        };
        // Adjust Y direction to maintain handedness
        if is_right_handed {
            result.vydir = NDir2d::new(-x_direction.y(), x_direction.x())?;
        } else {
            result.vydir = NDir2d::new(x_direction.y(), -x_direction.x())?;
        }
        Ok(result)
    }

    fn new_with_sense(location: NPoint2d, direction: NDir2d, is_right_handed: bool) -> Self {
        let vydir = if is_right_handed {
            NDir2d::new(-direction.y(), direction.x()).expect("Invalid Y direction")
        } else {
            NDir2d::new(direction.y(), -direction.x()).expect("Invalid Y direction")
        };
        NAx22d {
            point: location,
            vxdir: direction,
            vydir,
        }
    }

    fn from_ax2d_with_sense(axis: &NAx2d, is_right_handed: bool) -> Self {
        Self::new_with_sense(
            axis.location().clone(),
            axis.direction().clone(),
            is_right_handed,
        )
    }

    fn origin() -> Self {
        NAx22d {
            point: NPoint2d::new(0.0, 0.0),
            vxdir: NDir2d::new(1.0, 0.0).expect("Invalid X direction"),
            vydir: NDir2d::new(0.0, 1.0).expect("Invalid Y direction"),
        }
    }

    fn set_axis(&mut self, axis: &Self) {
        self.point = axis.location().clone();
        self.vxdir = axis.x_direction().clone();
        self.vydir = axis.y_direction().clone();
    }

    fn set_x_axis(&mut self, axis: &NAx2d) {
        let is_right_handed = self.vxdir.crossed(&self.vydir) >= 0.0;
        self.point = axis.location().clone();
        self.vxdir = axis.direction().clone();
        self.vydir = if is_right_handed {
            NDir2d::new(-self.vxdir.y(), self.vxdir.x()).expect("Invalid Y direction")
        } else {
            NDir2d::new(self.vxdir.y(), -self.vxdir.x()).expect("Invalid Y direction")
        };
    }

    fn set_y_axis(&mut self, axis: &NAx2d) {
        let is_right_handed = self.vxdir.crossed(&self.vydir) >= 0.0;
        self.point = axis.location().clone();
        self.vydir = axis.direction().clone();
        self.vxdir = if is_right_handed {
            NDir2d::new(self.vydir.y(), -self.vydir.x()).expect("Invalid X direction")
        } else {
            NDir2d::new(-self.vydir.y(), self.vydir.x()).expect("Invalid X direction")
        };
    }

    fn set_location(&mut self, location: NPoint2d) {
        self.point = location;
    }

    fn set_x_direction(&mut self, x_direction: NDir2d) {
        let is_right_handed = self.vxdir.crossed(&self.vydir) >= 0.0;
        self.vxdir = x_direction;
        self.vydir = if is_right_handed {
            NDir2d::new(-x_direction.y(), x_direction.x()).expect("Invalid Y direction")
        } else {
            NDir2d::new(x_direction.y(), -x_direction.x()).expect("Invalid Y direction")
        };
    }

    fn set_y_direction(&mut self, y_direction: NDir2d) {
        let is_right_handed = self.vxdir.crossed(&self.vydir) >= 0.0;
        self.vydir = y_direction;
        self.vxdir = if is_right_handed {
            NDir2d::new(y_direction.y(), -y_direction.x()).expect("Invalid X direction")
        } else {
            NDir2d::new(-y_direction.y(), y_direction.x()).expect("Invalid X direction")
        };
    }

    fn x_axis(&self) -> NAx2d {
        NAx2d::new(self.point.clone(), self.vxdir.clone())
    }

    fn y_axis(&self) -> NAx2d {
        NAx2d::new(self.point.clone(), self.vydir.clone())
    }

    fn location(&self) -> &NPoint2d {
        &self.point
    }

    fn x_direction(&self) -> &NDir2d {
        &self.vxdir
    }

    fn y_direction(&self) -> &NDir2d {
        &self.vydir
    }

    fn mirror_point3d(&mut self, point: &NPoint2d) {
        self.point.mirror_point3d(point);
        self.vxdir.reverse();
        self.vydir.reverse();
    }

    fn mirrored_point3d(&self, point: &NPoint2d) -> Self {
        let mut result = self.clone();
        result.mirror_point3d(point);
        result
    }

    fn mirror_ax2d(&mut self, axis: &NAx2d) {
        self.vydir.mirror_ax2d(axis);
        self.vxdir.mirror_ax2d(axis);
        self.point.mirror_ax2d(axis);
    }

    fn mirrored_ax2d(&self, axis: &NAx2d) -> Self {
        let mut result = self.clone();
        result.mirror_ax2d(axis);
        result
    }

    fn rotate(&mut self, point: &NPoint2d, angle: f64) {
        self.point.rotate(point, angle);
        self.vxdir.rotate(angle);
        self.vydir.rotate(angle);
    }

    fn rotated(&self, point: &NPoint2d, angle: f64) -> Self {
        let mut result = self.clone();
        result.rotate(point, angle);
        result
    }

    fn scale(&mut self, point: &NPoint2d, factor: f64) {
        self.point.scale(point, factor);
        if factor < 0.0 {
            self.vxdir.reverse();
            self.vydir.reverse();
        }
    }

    fn scaled(&self, point: &NPoint2d, factor: f64) -> Self {
        let mut result = self.clone();
        result.scale(point, factor);
        result
    }

    fn transform(&mut self, transformation: &NTrsf2d) {
        self.point.transform(transformation);
        self.vxdir.transform(transformation);
        self.vydir.transform(transformation);
    }

    fn transformed(&self, transformation: &NTrsf2d) -> Self {
        let mut result = self.clone();
        result.transform(transformation);
        result
    }

    fn translate_vec(&mut self, vector: &NVec2d) {
        self.point.translate_vec(vector);
    }

    fn translated_vec(&self, vector: &NVec2d) -> Self {
        let mut result = self.clone();
        result.translate_vec(vector);
        result
    }

    fn translate_point3d(&mut self, from: &NPoint2d, to: &NPoint2d) {
        self.point.translate_point3d(from, to);
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

    fn ax22d(loc: (f64, f64), x_dir: (f64, f64), y_dir: (f64, f64)) -> NAx22d {
        NAx22d::new(
            NPoint2d::new(loc.0, loc.1),
            NDir2d::new(x_dir.0, x_dir.1).expect("Invalid X direction"),
            NDir2d::new(y_dir.0, y_dir.1).expect("Invalid Y direction"),
        )
        .expect("Failed to create NAx22d")
    }

    #[test]
    fn test_origin() {
        let origin = NAx22d::origin();
        assert_eq!(origin.location(), &NPoint2d::new(0.0, 0.0));
        assert_eq!(origin.x_direction(), &NDir2d::new(1.0, 0.0).unwrap());
        assert_eq!(origin.y_direction(), &NDir2d::new(0.0, 1.0).unwrap());
    }

    #[test]
    fn test_new_with_sense() {
        let ax = NAx22d::new_with_sense(
            NPoint2d::new(1.0, 2.0),
            NDir2d::new(1.0, 0.0).unwrap(),
            true,
        );
        assert_eq!(ax.location(), &NPoint2d::new(1.0, 2.0));
        assert_eq!(ax.x_direction(), &NDir2d::new(1.0, 0.0).unwrap());
        assert_eq!(ax.y_direction(), &NDir2d::new(0.0, 1.0).unwrap());

        let ax_left = NAx22d::new_with_sense(
            NPoint2d::new(1.0, 2.0),
            NDir2d::new(1.0, 0.0).unwrap(),
            false,
        );
        assert_eq!(ax_left.y_direction(), &NDir2d::new(0.0, -1.0).unwrap());
    }

    #[test]
    fn test_from_ax2d_with_sense() {
        let ax2d = NAx2d::new(NPoint2d::new(1.0, 2.0), NDir2d::new(1.0, 0.0).unwrap());
        let ax = NAx22d::from_ax2d_with_sense(&ax2d, true);
        assert_eq!(ax.location(), &NPoint2d::new(1.0, 2.0));
        assert_eq!(ax.x_direction(), &NDir2d::new(1.0, 0.0).unwrap());
        assert_eq!(ax.y_direction(), &NDir2d::new(0.0, 1.0).unwrap());
    }

    #[test]
    fn test_set_axis() {
        let mut ax = NAx22d::origin();
        let new_ax = ax22d((1.0, 2.0), (0.0, 1.0), (1.0, 0.0));
        ax.set_axis(&new_ax);
        assert_eq!(ax.location(), &NPoint2d::new(1.0, 2.0));
        assert_eq!(ax.x_direction(), &NDir2d::new(0.0, 1.0).unwrap());
        assert_eq!(ax.y_direction(), &NDir2d::new(1.0, 0.0).unwrap());
    }

    #[test]
    fn test_set_x_axis() {
        let mut ax = NAx22d::origin();
        let new_axis = NAx2d::new(NPoint2d::new(1.0, 2.0), NDir2d::new(0.0, 1.0).unwrap());
        ax.set_x_axis(&new_axis);
        assert_eq!(ax.location(), &NPoint2d::new(1.0, 2.0));
        assert_eq!(ax.x_direction(), &NDir2d::new(0.0, 1.0).unwrap());
        assert_eq!(ax.y_direction(), &NDir2d::new(1.0, 0.0).unwrap());
    }

    #[test]
    fn test_set_y_axis() {
        let mut ax = NAx22d::origin();
        let new_axis = NAx2d::new(NPoint2d::new(1.0, 2.0), NDir2d::new(0.0, 1.0).unwrap());
        ax.set_y_axis(&new_axis);
        assert_eq!(ax.location(), &NPoint2d::new(1.0, 2.0));
        assert_eq!(ax.y_direction(), &NDir2d::new(0.0, 1.0).unwrap());
        assert_eq!(ax.x_direction(), &NDir2d::new(-1.0, 0.0).unwrap());
    }

    #[test]
    fn test_set_x_direction() {
        let mut ax = NAx22d::origin();
        let new_x_dir = NDir2d::new(0.0, 1.0).unwrap();
        ax.set_x_direction(new_x_dir.clone());
        assert_eq!(ax.x_direction(), &new_x_dir);
        assert_eq!(ax.y_direction(), &NDir2d::new(1.0, 0.0).unwrap());
    }

    #[test]
    fn test_set_y_direction() {
        let mut ax = NAx22d::origin();
        let new_y_dir = NDir2d::new(1.0, 0.0).unwrap();
        ax.set_y_direction(new_y_dir.clone());
        assert_eq!(ax.y_direction(), &new_y_dir);
        assert_eq!(ax.x_direction(), &NDir2d::new(0.0, -1.0).unwrap());
    }

    #[test]
    fn test_x_axis_y_axis() {
        let ax = ax22d((1.0, 2.0), (1.0, 0.0), (0.0, 1.0));
        let x_axis = ax.x_axis();
        let y_axis = ax.y_axis();
        assert_eq!(x_axis.location(), &NPoint2d::new(1.0, 2.0));
        assert_eq!(x_axis.direction(), &NDir2d::new(1.0, 0.0).unwrap());
        assert_eq!(y_axis.location(), &NPoint2d::new(1.0, 2.0));
        assert_eq!(y_axis.direction(), &NDir2d::new(0.0, 1.0).unwrap());
    }

    #[test]
    fn test_mirror_point3d() {
        let mut ax = ax22d((1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
        let point = NPoint2d::new(0.0, 0.0);
        ax.mirror_point3d(&point);
        assert_eq!(ax.location(), &NPoint2d::new(-1.0, 0.0));
        assert_eq!(ax.x_direction(), &NDir2d::new(-1.0, 0.0).unwrap());
        assert_eq!(ax.y_direction(), &NDir2d::new(0.0, -1.0).unwrap());
    }

    #[test]
    fn test_mirror_ax2d() {
        let mut ax = ax22d((1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
        let axis = NAx2d::new(NPoint2d::new(0.0, 0.0), NDir2d::new(0.0, 1.0).unwrap());
        ax.mirror_ax2d(&axis);
        assert_eq!(ax.location(), &NPoint2d::new(1.0, 0.0));
        assert_eq!(ax.x_direction(), &NDir2d::new(-1.0, 0.0).unwrap());
        assert_eq!(ax.y_direction(), &NDir2d::new(0.0, 1.0).unwrap());
    }

    #[test]
    fn test_rotate() {
        let mut ax = ax22d((1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
        let point = NPoint2d::new(0.0, 0.0);
        ax.rotate(&point, PI / 2.0);
        assert!((ax.location().x() - 0.0).abs() < 1e-5);
        assert!((ax.location().y() - 1.0).abs() < 1e-5);
        assert!((ax.x_direction().x() - 0.0).abs() < 1e-5);
        assert!((ax.x_direction().y() - 1.0).abs() < 1e-5);
        assert!((ax.y_direction().x() + 1.0).abs() < 1e-5);
        assert!((ax.y_direction().y() - 0.0).abs() < 1e-5);
    }

    #[test]
    fn test_scale() {
        let mut ax = ax22d((1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
        let point = NPoint2d::new(0.0, 0.0);
        ax.scale(&point, 2.0);
        assert_eq!(ax.location(), &NPoint2d::new(2.0, 0.0));
        assert_eq!(ax.x_direction(), &NDir2d::new(1.0, 0.0).unwrap());

        ax.scale(&point, -2.0);
        assert_eq!(ax.location(), &NPoint2d::new(-4.0, 0.0));
        assert_eq!(ax.x_direction(), &NDir2d::new(-1.0, 0.0).unwrap());
        assert_eq!(ax.y_direction(), &NDir2d::new(0.0, -1.0).unwrap());
    }

    #[test]
    fn test_translate_vec() {
        let mut ax = ax22d((1.0, 2.0), (1.0, 0.0), (0.0, 1.0));
        let vec = NVec2d::new(1.0, 1.0);
        ax.translate_vec(&vec);
        assert_eq!(ax.location(), &NPoint2d::new(2.0, 3.0));
        assert_eq!(ax.x_direction(), &NDir2d::new(1.0, 0.0).unwrap());
    }

    #[test]
    fn test_translate_point3d() {
        let mut ax = ax22d((1.0, 2.0), (1.0, 0.0), (0.0, 1.0));
        let p1 = NPoint2d::new(0.0, 0.0);
        let p2 = NPoint2d::new(1.0, 1.0);
        ax.translate_point3d(&p1, &p2);
        assert_eq!(ax.location(), &NPoint2d::new(2.0, 3.0));
        assert_eq!(ax.x_direction(), &NDir2d::new(1.0, 0.0).unwrap());
    }
}
