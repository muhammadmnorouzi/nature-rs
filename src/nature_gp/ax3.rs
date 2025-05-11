use serde::{Deserialize, Serialize};

use crate::{
    nature_errors::NErrors,
    nature_gp::{NAx1, NAx2, NDir, NGP, NPoint3d, NTrsf, NVec},
};

// Trait to define the behavior of a coordinate system in 3D space
pub trait Ax3 {
    fn new(location: NPoint3d, direction: NDir, x_direction: NDir) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn new_from_point_and_direction(location: NPoint3d, direction: NDir) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn from_ax2(ax2: &NAx2) -> Self;
    fn origin() -> Self;
    fn x_reverse(&mut self);
    fn y_reverse(&mut self);
    fn z_reverse(&mut self);
    fn set_axis(&mut self, axis: &NAx1) -> Result<(), NErrors>;
    fn set_direction(&mut self, direction: NDir) -> Result<(), NErrors>;
    fn set_location(&mut self, location: NPoint3d);
    fn set_x_direction(&mut self, x_direction: NDir) -> Result<(), NErrors>;
    fn set_y_direction(&mut self, y_direction: NDir) -> Result<(), NErrors>;
    fn angle(&self, other: &Self) -> f64;
    fn axis(&self) -> &NAx1;
    fn ax2(&self) -> NAx2;
    fn direction(&self) -> &NDir;
    fn location(&self) -> &NPoint3d;
    fn x_direction(&self) -> &NDir;
    fn y_direction(&self) -> &NDir;
    fn direct(&self) -> bool;
    fn is_coplanar_ax3(&self, other: &Self, linear_tolerance: f64, angular_tolerance: f64) -> bool;
    fn is_coplanar_ax1(&self, axis: &NAx1, linear_tolerance: f64, angular_tolerance: f64) -> bool;
    fn mirror_pnt(&mut self, point: &NPoint3d);
    fn mirrored_pnt(&self, point: &NPoint3d) -> Self;
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

// Struct representing a coordinate system in 3D space (right- or left-handed)
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NAx3 {
    axis: NAx1,
    vxdir: NDir,
    vydir: NDir,
}

impl Default for NAx3 {
    fn default() -> Self {
        Self::origin()
    }
}

impl Ax3 for NAx3 {
    fn new(location: NPoint3d, direction: NDir, x_direction: NDir) -> Result<Self, NErrors> {
        // Check if x_direction is not parallel to direction
        if direction.is_parallel(&x_direction, NGP::angular()) {
            return Err(NErrors::ParallelVectors);
        }
        let mut result = NAx3 {
            axis: NAx1::new(location, direction),
            vxdir: x_direction,
            vydir: direction, // Temporary, will be updated
        };
        result.vxdir = direction.cross_crossed(&x_direction, &direction)?;
        result.vydir = direction.crossed(&result.vxdir)?;
        Ok(result)
    }

    fn new_from_point_and_direction(location: NPoint3d, direction: NDir) -> Result<Self, NErrors> {
        let (x, y, z) = (direction.x(), direction.y(), direction.z());
        let (x_abs, y_abs, z_abs) = (x.abs(), y.abs(), z.abs());
        let mut x_dir = NDir::new(0.0, 0.0, 0.0).map_err(|_| NErrors::InvalidDirection)?;

        // Compute X direction orthogonal to direction
        if y_abs <= x_abs && y_abs <= z_abs {
            if x_abs > z_abs {
                x_dir = NDir::new(-z, 0.0, x).map_err(|_| NErrors::InvalidDirection)?;
            } else {
                x_dir = NDir::new(z, 0.0, -x).map_err(|_| NErrors::InvalidDirection)?;
            }
        } else if x_abs <= y_abs && x_abs <= z_abs {
            if y_abs > z_abs {
                x_dir = NDir::new(0.0, -z, y).map_err(|_| NErrors::InvalidDirection)?;
            } else {
                x_dir = NDir::new(0.0, z, -y).map_err(|_| NErrors::InvalidDirection)?;
            }
        } else {
            if x_abs > y_abs {
                x_dir = NDir::new(-y, x, 0.0).map_err(|_| NErrors::InvalidDirection)?;
            } else {
                x_dir = NDir::new(y, -x, 0.0).map_err(|_| NErrors::InvalidDirection)?;
            }
        }

        NAx3::new(location, direction, x_dir)
    }

    fn from_ax2(ax2: &NAx2) -> Self {
        NAx3 {
            axis: ax2.axis().clone(),
            vxdir: ax2.x_direction().clone(),
            vydir: ax2.y_direction().clone(),
        }
    }

    fn origin() -> Self {
        NAx3 {
            axis: NAx1::z_axis(),
            vxdir: NDir::new(1.0, 0.0, 0.0).expect("Invalid X direction"),
            vydir: NDir::new(0.0, 1.0, 0.0).expect("Invalid Y direction"),
        }
    }

    fn x_reverse(&mut self) {
        self.vxdir.reverse();
    }

    fn y_reverse(&mut self) {
        self.vydir.reverse();
    }

    fn z_reverse(&mut self) {
        self.axis.reverse();
    }

    fn set_axis(&mut self, axis: &NAx1) -> Result<(), NErrors> {
        self.axis.set_location(axis.location().clone());
        self.set_direction(axis.direction().clone())
    }

    fn set_direction(&mut self, direction: NDir) -> Result<(), NErrors> {
        let dot = direction.dot(&self.vxdir);
        if (1.0 - dot.abs()) <= NGP::angular() {
            if dot > 0.0 {
                self.vxdir = self.vydir.clone();
                self.vydir = self.axis.direction().clone();
            } else {
                self.vxdir = self.axis.direction().clone();
            }
            self.axis.set_direction(direction);
        } else {
            let direct = self.direct();
            self.axis.set_direction(direction);
            self.vxdir = self
                .axis
                .direction()
                .cross_crossed(&self.vxdir, &self.axis.direction())?;
            self.vydir = if direct {
                self.axis.direction().crossed(&self.vxdir)?
            } else {
                self.vxdir.crossed(&self.axis.direction())?
            };
        }
        Ok(())
    }

    fn set_location(&mut self, location: NPoint3d) {
        self.axis.set_location(location);
    }

    fn set_x_direction(&mut self, x_direction: NDir) -> Result<(), NErrors> {
        let dot = x_direction.dot(&self.axis.direction());
        if (1.0 - dot.abs()) <= NGP::angular() {
            if dot > 0.0 {
                self.axis.set_direction(self.vxdir.clone());
                self.vydir.reverse();
            } else {
                self.axis.set_direction(self.vxdir.clone());
            }
            self.vxdir = x_direction;
        } else {
            let direct = self.direct();
            self.vxdir = self
                .axis
                .direction()
                .cross_crossed(&x_direction, &self.axis.direction())?;
            self.vydir = if direct {
                self.axis.direction().crossed(&self.vxdir)?
            } else {
                self.vxdir.crossed(&self.axis.direction())?
            };
        }
        Ok(())
    }

    fn set_y_direction(&mut self, y_direction: NDir) -> Result<(), NErrors> {
        let dot = y_direction.dot(&self.axis.direction());
        if (1.0 - dot.abs()) <= NGP::angular() {
            if dot > 0.0 {
                self.axis.set_direction(self.vydir.clone());
                self.vxdir.reverse();
            } else {
                self.axis.set_direction(self.vydir.clone());
            }
            self.vydir = y_direction;
        } else {
            let direct = self.direct();
            self.vxdir = y_direction.crossed(&self.axis.direction())?;
            self.vydir = if direct {
                self.axis.direction().crossed(&self.vxdir)?
            } else {
                self.vxdir.crossed(&self.axis.direction())?
            };
            if !direct {
                self.vxdir.reverse();
            }
        }
        Ok(())
    }

    fn angle(&self, other: &Self) -> f64 {
        self.axis.angle(&other.axis)
    }

    fn axis(&self) -> &NAx1 {
        &self.axis
    }

    fn ax2(&self) -> NAx2 {
        let mut z_dir = self.axis.direction().clone();
        if !self.direct() {
            z_dir.reverse();
        }
        NAx2::new(self.axis.location().clone(), z_dir, self.vxdir.clone()).expect("Invalid NAx2")
    }

    fn direction(&self) -> &NDir {
        self.axis.direction()
    }

    fn location(&self) -> &NPoint3d {
        self.axis.location()
    }

    fn x_direction(&self) -> &NDir {
        &self.vxdir
    }

    fn y_direction(&self) -> &NDir {
        &self.vydir
    }

    fn direct(&self) -> bool {
        self.vxdir
            .crossed(&self.vydir)
            .unwrap()
            .dot(&self.axis.direction())
            > 0.0
    }

    fn is_coplanar_ax3(&self, other: &Self, linear_tolerance: f64, angular_tolerance: f64) -> bool {
        let vec = NVec::new_from_points(self.axis.location(), other.axis.location());
        let d1 = NVec::new_from_dir(self.axis.direction()).dot(&vec).abs();
        let d2 = NVec::new_from_dir(other.axis.direction()).dot(&vec).abs();
        d1 <= linear_tolerance
            && d2 <= linear_tolerance
            && self.axis.is_parallel(&other.axis, angular_tolerance)
    }

    fn is_coplanar_ax1(&self, axis: &NAx1, linear_tolerance: f64, angular_tolerance: f64) -> bool {
        let vec = NVec::new_from_points(self.axis.location(), axis.location());
        let d1 = NVec::new_from_dir(self.axis.direction()).dot(&vec).abs();
        let d2 = NVec::new_from_dir(axis.direction())
            .crossed(&vec)
            .magnitude()
            .abs();
        d1 <= linear_tolerance
            && d2 <= linear_tolerance
            && self.axis.is_normal(axis, angular_tolerance)
    }

    fn mirror_pnt(&mut self, point: &NPoint3d) {
        self.axis.mirror_pnt(point);
        self.vxdir.reverse();
        self.vydir.reverse();
    }

    fn mirrored_pnt(&self, point: &NPoint3d) -> Self {
        let mut result = self.clone();
        result.mirror_pnt(point);
        result
    }

    fn mirror_ax1(&mut self, axis: &NAx1) {
        self.vydir.mirror_ax1(axis);
        self.vxdir.mirror_ax1(axis);
        self.axis.mirror_ax1(axis);
    }

    fn mirrored_ax1(&self, axis: &NAx1) -> Self {
        let mut result = self.clone();
        result.mirror_ax1(axis);
        result
    }

    fn mirror_ax2(&mut self, plane: &NAx2) {
        self.vydir.mirror_ax2(plane);
        self.vxdir.mirror_ax2(plane);
        self.axis.mirror_ax2(plane);
    }

    fn mirrored_ax2(&self, plane: &NAx2) -> Self {
        let mut result = self.clone();
        result.mirror_ax2(plane);
        result
    }

    fn rotate(&mut self, axis: &NAx1, angle: f64) {
        self.axis.rotate(axis, angle);
        self.vxdir.rotate(axis, angle);
        self.vydir.rotate(axis, angle);
    }

    fn rotated(&self, axis: &NAx1, angle: f64) -> Self {
        let mut result = self.clone();
        result.rotate(axis, angle);
        result
    }

    fn scale(&mut self, point: &NPoint3d, factor: f64) {
        self.axis.scale(point, factor);
        if factor < 0.0 {
            self.vxdir.reverse();
            self.vydir.reverse();
        }
    }

    fn scaled(&self, point: &NPoint3d, factor: f64) -> Self {
        let mut result = self.clone();
        result.scale(point, factor);
        result
    }

    fn transform(&mut self, transformation: &NTrsf) {
        self.axis.transform(transformation);
        self.vxdir.transform(transformation);
        self.vydir.transform(transformation);
    }

    fn transformed(&self, transformation: &NTrsf) -> Self {
        let mut result = self.clone();
        result.transform(transformation);
        result
    }

    fn translate_vec(&mut self, vector: &NVec) {
        self.axis.translate_vec(vector);
    }

    fn translated_vec(&self, vector: &NVec) -> Self {
        let mut result = self.clone();
        result.translate_vec(vector);
        result
    }

    fn translate_point3d(&mut self, from: &NPoint3d, to: &NPoint3d) {
        self.axis.translate_point3d(from, to);
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

    fn ax3(loc: (f64, f64, f64), dir: (f64, f64, f64), x_dir: (f64, f64, f64)) -> NAx3 {
        NAx3::new(
            NPoint3d::new(loc.0, loc.1, loc.2),
            NDir::new(dir.0, dir.1, dir.2).expect("Invalid direction"),
            NDir::new(x_dir.0, x_dir.1, x_dir.2).expect("Invalid X direction"),
        )
        .expect("Failed to create NAx3")
    }

    #[test]
    fn test_origin() {
        let origin = NAx3::origin();
        assert_eq!(origin.location(), &NPoint3d::new(0.0, 0.0, 0.0));
        assert_eq!(origin.direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
        assert_eq!(origin.x_direction(), &NDir::new(1.0, 0.0, 0.0).unwrap());
        assert_eq!(origin.y_direction(), &NDir::new(0.0, 1.0, 0.0).unwrap());
        assert!(origin.direct());
    }

    #[test]
    fn test_new_from_point_and_direction() {
        let ax3 = NAx3::new_from_point_and_direction(
            NPoint3d::new(1.0, 2.0, 3.0),
            NDir::new(0.0, 0.0, 1.0).unwrap(),
        )
        .unwrap();
        assert_eq!(ax3.location(), &NPoint3d::new(1.0, 2.0, 3.0));
        assert_eq!(ax3.direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
        assert!(ax3.x_direction().is_normal(ax3.direction(), 1e-5));
        assert!(ax3.y_direction().is_normal(ax3.direction(), 1e-5));
        assert!(ax3.x_direction().is_normal(ax3.y_direction(), 1e-5));
    }

    #[test]
    fn test_from_ax2() {
        let ax2 = NAx2::new(
            NPoint3d::new(1.0, 2.0, 3.0),
            NDir::new(0.0, 0.0, 1.0).unwrap(),
            NDir::new(1.0, 0.0, 0.0).unwrap(),
        )
        .unwrap();
        let ax3 = NAx3::from_ax2(&ax2);
        assert_eq!(ax3.location(), ax2.location());
        assert_eq!(ax3.direction(), ax2.direction());
        assert_eq!(ax3.x_direction(), ax2.x_direction());
        assert_eq!(ax3.y_direction(), ax2.y_direction());
        assert!(ax3.direct());
    }

    #[test]
    fn test_reverse() {
        let mut ax3 = NAx3::origin();
        ax3.x_reverse();
        assert_eq!(ax3.x_direction(), &NDir::new(-1.0, 0.0, 0.0).unwrap());
        ax3.y_reverse();
        assert_eq!(ax3.y_direction(), &NDir::new(0.0, -1.0, 0.0).unwrap());
        ax3.z_reverse();
        assert_eq!(ax3.direction(), &NDir::new(0.0, 0.0, -1.0).unwrap());
    }

    #[test]
    fn test_set_axis() {
        let mut ax3 = NAx3::origin();
        let new_axis = NAx1::new(
            NPoint3d::new(1.0, 0.0, 0.0),
            NDir::new(1.0, 0.0, 0.0).unwrap(),
        );
        ax3.set_axis(&new_axis).unwrap();
        assert_eq!(ax3.axis(), &new_axis);
        assert!(ax3.x_direction().is_normal(ax3.direction(), 1e-5));
        assert!(ax3.y_direction().is_normal(ax3.direction(), 1e-5));
    }

    #[test]
    fn test_set_direction() {
        let mut ax3 = NAx3::origin();
        let new_dir = NDir::new(1.0, 0.0, 0.0).unwrap();
        ax3.set_direction(new_dir.clone()).unwrap();
        assert_eq!(ax3.direction(), &new_dir);
        assert!(ax3.x_direction().is_normal(&new_dir, 1e-5));
        assert!(ax3.y_direction().is_normal(&new_dir, 1e-5));
    }

    #[test]
    fn test_set_x_direction() {
        let mut ax3 = NAx3::origin();
        let new_x_dir = NDir::new(0.0, 1.0, 0.0).unwrap();
        ax3.set_x_direction(new_x_dir.clone()).unwrap();
        assert_eq!(ax3.x_direction(), &new_x_dir);
        assert!(ax3.y_direction().is_normal(ax3.direction(), 1e-5));
        assert!(ax3.x_direction().is_normal(ax3.y_direction(), 1e-5));
    }

    #[test]
    fn test_ax2() {
        let ax3 = ax3((0.0, 0.0, 0.0), (0.0, 0.0, -1.0), (-1.0, 0.0, 0.0));
        let ax2 = ax3.ax2();
        assert_eq!(ax2.direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
        assert_eq!(ax2.x_direction(), &NDir::new(-1.0, 0.0, 0.0).unwrap());
        assert!(ax2.y_direction().is_normal(ax2.direction(), 1e-5));
    }

    #[test]
    fn test_direct() {
        let ax3 = NAx3::origin();
        assert!(ax3.direct());
        let ax3_left = ax3((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (-1.0, 0.0, 0.0));
        assert!(!ax3_left.direct());
    }

    #[test]
    fn test_is_coplanar_ax3() {
        let ax3_1 = ax3((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        let ax3_2 = ax3((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        assert!(ax3_1.is_coplanar_ax3(&ax3_2, 1e-5, 1e-5));
    }

    #[test]
    fn test_is_coplanar_ax1() {
        let ax3 = ax3((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        let ax1 = NAx1::new(
            NPoint3d::new(1.0, 0.0, 0.0),
            NDir::new(1.0, 0.0, 0.0).unwrap(),
        );
        assert!(ax3.is_coplanar_ax1(&ax1, 1e-5, 1e-5));
    }

    #[test]
    fn test_mirror_pnt() {
        let mut ax3 = ax3((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        let point = NPoint3d::new(0.0, 0.0, 0.0);
        ax3.mirror_pnt(&point);
        assert_eq!(ax3.location(), &NPoint3d::new(-1.0, 0.0, 0.0));
        assert_eq!(ax3.x_direction(), &NDir::new(-1.0, 0.0, 0.0).unwrap());
        assert_eq!(ax3.y_direction(), &NDir::new(0.0, -1.0, 0.0).unwrap());
    }

    #[test]
    fn test_scale() {
        let mut ax3 = ax3((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        let point = NPoint3d::new(0.0, 0.0, 0.0);
        ax3.scale(&point, 2.0);
        assert_eq!(ax3.location(), &NPoint3d::new(2.0, 0.0, 0.0));
        assert_eq!(ax3.direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());

        ax3.scale(&point, -2.0);
        assert_eq!(ax3.location(), &NPoint3d::new(-4.0, 0.0, 0.0));
        assert_eq!(ax3.x_direction(), &NDir::new(-1.0, 0.0, 0.0).unwrap());
        assert_eq!(ax3.y_direction(), &NDir::new(0.0, -1.0, 0.0).unwrap());
    }

    #[test]
    fn test_translate_vec() {
        let mut ax3 = ax3((1.0, 2.0, 3.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        let vec = NVec::new(1.0, 1.0, 1.0);
        ax3.translate_vec(&vec);
        assert_eq!(ax3.location(), &NPoint3d::new(2.0, 3.0, 4.0));
        assert_eq!(ax3.direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
    }
}
