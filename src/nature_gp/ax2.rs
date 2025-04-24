use serde::{Deserialize, Serialize};

use crate::{
    nature_errors::NErrors,
    nature_gp::{NAx1, NDir, NGP, NPoint3d, NTrsf, NVec},
};

// Trait to define the behavior of a right-handed coordinate system in 3D space
pub trait Ax2 {
    fn new(location: NPoint3d, direction: NDir, x_direction: NDir) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn new_from_point_and_direction(location: NPoint3d, direction: NDir) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn origin() -> Self;
    fn set_axis(&mut self, axis: &NAx1) -> Result<(), NErrors>;
    fn set_direction(&mut self, direction: NDir) -> Result<(), NErrors>;
    fn set_location(&mut self, location: NPoint3d);
    fn set_x_direction(&mut self, x_direction: NDir) -> Result<(), NErrors>;
    fn set_y_direction(&mut self, y_direction: NDir) -> Result<(), NErrors>;
    fn angle(&self, other: &Self) -> f64;
    fn axis(&self) -> &NAx1;
    fn direction(&self) -> &NDir;
    fn location(&self) -> &NPoint3d;
    fn x_direction(&self) -> &NDir;
    fn y_direction(&self) -> &NDir;
    fn is_coplanar_ax2(&self, other: &Self, linear_tolerance: f64, angular_tolerance: f64) -> bool;
    fn is_coplanar_ax1(&self, axis: &NAx1, linear_tolerance: f64, angular_tolerance: f64) -> bool;
    fn mirror_pnt(&mut self, point: &NPoint3d);
    fn mirrored_pnt(&self, point: &NPoint3d) -> Self;
    fn mirror_ax1(&mut self, axis: &NAx1);
    fn mirrored_ax1(&self, axis: &NAx1) -> Self;
    fn mirror_ax2(&mut self, plane: &Self);
    fn mirrored_ax2(&self, plane: &Self) -> Self;
    fn rotate(&mut self, axis: &NAx1, angle: f64);
    fn rotated(&self, axis: &NAx1, angle: f64) -> Self;
    fn scale(&mut self, point: &NPoint3d, factor: f64);
    fn scaled(&self, point: &NPoint3d, factor: f64) -> Self;
    fn transform(&mut self, transformation: &NTrsf);
    fn transformed(&self, transformation: &NTrsf) -> Self;
    fn translate_vec(&mut self, vector: &NVec);
    fn translated_vec(&self, vector: &NVec) -> Self;
    fn translate_pnts(&mut self, from: &NPoint3d, to: &NPoint3d);
    fn translated_pnts(&self, from: &NPoint3d, to: &NPoint3d) -> Self;
}

// Struct representing a right-handed coordinate system in 3D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NAx2 {
    axis: NAx1,
    vxdir: NDir,
    vydir: NDir,
}

impl Default for NAx2 {
    fn default() -> Self {
        Self::origin()
    }
}

impl Ax2 for NAx2 {
    fn new(location: NPoint3d, direction: NDir, x_direction: NDir) -> Result<Self, NErrors> {
        // Check if x_direction is not parallel to direction
        if direction.is_parallel(&x_direction, NGP::angular()) {
            return Err(NErrors::ParallelVectors);
        }
        let mut result = NAx2 {
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

        NAx2::new(location, direction, x_dir)
    }

    fn origin() -> Self {
        NAx2 {
            axis: NAx1::z_axis(),
            vxdir: NDir::new(1.0, 0.0, 0.0).expect("Invalid X direction"),
            vydir: NDir::new(0.0, 1.0, 0.0).expect("Invalid Y direction"),
        }
    }

    fn set_axis(&mut self, axis: &NAx1) -> Result<(), NErrors> {
        let dot = self.axis.direction().dot(&self.vxdir);
        if (dot.abs() - 1.0).abs() <= NGP::angular() {
            if dot > 0.0 {
                self.vxdir = self.vydir.clone();
                self.vydir = self.axis.direction().clone();
                self.axis = axis.clone();
            } else {
                self.vxdir = self.axis.direction().clone();
                self.axis = axis.clone();
            }
        } else {
            self.axis = axis.clone();
            self.vxdir = self
                .axis
                .direction()
                .cross_crossed(&self.vxdir, &self.axis.direction())?;
            self.vydir = self.axis.direction().crossed(&self.vxdir)?;
        }
        Ok(())
    }

    fn set_direction(&mut self, direction: NDir) -> Result<(), NErrors> {
        let dot = direction.dot(&self.vxdir);
        if (dot.abs() - 1.0).abs() <= NGP::angular() {
            if dot > 0.0 {
                self.vxdir = self.vydir.clone();
                self.vydir = self.axis.direction().clone();
                self.axis.set_direction(direction);
            } else {
                self.vxdir = self.axis.direction().clone();
                self.axis.set_direction(direction);
            }
        } else {
            self.axis.set_direction(direction);
            self.vxdir = self
                .axis
                .direction()
                .cross_crossed(&self.vxdir, &self.axis.direction())?;
            self.vydir = self.axis.direction().crossed(&self.vxdir)?;
        }
        Ok(())
    }

    fn set_location(&mut self, location: NPoint3d) {
        self.axis.set_location(location);
    }

    fn set_x_direction(&mut self, x_direction: NDir) -> Result<(), NErrors> {
        if self
            .axis
            .direction()
            .is_parallel(&x_direction, NGP::angular())
        {
            return Err(NErrors::ParallelVectors);
        }
        self.vxdir = self
            .axis
            .direction()
            .cross_crossed(&x_direction, &self.axis.direction())?;
        self.vydir = self.axis.direction().crossed(&self.vxdir)?;
        Ok(())
    }

    fn set_y_direction(&mut self, y_direction: NDir) -> Result<(), NErrors> {
        if self
            .axis
            .direction()
            .is_parallel(&y_direction, NGP::angular())
        {
            return Err(NErrors::ParallelVectors);
        }
        self.vxdir = y_direction.crossed(&self.axis.direction())?;
        self.vydir = self.axis.direction().crossed(&self.vxdir)?;
        Ok(())
    }

    fn angle(&self, other: &Self) -> f64 {
        self.axis.angle(&other.axis)
    }

    fn axis(&self) -> &NAx1 {
        &self.axis
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

    fn is_coplanar_ax2(&self, other: &Self, linear_tolerance: f64, angular_tolerance: f64) -> bool {
        let dd = self.axis.direction();
        let pp = self.axis.location();
        let op = other.axis.location();
        let d1 =
            dd.x() * (op.x() - pp.x()) + dd.y() * (op.y() - pp.y()) + dd.z() * (op.z() - pp.z());
        d1.abs() <= linear_tolerance && self.axis.is_parallel(&other.axis, angular_tolerance)
    }

    fn is_coplanar_ax1(&self, axis: &NAx1, linear_tolerance: f64, angular_tolerance: f64) -> bool {
        let dd = self.axis.direction();
        let pp = self.axis.location();
        let ap = axis.location();
        let d1 =
            dd.x() * (ap.x() - pp.x()) + dd.y() * (ap.y() - pp.y()) + dd.z() * (ap.z() - pp.z());
        d1.abs() <= linear_tolerance && self.axis.is_normal(axis, angular_tolerance)
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
        self.axis
            .set_direction(self.vxdir.crossed(&self.vydir).expect("Invalid direction"));
    }

    fn mirrored_ax1(&self, axis: &NAx1) -> Self {
        let mut result = self.clone();
        result.mirror_ax1(axis);
        result
    }

    fn mirror_ax2(&mut self, plane: &Self) {
        self.vydir.mirror_ax2(plane);
        self.vxdir.mirror_ax2(plane);
        self.axis.mirror_ax2(plane);
        self.axis
            .set_direction(self.vxdir.crossed(&self.vydir).expect("Invalid direction"));
    }

    fn mirrored_ax2(&self, plane: &Self) -> Self {
        let mut result = self.clone();
        result.mirror_ax2(plane);
        result
    }

    fn rotate(&mut self, axis: &NAx1, angle: f64) {
        self.axis.rotate(axis, angle);
        self.vxdir.rotate(axis, angle);
        self.vydir.rotate(axis, angle);
        self.axis
            .set_direction(self.vxdir.crossed(&self.vydir).expect("Invalid direction"));
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
        self.axis
            .set_direction(self.vxdir.crossed(&self.vydir).expect("Invalid direction"));
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

    fn translate_pnts(&mut self, from: &NPoint3d, to: &NPoint3d) {
        self.axis.translate_pnts(from, to);
    }

    fn translated_pnts(&self, from: &NPoint3d, to: &NPoint3d) -> Self {
        let mut result = self.clone();
        result.translate_pnts(from, to);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn ax2(loc: (f64, f64, f64), dir: (f64, f64, f64), x_dir: (f64, f64, f64)) -> NAx2 {
        NAx2::new(
            NPoint3d::new(loc.0, loc.1, loc.2),
            NDir::new(dir.0, dir.1, dir.2).expect("Invalid direction"),
            NDir::new(x_dir.0, x_dir.1, x_dir.2).expect("Invalid X direction"),
        )
        .expect("Failed to create NAx2")
    }

    #[test]
    fn test_origin() {
        let origin = NAx2::origin();
        assert_eq!(origin.location(), &NPoint3d::new(0.0, 0.0, 0.0));
        assert_eq!(origin.direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
        assert_eq!(origin.x_direction(), &NDir::new(1.0, 0.0, 0.0).unwrap());
        assert_eq!(origin.y_direction(), &NDir::new(0.0, 1.0, 0.0).unwrap());
    }

    #[test]
    fn test_new_from_point_and_direction() {
        let ax2 = NAx2::new_from_point_and_direction(
            NPoint3d::new(1.0, 2.0, 3.0),
            NDir::new(0.0, 0.0, 1.0).unwrap(),
        )
        .unwrap();
        assert_eq!(ax2.location(), &NPoint3d::new(1.0, 2.0, 3.0));
        assert_eq!(ax2.direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
        assert!(ax2.x_direction().is_normal(ax2.direction(), 1e-5));
        assert!(ax2.y_direction().is_normal(ax2.direction(), 1e-5));
        assert!(ax2.x_direction().is_normal(ax2.y_direction(), 1e-5));
    }

    #[test]
    fn test_set_axis() {
        let mut ax2 = NAx2::origin();
        let new_axis = NAx1::new(
            NPoint3d::new(1.0, 0.0, 0.0),
            NDir::new(1.0, 0.0, 0.0).unwrap(),
        );
        ax2.set_axis(&new_axis).unwrap();
        assert_eq!(ax2.axis(), &new_axis);
        assert!(ax2.x_direction().is_normal(ax2.direction(), 1e-5));
        assert!(ax2.y_direction().is_normal(ax2.direction(), 1e-5));
    }

    #[test]
    fn test_set_direction() {
        let mut ax2 = NAx2::origin();
        let new_dir = NDir::new(1.0, 0.0, 0.0).unwrap();
        ax2.set_direction(new_dir.clone()).unwrap();
        assert_eq!(ax2.direction(), &new_dir);
        assert!(ax2.x_direction().is_normal(&new_dir, 1e-5));
        assert!(ax2.y_direction().is_normal(&new_dir, 1e-5));
    }

    #[test]
    fn test_set_x_direction() {
        let mut ax2 = NAx2::origin();
        let new_x_dir = NDir::new(0.0, 1.0, 0.0).unwrap();
        ax2.set_x_direction(new_x_dir.clone()).unwrap();
        assert_eq!(ax2.x_direction(), &new_x_dir);
        assert!(ax2.y_direction().is_normal(ax2.direction(), 1e-5));
        assert!(ax2.x_direction().is_normal(ax2.y_direction(), 1e-5));
    }

    #[test]
    fn test_angle() {
        let ax2_1 = ax2((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        let ax2_2 = ax2((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0));
        assert!((ax2_1.angle(&ax2_2) - PI / 2.0).abs() < 1e-5);
    }

    #[test]
    fn test_is_coplanar_ax2() {
        let ax2_1 = ax2((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        let ax2_2 = ax2((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        assert!(ax2_1.is_coplanar_ax2(&ax2_2, 1e-5, 1e-5));
    }

    #[test]
    fn test_is_coplanar_ax1() {
        let ax2 = ax2((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        let ax1 = NAx1::new(
            NPoint3d::new(1.0, 0.0, 0.0),
            NDir::new(1.0, 0.0, 0.0).unwrap(),
        );
        assert!(ax2.is_coplanar_ax1(&ax1, 1e-5, 1e-5));
    }

    #[test]
    fn test_mirror_pnt() {
        let mut ax2 = ax2((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        let point = NPoint3d::new(0.0, 0.0, 0.0);
        ax2.mirror_pnt(&point);
        assert_eq!(ax2.location(), &NPoint3d::new(-1.0, 0.0, 0.0));
        assert_eq!(ax2.x_direction(), &NDir::new(-1.0, 0.0, 0.0).unwrap());
        assert_eq!(ax2.y_direction(), &NDir::new(0.0, -1.0, 0.0).unwrap());
    }

    #[test]
    fn test_mirror_ax1() {
        let mut ax2 = ax2((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        let axis = NAx1::new(
            NPoint3d::new(0.0, 0.0, 0.0),
            NDir::new(0.0, 0.0, 1.0).unwrap(),
        );
        ax2.mirror_ax1(&axis);
        assert_eq!(ax2.location(), &NPoint3d::new(1.0, 0.0, 0.0));
        assert!(
            ax2.direction()
                .is_normal(&NDir::new(0.0, 0.0, 1.0).unwrap(), 1e-5)
        );
    }

    #[test]
    fn test_scale() {
        let mut ax2 = ax2((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        let point = NPoint3d::new(0.0, 0.0, 0.0);
        ax2.scale(&point, 2.0);
        assert_eq!(ax2.location(), &NPoint3d::new(2.0, 0.0, 0.0));
        assert_eq!(ax2.direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());

        ax2.scale(&point, -2.0);
        assert_eq!(ax2.location(), &NPoint3d::new(-4.0, 0. sited(&point)));
        assert_eq!(ax2.x_direction(), &NDir::new(-1.0, 0.0, 0.0).unwrap());
        assert_eq!(ax2.y_direction(), &NDir::new(0.0, -1.0, 0.0).unwrap());
    }

    #[test]
    fn test_translate_vec() {
        let mut ax2 = ax2((1.0, 2.0, 3.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0));
        let vec = NVec::new(1.0, 1.0, 1.0);
        ax2.translate_vec(&vec);
        assert_eq!(ax2.location(), &NPoint3d::new(2.0, 3.0, 4.0));
        assert_eq!(ax2.direction(), &NDir::new(0.0, 0.0, 1.0).unwrap());
    }
}
