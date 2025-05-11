use serde::{Deserialize, Serialize};

use crate::gp::{NErrors, NQuaternion, NTrsf, NXYZ};

mod gp {
    pub fn resolution() -> f64 {
        1e-12 // Consistent with trsf.rs
    }
}

/// Trait for linear interpolation between two values of type T.
pub trait Lerp {
    type Output;

    /// Initialize the lerp with start and end values.
    fn init(&mut self, start: Self::Output, end: Self::Output);

    /// Compute the interpolated value at parameter t in [0, 1].
    fn interpolate(&self, t: f64, result: &mut Self::Output) -> Result<(), NErrors>;
}

/// Linear interpolation for 3D coordinates (NXYZ).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NXYZLerp {
    start: NXYZ,
    end: NXYZ,
}

impl Lerp for NXYZLerp {
    type Output = NXYZ;

    fn init(&mut self, start: NXYZ, end: NXYZ) {
        self.start = start;
        self.end = end;
    }

    fn interpolate(&self, t: f64, result: &mut NXYZ) -> Result<(), NErrors> {
        if t < 0.0 || t > 1.0 {
            return Err(NErrors::InvalidParameter);
        }
        result.set_coord(
            self.start.x() + t * (self.end.x() - self.start.x()),
            self.start.y() + t * (self.end.y() - self.start.y()),
            self.start.z() + t * (self.end.z() - self.start.z()),
        );
        Ok(())
    }
}

impl NXYZLerp {
    pub fn new(start: NXYZ, end: NXYZ) -> Self {
        NXYZLerp { start, end }
    }
}

/// Linear interpolation for scale factor (f64).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct F64Lerp {
    start: f64,
    end: f64,
}

impl Lerp for F64Lerp {
    type Output = f64;

    fn init(&mut self, start: f64, end: f64) {
        self.start = start;
        self.end = end;
    }

    fn interpolate(&self, t: f64, result: &mut f64) -> Result<(), NErrors> {
        if t < 0.0 || t > 1.0 {
            return Err(NErrors::InvalidParameter);
        }
        *result = self.start + t * (self.end - self.start);
        Ok(())
    }
}

impl F64Lerp {
    pub fn new(start: f64, end: f64) -> Self {
        F64Lerp { start, end }
    }
}

/// Normalized linear interpolation for quaternions.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NQuaternionNLerp {
    start: NQuaternion,
    end: NQuaternion,
}

impl Lerp for NQuaternionNLerp {
    type Output = NQuaternion;

    fn init(&mut self, start: NQuaternion, end: NQuaternion) {
        self.start = start;
        self.end = end;
    }

    fn interpolate(&self, t: f64, result: &mut NQuaternion) -> Result<(), NErrors> {
        if t < 0.0 || t > 1.0 {
            return Err(NErrors::InvalidParameter);
        }
        // Compute dot product to determine shortest path
        let dot = self.start.dot(&self.end);
        let end_quat = if dot < 0.0 {
            // Negate end quaternion for shortest path
            self.end.negated()
        } else {
            self.end.clone()
        };

        // Linear interpolation
        let mut lerp = NQuaternion::new();
        lerp.set_coords(
            self.start.x() + t * (end_quat.x() - self.start.x()),
            self.start.y() + t * (end_quat.y() - self.start.y()),
            self.start.z() + t * (end_quat.z() - self.start.z()),
            self.start.w() + t * (end_quat.w() - self.start.w()),
        );

        // Normalize the result
        lerp.normalize()?;
        *result = lerp;
        Ok(())
    }
}

impl NQuaternionNLerp {
    pub fn new(start: NQuaternion, end: NQuaternion) -> Self {
        NQuaternionNLerp { start, end }
    }
}

/// Linear interpolation for 3D transformations (NTrsf).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NTrsfNLerp {
    trsf_start: NTrsf,
    trsf_end: NTrsf,
    loc_lerp: NXYZLerp,
    rot_lerp: NQuaternionNLerp,
    scale_lerp: F64Lerp,
}

impl Lerp for NTrsfNLerp {
    type Output = NTrsf;

    fn init(&mut self, start: NTrsf, end: NTrsf) {
        self.trsf_start = start.clone();
        self.trsf_end = end.clone();
        self.loc_lerp
            .init(start.translation_part(), end.translation_part());
        self.rot_lerp.init(start.get_rotation(), end.get_rotation());
        self.scale_lerp
            .init(start.scale_factor(), end.scale_factor());
    }

    fn interpolate(&self, t: f64, result: &mut NTrsf) -> Result<(), NErrors> {
        if t < 0.0 || t > 1.0 {
            return Err(NErrors::InvalidParameter);
        }
        if (t - 0.0).abs() < gp::resolution() {
            *result = self.trsf_start.clone();
            return Ok(());
        }
        if (t - 1.0).abs() < gp::resolution() {
            *result = self.trsf_end.clone();
            return Ok(());
        }

        let mut loc = NXYZ::new(0.0, 0.0, 0.0);
        let mut rot = NQuaternion::new();
        let mut scale = 1.0;

        self.loc_lerp.interpolate(t, &mut loc)?;
        self.rot_lerp.interpolate(t, &mut rot)?;
        self.scale_lerp.interpolate(t, &mut scale)?;

        *result = NTrsf::new();
        result.set_rotation_quat(&rot);
        result.set_translation_part(&loc.into());
        result.set_scale_factor(scale)?;
        Ok(())
    }
}

impl NTrsfNLerp {
    /// Creates a new lerp instance for transformations.
    pub fn new(start: &NTrsf, end: &NTrsf) -> Self {
        let mut lerp = NTrsfNLerp {
            trsf_start: start.clone(),
            trsf_end: end.clone(),
            loc_lerp: NXYZLerp::new(start.translation_part(), end.translation_part()),
            rot_lerp: NQuaternionNLerp::new(start.get_rotation(), end.get_rotation()),
            scale_lerp: F64Lerp::new(start.scale_factor(), end.scale_factor()),
        };
        lerp.init(start.clone(), end.clone());
        lerp
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::{NPoint3d, NVec};

    #[test]
    fn test_nxyz_lerp() {
        let start = NXYZ::new(0.0, 0.0, 0.0);
        let end = NXYZ::new(10.0, 20.0, 30.0);
        let mut lerp = NXYZLerp::new(start, end);
        let mut result = NXYZ::new(0.0, 0.0, 0.0);

        lerp.interpolate(0.0, &mut result).unwrap();
        assert_eq!(result, start);

        lerp.interpolate(1.0, &mut result).unwrap();
        assert_eq!(result, end);

        lerp.interpolate(0.5, &mut result).unwrap();
        assert_eq!(result, NXYZ::new(5.0, 10.0, 15.0));
    }

    #[test]
    fn test_f64_lerp() {
        let start = 1.0;
        let end = 3.0;
        let mut lerp = F64Lerp::new(start, end);
        let mut result = 0.0;

        lerp.interpolate(0.0, &mut result).unwrap();
        assert_eq!(result, 1.0);

        lerp.interpolate(1.0, &mut result).unwrap();
        assert_eq!(result, 3.0);

        lerp.interpolate(0.5, &mut result).unwrap();
        assert_eq!(result, 2.0);
    }

    #[test]
    fn test_trsf_nlerp() {
        let mut start = NTrsf::new();
        let mut end = NTrsf::new();
        start.set_translation_vec(&NVec::new(0.0, 0.0, 0.0));
        end.set_translation_vec(&NVec::new(10.0, 20.0, 30.0));
        start.set_scale_factor(1.0).unwrap();
        end.set_scale_factor(2.0).unwrap();

        let lerp = NTrsfNLerp::new(&start, &end);
        let mut result = NTrsf::new();

        lerp.interpolate(0.0, &mut result).unwrap();
        assert_eq!(result.translation_part(), start.translation_part());
        assert_eq!(result.scale_factor(), start.scale_factor());

        lerp.interpolate(1.0, &mut result).unwrap();
        assert_eq!(result.translation_part(), end.translation_part());
        assert_eq!(result.scale_factor(), end.scale_factor());

        lerp.interpolate(0.5, &mut result).unwrap();
        assert_eq!(result.translation_part(), NXYZ::new(5.0, 10.0, 15.0));
        assert_eq!(result.scale_factor(), 1.5);
    }

    #[test]
    fn test_quaternion_nlerp() {
        let start = NQuaternion::new_from_axis_angle(&NXYZ::new(0.0, 0.0, 1.0), 0.0).unwrap();
        let end =
            NQuaternion::new_from_axis_angle(&NXYZ::new(0.0, 0.0, 1.0), std::f64::consts::PI / 2.0)
                .unwrap();
        let mut lerp = NQuaternionNLerp::new(start.clone(), end.clone());
        let mut result = NQuaternion::new();

        lerp.interpolate(0.0, &mut result).unwrap();
        assert!((result.dot(&start) - 1.0).abs() < gp::resolution());

        lerp.interpolate(1.0, &mut result).unwrap();
        assert!((result.dot(&end) - 1.0).abs() < gp::resolution());

        lerp.interpolate(0.5, &mut result).unwrap();
        let expected =
            NQuaternion::new_from_axis_angle(&NXYZ::new(0.0, 0.0, 1.0), std::f64::consts::PI / 4.0)
                .unwrap();
        assert!((result.dot(&expected) - 1.0).abs() < gp::resolution());
    }
}
