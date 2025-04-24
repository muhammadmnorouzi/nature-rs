use std::f64::consts::{FRAC_1_SQRT_2, PI};
use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::gp::{NEulerSequence, NMat, NVec};

const RESOLUTION: f64 = f64::EPSILON; // Approximation of gp::Resolution()

// Trait to define the behavior of a quaternion for 3D rotations
pub trait Quaternion {
    fn new() -> Self;
    fn new_with_components(x: f64, y: f64, z: f64, w: f64) -> Self;
    fn new_with_vectors(from: &NVec, to: &NVec) -> Self;
    fn new_with_vectors_and_help(from: &NVec, to: &NVec, help: &NVec) -> Self;
    fn new_with_axis_angle(axis: &NVec, angle: f64) -> Self;
    fn new_with_matrix(mat: &NMat) -> Self;
    fn is_equal(&self, other: &Self) -> bool;
    fn set_rotation(&mut self, from: &NVec, to: &NVec);
    fn set_rotation_with_help(&mut self, from: &NVec, to: &NVec, help: &NVec);
    fn set_vector_and_angle(&mut self, axis: &NVec, angle: f64);
    fn get_vector_and_angle(&self, axis: &mut NVec, angle: &mut f64);
    fn set_matrix(&mut self, mat: &NMat);
    fn get_matrix(&self) -> NMat;
    fn set_euler_angles(&mut self, order: NEulerSequence, alpha: f64, beta: f64, gamma: f64);
    fn get_euler_angles(
        &self,
        order: NEulerSequence,
        alpha: &mut f64,
        beta: &mut f64,
        gamma: &mut f64,
    );
    fn set(&mut self, x: f64, y: f64, z: f64, w: f64);
    fn set_quaternion(&mut self, other: &Self);
    fn x(&self) -> f64;
    fn y(&self) -> f64;
    fn z(&self) -> f64;
    fn w(&self) -> f64;
    fn set_ident(&mut self);
    fn reverse(&mut self);
    fn reversed(&self) -> Self
    where
        Self: Sized;
    fn invert(&mut self);
    fn inverted(&self) -> Self
    where
        Self: Sized;
    fn square_norm(&self) -> f64;
    fn norm(&self) -> f64;
    fn scale(&mut self, scale: f64);
    fn scaled(&self, scale: f64) -> Self
    where
        Self: Sized;
    fn stabilize_length(&mut self);
    fn normalize(&mut self);
    fn normalized(&self) -> Self
    where
        Self: Sized;
    fn negated(&self) -> Self
    where
        Self: Sized;
    fn added(&self, other: &Self) -> Self
    where
        Self: Sized;
    fn subtracted(&self, other: &Self) -> Self
    where
        Self: Sized;
    fn multiplied(&self, other: &Self) -> Self
    where
        Self: Sized;
    fn add(&mut self, other: &Self);
    fn subtract(&mut self, other: &Self);
    fn multiply(&mut self, other: &Self);
    fn dot(&self, other: &Self) -> f64;
    fn get_rotation_angle(&self) -> f64;
    fn multiply_vec(&self, vec: &NVec) -> NVec;
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
}

// Struct representing a quaternion for 3D rotations
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NQuaternion {
    x: f64,
    y: f64,
    z: f64,
    w: f64,
}

impl Quaternion for NQuaternion {
    /// Creates an identity quaternion (no rotation).
    fn new() -> Self {
        NQuaternion {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            w: 1.0,
        }
    }

    /// Creates a quaternion from component values.
    fn new_with_components(x: f64, y: f64, z: f64, w: f64) -> Self {
        NQuaternion { x, y, z, w }
    }

    /// Creates a quaternion for shortest-arc rotation from one vector to another.
    fn new_with_vectors(from: &NVec, to: &NVec) -> Self {
        let mut q = NQuaternion::new();
        q.set_rotation(from, to);
        q
    }

    /// Creates a quaternion for shortest-arc rotation with a help vector.
    fn new_with_vectors_and_help(from: &NVec, to: &NVec, help: &NVec) -> Self {
        let mut q = NQuaternion::new();
        q.set_rotation_with_help(from, to, help);
        q
    }

    /// Creates a quaternion from axis and angle.
    fn new_with_axis_angle(axis: &NVec, angle: f64) -> Self {
        let mut q = NQuaternion::new();
        q.set_vector_and_angle(axis, angle);
        q
    }

    /// Creates a quaternion from a rotation matrix.
    fn new_with_matrix(mat: &NMat) -> Self {
        let mut q = NQuaternion::new();
        q.set_matrix(mat);
        q
    }

    /// Checks if two quaternions are equal within resolution.
    fn is_equal(&self, other: &Self) -> bool {
        (self.x - other.x).abs() <= RESOLUTION
            && (self.y - other.y).abs() <= RESOLUTION
            && (self.z - other.z).abs() <= RESOLUTION
            && (self.w - other.w).abs() <= RESOLUTION
    }

    /// Sets quaternion to shortest-arc rotation from one vector to another.
    fn set_rotation(&mut self, from: &NVec, to: &NVec) {
        let cross = from.crossed(to);
        self.set(cross.x(), cross.y(), cross.z(), from.dot(to));
        self.normalize();
        self.w += 1.0;
        if self.w <= RESOLUTION {
            if from.z() * from.z() > from.x() * from.x() {
                self.set(0.0, from.z(), -from.y(), self.w);
            } else {
                self.set(from.y(), -from.x(), 0.0, self.w);
            }
        }
        self.normalize();
    }

    /// Sets quaternion to shortest-arc rotation with a help vector.
    fn set_rotation_with_help(&mut self, from: &NVec, to: &NVec, help: &NVec) {
        let cross = from.crossed(to);
        self.set(cross.x(), cross.y(), cross.z(), from.dot(to));
        self.normalize();
        self.w += 1.0;
        if self.w <= RESOLUTION {
            let axis = from.crossed(help);
            self.set(axis.x(), axis.y(), axis.z(), self.w);
        }
        self.normalize();
    }

    /// Sets quaternion from axis and angle.
    fn set_vector_and_angle(&mut self, axis: &NVec, angle: f64) {
        let axis = axis.normalized();
        let half_angle = 0.5 * angle;
        let sin_a = half_angle.sin();
        self.set(
            axis.x() * sin_a,
            axis.y() * sin_a,
            axis.z() * sin_a,
            half_angle.cos(),
        );
    }

    /// Gets axis and angle of rotation.
    fn get_vector_and_angle(&self, axis: &mut NVec, angle: &mut f64) {
        let vl = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        if vl > RESOLUTION {
            let ivl = 1.0 / vl;
            axis.set_coords(self.x * ivl, self.y * ivl, self.z * ivl);
            *angle = if self.w < 0.0 {
                2.0 * (-vl).atan2(-self.w)
            } else {
                2.0 * vl.atan2(self.w)
            };
        } else {
            axis.set_coords(0.0, 0.0, 1.0);
            *angle = 0.0;
        }
    }

    /// Sets quaternion from a rotation matrix.
    fn set_matrix(&mut self, mat: &NMat) {
        let tr = mat.get(1, 1) + mat.get(2, 2) + mat.get(3, 3);
        if tr > 0.0 {
            self.set(
                mat.get(3, 2) - mat.get(2, 3),
                mat.get(1, 3) - mat.get(3, 1),
                mat.get(2, 1) - mat.get(1, 2),
                tr + 1.0,
            );
            self.scale(0.5 / self.w.sqrt());
        } else if mat.get(1, 1) > mat.get(2, 2) && mat.get(1, 1) > mat.get(3, 3) {
            self.set(
                1.0 + mat.get(1, 1) - mat.get(2, 2) - mat.get(3, 3),
                mat.get(1, 2) + mat.get(2, 1),
                mat.get(1, 3) + mat.get(3, 1),
                mat.get(3, 2) - mat.get(2, 3),
            );
            self.scale(0.5 / self.x.sqrt());
        } else if mat.get(2, 2) > mat.get(3, 3) {
            self.set(
                mat.get(1, 2) + mat.get(2, 1),
                1.0 + mat.get(2, 2) - mat.get(1, 1) - mat.get(3, 3),
                mat.get(2, 3) + mat.get(3, 2),
                mat.get(1, 3) - mat.get(3, 1),
            );
            self.scale(0.5 / self.y.sqrt());
        } else {
            self.set(
                mat.get(1, 3) + mat.get(3, 1),
                mat.get(2, 3) + mat.get(3, 2),
                1.0 + mat.get(3, 3) - mat.get(1, 1) - mat.get(2, 2),
                mat.get(2, 1) - mat.get(1, 2),
            );
            self.scale(0.5 / self.z.sqrt());
        }
    }

    /// Returns the rotation matrix.
    fn get_matrix(&self) -> NMat {
        let s = 2.0 / self.square_norm();
        let x2 = self.x * s;
        let y2 = self.y * s;
        let z2 = self.z * s;
        let xx = self.x * x2;
        let xy = self.x * y2;
        let xz = self.x * z2;
        let yy = self.y * y2;
        let yz = self.y * z2;
        let zz = self.z * z2;
        let wx = self.w * x2;
        let wy = self.w * y2;
        let wz = self.w * z2;

        let mut mat = NMat::new();
        mat.set(1, 1, 1.0 - (yy + zz));
        mat.set(1, 2, xy - wz);
        mat.set(1, 3, xz + wy);
        mat.set(2, 1, xy + wz);
        mat.set(2, 2, 1.0 - (xx + zz));
        mat.set(2, 3, yz - wx);
        mat.set(3, 1, xz - wy);
        mat.set(3, 2, yz + wx);
        mat.set(3, 3, 1.0 - (xx + yy));
        mat
    }

    /// Sets quaternion from Euler angles.
    fn set_euler_angles(&mut self, order: NEulerSequence, alpha: f64, beta: f64, gamma: f64) {
        let o = translate_euler_sequence(order);

        let (a, b, c) = if !o.is_extrinsic {
            (gamma, beta, alpha)
        } else {
            (alpha, beta, gamma)
        };
        let b = if o.is_odd { -b } else { b };

        let ti = 0.5 * a;
        let tj = 0.5 * b;
        let th = 0.5 * c;
        let ci = ti.cos();
        let cj = tj.cos();
        let ch = th.cos();
        let si = ti.sin();
        let sj = tj.sin();
        let sh = th.sin();
        let cc = ci * ch;
        let cs = ci * sh;
        let sc = si * ch;
        let ss = si * sh;

        let mut values = [0.0; 4]; // [w, x, y, z]
        if o.is_two_axes {
            values[o.i] = cj * (cs + sc);
            values[o.j] = sj * (cc + ss);
            values[o.k] = sj * (cs - sc);
            values[0] = cj * (cc - ss);
        } else {
            values[o.i] = cj * sc - sj * cs;
            values[o.j] = cj * ss + sj * cc;
            values[o.k] = cj * cs - sj * sc;
            values[0] = cj * cc + sj * ss;
        }
        if o.is_odd {
            values[o.j] = -values[o.j];
        }

        self.x = values[1];
        self.y = values[2];
        self.z = values[3];
        self.w = values[0];
    }

    /// Gets Euler angles describing the rotation.
    fn get_euler_angles(
        &self,
        order: NEulerSequence,
        alpha: &mut f64,
        beta: &mut f64,
        gamma: &mut f64,
    ) {
        let mat = self.get_matrix();
        let o = translate_euler_sequence(order);

        if o.is_two_axes {
            let sy = (mat.get(o.i, o.j).powi(2) + mat.get(o.i, o.k).powi(2)).sqrt();
            if sy > 16.0 * f64::EPSILON {
                *alpha = mat.get(o.i, o.j).atan2(mat.get(o.i, o.k));
                *gamma = mat.get(o.j, o.i).atan2(-mat.get(o.k, o.i));
            } else {
                *alpha = (-mat.get(o.j, o.k)).atan2(mat.get(o.j, o.j));
                *gamma = 0.0;
            }
            *beta = sy.atan2(mat.get(o.i, o.i));
        } else {
            let cy = (mat.get(o.i, o.i).powi(2) + mat.get(o.j, o.i).powi(2)).sqrt();
            if cy > 16.0 * f64::EPSILON {
                *alpha = mat.get(o.k, o.j).atan2(mat.get(o.k, o.k));
                *gamma = mat.get(o.j, o.i).atan2(mat.get(o.i, o.i));
            } else {
                *alpha = (-mat.get(o.j, o.k)).atan2(mat.get(o.j, o.j));
                *gamma = 0.0;
            }
            *beta = (-mat.get(o.k, o.i)).atan2(cy);
        }

        if o.is_odd {
            *alpha = -(*alpha);
            *beta = -(*beta);
            *gamma = -(*gamma);
        }
        if !o.is_extrinsic {
            let temp = *alpha;
            *alpha = *gamma;
            *gamma = temp;
        }
    }

    /// Sets quaternion components.
    fn set(&mut self, x: f64, y: f64, z: f64, w: f64) {
        self.x = x;
        self.y = y;
        self.z = z;
        self.w = w;
    }

    /// Sets quaternion from another quaternion.
    fn set_quaternion(&mut self, other: &Self) {
        self.x = other.x;
        self.y = other.y;
        self.z = other.z;
        self.w = other.w;
    }

    /// Returns the x component.
    fn x(&self) -> f64 {
        self.x
    }

    /// Returns the y component.
    fn y(&self) -> f64 {
        self.y
    }

    /// Returns the z component.
    fn z(&self) -> f64 {
        self.z
    }

    /// Returns the w component.
    fn w(&self) -> f64 {
        self.w
    }

    /// Sets to identity quaternion.
    fn set_ident(&mut self) {
        self.x = 0.0;
        self.y = 0.0;
        self.z = 0.0;
        self.w = 1.0;
    }

    /// Reverses the rotation direction (conjugate).
    fn reverse(&mut self) {
        self.x = -self.x;
        self.y = -self.y;
        self.z = -self.z;
    }

    /// Returns a quaternion with reversed rotation direction.
    fn reversed(&self) -> Self {
        NQuaternion {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: self.w,
        }
    }

    /// Inverts the quaternion (both direction and norm).
    fn invert(&mut self) {
        let inv_norm = 1.0 / self.square_norm();
        self.x = -self.x * inv_norm;
        self.y = -self.y * inv_norm;
        self.z = -self.z * inv_norm;
        self.w *= inv_norm;
    }

    /// Returns the inverted quaternion.
    fn inverted(&self) -> Self {
        let inv_norm = 1.0 / self.square_norm();
        NQuaternion {
            x: -self.x * inv_norm,
            y: -self.y * inv_norm,
            z: -self.z * inv_norm,
            w: self.w * inv_norm,
        }
    }

    /// Returns the square norm.
    fn square_norm(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z + self.w * self.w
    }

    /// Returns the norm.
    fn norm(&self) -> f64 {
        self.square_norm().sqrt()
    }

    /// Scales all components.
    fn scale(&mut self, scale: f64) {
        self.x *= scale;
        self.y *= scale;
        self.z *= scale;
        self.w *= scale;
    }

    /// Returns a scaled quaternion.
    fn scaled(&self, scale: f64) -> Self {
        NQuaternion {
            x: self.x * scale,
            y: self.y * scale,
            z: self.z * scale,
            w: self.w * scale,
        }
    }

    /// Stabilizes quaternion length to avoid zero or infinity.
    fn stabilize_length(&mut self) {
        let cs = self.x.abs() + self.y.abs() + self.z.abs() + self.w.abs();
        if cs > 0.0 {
            self.x /= cs;
            self.y /= cs;
            self.z /= cs;
            self.w /= cs;
        } else {
            self.set_ident();
        }
    }

    /// Normalizes the quaternion to unit length.
    fn normalize(&mut self) {
        let mut magn = self.norm();
        if magn < RESOLUTION {
            self.stabilize_length();
            magn = self.norm();
        }
        self.scale(1.0 / magn);
    }

    /// Returns a normalized quaternion.
    fn normalized(&self) -> Self {
        let mut q = self.clone();
        q.normalize();
        q
    }

    /// Returns a quaternion with negated components.
    fn negated(&self) -> Self {
        NQuaternion {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: -self.w,
        }
    }

    /// Returns the sum of two quaternions (rotation mix).
    fn added(&self, other: &Self) -> Self {
        NQuaternion {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
            w: self.w + other.w,
        }
    }

    /// Returns the difference of two quaternions (rotation mix).
    fn subtracted(&self, other: &Self) -> Self {
        NQuaternion {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
            w: self.w - other.w,
        }
    }

    /// Returns the product of two quaternions (rotation composition).
    fn multiplied(&self, other: &Self) -> Self {
        NQuaternion {
            x: self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y,
            y: self.w * other.y + self.y * other.w + self.z * other.x - self.x * other.z,
            z: self.w * other.z + self.z * other.w + self.x * other.y - self.y * other.x,
            w: self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z,
        }
    }

    /// Adds another quaternion (rotation mix).
    fn add(&mut self, other: &Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
        self.w += other.w;
    }

    /// Subtracts another quaternion (rotation mix).
    fn subtract(&mut self, other: &Self) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
        self.w -= other.w;
    }

    /// Multiplies by another quaternion (rotation composition).
    fn multiply(&mut self, other: &Self) {
        *self = self.multiplied(other);
    }

    /// Computes the dot product.
    fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z + self.w * other.w
    }

    /// Returns the rotation angle from -PI to PI.
    fn get_rotation_angle(&self) -> f64 {
        if self.w < 0.0 {
            2.0 * (-(self.x * self.x + self.y * self.y + self.z * self.z).sqrt()).atan2(-self.w)
        } else {
            2.0 * (self.x * self.x + self.y * self.y + self.z * self.z)
                .sqrt()
                .atan2(self.w)
        }
    }

    /// Rotates a vector by the quaternion.
    fn multiply_vec(&self, vec: &NVec) -> NVec {
        let q = NQuaternion {
            x: vec.x() * self.w + vec.z() * self.y - vec.y() * self.z,
            y: vec.y() * self.w + vec.x() * self.z - vec.z() * self.x,
            z: vec.z() * self.w + vec.y() * self.x - vec.x() * self.y,
            w: vec.x() * self.x + vec.y() * self.y + vec.z() * self.z,
        };
        let inv_norm = 1.0 / self.square_norm();
        NVec::new_with_coords(
            (self.w * q.x + self.x * q.w + self.y * q.z - self.z * q.y) * inv_norm,
            (self.w * q.y + self.y * q.w + self.z * q.x - self.x * q.z) * inv_norm,
            (self.w * q.z + self.z * q.w + self.x * q.y - self.y * q.x) * inv_norm,
        )
    }

    /// Dumps the quaternion as JSON.
    fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(out, "{} {{", indent).unwrap();
        writeln!(out, "{}   \"type\": \"NQuaternion\",", indent).unwrap();
        writeln!(
            out,
            "{}   \"components\": [{}, {}, {}, {}]",
            indent, self.x, self.y, self.z, self.w
        )
        .unwrap();
        writeln!(out, "{} }}", indent).unwrap();
    }
}

// Helper struct for Euler sequence parameters
#[derive(Clone, Copy)]
struct EulerSequenceParams {
    i: usize,
    j: usize,
    k: usize,
    is_odd: bool,
    is_two_axes: bool,
    is_extrinsic: bool,
}

impl EulerSequenceParams {
    fn new(ax1: usize, is_odd: bool, is_two_axes: bool, is_extrinsic: bool) -> Self {
        EulerSequenceParams {
            i: ax1,
            j: 1 + (ax1 + if is_odd { 1 } else { 0 }) % 3,
            k: 1 + (ax1 + if is_odd { 0 } else { 1 }) % 3,
            is_odd,
            is_two_axes,
            is_extrinsic,
        }
    }
}

// Translates Euler sequence to parameters
fn translate_euler_sequence(seq: NEulerSequence) -> EulerSequenceParams {
    use NEulerSequence::*;
    match seq {
        Extrinsic_XYZ => EulerSequenceParams::new(1, false, false, true),
        Extrinsic_XZY => EulerSequenceParams::new(1, true, false, true),
        Extrinsic_YZX => EulerSequenceParams::new(2, false, false, true),
        Extrinsic_YXZ => EulerSequenceParams::new(2, true, false, true),
        Extrinsic_ZXY => EulerSequenceParams::new(3, false, false, true),
        Extrinsic_ZYX => EulerSequenceParams::new(3, true, false, true),
        Intrinsic_XYZ => EulerSequenceParams::new(3, true, false, false),
        Intrinsic_XZY => EulerSequenceParams::new(2, false, false, false),
        Intrinsic_YZX => EulerSequenceParams::new(1, true, false, false),
        Intrinsic_YXZ => EulerSequenceParams::new(3, false, false, false),
        Intrinsic_ZXY => EulerSequenceParams::new(2, true, false, false),
        Intrinsic_ZYX => EulerSequenceParams::new(1, false, false, false),
        Extrinsic_XYX => EulerSequenceParams::new(1, false, true, true),
        Extrinsic_XZX => EulerSequenceParams::new(1, true, true, true),
        Extrinsic_YZY => EulerSequenceParams::new(2, false, true, true),
        Extrinsic_YXY => EulerSequenceParams::new(2, true, true, true),
        Extrinsic_ZXZ => EulerSequenceParams::new(3, false, true, true),
        Extrinsic_ZYZ => EulerSequenceParams::new(3, true, true, true),
        Intrinsic_XYX => EulerSequenceParams::new(1, false, true, false),
        Intrinsic_XZX => EulerSequenceParams::new(1, true, true, false),
        Intrinsic_YZY => EulerSequenceParams::new(2, false, true, false),
        Intrinsic_YXY => EulerSequenceParams::new(2, true, true, false),
        Intrinsic_ZXZ => EulerSequenceParams::new(3, false, true, false),
        Intrinsic_ZYZ => EulerSequenceParams::new(3, true, true, false),
        EulerAngles => EulerSequenceParams::new(3, false, true, false),
        YawPitchRoll => EulerSequenceParams::new(1, false, false, false),
    }
}

// Operator overloads
impl std::ops::Mul<f64> for NQuaternion {
    type Output = Self;
    fn mul(self, scale: f64) -> Self {
        self.scaled(scale)
    }
}

impl std::ops::Neg for NQuaternion {
    type Output = Self;
    fn neg(self) -> Self {
        self.negated()
    }
}

impl std::ops::Add for NQuaternion {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        self.added(&other)
    }
}

impl std::ops::Sub for NQuaternion {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self.subtracted(&other)
    }
}

impl std::ops::Mul for NQuaternion {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        self.multiplied(&other)
    }
}

impl std::ops::Mul<NVec> for NQuaternion {
    type Output = NVec;
    fn mul(self, vec: NVec) -> NVec {
        self.multiply_vec(&vec)
    }
}

impl std::ops::MulAssign<f64> for NQuaternion {
    fn mul_assign(&mut self, scale: f64) {
        self.scale(scale);
    }
}

impl std::ops::AddAssign for NQuaternion {
    fn add_assign(&mut self, other: Self) {
        self.add(&other);
    }
}

impl std::ops::SubAssign for NQuaternion {
    fn sub_assign(&mut self, other: Self) {
        self.subtract(&other);
    }
}

impl std::ops::MulAssign for NQuaternion {
    fn mul_assign(&mut self, other: Self) {
        self.multiply(&other);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_quaternion() -> NQuaternion {
        NQuaternion::new_with_components(0.5, 0.5, 0.5, 0.5)
    }

    #[test]
    fn test_new() {
        let q = NQuaternion::new();
        assert_eq!(q.x(), 0.0);
        assert_eq!(q.w(), 1.0);
    }

    #[test]
    fn test_is_equal() {
        let q1 = create_test_quaternion();
        let q2 = NQuaternion::new_with_components(0.5 + RESOLUTION / 2.0, 0.5, 0.5, 0.5);
        assert!(q1.is_equal(&q2));
    }

    #[test]
    fn test_normalize() {
        let mut q = create_test_quaternion();
        q.normalize();
        assert!((q.norm() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_set_vector_and_angle() {
        let mut q = NQuaternion::new();
        let axis = NVec::new_with_coords(0.0, 0.0, 1.0);
        q.set_vector_and_angle(&axis, PI);
        assert!((q.x().abs() - 0.0).abs() < 1e-9);
        assert!((q.w().abs() - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_get_vector_and_angle() {
        let q = NQuaternion::new_with_components(0.0, 0.0, FRAC_1_SQRT_2, FRAC_1_SQRT_2);
        let mut axis = NVec::new();
        let mut angle = 0.0;
        q.get_vector_and_angle(&mut axis, &mut angle);
        assert!((axis.z() - 1.0).abs() < 1e-9);
        assert!((angle - PI / 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_matrix_conversion() {
        let q = NQuaternion::new_with_components(0.0, 0.0, FRAC_1_SQRT_2, FRAC_1_SQRT_2);
        let mat = q.get_matrix();
        let mut q2 = NQuaternion::new();
        q2.set_matrix(&mat);
        assert!(q.is_equal(&q2));
    }

    #[test]
    fn test_multiply_vec() {
        let q = NQuaternion::new_with_components(0.0, 0.0, FRAC_1_SQRT_2, FRAC_1_SQRT_2);
        let vec = NVec::new_with_coords(1.0, 0.0, 0.0);
        let rotated = q.multiply_vec(&vec);
        assert!((rotated.y() - 1.0).abs() < 1e-9);
        assert!((rotated.x().abs() - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_operators() {
        let q1 = create_test_quaternion();
        let q2 = q1.scaled(2.0);
        assert!((q2.x() - 1.0).abs() < 1e-9);
        let q3 = q1 + q1;
        assert!(q3.is_equal(&q2));
        let q4 = q3 - q1;
        assert!(q4.is_equal(&q1));
    }

    #[test]
    fn test_dump_json() {
        let q = create_test_quaternion();
        let mut output = Vec::new();
        q.dump_json(&mut output, 0);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NQuaternion\""));
        assert!(json.contains("[0.5, 0.5, 0.5, 0.5]"));
    }
}
