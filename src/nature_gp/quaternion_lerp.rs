use crate::gp::NQuaternion;

// Linear interpolation of quaternions (non-unit result)
#[derive(Clone, PartialEq, Debug)]
pub struct NQuaternionNLerp {
    q_start: NQuaternion,
    q_end: NQuaternion,
}

impl NQuaternionNLerp {
    /// Creates an empty interpolator.
    pub fn new() -> Self {
        NQuaternionNLerp {
            q_start: NQuaternion::new(),
            q_end: NQuaternion::new(),
        }
    }

    /// Creates an interpolator with start and end quaternions.
    pub fn new_with_quaternions(q_start: &NQuaternion, q_end: &NQuaternion) -> Self {
        let mut lerp = NQuaternionNLerp::new();
        lerp.init(q_start, q_end);
        lerp
    }

    /// Initializes with start and end quaternions.
    pub fn init(&mut self, q_start: &NQuaternion, q_end: &NQuaternion) {
        self.init_from_unit(&q_start.normalized(), &q_end.normalized());
    }

    /// Initializes with unit start and end quaternions.
    pub fn init_from_unit(&mut self, q_start: &NQuaternion, q_end: &NQuaternion) {
        self.q_start = q_start.clone();
        self.q_end = q_end.clone();
        if self.q_start.dot(&self.q_end) < 0.0 {
            self.q_end = self.q_end.negated();
        }
        self.q_end = self.q_end.subtracted(&self.q_start);
    }

    /// Computes interpolated quaternion at t (0.0 to 1.0).
    pub fn interpolate(&self, t: f64) -> NQuaternion {
        let mut result = self.q_end.scaled(t);
        result.add(&self.q_start);
        result
    }

    /// Static method to interpolate between two quaternions.
    pub fn interpolate_static(q_start: &NQuaternion, q_end: &NQuaternion, t: f64) -> NQuaternion {
        let lerp = NQuaternionNLerp::new_with_quaternions(q_start, q_end);
        lerp.interpolate(t)
    }
}

// Spherical linear interpolation of quaternions (unit result)
#[derive(Clone, PartialEq, Debug)]
pub struct NQuaternionSLerp {
    q_start: NQuaternion,
    q_end: NQuaternion,
    omega: f64,
}

impl NQuaternionSLerp {
    /// Creates an empty interpolator.
    pub fn new() -> Self {
        NQuaternionSLerp {
            q_start: NQuaternion::new(),
            q_end: NQuaternion::new(),
            omega: 0.0,
        }
    }

    /// Creates an interpolator with start and end quaternions.
    pub fn new_with_quaternions(q_start: &NQuaternion, q_end: &NQuaternion) -> Self {
        let mut lerp = NQuaternionSLerp::new();
        lerp.init(q_start, q_end);
        lerp
    }

    /// Initializes with start and end quaternions.
    pub fn init(&mut self, q_start: &NQuaternion, q_end: &NQuaternion) {
        self.init_from_unit(&q_start.normalized(), &q_end.normalized());
    }

    /// Initializes with unit start and end quaternions.
    pub fn init_from_unit(&mut self, q_start: &NQuaternion, q_end: &NQuaternion) {
        self.q_start = q_start.clone();
        self.q_end = q_end.clone();
        let mut cos_omega = self.q_start.dot(&self.q_end);
        if cos_omega < 0.0 {
            cos_omega = -cos_omega;
            self.q_end = self.q_end.negated();
        }
        if cos_omega > 0.9999 {
            cos_omega = 0.9999;
        }
        self.omega = cos_omega.acos();
        let inv_sin_omega = 1.0 / self.omega.sin();
        self.q_start = self.q_start.scaled(inv_sin_omega);
        self.q_end = self.q_end.scaled(inv_sin_omega);
    }

    /// Computes interpolated quaternion at t (0.0 to 1.0).
    pub fn interpolate(&self, t: f64) -> NQuaternion {
        let q1 = self.q_start.scaled(((1.0 - t) * self.omega).sin());
        let q2 = self.q_end.scaled((t * self.omega).sin());
        q1.added(&q2)
    }

    /// Static method to interpolate between two quaternions.
    pub fn interpolate_static(q_start: &NQuaternion, q_end: &NQuaternion, t: f64) -> NQuaternion {
        let lerp = NQuaternionSLerp::new_with_quaternions(q_start, q_end);
        lerp.interpolate(t)
    }
}