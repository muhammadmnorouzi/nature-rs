use super::prelude::*;
use crate::nature_common::prelude::*;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

pub trait XYZ
where
    Self: Default + Clone + PartialEq + Eq + From<f64> + Debug,
{
    fn cross_square_magnitude(&self, other: &Self) -> f64;
    fn set_linear_form02(&mut self, xyz1: &Self, xyz2: &Self);
    fn set_linear_form12(&mut self, a1: f64, xyz1: &Self, xyz2: &Self);
    fn set_linear_form22(&mut self, a1: f64, xyz1: &Self, a2: f64, xyz2: &Self);
    fn set_linear_form23(&mut self, a1: f64, xyz1: &Self, a2: f64, xyz2: &Self, xyz3: &Self);
    fn set_linear_form33(
        &mut self,
        a1: f64,
        xyz1: &Self,
        a2: f64,
        xyz2: &Self,
        a3: f64,
        xyz3: &Self,
    );
    #[allow(clippy::too_many_arguments)]
    fn set_linear_form34(
        &mut self,
        a1: f64,
        xyz1: &Self,
        a2: f64,
        xyz2: &Self,
        a3: f64,
        xyz3: &Self,
        xyz4: &Self,
    );
    fn cross_crossed(&self, left: &Self, right: &Self) -> Self;
    fn cross_cross(&mut self, left: &Self, right: &Self);
    fn dot_cross(&self, left: &Self, right: &Self) -> f64;
    fn crossed(&self, other: &Self) -> Self;
    fn subtracted(&self, other: &Self) -> Self;
    fn multiplied(&self, scalar: f64) -> Self;
    fn divided(&self, scalar: f64) -> Self;
    fn reversed(&mut self) -> Self;
    fn reverse(&mut self);
    fn added(&self, other: &Self) -> Self;
    fn coords(&self) -> (f64, f64, f64);
    fn change_coord(&mut self, index: usize) -> Result<&mut f64, NErrors>;
    fn coord_by_index(&self, index: usize) -> Result<f64, NErrors>;
    fn set_coord_by_index(&mut self, index: usize, scalar: f64) -> Result<(), NErrors>;
    fn new(x: f64, y: f64, z: f64) -> Self;
    fn zero() -> Self;
    fn multiply_xyz(&mut self, other: &Self);
    fn x(&self) -> f64;
    fn y(&self) -> f64;
    fn z(&self) -> f64;
    fn set_x(&mut self, x: f64);
    fn set_y(&mut self, y: f64);
    fn set_z(&mut self, z: f64);
    fn set_coords(&mut self, x: f64, y: f64, z: f64);
    fn add(&mut self, other: &Self);
    fn subtract(&mut self, other: &Self);
    fn multiply(&mut self, scalar: f64);
    fn divide(&mut self, scalar: f64);
    fn dot(&self, other: &Self) -> f64;
    fn cross(&mut self, other: &Self);
    fn cross_magnitude(&self, other: &Self) -> f64;
    fn normalize(&mut self) -> Result<(), NErrors>;
    fn normalized(&self) -> Result<Self, NErrors>;
    fn modulus(&self) -> f64;
    fn square_modulus(&self) -> f64;
    fn is_equal(&self, other: &Self, tolerance: f64) -> bool;
}

#[derive(Serialize, Deserialize)]
pub struct NXYZ {
    x: f64,
    y: f64,
    z: f64,
}

impl Default for NXYZ {
    fn default() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
}

impl From<f64> for NXYZ {
    fn from(scalar: f64) -> Self {
        Self::new(scalar, scalar, scalar)
    }
}

impl Clone for NXYZ {
    fn clone(&self) -> Self {
        Self {
            x: self.x,
            y: self.y,
            z: self.z,
        }
    }
}

impl PartialEq for NXYZ {
    fn eq(&self, other: &Self) -> bool {
        self.is_equal(other, NGP::resolution())
    }
}

impl Eq for NXYZ {}

impl Debug for NXYZ {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NXYZ")
            .field("x", &self.x)
            .field("y", &self.y)
            .field("z", &self.z)
            .finish()
    }
}

impl XYZ for NXYZ {
    fn cross_square_magnitude(&self, other: &Self) -> f64 {
        let x_result = self.y * other.z - self.z * other.y;
        let y_result = self.z * other.x - self.x * other.z;
        let z_result = self.x * other.y - self.y * other.x;
        x_result.powi(2) + y_result.powi(2) + z_result.powi(2)
    }

    fn set_linear_form02(&mut self, xyz1: &Self, xyz2: &Self) {
        self.x = xyz1.x + xyz2.x;
        self.y = xyz1.y + xyz2.y;
        self.z = xyz1.z + xyz2.z;
    }

    fn set_linear_form12(&mut self, a1: f64, xyz1: &Self, xyz2: &Self) {
        self.x = a1 * xyz1.x + xyz2.x;
        self.y = a1 * xyz1.y + xyz2.y;
        self.z = a1 * xyz1.z + xyz2.z;
    }

    fn set_linear_form22(&mut self, a1: f64, xyz1: &Self, a2: f64, xyz2: &Self) {
        self.x = a1 * xyz1.x + a2 * xyz2.x;
        self.y = a1 * xyz1.y + a2 * xyz2.y;
        self.z = a1 * xyz1.z + a2 * xyz2.z;
    }

    fn set_linear_form23(&mut self, a1: f64, xyz1: &Self, a2: f64, xyz2: &Self, xyz3: &Self) {
        self.x = a1 * xyz1.x + a2 * xyz2.x + xyz3.x;
        self.y = a1 * xyz1.y + a2 * xyz2.y + xyz3.y;
        self.z = a1 * xyz1.z + a2 * xyz2.z + xyz3.z;
    }

    fn set_linear_form33(
        &mut self,
        a1: f64,
        xyz1: &Self,
        a2: f64,
        xyz2: &Self,
        a3: f64,
        xyz3: &Self,
    ) {
        self.x = a1 * xyz1.x + a2 * xyz2.x + a3 * xyz3.x;
        self.y = a1 * xyz1.y + a2 * xyz2.y + a3 * xyz3.y;
        self.z = a1 * xyz1.z + a2 * xyz2.z + a3 * xyz3.z;
    }

    fn set_linear_form34(
        &mut self,
        a1: f64,
        xyz1: &Self,
        a2: f64,
        xyz2: &Self,
        a3: f64,
        xyz3: &Self,
        xyz4: &Self,
    ) {
        self.x = a1 * xyz1.x + a2 * xyz2.x + a3 * xyz3.x + xyz4.x;
        self.y = a1 * xyz1.y + a2 * xyz2.y + a3 * xyz3.y + xyz4.y;
        self.z = a1 * xyz1.z + a2 * xyz2.z + a3 * xyz3.z + xyz4.z;
    }

    fn cross_crossed(&self, left: &Self, right: &Self) -> Self {
        let mut self_clone = self.clone();
        self_clone.cross_cross(left, right);
        self_clone
    }

    fn cross_cross(&mut self, left: &Self, right: &Self) {
        let x_result = self.y * (left.x * right.y - left.y * right.x)
            - self.z * (left.z * right.x - left.x * right.z);

        let y_result = self.z * (left.y * right.z - left.z * right.y)
            - self.x * (left.x * right.y - left.y * right.x);

        self.z = self.x * (left.z * right.x - left.x * right.z)
            - self.y * (left.y * right.z - left.z * right.y);

        self.x = x_result;
        self.y = y_result;
    }

    fn dot_cross(&self, left: &Self, right: &Self) -> f64 {
        let x_result = left.y * right.z - left.z * right.y;
        let y_result = left.z * right.x - left.x * right.z;
        let z_result = left.x * right.y - left.y * right.x;

        self.x * x_result + self.y * y_result + self.z * z_result
    }

    fn crossed(&self, other: &Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    fn subtracted(&self, other: &Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    fn multiplied(&self, scalar: f64) -> Self {
        Self {
            x: self.x * scalar,
            y: self.y * scalar,
            z: self.z * scalar,
        }
    }

    fn divided(&self, scalar: f64) -> Self {
        Self {
            x: self.x / scalar,
            y: self.y / scalar,
            z: self.z / scalar,
        }
    }

    fn reversed(&mut self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }

    fn reverse(&mut self) {
        self.x = -self.x;
        self.y = -self.y;
        self.z = -self.z;
    }

    fn added(&self, other: &Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }

    fn coords(&self) -> (f64, f64, f64) {
        (self.x(), self.y(), self.z())
    }

    fn change_coord(&mut self, index: usize) -> Result<&mut f64, NErrors> {
        match index {
            1 => Ok(&mut self.x),
            2 => Ok(&mut self.y),
            3 => Ok(&mut self.z),
            _ => Err(NErrors::IndexOutOfRange),
        }
    }

    fn coord_by_index(&self, index: usize) -> Result<f64, NErrors> {
        match index {
            1 => Ok(self.x),
            2 => Ok(self.y),
            3 => Ok(self.z),
            _ => Err(NErrors::IndexOutOfRange),
        }
    }

    fn set_coord_by_index(&mut self, index: usize, scalar: f64) -> Result<(), NErrors> {
        match index {
            1 => self.x = scalar,
            2 => self.y = scalar,
            3 => self.z = scalar,
            _ => return Err(NErrors::IndexOutOfRange),
        }

        Ok(())
    }

    fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    fn zero() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    fn multiply_xyz(&mut self, other: &Self) {
        self.x *= other.x;
        self.y *= other.y;
        self.z *= other.z;
    }

    fn x(&self) -> f64 {
        self.x
    }

    fn y(&self) -> f64 {
        self.y
    }

    fn z(&self) -> f64 {
        self.z
    }

    fn set_x(&mut self, x: f64) {
        self.x = x;
    }

    fn set_y(&mut self, y: f64) {
        self.y = y;
    }

    fn set_z(&mut self, z: f64) {
        self.z = z;
    }

    fn set_coords(&mut self, x: f64, y: f64, z: f64) {
        self.x = x;
        self.y = y;
        self.z = z;
    }

    fn add(&mut self, other: &Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }

    fn subtract(&mut self, other: &Self) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }

    fn multiply(&mut self, scalar: f64) {
        self.x *= scalar;
        self.y *= scalar;
        self.z *= scalar;
    }

    fn divide(&mut self, scalar: f64) {
        self.x /= scalar;
        self.y /= scalar;
        self.z /= scalar;
    }

    fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn cross(&mut self, other: &Self) {
        let x_result = self.y * other.z - self.z * other.y;
        let y_result = self.z * other.x - self.x * other.z;

        self.z = self.x * other.y - self.y * other.x;
        self.x = x_result;
        self.y = y_result;
    }

    fn cross_magnitude(&self, other: &Self) -> f64 {
        let x_result = self.y * other.z - self.z * other.y;
        let y_result = self.z * other.x - self.x * other.z;
        let z_result = self.x * other.y - self.y * other.x;
        f64::sqrt(x_result.powi(2) + y_result.powi(2) + z_result.powi(2))
    }

    fn normalize(&mut self) -> Result<(), NErrors> {
        let m = self.modulus();

        if m <= NGP::resolution() {
            return Err(NErrors::DivisionByZero);
        }

        self.x /= m;
        self.y /= m;
        self.z /= m;

        Ok(())
    }

    fn normalized(&self) -> Result<Self, NErrors> {
        let m = self.modulus();

        if m <= NGP::resolution() {
            return Err(NErrors::DivisionByZero);
        }

        Ok(Self {
            x: self.x / m,
            y: self.y / m,
            z: self.z / m,
        })
    }

    fn modulus(&self) -> f64 {
        f64::sqrt(self.square_modulus())
    }

    fn square_modulus(&self) -> f64 {
        self.x.powi(2) + self.y.powi(2) + self.z.powi(2)
    }

    fn is_equal(&self, other: &Self, tolerance: f64) -> bool {
        let val_x = self.x - other.x;
        let val_y = self.y - other.y;
        let val_z = self.z - other.z;

        val_x.abs() <= tolerance && val_y.abs() <= tolerance && val_z.abs() <= tolerance
    }
}

#[cfg(test)]
mod tests_gp_xyz {
    use super::*;

    fn vec(x: f64, y: f64, z: f64) -> NXYZ {
        NXYZ::new(x, y, z)
    }

    #[test]
    fn test_clone() {
        let a = vec(1.0, 2.0, 3.0);
        let b = a.clone();
        assert_eq!(a.x, b.x);
        assert_eq!(a.y, b.y);
        assert_eq!(a.z, b.z);
    }

    #[test]
    fn test_add_subtract() {
        let a = vec(1.0, 2.0, 3.0);
        let b = vec(4.0, 5.0, 6.0);
        let mut a_cloned = a.clone();
        let mut b_cloned = b.clone();
        a_cloned.add(&b);
        b_cloned.subtract(&a);

        assert_eq!(a_cloned, vec(5.0, 7.0, 9.0));
        assert_eq!(b_cloned, vec(3.0, 3.0, 3.0));
    }

    #[test]
    fn test_added_subtracted() {
        let a = vec(1.0, 2.0, 3.0);
        let b = vec(4.0, 5.0, 6.0);
        assert_eq!(a.added(&b), vec(5.0, 7.0, 9.0));
        assert_eq!(b.subtracted(&a), vec(3.0, 3.0, 3.0));
    }

    #[test]
    fn test_multiply_divide() {
        let mut a = vec(2.0, 4.0, 6.0);
        a.multiply(2.0);
        assert_eq!(a, vec(4.0, 8.0, 12.0));
        a.divide(4.0);
        assert_eq!(a, vec(1.0, 2.0, 3.0));
    }

    #[test]
    fn test_multiplied_divided() {
        let a = vec(2.0, 4.0, 6.0);
        assert_eq!(a.multiplied(2.0), vec(4.0, 8.0, 12.0));
        assert_eq!(a.divided(2.0), vec(1.0, 2.0, 3.0));
    }

    #[test]
    fn test_reverse_reversed() {
        let mut a = vec(1.0, -2.0, 3.0);
        let reversed = a.reversed();
        assert_eq!(reversed, vec(-1.0, 2.0, -3.0));
        a.reverse();
        assert_eq!(a, vec(-1.0, 2.0, -3.0));
    }

    #[test]
    fn test_cross_crossed_and_dot() {
        let a = vec(1.0, 0.0, 0.0);
        let b = vec(0.0, 1.0, 0.0);
        assert_eq!(a.crossed(&b), vec(0.0, 0.0, 1.0));
        let mut a_cloned = a.clone();
        a_cloned.cross(&b);
        assert_eq!(a_cloned, vec(0.0, 0.0, 1.0));
        assert_eq!(a.dot(&b), 0.0);
        assert_eq!(a.dot(&a), 1.0);
    }

    #[test]
    fn test_cross_cross() {
        let a = vec(1.0, 2.0, 3.0);
        let b = vec(4.0, 5.0, 6.0);
        let c = vec(7.0, 8.0, 9.0);
        let mut clone = a.clone();
        clone.cross_cross(&b, &c);
        let value = a.cross_crossed(&b, &c);
        assert_eq!(clone, value);
    }

    #[test]
    fn test_dot_cross() {
        let a = vec(1.0, 2.0, 3.0);
        let b = vec(4.0, 5.0, 6.0);
        let c = vec(7.0, 8.0, 9.0);
        let result = a.dot_cross(&b, &c);
        // Manually calculated result should be zero (b x c is parallel)
        assert_eq!(result, 0.0);
    }

    #[test]
    fn test_normalize_normalize_success() {
        let mut a = vec(3.0, 4.0, 0.0);
        assert!(a.normalize().is_ok());
        assert_eq!(a, vec(0.6, 0.8, 0.0));
        assert!((a.modulus() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalize_failure() {
        let mut a = NXYZ::zero();
        assert_eq!(a.normalize(), Err(NErrors::DivisionByZero));
        assert_eq!(a.normalized(), Err(NErrors::DivisionByZero));
    }

    #[test]
    fn test_coord_access() {
        let mut a = vec(1.0, 2.0, 3.0);
        assert_eq!(a.coord_by_index(1).unwrap(), 1.0);
        assert_eq!(a.coord_by_index(4), Err(NErrors::IndexOutOfRange));

        *a.change_coord(1).unwrap() = 5.0;
        assert_eq!(a.x, 5.0);
        assert_eq!(a.change_coord(0), Err(NErrors::IndexOutOfRange));
    }

    #[test]
    fn test_set_coord_by_index() {
        let mut a = NXYZ::zero();
        assert!(a.set_coord_by_index(1, 7.0).is_ok());
        assert_eq!(a.x, 7.0);
        assert_eq!(a.set_coord_by_index(4, 1.0), Err(NErrors::IndexOutOfRange));
    }

    #[test]
    fn test_cross_magnitude() {
        let a = vec(1.0, 0.0, 0.0);
        let b = vec(0.0, 1.0, 0.0);
        assert_eq!(a.cross_magnitude(&b), 1.0);
        assert_eq!(a.cross_square_magnitude(&b), 1.0);
    }

    #[test]
    fn test_is_equal() {
        let a = vec(1.0, 2.0, 3.0);
        let b = vec(1.000001, 2.0, 3.0);
        assert!(a.is_equal(&b, 1e-3));
        assert!(!a.is_equal(&b, 1e-7));
    }

    #[test]
    fn test_set_linear_forms() {
        let a = vec(1.0, 2.0, 3.0);
        let b = vec(4.0, 5.0, 6.0);
        let c = vec(7.0, 8.0, 9.0);
        let d = vec(1.0, 1.0, 1.0);

        let mut target = NXYZ::zero();
        target.set_linear_form34(1.0, &a, 2.0, &b, 3.0, &c, &d);

        assert_eq!(
            target,
            vec(
                1.0 * 1.0 + 2.0 * 4.0 + 3.0 * 7.0 + 1.0,
                1.0 * 2.0 + 2.0 * 5.0 + 3.0 * 8.0 + 1.0,
                1.0 * 3.0 + 2.0 * 6.0 + 3.0 * 9.0 + 1.0
            )
        );

        target.set_linear_form33(1.0, &a, 2.0, &b, 3.0, &c);

        assert_eq!(
            target,
            vec(
                1.0 * 1.0 + 2.0 * 4.0 + 3.0 * 7.0,
                1.0 * 2.0 + 2.0 * 5.0 + 3.0 * 8.0,
                1.0 * 3.0 + 2.0 * 6.0 + 3.0 * 9.0
            )
        );

        target.set_linear_form23(1.0, &a, 2.0, &b, &c);

        assert_eq!(
            target,
            vec(
                1.0 * 1.0 + 2.0 * 4.0 + 7.0,
                1.0 * 2.0 + 2.0 * 5.0 + 8.0,
                1.0 * 3.0 + 2.0 * 6.0 + 9.0
            )
        );

        target.set_linear_form22(1.0, &a, 2.0, &b);

        assert_eq!(
            target,
            vec(
                1.0 * 1.0 + 2.0 * 4.0,
                1.0 * 2.0 + 2.0 * 5.0,
                1.0 * 3.0 + 2.0 * 6.0
            )
        );

        target.set_linear_form12(1.0, &a, &b);

        assert_eq!(
            target,
            vec(1.0 * 1.0 + 4.0, 1.0 * 2.0 + 5.0, 1.0 * 3.0 + 6.0)
        );

        target.set_linear_form02(&a, &b);

        assert_eq!(target, vec(1.0 + 4.0, 2.0 + 5.0, 3.0 + 6.0));
    }

    #[test]
    fn test_coords() {
        let a = vec(1.1, 2.2, 3.3);
        let (x, y, z) = a.coords();
        assert_eq!(vec(x, y, z), vec(1.1, 2.2, 3.3));
    }

    #[test]
    fn test_mutating_operations() {
        let mut a = vec(1.0, 2.0, 3.0);
        let b = vec(2.0, 3.0, 4.0);

        a.add(&b);
        assert_eq!(a, vec(3.0, 5.0, 7.0));

        a.subtract(&b);
        assert_eq!(a, vec(1.0, 2.0, 3.0));

        a.multiply(2.0);
        assert_eq!(a, vec(2.0, 4.0, 6.0));

        a.divide(2.0);
        assert_eq!(a, vec(1.0, 2.0, 3.0));

        a.multiply_xyz(&b);
        assert_eq!(a, vec(2.0, 6.0, 12.0));
    }

    #[test]
    fn test_setters() {
        let mut a = NXYZ::zero();
        a.set_x(1.0);
        a.set_y(2.0);
        a.set_z(3.0);
        assert_eq!(a, vec(1.0, 2.0, 3.0));
        a.set_coords(4.0, 5.0, 6.0);
        assert_eq!(a, vec(4.0, 5.0, 6.0));
    }

    #[test]
    fn test_from() {
        let a = NXYZ::from(5.0);
        assert_eq!(a, vec(5.0, 5.0, 5.0));
    }

    #[test]
    fn test_new() {
        let a = NXYZ::new(1.0, 2.0, 3.0);
        assert_eq!(a, vec(1.0, 2.0, 3.0));
    }

    #[test]
    fn test_zero() {
        let a = NXYZ::zero();
        assert_eq!(a, vec(0.0, 0.0, 0.0));
    }

    #[test]
    fn test_modulus_square_modulus() {
        let a = vec(3.0, 4.0, 0.0);
        let b = NXYZ::from(1.0);
        assert_eq!(a.modulus(), 5.0);
        assert_eq!(a.square_modulus(), 25.0);
        assert_eq!(b.modulus(), f64::sqrt(3.0));
        assert_eq!(b.square_modulus(), 3.0);
    }
}
