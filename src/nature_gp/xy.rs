use std::fmt::Debug;

use crate::nature_errors::NErrors;

use super::gp::{GP, NGP};

pub trait XY
where
    Self: Default + Clone + PartialEq + Eq + From<f64> + Debug,
{
    fn cross_square_magnitude(&self, other: &Self) -> f64;
    fn set_linear_form02(&mut self, xyz1: &Self, xyz2: &Self);
    fn set_linear_form12(&mut self, a1: f64, xyz1: &Self, xyz2: &Self);
    fn set_linear_form22(&mut self, a1: f64, xyz1: &Self, a2: f64, xyz2: &Self);
    fn set_linear_form23(&mut self, a1: f64, xyz1: &Self, a2: f64, xyz2: &Self, xyz3: &Self);
    fn crossed(&self, other: &Self) -> f64;
    fn subtracted(&self, other: &Self) -> Self;
    fn multiplied(&self, scalar: f64) -> Self;
    fn divided(&self, scalar: f64) -> Self;
    fn reversed(&mut self) -> Self;
    fn reverse(&mut self);
    fn added(&self, other: &Self) -> Self;
    fn coords(&self, x: &mut f64, y: &mut f64);
    fn change_coord(&mut self, index: usize) -> Result<&mut f64, NErrors>;
    fn coord_by_index(&mut self, index: usize) -> Result<f64, NErrors>;
    fn set_coord_by_index(&mut self, index: usize, scalar: f64) -> Result<f64, NErrors>;
    fn new(x: f64, y: f64) -> Self;
    fn zero() -> Self;
    fn multiply_xy(&mut self, other: &Self);
    fn x(&self) -> f64;
    fn y(&self) -> f64;
    fn set_x(&mut self, x: f64);
    fn set_y(&mut self, y: f64);
    fn set_coords(&mut self, x: f64, y: f64);
    fn add(&mut self, other: &Self);
    fn subtract(&mut self, other: &Self);
    fn multiply(&mut self, scalar: f64);
    fn divide(&mut self, scalar: f64);
    fn dot(&self, other: &Self) -> f64;
    fn cross_magnitude(&self, other: &Self) -> f64;
    fn normalize(&mut self) -> Result<(), NErrors>;
    fn normalized(&self) -> Result<Self, NErrors>;
    fn modulus(&self) -> f64;
    fn square_modulus(&self) -> f64;
    fn is_equal(&self, other: &Self, tolerance: f64) -> bool;
}

pub struct NXY {
    x: f64,
    y: f64,
}

impl Default for NXY {
    fn default() -> Self {
        Self { x: 0.0, y: 0.0 }
    }
}

impl From<f64> for NXY {
    fn from(scalar: f64) -> Self {
        Self::new(scalar, scalar)
    }
}

impl Clone for NXY {
    fn clone(&self) -> Self {
        Self {
            x: self.x,
            y: self.y,
        }
    }
}

impl PartialEq for NXY {
    fn eq(&self, other: &Self) -> bool {
        self.is_equal(other, NGP::resolution())
    }
}

impl Eq for NXY {}

impl Debug for NXY {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NXY")
            .field("x", &self.x)
            .field("y", &self.y)
            .finish()
    }
}

impl XY for NXY {
    fn cross_square_magnitude(&self, other: &Self) -> f64 {
        let a = self.x * other.y - self.y * other.x;
        a * a
    }

    fn set_linear_form02(&mut self, xyz1: &Self, xyz2: &Self) {
        self.x = xyz1.x + xyz2.x;
        self.y = xyz1.y + xyz2.y;
    }

    fn set_linear_form12(&mut self, a1: f64, xyz1: &Self, xyz2: &Self) {
        self.x = a1 * xyz1.x + xyz2.x;
        self.y = a1 * xyz1.y + xyz2.y;
    }

    fn set_linear_form22(&mut self, a1: f64, xyz1: &Self, a2: f64, xyz2: &Self) {
        self.x = a1 * xyz1.x + a2 * xyz2.x;
        self.y = a1 * xyz1.y + a2 * xyz2.y;
    }

    fn set_linear_form23(&mut self, a1: f64, xyz1: &Self, a2: f64, xyz2: &Self, xyz3: &Self) {
        self.x = a1 * xyz1.x + a2 * xyz2.x + xyz3.x;
        self.y = a1 * xyz1.y + a2 * xyz2.y + xyz3.y;
    }

    fn crossed(&self, other: &Self) -> f64 {
        self.x * other.y - self.y * other.x
    }

    fn subtracted(&self, other: &Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }

    fn multiplied(&self, scalar: f64) -> Self {
        Self {
            x: self.x * scalar,
            y: self.y * scalar,
        }
    }

    fn divided(&self, scalar: f64) -> Self {
        Self {
            x: self.x / scalar,
            y: self.y / scalar,
        }
    }

    fn reversed(&mut self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
        }
    }

    fn reverse(&mut self) {
        self.x = -self.x;
        self.y = -self.y;
    }

    fn added(&self, other: &Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }

    fn coords(&self, x: &mut f64, y: &mut f64) {
        *x = self.x;
        *y = self.y;
    }

    fn change_coord(&mut self, index: usize) -> Result<&mut f64, NErrors> {
        match index {
            1 => Ok(&mut self.x),
            2 => Ok(&mut self.y),
            _ => Err(NErrors::IndexOutOfRange),
        }
    }

    fn coord_by_index(&mut self, index: usize) -> Result<f64, NErrors> {
        match index {
            1 => Ok(self.x),
            2 => Ok(self.y),
            _ => Err(NErrors::IndexOutOfRange),
        }
    }

    fn set_coord_by_index(&mut self, index: usize, scalar: f64) -> Result<f64, NErrors> {
        if !(1..=2).contains(&index) {
            return Err(NErrors::IndexOutOfRange);
        }

        match index {
            1 => self.x = scalar,
            2 => self.y = scalar,
            _ => (),
        }

        Ok(scalar)
    }

    fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    fn zero() -> Self {
        Self { x: 0.0, y: 0.0 }
    }

    fn multiply_xy(&mut self, other: &Self) {
        self.x *= other.x;
        self.y *= other.y;
    }

    fn x(&self) -> f64 {
        self.x
    }

    fn y(&self) -> f64 {
        self.y
    }

    fn set_x(&mut self, x: f64) {
        self.x = x;
    }

    fn set_y(&mut self, y: f64) {
        self.y = y;
    }

    fn set_coords(&mut self, x: f64, y: f64) {
        self.x = x;
        self.y = y;
    }

    fn add(&mut self, other: &Self) {
        self.x += other.x;
        self.y += other.y;
    }

    fn subtract(&mut self, other: &Self) {
        self.x -= other.x;
        self.y -= other.y;
    }

    fn multiply(&mut self, scalar: f64) {
        self.x *= scalar;
        self.y *= scalar;
    }

    fn divide(&mut self, scalar: f64) {
        self.x /= scalar;
        self.y /= scalar;
    }

    fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y
    }

    fn cross_magnitude(&self, other: &Self) -> f64 {
        let a = self.x * other.y - self.y * other.x;
        a.abs()
    }

    fn normalize(&mut self) -> Result<(), NErrors> {
        let m = self.modulus();

        if m <= NGP::resolution() {
            return Err(NErrors::DivisionByZero);
        }

        self.x /= m;
        self.y /= m;

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
        })
    }

    fn modulus(&self) -> f64 {
        f64::sqrt(self.square_modulus())
    }

    fn square_modulus(&self) -> f64 {
        self.x.powi(2) + self.y.powi(2)
    }

    fn is_equal(&self, other: &Self, tolerance: f64) -> bool {
        let val_x = self.x - other.x;
        let val_y = self.y - other.y;

        val_x.abs() <= tolerance && val_y.abs() <= tolerance
    }
}

#[cfg(test)]
mod tests_gp_xyz {
    use super::*;

    fn vec(x: f64, y: f64) -> NXY {
        NXY::new(x, y)
    }

    #[test]
    fn test_clone() {
        let a = vec(1.0, 2.0);
        let b = a.clone();
        assert_eq!(a.x, b.x);
        assert_eq!(a.y, b.y);
    }

    #[test]
    fn test_add_subtract() {
        let a = vec(1.0, 2.0);
        let b = vec(4.0, 5.0);
        let mut a_cloned = a.clone();
        let mut b_cloned = b.clone();
        a_cloned.add(&b);
        b_cloned.subtract(&a);

        assert_eq!(a_cloned, vec(5.0, 7.0));
        assert_eq!(b_cloned, vec(3.0, 3.0));
    }

    #[test]
    fn test_added_subtracted() {
        let a = vec(1.0, 2.0);
        let b = vec(4.0, 5.0);
        assert_eq!(a.added(&b), vec(5.0, 7.0));
        assert_eq!(b.subtracted(&a), vec(3.0, 3.0));
    }

    #[test]
    fn test_multiply_divide() {
        let mut a = vec(2.0, 4.0);
        a.multiply(2.0);
        assert_eq!(a, vec(4.0, 8.0));
        a.divide(4.0);
        assert_eq!(a, vec(1.0, 2.0));
    }

    #[test]
    fn test_multiplied_divided() {
        let a = vec(2.0, 4.0);
        assert_eq!(a.multiplied(2.0), vec(4.0, 8.0));
        assert_eq!(a.divided(2.0), vec(1.0, 2.0));
    }

    #[test]
    fn test_reverse_reversed() {
        let mut a = vec(1.0, -2.0);
        let reversed = a.reversed();
        assert_eq!(reversed, vec(-1.0, 2.0));
        a.reverse();
        assert_eq!(a, vec(-1.0, 2.0));
    }

    #[test]
    fn test_crossed_and_dot() {
        let a = vec(1.0, 0.0);
        let b = vec(0.0, 1.0);
        assert_eq!(a.crossed(&b), 1.0);
        assert_eq!(a.dot(&b), 0.0);
        assert_eq!(a.dot(&a), 1.0);
    }

    #[test]
    fn test_normalize_normalize_success() {
        let mut a = vec(3.0, 4.0);
        assert!(a.normalize().is_ok());
        assert_eq!(a, vec(0.6, 0.8));
        assert!((a.modulus() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalize_failure() {
        let mut a = NXY::zero();
        assert_eq!(a.normalize(), Err(NErrors::DivisionByZero));
        assert_eq!(a.normalized(), Err(NErrors::DivisionByZero));
    }

    #[test]
    fn test_coord_access() {
        let mut a = vec(1.0, 2.0);
        assert_eq!(a.coord_by_index(1).unwrap(), 1.0);
        assert_eq!(a.coord_by_index(4), Err(NErrors::IndexOutOfRange));

        *a.change_coord(1).unwrap() = 5.0;
        assert_eq!(a.x, 5.0);
        assert_eq!(a.change_coord(0), Err(NErrors::IndexOutOfRange));
    }

    #[test]
    fn test_set_coord_by_index() {
        let mut a = NXY::zero();
        assert!(a.set_coord_by_index(1, 7.0).is_ok());
        assert_eq!(a.x, 7.0);
        assert_eq!(a.set_coord_by_index(4, 1.0), Err(NErrors::IndexOutOfRange));
    }

    #[test]
    fn test_cross_magnitude() {
        let a = vec(1.0, 0.0);
        let b = vec(0.0, 1.0);
        assert_eq!(a.cross_magnitude(&b), 1.0);
        assert_eq!(a.cross_square_magnitude(&b), 1.0);
    }

    #[test]
    fn test_is_equal() {
        let a = vec(1.0, 2.0);
        let b = vec(1.000001, 2.0);
        assert!(a.is_equal(&b, 1e-3));
        assert!(!a.is_equal(&b, 1e-7));
    }

    #[test]
    fn test_set_linear_forms() {
        let a = vec(1.0, 2.0);
        let b = vec(4.0, 5.0);
        let c = vec(7.0, 8.0);
        let d = vec(1.0, 1.0);

        let mut target = NXY::zero();

        target.set_linear_form23(1.0, &a, 2.0, &b, &c);

        assert_eq!(
            target,
            vec(1.0 * 1.0 + 2.0 * 4.0 + 7.0, 1.0 * 2.0 + 2.0 * 5.0 + 8.0)
        );

        target.set_linear_form22(1.0, &a, 2.0, &b);
        assert_eq!(target, vec(1.0 * 1.0 + 2.0 * 4.0, 1.0 * 2.0 + 2.0 * 5.0));

        target.set_linear_form12(1.0, &a, &b);
        assert_eq!(target, vec(1.0 * 1.0 + 4.0, 1.0 * 2.0 + 5.0));

        target.set_linear_form02(&a, &b);
        assert_eq!(target, vec(1.0 + 4.0, 2.0 + 5.0));
    }

    #[test]
    fn test_coords() {
        let a = vec(1.1, 2.2);
        let (mut x, mut y) = (0.0, 0.0);
        a.coords(&mut x, &mut y);
        assert_eq!(vec(x, y), vec(1.1, 2.2));
    }

    #[test]
    fn test_mutating_operations() {
        let mut a = vec(1.0, 2.0);
        let b = vec(2.0, 3.0);

        a.add(&b);
        assert_eq!(a, vec(3.0, 5.0));

        a.subtract(&b);
        assert_eq!(a, vec(1.0, 2.0));

        a.multiply(2.0);
        assert_eq!(a, vec(2.0, 4.0));

        a.divide(2.0);
        assert_eq!(a, vec(1.0, 2.0,));

        a.multiply_xy(&b);
        assert_eq!(a, vec(2.0, 6.0));
    }

    #[test]
    fn test_setters() {
        let mut a = NXY::zero();
        a.set_x(1.0);
        a.set_y(2.0);
        assert_eq!(a, vec(1.0, 2.0));
        a.set_coords(4.0, 5.0);
        assert_eq!(a, vec(4.0, 5.0));
    }

    #[test]
    fn test_from() {
        let a = NXY::from(5.0);
        assert_eq!(a, vec(5.0, 5.0));
    }

    #[test]
    fn test_new() {
        let a = NXY::new(1.0, 2.0);
        assert_eq!(a, vec(1.0, 2.0));
    }

    #[test]
    fn test_zero() {
        let a = NXY::zero();
        assert_eq!(a, vec(0.0, 0.0));
    }

    #[test]
    fn test_modulus_square_modulus() {
        let a = vec(3.0, 4.0);
        let b = NXY::from(1.0);
        assert_eq!(a.modulus(), 5.0);
        assert_eq!(a.square_modulus(), 25.0);
        assert_eq!(b.modulus(), f64::sqrt(2.0));
        assert_eq!(b.square_modulus(), 2.0);
    }
}
