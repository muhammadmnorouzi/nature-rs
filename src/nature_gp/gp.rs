use std::f64;

use super::{
    ax1::{Ax1, NAx1},
    ax2::{Ax2, NAx2},
    ax2d::{Ax2d, NAx2d},
    dir::{Dir, NDir},
    dir2d::{Dir2d, NDir2d},
    point2d::{NPoint2d, Pnt2d},
    point3d::{NPoint3d, Point3d},
};

pub trait GP {
    fn resolution() -> f64;

    /// Identifies a Cartesian point with coordinates X = Y = Z = 0.0.
    fn origin() -> NPoint3d;

    /// Returns a unit vector with the combination (1,0,0).
    fn dx() -> NDir;

    /// Returns a unit vector with the combination (0,1,0).
    fn dy() -> NDir;

    /// Returns a unit vector with the combination (0,0,1).
    fn dz() -> NDir;

    /// Identifies an axis where its origin is Origin and its unit vector coordinates X = 1.0, Y = Z = 0.0.
    fn ox() -> NAx1;

    /// Identifies an axis where its origin is Origin and its unit vector coordinates Y = 1.0, X = Z = 0.0.
    fn oy() -> NAx1;

    /// Identifies an axis where its origin is Origin and its unit vector coordinates Z = 1.0, Y = X = 0.0.
    fn oz() -> NAx1;

    /// Identifies a coordinate system where its origin is Origin,
    /// and its main direction and X direction coordinates Z = 1.0, X = Y = 0.0
    /// and X direction coordinates X = 1.0, Y = Z = 0.0.
    fn xoy() -> NAx2;

    /// Identifies a coordinate system where its origin is Origin,
    /// and its main direction and X direction coordinates Y = 1.0, X = Z = 0.0
    /// and X direction coordinates Z = 1.0, X = Y = 0.0.
    fn zox() -> NAx2;

    /// Identifies a coordinate system where its origin is Origin,
    /// and its main direction and X direction coordinates X = 1.0, Z = Y = 0.0
    /// and X direction coordinates Y = 1.0, X = Z = 0.0.
    fn yoz() -> NAx2;

    /// Identifies a Cartesian point with coordinates X = Y = 0.0 in 2D.
    fn origin2d() -> NPoint2d;

    /// Returns a unit vector with the combinations (1,0) in 2D.
    fn dx2d() -> NDir2d;

    /// Returns a unit vector with the combinations (0,1) in 2D.
    fn dy2d() -> NDir2d;

    /// Identifies an axis where its origin is Origin2d
    /// and its unit vector coordinates are X = 1.0, Y = 0.0.
    fn ox2d() -> NAx2d;

    /// Identifies an axis where its origin is Origin2d
    /// and its unit vector coordinates are Y = 1.0, X = 0.0.
    fn oy2d() -> NAx2d;
}

/// Geometric processor struct, providing access to standard geometric entities.
#[derive(Clone, Copy, Debug)]
pub struct NGP;

impl GP for NGP {
    fn resolution() -> f64 {
        f64::EPSILON
    }

    fn origin() -> NPoint3d {
        NPoint3d::new(0.0, 0.0, 0.0)
    }

    fn dx() -> NDir {
        NDir::new(1.0, 0.0, 0.0).unwrap()
    }

    fn dy() -> NDir {
        NDir::new(0.0, 1.0, 0.0).unwrap()
    }

    fn dz() -> NDir {
        NDir::new(0.0, 0.0, 1.0).unwrap()
    }

    fn ox() -> NAx1 {
        NAx1::new(NGP::origin(), NGP::dx())
    }

    fn oy() -> NAx1 {
        NAx1::new(NGP::origin(), NGP::dy())
    }

    fn oz() -> NAx1 {
        NAx1::new(NGP::origin(), NGP::dz())
    }

    fn xoy() -> NAx2 {
        NAx2::new(NGP::origin(), NGP::dz(), NGP::dx()).unwrap()
    }

    fn zox() -> NAx2 {
        NAx2::new(NGP::origin(), NGP::dy(), NGP::dz()).unwrap()
    }

    fn yoz() -> NAx2 {
        NAx2::new(NGP::origin(), NGP::dx(), NGP::dy()).unwrap()
    }

    fn origin2d() -> NPoint2d {
        NPoint2d::new_with_coords(0.0, 0.0)
    }

    fn dx2d() -> NDir2d {
        NDir2d::new(1.0, 0.0).unwrap()
    }

    fn dy2d() -> NDir2d {
        NDir2d::new(0.0, 1.0).unwrap()
    }

    fn ox2d() -> NAx2d {
        NAx2d::new(NGP::origin2d(), NGP::dx2d())
    }

    fn oy2d() -> NAx2d {
        NAx2d::new(NGP::origin2d(), NGP::dy2d())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resolution() {
        assert_eq!(NGP::resolution(), f64::EPSILON);
    }

    #[test]
    fn test_origin() {
        assert_eq!(NGP::origin(), NPoint3d::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn test_dx() {
        assert_eq!(NGP::dx(), NDir::new(1.0, 0.0, 0.0).unwrap());
    }

    #[test]
    fn test_dy() {
        assert_eq!(NGP::dy(), NDir::new(0.0, 1.0, 0.0).unwrap());
    }

    #[test]
    fn test_dz() {
        assert_eq!(NGP::dz(), NDir::new(0.0, 0.0, 1.0).unwrap());
    }

    #[test]
    fn test_ox() {
        assert_eq!(NGP::ox(), NAx1::new(NGP::origin(), NGP::dx()));
    }

    #[test]
    fn test_oy() {
        assert_eq!(NGP::oy(), NAx1::new(NGP::origin(), NGP::dy()));
    }

    #[test]
    fn test_oz() {
        assert_eq!(NGP::oz(), NAx1::new(NGP::origin(), NGP::dz()));
    }

    #[test]
    fn test_xoy() {
        assert_eq!(
            NGP::xoy(),
            NAx2::new(NGP::origin(), NGP::dz(), NGP::dx()).unwrap()
        );
    }

    #[test]
    fn test_zox() {
        assert_eq!(
            NGP::zox(),
            NAx2::new(NGP::origin(), NGP::dy(), NGP::dz()).unwrap()
        );
    }

    #[test]
    fn test_yoz() {
        assert_eq!(
            NGP::yoz(),
            NAx2::new(NGP::origin(), NGP::dx(), NGP::dy()).unwrap()
        );
    }

    #[test]
    fn test_origin2d() {
        assert_eq!(NGP::origin2d(), NPoint2d::new_with_coords(0.0, 0.0));
    }

    #[test]
    fn test_dx2d() {
        assert_eq!(NGP::dx2d(), NDir2d::new(1.0, 0.0).unwrap());
    }

    #[test]
    fn test_dy2d() {
        assert_eq!(NGP::dy2d(), NDir2d::new(0.0, 1.0).unwrap());
    }

    #[test]
    fn test_ox2d() {
        assert_eq!(NGP::ox2d(), NAx2d::new(NGP::origin2d(), NGP::dx2d()));
    }

    #[test]
    fn test_oy2d() {
        assert_eq!(NGP::oy2d(), NAx2d::new(NGP::origin2d(), NGP::dy2d()));
    }
}
