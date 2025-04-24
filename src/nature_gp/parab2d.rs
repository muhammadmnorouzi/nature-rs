use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::{
    nature_errors::NErrors,
    nature_gp::{NAx2d, NAx22d, NDir2d, NPnt2d, NTrsf2d, NVec2d},
};

// Trait to define the behavior of a parabola in 2D space
pub trait Parab2d {
    fn new() -> Self;
    fn new_with_mirror_axis_focal(
        mirror_axis: &NAx2d,
        focal_length: f64,
        sense: bool,
    ) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn new_with_axes_focal(axes: &NAx22d, focal_length: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn new_with_directrix_focus(directrix: &NAx2d, focus: &NPnt2d, sense: bool) -> Self;
    fn set_focal(&mut self, focal: f64) -> Result<(), NErrors>;
    fn set_location(&mut self, p: &NPnt2d);
    fn set_mirror_axis(&mut self, a: &NAx2d);
    fn set_axis(&mut self, a: &NAx22d);
    fn coefficients(&self) -> (f64, f64, f64, f64, f64, f64);
    fn directrix(&self) -> NAx2d;
    fn focal(&self) -> f64;
    fn focus(&self) -> NPnt2d;
    fn location(&self) -> NPnt2d;
    fn mirror_axis(&self) -> NAx2d;
    fn axis(&self) -> NAx22d;
    fn parameter(&self) -> f64;
    fn reverse(&mut self);
    fn reversed(&self) -> Self
    where
        Self: Sized;
    fn is_direct(&self) -> bool;
    fn mirror_pnt(&mut self, p: &NPnt2d);
    fn mirrored_pnt(&self, p: &NPnt2d) -> Self
    where
        Self: Sized;
    fn mirror_ax2d(&mut self, a: &NAx2d);
    fn mirrored_ax2d(&self, a: &NAx2d) -> Self
    where
        Self: Sized;
    fn rotate(&mut self, p: &NPnt2d, ang: f64);
    fn rotated(&self, p: &NPnt2d, ang: f64) -> Self
    where
        Self: Sized;
    fn scale(&mut self, p: &NPnt2d, s: f64);
    fn scaled(&self, p: &NPnt2d, s: f64) -> Self
    where
        Self: Sized;
    fn transform(&mut self, t: &NTrsf2d);
    fn transformed(&self, t: &NTrsf2d) -> Self
    where
        Self: Sized;
    fn translate_vec(&mut self, v: &NVec2d);
    fn translated_vec(&self, v: &NVec2d) -> Self
    where
        Self: Sized;
    fn translate_pnts(&mut self, p1: &NPnt2d, p2: &NPnt2d);
    fn translated_pnts(&self, p1: &NPnt2d, p2: &NPnt2d) -> Self
    where
        Self: Sized;
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
}

// Struct representing a parabola in 2D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NParab2d {
    pos: NAx22d,
    focal_length: f64,
}

impl Parab2d for NParab2d {
    /// Creates an indefinite parabola with a very large focal length.
    fn new() -> Self {
        NParab2d {
            pos: NAx22d::new(
                &NPnt2d::new(f64::MAX, f64::MAX),
                &NDir2d::new(1.0, 0.0),
                &NDir2d::new(0.0, 1.0),
            )
            .unwrap(),
            focal_length: f64::MAX,
        }
    }

    /// Creates a parabola with its mirror axis, focal length, and sense.
    fn new_with_mirror_axis_focal(
        mirror_axis: &NAx2d,
        focal_length: f64,
        sense: bool,
    ) -> Result<Self, NErrors> {
        if focal_length < 0.0 {
            return Err(NErrors::ConstructionError);
        }
        Ok(NParab2d {
            pos: NAx22d::new_with_axis_sense(mirror_axis, sense),
            focal_length,
        })
    }

    /// Creates a parabola with its axes and focal length.
    fn new_with_axes_focal(axes: &NAx22d, focal_length: f64) -> Result<Self, NErrors> {
        if focal_length < 0.0 {
            return Err(NErrors::ConstructionError);
        }
        Ok(NParab2d {
            pos: axes.clone(),
            focal_length,
        })
    }

    /// Creates a parabola from a directrix and focus point.
    fn new_with_directrix_focus(directrix: &NAx2d, focus: &NPnt2d, sense: bool) -> Self {
        let dir_loc = directrix.location();
        let dir_vec = directrix.direction();
        let f_vec = NVec2d::new_with_points(&dir_loc, focus);
        let origin = NPnt2d::new(
            dir_loc.x() + dir_vec.x() * f_vec.dot(&NVec2d::new_with_dir(&dir_vec)),
            dir_loc.y() + dir_vec.y() * f_vec.dot(&NVec2d::new_with_dir(&dir_vec)),
        );
        let apex = NPnt2d::new(
            0.5 * (origin.x() + focus.x()),
            0.5 * (origin.y() + focus.y()),
        );
        let focal_length = 0.5 * origin.distance(focus);
        let x_dir = if focal_length > 0.0 {
            NDir2d::new(focus.x() - origin.x(), focus.y() - origin.y())
        } else {
            let angle = if sense {
                -std::f64::consts::FRAC_PI_2
            } else {
                std::f64::consts::FRAC_PI_2
            };
            directrix.rotated(&dir_loc, angle).direction()
        };
        let pos = NAx22d::new(&apex, &x_dir, &dir_vec).unwrap();
        NParab2d { pos, focal_length }
    }

    /// Changes the focal distance of the parabola.
    fn set_focal(&mut self, focal: f64) -> Result<(), NErrors> {
        if focal < 0.0 {
            return Err(NErrors::ConstructionError);
        }
        self.focal_length = focal;
        Ok(())
    }

    /// Changes the vertex of the parabola.
    fn set_location(&mut self, p: &NPnt2d) {
        self.pos.set_location(p);
    }

    /// Sets the mirror axis (X-axis) of the parabola.
    fn set_mirror_axis(&mut self, a: &NAx2d) {
        self.pos.set_x_axis(a);
    }

    /// Sets the local coordinate system of the parabola.
    fn set_axis(&mut self, a: &NAx22d) {
        self.pos.set_axis(a);
    }

    /// Computes the coefficients of the implicit equation A*X^2 + B*Y^2 + 2*C*X*Y + 2*D*X + 2*E*Y + F = 0.
    fn coefficients(&self) -> (f64, f64, f64, f64, f64, f64) {
        let p = 2.0 * self.focal_length;
        let t = NTrsf2d::new_transformation(&self.pos.x_axis());
        let t11 = t.value(1, 1).unwrap();
        let t12 = t.value(1, 2).unwrap();
        let t13 = t.value(1, 3).unwrap();
        let t21 = t.value(2, 1).unwrap();
        let t22 = t.value(2, 2).unwrap();
        let t23 = t.value(2, 3).unwrap();
        let a = t21 * t21;
        let b = t22 * t22;
        let c = t21 * t22;
        let d = t21 * t23 - p * t11;
        let e = t22 * t23 - p * t12;
        let f = t23 * t23 - 2.0 * p * t13;
        (a, b, c, d, e, f)
    }

    /// Computes the directrix of the parabola.
    fn directrix(&self) -> NAx2d {
        let p = NPnt2d::new(
            self.pos.location().x() - self.focal_length * self.pos.x_direction().x(),
            self.pos.location().y() - self.focal_length * self.pos.x_direction().y(),
        );
        NAx2d::new(&p, &self.pos.y_direction())
    }

    /// Returns the focal length of the parabola.
    fn focal(&self) -> f64 {
        self.focal_length
    }

    /// Returns the focus of the parabola.
    fn focus(&self) -> NPnt2d {
        NPnt2d::new(
            self.pos.location().x() + self.focal_length * self.pos.x_direction().x(),
            self.pos.location().y() + self.focal_length * self.pos.x_direction().y(),
        )
    }

    /// Returns the vertex of the parabola.
    fn location(&self) -> NPnt2d {
        self.pos.location()
    }

    /// Returns the symmetry axis of the parabola.
    fn mirror_axis(&self) -> NAx2d {
        self.pos.x_axis()
    }

    /// Returns the local coordinate system of the parabola.
    fn axis(&self) -> NAx22d {
        self.pos.clone()
    }

    /// Returns the parameter (distance between focus and directrix).
    fn parameter(&self) -> f64 {
        2.0 * self.focal_length
    }

    /// Reverses the Y-direction of the parabola.
    fn reverse(&mut self) {
        let mut y_dir = self.pos.y_direction();
        y_dir.reverse();
        self.pos
            .set_axis(&NAx22d::new(&self.pos.location(), &self.pos.x_direction(), &y_dir).unwrap());
    }

    /// Returns a parabola with reversed Y-direction.
    fn reversed(&self) -> Self {
        let mut prb = self.clone();
        prb.reverse();
        prb
    }

    /// Returns true if the coordinate system is direct (right-handed).
    fn is_direct(&self) -> bool {
        self.pos.x_direction().crossed(&self.pos.y_direction()) >= 0.0
    }

    /// Mirrors the parabola with respect to a point.
    fn mirror_pnt(&mut self, p: &NPnt2d) {
        self.pos.mirror_pnt(p);
    }

    /// Returns the parabola mirrored with respect to a point.
    fn mirrored_pnt(&self, p: &NPnt2d) -> Self {
        let mut prb = self.clone();
        prb.mirror_pnt(p);
        prb
    }

    /// Mirrors the parabola with respect to an axis.
    fn mirror_ax2d(&mut self, a: &NAx2d) {
        self.pos.mirror_ax2d(a);
    }

    /// Returns the parabola mirrored with respect to an axis.
    fn mirrored_ax2d(&self, a: &NAx2d) -> Self {
        let mut prb = self.clone();
        prb.mirror_ax2d(a);
        prb
    }

    /// Rotates the parabola around a point by an angle.
    fn rotate(&mut self, p: &NPnt2d, ang: f64) {
        self.pos.rotate(p, ang);
    }

    /// Returns the parabola rotated around a point by an angle.
    fn rotated(&self, p: &NPnt2d, ang: f64) -> Self {
        let mut prb = self.clone();
        prb.rotate(p, ang);
        prb
    }

    /// Scales the parabola with respect to a point.
    fn scale(&mut self, p: &NPnt2d, s: f64) {
        self.focal_length *= s.abs();
        self.pos.scale(p, s);
    }

    /// Returns the parabola scaled with respect to a point.
    fn scaled(&self, p: &NPnt2d, s: f64) -> Self {
        let mut prb = self.clone();
        prb.scale(p, s);
        prb
    }

    /// Transforms the parabola with a transformation.
    fn transform(&mut self, t: &NTrsf2d) {
        self.focal_length *= t.scale_factor().abs();
        self.pos.transform(t);
    }

    /// Returns the parabola transformed with a transformation.
    fn transformed(&self, t: &NTrsf2d) -> Self {
        let mut prb = self.clone();
        prb.transform(t);
        prb
    }

    /// Translates the parabola by a vector.
    fn translate_vec(&mut self, v: &NVec2d) {
        self.pos.translate_vec(v);
    }

    /// Returns the parabola translated by a vector.
    fn translated_vec(&self, v: &NVec2d) -> Self {
        let mut prb = self.clone();
        prb.translate_vec(v);
        prb
    }

    /// Translates the parabola from one point to another.
    fn translate_pnts(&mut self, p1: &NPnt2d, p2: &NPnt2d) {
        self.pos.translate_pnts(p1, p2);
    }

    /// Returns the parabola translated from one point to another.
    fn translated_pnts(&self, p1: &NPnt2d, p2: &NPnt2d) -> Self {
        let mut prb = self.clone();
        prb.translate_pnts(p1, p2);
        prb
    }

    /// Dumps the parabola as JSON.
    fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(out, "{} {{", indent).unwrap();
        writeln!(out, "{}   \"type\": \"NParab2d\",", indent).unwrap();
        writeln!(out, "{}   \"focal_length\": {},", indent, self.focal_length).unwrap();
        writeln!(out, "{}   \"position\":", indent).unwrap();
        self.pos.dump_json(out, depth + 2);
        writeln!(out, "{} }}", indent).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_parabola() -> NParab2d {
        let mirror_axis = NAx2d::new(&NPnt2d::new(0.0, 0.0), &NDir2d::new(1.0, 0.0));
        NParab2d::new_with_mirror_axis_focal(&mirror_axis, 1.0, true).unwrap()
    }

    #[test]
    fn test_new() {
        let parab = NParab2d::new();
        assert_eq!(parab.focal_length, f64::MAX);
    }

    #[test]
    fn test_new_with_mirror_axis_focal() {
        let mirror_axis = NAx2d::new(&NPnt2d::new(0.0, 0.0), &NDir2d::new(1.0, 0.0));
        let parab = NParab2d::new_with_mirror_axis_focal(&mirror_axis, 1.0, true).unwrap();
        assert_eq!(parab.focal_length, 1.0);
        assert_eq!(parab.pos.x_axis(), mirror_axis);
        assert!(matches!(
            NParab2d::new_with_mirror_axis_focal(&mirror_axis, -1.0, true),
            Err(NErrors::ConstructionError)
        ));
    }

    #[test]
    fn test_new_with_directrix_focus() {
        let directrix = NAx2d::new(&NPnt2d::new(0.0, -1.0), &NDir2d::new(0.0, 1.0));
        let focus = NPnt2d::new(0.0, 1.0);
        let parab = NParab2d::new_with_directrix_focus(&directrix, &focus, true);
        assert!((parab.focal_length - 1.0).abs() < 1e-9);
        assert_eq!(parab.location(), NPnt2d::new(0.0, 0.0));
        assert_eq!(parab.focus(), NPnt2d::new(0.0, 1.0));
    }

    #[test]
    fn test_set_focal() {
        let mut parab = create_test_parabola();
        parab.set_focal(2.0).unwrap();
        assert_eq!(parab.focal_length, 2.0);
        assert!(matches!(
            parab.set_focal(-1.0),
            Err(NErrors::ConstructionError)
        ));
    }

    #[test]
    fn test_directrix() {
        let parab = create_test_parabola();
        let directrix = parab.directrix();
        assert_eq!(directrix.location(), NPnt2d::new(-1.0, 0.0));
        assert_eq!(directrix.direction(), NDir2d::new(0.0, 1.0));
    }

    #[test]
    fn test_focus() {
        let parab = create_test_parabola();
        let focus = parab.focus();
        assert_eq!(focus, NPnt2d::new(1.0, 0.0));
    }

    #[test]
    fn test_coefficients() {
        let parab = create_test_parabola();
        let (a, b, c, d, e, f) = parab.coefficients();
        assert!((a - 0.0).abs() < 1e-9);
        assert!((b - 1.0).abs() < 1e-9);
        assert!((c - 0.0).abs() < 1e-9);
        assert!((d + 2.0).abs() < 1e-9);
        assert!((e - 0.0).abs() < 1e-9);
        assert!((f - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_reverse() {
        let mut parab = create_test_parabola();
        let y_dir_before = parab.pos.y_direction();
        parab.reverse();
        let y_dir_after = parab.pos.y_direction();
        assert_eq!(y_dir_after, y_dir_before.reversed());
        assert!(parab.is_direct() != create_test_parabola().is_direct());
    }

    #[test]
    fn test_is_direct() {
        let parab = create_test_parabola();
        assert!(parab.is_direct());
        let parab = NParab2d::new_with_mirror_axis_focal(
            &NAx2d::new(&NPnt2d::new(0.0, 0.0), &NDir2d::new(1.0, 0.0)),
            1.0,
            false,
        )
        .unwrap();
        assert!(!parab.is_direct());
    }

    #[test]
    fn test_scale() {
        let mut parab = create_test_parabola();
        parab.scale(&NPnt2d::new(0.0, 0.0), 2.0);
        assert_eq!(parab.focal_length, 2.0);
        parab.scale(&NPnt2d::new(0.0, 0.0), -2.0);
        assert_eq!(parab.focal_length, 2.0); // Abs taken
    }

    #[test]
    fn test_dump_json() {
        let parab = create_test_parabola();
        let mut output = Vec::new();
        parab.dump_json(&mut output, 0);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NParab2d\""));
        assert!(json.contains("\"focal_length\": 1"));
    }
}
