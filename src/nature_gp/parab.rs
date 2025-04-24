use std::io::Write;

use serde::{Deserialize, Serialize};

use crate::{
    gp::{NAx1, NAx2, NDir, NLin, NPnt, NTrsf, NVec},
    nature_errors::NErrors,
};

// Trait to define the behavior of a parabola in 3D space
pub trait Parab {
    fn new() -> Self;
    fn new_with_ax2_focal(a2: &NAx2, focal: f64) -> Result<Self, NErrors>
    where
        Self: Sized;
    fn new_with_directrix_focus(d: &NAx1, f: &NPnt) -> Self;
    fn set_axis(&mut self, a1: &NAx1) -> Result<(), NErrors>;
    fn set_focal(&mut self, focal: f64) -> Result<(), NErrors>;
    fn set_location(&mut self, p: &NPnt);
    fn set_position(&mut self, a2: &NAx2);
    fn axis(&self) -> NAx1;
    fn directrix(&self) -> NAx1;
    fn focal(&self) -> f64;
    fn focus(&self) -> NPnt;
    fn location(&self) -> NPnt;
    fn parameter(&self) -> f64;
    fn position(&self) -> NAx2;
    fn x_axis(&self) -> NAx1;
    fn y_axis(&self) -> NAx1;
    fn mirror_pnt(&mut self, p: &NPnt);
    fn mirrored_pnt(&self, p: &NPnt) -> Self
    where
        Self: Sized;
    fn mirror_ax1(&mut self, a1: &NAx1);
    fn mirrored_ax1(&self, a1: &NAx1) -> Self
    where
        Self: Sized;
    fn mirror_ax2(&mut self, a2: &NAx2);
    fn mirrored_ax2(&self, a2: &NAx2) -> Self
    where
        Self: Sized;
    fn rotate(&mut self, a1: &NAx1, ang: f64);
    fn rotated(&self, a1: &NAx1, ang: f64) -> Self
    where
        Self: Sized;
    fn scale(&mut self, p: &NPnt, s: f64);
    fn scaled(&self, p: &NPnt, s: f64) -> Self
    where
        Self: Sized;
    fn transform(&mut self, t: &NTrsf);
    fn transformed(&self, t: &NTrsf) -> Self
    where
        Self: Sized;
    fn translate_vec(&mut self, v: &NVec);
    fn translated_vec(&self, v: &NVec) -> Self
    where
        Self: Sized;
    fn translate_pnts(&mut self, p1: &NPnt, p2: &NPnt);
    fn translated_pnts(&self, p1: &NPnt, p2: &NPnt) -> Self
    where
        Self: Sized;
    fn dump_json(&self, out: &mut dyn Write, depth: i32);
}

// Struct representing a parabola in 3D space
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct NParab {
    pos: NAx2,
    focal_length: f64,
}

impl Parab for NParab {
    /// Creates an indefinite parabola with a very large focal length.
    fn new() -> Self {
        NParab {
            pos: NAx2::new(
                &NPnt::new(f64::MAX, f64::MAX, f64::MAX),
                &NDir::new(1.0, 0.0, 0.0),
                &NDir::new(0.0, 1.0, 0.0),
            )
            .unwrap(),
            focal_length: f64::MAX,
        }
    }

    /// Creates a parabola with its local coordinate system and focal length.
    fn new_with_ax2_focal(a2: &NAx2, focal: f64) -> Result<Self, NErrors> {
        if focal < 0.0 {
            return Err(NErrors::ConstructionError);
        }
        Ok(NParab {
            pos: a2.clone(),
            focal_length: focal,
        })
    }

    /// Creates a parabola from a directrix and focus point.
    fn new_with_directrix_focus(d: &NAx1, f: &NPnt) -> Self {
        let lin = NLin::new_with_ax1(d);
        let focal_length = lin.distance(f) / 2.0;
        let ax = lin.normal(f).position();
        let ay = lin.position();
        let dd = ax.direction();
        let vertex = NPnt::new(
            f.x() - focal_length * dd.x(),
            f.y() - focal_length * dd.y(),
            f.z() - focal_length * dd.z(),
        );
        let normal_dir = dd.crossed(&ay.direction());
        let pos = NAx2::new(&vertex, &normal_dir, &ax.direction()).unwrap();
        NParab { pos, focal_length }
    }

    /// Sets the main axis of the parabola.
    fn set_axis(&mut self, a1: &NAx1) -> Result<(), NErrors> {
        self.pos.set_axis(a1)
    }

    /// Changes the focal distance of the parabola.
    fn set_focal(&mut self, focal: f64) -> Result<(), NErrors> {
        if focal < 0.0 {
            return Err(NErrors::ConstructionError);
        }
        self.focal_length = focal;
        Ok(())
    }

    /// Changes the location (vertex) of the parabola.
    fn set_location(&mut self, p: &NPnt) {
        self.pos.set_location(p);
    }

    /// Changes the local coordinate system of the parabola.
    fn set_position(&mut self, a2: &NAx2) {
        self.pos = a2.clone();
    }

    /// Returns the main axis of the parabola.
    fn axis(&self) -> NAx1 {
        self.pos.axis()
    }

    /// Computes the directrix of the parabola.
    fn directrix(&self) -> NAx1 {
        let p = self.pos.location();
        let d = self.pos.x_direction();
        let directrix_p = NPnt::new(
            p.x() - self.focal_length * d.x(),
            p.y() - self.focal_length * d.y(),
            p.z() - self.focal_length * d.z(),
        );
        NAx1::new(&directrix_p, &self.pos.y_direction())
    }

    /// Returns the focal length of the parabola.
    fn focal(&self) -> f64 {
        self.focal_length
    }

    /// Computes the focus of the parabola.
    fn focus(&self) -> NPnt {
        let p = self.pos.location();
        let d = self.pos.x_direction();
        NPnt::new(
            p.x() + self.focal_length * d.x(),
            p.y() + self.focal_length * d.y(),
            p.z() + self.focal_length * d.z(),
        )
    }

    /// Returns the vertex of the parabola.
    fn location(&self) -> NPnt {
        self.pos.location()
    }

    /// Computes the parameter of the parabola (twice the focal length).
    fn parameter(&self) -> f64 {
        2.0 * self.focal_length
    }

    /// Returns the local coordinate system of the parabola.
    fn position(&self) -> NAx2 {
        self.pos.clone()
    }

    /// Returns the symmetry axis of the parabola.
    fn x_axis(&self) -> NAx1 {
        NAx1::new(&self.pos.location(), &self.pos.x_direction())
    }

    /// Returns an axis parallel to the directrix.
    fn y_axis(&self) -> NAx1 {
        NAx1::new(&self.pos.location(), &self.pos.y_direction())
    }

    /// Mirrors the parabola with respect to a point.
    fn mirror_pnt(&mut self, p: &NPnt) {
        self.pos.mirror_pnt(p);
    }

    /// Returns the parabola mirrored with respect to a point.
    fn mirrored_pnt(&self, p: &NPnt) -> Self {
        let mut prb = self.clone();
        prb.mirror_pnt(p);
        prb
    }

    /// Mirrors the parabola with respect to an axis.
    fn mirror_ax1(&mut self, a1: &NAx1) {
        self.pos.mirror_ax1(a1);
    }

    /// Returns the parabola mirrored with respect to an axis.
    fn mirrored_ax1(&self, a1: &NAx1) -> Self {
        let mut prb = self.clone();
        prb.mirror_ax1(a1);
        prb
    }

    /// Mirrors the parabola with respect to a plane.
    fn mirror_ax2(&mut self, a2: &NAx2) {
        self.pos.mirror_ax2(a2);
    }

    /// Returns the parabola mirrored with respect to a plane.
    fn mirrored_ax2(&self, a2: &NAx2) -> Self {
        let mut prb = self.clone();
        prb.mirror_ax2(a2);
        prb
    }

    /// Rotates the parabola around an axis by an angle.
    fn rotate(&mut self, a1: &NAx1, ang: f64) {
        self.pos.rotate(a1, ang);
    }

    /// Returns the parabola rotated around an axis by an angle.
    fn rotated(&self, a1: &NAx1, ang: f64) -> Self {
        let mut prb = self.clone();
        prb.rotate(a1, ang);
        prb
    }

    /// Scales the parabola with respect to a point.
    fn scale(&mut self, p: &NPnt, s: f64) {
        self.focal_length *= s.abs();
        self.pos.scale(p, s);
    }

    /// Returns the parabola scaled with respect to a point.
    fn scaled(&self, p: &NPnt, s: f64) -> Self {
        let mut prb = self.clone();
        prb.scale(p, s);
        prb
    }

    /// Transforms the parabola with a transformation.
    fn transform(&mut self, t: &NTrsf) {
        self.focal_length *= t.scale_factor().abs();
        self.pos.transform(t);
    }

    /// Returns the parabola transformed with a transformation.
    fn transformed(&self, t: &NTrsf) -> Self {
        let mut prb = self.clone();
        prb.transform(t);
        prb
    }

    /// Translates the parabola by a vector.
    fn translate_vec(&mut self, v: &NVec) {
        self.pos.translate_vec(v);
    }

    /// Returns the parabola translated by a vector.
    fn translated_vec(&self, v: &NVec) -> Self {
        let mut prb = self.clone();
        prb.translate_vec(v);
        prb
    }

    /// Translates the parabola from one point to another.
    fn translate_pnts(&mut self, p1: &NPnt, p2: &NPnt) {
        self.pos.translate_pnts(p1, p2);
    }

    /// Returns the parabola translated from one point to another.
    fn translated_pnts(&self, p1: &NPnt, p2: &NPnt) -> Self {
        let mut prb = self.clone();
        prb.translate_pnts(p1, p2);
        prb
    }

    /// Dumps the parabola as JSON.
    fn dump_json(&self, out: &mut dyn Write, depth: i32) {
        let indent = " ".repeat((depth * 2) as usize);
        writeln!(out, "{} {{", indent).unwrap();
        writeln!(out, "{}   \"type\": \"NParab\",", indent).unwrap();
        writeln!(out, "{}   \"focal_length\": {},", indent, self.focal_length).unwrap();
        writeln!(out, "{}   \"position\":", indent).unwrap();
        self.pos.dump_json(out, depth + 2);
        writeln!(out, "{} }}", indent).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_parabola() -> NParab {
        let a2 = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0),
            &NDir::new(1.0, 0.0, 0.0),
        )
        .unwrap();
        NParab::new_with_ax2_focal(&a2, 1.0).unwrap()
    }

    #[test]
    fn test_new() {
        let parab = NParab::new();
        assert_eq!(parab.focal_length, f64::MAX);
    }

    #[test]
    fn test_new_with_ax2_focal() {
        let a2 = NAx2::new(
            &NPnt::new(0.0, 0.0, 0.0),
            &NDir::new(0.0, 0.0, 1.0),
            &NDir::new(1.0, 0.0, 0.0),
        )
        .unwrap();
        let parab = NParab::new_with_ax2_focal(&a2, 1.0).unwrap();
        assert_eq!(parab.focal_length, 1.0);
        assert_eq!(parab.pos, a2);
        assert!(matches!(
            NParab::new_with_ax2_focal(&a2, -1.0),
            Err(NErrors::ConstructionError)
        ));
    }

    #[test]
    fn test_new_with_directrix_focus() {
        let d = NAx1::new(&NPnt::new(0.0, -1.0, 0.0), &NDir::new(0.0, 1.0, 0.0));
        let f = NPnt::new(0.0, 1.0, 0.0);
        let parab = NParab::new_with_directrix_focus(&d, &f);
        assert!((parab.focal_length - 1.0).abs() < 1e-9);
        assert_eq!(parab.location(), NPnt::new(0.0, 0.0, 0.0));
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
        assert_eq!(directrix.location(), NPnt::new(-1.0, 0.0, 0.0));
        assert_eq!(directrix.direction(), NDir::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn test_focus() {
        let parab = create_test_parabola();
        let focus = parab.focus();
        assert_eq!(focus, NPnt::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn test_parameter() {
        let parab = create_test_parabola();
        assert_eq!(parab.parameter(), 2.0);
    }

    #[test]
    fn test_scale() {
        let mut parab = create_test_parabola();
        parab.scale(&NPnt::new(0.0, 0.0, 0.0), 2.0);
        assert_eq!(parab.focal_length, 2.0);
        parab.scale(&NPnt::new(0.0, 0.0, 0.0), -2.0);
        assert_eq!(parab.focal_length, 2.0); // Abs taken
    }

    #[test]
    fn test_transform() {
        let mut parab = create_test_parabola();
        let t = NTrsf::new_scale(&NPnt::new(0.0, 0.0, 0.0), 2.0);
        parab.transform(&t);
        assert_eq!(parab.focal_length, 2.0);
    }

    #[test]
    fn test_dump_json() {
        let parab = create_test_parabola();
        let mut output = Vec::new();
        parab.dump_json(&mut output, 0);
        let json = String::from_utf8(output).unwrap();
        assert!(json.contains("\"type\": \"NParab\""));
        assert!(json.contains("\"focal_length\": 1"));
    }
}
