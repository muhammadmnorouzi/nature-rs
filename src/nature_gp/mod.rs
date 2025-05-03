pub mod ax1;
pub mod ax2;
pub mod ax22d;
pub mod ax2d;
pub mod ax3;
pub mod circ;
pub mod circ2d;
pub mod cone;
pub mod cylinder;
pub mod dir;
pub mod dir2d;
pub mod elips;
pub mod elips2d;
pub mod euler_sequence;
pub mod gp;
pub mod gtrsf;
pub mod gtrsf2d;
pub mod hypr;
pub mod hypr2d;
pub mod line;
pub mod line2d;
pub mod mat;
pub mod mat2d;
pub mod parab;
pub mod parab2d;
pub mod plane;
pub mod point2d;
pub mod point3d;
pub mod quaternion;
pub mod quaternion_lerp;
pub mod sphere;
pub mod torus;
pub mod trsf;
pub mod trsf2d;
pub mod trsf_form;
pub mod trsf_nlerp;
pub mod vec;
pub mod vec2d;
pub mod vec2f;
pub mod vec3f;
pub mod xy;
pub mod xyz;

pub mod prelude {
    use super::*;

    pub use gp::{GP, NGP};
    pub use point2d::{NPoint2d, Point2d};
    pub use point3d::{NPoint3d, Point3d};
    pub use xy::{NXY, XY};
    pub use xyz::{NXYZ, XYZ};

    pub use ax2d::{Ax2d, NAx2d};
    pub use dir2d::{Dir2d, NDir2d};
    pub use trsf::{NTrsf, Trsf};
    pub use trsf_form::NTrsfForm;
    pub use trsf2d::{NTrsf2d, Trsf2d};
    pub use vec2d::{NVec2d, Vec2d};
}
