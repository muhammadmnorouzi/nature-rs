use serde::{Deserialize, Serialize};

/// Enumerates all 24 possible variants of generalized Euler angles,
/// defining general 3D rotation by three rotations around main axes
/// of a coordinate system, in different possible orders.
///
/// The name of each variant corresponds to the order of rotations,
/// prefixed by the type of coordinate system used:
/// - `Intrinsic`: Rotations are made around axes of a rotating
///   coordinate system associated with the object.
/// - `Extrinsic`: Rotations are made around axes of a fixed
///   (static) coordinate system.
///
/// Two specific variants are provided for commonly used conventions:
/// `EulerAngles` (alias for `Intrinsic_ZXZ`) and `YawPitchRoll`
/// (alias for `Intrinsic_ZYX`).
#[derive(Clone, Copy, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub enum NEulerSequence {
    /// Classic Euler angles, alias to `Intrinsic_ZXZ`.
    EulerAngles,

    /// Yaw-Pitch-Roll (or nautical) angles, alias to `Intrinsic_ZYX`.
    YawPitchRoll,

    // Tait-Bryan angles (using three different axes)
    ExtrinsicXYZ,
    ExtrinsicXZY,
    ExtrinsicYZX,
    ExtrinsicYXZ,
    ExtrinsicZXY,
    ExtrinsicZYX,

    IntrinsicXYZ,
    IntrinsicXZY,
    IntrinsicYZX,
    IntrinsicYXZ,
    IntrinsicZXY,
    IntrinsicZYX,

    // Proper Euler angles (using two different axes, first and third the same)
    ExtrinsicXYX,
    ExtrinsicXZX,
    ExtrinsicYZY,
    ExtrinsicYXY,
    ExtrinsicZYZ,
    ExtrinsicZXZ,

    IntrinsicXYX,
    IntrinsicXZX,
    IntrinsicYZY,
    IntrinsicYXY,
    IntrinsicZXZ,
    IntrinsicZYZ,
}

impl NEulerSequence {
    /// Returns true if the sequence is intrinsic (rotations around axes of a rotating coordinate system).
    pub fn is_intrinsic(&self) -> bool {
        matches!(
            self,
            NEulerSequence::EulerAngles
                | NEulerSequence::YawPitchRoll
                | NEulerSequence::IntrinsicXYZ
                | NEulerSequence::IntrinsicXZY
                | NEulerSequence::IntrinsicYZX
                | NEulerSequence::IntrinsicYXZ
                | NEulerSequence::IntrinsicZXY
                | NEulerSequence::IntrinsicZYX
                | NEulerSequence::IntrinsicXYX
                | NEulerSequence::IntrinsicXZX
                | NEulerSequence::IntrinsicYZY
                | NEulerSequence::IntrinsicYXY
                | NEulerSequence::IntrinsicZXZ
                | NEulerSequence::IntrinsicZYZ
        )
    }

    /// Returns true if the sequence is extrinsic (rotations around axes of a fixed coordinate system).
    pub fn is_extrinsic(&self) -> bool {
        !self.is_intrinsic()
    }

    /// Returns the canonical sequence, resolving aliases (`EulerAngles` -> `Intrinsic_ZXZ`, `YawPitchRoll` -> `Intrinsic_ZYX`).
    pub fn canonical(&self) -> Self {
        match self {
            NEulerSequence::EulerAngles => NEulerSequence::IntrinsicZXZ,
            NEulerSequence::YawPitchRoll => NEulerSequence::IntrinsicZYX,
            other => *other,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_intrinsic() {
        assert!(NEulerSequence::EulerAngles.is_intrinsic());
        assert!(NEulerSequence::YawPitchRoll.is_intrinsic());
        assert!(NEulerSequence::IntrinsicZXZ.is_intrinsic());
        assert!(NEulerSequence::IntrinsicZYX.is_intrinsic());
        assert!(!NEulerSequence::ExtrinsicXYZ.is_intrinsic());
        assert!(!NEulerSequence::ExtrinsicZXZ.is_intrinsic());
    }

    #[test]
    fn test_is_extrinsic() {
        assert!(!NEulerSequence::EulerAngles.is_extrinsic());
        assert!(!NEulerSequence::YawPitchRoll.is_extrinsic());
        assert!(!NEulerSequence::IntrinsicZXZ.is_extrinsic());
        assert!(!NEulerSequence::IntrinsicZYX.is_extrinsic());
        assert!(NEulerSequence::ExtrinsicXYZ.is_extrinsic());
        assert!(NEulerSequence::ExtrinsicZXZ.is_extrinsic());
    }

    #[test]
    fn test_canonical() {
        assert_eq!(
            NEulerSequence::EulerAngles.canonical(),
            NEulerSequence::IntrinsicZXZ
        );
        assert_eq!(
            NEulerSequence::YawPitchRoll.canonical(),
            NEulerSequence::IntrinsicZYX
        );
        assert_eq!(
            NEulerSequence::IntrinsicXYZ.canonical(),
            NEulerSequence::IntrinsicXYZ
        );
        assert_eq!(
            NEulerSequence::ExtrinsicZXZ.canonical(),
            NEulerSequence::ExtrinsicZXZ
        );
    }

    #[test]
    fn test_serialization() {
        let seq = NEulerSequence::EulerAngles;
        let serialized = serde_json::to_string(&seq).unwrap();
        assert_eq!(serialized, "\"EulerAngles\"");
        let deserialized: NEulerSequence = serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized, NEulerSequence::EulerAngles);

        let seq = NEulerSequence::IntrinsicZYX;
        let serialized = serde_json::to_string(&seq).unwrap();
        assert_eq!(serialized, "\"Intrinsic_ZYX\"");
        let deserialized: NEulerSequence = serde_json::from_str(&serialized).unwrap();
        assert_eq!(deserialized, NEulerSequence::IntrinsicZYX);
    }
}
