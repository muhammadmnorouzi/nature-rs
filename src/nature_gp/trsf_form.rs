use serde::{Deserialize, Serialize};

/// Identifies the type of a geometric transformation.
#[derive(Clone, Copy, PartialEq, Debug, Serialize, Deserialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub enum NTrsfForm {
    /// No transformation (matrix is identity)
    Identity = 0,
    /// Rotation
    Rotation = 1,
    /// Translation
    Translation = 2,
    /// Central symmetry
    PntMirror = 3,
    /// Rotational symmetry
    Ax1Mirror = 4,
    /// Bilateral symmetry
    Ax2Mirror = 5,
    /// Scale
    Scale = 6,
    /// Combination of the above transformations
    CompoundTrsf = 7,
    /// Transformation with not-orthogonal matrix
    Other = 8,
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json;

    #[test]
    fn test_trsf_form_serialization() {
        let form = NTrsfForm::Identity;
        let serialized = serde_json::to_string(&form).unwrap();
        assert_eq!(serialized, "\"IDENTITY\"");

        let form = NTrsfForm::Rotation;
        let serialized = serde_json::to_string(&form).unwrap();
        assert_eq!(serialized, "\"ROTATION\"");
    }

    #[test]
    fn test_trsf_form_deserialization() {
        let json = "\"TRANSLATION\"";
        let form: NTrsfForm = serde_json::from_str(json).unwrap();
        assert_eq!(form, NTrsfForm::Translation);

        let json = "\"COMPOUND_TRSF\"";
        let form: NTrsfForm = serde_json::from_str(json).unwrap();
        assert_eq!(form, NTrsfForm::CompoundTrsf);
    }
}