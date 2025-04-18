use std::fmt::Debug;

pub enum NErrors {
    IndexOutOfRange,
    DivisionByZero,
}

impl PartialEq for NErrors {
    fn eq(&self, other: &Self) -> bool {
        core::mem::discriminant(self) == core::mem::discriminant(other)
    }
}

impl Debug for NErrors {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IndexOutOfRange => write!(f, "IndexOutOfRange"),
            Self::DivisionByZero => write!(f, "DivisionByZero"),
        }
    }
}
