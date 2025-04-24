use std::fmt::Debug;

pub enum NErrors {
    IndexOutOfRange,
    DivisionByZero,
    InvalidDirection,
    ParallelVectors,
}

impl PartialEq for NErrors {
    fn eq(&self, other: &Self) -> bool {
        core::mem::discriminant(self) == core::mem::discriminant(other)
    }
}

impl Debug for NErrors {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NErrors::IndexOutOfRange => write!(f, "IndexOutOfRange"),
            NErrors::DivisionByZero => write!(f, "DivisionByZero"),
            NErrors::InvalidDirection => write!(f, "InvalidDirection"),
            NErrors::ParallelVectors => write!(f, "ParallelVectors"),
        }
    }
}
