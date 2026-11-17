//! Defines `UcmdError`, the library's invalid-input error type.
//!
//! The type carries no data and no units. It identifies the class of
//! impossible input so callers can propagate errors instead of substituting
//! zero, stale, or default values.

/// Library error for invalid input. Units: none.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UcmdError {
    /// A position coordinate is NaN or infinite. Units: none.
    NonFiniteCoordinate,
    /// Two atoms coincide where a direction is needed. Units: none.
    CoincidentAtoms,
    /// An index is outside its valid range, or two arrays that must have the same length do
    /// not. Units: none.
    IndexOutOfRange,
    /// A periodic cell is invalid. Units: none.
    InvalidCell,
    /// A parameter is outside its valid range or not finite (for example a non-positive Ewald
    /// splitting parameter or cutoff, a negative force constant, an unknown functional form).
    /// Units: none.
    InvalidParameter,
    /// The total charge of a periodic system differs from zero by more than 1e-6 e (R4.1).
    /// Units: none.
    NonNeutralSystem,
    /// In a periodic system two bonded atoms are at least half the smallest perpendicular cell
    /// height apart: the molecule is not whole (spec/system-package.md §7, E_MOLECULE_SPLIT).
    /// Units: none.
    MoleculeSplit,
}

impl std::fmt::Display for UcmdError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let text = match self {
            UcmdError::NonFiniteCoordinate => "a position coordinate is NaN or infinite",
            UcmdError::CoincidentAtoms => "two interacting atoms are at the same position",
            UcmdError::IndexOutOfRange => {
                "an index is out of range, or the number of positions does not match the system"
            }
            UcmdError::InvalidCell => {
                "invalid periodic cell, or a cutoff not below half the smallest cell height"
            }
            UcmdError::InvalidParameter => "a parameter is out of range or not finite",
            UcmdError::NonNeutralSystem => {
                "the total charge of the periodic system differs from zero by more than 1e-6 e"
            }
            UcmdError::MoleculeSplit => {
                "a bond is at least half the smallest cell height long: the molecule is not whole"
            }
        };
        f.write_str(text)
    }
}

impl std::error::Error for UcmdError {}
