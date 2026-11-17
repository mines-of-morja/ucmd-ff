//! Computes the minimum-image displacement vector and its Euclidean magnitude
//! for a pair of atoms in a periodic cell.
//!
//! Units: positions in nm; returned displacement vector in nm; returned magnitude in nm.

use crate::cell::{self, Cell};
use crate::errors::UcmdError;
use crate::vec3;

/// Computes the minimum-image displacement vector and its magnitude
/// for a pair of atoms in a periodic cell.
///
/// Units: `r_i`, `r_j` in nm; returned vector in nm; returned scalar in nm.
///
/// The displacement is `d_min = r_i - r_j` corrected by integer lattice
/// translations so that each fractional component lies in [−½, ½].
/// The scalar is the Euclidean norm of `d_min`.
///
/// # Errors
/// * `NonFiniteCoordinate` — any coordinate of `r_i` or `r_j` is NaN or infinite.
/// * `InvalidCell` — the cell fails validation (non-positive diagonal, zero volume).
/// * `CoincidentAtoms` — the minimum-image displacement has zero magnitude.
pub fn minimum_image_displacement(
    r_i: [f64; 3],
    r_j: [f64; 3],
    cell: &Cell,
) -> Result<([f64; 3], f64), UcmdError> {
    let d_min = cell::minimum_image(r_i, r_j, cell)?;
    let dist = vec3::norm(d_min);
    if dist == 0.0 {
        return Err(UcmdError::CoincidentAtoms);
    }
    Ok((d_min, dist))
}
