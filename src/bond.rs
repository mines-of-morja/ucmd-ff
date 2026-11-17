//! Module `bond`: harmonic bond term for UnitcellMD.
//!
//! Computes the energy, forces, and virial for a single harmonic bond
//! `E = 0.5 * k * (r - r0)^2`.
//!
//! Units: positions `r_i` and `r_j` in nm; `r0_nm` in nm;
//! `k_kjmol_nm2` in kJ mol^-1 nm^-2; energy in kJ/mol;
//! forces in kJ mol^-1 nm^-1; virial in kJ/mol.

use crate::errors::UcmdError;
use crate::vec3::{norm, scale, sub};

/// Checks that all supplied scalar inputs are finite.
///
/// Returns `Ok(())` if all values are finite,
/// `Err(UcmdError::NonFiniteCoordinate)` otherwise.
/// Units: none.
fn validate_inputs(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r0_nm: f64,
    k_kjmol_nm2: f64,
) -> Result<(), UcmdError> {
    let values: [f64; 8] = [
        r_i[0],
        r_i[1],
        r_i[2],
        r_j[0],
        r_j[1],
        r_j[2],
        r0_nm,
        k_kjmol_nm2,
    ];
    if values.iter().all(|v| v.is_finite()) {
        Ok(())
    } else {
        Err(UcmdError::NonFiniteCoordinate)
    }
}

/// Computes the energy of a harmonic bond.
/// Units: r_i, r_j in nm; r0_nm in nm; k_kjmol_nm2 in kJ mol^-1 nm^-2; return in kJ/mol.
pub fn bond_energy(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r0_nm: f64,
    k_kjmol_nm2: f64,
) -> Result<f64, UcmdError> {
    validate_inputs(r_i, r_j, r0_nm, k_kjmol_nm2)?;

    let d = sub(r_i, r_j);
    let r = norm(d);

    if !r.is_finite() {
        return Err(UcmdError::NonFiniteCoordinate);
    }

    if r == 0.0 {
        return Err(UcmdError::CoincidentAtoms);
    }

    let dr = r - r0_nm;
    let energy = 0.5 * k_kjmol_nm2 * dr * dr;
    Ok(energy)
}

/// Computes the forces on atoms i and j for a harmonic bond.
/// Units: r_i, r_j in nm; r0_nm in nm; k_kjmol_nm2 in kJ mol^-1 nm^-2; return forces in kJ mol^-1 nm^-1.
pub fn bond_forces(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r0_nm: f64,
    k_kjmol_nm2: f64,
) -> Result<([f64; 3], [f64; 3]), UcmdError> {
    validate_inputs(r_i, r_j, r0_nm, k_kjmol_nm2)?;

    let d = sub(r_i, r_j);
    let r = norm(d);

    if !r.is_finite() {
        return Err(UcmdError::NonFiniteCoordinate);
    }

    if r == 0.0 {
        return Err(UcmdError::CoincidentAtoms);
    }

    let factor = -k_kjmol_nm2 * (r - r0_nm) / r;
    let f_i = scale(d, factor);
    let f_j = scale(f_i, -1.0);

    Ok((f_i, f_j))
}

/// Computes the virial contribution of a harmonic bond.
/// Units: r_i, r_j in nm; r0_nm in nm; k_kjmol_nm2 in kJ mol^-1 nm^-2; return in kJ/mol.
pub fn bond_virial(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r0_nm: f64,
    k_kjmol_nm2: f64,
) -> Result<[[f64; 3]; 3], UcmdError> {
    validate_inputs(r_i, r_j, r0_nm, k_kjmol_nm2)?;

    let d = sub(r_i, r_j);
    let r = norm(d);

    if !r.is_finite() {
        return Err(UcmdError::NonFiniteCoordinate);
    }

    if r == 0.0 {
        return Err(UcmdError::CoincidentAtoms);
    }

    let factor = -k_kjmol_nm2 * (r - r0_nm) / r;
    let f_i = scale(d, factor);

    let mut w = [[0.0_f64; 3]; 3];
    for a in 0..3 {
        w[a] = scale(f_i, d[a]);
    }

    Ok(w)
}
