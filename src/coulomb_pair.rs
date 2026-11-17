//! Coulomb pair interaction: energy, forces, and virial for a single pair of point charges.
//!
//! Units: charges in e, positions in nm, energy in kJ/mol, forces in kJ mol⁻¹ nm⁻¹,
//! virial in kJ/mol.

use crate::constants::COULOMB_KJ_MOL_NM_E2;
use crate::errors::UcmdError;
use crate::vec3;

/// Computes the Coulomb pair energy.
/// Units: q_i, q_j in e; r_i, r_j in nm; return in kJ/mol.
///
/// E = k_e * q_i * q_j / r, where r = |r_i - r_j| and
/// k_e = COULOMB_KJ_MOL_NM_E2 (kJ mol⁻¹ nm e⁻²).
///
/// # Errors
/// * `NonFiniteCoordinate` — any coordinate of r_i or r_j is NaN or infinite.
/// * `CoincidentAtoms` — r_i and r_j are identical (r = 0).
pub fn coulomb_pair_energy(
    q_i: f64,
    q_j: f64,
    r_i: [f64; 3],
    r_j: [f64; 3],
) -> Result<f64, UcmdError> {
    check_finite(r_i)?;
    check_finite(r_j)?;
    let r_ij = vec3::sub(r_i, r_j);
    let r = vec3::norm(r_ij);
    if r == 0.0 {
        return Err(UcmdError::CoincidentAtoms);
    }
    Ok(COULOMB_KJ_MOL_NM_E2 * q_i * q_j / r)
}

/// Computes the forces on atoms i and j for a Coulomb pair.
/// Units: q_i, q_j in e; r_i, r_j in nm; return forces in kJ mol⁻¹ nm⁻¹.
///
/// F_i =  k_e * q_i * q_j * r_ij / r³
/// F_j = -k_e * q_i * q_j * r_ij / r³
/// where r_ij = r_i - r_j and r = |r_ij|.
///
/// # Errors
/// * `NonFiniteCoordinate` — any coordinate of r_i or r_j is NaN or infinite.
/// * `CoincidentAtoms` — r_i and r_j are identical (r = 0).
pub fn coulomb_pair_forces(
    q_i: f64,
    q_j: f64,
    r_i: [f64; 3],
    r_j: [f64; 3],
) -> Result<([f64; 3], [f64; 3]), UcmdError> {
    check_finite(r_i)?;
    check_finite(r_j)?;
    let r_ij = vec3::sub(r_i, r_j);
    let r = vec3::norm(r_ij);
    if r == 0.0 {
        return Err(UcmdError::CoincidentAtoms);
    }
    let factor = COULOMB_KJ_MOL_NM_E2 * q_i * q_j / (r * r * r);
    let f_i = vec3::scale(r_ij, factor);
    let f_j = vec3::scale(r_ij, -factor);
    Ok((f_i, f_j))
}

/// Computes the 3×3 virial matrix for a Coulomb pair.
/// Units: q_i, q_j in e; r_i, r_j in nm; return in kJ/mol.
///
/// W = r_ij ⊗ F_ij (outer product), where r_ij = r_i - r_j and
/// F_ij is the force on i due to j.
///
/// # Errors
/// * `NonFiniteCoordinate` — any coordinate of r_i or r_j is NaN or infinite.
/// * `CoincidentAtoms` — r_i and r_j are identical (r = 0).
pub fn coulomb_pair_virial(
    q_i: f64,
    q_j: f64,
    r_i: [f64; 3],
    r_j: [f64; 3],
) -> Result<[[f64; 3]; 3], UcmdError> {
    check_finite(r_i)?;
    check_finite(r_j)?;
    let r_ij = vec3::sub(r_i, r_j);
    let r = vec3::norm(r_ij);
    if r == 0.0 {
        return Err(UcmdError::CoincidentAtoms);
    }
    let factor = COULOMB_KJ_MOL_NM_E2 * q_i * q_j / (r * r * r);
    let f_i = vec3::scale(r_ij, factor);
    let mut w = [[0.0_f64; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            w[a][b] = r_ij[a] * f_i[b];
        }
    }
    Ok(w)
}

/// Checks that all coordinates in `v` are finite.
/// Units: coordinates in nm (or any length unit); return is unitless.
fn check_finite(v: [f64; 3]) -> Result<(), UcmdError> {
    for &x in v.iter() {
        if !x.is_finite() {
            return Err(UcmdError::NonFiniteCoordinate);
        }
    }
    Ok(())
}
