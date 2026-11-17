//! Lennard-Jones 12-6 pair potential: energy, forces, and virial for a single pair.
//!
//! Units: positions in nm, energy in kJ/mol, forces in kJ mol^-1 nm^-1,
//! virial in kJ/mol, sigma in nm, epsilon in kJ/mol.

use crate::errors::UcmdError;
use crate::vec3;

/// Computes the Lennard-Jones 12-6 pair energy.
/// Units: sigma_nm in nm; eps_kjmol in kJ/mol; r_i, r_j in nm; return in kJ/mol.
///
/// E = 4 * eps * [(sigma/r)^12 - (sigma/r)^6], where r = |r_i - r_j|.
///
/// # Errors
/// * `NonFiniteCoordinate` — any coordinate of r_i or r_j is NaN or infinite.
/// * `CoincidentAtoms` — r_i == r_j (distance is zero).
pub fn lj_pair_energy(
    sigma_nm: f64,
    eps_kjmol: f64,
    r_i: [f64; 3],
    r_j: [f64; 3],
) -> Result<f64, UcmdError> {
    check_inputs(r_i, r_j)?;
    if sigma_nm == 0.0 || eps_kjmol == 0.0 {
        return Ok(0.0);
    }
    let r_ij = vec3::sub(r_i, r_j);
    let r = vec3::norm(r_ij);
    let sr = sigma_nm / r;
    let sr6 = sr.powi(6);
    let sr12 = sr6 * sr6;
    let e = 4.0 * eps_kjmol * (sr12 - sr6);
    Ok(e)
}

/// Computes the Lennard-Jones 12-6 pair forces on atoms i and j.
/// Units: sigma_nm in nm; eps_kjmol in kJ/mol; r_i, r_j in nm;
///        returned forces in kJ mol^-1 nm^-1.
///
/// Returns `(F_i, F_j)` where `F_i` is the force on atom i due to atom j
/// and `F_j = -F_i` (Newton's third law).
///
/// # Errors
/// * `NonFiniteCoordinate` — any coordinate of r_i or r_j is NaN or infinite.
/// * `CoincidentAtoms` — r_i == r_j (distance is zero).
pub fn lj_pair_forces(
    sigma_nm: f64,
    eps_kjmol: f64,
    r_i: [f64; 3],
    r_j: [f64; 3],
) -> Result<([f64; 3], [f64; 3]), UcmdError> {
    check_inputs(r_i, r_j)?;
    if sigma_nm == 0.0 || eps_kjmol == 0.0 {
        return Ok(([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]));
    }
    let r_ij = vec3::sub(r_i, r_j);
    let r = vec3::norm(r_ij);
    let sr = sigma_nm / r;
    let sr6 = sr.powi(6);
    let sr12 = sr6 * sr6;
    let scalar = (24.0 * eps_kjmol / (r * r)) * (2.0 * sr12 - sr6);
    let f_i = vec3::scale(r_ij, scalar);
    let f_j = vec3::scale(f_i, -1.0);
    Ok((f_i, f_j))
}

/// Computes the Lennard-Jones 12-6 pair virial (3×3 matrix).
/// Units: sigma_nm in nm; eps_kjmol in kJ/mol; r_i, r_j in nm;
///        returned virial in kJ/mol.
///
/// W = r_ij ⊗ F_i, where r_ij = r_i - r_j and F_i is the force on i due to j.
/// W[a][b] = r_ij[a] * F_i[b].
///
/// # Errors
/// * `NonFiniteCoordinate` — any coordinate of r_i or r_j is NaN or infinite.
/// * `CoincidentAtoms` — r_i == r_j (distance is zero).
pub fn lj_pair_virial(
    sigma_nm: f64,
    eps_kjmol: f64,
    r_i: [f64; 3],
    r_j: [f64; 3],
) -> Result<[[f64; 3]; 3], UcmdError> {
    check_inputs(r_i, r_j)?;
    if sigma_nm == 0.0 || eps_kjmol == 0.0 {
        return Ok([[0.0; 3]; 3]);
    }
    let r_ij = vec3::sub(r_i, r_j);
    let r = vec3::norm(r_ij);
    let sr = sigma_nm / r;
    let sr6 = sr.powi(6);
    let sr12 = sr6 * sr6;
    let scalar = (24.0 * eps_kjmol / (r * r)) * (2.0 * sr12 - sr6);
    let f_i = vec3::scale(r_ij, scalar);
    let mut w: [[f64; 3]; 3] = [[0.0; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            w[a][b] = r_ij[a] * f_i[b];
        }
    }
    Ok(w)
}

/// Validates that all coordinates are finite and that the two positions are not coincident.
///
/// # Errors
/// * `NonFiniteCoordinate` — any coordinate is NaN or infinite.
/// * `CoincidentAtoms` — r_i and r_j are identical.
fn check_inputs(r_i: [f64; 3], r_j: [f64; 3]) -> Result<(), UcmdError> {
    for a in 0..3 {
        if !r_i[a].is_finite() || !r_j[a].is_finite() {
            return Err(UcmdError::NonFiniteCoordinate);
        }
    }
    if r_i[0] == r_j[0] && r_i[1] == r_j[1] && r_i[2] == r_j[2] {
        return Err(UcmdError::CoincidentAtoms);
    }
    Ok(())
}
