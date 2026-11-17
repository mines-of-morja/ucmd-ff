//! Periodic dihedral term: energy, forces, and virial.
//!
//! Implements $E = k [1 + \cos(n\phi - \phi_0)]$ where $\phi$ is the IUPAC dihedral angle.
//!
//! Units: positions in nm; n dimensionless (integer 0..6); phi0 in rad; k in kJ/mol;
//! energy in kJ/mol; forces in kJ mol^-1 nm^-1; virial in kJ/mol.

use crate::dihedral_angle::{dihedral_angle, dihedral_gradient};
use crate::errors::UcmdError;
use crate::vec3::PerAtom4;
use crate::vec3::scale;

/// Computes the energy of a periodic dihedral term.
/// Units: r_i, r_j, r_k, r_l in nm; n dimensionless (integer 0..6); phi0 in rad; k_kjmol in kJ/mol; return in kJ/mol.
pub fn periodic_dihedral_energy(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
    n: i32,
    phi0: f64,
    k_kjmol: f64,
) -> Result<f64, UcmdError> {
    if !(0..=6).contains(&n) {
        return Err(UcmdError::IndexOutOfRange);
    }
    let phi = dihedral_angle(r_i, r_j, r_k, r_l)?;
    let n_f64 = n as f64;
    let arg = n_f64 * phi - phi0;
    let energy = k_kjmol * (1.0 + arg.cos());
    Ok(energy)
}

/// Computes the forces on atoms i, j, k, l for a periodic dihedral term.
/// Units: r_i, r_j, r_k, r_l in nm; n dimensionless (integer 0..6); phi0 in rad; k_kjmol in kJ/mol; return forces in kJ mol^-1 nm^-1.
pub fn periodic_dihedral_forces(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
    n: i32,
    phi0: f64,
    k_kjmol: f64,
) -> Result<PerAtom4, UcmdError> {
    if !(0..=6).contains(&n) {
        return Err(UcmdError::IndexOutOfRange);
    }
    let phi = dihedral_angle(r_i, r_j, r_k, r_l)?;
    let (g_i, g_j, g_k, g_l) = dihedral_gradient(r_i, r_j, r_k, r_l)?;
    let n_f64 = n as f64;
    let arg = n_f64 * phi - phi0;
    let factor = k_kjmol * n_f64 * arg.sin();
    let f_i = scale(g_i, factor);
    let f_j = scale(g_j, factor);
    let f_k = scale(g_k, factor);
    let f_l = scale(g_l, factor);
    Ok((f_i, f_j, f_k, f_l))
}

/// Computes the virial tensor for a periodic dihedral term.
/// Units: r_i, r_j, r_k, r_l in nm; n dimensionless (integer 0..6); phi0 in rad; k_kjmol in kJ/mol; return virial in kJ/mol.
pub fn periodic_dihedral_virial(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
    n: i32,
    phi0: f64,
    k_kjmol: f64,
) -> Result<[[f64; 3]; 3], UcmdError> {
    if !(0..=6).contains(&n) {
        return Err(UcmdError::IndexOutOfRange);
    }
    let (f_i, f_j, f_k, f_l) = periodic_dihedral_forces(r_i, r_j, r_k, r_l, n, phi0, k_kjmol)?;
    let atoms = [r_i, r_j, r_k, r_l];
    let forces = [f_i, f_j, f_k, f_l];
    let mut virial: [[f64; 3]; 3] = [[0.0; 3]; 3];
    for m in 0..4 {
        for a in 0..3 {
            for b in 0..3 {
                virial[a][b] += atoms[m][a] * forces[m][b];
            }
        }
    }
    Ok(virial)
}
