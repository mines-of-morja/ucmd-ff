//! OPLS-Fourier dihedral term for four atoms.
//!
//! Implements the energy, Cartesian forces, and virial tensor for the
//! OPLS-Fourier dihedral potential:
//!
//! ```text
//! E = 0.5 * [V1*(1+cos(phi)) + V2*(1-cos(2*phi))
//!           + V3*(1+cos(3*phi)) + V4*(1-cos(4*phi))]
//! ```
//!
//! Units: positions in nm; V1..V4 in kJ/mol; energy in kJ/mol;
//! forces in kJ mol^-1 nm^-1; virial in kJ/mol.

use crate::dihedral_angle::{dihedral_angle, dihedral_gradient};
use crate::errors::UcmdError;
use crate::vec3::PerAtom4;
use crate::vec3::scale;

/// Computes the energy of an OPLS-Fourier dihedral term.
/// Units: r_i, r_j, r_k, r_l in nm; v1, v2, v3, v4 in kJ/mol; return in kJ/mol.
#[allow(clippy::too_many_arguments)] // the four positions and V1..V4 of the formula
pub fn opls_fourier_dihedral_energy(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
    v1: f64,
    v2: f64,
    v3: f64,
    v4: f64,
) -> Result<f64, UcmdError> {
    let phi = dihedral_angle(r_i, r_j, r_k, r_l)?;

    let phi2: f64 = 2.0 * phi;
    let phi3: f64 = 3.0 * phi;
    let phi4: f64 = 4.0 * phi;

    let e = 0.5
        * (v1 * (1.0 + phi.cos())
            + v2 * (1.0 - phi2.cos())
            + v3 * (1.0 + phi3.cos())
            + v4 * (1.0 - phi4.cos()));

    Ok(e)
}

/// Computes the forces on atoms i, j, k, l for an OPLS-Fourier dihedral term.
/// Units: r_i, r_j, r_k, r_l in nm; v1, v2, v3, v4 in kJ/mol; return forces in kJ mol^-1 nm^-1.
#[allow(clippy::too_many_arguments)] // the four positions and V1..V4 of the formula
pub fn opls_fourier_dihedral_forces(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
    v1: f64,
    v2: f64,
    v3: f64,
    v4: f64,
) -> Result<PerAtom4, UcmdError> {
    let (g_i, g_j, g_k, g_l) = dihedral_gradient(r_i, r_j, r_k, r_l)?;

    let phi = dihedral_angle(r_i, r_j, r_k, r_l)?;

    let phi2: f64 = 2.0 * phi;
    let phi3: f64 = 3.0 * phi;
    let phi4: f64 = 4.0 * phi;

    let d_e_d_phi = 0.5
        * (-v1 * phi.sin() + 2.0 * v2 * phi2.sin() - 3.0 * v3 * phi3.sin() + 4.0 * v4 * phi4.sin());

    let f_i = scale(g_i, -d_e_d_phi);
    let f_j = scale(g_j, -d_e_d_phi);
    let f_k = scale(g_k, -d_e_d_phi);
    let f_l = scale(g_l, -d_e_d_phi);

    Ok((f_i, f_j, f_k, f_l))
}

/// Computes the virial tensor for an OPLS-Fourier dihedral term.
/// Units: r_i, r_j, r_k, r_l in nm; v1, v2, v3, v4 in kJ/mol; return virial in kJ/mol.
#[allow(clippy::too_many_arguments)] // the four positions and V1..V4 of the formula
pub fn opls_fourier_dihedral_virial(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
    v1: f64,
    v2: f64,
    v3: f64,
    v4: f64,
) -> Result<[[f64; 3]; 3], UcmdError> {
    let (f_i, f_j, f_k, f_l) = opls_fourier_dihedral_forces(r_i, r_j, r_k, r_l, v1, v2, v3, v4)?;

    let mut w = [[0.0_f64; 3]; 3];

    let atoms = [r_i, r_j, r_k, r_l];
    let forces = [f_i, f_j, f_k, f_l];

    for idx in 0..4 {
        for a in 0..3 {
            for b in 0..3 {
                w[a][b] += forces[idx][a] * atoms[idx][b];
            }
        }
    }

    Ok(w)
}
