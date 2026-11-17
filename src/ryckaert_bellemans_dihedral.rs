//! Ryckaert-Bellemans dihedral term: energy, forces, and virial.
//! Units: positions in nm; coefficients in kJ/mol; energy in kJ/mol;
//! forces in kJ mol^-1 nm^-1; virial in kJ/mol.

use crate::dihedral_angle::{dihedral_angle, dihedral_gradient};
use crate::errors::UcmdError;
use crate::opls_to_rb::opls_to_rb as opls_to_rb_impl;
use crate::rb_dihedral_math::{rb_denergy_dphi, rb_energy};
use crate::vec3::PerAtom4;
use crate::vec3::scale;

/// Computes the energy of a Ryckaert-Bellemans dihedral term.
/// Units: r_i, r_j, r_k, r_l in nm; c in kJ/mol; return in kJ/mol.
/// E = sum_{n=0..5} c[n] * cos^n(psi), where psi = phi - pi and phi is the IUPAC dihedral.
pub fn ryckaert_bellemans_dihedral_energy(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
    c: [f64; 6],
) -> Result<f64, UcmdError> {
    let phi = dihedral_angle(r_i, r_j, r_k, r_l)?;
    Ok(rb_energy(phi, c))
}

/// Computes the forces on atoms i, j, k, l for a Ryckaert-Bellemans dihedral term.
/// Units: r_i, r_j, r_k, r_l in nm; c in kJ/mol; return forces in kJ mol^-1 nm^-1.
/// F_a = -(dE/dphi) * grad_a(phi), where grad_a is the position gradient of phi.
pub fn ryckaert_bellemans_dihedral_forces(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
    c: [f64; 6],
) -> Result<PerAtom4, UcmdError> {
    let phi = dihedral_angle(r_i, r_j, r_k, r_l)?;
    let (g_i, g_j, g_k, g_l) = dihedral_gradient(r_i, r_j, r_k, r_l)?;
    let factor = -rb_denergy_dphi(phi, c);
    Ok((
        scale(g_i, factor),
        scale(g_j, factor),
        scale(g_k, factor),
        scale(g_l, factor),
    ))
}

/// Computes the virial tensor for a Ryckaert-Bellemans dihedral term.
/// Units: r_i, r_j, r_k, r_l in nm; c in kJ/mol; return W in kJ/mol.
/// W[a][b] = r_i[a]*F_i[b] + r_j[a]*F_j[b] + r_k[a]*F_k[b] + r_l[a]*F_l[b].
pub fn ryckaert_bellemans_dihedral_virial(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
    c: [f64; 6],
) -> Result<[[f64; 3]; 3], UcmdError> {
    let (f_i, f_j, f_k, f_l) = ryckaert_bellemans_dihedral_forces(r_i, r_j, r_k, r_l, c)?;

    let mut w: [[f64; 3]; 3] = [[0.0; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            w[a][b] = r_i[a] * f_i[b] + r_j[a] * f_j[b] + r_k[a] * f_k[b] + r_l[a] * f_l[b];
        }
    }
    Ok(w)
}

/// Converts OPLS-Fourier torsion coefficients to Ryckaert-Bellemans coefficients.
/// Units: input v and output c are both in kJ/mol.
/// Re-export of `crate::opls_to_rb::opls_to_rb`.
pub fn opls_to_rb(v: [f64; 4]) -> [f64; 6] {
    opls_to_rb_impl(v)
}
