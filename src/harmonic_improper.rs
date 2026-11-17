//! Harmonic improper torsion term.
//!
//! The improper angle ξ is the IUPAC dihedral angle of atoms (i, j, k, l).
//! The energy is E = ½ k (ξ − ξ₀)² where the deviation (ξ − ξ₀) is wrapped
//! to (−π, π] before squaring, making the potential periodic with period 2π.
//!
//! Units: positions in nm; angles in rad; force constant in kJ mol⁻¹ rad⁻²;
//! energy in kJ/mol; forces in kJ mol⁻¹ nm⁻¹; virial in kJ/mol.

use crate::dihedral_angle::{dihedral_angle, dihedral_gradient};
use crate::errors::UcmdError;
use crate::vec3::PerAtom4;
use crate::vec3::scale;

/// Computes the energy of a harmonic improper term.
/// Units: r_i, r_j, r_k, r_l in nm; xi0 in rad; k_kjmol_rad2 in kJ mol⁻¹ rad⁻²; return in kJ/mol.
pub fn harmonic_improper_energy(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
    xi0: f64,
    k_kjmol_rad2: f64,
) -> Result<f64, UcmdError> {
    let xi = dihedral_angle(r_i, r_j, r_k, r_l)?;
    let delta = wrap_delta(xi - xi0);
    Ok(0.5 * k_kjmol_rad2 * delta * delta)
}

/// Computes the forces on atoms i, j, k, l for a harmonic improper term.
/// Units: r_i, r_j, r_k, r_l in nm; xi0 in rad; k_kjmol_rad2 in kJ mol⁻¹ rad⁻²;
/// return forces in kJ mol⁻¹ nm⁻¹.
pub fn harmonic_improper_forces(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
    xi0: f64,
    k_kjmol_rad2: f64,
) -> Result<PerAtom4, UcmdError> {
    let xi = dihedral_angle(r_i, r_j, r_k, r_l)?;
    let delta = wrap_delta(xi - xi0);

    if k_kjmol_rad2 == 0.0 {
        return Ok(([0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3]));
    }

    let (g_i, g_j, g_k, g_l) = dihedral_gradient(r_i, r_j, r_k, r_l)?;
    let factor = -k_kjmol_rad2 * delta;

    let f_i = scale(g_i, factor);
    let f_j = scale(g_j, factor);
    let f_k = scale(g_k, factor);
    let f_l = scale(g_l, factor);

    Ok((f_i, f_j, f_k, f_l))
}

/// Computes the virial tensor for a harmonic improper term.
/// Units: r_i, r_j, r_k, r_l in nm; xi0 in rad; k_kjmol_rad2 in kJ mol⁻¹ rad⁻²;
/// return W as a 3×3 matrix in kJ/mol.
pub fn harmonic_improper_virial(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
    xi0: f64,
    k_kjmol_rad2: f64,
) -> Result<[[f64; 3]; 3], UcmdError> {
    let xi = dihedral_angle(r_i, r_j, r_k, r_l)?;
    let delta = wrap_delta(xi - xi0);

    if k_kjmol_rad2 == 0.0 {
        return Ok([[0.0; 3]; 3]);
    }

    let (g_i, g_j, g_k, g_l) = dihedral_gradient(r_i, r_j, r_k, r_l)?;
    let factor = -k_kjmol_rad2 * delta;

    let f_i = scale(g_i, factor);
    let f_j = scale(g_j, factor);
    let f_k = scale(g_k, factor);
    let f_l = scale(g_l, factor);

    let mut w = [[0.0; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            w[a][b] = r_i[a] * f_i[b] + r_j[a] * f_j[b] + r_k[a] * f_k[b] + r_l[a] * f_l[b];
        }
    }

    Ok(w)
}

/// Wraps a deviation to (−π, π].
/// Units: input and output in rad.
fn wrap_delta(delta: f64) -> f64 {
    let pi: f64 = std::f64::consts::PI;
    let two_pi: f64 = 2.0 * pi;
    if delta > pi {
        delta - two_pi
    } else if delta <= -pi {
        delta + two_pi
    } else {
        delta
    }
}
