//! Reciprocal-space part of the Ewald sum (M4, R4.1, R4.2).
//!
//! E_rec = (2π k_e/V) Σ_{0<|k|<k_c} f(k) |S(k)|²,  f(k) = e^(−k²/4α²)/k²,  S(k) = Σ_j q_j e^(i k·r_j),
//! with k = 2π (m0 a* + m1 b* + m2 c*) (a*·a = 1). The sum runs over half of k-space, k and −k
//! giving equal terms.
//!
//! - Force: F_j = (4π k_e q_j/V) Σ_k f(k) k [C(k) sin(k·r_j) − S_s(k) cos(k·r_j)], where
//!   C = Σ q cos(k·r) and S_s = Σ q sin(k·r).
//! - Virial: W_ab = (2π k_e/V) Σ_k f(k) |S(k)|² [δ_ab − 2 k_a k_b (1/k² + 1/(4α²))]. Derivation:
//!   under r_b → r_b + ε_ab r_a (positions and lattice vectors) every k·r_j is unchanged,
//!   ∂k_c/∂ε_ab = −δ_ca k_b, so ∂k²/∂ε_ab = −2 k_a k_b, and ∂V/∂ε_ab = V δ_ab.
//!
//! The set of k vectors is the one of the current cell; the virial treats it as fixed, which is
//! exact unless a |k| lies exactly at the cutoff.
//!
//! Units: positions and cell in nm; charges in e; α and k in nm⁻¹; energy in kJ/mol; forces in
//! kJ mol⁻¹ nm⁻¹; virial in kJ/mol.

use crate::cell::{Cell, cell_volume};
use crate::constants::COULOMB_KJ_MOL_NM_E2;
use crate::errors::UcmdError;
use crate::ewald::{EwaldParameters, EwaldPart, validate_system};

/// Reciprocal vectors with the factor 2π: g_i · a_j = 2π δ_ij for the lattice rows a_j of the
/// lower-triangular cell. Units: nm⁻¹.
pub(crate) fn reciprocal_vectors(cell: &Cell) -> [[f64; 3]; 3] {
    let two_pi = 2.0 * std::f64::consts::PI;
    let (ax, bx, by, cx, cy, cz) = (cell.ax, cell.bx, cell.by, cell.cx, cell.cy, cell.cz);
    // Columns of the inverse of the lattice matrix (rows a, b, c).
    let g0 = [
        1.0 / ax,
        -bx / (ax * by),
        (bx * cy - by * cx) / (ax * by * cz),
    ];
    let g1 = [0.0, 1.0 / by, -cy / (by * cz)];
    let g2 = [0.0, 0.0, 1.0 / cz];
    let scale = |g: [f64; 3]| [two_pi * g[0], two_pi * g[1], two_pi * g[2]];
    [scale(g0), scale(g1), scale(g2)]
}

/// Reciprocal-space part of the Ewald sum with forces and virial (see the module documentation).
/// Units: positions and cell in nm; charges in e; results in kJ/mol, kJ mol⁻¹ nm⁻¹, kJ/mol.
/// Errors: as `ewald::validate_system` and `EwaldParameters::validate`.
pub fn ewald_reciprocal(
    positions: &[[f64; 3]],
    charges: &[f64],
    cell: &Cell,
    params: &EwaldParameters,
) -> Result<EwaldPart, UcmdError> {
    validate_system(positions, charges, cell)?;
    params.validate()?;
    let n = positions.len();
    let volume = cell_volume(cell)?;
    let g = reciprocal_vectors(cell);
    let alpha = params.alpha_per_nm;
    let kc = params.recip_cutoff_per_nm;
    let kc2 = kc * kc;
    let inv_4a2 = 1.0 / (4.0 * alpha * alpha);
    let two_pi = 2.0 * std::f64::consts::PI;
    // |m_i| = |k . a_i| / (2 pi) <= k_c |a_i| / (2 pi).
    let lengths = [
        cell.ax,
        (cell.bx * cell.bx + cell.by * cell.by).sqrt(),
        (cell.cx * cell.cx + cell.cy * cell.cy + cell.cz * cell.cz).sqrt(),
    ];
    let m_max = lengths.map(|len| (kc * len / two_pi).floor() as i64);
    // Half space: factor 2 on every term.
    let pref = 2.0 * two_pi * COULOMB_KJ_MOL_NM_E2 / volume;

    let mut part = EwaldPart::zero(n);
    let mut cos_t = vec![0.0_f64; n];
    let mut sin_t = vec![0.0_f64; n];
    for m0 in 0..=m_max[0] {
        for m1 in -m_max[1]..=m_max[1] {
            for m2 in -m_max[2]..=m_max[2] {
                let in_half = m0 > 0 || (m0 == 0 && (m1 > 0 || (m1 == 0 && m2 > 0)));
                if !in_half {
                    continue;
                }
                let (f0, f1, f2) = (m0 as f64, m1 as f64, m2 as f64);
                let k = [
                    f0 * g[0][0] + f1 * g[1][0] + f2 * g[2][0],
                    f0 * g[0][1] + f1 * g[1][1] + f2 * g[2][1],
                    f0 * g[0][2] + f1 * g[1][2] + f2 * g[2][2],
                ];
                let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                if k2 >= kc2 {
                    continue;
                }
                let f = (-k2 * inv_4a2).exp() / k2;
                let (mut c_sum, mut s_sum) = (0.0_f64, 0.0_f64);
                for j in 0..n {
                    let r = positions[j];
                    let theta = k[0] * r[0] + k[1] * r[1] + k[2] * r[2];
                    let (s, c) = theta.sin_cos();
                    cos_t[j] = c;
                    sin_t[j] = s;
                    c_sum += charges[j] * c;
                    s_sum += charges[j] * s;
                }
                let s2 = c_sum * c_sum + s_sum * s_sum;
                let e = pref * f * s2;
                part.energy_kjmol += e;
                for j in 0..n {
                    let w = 2.0 * pref * f * charges[j] * (c_sum * sin_t[j] - s_sum * cos_t[j]);
                    for c in 0..3 {
                        part.forces[j][c] += w * k[c];
                    }
                }
                let b_fac = 2.0 * (1.0 / k2 + inv_4a2);
                for a in 0..3 {
                    for b in 0..3 {
                        let delta = if a == b { 1.0 } else { 0.0 };
                        part.virial[a][b] += e * (delta - b_fac * k[a] * k[b]);
                    }
                }
            }
        }
    }
    Ok(part)
}
