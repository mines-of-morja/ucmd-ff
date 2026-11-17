//! Real-space part of the Ewald sum (M4, R4.1, R4.2).
//!
//! E_real = k_e Σ_{i<j} Σ_n q_i q_j erfc(α|r_ij+n|)/|r_ij+n|  +  ½ k_e Σ_i q_i² Σ_{n≠0} erfc(α|n|)/|n|,
//! over the lattice vectors n with |r_ij+n| (or |n|) below the real-space cutoff. A scaled pair
//! (see `ewald::ScaledPair`) leaves out its nearest image, which `ewald::ewald_correction` treats.
//! The image search is complete for any cutoff, also one above half the cell height.
//!
//! Pair force on i: F_i = k_e q_i q_j [erfc(αr)/r + (2α/√π) e^(−α²r²)] d/r², d the image
//! displacement; virial W = Σ d ⊗ F_i, and for the self images ½ Σ_n n ⊗ F(n).
//!
//! Units: positions and cell in nm; charges in e; α in nm⁻¹; energy in kJ/mol; forces in
//! kJ mol⁻¹ nm⁻¹; virial in kJ/mol.

use crate::cell::{Cell, minimum_image};
use crate::constants::COULOMB_KJ_MOL_NM_E2;
use crate::erf::erfc;
use crate::errors::UcmdError;
use crate::ewald::{
    EwaldParameters, EwaldPart, ScaledPair, add_pair, scaled_pair_set, validate_system,
};
use crate::minimum_image::perpendicular_heights;

/// Lattice vectors of the cell as rows. Units: nm.
pub(crate) fn lattice(cell: &Cell) -> [[f64; 3]; 3] {
    [
        [cell.ax, 0.0, 0.0],
        [cell.bx, cell.by, 0.0],
        [cell.cx, cell.cy, cell.cz],
    ]
}

/// Fractional coordinates s of a vector d = s0 a + s1 b + s2 c (lower-triangular cell).
/// Units: d in nm; result dimensionless.
fn fractional(d: [f64; 3], cell: &Cell) -> [f64; 3] {
    let s2 = d[2] / cell.cz;
    let s1 = (d[1] - s2 * cell.cy) / cell.by;
    let s0 = (d[0] - s1 * cell.bx - s2 * cell.cx) / cell.ax;
    [s0, s1, s2]
}

/// Integer range of the lattice index along one direction for images closer than r_c: the
/// fractional coordinate s + n of an image at distance below r_c satisfies |s + n| < r_c/h.
/// Units: s dimensionless; reach = r_c/h dimensionless.
fn index_range(s: f64, reach: f64) -> std::ops::RangeInclusive<i64> {
    let lo = (-reach - s).ceil() as i64;
    let hi = (reach - s).floor() as i64;
    lo..=hi
}

/// Real-space part of the Ewald sum with forces and virial (see the module documentation).
/// Units: positions and cell in nm; charges in e; results in kJ/mol, kJ mol⁻¹ nm⁻¹, kJ/mol.
/// Errors: as `ewald::validate_system`, `EwaldParameters::validate` and the scaled-pair checks
/// (`IndexOutOfRange`, `InvalidParameter`); `CoincidentAtoms` if two charged atoms coincide
/// (modulo a lattice vector) in an interaction that is summed.
pub fn ewald_real(
    positions: &[[f64; 3]],
    charges: &[f64],
    cell: &Cell,
    params: &EwaldParameters,
    scaled: &[ScaledPair],
) -> Result<EwaldPart, UcmdError> {
    validate_system(positions, charges, cell)?;
    params.validate()?;
    let special = scaled_pair_set(positions.len(), scaled)?;
    let n = positions.len();
    let lat = lattice(cell);
    let heights = perpendicular_heights(cell)?;
    let rc = params.real_cutoff_nm;
    let rc2 = rc * rc;
    let alpha = params.alpha_per_nm;
    let two_alpha_over_sqrt_pi = 2.0 * alpha / std::f64::consts::PI.sqrt();
    let reach = [rc / heights[0], rc / heights[1], rc / heights[2]];

    // Pair factor g(r) with F = g d: k_e q q [erfc(ar)/r + (2 alpha/sqrt(pi)) exp(-(ar)^2)] / r^2.
    let pair = |qq: f64, r2: f64| -> (f64, f64) {
        let r = r2.sqrt();
        let ar = alpha * r;
        let ec = erfc(ar);
        let energy = qq * ec / r;
        let g = qq * (ec / r + two_alpha_over_sqrt_pi * (-ar * ar).exp()) / r2;
        (energy, g)
    };

    let mut part = EwaldPart::zero(n);
    for i in 0..n {
        for j in (i + 1)..n {
            let qq = COULOMB_KJ_MOL_NM_E2 * charges[i] * charges[j];
            if qq == 0.0 {
                continue;
            }
            let d0 = minimum_image(positions[i], positions[j], cell)?;
            let skip_nearest = special.contains(&(i, j));
            let s = fractional(d0, cell);
            for na in index_range(s[0], reach[0]) {
                for nb in index_range(s[1], reach[1]) {
                    for nc in index_range(s[2], reach[2]) {
                        if skip_nearest && na == 0 && nb == 0 && nc == 0 {
                            continue;
                        }
                        let (fa, fb, fc) = (na as f64, nb as f64, nc as f64);
                        let d = [
                            d0[0] + fa * lat[0][0] + fb * lat[1][0] + fc * lat[2][0],
                            d0[1] + fa * lat[0][1] + fb * lat[1][1] + fc * lat[2][1],
                            d0[2] + fa * lat[0][2] + fb * lat[1][2] + fc * lat[2][2],
                        ];
                        let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                        if r2 >= rc2 {
                            continue;
                        }
                        if r2 == 0.0 {
                            return Err(UcmdError::CoincidentAtoms);
                        }
                        let (e, g) = pair(qq, r2);
                        part.energy_kjmol += e;
                        add_pair(&mut part, i, j, d, g);
                    }
                }
            }
        }
    }

    // Every atom with its own images: 1/2 k_e sum_i q_i^2 sum_{n != 0} erfc(a|n|)/|n|.
    let sum_q2: f64 = charges.iter().map(|q| q * q).sum();
    let qq_self = COULOMB_KJ_MOL_NM_E2 * sum_q2;
    if qq_self != 0.0 {
        for na in index_range(0.0, reach[0]) {
            for nb in index_range(0.0, reach[1]) {
                for nc in index_range(0.0, reach[2]) {
                    if na == 0 && nb == 0 && nc == 0 {
                        continue;
                    }
                    let (fa, fb, fc) = (na as f64, nb as f64, nc as f64);
                    let v = [
                        fa * lat[0][0] + fb * lat[1][0] + fc * lat[2][0],
                        fa * lat[0][1] + fb * lat[1][1] + fc * lat[2][1],
                        fa * lat[0][2] + fb * lat[1][2] + fc * lat[2][2],
                    ];
                    let r2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
                    if r2 >= rc2 {
                        continue;
                    }
                    let (e, g) = pair(qq_self, r2);
                    part.energy_kjmol += 0.5 * e;
                    for a in 0..3 {
                        for b in 0..3 {
                            part.virial[a][b] += 0.5 * g * v[a] * v[b];
                        }
                    }
                }
            }
        }
    }
    Ok(part)
}
