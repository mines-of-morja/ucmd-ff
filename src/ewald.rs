//! Classic Ewald summation of the Coulomb energy of point charges in a periodic cell (M4).
//!
//! E = E_real + E_reciprocal + E_self + E_correction, tinfoil boundary conditions, for a neutral
//! system (R4.1). With k_e the Coulomb constant, α the splitting parameter, r_ij = r_i − r_j,
//! n a lattice vector, k a reciprocal-lattice vector and V the cell volume:
//!
//! - real space (`ewald_real`): k_e Σ_{i<j} Σ_n q_i q_j erfc(α|r_ij+n|)/|r_ij+n| over the images
//!   closer than the real-space cutoff, plus the interaction of every atom with its own images,
//!   ½ k_e Σ_i q_i² Σ_{n≠0} erfc(α|n|)/|n|. A scaled pair (below) leaves out its nearest image.
//! - reciprocal space (`ewald_reciprocal`): (2π k_e/V) Σ_{0<|k|<k_c} e^(−k²/4α²)/k² |S(k)|²,
//!   S(k) = Σ_j q_j e^(i k·r_j).
//! - self term (`ewald_self`): −k_e (α/√π) Σ_i q_i².
//! - correction (`ewald_correction`): for every scaled pair (i, j, s) with the minimum-image
//!   distance r, k_e q_i q_j [s − erf(αr)]/r. The reciprocal sum contains the pair's nearest
//!   interaction k_e q_i q_j erf(αr)/r all the same; the correction removes it and adds the
//!   scaled direct interaction s k_e q_i q_j / r. An excluded pair has s = 0: it has no Coulomb
//!   interaction at all. A 1-4 pair has s = the 1-4 Coulomb scale factor.
//!
//! So the converged total changes by (s − 1) k_e q_i q_j / r when a pair gets the scale s,
//! independent of α (R4.3). Every part returns its energy, forces F = −∂E/∂r and virial
//! W_ab = −∂E/∂ε_ab (rule 2, spec/units.md §4) (R4.2).
//!
//! Units: positions and cell in nm; charges in e; α and k in nm⁻¹; energies in kJ/mol; forces in
//! kJ mol⁻¹ nm⁻¹; virial in kJ/mol.

use std::collections::HashSet;

use crate::cell::{Cell, validate_cell};
use crate::constants::COULOMB_KJ_MOL_NM_E2;
use crate::erf::erf;
use crate::errors::UcmdError;
use crate::ewald_real::ewald_real;
use crate::ewald_reciprocal::ewald_reciprocal;
use crate::exclusions::PairClass;
use crate::pair_class_map;

/// Largest total charge, in e, that counts as neutral (R4.1).
pub const NEUTRALITY_TOLERANCE_E: f64 = 1e-6;

/// Parameters of an Ewald sum.
/// Units: alpha_per_nm in nm⁻¹; real_cutoff_nm in nm; recip_cutoff_per_nm in nm⁻¹.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EwaldParameters {
    /// Splitting parameter α (the real-space interaction is damped by erfc(αr)). Units: nm⁻¹.
    pub alpha_per_nm: f64,
    /// Real-space cutoff: images at a distance r < real_cutoff_nm are summed. It may exceed half
    /// the cell height; all images within it are found. Units: nm.
    pub real_cutoff_nm: f64,
    /// Reciprocal-space cutoff: vectors with 0 < |k| < recip_cutoff_per_nm are summed. Units: nm⁻¹.
    pub recip_cutoff_per_nm: f64,
}

impl EwaldParameters {
    /// Checks that all three parameters are finite and positive.
    /// Units: none. Errors: `InvalidParameter`.
    pub fn validate(&self) -> Result<(), UcmdError> {
        for v in [
            self.alpha_per_nm,
            self.real_cutoff_nm,
            self.recip_cutoff_per_nm,
        ] {
            if !v.is_finite() || v <= 0.0 {
                return Err(UcmdError::InvalidParameter);
            }
        }
        Ok(())
    }

    /// Parameters for a given real-space cutoff and tolerance δ: α solves erfc(α r_c) = δ (by
    /// bisection), and k_c = 2α √(−ln δ), so that e^(−k_c²/4α²) = δ. The neglected terms of both
    /// sums are then of relative size δ.
    /// Units: real_cutoff_nm in nm; tolerance dimensionless, 0 < δ < 1.
    /// Errors: `InvalidParameter` for a non-positive or non-finite cutoff or δ outside (0, 1).
    pub fn from_tolerance(real_cutoff_nm: f64, tolerance: f64) -> Result<Self, UcmdError> {
        if !real_cutoff_nm.is_finite() || real_cutoff_nm <= 0.0 {
            return Err(UcmdError::InvalidParameter);
        }
        if !tolerance.is_finite() || tolerance <= 0.0 || tolerance >= 1.0 {
            return Err(UcmdError::InvalidParameter);
        }
        // erfc(x) decreases from 1 at x = 0; erfc(27) ≈ 5e-319, so x lies in (0, 27] for δ > 1e-300.
        let (mut lo, mut hi) = (0.0_f64, 27.0_f64);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if crate::erf::erfc(mid) > tolerance {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let alpha = 0.5 * (lo + hi) / real_cutoff_nm;
        let params = EwaldParameters {
            alpha_per_nm: alpha,
            real_cutoff_nm,
            recip_cutoff_per_nm: 2.0 * alpha * (-tolerance.ln()).sqrt(),
        };
        params.validate()?;
        Ok(params)
    }
}

/// A pair whose direct Coulomb interaction is scaled: 0 for an excluded pair, the 1-4 Coulomb
/// scale factor for a 1-4 pair. It applies to the nearest (minimum-image) interaction of the pair
/// only; its periodic images interact fully.
/// Units: i, j atom indices; scale dimensionless, in [0, 1].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScaledPair {
    /// First atom index. Units: none.
    pub i: usize,
    /// Second atom index (different from i). Units: none.
    pub j: usize,
    /// Factor of the direct Coulomb interaction. Units: none.
    pub scale: f64,
}

/// Energy, forces and virial of one part of the Ewald sum (or of the total).
/// Units: energy_kjmol in kJ/mol; forces in kJ mol⁻¹ nm⁻¹; virial in kJ/mol.
#[derive(Debug, Clone, PartialEq)]
pub struct EwaldPart {
    /// Energy. Units: kJ/mol.
    pub energy_kjmol: f64,
    /// Force on each atom, F = −∂E/∂r. Units: kJ mol⁻¹ nm⁻¹.
    pub forces: Vec<[f64; 3]>,
    /// Virial W_ab = −∂E/∂ε_ab. Units: kJ/mol.
    pub virial: [[f64; 3]; 3],
}

impl EwaldPart {
    /// A part with zero energy, forces and virial for n atoms. Units: as `EwaldPart`.
    pub(crate) fn zero(n: usize) -> Self {
        EwaldPart {
            energy_kjmol: 0.0,
            forces: vec![[0.0; 3]; n],
            virial: [[0.0; 3]; 3],
        }
    }
}

/// The four parts of an Ewald sum and their total.
/// Units: as `EwaldPart`.
#[derive(Debug, Clone, PartialEq)]
pub struct EwaldResult {
    /// Real-space part. Units: as `EwaldPart`.
    pub real: EwaldPart,
    /// Reciprocal-space part. Units: as `EwaldPart`.
    pub reciprocal: EwaldPart,
    /// Self term (no forces, no virial). Units: as `EwaldPart`.
    pub self_term: EwaldPart,
    /// Correction for scaled (excluded and 1-4) pairs. Units: as `EwaldPart`.
    pub correction: EwaldPart,
    /// Sum of the four parts. Units: as `EwaldPart`.
    pub total: EwaldPart,
}

/// Checks a periodic system of point charges: equal lengths, finite positions and charges, a
/// valid cell, and a total charge within 1e-6 e of zero.
/// Units: positions in nm; charges in e.
/// Errors: `IndexOutOfRange` (lengths differ), `NonFiniteCoordinate`, `InvalidParameter`
/// (non-finite charge), `InvalidCell`, `NonNeutralSystem`.
pub fn validate_system(
    positions: &[[f64; 3]],
    charges: &[f64],
    cell: &Cell,
) -> Result<(), UcmdError> {
    if positions.len() != charges.len() {
        return Err(UcmdError::IndexOutOfRange);
    }
    if positions.iter().flatten().any(|c| !c.is_finite()) {
        return Err(UcmdError::NonFiniteCoordinate);
    }
    if charges.iter().any(|q| !q.is_finite()) {
        return Err(UcmdError::InvalidParameter);
    }
    validate_cell(cell)?;
    let total: f64 = charges.iter().sum();
    if total.abs() > NEUTRALITY_TOLERANCE_E {
        return Err(UcmdError::NonNeutralSystem);
    }
    Ok(())
}

/// Checks the scaled pairs of an n-atom system and returns them as a set of (min, max) index
/// pairs. Units: none.
/// Errors: `IndexOutOfRange` (index ≥ n), `InvalidParameter` (i = j, a scale outside [0, 1] or
/// not finite, or a pair given twice).
pub(crate) fn scaled_pair_set(
    n: usize,
    scaled: &[ScaledPair],
) -> Result<HashSet<(usize, usize)>, UcmdError> {
    let mut set = HashSet::with_capacity(scaled.len());
    for p in scaled {
        if p.i >= n || p.j >= n {
            return Err(UcmdError::IndexOutOfRange);
        }
        if p.i == p.j || !p.scale.is_finite() || !(0.0..=1.0).contains(&p.scale) {
            return Err(UcmdError::InvalidParameter);
        }
        if !set.insert((p.i.min(p.j), p.i.max(p.j))) {
            return Err(UcmdError::InvalidParameter);
        }
    }
    Ok(set)
}

/// The scaled pairs that the bond graph defines (spec/system-package.md §5.2): pairs at shortest
/// path length 1 or 2 get scale 0, 1-4 pairs (path length 3) get `coulomb_scale_14`. A 1-4 scale
/// of exactly 1 adds no pair. Pairs are returned with i < j, in increasing order.
/// Units: n_atoms and bond indices none; coulomb_scale_14 dimensionless, in [0, 1].
/// Errors: `IndexOutOfRange` (bond index ≥ n_atoms), `InvalidParameter` (scale outside [0, 1]).
pub fn scaled_pairs_from_bonds(
    n_atoms: usize,
    bonds: &[(usize, usize)],
    coulomb_scale_14: f64,
) -> Result<Vec<ScaledPair>, UcmdError> {
    if !coulomb_scale_14.is_finite() || !(0.0..=1.0).contains(&coulomb_scale_14) {
        return Err(UcmdError::InvalidParameter);
    }
    let map = pair_class_map::build(n_atoms, bonds)?;
    let mut pairs = Vec::new();
    for i in 0..n_atoms {
        for j in (i + 1)..n_atoms {
            match map.class_at(i, j)? {
                PairClass::Excluded => pairs.push(ScaledPair { i, j, scale: 0.0 }),
                PairClass::OneFour if coulomb_scale_14 != 1.0 => pairs.push(ScaledPair {
                    i,
                    j,
                    scale: coulomb_scale_14,
                }),
                _ => {}
            }
        }
    }
    Ok(pairs)
}

/// Self term −k_e (α/√π) Σ_i q_i². It depends on neither positions nor cell: its forces and
/// virial are zero.
/// Units: charges in e; returned part in kJ/mol (forces kJ mol⁻¹ nm⁻¹, virial kJ/mol).
/// Errors: as `validate_system` and `EwaldParameters::validate`.
pub fn ewald_self(
    positions: &[[f64; 3]],
    charges: &[f64],
    cell: &Cell,
    params: &EwaldParameters,
) -> Result<EwaldPart, UcmdError> {
    validate_system(positions, charges, cell)?;
    params.validate()?;
    let sum_q2: f64 = charges.iter().map(|q| q * q).sum();
    let mut part = EwaldPart::zero(charges.len());
    part.energy_kjmol =
        -COULOMB_KJ_MOL_NM_E2 * params.alpha_per_nm / std::f64::consts::PI.sqrt() * sum_q2;
    Ok(part)
}

/// Correction for scaled pairs: Σ k_e q_i q_j [s − erf(αr)]/r over the scaled pairs, with r the
/// minimum-image distance (see the module documentation), with its forces and virial
/// (pair form W = Σ r_ij ⊗ F_ij).
/// Units: positions and cell in nm; charges in e; returned part in kJ/mol, kJ mol⁻¹ nm⁻¹, kJ/mol.
/// Errors: as `validate_system`, `EwaldParameters::validate` and the scaled-pair checks;
/// `CoincidentAtoms` if a scaled pair is at distance zero.
pub fn ewald_correction(
    positions: &[[f64; 3]],
    charges: &[f64],
    cell: &Cell,
    params: &EwaldParameters,
    scaled: &[ScaledPair],
) -> Result<EwaldPart, UcmdError> {
    validate_system(positions, charges, cell)?;
    params.validate()?;
    scaled_pair_set(positions.len(), scaled)?;
    let alpha = params.alpha_per_nm;
    let two_alpha_over_sqrt_pi = 2.0 * alpha / std::f64::consts::PI.sqrt();
    let mut part = EwaldPart::zero(positions.len());
    for p in scaled {
        let qq = COULOMB_KJ_MOL_NM_E2 * charges[p.i] * charges[p.j];
        if qq == 0.0 {
            continue;
        }
        let d = crate::cell::minimum_image(positions[p.i], positions[p.j], cell)?;
        let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        if r2 == 0.0 {
            return Err(UcmdError::CoincidentAtoms);
        }
        let r = r2.sqrt();
        let ar = alpha * r;
        let bracket = p.scale - erf(ar);
        part.energy_kjmol += qq * bracket / r;
        // dE/dr = qq [ -(s - erf(ar))/r^2 - (2 alpha/sqrt(pi)) exp(-(ar)^2)/r ]
        let de_dr = qq * (-bracket / r2 - two_alpha_over_sqrt_pi * (-ar * ar).exp() / r);
        let g = -de_dr / r; // F_i = g * d
        add_pair(&mut part, p.i, p.j, d, g);
    }
    Ok(part)
}

/// Adds the pair force F_i = g d (and F_j = −F_i) and its virial d ⊗ F_i to a part.
/// Units: d in nm; g in kJ mol⁻¹ nm⁻².
pub(crate) fn add_pair(part: &mut EwaldPart, i: usize, j: usize, d: [f64; 3], g: f64) {
    let f = [g * d[0], g * d[1], g * d[2]];
    for c in 0..3 {
        part.forces[i][c] += f[c];
        part.forces[j][c] -= f[c];
    }
    for a in 0..3 {
        for b in 0..3 {
            part.virial[a][b] += d[a] * f[b];
        }
    }
}

/// The complete Ewald sum: real space, reciprocal space, self term and the correction for scaled
/// pairs, each with energy, forces and virial, and their total.
/// Units: positions and cell in nm; charges in e; parameters as `EwaldParameters`; results as
/// `EwaldPart`.
/// Errors: `NonNeutralSystem` when |Σq| > 1e-6 e, and the errors of the parts.
pub fn ewald(
    positions: &[[f64; 3]],
    charges: &[f64],
    cell: &Cell,
    params: &EwaldParameters,
    scaled: &[ScaledPair],
) -> Result<EwaldResult, UcmdError> {
    let real = ewald_real(positions, charges, cell, params, scaled)?;
    let reciprocal = ewald_reciprocal(positions, charges, cell, params)?;
    let self_term = ewald_self(positions, charges, cell, params)?;
    let correction = ewald_correction(positions, charges, cell, params, scaled)?;
    let mut total = EwaldPart::zero(positions.len());
    for part in [&real, &reciprocal, &self_term, &correction] {
        total.energy_kjmol += part.energy_kjmol;
        for (t, f) in total.forces.iter_mut().zip(&part.forces) {
            for c in 0..3 {
                t[c] += f[c];
            }
        }
        for a in 0..3 {
            for b in 0..3 {
                total.virial[a][b] += part.virial[a][b];
            }
        }
    }
    Ok(EwaldResult {
        real,
        reciprocal,
        self_term,
        correction,
        total,
    })
}
