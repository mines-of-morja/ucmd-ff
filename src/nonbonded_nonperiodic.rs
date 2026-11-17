//! Evaluates all non-bonded pairs in a non-periodic system.
//! Units: positions in nm; charges in e; LJ sigma in nm; LJ epsilon in kJ/mol;
//!        energy in kJ/mol; forces in kJ mol⁻¹ nm⁻¹; virial in kJ/mol.

use crate::combination_rule::{CombinationRule, combine_eps, combine_sigma};
use crate::coulomb_pair::{coulomb_pair_energy, coulomb_pair_forces};
use crate::errors::UcmdError;
use crate::exclusions::{PairClass, build_exclusion_map};
use crate::lj_pair::{lj_pair_energy, lj_pair_forces};
use crate::vec3::{add, scale, sub};

/// Input data for a non-periodic non-bonded evaluation.
/// Units: positions in nm; charges in e; lj_sigma in nm; lj_eps in kJ/mol;
///        lj_scale_14 and coul_scale_14 are dimensionless.
#[derive(Debug, Clone, PartialEq)]
pub struct NonbondedInput {
    /// Atomic positions, unwrapped Cartesian. Units: nm.
    pub positions: Vec<[f64; 3]>,
    /// Atom-type index for each atom (index into lj_sigma, lj_eps). Units: none.
    pub types: Vec<usize>,
    /// Point charge on each atom. Units: e.
    pub charges: Vec<f64>,
    /// Per-type Lennard-Jones sigma. Units: nm.
    pub lj_sigma: Vec<f64>,
    /// Per-type Lennard-Jones epsilon. Units: kJ/mol.
    pub lj_eps: Vec<f64>,
    /// Combination rule for deriving pair sigma from per-type sigmas.
    pub rule: CombinationRule,
    /// Bonded pairs as (i, j) index tuples. Units: none.
    pub bonds: Vec<(usize, usize)>,
    /// Scale factor applied to LJ energy, force, and virial for 1-4 pairs. Dimensionless.
    pub lj_scale_14: f64,
    /// Scale factor applied to Coulomb energy, force, and virial for 1-4 pairs. Dimensionless.
    pub coul_scale_14: f64,
}

/// Result of a non-periodic non-bonded evaluation.
/// Units: energy_lj and energy_coulomb in kJ/mol; forces in kJ mol⁻¹ nm⁻¹;
///        virial in kJ/mol.
#[derive(Debug, Clone, PartialEq)]
pub struct NonbondedResult {
    /// Total Lennard-Jones energy (including 1-4 scaling). Units: kJ/mol.
    pub energy_lj: f64,
    /// Total Coulomb energy (including 1-4 scaling). Units: kJ/mol.
    pub energy_coulomb: f64,
    /// Force on each atom. Units: kJ mol⁻¹ nm⁻¹.
    pub forces: Vec<[f64; 3]>,
    /// Virial tensor W where W_ab = -∂E/∂ε_ab. Units: kJ/mol.
    pub virial: [[f64; 3]; 3],
}

/// Evaluates all non-bonded pairs in a non-periodic system.
/// Units: see `NonbondedInput` and `NonbondedResult`.
///
/// # Errors
/// * `NonFiniteCoordinate` — any position coordinate is NaN or infinite.
/// * `CoincidentAtoms` — two interacting atoms are at the same position.
/// * `IndexOutOfRange` — a bond index or type index is out of range.
pub fn evaluate_nonperiodic(input: &NonbondedInput) -> Result<NonbondedResult, UcmdError> {
    let n = input.positions.len();

    // Validate that all per-atom vectors have the same length.
    if input.types.len() != n || input.charges.len() != n {
        return Err(UcmdError::IndexOutOfRange);
    }

    // Validate that all position coordinates are finite.
    for position in &input.positions {
        for coordinate in position {
            if !coordinate.is_finite() {
                return Err(UcmdError::NonFiniteCoordinate);
            }
        }
    }

    // Validate bond indices.
    for &(bi, bj) in &input.bonds {
        if bi >= n || bj >= n {
            return Err(UcmdError::IndexOutOfRange);
        }
    }

    // Validate type indices.
    for i in 0..n {
        if input.types[i] >= input.lj_sigma.len() || input.types[i] >= input.lj_eps.len() {
            return Err(UcmdError::IndexOutOfRange);
        }
    }

    // Build the exclusion map.
    let exclusion_map = build_exclusion_map(n, &input.bonds)?;

    // Initialize accumulators.
    let mut energy_lj: f64 = 0.0;
    let mut energy_coulomb: f64 = 0.0;
    let mut forces: Vec<[f64; 3]> = vec![[0.0; 3]; n];
    let mut virial: [[f64; 3]; 3] = [[0.0; 3]; 3];

    // Pair loop.
    for i in 0..n {
        for j in (i + 1)..n {
            let cls = exclusion_map[i][j];
            if cls == PairClass::Excluded {
                continue;
            }

            let pos_i = input.positions[i];
            let pos_j = input.positions[j];
            let r_ij = sub(pos_i, pos_j);

            // Determine scale factors.
            let lj_scale = if cls == PairClass::OneFour {
                input.lj_scale_14
            } else {
                1.0
            };
            let coul_scale = if cls == PairClass::OneFour {
                input.coul_scale_14
            } else {
                1.0
            };

            // --- LJ contribution ---
            let sigma_ij = combine_sigma(
                input.rule,
                input.lj_sigma[input.types[i]],
                input.lj_sigma[input.types[j]],
            );
            let eps_ij = combine_eps(input.lj_eps[input.types[i]], input.lj_eps[input.types[j]]);

            let (e_lj, fi_lj, fj_lj, w_lj) = if sigma_ij == 0.0 || eps_ij == 0.0 {
                (0.0, [0.0; 3], [0.0; 3], [[0.0; 3]; 3])
            } else {
                let e = lj_pair_energy(sigma_ij, eps_ij, pos_i, pos_j)?;
                let (fi, fj) = lj_pair_forces(sigma_ij, eps_ij, pos_i, pos_j)?;
                let mut w = [[0.0; 3]; 3];
                for a in 0..3 {
                    for b in 0..3 {
                        w[a][b] = r_ij[a] * fi[b];
                    }
                }
                (e, fi, fj, w)
            };

            let e_lj = e_lj * lj_scale;
            let fi_lj = scale(fi_lj, lj_scale);
            let fj_lj = scale(fj_lj, lj_scale);
            let w_lj = scale_matrix(w_lj, lj_scale);

            // --- Coulomb contribution ---
            let (e_coul, fi_coul, fj_coul, w_coul) =
                if input.charges[i] == 0.0 && input.charges[j] == 0.0 {
                    (0.0, [0.0; 3], [0.0; 3], [[0.0; 3]; 3])
                } else {
                    let e = coulomb_pair_energy(input.charges[i], input.charges[j], pos_i, pos_j)?;
                    let (fi, fj) =
                        coulomb_pair_forces(input.charges[i], input.charges[j], pos_i, pos_j)?;
                    let mut w = [[0.0; 3]; 3];
                    for a in 0..3 {
                        for b in 0..3 {
                            w[a][b] = r_ij[a] * fi[b];
                        }
                    }
                    (e, fi, fj, w)
                };

            let e_coul = e_coul * coul_scale;
            let fi_coul = scale(fi_coul, coul_scale);
            let fj_coul = scale(fj_coul, coul_scale);
            let w_coul = scale_matrix(w_coul, coul_scale);

            // Accumulate.
            energy_lj += e_lj;
            energy_coulomb += e_coul;
            forces[i] = add(forces[i], add(fi_lj, fi_coul));
            forces[j] = add(forces[j], add(fj_lj, fj_coul));
            for a in 0..3 {
                for b in 0..3 {
                    virial[a][b] += w_lj[a][b] + w_coul[a][b];
                }
            }
        }
    }

    Ok(NonbondedResult {
        energy_lj,
        energy_coulomb,
        forces,
        virial,
    })
}

/// Scales a 3×3 matrix by a scalar.
/// Units: output units = input units × scalar units.
fn scale_matrix(m: [[f64; 3]; 3], s: f64) -> [[f64; 3]; 3] {
    let mut out = [[0.0; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            out[a][b] = m[a][b] * s;
        }
    }
    out
}
