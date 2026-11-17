//! Evaluates non-bonded interactions (Lennard-Jones and Coulomb) in a periodic system.
//!
//! Units: positions in nm; cell in nm; lj_cutoff_nm in nm;
//!        energies in kJ/mol; forces in kJ mol⁻¹ nm⁻¹; virial in kJ/mol.
//!
//! For each pair (i, j) with i < j, the minimum-image displacement is computed.
//! LJ pairs contribute only if the minimum-image distance is strictly less than
//! `lj_cutoff_nm`; Coulomb pairs always contribute. Pair parameters are derived
//! via the combination rule. Exclusions and 1-4 scaling are handled via the
//! pair class map.

use crate::cell::{Cell, validate_cell};
use crate::combination_rule::{combine_eps, combine_sigma};
use crate::coulomb_pair::{coulomb_pair_energy, coulomb_pair_forces};
use crate::errors::UcmdError;
use crate::exclusions::PairClass;
use crate::lj_pair::{lj_pair_energy, lj_pair_forces};
use crate::min_image_disp::minimum_image_displacement;
use crate::minimum_image::{perpendicular_heights, validate_cutoff};
use crate::nonbonded_nonperiodic::NonbondedInput;
use crate::pair_class_map::build as build_pair_class_map;
use crate::vec3::add;

/// Result of a periodic non-bonded evaluation.
/// Units: energies in kJ/mol; forces in kJ mol⁻¹ nm⁻¹; virial in kJ/mol.
#[derive(Debug, Clone, PartialEq)]
pub struct NonbondedPeriodicResult {
    /// Total Lennard-Jones energy. Units: kJ/mol.
    pub energy_lj_kjmol: f64,
    /// Total Coulomb energy. Units: kJ/mol.
    pub energy_coulomb_kjmol: f64,
    /// Forces on each atom. Units: kJ mol⁻¹ nm⁻¹.
    pub forces: Vec<[f64; 3]>,
    /// Virial tensor W[a][b] = Σ_{i<j} d_ij[a] * F_ij[b]. Units: kJ/mol.
    pub virial: [[f64; 3]; 3],
}

/// Evaluates all non-bonded pairs in a periodic system.
/// Units: positions in nm; cell in nm; lj_cutoff_nm in nm;
///        energies in kJ/mol; forces in kJ mol⁻¹ nm⁻¹; virial in kJ/mol.
///
/// # Errors
/// * `InvalidCell` — cell invalid or cutoff ≥ min(perpendicular_heights)/2.
/// * `NonFiniteCoordinate` — any position is NaN or infinite.
/// * `CoincidentAtoms` — atoms coincide after minimum-image correction.
pub fn evaluate_periodic(
    input: &NonbondedInput,
    cell: &Cell,
    lj_cutoff_nm: f64,
) -> Result<NonbondedPeriodicResult, UcmdError> {
    // 1. Validation
    validate_cell(cell)?;
    perpendicular_heights(cell)?;
    validate_cutoff(cell, lj_cutoff_nm)?;

    let n = input.positions.len();

    // 2. Exclusions
    let pair_map = build_pair_class_map(n, &input.bonds)?;

    // 3. Loop over pairs
    let mut energy_lj: f64 = 0.0;
    let mut energy_coul: f64 = 0.0;
    let mut forces: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0]; n];
    let mut virial: [[f64; 3]; 3] = [[0.0; 3]; 3];

    for i in 0..n {
        for j in (i + 1)..n {
            let class = pair_map.classes[pair_index(n, i, j)];
            if class == PairClass::Excluded {
                continue;
            }

            let (d_min, r) =
                minimum_image_displacement(input.positions[i], input.positions[j], cell)?;

            let sigma_ij = combine_sigma(
                input.rule,
                input.lj_sigma[input.types[i]],
                input.lj_sigma[input.types[j]],
            );
            let eps_ij = combine_eps(input.lj_eps[input.types[i]], input.lj_eps[input.types[j]]);

            // Construct positions so that pair functions see the minimum-image displacement.
            // Pass r_i and r_j' where r_j' = r_i - d_min, so that r_i - r_j' = d_min.
            let r_j_corrected = [
                input.positions[i][0] - d_min[0],
                input.positions[i][1] - d_min[1],
                input.positions[i][2] - d_min[2],
            ];

            let mut f_i_total: [f64; 3] = [0.0, 0.0, 0.0];

            // LJ contribution
            if r < lj_cutoff_nm {
                let e_lj = lj_pair_energy(sigma_ij, eps_ij, input.positions[i], r_j_corrected)?;
                let (f_i, f_j) =
                    lj_pair_forces(sigma_ij, eps_ij, input.positions[i], r_j_corrected)?;

                let (e_lj_scaled, f_i_scaled, f_j_scaled) = if class == PairClass::OneFour {
                    (
                        e_lj * input.lj_scale_14,
                        [
                            f_i[0] * input.lj_scale_14,
                            f_i[1] * input.lj_scale_14,
                            f_i[2] * input.lj_scale_14,
                        ],
                        [
                            f_j[0] * input.lj_scale_14,
                            f_j[1] * input.lj_scale_14,
                            f_j[2] * input.lj_scale_14,
                        ],
                    )
                } else {
                    (e_lj, f_i, f_j)
                };

                energy_lj += e_lj_scaled;
                f_i_total = add(f_i_total, f_i_scaled);
                forces[j] = add(forces[j], f_j_scaled);
            }

            // Coulomb contribution
            let e_coul = coulomb_pair_energy(
                input.charges[i],
                input.charges[j],
                input.positions[i],
                r_j_corrected,
            )?;
            let (f_c_i, f_c_j) = coulomb_pair_forces(
                input.charges[i],
                input.charges[j],
                input.positions[i],
                r_j_corrected,
            )?;

            let (e_coul_scaled, f_c_i_scaled, f_c_j_scaled) = if class == PairClass::OneFour {
                (
                    e_coul * input.coul_scale_14,
                    [
                        f_c_i[0] * input.coul_scale_14,
                        f_c_i[1] * input.coul_scale_14,
                        f_c_i[2] * input.coul_scale_14,
                    ],
                    [
                        f_c_j[0] * input.coul_scale_14,
                        f_c_j[1] * input.coul_scale_14,
                        f_c_j[2] * input.coul_scale_14,
                    ],
                )
            } else {
                (e_coul, f_c_i, f_c_j)
            };

            energy_coul += e_coul_scaled;
            f_i_total = add(f_i_total, f_c_i_scaled);
            forces[i] = add(forces[i], f_i_total);
            forces[j] = add(forces[j], f_c_j_scaled);

            // Virial: W[a][b] += d_min[a] * F_i[b]
            for a in 0..3 {
                for b in 0..3 {
                    virial[a][b] += d_min[a] * f_i_total[b];
                }
            }
        }
    }

    // 4. Return
    Ok(NonbondedPeriodicResult {
        energy_lj_kjmol: energy_lj,
        energy_coulomb_kjmol: energy_coul,
        forces,
        virial,
    })
}

/// Computes the flat index into the pair class vector for pair (i, j) with i < j.
/// Units: none (integer index).
fn pair_index(n: usize, i: usize, j: usize) -> usize {
    let before = if i == 0 {
        0
    } else {
        i * (n - 1) - i * (i - 1) / 2
    };
    before + (j - i - 1)
}
