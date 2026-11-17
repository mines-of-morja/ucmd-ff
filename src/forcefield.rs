//! Total evaluation of a parameterized system: energy, forces and virial as the sum of all
//! terms, with named energy terms that sum to the energy (R5.2).
//!
//! Terms, in this order:
//! - `bond`, `angle`, `dihedral`, `improper`: the bonded terms of the topology. They use the
//!   unwrapped positions directly, because molecules are whole (spec/units.md §3).
//! - `lj`: Lennard-Jones with exclusions and 1-4 scaling. A non-periodic system uses all pairs;
//!   a periodic one uses the minimum image and the plain cutoff of `PeriodicSettings`.
//! - Coulomb, with exclusions and 1-4 scaling:
//!   - non-periodic: `coulomb`, the direct sum over all pairs;
//!   - periodic: the Ewald sum, as `coulomb_real`, `coulomb_reciprocal`, `coulomb_self` and
//!     `coulomb_correction`.
//!
//! The energy is the sum of the terms in this order, so the terms add up to it exactly. The virial
//! is W_ab = −∂E/∂ε_ab (spec/units.md §4). For a periodic system the strain acts on the positions
//! and the lattice vectors. For a non-periodic one it acts on the positions only. Every term's
//! virial is the sum over its atoms of (r_k − r_ref) ⊗ F_k, which is equivalent for
//! translation-invariant terms.
//!
//! Units: positions and cell in nm; energies in kJ/mol; forces in kJ mol⁻¹ nm⁻¹; virial in kJ/mol.

use crate::angle::{angle_energy, angle_forces};
use crate::bond::{bond_energy, bond_forces};
use crate::cell::{Cell, validate_cell};
use crate::errors::UcmdError;
use crate::ewald::{EwaldParameters, ScaledPair, ewald, scaled_pairs_from_bonds};
use crate::harmonic_improper::{harmonic_improper_energy, harmonic_improper_forces};
use crate::minimum_image::perpendicular_heights;
use crate::nonbonded_nonperiodic::{NonbondedInput, evaluate_nonperiodic};
use crate::nonbonded_periodic::evaluate_periodic;
use crate::opls_fourier_dihedral::{opls_fourier_dihedral_energy, opls_fourier_dihedral_forces};
use crate::periodic_dihedral::{periodic_dihedral_energy, periodic_dihedral_forces};
use crate::ryckaert_bellemans_dihedral::{
    ryckaert_bellemans_dihedral_energy, ryckaert_bellemans_dihedral_forces,
};
use crate::system::{DihedralForm, ImproperForm, System, Topology};

/// Run settings of a periodic evaluation (cutoffs and electrostatics are run settings, not part
/// of the system, spec/system-package.md §5.2).
/// Units: lj_cutoff_nm in nm; ewald as `EwaldParameters`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PeriodicSettings {
    /// Plain Lennard-Jones cutoff, below half the smallest cell height (R3.5). Units: nm.
    pub lj_cutoff_nm: f64,
    /// Parameters of the Ewald sum for the Coulomb interaction. Units: as `EwaldParameters`.
    pub ewald: EwaldParameters,
}

/// A named part of the energy. Units: energy_kjmol in kJ/mol.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnergyTerm {
    /// Name of the term (see the module documentation). Units: none.
    pub name: &'static str,
    /// Energy of the term. Units: kJ/mol.
    pub energy_kjmol: f64,
}

/// Result of one evaluation.
/// Units: energy_kjmol in kJ/mol; forces in kJ mol⁻¹ nm⁻¹; virial in kJ/mol.
#[derive(Debug, Clone, PartialEq)]
pub struct Evaluation {
    /// Potential energy, the sum of `energy_terms`. Units: kJ/mol.
    pub energy_kjmol: f64,
    /// Force on every atom in canonical order, F = −∂E/∂r. Units: kJ mol⁻¹ nm⁻¹.
    pub forces: Vec<[f64; 3]>,
    /// Virial W_ab = −∂E/∂ε_ab. Units: kJ/mol.
    pub virial: [[f64; 3]; 3],
    /// Named energy terms. Units: kJ/mol.
    pub energy_terms: Vec<EnergyTerm>,
}

/// A force field ready to evaluate one system: the expanded topology and the pair lists derived
/// from its bonds. Units: none.
#[derive(Debug, Clone, PartialEq)]
pub struct ForceField {
    topology: Topology,
    bond_pairs: Vec<(usize, usize)>,
    scaled_pairs: Vec<ScaledPair>,
}

/// Forces and virial being summed over terms.
struct Accumulator<'a> {
    positions: &'a [[f64; 3]],
    forces: Vec<[f64; 3]>,
    virial: [[f64; 3]; 3],
}

impl Accumulator<'_> {
    /// Adds the forces of one bonded term on its atoms and its virial Σ (r_k − r_ref) ⊗ F_k.
    fn add_term(&mut self, atoms: &[usize], forces: &[[f64; 3]]) {
        let r_ref = self.positions[atoms[0]];
        for (&a, f) in atoms.iter().zip(forces) {
            let r = self.positions[a];
            for c in 0..3 {
                self.forces[a][c] += f[c];
            }
            for x in 0..3 {
                for y in 0..3 {
                    self.virial[x][y] += (r[x] - r_ref[x]) * f[y];
                }
            }
        }
    }

    /// Adds forces and virial that are already complete for the whole system.
    fn add_all(&mut self, forces: &[[f64; 3]], virial: &[[f64; 3]; 3]) {
        for (t, f) in self.forces.iter_mut().zip(forces) {
            for c in 0..3 {
                t[c] += f[c];
            }
        }
        for x in 0..3 {
            for y in 0..3 {
                self.virial[x][y] += virial[x][y];
            }
        }
    }
}

impl ForceField {
    /// Validates the system and prepares its evaluation.
    /// Units: none. Errors: as `System::topology`.
    pub fn new(system: &System) -> Result<Self, UcmdError> {
        let topology = system.topology()?;
        let bond_pairs: Vec<(usize, usize)> = topology
            .bonds
            .iter()
            .map(|b| (b.atoms[0], b.atoms[1]))
            .collect();
        let scaled_pairs = scaled_pairs_from_bonds(
            topology.atom_count(),
            &bond_pairs,
            topology.forcefield.scale_14_coulomb,
        )?;
        Ok(ForceField {
            topology,
            bond_pairs,
            scaled_pairs,
        })
    }

    /// The expanded topology (canonical atom order). Units: as `Topology`.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// Energy, forces, virial and energy terms of a configuration. `periodic` is `None` for a
    /// non-periodic system, or the cell with the run settings.
    /// Units: positions and cell in nm; results as `Evaluation`.
    /// Errors: `IndexOutOfRange` (wrong number of positions), `NonFiniteCoordinate`,
    /// `InvalidCell` (invalid cell, or an LJ cutoff ≥ h_min/2), `MoleculeSplit` (a bond longer
    /// than half the smallest cell height in a periodic system), `NonNeutralSystem` (periodic),
    /// `InvalidParameter`, and the errors of the terms (`CoincidentAtoms`, collinear dihedrals).
    pub fn evaluate(
        &self,
        positions: &[[f64; 3]],
        periodic: Option<(&Cell, &PeriodicSettings)>,
    ) -> Result<Evaluation, UcmdError> {
        let t = &self.topology;
        let n = t.atom_count();
        if positions.len() != n {
            return Err(UcmdError::IndexOutOfRange);
        }
        if positions.iter().flatten().any(|c| !c.is_finite()) {
            return Err(UcmdError::NonFiniteCoordinate);
        }
        if let Some((cell, _)) = periodic {
            validate_cell(cell)?;
            let h = perpendicular_heights(cell)?;
            let half = 0.5 * h[0].min(h[1]).min(h[2]);
            for b in &t.bonds {
                let (p, q) = (positions[b.atoms[0]], positions[b.atoms[1]]);
                let d2 = (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2);
                if d2.sqrt() >= half {
                    return Err(UcmdError::MoleculeSplit);
                }
            }
        }
        let mut acc = Accumulator {
            positions,
            forces: vec![[0.0; 3]; n],
            virial: [[0.0; 3]; 3],
        };
        let mut terms = Vec::new();
        let p = |i: usize| positions[i];

        let mut e = 0.0;
        for b in &t.bonds {
            let [i, j] = b.atoms;
            e += bond_energy(p(i), p(j), b.r0_nm, b.k_kjmol_nm2)?;
            let (fi, fj) = bond_forces(p(i), p(j), b.r0_nm, b.k_kjmol_nm2)?;
            acc.add_term(&b.atoms, &[fi, fj]);
        }
        terms.push(EnergyTerm {
            name: "bond",
            energy_kjmol: e,
        });

        let mut e = 0.0;
        for a in &t.angles {
            let [i, j, k] = a.atoms;
            e += angle_energy(p(i), p(j), p(k), a.theta0_rad, a.k_kjmol_rad2)?;
            let (fi, fj, fk) = angle_forces(p(i), p(j), p(k), a.theta0_rad, a.k_kjmol_rad2)?;
            acc.add_term(&a.atoms, &[fi, fj, fk]);
        }
        terms.push(EnergyTerm {
            name: "angle",
            energy_kjmol: e,
        });

        let mut e = 0.0;
        for d in &t.dihedrals {
            let [i, j, k, l] = d.atoms;
            let (ri, rj, rk, rl) = (p(i), p(j), p(k), p(l));
            let (energy, f) = match d.form {
                DihedralForm::Periodic {
                    k_kjmol,
                    n,
                    phi0_rad,
                } => (
                    periodic_dihedral_energy(ri, rj, rk, rl, n, phi0_rad, k_kjmol)?,
                    periodic_dihedral_forces(ri, rj, rk, rl, n, phi0_rad, k_kjmol)?,
                ),
                DihedralForm::OplsFourier { v_kjmol: v } => (
                    opls_fourier_dihedral_energy(ri, rj, rk, rl, v[0], v[1], v[2], v[3])?,
                    opls_fourier_dihedral_forces(ri, rj, rk, rl, v[0], v[1], v[2], v[3])?,
                ),
                DihedralForm::RyckaertBellemans { c_kjmol } => (
                    ryckaert_bellemans_dihedral_energy(ri, rj, rk, rl, c_kjmol)?,
                    ryckaert_bellemans_dihedral_forces(ri, rj, rk, rl, c_kjmol)?,
                ),
            };
            e += energy;
            acc.add_term(&d.atoms, &[f.0, f.1, f.2, f.3]);
        }
        terms.push(EnergyTerm {
            name: "dihedral",
            energy_kjmol: e,
        });

        let mut e = 0.0;
        for im in &t.impropers {
            let [i, j, k, l] = im.atoms;
            let (ri, rj, rk, rl) = (p(i), p(j), p(k), p(l));
            let (energy, f) = match im.form {
                ImproperForm::Periodic {
                    k_kjmol,
                    n,
                    phi0_rad,
                } => (
                    periodic_dihedral_energy(ri, rj, rk, rl, n, phi0_rad, k_kjmol)?,
                    periodic_dihedral_forces(ri, rj, rk, rl, n, phi0_rad, k_kjmol)?,
                ),
                ImproperForm::Harmonic {
                    xi0_rad,
                    k_kjmol_rad2,
                } => (
                    harmonic_improper_energy(ri, rj, rk, rl, xi0_rad, k_kjmol_rad2)?,
                    harmonic_improper_forces(ri, rj, rk, rl, xi0_rad, k_kjmol_rad2)?,
                ),
            };
            e += energy;
            acc.add_term(&im.atoms, &[f.0, f.1, f.2, f.3]);
        }
        terms.push(EnergyTerm {
            name: "improper",
            energy_kjmol: e,
        });

        let ff = &t.forcefield;
        let mut input = NonbondedInput {
            positions: positions.to_vec(),
            types: t.atom_types.clone(),
            charges: t.charges.clone(),
            lj_sigma: ff.atom_types.iter().map(|a| a.sigma_nm).collect(),
            lj_eps: ff.atom_types.iter().map(|a| a.epsilon_kjmol).collect(),
            rule: ff.combination_rule,
            bonds: self.bond_pairs.clone(),
            lj_scale_14: ff.scale_14_lj,
            coul_scale_14: ff.scale_14_coulomb,
        };
        match periodic {
            None => {
                let r = evaluate_nonperiodic(&input)?;
                terms.push(EnergyTerm {
                    name: "lj",
                    energy_kjmol: r.energy_lj,
                });
                terms.push(EnergyTerm {
                    name: "coulomb",
                    energy_kjmol: r.energy_coulomb,
                });
                acc.add_all(&r.forces, &r.virial);
            }
            Some((cell, settings)) => {
                // Lennard-Jones only: the Coulomb part of the pair evaluation is switched off.
                input.charges = vec![0.0; n];
                let lj = evaluate_periodic(&input, cell, settings.lj_cutoff_nm)?;
                terms.push(EnergyTerm {
                    name: "lj",
                    energy_kjmol: lj.energy_lj_kjmol,
                });
                acc.add_all(&lj.forces, &lj.virial);
                let ew = ewald(
                    positions,
                    &t.charges,
                    cell,
                    &settings.ewald,
                    &self.scaled_pairs,
                )?;
                for (name, part) in [
                    ("coulomb_real", &ew.real),
                    ("coulomb_reciprocal", &ew.reciprocal),
                    ("coulomb_self", &ew.self_term),
                    ("coulomb_correction", &ew.correction),
                ] {
                    terms.push(EnergyTerm {
                        name,
                        energy_kjmol: part.energy_kjmol,
                    });
                }
                acc.add_all(&ew.total.forces, &ew.total.virial);
            }
        }
        let energy_kjmol = terms.iter().map(|t| t.energy_kjmol).sum();
        Ok(Evaluation {
            energy_kjmol,
            forces: acc.forces,
            virial: acc.virial,
            energy_terms: terms,
        })
    }
}
