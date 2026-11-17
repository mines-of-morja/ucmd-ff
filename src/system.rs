//! Data model of a parameterized system (R5.1, spec/system-package.md §3-§5).
//!
//! A `System` holds molecule types (atoms with mass, charge and atom type; bonds, angles,
//! dihedrals and impropers with their functional form and parameters), the components (molecule
//! type × count), and the force field (atom types with σ and ε, the combination rule and the 1-4
//! scale factors). `System::topology` expands it into one list of atoms and interactions in the
//! canonical order of spec/system-package.md §4: components in list order, then instances
//! 0…count−1, then atoms in molecule-type order.
//!
//! Angles are radians in this data model (rule 1); the package files give degrees, so a reader
//! converts them (`f64::to_radians`). Constraints of the package are not part of the energy and
//! are not represented here. Reading `system.json` is the engine's task (ucmd-io).
//!
//! Units: masses in amu; charges in e; lengths in nm; angles in rad; energies in kJ/mol; force
//! constants in kJ mol⁻¹ nm⁻² (bonds) and kJ mol⁻¹ rad⁻² (angles, harmonic impropers).

use crate::combination_rule::CombinationRule;
use crate::errors::UcmdError;

/// A Lennard-Jones atom type (spec/system-package.md §5).
/// Units: sigma_nm in nm; epsilon_kjmol in kJ/mol.
#[derive(Debug, Clone, PartialEq)]
pub struct AtomType {
    /// Name of the type (a key of `atom_types` in the package). Units: none.
    pub name: String,
    /// σ of the type, ≥ 0. Units: nm.
    pub sigma_nm: f64,
    /// ε of the type, ≥ 0. Units: kJ/mol.
    pub epsilon_kjmol: f64,
}

/// An atom of a molecule type.
/// Units: mass_amu in amu; charge_e in e.
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    /// Name of the atom within its molecule type. Units: none.
    pub name: String,
    /// Mass, > 0. Units: amu.
    pub mass_amu: f64,
    /// Partial charge. Units: e.
    pub charge_e: f64,
    /// Index of the atom type in `ForceFieldParameters::atom_types`. Units: none.
    pub atom_type: usize,
}

/// Harmonic bond E = ½ k (r − r0)², atoms indexed within the molecule type.
/// Units: r0_nm in nm (> 0); k_kjmol_nm2 in kJ mol⁻¹ nm⁻² (≥ 0).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bond {
    /// The two atoms. Units: none.
    pub atoms: [usize; 2],
    /// Equilibrium length. Units: nm.
    pub r0_nm: f64,
    /// Force constant. Units: kJ mol⁻¹ nm⁻².
    pub k_kjmol_nm2: f64,
}

/// Harmonic angle E = ½ k (θ − θ0)²; `atoms[1]` is the vertex.
/// Units: theta0_rad in rad, in [0, π]; k_kjmol_rad2 in kJ mol⁻¹ rad⁻² (≥ 0).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Angle {
    /// The three atoms i, j (vertex), k. Units: none.
    pub atoms: [usize; 3],
    /// Equilibrium angle. Units: rad.
    pub theta0_rad: f64,
    /// Force constant. Units: kJ mol⁻¹ rad⁻².
    pub k_kjmol_rad2: f64,
}

/// Functional form of a dihedral (spec/system-package.md §5.3).
/// Units: energies in kJ/mol; angles in rad.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DihedralForm {
    /// E = k [1 + cos(nφ − φ0)], n an integer in 0…6. Units: k in kJ/mol; phi0 in rad.
    Periodic {
        /// Barrier constant. Units: kJ/mol.
        k_kjmol: f64,
        /// Multiplicity, 0…6. Units: none.
        n: i32,
        /// Phase. Units: rad.
        phi0_rad: f64,
    },
    /// E = ½[V1(1+cos φ) + V2(1−cos 2φ) + V3(1+cos 3φ) + V4(1−cos 4φ)]. Units: V in kJ/mol.
    OplsFourier {
        /// V1…V4. Units: kJ/mol.
        v_kjmol: [f64; 4],
    },
    /// E = Σₙ₌₀⁵ Cₙ cosⁿψ, ψ = φ − 180°. Units: C in kJ/mol.
    RyckaertBellemans {
        /// C0…C5. Units: kJ/mol.
        c_kjmol: [f64; 6],
    },
}

/// Functional form of an improper (spec/system-package.md §5.3).
/// Units: energies in kJ/mol; angles in rad.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ImproperForm {
    /// E = k [1 + cos(nξ − ξ0)], as the periodic dihedral. Units: k in kJ/mol; phi0 in rad.
    Periodic {
        /// Barrier constant. Units: kJ/mol.
        k_kjmol: f64,
        /// Multiplicity, 0…6. Units: none.
        n: i32,
        /// Phase. Units: rad.
        phi0_rad: f64,
    },
    /// E = ½ k (ξ − ξ0)², with ξ − ξ0 wrapped to (−π, π]. Units: xi0 in rad; k in kJ mol⁻¹ rad⁻².
    Harmonic {
        /// Equilibrium angle. Units: rad.
        xi0_rad: f64,
        /// Force constant, ≥ 0. Units: kJ mol⁻¹ rad⁻².
        k_kjmol_rad2: f64,
    },
}

/// A dihedral i-j-k-l with its form. Units: as `DihedralForm`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dihedral {
    /// The four atoms. Units: none.
    pub atoms: [usize; 4],
    /// Form and parameters. Units: as `DihedralForm`.
    pub form: DihedralForm,
}

/// An improper i-j-k-l (the angle of the planes (i,j,k) and (j,k,l)) with its form.
/// Units: as `ImproperForm`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Improper {
    /// The four atoms. Units: none.
    pub atoms: [usize; 4],
    /// Form and parameters. Units: as `ImproperForm`.
    pub form: ImproperForm,
}

/// A molecule type: atoms and bonded interactions, indices 0-based within the type.
/// Units: as the member types.
#[derive(Debug, Clone, PartialEq)]
pub struct MoleculeType {
    /// Name. Units: none.
    pub name: String,
    /// Ionic charge number of the molecule (independent of the partial charges). Units: e.
    pub formal_charge_e: i32,
    /// Atoms, at least one. Units: as `Atom`.
    pub atoms: Vec<Atom>,
    /// Bonds; they also define the exclusions (spec/system-package.md §5.2). Units: as `Bond`.
    pub bonds: Vec<Bond>,
    /// Angles. Units: as `Angle`.
    pub angles: Vec<Angle>,
    /// Dihedrals. Units: as `Dihedral`.
    pub dihedrals: Vec<Dihedral>,
    /// Impropers. Units: as `Improper`.
    pub impropers: Vec<Improper>,
}

/// A component: `count` instances of a molecule type. Units: none.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Component {
    /// Index into `System::molecule_types`. Units: none.
    pub molecule_type: usize,
    /// Number of instances, ≥ 1. Units: none.
    pub count: usize,
}

/// Force-field parameters of the package (spec/system-package.md §5).
/// Units: scale factors dimensionless, in [0, 1]; atom types as `AtomType`.
#[derive(Debug, Clone, PartialEq)]
pub struct ForceFieldParameters {
    /// Combination rule for σ_ij (ε_ij is always the geometric mean). Units: none.
    pub combination_rule: CombinationRule,
    /// Factor of the Lennard-Jones interaction of 1-4 pairs. Units: none.
    pub scale_14_lj: f64,
    /// Factor of the Coulomb interaction of 1-4 pairs. Units: none.
    pub scale_14_coulomb: f64,
    /// Atom types. Units: as `AtomType`.
    pub atom_types: Vec<AtomType>,
}

/// A parameterized system: molecule types, components and force field (R5.1).
/// Units: as the member types.
#[derive(Debug, Clone, PartialEq)]
pub struct System {
    /// Molecule types. Units: as `MoleculeType`.
    pub molecule_types: Vec<MoleculeType>,
    /// Components in canonical order. Units: as `Component`.
    pub components: Vec<Component>,
    /// Force field. Units: as `ForceFieldParameters`.
    pub forcefield: ForceFieldParameters,
}

/// The system expanded into global atoms and interactions in canonical order. Interaction
/// indices are global atom indices.
/// Units: as the member types.
#[derive(Debug, Clone, PartialEq)]
pub struct Topology {
    /// Mass of every atom. Units: amu.
    pub masses: Vec<f64>,
    /// Charge of every atom. Units: e.
    pub charges: Vec<f64>,
    /// Atom-type index of every atom. Units: none.
    pub atom_types: Vec<usize>,
    /// Global molecule index of every atom. Units: none.
    pub molecule: Vec<usize>,
    /// All bonds. Units: as `Bond`.
    pub bonds: Vec<Bond>,
    /// All angles. Units: as `Angle`.
    pub angles: Vec<Angle>,
    /// All dihedrals. Units: as `Dihedral`.
    pub dihedrals: Vec<Dihedral>,
    /// All impropers. Units: as `Improper`.
    pub impropers: Vec<Improper>,
    /// Force-field parameters. Units: as `ForceFieldParameters`.
    pub forcefield: ForceFieldParameters,
}

impl Topology {
    /// Number of atoms. Units: none.
    pub fn atom_count(&self) -> usize {
        self.masses.len()
    }
}

fn finite(v: f64) -> Result<f64, UcmdError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(UcmdError::InvalidParameter)
    }
}

fn non_negative(v: f64) -> Result<(), UcmdError> {
    if finite(v)? >= 0.0 {
        Ok(())
    } else {
        Err(UcmdError::InvalidParameter)
    }
}

fn unit_interval(v: f64) -> Result<(), UcmdError> {
    if (0.0..=1.0).contains(&finite(v)?) {
        Ok(())
    } else {
        Err(UcmdError::InvalidParameter)
    }
}

fn multiplicity(n: i32) -> Result<(), UcmdError> {
    if (0..=6).contains(&n) {
        Ok(())
    } else {
        Err(UcmdError::InvalidParameter)
    }
}

/// Atom indices of one interaction: all below n and pairwise different (E_INDEX).
fn check_indices(atoms: &[usize], n: usize) -> Result<(), UcmdError> {
    for (k, &a) in atoms.iter().enumerate() {
        if a >= n || atoms[..k].contains(&a) {
            return Err(UcmdError::IndexOutOfRange);
        }
    }
    Ok(())
}

impl MoleculeType {
    /// Checks the molecule type against the rules of spec/system-package.md §3 and §5.3.
    /// Units: none. Errors: `IndexOutOfRange` (an interaction index out of range or repeated,
    /// an unknown atom type), `InvalidParameter` (no atoms, a mass ≤ 0, a non-finite charge, a
    /// parameter out of range).
    pub fn validate(&self, n_atom_types: usize) -> Result<(), UcmdError> {
        let n = self.atoms.len();
        if n == 0 {
            return Err(UcmdError::InvalidParameter);
        }
        for a in &self.atoms {
            if finite(a.mass_amu)? <= 0.0 {
                return Err(UcmdError::InvalidParameter);
            }
            finite(a.charge_e)?;
            if a.atom_type >= n_atom_types {
                return Err(UcmdError::IndexOutOfRange);
            }
        }
        for b in &self.bonds {
            check_indices(&b.atoms, n)?;
            if finite(b.r0_nm)? <= 0.0 {
                return Err(UcmdError::InvalidParameter);
            }
            non_negative(b.k_kjmol_nm2)?;
        }
        for a in &self.angles {
            check_indices(&a.atoms, n)?;
            if !(0.0..=std::f64::consts::PI).contains(&finite(a.theta0_rad)?) {
                return Err(UcmdError::InvalidParameter);
            }
            non_negative(a.k_kjmol_rad2)?;
        }
        for d in &self.dihedrals {
            check_indices(&d.atoms, n)?;
            match d.form {
                DihedralForm::Periodic {
                    k_kjmol,
                    n,
                    phi0_rad,
                } => {
                    finite(k_kjmol)?;
                    finite(phi0_rad)?;
                    multiplicity(n)?;
                }
                DihedralForm::OplsFourier { v_kjmol } => {
                    for v in v_kjmol {
                        finite(v)?;
                    }
                }
                DihedralForm::RyckaertBellemans { c_kjmol } => {
                    for c in c_kjmol {
                        finite(c)?;
                    }
                }
            }
        }
        for im in &self.impropers {
            check_indices(&im.atoms, n)?;
            match im.form {
                ImproperForm::Periodic {
                    k_kjmol,
                    n,
                    phi0_rad,
                } => {
                    finite(k_kjmol)?;
                    finite(phi0_rad)?;
                    multiplicity(n)?;
                }
                ImproperForm::Harmonic {
                    xi0_rad,
                    k_kjmol_rad2,
                } => {
                    finite(xi0_rad)?;
                    non_negative(k_kjmol_rad2)?;
                }
            }
        }
        Ok(())
    }
}

impl System {
    /// Checks the whole system: every molecule type, every component (known type, count ≥ 1),
    /// every atom type (σ, ε finite and ≥ 0) and the 1-4 scale factors (in [0, 1]).
    /// Units: none. Errors: `IndexOutOfRange`, `InvalidParameter`.
    pub fn validate(&self) -> Result<(), UcmdError> {
        let ff = &self.forcefield;
        for t in &ff.atom_types {
            non_negative(t.sigma_nm)?;
            non_negative(t.epsilon_kjmol)?;
        }
        unit_interval(ff.scale_14_lj)?;
        unit_interval(ff.scale_14_coulomb)?;
        for m in &self.molecule_types {
            m.validate(ff.atom_types.len())?;
        }
        if self.components.is_empty() {
            return Err(UcmdError::InvalidParameter);
        }
        for c in &self.components {
            if c.molecule_type >= self.molecule_types.len() {
                return Err(UcmdError::IndexOutOfRange);
            }
            if c.count == 0 {
                return Err(UcmdError::InvalidParameter);
            }
        }
        Ok(())
    }

    /// Number of atoms: Σ over components of count × atoms per molecule type.
    /// Units: none. Errors: `IndexOutOfRange` for an unknown molecule type.
    pub fn atom_count(&self) -> Result<usize, UcmdError> {
        let mut n = 0;
        for c in &self.components {
            let m = self
                .molecule_types
                .get(c.molecule_type)
                .ok_or(UcmdError::IndexOutOfRange)?;
            n += c.count * m.atoms.len();
        }
        Ok(n)
    }

    /// Validates the system and expands it into global atoms and interactions in canonical
    /// order (spec/system-package.md §4).
    /// Units: as `Topology`. Errors: as `System::validate`.
    pub fn topology(&self) -> Result<Topology, UcmdError> {
        self.validate()?;
        let mut t = Topology {
            masses: Vec::new(),
            charges: Vec::new(),
            atom_types: Vec::new(),
            molecule: Vec::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
            dihedrals: Vec::new(),
            impropers: Vec::new(),
            forcefield: self.forcefield.clone(),
        };
        let mut molecule = 0;
        for c in &self.components {
            let m = &self.molecule_types[c.molecule_type];
            for _ in 0..c.count {
                let off = t.masses.len();
                for a in &m.atoms {
                    t.masses.push(a.mass_amu);
                    t.charges.push(a.charge_e);
                    t.atom_types.push(a.atom_type);
                    t.molecule.push(molecule);
                }
                let shift = |idx: &[usize]| -> Vec<usize> { idx.iter().map(|i| i + off).collect() };
                for b in &m.bonds {
                    let g = shift(&b.atoms);
                    t.bonds.push(Bond {
                        atoms: [g[0], g[1]],
                        ..*b
                    });
                }
                for a in &m.angles {
                    let g = shift(&a.atoms);
                    t.angles.push(Angle {
                        atoms: [g[0], g[1], g[2]],
                        ..*a
                    });
                }
                for d in &m.dihedrals {
                    let g = shift(&d.atoms);
                    t.dihedrals.push(Dihedral {
                        atoms: [g[0], g[1], g[2], g[3]],
                        ..*d
                    });
                }
                for im in &m.impropers {
                    let g = shift(&im.atoms);
                    t.impropers.push(Improper {
                        atoms: [g[0], g[1], g[2], g[3]],
                        ..*im
                    });
                }
                molecule += 1;
            }
        }
        Ok(t)
    }
}
