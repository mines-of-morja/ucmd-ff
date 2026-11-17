//! Tests of the data model and the total evaluation (M5: R5.1-R5.4).
//!
//! The conformance tests follow spec/force-provider.md §6 (R5.3). The systems:
//! - R5.4: the UnitcellMD fixture `water-2`, two flexible SPC/Fw waters (parameters and positions
//!   from spec/fixtures/valid/water-2, pre-approved). It is checked non-periodic, and periodic
//!   in an arbitrary triclinic test cell with Ewald electrostatics.
//! - A test molecule with every bonded form and 1-4 pairs, two copies, non-periodic and
//!   periodic. All its parameters are arbitrary test values, not a physical force field.
//! - The UnitcellMD fixture `nacl-crystal`: rock salt, 2×2×2 conventional cells, a = 0.564 nm.
//!
//! The water-2 reference energies, forces and virial come from an independent Python
//! evaluation. That script reads the fixture's system.json and positions.npy directly and
//! evaluates each term from its formula, with CODATA 2018 k_e. Its forces and virial are
//! central differences (h = 1e-6 nm, strain 1e-6).
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::cell::Cell;
use ucmd_ff::combination_rule::CombinationRule;
use ucmd_ff::constants::COULOMB_KJ_MOL_NM_E2;
use ucmd_ff::errors::UcmdError;
use ucmd_ff::ewald::EwaldParameters;
use ucmd_ff::fd_test::{numerical_forces, numerical_virial, numerical_virial_periodic};
use ucmd_ff::forcefield::{Evaluation, ForceField, PeriodicSettings};
use ucmd_ff::system::{
    Angle, Atom, AtomType, Bond, Component, Dihedral, DihedralForm, ForceFieldParameters, Improper,
    ImproperForm, MoleculeType, System,
};

// ---------------------------------------------------------------- systems

fn atom(name: &str, mass_amu: f64, charge_e: f64, atom_type: usize) -> Atom {
    Atom {
        name: name.to_string(),
        mass_amu,
        charge_e,
        atom_type,
    }
}

/// The fixture water-2 (spec/fixtures/valid/water-2/system.json): SPC/Fw, Wu, Tepper and Voth,
/// J. Chem. Phys. 124, 024503 (2006), as transcribed by the fixture.
fn water2() -> System {
    let bond = |i, j| Bond {
        atoms: [i, j],
        r0_nm: 0.1012,
        k_kjmol_nm2: 443153.3808,
    };
    let water = MoleculeType {
        name: "SPCFw".to_string(),
        formal_charge_e: 0,
        atoms: vec![
            atom("OW", 15.9994, -0.82, 0),
            atom("HW1", 1.008, 0.41, 1),
            atom("HW2", 1.008, 0.41, 1),
        ],
        bonds: vec![bond(0, 1), bond(0, 2)],
        angles: vec![Angle {
            atoms: [1, 0, 2],
            theta0_rad: 113.24_f64.to_radians(),
            k_kjmol_rad2: 317.5656,
        }],
        dihedrals: vec![],
        impropers: vec![],
    };
    System {
        molecule_types: vec![water],
        components: vec![Component {
            molecule_type: 0,
            count: 2,
        }],
        forcefield: ForceFieldParameters {
            combination_rule: CombinationRule::LorentzBerthelot,
            scale_14_lj: 0.5,
            scale_14_coulomb: 0.5,
            atom_types: vec![
                AtomType {
                    name: "OW".to_string(),
                    sigma_nm: 0.3165492,
                    epsilon_kjmol: 0.6502994552,
                },
                AtomType {
                    name: "HW".to_string(),
                    sigma_nm: 0.0,
                    epsilon_kjmol: 0.0,
                },
            ],
        },
    }
}

/// spec/fixtures/valid/water-2/positions.npy (nm).
fn water2_positions() -> Vec<[f64; 3]> {
    vec![
        [
            0.0014879054633459199,
            0.00023628648429753563,
            -0.000789356853724116,
        ],
        [
            0.08393883634235291,
            0.05485979869227756,
            0.0018390991051459836,
        ],
        [
            -0.08620544326799455,
            0.0553171257097342,
            0.001786700948775454,
        ],
        [
            0.29147061232497856,
            0.029893536587790245,
            0.008219852201468345,
        ],
        [
            0.37580409511058643,
            0.059451909021630854,
            0.057423091208704626,
        ],
        [
            0.20354470422243195,
            0.05944960237012237,
            0.05865817528988138,
        ],
    ]
}

/// Arbitrary test molecule: chain 0-1-2-3 with side atoms 4 (on 1) and 5 (on 2); 1-4 pairs
/// 0-3, 0-5, 4-3, 4-5. Every dihedral and improper form appears once or twice.
fn test_molecule() -> MoleculeType {
    let bond = |i, j, r0| Bond {
        atoms: [i, j],
        r0_nm: r0,
        k_kjmol_nm2: 2.0e5,
    };
    let angle = |i, j, k| Angle {
        atoms: [i, j, k],
        theta0_rad: 1.9,
        k_kjmol_rad2: 300.0,
    };
    MoleculeType {
        name: "test".to_string(),
        formal_charge_e: 0,
        atoms: vec![
            atom("C0", 12.0, -0.3, 0),
            atom("C1", 12.0, 0.2, 0),
            atom("C2", 12.0, 0.1, 0),
            atom("C3", 12.0, -0.25, 0),
            atom("H4", 1.0, 0.15, 1),
            atom("H5", 1.0, 0.1, 1),
        ],
        bonds: vec![
            bond(0, 1, 0.15),
            bond(1, 2, 0.153),
            bond(2, 3, 0.152),
            bond(1, 4, 0.11),
            bond(2, 5, 0.109),
        ],
        angles: vec![
            angle(0, 1, 2),
            angle(1, 2, 3),
            angle(0, 1, 4),
            angle(4, 1, 2),
            angle(1, 2, 5),
            angle(5, 2, 3),
        ],
        dihedrals: vec![
            Dihedral {
                atoms: [0, 1, 2, 3],
                form: DihedralForm::Periodic {
                    k_kjmol: 5.0,
                    n: 3,
                    phi0_rad: 0.2,
                },
            },
            Dihedral {
                atoms: [4, 1, 2, 5],
                form: DihedralForm::OplsFourier {
                    v_kjmol: [3.0, -1.5, 2.0, 0.5],
                },
            },
            Dihedral {
                atoms: [0, 1, 2, 5],
                form: DihedralForm::RyckaertBellemans {
                    c_kjmol: [1.0, 2.0, -1.0, 0.5, 0.3, -0.2],
                },
            },
        ],
        impropers: vec![
            Improper {
                atoms: [1, 0, 2, 4],
                form: ImproperForm::Harmonic {
                    xi0_rad: 0.3,
                    k_kjmol_rad2: 40.0,
                },
            },
            Improper {
                atoms: [2, 1, 3, 5],
                form: ImproperForm::Periodic {
                    k_kjmol: 4.0,
                    n: 2,
                    phi0_rad: std::f64::consts::PI,
                },
            },
        ],
    }
}

fn test_system() -> System {
    System {
        molecule_types: vec![test_molecule()],
        components: vec![Component {
            molecule_type: 0,
            count: 2,
        }],
        forcefield: ForceFieldParameters {
            combination_rule: CombinationRule::Geometric,
            scale_14_lj: 0.5,
            scale_14_coulomb: 0.8333333333333334,
            atom_types: vec![
                AtomType {
                    name: "A".to_string(),
                    sigma_nm: 0.3,
                    epsilon_kjmol: 0.5,
                },
                AtomType {
                    name: "B".to_string(),
                    sigma_nm: 0.25,
                    epsilon_kjmol: 0.2,
                },
            ],
        },
    }
}

/// Arbitrary non-planar geometry of two copies of the test molecule (nm).
fn test_positions() -> Vec<[f64; 3]> {
    let one = [
        [0.000, 0.000, 0.000],
        [0.150, 0.020, 0.010],
        [0.210, 0.160, 0.030],
        [0.360, 0.180, 0.090],
        [0.170, -0.060, 0.120],
        [0.190, 0.230, -0.090],
    ];
    let mut p: Vec<[f64; 3]> = one.to_vec();
    for (k, r) in one.iter().enumerate() {
        let wiggle = 0.01 * ((k * 7) % 5) as f64 - 0.02;
        p.push([
            r[0] + 0.45 + wiggle,
            r[1] + 0.35 - wiggle,
            r[2] + 0.30 + 0.5 * wiggle,
        ]);
    }
    p
}

/// Arbitrary triclinic test cell in reduced form (nm); smallest height about 1.08 nm.
fn test_cell() -> Cell {
    Cell {
        ax: 1.3,
        bx: 0.2,
        by: 1.2,
        cx: -0.3,
        cy: 0.1,
        cz: 1.1,
    }
}

fn periodic_settings() -> PeriodicSettings {
    PeriodicSettings {
        lj_cutoff_nm: 0.48,
        ewald: EwaldParameters::from_tolerance(0.9, 1e-14).unwrap(),
    }
}

// ---------------------------------------------------------------- conformance (R5.3)

const STEPS: [f64; 5] = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
const LOWER: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (1, 0), (2, 0), (2, 1)];

/// The provider test kit of spec/force-provider.md §6, items 1-6.
fn conformance(
    ff: &ForceField,
    pos: &[[f64; 3]],
    periodic: Option<(&Cell, &PeriodicSettings)>,
    label: &str,
) {
    let ev = ff.evaluate(pos, periodic).unwrap();
    let energy = |p: &[[f64; 3]]| ff.evaluate(p, periodic).map(|e| e.energy_kjmol);

    // 5. Finiteness and shape.
    assert_eq!(ev.forces.len(), pos.len(), "{label}: forces length");
    assert!(ev.energy_kjmol.is_finite(), "{label}: energy");
    assert!(
        ev.forces.iter().flatten().all(|f| f.is_finite()),
        "{label}: forces finite"
    );
    assert!(
        ev.virial.iter().flatten().all(|w| w.is_finite()),
        "{label}: virial finite"
    );

    // 1. Finite-difference forces: at the best step, relative error <= 1e-5 on |F| > 1e-3.
    let mut best = f64::INFINITY;
    for h in STEPS {
        let num = numerical_forces(&energy, pos, h).unwrap();
        let mut worst = 0.0_f64;
        for (fa, fn_) in ev.forces.iter().zip(&num) {
            for c in 0..3 {
                if fa[c].abs() > 1e-3 {
                    worst = worst.max((fa[c] - fn_[c]).abs() / fa[c].abs());
                }
            }
        }
        best = best.min(worst);
    }
    assert!(
        best <= 1e-5,
        "{label}: finite-difference forces, best relative error {best:e}"
    );

    // 2. Finite-difference virial (relative 1e-5 of the largest component), at the best strain.
    let w_scale = ev
        .virial
        .iter()
        .flatten()
        .fold(0.0_f64, |m, w| m.max(w.abs()));
    assert!(w_scale > 0.0, "{label}: virial is zero");
    let mut best = f64::INFINITY;
    for eps in STEPS {
        let mut worst = 0.0_f64;
        match periodic {
            None => {
                let num = numerical_virial(&energy, pos, eps).unwrap();
                for a in 0..3 {
                    for b in 0..3 {
                        worst = worst.max((ev.virial[a][b] - num[a][b]).abs() / w_scale);
                    }
                }
            }
            Some((cell, settings)) => {
                let e_cell = |p: &[[f64; 3]], c: &Cell| {
                    ff.evaluate(p, Some((c, settings))).map(|e| e.energy_kjmol)
                };
                for (a, b) in LOWER {
                    let num = numerical_virial_periodic(&e_cell, pos, cell, a, b, eps).unwrap();
                    worst = worst.max((ev.virial[a][b] - num).abs() / w_scale);
                }
            }
        }
        best = best.min(worst);
    }
    assert!(
        best <= 1e-5,
        "{label}: finite-difference virial, best relative error {best:e}"
    );

    // 3. Energy terms sum to the energy (relative 1e-10).
    let sum: f64 = ev.energy_terms.iter().map(|t| t.energy_kjmol).sum();
    assert!(
        (sum - ev.energy_kjmol).abs() <= 1e-10 * ev.energy_kjmol.abs().max(1e-300),
        "{label}: terms sum {sum} vs energy {}",
        ev.energy_kjmol
    );

    // 4. Momentum conservation: |sum F| <= 1e-8 sum |F|.
    let total: f64 = ev
        .forces
        .iter()
        .map(|f| (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt())
        .sum();
    let mut s = [0.0_f64; 3];
    for f in &ev.forces {
        for c in 0..3 {
            s[c] += f[c];
        }
    }
    let s_norm = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
    assert!(
        s_norm <= 1e-8 * total,
        "{label}: |sum F| = {s_norm:e}, sum |F| = {total}"
    );

    // 6. Determinism: a second evaluation is bitwise identical.
    let again = ff.evaluate(pos, periodic).unwrap();
    assert_eq!(
        bits(&ev),
        bits(&again),
        "{label}: not bitwise deterministic"
    );
}

fn bits(e: &Evaluation) -> Vec<u64> {
    let mut v = vec![e.energy_kjmol.to_bits()];
    v.extend(e.forces.iter().flatten().map(|x| x.to_bits()));
    v.extend(e.virial.iter().flatten().map(|x| x.to_bits()));
    v.extend(e.energy_terms.iter().map(|t| t.energy_kjmol.to_bits()));
    v
}

#[test]
fn conformance_water2_nonperiodic() {
    let ff = ForceField::new(&water2()).unwrap();
    conformance(&ff, &water2_positions(), None, "water-2");
}

#[test]
fn conformance_water2_periodic_ewald() {
    let ff = ForceField::new(&water2()).unwrap();
    let cell = test_cell();
    conformance(
        &ff,
        &water2_positions(),
        Some((&cell, &periodic_settings())),
        "water-2 periodic",
    );
}

#[test]
fn conformance_test_molecule_nonperiodic() {
    let ff = ForceField::new(&test_system()).unwrap();
    conformance(&ff, &test_positions(), None, "test molecule");
}

#[test]
fn conformance_test_molecule_periodic_ewald() {
    let ff = ForceField::new(&test_system()).unwrap();
    let cell = test_cell();
    conformance(
        &ff,
        &test_positions(),
        Some((&cell, &periodic_settings())),
        "test molecule periodic",
    );
}

// ---------------------------------------------------------------- R5.2 and R5.4 values

#[test]
fn water2_terms_match_independent_evaluation() {
    let ev = ForceField::new(&water2())
        .unwrap()
        .evaluate(&water2_positions(), None)
        .unwrap();
    let names: Vec<&str> = ev.energy_terms.iter().map(|t| t.name).collect();
    assert_eq!(
        names,
        ["bond", "angle", "dihedral", "improper", "lj", "coulomb"]
    );
    let want = [
        ("bond", 6.808423179449083),
        ("angle", 0.11099804198932818),
        ("dihedral", 0.0),
        ("improper", 0.0),
        ("lj", 2.70276461066953),
        ("coulomb", 22.412132333352886),
    ];
    for (t, (name, value)) in ev.energy_terms.iter().zip(want) {
        assert_eq!(t.name, name);
        assert!(
            (t.energy_kjmol - value).abs() <= 1e-10 * value.abs(),
            "{name}: {} vs {value}",
            t.energy_kjmol
        );
    }
    let total = 32.034318165460824;
    assert!((ev.energy_kjmol - total).abs() <= 1e-10 * total);
}

#[test]
fn water2_forces_and_virial_match_independent_reference() {
    let ev = ForceField::new(&water2())
        .unwrap()
        .evaluate(&water2_positions(), None)
        .unwrap();
    let want_f = [
        [-1807.230602274501, 94.9329765234097, 268.4669331642908],
        [402.0786526837128, 424.3359381312928, -558.4895816994617],
        [881.2019737582943, -542.2786902613552, -83.0039597907728],
        [-1417.6336072004858, 979.6046004595382, 1188.1840621867923],
        [
            -207.85477023466115,
            -173.56808417545722,
            -202.26153028346516,
        ],
        [2149.4383532996153, -783.0267407165081, -612.8959235525144],
    ];
    let want_w = [
        [-98.7077961696059, 143.42430056402122, 106.2341862940741],
        [143.42430057112665, -34.28161916119166, -48.108960587711636],
        [106.2341862905214, -48.10896065876591, -39.18647852785284],
    ];
    for i in 0..6 {
        for c in 0..3 {
            let (got, want) = (ev.forces[i][c], want_f[i][c]);
            assert!(
                (got - want).abs() < 1e-6 * 2149.4,
                "F[{i}][{c}]: {got} vs {want}"
            );
        }
    }
    for a in 0..3 {
        for b in 0..3 {
            let (got, want) = (ev.virial[a][b], want_w[a][b]);
            assert!(
                (got - want).abs() < 1e-6 * 143.4,
                "W[{a}][{b}]: {got} vs {want}"
            );
        }
    }
}

#[test]
fn nacl_crystal_gives_the_madelung_constant() {
    // Fixture nacl-crystal: 2x2x2 conventional rock-salt cells, a = 0.564 nm (cell 1.128 nm),
    // charges +-1 e, no Lennard-Jones. E = -32 M k_e / (a/2), M = 1.7475645946 (R4.3).
    let a = 0.564;
    let ion = |name: &str, q: f64| MoleculeType {
        name: name.to_string(),
        formal_charge_e: q as i32,
        atoms: vec![atom(name, 1.0, q, 0)], // the mass does not enter the energy
        bonds: vec![],
        angles: vec![],
        dihedrals: vec![],
        impropers: vec![],
    };
    let system = System {
        molecule_types: vec![ion("Na", 1.0), ion("Cl", -1.0)],
        components: vec![
            Component {
                molecule_type: 0,
                count: 32,
            },
            Component {
                molecule_type: 1,
                count: 32,
            },
        ],
        forcefield: ForceFieldParameters {
            combination_rule: CombinationRule::Geometric,
            scale_14_lj: 1.0,
            scale_14_coulomb: 1.0,
            atom_types: vec![AtomType {
                name: "none".to_string(),
                sigma_nm: 0.0,
                epsilon_kjmol: 0.0,
            }],
        },
    };
    let fcc = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ];
    let mut na = Vec::new();
    let mut cl = Vec::new();
    for ix in 0..2 {
        for iy in 0..2 {
            for iz in 0..2 {
                for s in fcc {
                    let r = [
                        a * (ix as f64 + s[0]),
                        a * (iy as f64 + s[1]),
                        a * (iz as f64 + s[2]),
                    ];
                    na.push(r);
                    cl.push([r[0] + 0.5 * a, r[1], r[2]]);
                }
            }
        }
    }
    let pos: Vec<[f64; 3]> = na.into_iter().chain(cl).collect();
    let cell = Cell {
        ax: 2.0 * a,
        bx: 0.0,
        by: 2.0 * a,
        cx: 0.0,
        cy: 0.0,
        cz: 2.0 * a,
    };
    let settings = PeriodicSettings {
        lj_cutoff_nm: 0.5,
        ewald: EwaldParameters::from_tolerance(0.55, 1e-15).unwrap(),
    };
    let ev = ForceField::new(&system)
        .unwrap()
        .evaluate(&pos, Some((&cell, &settings)))
        .unwrap();
    let m = -ev.energy_kjmol * (0.5 * a) / (32.0 * COULOMB_KJ_MOL_NM_E2);
    assert!(
        ((m - 1.7475645946) / 1.7475645946).abs() < 1e-8,
        "Madelung constant {m}"
    );
    let lj = ev.energy_terms.iter().find(|t| t.name == "lj").unwrap();
    assert_eq!(lj.energy_kjmol, 0.0);
}

// ---------------------------------------------------------------- R5.1 data model

#[test]
fn canonical_order_of_atoms_and_interactions() {
    // Components: water x1, test molecule x2, water x1.
    let w = water2().molecule_types.remove(0);
    let mut sys = test_system();
    sys.molecule_types.push(w);
    sys.components = vec![
        Component {
            molecule_type: 1,
            count: 1,
        },
        Component {
            molecule_type: 0,
            count: 2,
        },
        Component {
            molecule_type: 1,
            count: 1,
        },
    ];
    sys.forcefield.atom_types.push(AtomType {
        name: "HW".to_string(),
        sigma_nm: 0.0,
        epsilon_kjmol: 0.0,
    });
    assert_eq!(sys.atom_count().unwrap(), 3 + 12 + 3);
    let t = sys.topology().unwrap();
    assert_eq!(t.atom_count(), 18);
    assert_eq!(
        t.molecule,
        [0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3]
    );
    assert_eq!(&t.charges[..4], &[-0.82, 0.41, 0.41, -0.3]);
    assert_eq!(t.charges[15], -0.82);
    // Bonds: water (2), test molecule (5) twice, water (2), with global indices.
    let b: Vec<[usize; 2]> = t.bonds.iter().map(|b| b.atoms).collect();
    assert_eq!(b[0], [0, 1]);
    assert_eq!(b[2], [3, 4]);
    assert_eq!(b[7], [9, 10]);
    assert_eq!(b[12], [15, 16]);
    assert_eq!(t.bonds.len(), 2 + 5 + 5 + 2);
    assert_eq!(t.angles[1].atoms, [3, 4, 5]);
    assert_eq!(t.dihedrals.len(), 6);
    assert_eq!(t.dihedrals[3].atoms, [9, 10, 11, 12]);
    assert_eq!(t.impropers[2].atoms, [10, 9, 11, 13]);
}

#[test]
fn invalid_systems_are_errors() {
    let check = |f: &dyn Fn(&mut System), err: UcmdError, what: &str| {
        let mut s = test_system();
        f(&mut s);
        assert_eq!(s.validate(), Err(err), "{what}");
        assert_eq!(ForceField::new(&s).err(), Some(err), "{what}");
    };
    check(
        &|s| s.molecule_types[0].bonds[0].atoms = [0, 6],
        UcmdError::IndexOutOfRange,
        "bond index",
    );
    check(
        &|s| s.molecule_types[0].angles[0].atoms = [0, 1, 0],
        UcmdError::IndexOutOfRange,
        "repeated index",
    );
    check(
        &|s| s.molecule_types[0].atoms[2].atom_type = 2,
        UcmdError::IndexOutOfRange,
        "atom type",
    );
    check(
        &|s| s.components[0].molecule_type = 1,
        UcmdError::IndexOutOfRange,
        "molecule type",
    );
    check(
        &|s| s.components[0].count = 0,
        UcmdError::InvalidParameter,
        "count 0",
    );
    check(
        &|s| s.components.clear(),
        UcmdError::InvalidParameter,
        "no components",
    );
    check(
        &|s| s.molecule_types[0].bonds[1].r0_nm = 0.0,
        UcmdError::InvalidParameter,
        "r0 = 0",
    );
    check(
        &|s| s.molecule_types[0].bonds[1].k_kjmol_nm2 = -1.0,
        UcmdError::InvalidParameter,
        "k < 0",
    );
    check(
        &|s| s.molecule_types[0].angles[1].theta0_rad = 3.2,
        UcmdError::InvalidParameter,
        "theta0 > pi",
    );
    check(
        &|s| s.molecule_types[0].atoms[0].mass_amu = 0.0,
        UcmdError::InvalidParameter,
        "mass 0",
    );
    check(
        &|s| s.molecule_types[0].atoms[0].charge_e = f64::NAN,
        UcmdError::InvalidParameter,
        "charge NaN",
    );
    check(
        &|s| s.forcefield.scale_14_lj = 1.5,
        UcmdError::InvalidParameter,
        "scale > 1",
    );
    check(
        &|s| s.forcefield.atom_types[0].sigma_nm = -0.1,
        UcmdError::InvalidParameter,
        "sigma < 0",
    );
    check(
        &|s| {
            s.molecule_types[0].dihedrals[0].form = DihedralForm::Periodic {
                k_kjmol: 1.0,
                n: 7,
                phi0_rad: 0.0,
            }
        },
        UcmdError::InvalidParameter,
        "n = 7",
    );
    check(
        &|s| {
            s.molecule_types[0].impropers[0].form = ImproperForm::Harmonic {
                xi0_rad: 0.0,
                k_kjmol_rad2: -2.0,
            }
        },
        UcmdError::InvalidParameter,
        "improper k < 0",
    );
    check(
        &|s| s.molecule_types[0].atoms.clear(),
        UcmdError::InvalidParameter,
        "no atoms",
    );
}

#[test]
fn invalid_evaluations_are_errors() {
    let ff = ForceField::new(&water2()).unwrap();
    let pos = water2_positions();
    assert_eq!(
        ff.evaluate(&pos[..5], None).err(),
        Some(UcmdError::IndexOutOfRange)
    );
    let mut bad = pos.clone();
    bad[3][0] = f64::INFINITY;
    assert_eq!(
        ff.evaluate(&bad, None).err(),
        Some(UcmdError::NonFiniteCoordinate)
    );
    // A molecule split across the cell: H moved by a lattice vector away from its O.
    let cell = test_cell();
    let mut split = pos.clone();
    split[1][0] += cell.ax;
    assert_eq!(
        ff.evaluate(&split, Some((&cell, &periodic_settings())))
            .err(),
        Some(UcmdError::MoleculeSplit)
    );
    // LJ cutoff not below half the smallest height.
    let too_long = PeriodicSettings {
        lj_cutoff_nm: 0.6,
        ..periodic_settings()
    };
    assert_eq!(
        ff.evaluate(&pos, Some((&cell, &too_long))).err(),
        Some(UcmdError::InvalidCell)
    );
    // A non-neutral periodic system.
    let mut charged = water2();
    charged.molecule_types[0].atoms[0].charge_e = -0.8;
    let ff_charged = ForceField::new(&charged).unwrap();
    assert!(ff_charged.evaluate(&pos, None).is_ok());
    assert_eq!(
        ff_charged
            .evaluate(&pos, Some((&cell, &periodic_settings())))
            .err(),
        Some(UcmdError::NonNeutralSystem)
    );
}

#[test]
fn periodic_and_nonperiodic_agree_on_the_short_range_terms() {
    // Bonded and LJ terms of water-2 are the same in a large cell (all pairs within the cutoff).
    let ff = ForceField::new(&water2()).unwrap();
    let pos = water2_positions();
    let open = ff.evaluate(&pos, None).unwrap();
    let cell = test_cell();
    let per = ff
        .evaluate(&pos, Some((&cell, &periodic_settings())))
        .unwrap();
    for name in ["bond", "angle", "dihedral", "improper", "lj"] {
        let a = open
            .energy_terms
            .iter()
            .find(|t| t.name == name)
            .unwrap()
            .energy_kjmol;
        let b = per
            .energy_terms
            .iter()
            .find(|t| t.name == name)
            .unwrap()
            .energy_kjmol;
        assert!(
            (a - b).abs() <= 1e-12 * a.abs().max(1.0),
            "{name}: {a} vs {b}"
        );
    }
    let names: Vec<&str> = per.energy_terms.iter().map(|t| t.name).collect();
    assert_eq!(
        names,
        [
            "bond",
            "angle",
            "dihedral",
            "improper",
            "lj",
            "coulomb_real",
            "coulomb_reciprocal",
            "coulomb_self",
            "coulomb_correction"
        ]
    );
}

#[test]
fn test_molecule_terms_match_independent_evaluation() {
    // Reference: every term evaluated from its formula in Python (spec/system-package.md §5.3,
    // IUPAC dihedral by atan2), exclusions and 1-4 pairs (8 of them) from the bond graph. The
    // periodic Coulomb energy is the Python converged Ewald sum with excluded and 1-4 pairs
    // removed, plus s14 k_e q_i q_j / r (minimum image) for the 1-4 pairs.
    let ff = ForceField::new(&test_system()).unwrap();
    let pos = test_positions();
    let ev = ff.evaluate(&pos, None).unwrap();
    let want = [
        ("bond", 1029.6503098107414),
        ("angle", 104.56399447521024),
        ("dihedral", 5.9279377138501115),
        ("improper", 42.44541588683559),
        ("lj", -1.1481105394928768),
        ("coulomb", 3.785832682658132),
    ];
    for (t, (name, value)) in ev.energy_terms.iter().zip(want) {
        assert_eq!(t.name, name);
        assert!(
            (t.energy_kjmol - value).abs() <= 1e-10 * value.abs(),
            "{name}: {} vs {value}",
            t.energy_kjmol
        );
    }
    let cell = test_cell();
    let per = ff
        .evaluate(&pos, Some((&cell, &periodic_settings())))
        .unwrap();
    let term = |name: &str| {
        per.energy_terms
            .iter()
            .find(|t| t.name == name)
            .unwrap()
            .energy_kjmol
    };
    let lj = term("lj");
    assert!(
        (lj - -0.7049477585416924).abs() <= 1e-10 * 0.705,
        "periodic lj {lj}"
    );
    let coulomb: f64 = [
        "coulomb_real",
        "coulomb_reciprocal",
        "coulomb_self",
        "coulomb_correction",
    ]
    .iter()
    .map(|n| term(n))
    .sum();
    assert!(
        (coulomb - 3.350465861155718).abs() <= 1e-9,
        "periodic coulomb {coulomb}"
    );
}
