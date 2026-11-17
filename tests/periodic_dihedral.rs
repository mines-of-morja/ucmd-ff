//! Tests of the `periodic_dihedral` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::errors::UcmdError;
use ucmd_ff::periodic_dihedral::*;

// corrected after review: IUPAC dihedral angle for this geometry is +pi/2, not -pi/2
#[test]
fn test_energy_at_equilibrium() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];
    let r_k: [f64; 3] = [1.0, 1.0, 0.0];
    let r_l: [f64; 3] = [1.0, 1.0, 1.0];
    let n: i32 = 1;
    let phi0: f64 = std::f64::consts::FRAC_PI_2;
    let k: f64 = 10.0;

    let got = periodic_dihedral_energy(r_i, r_j, r_k, r_l, n, phi0, k).unwrap();
    let expected = 2.0 * k;
    let tol = 1e-10;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

// corrected after review: IUPAC dihedral angle for this geometry is +pi/2, so antiphase requires phi0 = phi - pi = -pi/2
#[test]
fn test_energy_at_antiphase() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];
    let r_k: [f64; 3] = [1.0, 1.0, 0.0];
    let r_l: [f64; 3] = [1.0, 1.0, 1.0];
    let n: i32 = 1;
    let phi0: f64 = -std::f64::consts::FRAC_PI_2;
    let k: f64 = 10.0;

    let got = periodic_dihedral_energy(r_i, r_j, r_k, r_l, n, phi0, k).unwrap();
    let expected = 0.0;
    let tol = 1e-10;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

#[test]
fn test_forces_finite_difference() {
    let r_i: [f64; 3] = [0.1, 0.2, 0.3];
    let r_j: [f64; 3] = [1.1, 0.2, 0.3];
    let r_k: [f64; 3] = [1.1, 1.2, 0.3];
    let r_l: [f64; 3] = [1.1, 1.2, 1.3];
    let n: i32 = 2;
    let phi0: f64 = 0.5;
    let k: f64 = 5.0;
    let h: f64 = 1e-6;

    let (f_i, f_j, f_k, f_l) = periodic_dihedral_forces(r_i, r_j, r_k, r_l, n, phi0, k).unwrap();

    let atoms = [r_i, r_j, r_k, r_l];
    let forces = [f_i, f_j, f_k, f_l];

    for a in 0..4 {
        for c in 0..3 {
            let mut r_plus = atoms[a];
            let mut r_minus = atoms[a];
            r_plus[c] += h;
            r_minus[c] -= h;

            let mut atoms_plus = atoms;
            let mut atoms_minus = atoms;
            atoms_plus[a] = r_plus;
            atoms_minus[a] = r_minus;

            let e_plus = periodic_dihedral_energy(
                atoms_plus[0],
                atoms_plus[1],
                atoms_plus[2],
                atoms_plus[3],
                n,
                phi0,
                k,
            )
            .unwrap();
            let e_minus = periodic_dihedral_energy(
                atoms_minus[0],
                atoms_minus[1],
                atoms_minus[2],
                atoms_minus[3],
                n,
                phi0,
                k,
            )
            .unwrap();

            let expected = -(e_plus - e_minus) / (2.0 * h);
            let got = forces[a][c];
            let tol = 1e-8;
            assert!(
                (got - expected).abs() <= tol,
                "atom {a}, coord {c}: got {got}, expected {expected}"
            );
        }
    }
}

#[test]
fn test_virial_finite_difference() {
    let r_i: [f64; 3] = [0.1, 0.2, 0.3];
    let r_j: [f64; 3] = [1.1, 0.2, 0.3];
    let r_k: [f64; 3] = [1.1, 1.2, 0.3];
    let r_l: [f64; 3] = [1.1, 1.2, 1.3];
    let n: i32 = 2;
    let phi0: f64 = 0.5;
    let k: f64 = 5.0;
    let eps: f64 = 1e-6;

    let got_virial = periodic_dihedral_virial(r_i, r_j, r_k, r_l, n, phi0, k).unwrap();

    let atoms = [r_i, r_j, r_k, r_l];

    for a in 0..3 {
        for b in 0..3 {
            let mut atoms_plus = atoms;
            let mut atoms_minus = atoms;

            for m in 0..4 {
                atoms_plus[m][b] += eps * atoms[m][a];
                atoms_minus[m][b] -= eps * atoms[m][a];
            }

            let e_plus = periodic_dihedral_energy(
                atoms_plus[0],
                atoms_plus[1],
                atoms_plus[2],
                atoms_plus[3],
                n,
                phi0,
                k,
            )
            .unwrap();
            let e_minus = periodic_dihedral_energy(
                atoms_minus[0],
                atoms_minus[1],
                atoms_minus[2],
                atoms_minus[3],
                n,
                phi0,
                k,
            )
            .unwrap();

            let expected = -(e_plus - e_minus) / (2.0 * eps);
            let got = got_virial[a][b];
            let tol = 1e-8;
            assert!(
                (got - expected).abs() <= tol,
                "virial[{a}][{b}]: got {got}, expected {expected}"
            );
        }
    }
}

// corrected after review: IUPAC dihedral angle for this geometry is +pi/2, so phi0 must be +pi/2 for equilibrium
#[test]
fn test_zero_force_at_equilibrium() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];
    let r_k: [f64; 3] = [1.0, 1.0, 0.0];
    let r_l: [f64; 3] = [1.0, 1.0, 1.0];
    let n: i32 = 1;
    let phi0: f64 = std::f64::consts::FRAC_PI_2;
    let k: f64 = 10.0;

    let (f_i, f_j, f_k, f_l) = periodic_dihedral_forces(r_i, r_j, r_k, r_l, n, phi0, k).unwrap();
    let forces = [f_i, f_j, f_k, f_l];

    for (a, force) in forces.iter().enumerate() {
        for (c, &got) in force.iter().enumerate() {
            let expected = 0.0;
            let tol = 1e-10;
            assert!(
                (got - expected).abs() <= tol,
                "atom {a}, coord {c}: got {got}, expected {expected}"
            );
        }
    }
}

#[test]
fn test_translation_invariance() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];
    let r_k: [f64; 3] = [1.0, 1.0, 0.0];
    let r_l: [f64; 3] = [1.0, 1.0, 1.0];
    let n: i32 = 1;
    let phi0: f64 = 0.0;
    let k: f64 = 1.0;

    let e_base = periodic_dihedral_energy(r_i, r_j, r_k, r_l, n, phi0, k).unwrap();

    let shift: [f64; 3] = [5.0, 3.0, 2.0];
    let r_i_t: [f64; 3] = [r_i[0] + shift[0], r_i[1] + shift[1], r_i[2] + shift[2]];
    let r_j_t: [f64; 3] = [r_j[0] + shift[0], r_j[1] + shift[1], r_j[2] + shift[2]];
    let r_k_t: [f64; 3] = [r_k[0] + shift[0], r_k[1] + shift[1], r_k[2] + shift[2]];
    let r_l_t: [f64; 3] = [r_l[0] + shift[0], r_l[1] + shift[1], r_l[2] + shift[2]];

    let e_translated = periodic_dihedral_energy(r_i_t, r_j_t, r_k_t, r_l_t, n, phi0, k).unwrap();

    let tol = 1e-12;
    assert!(
        (e_base - e_translated).abs() <= tol,
        "base {e_base}, translated {e_translated}"
    );
}

#[test]
fn test_collinear_error() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];
    let r_k: [f64; 3] = [2.0, 0.0, 0.0];
    let r_l: [f64; 3] = [3.0, 0.0, 0.0];
    let n: i32 = 1;
    let phi0: f64 = 0.0;
    let k: f64 = 1.0;

    let result = periodic_dihedral_energy(r_i, r_j, r_k, r_l, n, phi0, k);
    assert!(
        matches!(result, Err(UcmdError::CoincidentAtoms)),
        "expected CoincidentAtoms, got {result:?}"
    );
}
