//! Tests of the `lj_pair` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::constants::*;
use ucmd_ff::errors::*;
use ucmd_ff::lj_pair::*;

#[test]
fn test_lj_energy_at_minimum() {
    let sigma_nm: f64 = 0.34;
    let eps_kjmol: f64 = 120.0 * K_B_J_K * AVOGADRO_MOL / 1000.0;
    let r_min: f64 = 2.0_f64.powf(1.0 / 6.0) * sigma_nm;
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [r_min, 0.0, 0.0];

    let got = lj_pair_energy(sigma_nm, eps_kjmol, r_i, r_j).unwrap();
    let expected: f64 = -eps_kjmol;
    assert!(
        (got - expected).abs() <= 1e-12,
        "got {got}, expected {expected}"
    );
}

// corrected after review: tolerance 1e-15 is below the unavoidable f64 rounding residual (~7e-15) in the force bracket at the LJ minimum
#[test]
fn test_lj_force_at_minimum() {
    let sigma_nm: f64 = 0.34;
    let eps_kjmol: f64 = 120.0 * K_B_J_K * AVOGADRO_MOL / 1000.0;
    let r_min: f64 = 2.0_f64.powf(1.0 / 6.0) * sigma_nm;
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [r_min, 0.0, 0.0];

    let (f_i, f_j) = lj_pair_forces(sigma_nm, eps_kjmol, r_i, r_j).unwrap();
    for a in 0..3 {
        assert!(f_i[a].abs() <= 1e-14, "F_i[{a}] = {f_i:?}, expected 0");
        assert!(f_j[a].abs() <= 1e-14, "F_j[{a}] = {f_j:?}, expected 0");
    }
}

#[test]
fn test_lj_force_finite_difference() {
    let sigma_nm: f64 = 0.34;
    let eps_kjmol: f64 = 120.0 * K_B_J_K * AVOGADRO_MOL / 1000.0;
    let r_i: [f64; 3] = [0.5, 0.3, 0.1];
    let r_j: [f64; 3] = [0.0, 0.0, 0.0];
    let h: f64 = 1e-8;

    let (f_i, _) = lj_pair_forces(sigma_nm, eps_kjmol, r_i, r_j).unwrap();

    for a in 0..3 {
        let mut r_plus = r_i;
        let mut r_minus = r_i;
        r_plus[a] += h;
        r_minus[a] -= h;

        let e_plus = lj_pair_energy(sigma_nm, eps_kjmol, r_plus, r_j).unwrap();
        let e_minus = lj_pair_energy(sigma_nm, eps_kjmol, r_minus, r_j).unwrap();
        let fd: f64 = -(e_plus - e_minus) / (2.0 * h);

        assert!(
            (f_i[a] - fd).abs() <= 1e-6,
            "axis {a}: analytic F_i[{a}] = {f_i:?}, finite-diff = {fd}"
        );
    }
}

#[test]
fn test_lj_virial_finite_difference() {
    let sigma_nm: f64 = 0.34;
    let eps_kjmol: f64 = 120.0 * K_B_J_K * AVOGADRO_MOL / 1000.0;
    let r_i: [f64; 3] = [0.5, 0.0, 0.0];
    let r_j: [f64; 3] = [0.0, 0.0, 0.0];
    let eps: f64 = 1e-8;

    let w = lj_pair_virial(sigma_nm, eps_kjmol, r_i, r_j).unwrap();

    // Strained positions: r_i'[0] = r_i[0] * (1 + eps), r_j unchanged (r_j[0] = 0)
    let r_i_plus: [f64; 3] = [r_i[0] * (1.0 + eps), r_i[1], r_i[2]];
    let r_i_minus: [f64; 3] = [r_i[0] * (1.0 - eps), r_i[1], r_i[2]];

    let e_plus = lj_pair_energy(sigma_nm, eps_kjmol, r_i_plus, r_j).unwrap();
    let e_minus = lj_pair_energy(sigma_nm, eps_kjmol, r_i_minus, r_j).unwrap();
    let fd: f64 = -(e_plus - e_minus) / (2.0 * eps);

    assert!(
        (w[0][0] - fd).abs() <= 1e-6,
        "W[0][0] = {w:?}, finite-diff = {fd}"
    );
}

#[test]
fn test_lj_newtons_third_law() {
    let sigma_nm: f64 = 0.34;
    let eps_kjmol: f64 = 0.5;
    let r_i: [f64; 3] = [1.2, 0.7, 0.3];
    let r_j: [f64; 3] = [0.4, 1.1, 0.9];

    let (f_i, f_j) = lj_pair_forces(sigma_nm, eps_kjmol, r_i, r_j).unwrap();

    for a in 0..3 {
        let sum: f64 = f_i[a] + f_j[a];
        assert!(
            sum.abs() <= 1e-15,
            "component {a}: F_i[{a}] + F_j[{a}] = {sum}, expected 0"
        );
    }
}

#[test]
fn test_lj_zero_sigma() {
    let sigma_nm: f64 = 0.0;
    let eps_kjmol: f64 = 0.6502994552;
    let r_i: [f64; 3] = [0.1, 0.2, 0.3];
    let r_j: [f64; 3] = [0.4, 0.5, 0.6];

    let e = lj_pair_energy(sigma_nm, eps_kjmol, r_i, r_j).unwrap();
    assert_eq!(e, 0.0, "energy with sigma=0 should be exactly 0.0, got {e}");

    let (f_i, f_j) = lj_pair_forces(sigma_nm, eps_kjmol, r_i, r_j).unwrap();
    for a in 0..3 {
        assert_eq!(
            f_i[a], 0.0,
            "F_i[{a}] with sigma=0 should be exactly 0.0, got {f_i:?}"
        );
        assert_eq!(
            f_j[a], 0.0,
            "F_j[{a}] with sigma=0 should be exactly 0.0, got {f_j:?}"
        );
    }

    let w = lj_pair_virial(sigma_nm, eps_kjmol, r_i, r_j).unwrap();
    for a in 0..3 {
        for b in 0..3 {
            assert_eq!(
                w[a][b], 0.0,
                "W[{a}][{b}] with sigma=0 should be exactly 0.0, got {w:?}"
            );
        }
    }
}

#[test]
fn test_lj_coincident_atoms() {
    let sigma_nm: f64 = 0.34;
    let eps_kjmol: f64 = 0.5;
    let r_i: [f64; 3] = [1.0, 1.0, 1.0];
    let r_j: [f64; 3] = [1.0, 1.0, 1.0];

    let e_result = lj_pair_energy(sigma_nm, eps_kjmol, r_i, r_j);
    assert!(
        matches!(e_result, Err(UcmdError::CoincidentAtoms)),
        "lj_pair_energy should return Err(CoincidentAtoms), got {e_result:?}"
    );

    let f_result = lj_pair_forces(sigma_nm, eps_kjmol, r_i, r_j);
    assert!(
        matches!(f_result, Err(UcmdError::CoincidentAtoms)),
        "lj_pair_forces should return Err(CoincidentAtoms), got {f_result:?}"
    );

    let w_result = lj_pair_virial(sigma_nm, eps_kjmol, r_i, r_j);
    assert!(
        matches!(w_result, Err(UcmdError::CoincidentAtoms)),
        "lj_pair_virial should return Err(CoincidentAtoms), got {w_result:?}"
    );
}
