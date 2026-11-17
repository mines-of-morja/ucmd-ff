//! Tests of the `dihedral_angle` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::dihedral_angle::*;
use ucmd_ff::errors::UcmdError;

#[test]
fn test_trans_angle() {
    let r_i: [f64; 3] = [0.0, 1.0, 0.0];
    let r_j: [f64; 3] = [0.0, 0.0, 0.0];
    let r_k: [f64; 3] = [1.0, 0.0, 0.0];
    let r_l: [f64; 3] = [1.0, -1.0, 0.0];
    let phi = dihedral_angle(r_i, r_j, r_k, r_l).unwrap();
    let expected = std::f64::consts::PI;
    assert!(
        (phi - expected).abs() <= 1e-12,
        "got {phi}, expected {expected}"
    );
}

#[test]
fn test_cis_angle() {
    let r_i: [f64; 3] = [0.0, 1.0, 0.0];
    let r_j: [f64; 3] = [0.0, 0.0, 0.0];
    let r_k: [f64; 3] = [1.0, 0.0, 0.0];
    let r_l: [f64; 3] = [1.0, 1.0, 0.0];
    let phi = dihedral_angle(r_i, r_j, r_k, r_l).unwrap();
    let expected = 0.0;
    assert!(
        (phi - expected).abs() <= 1e-12,
        "got {phi}, expected {expected}"
    );
}

#[test]
fn test_gradient_finite_difference() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];
    let r_k: [f64; 3] = [1.0, 1.0, 0.0];
    let r_l: [f64; 3] = [2.0, 0.0, 1.0];
    let h = 1e-6;
    let tol = 1e-8;
    let (g_i, g_j, g_k, g_l) = dihedral_gradient(r_i, r_j, r_k, r_l).unwrap();
    let atoms = [r_i, r_j, r_k, r_l];
    let grads = [g_i, g_j, g_k, g_l];
    for (m, &r_m) in atoms.iter().enumerate() {
        for c in 0..3 {
            let mut r_plus = r_m;
            let mut r_minus = r_m;
            r_plus[c] += h;
            r_minus[c] -= h;
            let a = atoms[0];
            let b = atoms[1];
            let c_ = atoms[2];
            let d = atoms[3];
            let phi_plus = dihedral_angle(
                if m == 0 { r_plus } else { a },
                if m == 1 { r_plus } else { b },
                if m == 2 { r_plus } else { c_ },
                if m == 3 { r_plus } else { d },
            )
            .unwrap();
            let phi_minus = dihedral_angle(
                if m == 0 { r_minus } else { a },
                if m == 1 { r_minus } else { b },
                if m == 2 { r_minus } else { c_ },
                if m == 3 { r_minus } else { d },
            )
            .unwrap();
            let fd = (phi_plus - phi_minus) / (2.0 * h);
            let got = grads[m][c];
            assert!(
                (got - fd).abs() <= tol,
                "atom {m} coord {c}: got {got}, expected {fd}"
            );
        }
    }
}

#[test]
fn test_collinear_error() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];
    let r_k: [f64; 3] = [2.0, 0.0, 0.0];
    let r_l: [f64; 3] = [3.0, 0.0, 0.0];
    let result = dihedral_angle(r_i, r_j, r_k, r_l);
    assert!(
        matches!(result, Err(UcmdError::CoincidentAtoms)),
        "expected Err(CoincidentAtoms), got {result:?}"
    );
}

#[test]
fn test_translation_invariance() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];
    let r_k: [f64; 3] = [1.0, 1.0, 0.0];
    let r_l: [f64; 3] = [2.0, 0.0, 1.0];
    let shift: [f64; 3] = [10.0, 20.0, 30.0];
    let r_i_t = [r_i[0] + shift[0], r_i[1] + shift[1], r_i[2] + shift[2]];
    let r_j_t = [r_j[0] + shift[0], r_j[1] + shift[1], r_j[2] + shift[2]];
    let r_k_t = [r_k[0] + shift[0], r_k[1] + shift[1], r_k[2] + shift[2]];
    let r_l_t = [r_l[0] + shift[0], r_l[1] + shift[1], r_l[2] + shift[2]];
    let phi = dihedral_angle(r_i, r_j, r_k, r_l).unwrap();
    let phi_t = dihedral_angle(r_i_t, r_j_t, r_k_t, r_l_t).unwrap();
    assert!((phi - phi_t).abs() <= 1e-12, "got {phi}, expected {phi_t}");
    let (g_i, g_j, g_k, g_l) = dihedral_gradient(r_i, r_j, r_k, r_l).unwrap();
    let (g_i_t, g_j_t, g_k_t, g_l_t) = dihedral_gradient(r_i_t, r_j_t, r_k_t, r_l_t).unwrap();
    let grads = [g_i, g_j, g_k, g_l];
    let grads_t = [g_i_t, g_j_t, g_k_t, g_l_t];
    for m in 0..4 {
        for c in 0..3 {
            assert!(
                (grads[m][c] - grads_t[m][c]).abs() <= 1e-12,
                "atom {m} coord {c}: got {}, expected {}",
                grads[m][c],
                grads_t[m][c]
            );
        }
    }
}

// revised: added test_quadrant_gradient_finite_difference per revised specification (replaces the old test_trans_gradient_finite_difference which hit the atan2 branch discontinuity at phi = pi)
#[test]
fn test_quadrant_gradient_finite_difference() {
    let r_i: [f64; 3] = [0.0, 1.0, 0.0];
    let r_j: [f64; 3] = [0.0, 0.0, 0.0];
    let r_k: [f64; 3] = [1.0, 0.0, 0.0];
    let r_l: [f64; 3] = [1.0, 0.0, 1.0];
    let h = 1e-6;
    let tol = 1e-8;
    let (g_i, g_j, g_k, g_l) = dihedral_gradient(r_i, r_j, r_k, r_l).unwrap();
    let atoms = [r_i, r_j, r_k, r_l];
    let grads = [g_i, g_j, g_k, g_l];
    for (m, &r_m) in atoms.iter().enumerate() {
        for c in 0..3 {
            let mut r_plus = r_m;
            let mut r_minus = r_m;
            r_plus[c] += h;
            r_minus[c] -= h;
            let a = atoms[0];
            let b = atoms[1];
            let c_ = atoms[2];
            let d = atoms[3];
            let phi_plus = dihedral_angle(
                if m == 0 { r_plus } else { a },
                if m == 1 { r_plus } else { b },
                if m == 2 { r_plus } else { c_ },
                if m == 3 { r_plus } else { d },
            )
            .unwrap();
            let phi_minus = dihedral_angle(
                if m == 0 { r_minus } else { a },
                if m == 1 { r_minus } else { b },
                if m == 2 { r_minus } else { c_ },
                if m == 3 { r_minus } else { d },
            )
            .unwrap();
            let fd = (phi_plus - phi_minus) / (2.0 * h);
            let got = grads[m][c];
            assert!(
                (got - fd).abs() <= tol,
                "atom {m} coord {c}: got {got}, expected {fd}"
            );
        }
    }
}
