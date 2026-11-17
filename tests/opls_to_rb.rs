//! Tests of the `opls_to_rb` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::opls_to_rb::*;
use ucmd_ff::rb_dihedral_math::*;

#[test]
fn test_v1_only() {
    // revised: C1 and C3 formulas corrected to match the revised specification
    let v: [f64; 4] = [2.0, 0.0, 0.0, 0.0];
    let c = opls_to_rb(v);

    let expected: [f64; 6] = [
        0.5 * (v[0] + 2.0 * v[1] + v[2]),
        0.5 * (3.0 * v[2] - v[0]),
        4.0 * v[3] - v[1],
        -2.0 * v[2],
        -4.0 * v[3],
        0.0,
    ];

    for i in 0..6 {
        assert!(
            (c[i] - expected[i]).abs() <= 0.0,
            "C{}: got {c:?}, expected {expected:?}",
            i
        );
    }
}

#[test]
fn test_v2_only() {
    let v: [f64; 4] = [0.0, 2.0, 0.0, 0.0];
    let c = opls_to_rb(v);

    let expected: [f64; 6] = [
        0.5 * (v[0] + 2.0 * v[1] + v[2]),
        0.5 * (3.0 * v[2] - v[0]),
        4.0 * v[3] - v[1],
        -2.0 * v[2],
        -4.0 * v[3],
        0.0,
    ];

    for i in 0..6 {
        assert!(
            (c[i] - expected[i]).abs() <= 0.0,
            "C{}: got {c:?}, expected {expected:?}",
            i
        );
    }
}

#[test]
fn test_v3_only() {
    // revised: C1 and C3 formulas corrected to match the revised specification
    let v: [f64; 4] = [0.0, 0.0, 2.0, 0.0];
    let c = opls_to_rb(v);

    let expected: [f64; 6] = [
        0.5 * (v[0] + 2.0 * v[1] + v[2]),
        0.5 * (3.0 * v[2] - v[0]),
        4.0 * v[3] - v[1],
        -2.0 * v[2],
        -4.0 * v[3],
        0.0,
    ];

    for i in 0..6 {
        assert!(
            (c[i] - expected[i]).abs() <= 0.0,
            "C{}: got {c:?}, expected {expected:?}",
            i
        );
    }
}

#[test]
fn test_v4_only() {
    let v: [f64; 4] = [0.0, 0.0, 0.0, 2.0];
    let c = opls_to_rb(v);

    let expected: [f64; 6] = [
        0.5 * (v[0] + 2.0 * v[1] + v[2]),
        0.5 * (3.0 * v[2] - v[0]),
        4.0 * v[3] - v[1],
        -2.0 * v[2],
        -4.0 * v[3],
        0.0,
    ];

    for i in 0..6 {
        assert!(
            (c[i] - expected[i]).abs() <= 0.0,
            "C{}: got {c:?}, expected {expected:?}",
            i
        );
    }
}

#[test]
fn test_polynomial_identity_scan() {
    // revised: rb_energy is called with phi (not phi - PI) because it internally
    // computes (-cos phi)^n, matching the project definition sum C_n cos^n(psi)
    // with psi = phi - pi
    let v: [f64; 4] = [1.7, -3.2, 0.5, 2.1];
    let c = opls_to_rb(v);

    let n_points: u32 = 200;
    let two_pi: f64 = 2.0 * std::f64::consts::PI;

    for i in 0..n_points {
        let phi: f64 = (i as f64) * two_pi / (n_points as f64);

        let e_opls: f64 = 0.5
            * (v[0] * (1.0 + phi.cos())
                + v[1] * (1.0 - (2.0 * phi).cos())
                + v[2] * (1.0 + (3.0 * phi).cos())
                + v[3] * (1.0 - (4.0 * phi).cos()));

        let e_rb: f64 = rb_energy(phi, c);

        assert!(
            (e_opls - e_rb).abs() <= 1e-12,
            "phi = {phi}: E_OPLS = {e_opls}, E_RB = {e_rb}, diff = {}",
            (e_opls - e_rb).abs()
        );
    }
}

#[test]
fn test_zero_coefficients() {
    let v: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
    let c = opls_to_rb(v);

    let expected: [f64; 6] = [
        0.5 * (v[0] + 2.0 * v[1] + v[2]),
        0.5 * (3.0 * v[2] - v[0]),
        4.0 * v[3] - v[1],
        -2.0 * v[2],
        -4.0 * v[3],
        0.0,
    ];

    for i in 0..6 {
        assert!(
            (c[i] - expected[i]).abs() <= 0.0,
            "C{}: got {c:?}, expected {expected:?}",
            i
        );
    }
}

#[test]
fn test_c5_always_zero() {
    let v: [f64; 4] = [5.0, -1.0, 3.0, 0.7];
    let c = opls_to_rb(v);

    assert!(c[5] == 0.0, "C5: got {c:?}, expected 0.0");
}

#[test]
fn test_mixed_coefficients_known_values() {
    let v: [f64; 4] = [1.0, 2.0, 3.0, 4.0];
    let c = opls_to_rb(v);

    let expected: [f64; 6] = [
        0.5 * (v[0] + 2.0 * v[1] + v[2]),
        0.5 * (3.0 * v[2] - v[0]),
        4.0 * v[3] - v[1],
        -2.0 * v[2],
        -4.0 * v[3],
        0.0,
    ];

    for i in 0..6 {
        assert!(
            (c[i] - expected[i]).abs() <= 0.0,
            "C{}: got {c:?}, expected {expected:?}",
            i
        );
    }
}
