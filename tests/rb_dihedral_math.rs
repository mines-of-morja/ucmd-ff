//! Tests of the `rb_dihedral_math` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::rb_dihedral_math::*;

#[test]
fn rb_energy_constant_term() {
    let phi = 0.7_f64;
    let c: [f64; 6] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let got = rb_energy(phi, c);
    let expected = 1.0_f64;
    assert!(
        (got - expected).abs() <= 0.0,
        "got {got}, expected {expected}"
    );
}

#[test]
fn rb_energy_linear_term() {
    let phi = 1.0_f64;
    let c: [f64; 6] = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0];

    let got_e = rb_energy(phi, c);
    let expected_e = -(phi.cos());
    assert!(
        (got_e - expected_e).abs() <= 1e-15,
        "energy: got {got_e}, expected {expected_e}"
    );

    let got_d = rb_denergy_dphi(phi, c);
    let expected_d = phi.sin();
    assert!(
        (got_d - expected_d).abs() <= 1e-15,
        "derivative: got {got_d}, expected {expected_d}"
    );
}

#[test]
fn rb_energy_quadratic_term() {
    let phi = 0.5_f64;
    let c: [f64; 6] = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0];

    let got_e = rb_energy(phi, c);
    let expected_e = phi.cos() * phi.cos();
    assert!(
        (got_e - expected_e).abs() <= 1e-15,
        "energy: got {got_e}, expected {expected_e}"
    );

    let got_d = rb_denergy_dphi(phi, c);
    let expected_d = -2.0 * phi.sin() * phi.cos();
    assert!(
        (got_d - expected_d).abs() <= 1e-15,
        "derivative: got {got_d}, expected {expected_d}"
    );
}

// corrected after review: central-difference truncation error ~2.7e-10 exceeds 1e-10 tolerance
#[test]
fn rb_derivative_finite_difference() {
    let phi = 0.37_f64;
    let c: [f64; 6] = [2.0, -1.5, 3.0, 0.5, -0.25, 0.1];
    let h = 1e-6_f64;

    let fd = (rb_energy(phi + h, c) - rb_energy(phi - h, c)) / (2.0 * h);
    let analytic = rb_denergy_dphi(phi, c);

    assert!(
        (fd - analytic).abs() <= 1e-9,
        "finite difference {fd}, analytic {analytic}"
    );
}

#[test]
fn rb_derivative_zero_at_phi_zero() {
    let phi = 0.0_f64;
    let c: [f64; 6] = [1.0, -2.0, 3.0, -4.0, 5.0, -6.0];
    let got = rb_denergy_dphi(phi, c);
    let expected = 0.0_f64;
    assert!(
        (got - expected).abs() <= 0.0,
        "got {got}, expected {expected}"
    );
}

#[test]
fn rb_energy_at_phi_pi() {
    let phi = std::f64::consts::PI;
    let c: [f64; 6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];

    let got_e = rb_energy(phi, c);
    let expected_e = 1.0 + 2.0 + 3.0 + 4.0 + 5.0 + 6.0;
    assert!(
        (got_e - expected_e).abs() <= 1e-14,
        "energy: got {got_e}, expected {expected_e}"
    );

    let got_d = rb_denergy_dphi(phi, c);
    let expected_d = 0.0_f64;
    assert!(
        (got_d - expected_d).abs() <= 1e-14,
        "derivative: got {got_d}, expected {expected_d}"
    );
}

#[test]
fn rb_agrees_with_opls_scan() {
    // OPLS-Fourier coefficients (arbitrary test parameters)
    let v1 = 1.0_f64;
    let v2 = 2.0_f64;
    let v3 = 3.0_f64;
    let v4 = 4.0_f64;

    // Conversion identities: C_n in terms of V_n
    let c0 = 0.5 * (v1 + 2.0 * v2 + v3);
    let c1 = 0.5 * (3.0 * v3 - v1);
    let c2 = 4.0 * v4 - v2;
    let c3 = -2.0 * v3;
    let c4 = -4.0 * v4;
    let c5 = 0.0_f64;
    let c: [f64; 6] = [c0, c1, c2, c3, c4, c5];

    // OPLS energy: 0.5 * [V1*(1+cos phi) + V2*(1-cos 2phi) + V3*(1+cos 3phi) + V4*(1-cos 4phi)]
    fn opls_energy(phi: f64, v1: f64, v2: f64, v3: f64, v4: f64) -> f64 {
        0.5 * (v1 * (1.0 + phi.cos())
            + v2 * (1.0 - (2.0 * phi).cos())
            + v3 * (1.0 + (3.0 * phi).cos())
            + v4 * (1.0 - (4.0 * phi).cos()))
    }

    let n_points = 200;
    let pi = std::f64::consts::PI;
    let step = 2.0 * pi / n_points as f64;

    for k in 0..n_points {
        let phi = -pi + (k as f64 + 0.5) * step;
        let got = rb_energy(phi, c);
        let expected = opls_energy(phi, v1, v2, v3, v4);
        assert!(
            (got - expected).abs() <= 1e-10,
            "phi index {k}: got {got}, expected {expected}"
        );
    }
}
