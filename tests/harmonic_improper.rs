//! Tests of the `harmonic_improper` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::dihedral_angle::*;
use ucmd_ff::harmonic_improper::*;
use ucmd_ff::vec3::*;

#[test]
fn test_energy_zero_at_equilibrium() {
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.5, 0.5, 0.0];
    let r_l = [1.5, 1.0, 0.0];
    let xi_val = dihedral_angle(r_i, r_j, r_k, r_l).unwrap();
    let k: f64 = 50.0;
    let e = harmonic_improper_energy(r_i, r_j, r_k, r_l, xi_val, k).unwrap();
    assert!(e == 0.0, "got {e}, expected 0.0");
}

#[test]
fn test_energy_known_deviation() {
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.5, 0.5, 0.0];
    let r_l = [1.5, 1.0, 0.0];
    let xi_val = dihedral_angle(r_i, r_j, r_k, r_l).unwrap();
    let xi0 = xi_val - 0.5;
    let k: f64 = 40.0;
    let dr: f64 = 0.5;
    let expected = 0.5 * k * dr * dr;
    let e = harmonic_improper_energy(r_i, r_j, r_k, r_l, xi0, k).unwrap();
    assert!(
        (e - expected).abs() <= 1e-12,
        "got {e}, expected {expected}"
    );
}

#[test]
fn test_forces_finite_difference() {
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.5, 0.5, 0.0];
    let r_l = [1.5, 1.0, 0.0];
    let xi0: f64 = 0.0;
    let k: f64 = 30.0;
    let h: f64 = 1e-6;

    let (fi, fj, fk, fl) = harmonic_improper_forces(r_i, r_j, r_k, r_l, xi0, k).unwrap();

    let positions: [[f64; 3]; 4] = [r_i, r_j, r_k, r_l];
    let forces: [[f64; 3]; 4] = [fi, fj, fk, fl];

    for p in 0..4 {
        for c in 0..3 {
            let mut pos_plus = positions;
            let mut pos_minus = positions;
            pos_plus[p][c] += h;
            pos_minus[p][c] -= h;

            let e_plus = harmonic_improper_energy(
                pos_plus[0],
                pos_plus[1],
                pos_plus[2],
                pos_plus[3],
                xi0,
                k,
            )
            .unwrap();
            let e_minus = harmonic_improper_energy(
                pos_minus[0],
                pos_minus[1],
                pos_minus[2],
                pos_minus[3],
                xi0,
                k,
            )
            .unwrap();

            let fd = -(e_plus - e_minus) / (2.0 * h);
            let analytic = forces[p][c];

            let tol = 1e-5 * fd.abs().max(analytic.abs()).max(1e-30);
            assert!(
                (analytic - fd).abs() <= tol,
                "atom {p} component {c}: analytic {analytic}, fd {fd}"
            );
        }
    }
}

#[test]
fn test_virial_finite_difference() {
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.5, 0.5, 0.0];
    let r_l = [1.5, 1.0, 0.0];
    let xi0: f64 = 0.0;
    let k: f64 = 30.0;
    let eps: f64 = 1e-6;

    let w = harmonic_improper_virial(r_i, r_j, r_k, r_l, xi0, k).unwrap();

    let positions: [[f64; 3]; 4] = [r_i, r_j, r_k, r_l];

    let pairs: [(usize, usize); 6] = [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)];

    for (a, b) in pairs {
        let mut pos_plus = positions;
        let mut pos_minus = positions;
        for p in 0..4 {
            pos_plus[p][b] += eps * positions[p][a];
            pos_minus[p][b] -= eps * positions[p][a];
        }

        let e_plus =
            harmonic_improper_energy(pos_plus[0], pos_plus[1], pos_plus[2], pos_plus[3], xi0, k)
                .unwrap();
        let e_minus = harmonic_improper_energy(
            pos_minus[0],
            pos_minus[1],
            pos_minus[2],
            pos_minus[3],
            xi0,
            k,
        )
        .unwrap();

        let fd = -(e_plus - e_minus) / (2.0 * eps);
        let analytic = w[a][b];

        let tol = 1e-5 * fd.abs().max(analytic.abs()).max(1e-30);
        assert!(
            (analytic - fd).abs() <= tol,
            "W[{a}][{b}]: analytic {analytic}, fd {fd}"
        );
    }
}

#[test]
fn test_wrapping_at_boundary() {
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.5, 0.5, 0.0];
    let r_l = [1.5, 1.0, 0.0];
    let xi_val = dihedral_angle(r_i, r_j, r_k, r_l).unwrap();
    let k: f64 = 10.0;
    let pi: f64 = std::f64::consts::PI;

    // Case 1: xi0 = xi_val - pi, so raw delta = pi
    let xi0_1 = xi_val - pi;
    let e1 = harmonic_improper_energy(r_i, r_j, r_k, r_l, xi0_1, k).unwrap();
    let expected = 0.5 * k * pi * pi;
    assert!(
        (e1 - expected).abs() <= 1e-10,
        "case 1: got {e1}, expected {expected}"
    );

    // Case 2: xi0 = xi_val + pi, so raw delta = -pi, wraps to pi
    let xi0_2 = xi_val + pi;
    let e2 = harmonic_improper_energy(r_i, r_j, r_k, r_l, xi0_2, k).unwrap();
    assert!(
        (e2 - expected).abs() <= 1e-10,
        "case 2: got {e2}, expected {expected}"
    );
}

#[test]
fn test_momentum_conservation() {
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.5, 0.5, 0.0];
    let r_l = [2.0, 0.3, 0.1];
    let xi0: f64 = 0.5;
    let k: f64 = 25.0;

    let (fi, fj, fk, fl) = harmonic_improper_forces(r_i, r_j, r_k, r_l, xi0, k).unwrap();

    let sum = [
        fi[0] + fj[0] + fk[0] + fl[0],
        fi[1] + fj[1] + fk[1] + fl[1],
        fi[2] + fj[2] + fk[2] + fl[2],
    ];

    let total_norm = norm(fi) + norm(fj) + norm(fk) + norm(fl);
    let sum_norm = norm(sum);

    assert!(
        sum_norm <= 1e-8 * total_norm,
        "momentum not conserved: |sum| = {sum_norm}, total = {total_norm}"
    );
}

#[test]
fn test_translation_invariance() {
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.5, 0.5, 0.0];
    let r_l = [1.5, 1.0, 0.0];
    let xi0: f64 = 0.3;
    let k: f64 = 20.0;
    let t = [3.7, -1.2, 0.8];

    let e_orig = harmonic_improper_energy(r_i, r_j, r_k, r_l, xi0, k).unwrap();

    let ti = add(r_i, t);
    let tj = add(r_j, t);
    let tk = add(r_k, t);
    let tl = add(r_l, t);

    let e_trans = harmonic_improper_energy(ti, tj, tk, tl, xi0, k).unwrap();

    assert!(
        (e_trans - e_orig).abs() <= 1e-12,
        "got {e_trans}, expected {e_orig}"
    );
}

#[test]
fn test_forces_sign_with_nonzero_xi0() {
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.5, 0.5, 0.0];
    let r_l = [1.5, 1.0, 0.0];
    let xi0: f64 = 1.0;
    let k: f64 = 10.0;
    let h: f64 = 1e-6;

    let (fi, fj, fk, fl) = harmonic_improper_forces(r_i, r_j, r_k, r_l, xi0, k).unwrap();

    let positions: [[f64; 3]; 4] = [r_i, r_j, r_k, r_l];
    let forces: [[f64; 3]; 4] = [fi, fj, fk, fl];

    for p in 0..4 {
        for c in 0..3 {
            let mut pos_plus = positions;
            let mut pos_minus = positions;
            pos_plus[p][c] += h;
            pos_minus[p][c] -= h;

            let e_plus = harmonic_improper_energy(
                pos_plus[0],
                pos_plus[1],
                pos_plus[2],
                pos_plus[3],
                xi0,
                k,
            )
            .unwrap();
            let e_minus = harmonic_improper_energy(
                pos_minus[0],
                pos_minus[1],
                pos_minus[2],
                pos_minus[3],
                xi0,
                k,
            )
            .unwrap();

            let fd = -(e_plus - e_minus) / (2.0 * h);
            let analytic = forces[p][c];

            let tol = 1e-5 * fd.abs().max(analytic.abs()).max(1e-30);
            assert!(
                (analytic - fd).abs() <= tol,
                "atom {p} component {c}: analytic {analytic}, fd {fd}"
            );
        }
    }
}
