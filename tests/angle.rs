//! Tests of the `angle` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code
// Index loops over atoms and axes mirror the formulas they test.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

use ucmd_ff::angle::*;
use ucmd_ff::errors::UcmdError;

#[test]
fn test_angle_energy_equilibrium() {
    let r1: [f64; 3] = [0.0, 0.0, 0.0];
    let r2: [f64; 3] = [1.0, 0.0, 0.0];
    let r3: [f64; 3] = [1.0, 1.0, 0.0];
    let theta0: f64 = std::f64::consts::PI / 2.0;
    let k: f64 = 1.0;

    let got = angle_energy(r1, r2, r3, theta0, k).unwrap();
    let expected: f64 = 0.0;
    let tol: f64 = 1e-12;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

#[test]
fn test_angle_energy_known_value() {
    let r1: [f64; 3] = [0.0, 0.0, 0.0];
    let r2: [f64; 3] = [1.0, 0.0, 0.0];
    let r3: [f64; 3] = [1.0, 1.0, 0.0];
    let theta0: f64 = std::f64::consts::PI / 4.0;
    let k: f64 = 2.0;

    let got = angle_energy(r1, r2, r3, theta0, k).unwrap();
    let delta_theta: f64 = std::f64::consts::PI / 2.0 - theta0;
    let expected: f64 = 0.5 * k * delta_theta * delta_theta;
    let tol: f64 = 1e-10;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

#[test]
fn test_angle_forces_finite_difference() {
    let r1: [f64; 3] = [0.1, 0.2, 0.3];
    let r2: [f64; 3] = [1.0, 0.5, 0.1];
    let r3: [f64; 3] = [0.8, 1.2, 0.4];
    let theta0: f64 = 1.0;
    let k: f64 = 100.0;
    let h: f64 = 1e-5;

    let (f1, f2, f3) = angle_forces(r1, r2, r3, theta0, k).unwrap();

    let mut fd1: [f64; 3] = [0.0; 3];
    let mut fd2: [f64; 3] = [0.0; 3];
    let mut fd3: [f64; 3] = [0.0; 3];

    for dim in 0..3 {
        // Perturb r1
        let mut r1p = r1;
        let mut r1m = r1;
        r1p[dim] += h;
        r1m[dim] -= h;
        let ep = angle_energy(r1p, r2, r3, theta0, k).unwrap();
        let em = angle_energy(r1m, r2, r3, theta0, k).unwrap();
        fd1[dim] = (ep - em) / (2.0 * h);

        // Perturb r2
        let mut r2p = r2;
        let mut r2m = r2;
        r2p[dim] += h;
        r2m[dim] -= h;
        let ep = angle_energy(r1, r2p, r3, theta0, k).unwrap();
        let em = angle_energy(r1, r2m, r3, theta0, k).unwrap();
        fd2[dim] = (ep - em) / (2.0 * h);

        // Perturb r3
        let mut r3p = r3;
        let mut r3m = r3;
        r3p[dim] += h;
        r3m[dim] -= h;
        let ep = angle_energy(r1, r2, r3p, theta0, k).unwrap();
        let em = angle_energy(r1, r2, r3m, theta0, k).unwrap();
        fd3[dim] = (ep - em) / (2.0 * h);
    }

    // F = -dE/dr, so F_fd = -fd
    let tol: f64 = 1e-5;
    for dim in 0..3 {
        let f1_fd = -fd1[dim];
        let f2_fd = -fd2[dim];
        let f3_fd = -fd3[dim];

        let scale1 = f1_fd.abs().max(1.0);
        let scale2 = f2_fd.abs().max(1.0);
        let scale3 = f3_fd.abs().max(1.0);

        assert!(
            (f1[dim] - f1_fd).abs() <= tol * scale1,
            "F1[{}]: got {}, expected {}",
            dim,
            f1[dim],
            f1_fd
        );
        assert!(
            (f2[dim] - f2_fd).abs() <= tol * scale2,
            "F2[{}]: got {}, expected {}",
            dim,
            f2[dim],
            f2_fd
        );
        assert!(
            (f3[dim] - f3_fd).abs() <= tol * scale3,
            "F3[{}]: got {}, expected {}",
            dim,
            f3[dim],
            f3_fd
        );
    }
}

#[test]
fn test_angle_forces_momentum_conservation() {
    let r1: [f64; 3] = [0.1, 0.2, 0.3];
    let r2: [f64; 3] = [1.0, 0.5, 0.1];
    let r3: [f64; 3] = [0.8, 1.2, 0.4];
    let theta0: f64 = 1.0;
    let k: f64 = 100.0;

    let (f1, f2, f3) = angle_forces(r1, r2, r3, theta0, k).unwrap();

    let sum: [f64; 3] = [
        f1[0] + f2[0] + f3[0],
        f1[1] + f2[1] + f3[1],
        f1[2] + f2[2] + f3[2],
    ];
    let tol: f64 = 1e-10;
    for dim in 0..3 {
        assert!(
            sum[dim].abs() <= tol,
            "Sum of forces component {}: got {}",
            dim,
            sum[dim]
        );
    }
}

#[test]
fn test_angle_virial_finite_difference() {
    let r1: [f64; 3] = [0.1, 0.2, 0.3];
    let r2: [f64; 3] = [1.0, 0.5, 0.1];
    let r3: [f64; 3] = [0.8, 1.2, 0.4];
    let theta0: f64 = 1.0;
    let k: f64 = 100.0;
    let eps: f64 = 1e-5;

    let w = angle_virial(r1, r2, r3, theta0, k).unwrap();

    // Virial W_ab = -dE/d(eps_ab)
    // Deformation: r'_a = r_a + eps_ab * r_b
    // dE/d(eps_ab) = sum_i (dE/dr_i)_a * r_i_b = - sum_i F_i_a * r_i_b
    // W_ab = sum_i F_i_a * r_i_b
    // Since F1+F2+F3=0, W_ab = F1_a*r1_b + F2_a*r2_b + F3_a*r3_b
    // Or using pair-like form: W_ab = r12_a * F1_b + r32_a * F3_b (since F2 = -F1-F3)
    // Let's verify with finite difference.

    let mut w_fd: [[f64; 3]; 3] = [[0.0; 3]; 3];

    for a in 0..3 {
        for b in 0..3 {
            // Perturb coordinate a of all atoms by eps * coordinate b
            let mut r1p = r1;
            let mut r2p = r2;
            let mut r3p = r3;
            let mut r1m = r1;
            let mut r2m = r2;
            let mut r3m = r3;

            r1p[a] += eps * r1[b];
            r2p[a] += eps * r2[b];
            r3p[a] += eps * r3[b];

            r1m[a] -= eps * r1[b];
            r2m[a] -= eps * r2[b];
            r3m[a] -= eps * r3[b];

            let ep = angle_energy(r1p, r2p, r3p, theta0, k).unwrap();
            let em = angle_energy(r1m, r2m, r3m, theta0, k).unwrap();

            // W_ab = -dE/d(eps_ab) = -(ep - em) / (2 * eps)
            w_fd[a][b] = -(ep - em) / (2.0 * eps);
        }
    }

    let tol: f64 = 1e-5;
    for a in 0..3 {
        for b in 0..3 {
            let scale = w_fd[a][b].abs().max(1.0);
            assert!(
                (w[a][b] - w_fd[a][b]).abs() <= tol * scale,
                "W[{}][{}]: got {}, expected {}",
                a,
                b,
                w[a][b],
                w_fd[a][b]
            );
        }
    }
}

#[test]
fn test_angle_error_coincident() {
    let r1: [f64; 3] = [1.0, 0.0, 0.0];
    let r2: [f64; 3] = [1.0, 0.0, 0.0];
    let r3: [f64; 3] = [2.0, 0.0, 0.0];
    let theta0: f64 = 1.0;
    let k: f64 = 1.0;

    let result = angle_energy(r1, r2, r3, theta0, k);
    assert!(
        matches!(result, Err(UcmdError::CoincidentAtoms)),
        "Expected CoincidentAtoms, got {:?}",
        result
    );
}

#[test]
fn test_angle_error_straight() {
    let r1: [f64; 3] = [0.0, 0.0, 0.0];
    let r2: [f64; 3] = [1.0, 0.0, 0.0];
    let r3: [f64; 3] = [2.0, 0.0, 0.0];
    let theta0: f64 = 1.0;
    let k: f64 = 1.0;

    let result = angle_energy(r1, r2, r3, theta0, k);
    assert!(
        matches!(result, Err(UcmdError::CoincidentAtoms)),
        "Expected CoincidentAtoms, got {:?}",
        result
    );
}

#[test]
fn test_angle_energy_non_unit_lengths() {
    let r1: [f64; 3] = [0.0, 0.0, 0.0];
    let r2: [f64; 3] = [2.0, 0.0, 0.0];
    let r3: [f64; 3] = [2.0, 3.0, 0.0];
    let theta0: f64 = std::f64::consts::PI / 4.0;
    let k: f64 = 2.0;

    let got = angle_energy(r1, r2, r3, theta0, k).unwrap();
    let delta_theta: f64 = std::f64::consts::PI / 2.0 - theta0;
    let expected: f64 = 0.5 * k * delta_theta * delta_theta;
    let tol: f64 = 1e-10;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}
