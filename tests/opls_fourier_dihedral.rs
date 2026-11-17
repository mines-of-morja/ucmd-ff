//! Tests of the `opls_fourier_dihedral` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code
// Index loops over atoms and axes mirror the formulas they test.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

use ucmd_ff::opls_fourier_dihedral::*;

fn energy_value(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
    v1: f64,
    v2: f64,
    v3: f64,
    v4: f64,
) -> f64 {
    match opls_fourier_dihedral_energy(r_i, r_j, r_k, r_l, v1, v2, v3, v4) {
        Ok(e) => e,
        Err(e) => panic!("energy error: {:?}", e),
    }
}

fn perturb_coord(pos: [f64; 3], axis: usize, delta: f64) -> [f64; 3] {
    let mut p = pos;
    p[axis] += delta;
    p
}

// revised: geometry corrected to match specification (r_l = [0,1,0] gives phi=0, not [2,1,0] which gives phi=pi)
#[test]
fn test_energy_at_phi_zero() {
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.0, 1.0, 0.0];
    let r_l = [0.0, 1.0, 0.0];
    let v1 = 1.0;
    let v2 = 0.0;
    let v3 = 0.0;
    let v4 = 0.0;

    let got = energy_value(r_i, r_j, r_k, r_l, v1, v2, v3, v4);
    // phi = 0, E = 0.5 * v1 * (1 + cos(0)) = 0.5 * 1.0 * 2.0 = 1.0
    let expected = 0.5 * v1 * (1.0 + (0.0_f64).cos());
    let tol = 1e-12;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

#[test]
fn test_energy_at_phi_pi() {
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.0, 1.0, 0.0];
    let r_l = [2.0, 1.0, 0.0];
    let v1 = 1.0;
    let v2 = 0.0;
    let v3 = 0.0;
    let v4 = 0.0;

    let got = energy_value(r_i, r_j, r_k, r_l, v1, v2, v3, v4);
    // phi = pi, E = 0.5 * v1 * (1 + cos(pi)) = 0.5 * 1.0 * (1 + (-1)) = 0.0
    let expected = 0.5 * v1 * (1.0 + std::f64::consts::PI.cos());
    let tol = 1e-12;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

#[test]
fn test_forces_finite_difference() {
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.0, 1.0, 0.0];
    let r_l = [2.0, 1.0, 0.0];
    let v1 = 1.0;
    let v2 = 2.0;
    let v3 = 3.0;
    let v4 = 4.0;
    let h = 1e-5;

    let (f_i, f_j, f_k, f_l) =
        opls_fourier_dihedral_forces(r_i, r_j, r_k, r_l, v1, v2, v3, v4).expect("forces failed");

    let atoms = [r_i, r_j, r_k, r_l];
    let forces = [f_i, f_j, f_k, f_l];

    for a in 0..4 {
        for axis in 0..3 {
            let pos_plus = perturb_coord(atoms[a], axis, h);
            let pos_minus = perturb_coord(atoms[a], axis, -h);

            let (ri, rj, rk, rl) = match a {
                0 => (pos_plus, r_j, r_k, r_l),
                1 => (r_i, pos_plus, r_k, r_l),
                2 => (r_i, r_j, pos_plus, r_l),
                _ => (r_i, r_j, r_k, pos_plus),
            };
            let (ri_m, rj_m, rk_m, rl_m) = match a {
                0 => (pos_minus, r_j, r_k, r_l),
                1 => (r_i, pos_minus, r_k, r_l),
                2 => (r_i, r_j, pos_minus, r_l),
                _ => (r_i, r_j, r_k, pos_minus),
            };

            let e_plus = energy_value(ri, rj, rk, rl, v1, v2, v3, v4);
            let e_minus = energy_value(ri_m, rj_m, rk_m, rl_m, v1, v2, v3, v4);
            let expected = -(e_plus - e_minus) / (2.0 * h);

            let got = forces[a][axis];
            let tol = 1e-6;
            assert!(
                (got - expected).abs() <= tol,
                "atom {a} axis {axis}: got {got}, expected {expected}"
            );
        }
    }
}

#[test]
fn test_virial_finite_difference() {
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.0, 1.0, 0.0];
    let r_l = [2.0, 1.0, 0.0];
    let v1 = 1.0;
    let v2 = 2.0;
    let v3 = 3.0;
    let v4 = 4.0;
    let eps = 1e-5;

    let got_virial =
        opls_fourier_dihedral_virial(r_i, r_j, r_k, r_l, v1, v2, v3, v4).expect("virial failed");

    let atoms = [r_i, r_j, r_k, r_l];

    for a in 0..3 {
        for b in 0..3 {
            let mut pos_plus = [[0.0; 3]; 4];
            let mut pos_minus = [[0.0; 3]; 4];
            for idx in 0..4 {
                let mut p_plus = atoms[idx];
                let mut p_minus = atoms[idx];
                p_plus[a] += eps * atoms[idx][b];
                p_minus[a] -= eps * atoms[idx][b];
                pos_plus[idx] = p_plus;
                pos_minus[idx] = p_minus;
            }

            let e_plus = energy_value(
                pos_plus[0],
                pos_plus[1],
                pos_plus[2],
                pos_plus[3],
                v1,
                v2,
                v3,
                v4,
            );
            let e_minus = energy_value(
                pos_minus[0],
                pos_minus[1],
                pos_minus[2],
                pos_minus[3],
                v1,
                v2,
                v3,
                v4,
            );
            let expected = -(e_plus - e_minus) / (2.0 * eps);

            let got = got_virial[a][b];
            let tol = 1e-6;
            assert!(
                (got - expected).abs() <= tol,
                "virial [{a}][{b}]: got {got}, expected {expected}"
            );
        }
    }
}

#[test]
fn test_translation_invariance() {
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.0, 1.0, 0.0];
    let r_l = [2.0, 1.0, 0.0];
    let v1 = 1.0;
    let v2 = 2.0;
    let v3 = 3.0;
    let v4 = 4.0;

    let shift = [10.0, 20.0, 30.0];

    let r_i_t = [r_i[0] + shift[0], r_i[1] + shift[1], r_i[2] + shift[2]];
    let r_j_t = [r_j[0] + shift[0], r_j[1] + shift[1], r_j[2] + shift[2]];
    let r_k_t = [r_k[0] + shift[0], r_k[1] + shift[1], r_k[2] + shift[2]];
    let r_l_t = [r_l[0] + shift[0], r_l[1] + shift[1], r_l[2] + shift[2]];

    let e_orig = energy_value(r_i, r_j, r_k, r_l, v1, v2, v3, v4);
    let e_trans = energy_value(r_i_t, r_j_t, r_k_t, r_l_t, v1, v2, v3, v4);
    let tol = 1e-12;
    assert!(
        (e_orig - e_trans).abs() <= tol,
        "energy: got {e_trans}, expected {e_orig}"
    );

    let (f_i, f_j, f_k, f_l) =
        opls_fourier_dihedral_forces(r_i, r_j, r_k, r_l, v1, v2, v3, v4).expect("forces failed");
    let (f_i_t, f_j_t, f_k_t, f_l_t) =
        opls_fourier_dihedral_forces(r_i_t, r_j_t, r_k_t, r_l_t, v1, v2, v3, v4)
            .expect("forces failed");

    let forces_orig = [f_i, f_j, f_k, f_l];
    let forces_trans = [f_i_t, f_j_t, f_k_t, f_l_t];

    for a in 0..4 {
        for axis in 0..3 {
            let got = forces_trans[a][axis];
            let expected = forces_orig[a][axis];
            assert!(
                (got - expected).abs() <= tol,
                "force atom {a} axis {axis}: got {got}, expected {expected}"
            );
        }
    }
}

#[test]
fn test_zero_parameters() {
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.0, 1.0, 0.0];
    let r_l = [2.0, 1.0, 0.0];
    let v1 = 0.0;
    let v2 = 0.0;
    let v3 = 0.0;
    let v4 = 0.0;

    let e = energy_value(r_i, r_j, r_k, r_l, v1, v2, v3, v4);
    let tol = 1e-12;
    assert!(e.abs() <= tol, "energy: got {e}, expected 0.0");

    let (f_i, f_j, f_k, f_l) =
        opls_fourier_dihedral_forces(r_i, r_j, r_k, r_l, v1, v2, v3, v4).expect("forces failed");
    let forces = [f_i, f_j, f_k, f_l];
    for a in 0..4 {
        for axis in 0..3 {
            let got = forces[a][axis];
            assert!(
                got.abs() <= tol,
                "force atom {a} axis {axis}: got {got}, expected 0.0"
            );
        }
    }

    let w =
        opls_fourier_dihedral_virial(r_i, r_j, r_k, r_l, v1, v2, v3, v4).expect("virial failed");
    for a in 0..3 {
        for b in 0..3 {
            let got = w[a][b];
            assert!(
                got.abs() <= tol,
                "virial [{a}][{b}]: got {got}, expected 0.0"
            );
        }
    }
}

#[test]
fn test_collinear_error() {
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [2.0, 0.0, 0.0];
    let r_l = [3.0, 0.0, 0.0];
    let v1 = 1.0;
    let v2 = 0.0;
    let v3 = 0.0;
    let v4 = 0.0;

    let res = opls_fourier_dihedral_energy(r_i, r_j, r_k, r_l, v1, v2, v3, v4);
    assert!(
        matches!(res, Err(ucmd_ff::errors::UcmdError::CoincidentAtoms)),
        "expected CoincidentAtoms, got {:?}",
        res
    );
}

#[test]
fn test_energy_v2_term_at_phi_pi() {
    // Geometry from the analysis: phi = pi.
    // Correct energy: 0.5*(1*(1+cos(pi)) + 2*(1-cos(2*pi)) + 3*(1+cos(3*pi)) + 4*(1-cos(4*pi)))
    //                = 0.5*(0 + 0 + 3*(-1) + 4*(0)) = -1.5
    // Wrong version (phi2 = 1.0*phi instead of 2.0*phi):
    //                = 0.5*(0 + 2*(1-cos(pi)) + 3*(-1) + 4*(0)) = 0.5
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.0, 1.0, 0.0];
    let r_l = [2.0, 1.0, 0.0];
    let v1 = 1.0;
    let v2 = 2.0;
    let v3 = 3.0;
    let v4 = 4.0;

    let got = energy_value(r_i, r_j, r_k, r_l, v1, v2, v3, v4);

    // Compute expected from the specification formula with phi = pi.
    let phi: f64 = std::f64::consts::PI;
    let expected = 0.5
        * (v1 * (1.0 + phi.cos())
            + v2 * (1.0 - (2.0 * phi).cos())
            + v3 * (1.0 + (3.0 * phi).cos())
            + v4 * (1.0 - (4.0 * phi).cos()));

    let tol = 1e-12;
    assert!(
        (got - expected).abs() <= tol,
        "energy: got {got}, expected {expected}"
    );
}

#[test]
fn test_v4_energy_at_phi_pi_over_4() {
    // phi = pi/4, so 4*phi = pi, and 1 - cos(pi) = 2.
    // E = 0.5 * v4 * (1 - cos(4*phi)) = 0.5 * 1.0 * 2.0 = 1.0
    // This test fails if the V4 term uses cos(3*phi) instead of cos(4*phi),
    // or uses (1+cos) instead of (1-cos).
    let sqrt2_over_2: f64 = std::f64::consts::SQRT_2 / 2.0;
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.0, 1.0, 0.0];
    let r_l = [1.0 - sqrt2_over_2, 1.0, sqrt2_over_2];
    let v1 = 0.0;
    let v2 = 0.0;
    let v3 = 0.0;
    let v4 = 1.0;

    let got = energy_value(r_i, r_j, r_k, r_l, v1, v2, v3, v4);

    let phi: f64 = std::f64::consts::PI / 4.0;
    let expected = 0.5 * v4 * (1.0 - (4.0 * phi).cos());

    let tol = 1e-12;
    assert!(
        (got - expected).abs() <= tol,
        "energy: got {got}, expected {expected}"
    );
}

#[test]
fn test_v4_energy_at_phi_pi_over_3() {
    // phi = pi/3, so 4*phi = 4*pi/3, and 1 - cos(4*pi/3) = 1 - (-1/2) = 3/2.
    // E = 0.5 * v4 * (1 - cos(4*phi)) = 0.5 * 1.0 * 1.5 = 0.75
    // This test fails if the V4 term uses cos(3*phi) instead of cos(4*phi),
    // or uses (1+cos) instead of (1-cos).
    let sqrt3_over_2: f64 = (3.0_f64).sqrt() / 2.0;
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.0, 1.0, 0.0];
    let r_l = [1.0 - sqrt3_over_2, 1.0, 0.5];
    let v1 = 0.0;
    let v2 = 0.0;
    let v3 = 0.0;
    let v4 = 1.0;

    let got = energy_value(r_i, r_j, r_k, r_l, v1, v2, v3, v4);

    let phi: f64 = std::f64::consts::PI / 3.0;
    let expected = 0.5 * v4 * (1.0 - (4.0 * phi).cos());

    let tol = 1e-12;
    assert!(
        (got - expected).abs() <= tol,
        "energy: got {got}, expected {expected}"
    );
}

#[test]
fn test_forces_v4_at_phi_pi_over_8() {
    // phi = pi/8, so 4*phi = pi/2, sin(4*phi) = 1.
    // dE_dphi = 0.5 * 4 * v4 * sin(4*phi) = 0.5 * 4 * 1.0 * 1.0 = 2.0
    // The wrong version (factor 1.0 instead of 0.5) gives dE_dphi = 4.0.
    // This test checks forces via finite difference of energy, so it
    // catches the factor-of-2 error in dE_dphi.
    let cos_pi_8: f64 = (2.0 + (2.0_f64).sqrt()) / 2.0;
    let sin_pi_8: f64 = (2.0 - (2.0_f64).sqrt()) / 2.0;
    let r_i = [0.0, 0.0, 0.0];
    let r_j = [1.0, 0.0, 0.0];
    let r_k = [1.0, 1.0, 0.0];
    let r_l = [1.0 - cos_pi_8, 1.0, sin_pi_8];
    let v1 = 0.0;
    let v2 = 0.0;
    let v3 = 0.0;
    let v4 = 1.0;
    let h = 1e-5;

    let (f_i, f_j, f_k, f_l) =
        opls_fourier_dihedral_forces(r_i, r_j, r_k, r_l, v1, v2, v3, v4).expect("forces failed");

    let atoms = [r_i, r_j, r_k, r_l];
    let forces = [f_i, f_j, f_k, f_l];

    for a in 0..4 {
        for axis in 0..3 {
            let pos_plus = perturb_coord(atoms[a], axis, h);
            let pos_minus = perturb_coord(atoms[a], axis, -h);

            let (ri, rj, rk, rl) = match a {
                0 => (pos_plus, r_j, r_k, r_l),
                1 => (r_i, pos_plus, r_k, r_l),
                2 => (r_i, r_j, pos_plus, r_l),
                _ => (r_i, r_j, r_k, pos_plus),
            };
            let (ri_m, rj_m, rk_m, rl_m) = match a {
                0 => (pos_minus, r_j, r_k, r_l),
                1 => (r_i, pos_minus, r_k, r_l),
                2 => (r_i, r_j, pos_minus, r_l),
                _ => (r_i, r_j, r_k, pos_minus),
            };

            let e_plus = energy_value(ri, rj, rk, rl, v1, v2, v3, v4);
            let e_minus = energy_value(ri_m, rj_m, rk_m, rl_m, v1, v2, v3, v4);
            let expected = -(e_plus - e_minus) / (2.0 * h);

            let got = forces[a][axis];
            let tol = 1e-6;
            assert!(
                (got - expected).abs() <= tol,
                "atom {a} axis {axis}: got {got}, expected {expected}"
            );
        }
    }
}
