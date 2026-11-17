//! Tests of the `coulomb_pair` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::constants::COULOMB_KJ_MOL_NM_E2;
use ucmd_ff::coulomb_pair::*;
use ucmd_ff::errors::UcmdError;

#[test]
fn test_energy_unit_charges_one_nm() {
    let q_i = 1.0_f64;
    let q_j = 1.0_f64;
    let r_i: [f64; 3] = [1.0, 0.0, 0.0];
    let r_j: [f64; 3] = [0.0, 0.0, 0.0];

    let got = coulomb_pair_energy(q_i, q_j, r_i, r_j).unwrap();
    let expected = COULOMB_KJ_MOL_NM_E2 * q_i * q_j / 1.0;
    assert!(
        (got - expected).abs() <= 1e-15,
        "got {got}, expected {expected}"
    );
}

#[test]
fn test_force_finite_difference() {
    let q_i = 0.5_f64;
    let q_j = -0.3_f64;
    let r_i: [f64; 3] = [0.7, 0.2, -0.4];
    let r_j: [f64; 3] = [0.1, 0.5, 0.3];
    let h = 1e-8_f64;
    let tol = 1e-6_f64;

    let (f_i, f_j) = coulomb_pair_forces(q_i, q_j, r_i, r_j).unwrap();

    for a in 0..3 {
        // Finite difference for F_i
        let r_i_plus = perturb(r_i, a, h);
        let r_i_minus = perturb(r_i, a, -h);
        let e_plus = coulomb_pair_energy(q_i, q_j, r_i_plus, r_j).unwrap();
        let e_minus = coulomb_pair_energy(q_i, q_j, r_i_minus, r_j).unwrap();
        let fd_i = -(e_plus - e_minus) / (2.0 * h);
        assert!(
            (f_i[a] - fd_i).abs() <= tol,
            "F_i[{a}]: got {f_i:?}, fd {fd_i}"
        );

        // Finite difference for F_j
        let r_j_plus = perturb(r_j, a, h);
        let r_j_minus = perturb(r_j, a, -h);
        let e_plus_j = coulomb_pair_energy(q_i, q_j, r_i, r_j_plus).unwrap();
        let e_minus_j = coulomb_pair_energy(q_i, q_j, r_i, r_j_minus).unwrap();
        let fd_j = -(e_plus_j - e_minus_j) / (2.0 * h);
        assert!(
            (f_j[a] - fd_j).abs() <= tol,
            "F_j[{a}]: got {f_j:?}, fd {fd_j}"
        );
    }
}

// corrected after review: off-diagonal virial elements are non-zero for central forces; diagonal strain only probes diagonal elements
#[test]
fn test_virial_finite_difference_strain() {
    let q_i = 1.0_f64;
    let q_j = 1.0_f64;
    let r_i: [f64; 3] = [1.2, 0.5, 0.3];
    let r_j: [f64; 3] = [0.4, 0.1, 0.8];
    let eps = 1e-8_f64;
    let tol = 1e-6_f64;

    let w = coulomb_pair_virial(q_i, q_j, r_i, r_j).unwrap();

    for a in 0..3 {
        // Diagonal: finite difference under strain along axis a
        let r_i_strain_plus = strain_axis(r_i, a, eps);
        let r_j_strain_plus = strain_axis(r_j, a, eps);
        let r_i_strain_minus = strain_axis(r_i, a, -eps);
        let r_j_strain_minus = strain_axis(r_j, a, -eps);

        let e_plus = coulomb_pair_energy(q_i, q_j, r_i_strain_plus, r_j_strain_plus).unwrap();
        let e_minus = coulomb_pair_energy(q_i, q_j, r_i_strain_minus, r_j_strain_minus).unwrap();
        let fd = -(e_plus - e_minus) / (2.0 * eps);

        assert!(
            (w[a][a] - fd).abs() <= tol,
            "W[{a}][{a}]: got {w:?}, fd {fd}"
        );
    }

    // Symmetry check: W[a][b] == W[b][a] for all a, b
    for a in 0..3 {
        for b in 0..3 {
            assert!(
                (w[a][b] - w[b][a]).abs() <= tol,
                "W[{}][{}] = {}, W[{}][{}] = {}, expected equal",
                a,
                b,
                w[a][b],
                b,
                a,
                w[b][a]
            );
        }
    }
}

#[test]
fn test_newtons_third_law() {
    let q_i = 2.0_f64;
    let q_j = -1.5_f64;
    let r_i: [f64; 3] = [0.3, 1.1, -0.7];
    let r_j: [f64; 3] = [1.9, 0.4, 0.2];
    let tol = 1e-15_f64;

    let (f_i, f_j) = coulomb_pair_forces(q_i, q_j, r_i, r_j).unwrap();

    for a in 0..3 {
        let sum = f_i[a] + f_j[a];
        assert!(sum.abs() <= tol, "F_i[{a}] + F_j[{a}] = {sum}, expected 0");
    }
}

#[test]
fn test_zero_charge_gives_zero() {
    let q_i = 0.0_f64;
    let q_j = 1.0_f64;
    let r_i: [f64; 3] = [0.5, 0.5, 0.5];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];

    // Energy
    let e = coulomb_pair_energy(q_i, q_j, r_i, r_j).unwrap();
    assert!(e.abs() <= 0.0, "energy: got {e}, expected 0.0");

    // Forces
    let (f_i, f_j) = coulomb_pair_forces(q_i, q_j, r_i, r_j).unwrap();
    for a in 0..3 {
        assert!(f_i[a].abs() <= 0.0, "F_i[{a}]: got {f_i:?}, expected 0.0");
        assert!(f_j[a].abs() <= 0.0, "F_j[{a}]: got {f_j:?}, expected 0.0");
    }

    // Virial
    let w = coulomb_pair_virial(q_i, q_j, r_i, r_j).unwrap();
    for a in 0..3 {
        for b in 0..3 {
            assert!(w[a][b].abs() <= 0.0, "W[{a}][{b}]: got {w:?}, expected 0.0");
        }
    }
}

#[test]
fn test_coincident_atoms_error() {
    let q_i = 1.0_f64;
    let q_j = 1.0_f64;
    let r_i: [f64; 3] = [0.5, 0.5, 0.5];
    let r_j: [f64; 3] = [0.5, 0.5, 0.5];

    let e = coulomb_pair_energy(q_i, q_j, r_i, r_j);
    assert!(
        matches!(e, Err(UcmdError::CoincidentAtoms)),
        "energy: {e:?}"
    );

    let f = coulomb_pair_forces(q_i, q_j, r_i, r_j);
    assert!(
        matches!(f, Err(UcmdError::CoincidentAtoms)),
        "forces: {f:?}"
    );

    let w = coulomb_pair_virial(q_i, q_j, r_i, r_j);
    assert!(
        matches!(w, Err(UcmdError::CoincidentAtoms)),
        "virial: {w:?}"
    );
}

#[test]
fn test_nonfinite_coordinate_error() {
    let q_i = 1.0_f64;
    let q_j = 1.0_f64;
    let r_i: [f64; 3] = [f64::NAN, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];

    let e = coulomb_pair_energy(q_i, q_j, r_i, r_j);
    assert!(
        matches!(e, Err(UcmdError::NonFiniteCoordinate)),
        "energy: {e:?}"
    );

    let f = coulomb_pair_forces(q_i, q_j, r_i, r_j);
    assert!(
        matches!(f, Err(UcmdError::NonFiniteCoordinate)),
        "forces: {f:?}"
    );

    let w = coulomb_pair_virial(q_i, q_j, r_i, r_j);
    assert!(
        matches!(w, Err(UcmdError::NonFiniteCoordinate)),
        "virial: {w:?}"
    );
}

/// Perturbs component `axis` of vector `v` by `delta`.
fn perturb(v: [f64; 3], axis: usize, delta: f64) -> [f64; 3] {
    let mut out = v;
    out[axis] += delta;
    out
}

/// Scales component `axis` of vector `v` by `(1.0 + eps)`.
fn strain_axis(v: [f64; 3], axis: usize, eps: f64) -> [f64; 3] {
    let mut out = v;
    out[axis] *= 1.0 + eps;
    out
}
