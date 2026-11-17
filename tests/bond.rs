//! Tests of the `bond` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::bond::*;
use ucmd_ff::errors::*;

/// Helper: add two 3D vectors.
fn vec_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Helper: subtract vector b from vector a.
fn vec_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Helper: dot product of two vectors.
fn vec_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Helper: Euclidean norm of a vector.
fn vec_norm(v: [f64; 3]) -> f64 {
    (vec_dot(v, v)).sqrt()
}

#[test]
fn test_bond_energy_equilibrium() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [0.1, 0.0, 0.0];
    let r0_nm: f64 = 0.1;
    let k_kjmol_nm2: f64 = 100.0;

    let d = vec_sub(r_i, r_j);
    let r = vec_norm(d);
    let expected = 0.5 * k_kjmol_nm2 * (r - r0_nm) * (r - r0_nm);

    let result = bond_energy(r_i, r_j, r0_nm, k_kjmol_nm2);
    let got = match result {
        Ok(val) => val,
        Err(e) => panic!("Expected Ok, got Err: {:?}", e),
    };

    let tol: f64 = 1e-15;
    assert!(
        (got - expected).abs() <= tol,
        "got {}, expected {}",
        got,
        expected
    );
}

#[test]
fn test_bond_energy_stretched() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [0.2, 0.0, 0.0];
    let r0_nm: f64 = 0.1;
    let k_kjmol_nm2: f64 = 100.0;

    let d = vec_sub(r_i, r_j);
    let r = vec_norm(d);
    let expected = 0.5 * k_kjmol_nm2 * (r - r0_nm) * (r - r0_nm);

    let result = bond_energy(r_i, r_j, r0_nm, k_kjmol_nm2);
    let got = match result {
        Ok(val) => val,
        Err(e) => panic!("Expected Ok, got Err: {:?}", e),
    };

    let tol: f64 = 1e-12;
    assert!(
        (got - expected).abs() <= tol,
        "got {}, expected {}",
        got,
        expected
    );
}

#[test]
fn test_bond_forces_finite_difference() {
    let r_i: [f64; 3] = [0.1, 0.2, 0.3];
    let r_j: [f64; 3] = [0.4, 0.1, 0.2];
    let r0_nm: f64 = 0.5;
    let k_kjmol_nm2: f64 = 200.0;
    let h: f64 = 1e-5;

    let result = bond_forces(r_i, r_j, r0_nm, k_kjmol_nm2);
    let (f_i, f_j) = match result {
        Ok(val) => val,
        Err(e) => panic!("Expected Ok, got Err: {:?}", e),
    };

    // Check forces on r_i
    for axis in 0..3 {
        let mut r_i_plus = r_i;
        let mut r_i_minus = r_i;
        r_i_plus[axis] += h;
        r_i_minus[axis] -= h;

        let e_plus = match bond_energy(r_i_plus, r_j, r0_nm, k_kjmol_nm2) {
            Ok(val) => val,
            Err(e) => panic!("Expected Ok, got Err: {:?}", e),
        };
        let e_minus = match bond_energy(r_i_minus, r_j, r0_nm, k_kjmol_nm2) {
            Ok(val) => val,
            Err(e) => panic!("Expected Ok, got Err: {:?}", e),
        };

        let expected_fd = -(e_plus - e_minus) / (2.0 * h);
        let got = f_i[axis];
        let tol: f64 = 1e-5;

        // Use relative tolerance, but if expected is close to zero, use absolute
        let tol_eff = if expected_fd.abs() < 1e-10 {
            tol
        } else {
            tol * expected_fd.abs()
        };

        assert!(
            (got - expected_fd).abs() <= tol_eff,
            "Force on i axis {}: got {}, expected {}",
            axis,
            got,
            expected_fd
        );
    }

    // Check forces on r_j
    for axis in 0..3 {
        let mut r_j_plus = r_j;
        let mut r_j_minus = r_j;
        r_j_plus[axis] += h;
        r_j_minus[axis] -= h;

        let e_plus = match bond_energy(r_i, r_j_plus, r0_nm, k_kjmol_nm2) {
            Ok(val) => val,
            Err(e) => panic!("Expected Ok, got Err: {:?}", e),
        };
        let e_minus = match bond_energy(r_i, r_j_minus, r0_nm, k_kjmol_nm2) {
            Ok(val) => val,
            Err(e) => panic!("Expected Ok, got Err: {:?}", e),
        };

        let expected_fd = -(e_plus - e_minus) / (2.0 * h);
        let got = f_j[axis];
        let tol: f64 = 1e-5;

        let tol_eff = if expected_fd.abs() < 1e-10 {
            tol
        } else {
            tol * expected_fd.abs()
        };

        assert!(
            (got - expected_fd).abs() <= tol_eff,
            "Force on j axis {}: got {}, expected {}",
            axis,
            got,
            expected_fd
        );
    }
}

#[test]
fn test_bond_virial_strain() {
    let r_i: [f64; 3] = [0.1, 0.2, 0.3];
    let r_j: [f64; 3] = [0.4, 0.1, 0.2];
    let r0_nm: f64 = 0.5;
    let k_kjmol_nm2: f64 = 200.0;
    let eps: f64 = 1e-5;

    let result = bond_virial(r_i, r_j, r0_nm, k_kjmol_nm2);
    let w = match result {
        Ok(val) => val,
        Err(e) => panic!("Expected Ok, got Err: {:?}", e),
    };

    // Components to check: (0,0), (1,0), (2,0), (2,1)
    let components = [(0, 0), (1, 0), (2, 0), (2, 1)];

    for (a, b) in components {
        // Deform positions: r'_a = r_a + eps * r_b (if a != b) or r_a * (1 + eps) (if a == b)
        // This corresponds to r' = (I + eps) r where eps is a matrix with eps_ab = eps

        let mut r_i_plus = r_i;
        let mut r_j_plus = r_j;
        let mut r_i_minus = r_i;
        let mut r_j_minus = r_j;

        // Apply deformation for +eps
        if a == b {
            r_i_plus[a] *= 1.0 + eps;
            r_j_plus[a] *= 1.0 + eps;
        } else {
            r_i_plus[a] += eps * r_i[b];
            r_j_plus[a] += eps * r_j[b];
        }

        // Apply deformation for -eps
        if a == b {
            r_i_minus[a] *= 1.0 - eps;
            r_j_minus[a] *= 1.0 - eps;
        } else {
            r_i_minus[a] -= eps * r_i[b];
            r_j_minus[a] -= eps * r_j[b];
        }

        let e_plus = match bond_energy(r_i_plus, r_j_plus, r0_nm, k_kjmol_nm2) {
            Ok(val) => val,
            Err(e) => panic!("Expected Ok, got Err: {:?}", e),
        };
        let e_minus = match bond_energy(r_i_minus, r_j_minus, r0_nm, k_kjmol_nm2) {
            Ok(val) => val,
            Err(e) => panic!("Expected Ok, got Err: {:?}", e),
        };

        let expected_fd = -(e_plus - e_minus) / (2.0 * eps);
        let got = w[a][b];
        let tol: f64 = 1e-5;

        let tol_eff = if expected_fd.abs() < 1e-10 {
            tol
        } else {
            tol * expected_fd.abs()
        };

        assert!(
            (got - expected_fd).abs() <= tol_eff,
            "Virial component ({}, {}): got {}, expected {}",
            a,
            b,
            got,
            expected_fd
        );
    }
}

#[test]
fn test_bond_momentum_conservation() {
    let r_i: [f64; 3] = [0.1, 0.2, 0.3];
    let r_j: [f64; 3] = [0.4, 0.1, 0.2];
    let r0_nm: f64 = 0.5;
    let k_kjmol_nm2: f64 = 200.0;

    let result = bond_forces(r_i, r_j, r0_nm, k_kjmol_nm2);
    let (f_i, f_j) = match result {
        Ok(val) => val,
        Err(e) => panic!("Expected Ok, got Err: {:?}", e),
    };

    let sum = vec_add(f_i, f_j);
    let norm_sum = vec_norm(sum);
    let tol: f64 = 1e-10;

    assert!(
        norm_sum <= tol,
        "Sum of forces norm: got {}, expected <= {}",
        norm_sum,
        tol
    );
}

#[test]
fn test_bond_coincident_error() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [0.0, 0.0, 0.0];
    let r0_nm: f64 = 0.1;
    let k_kjmol_nm2: f64 = 100.0;

    let result = bond_energy(r_i, r_j, r0_nm, k_kjmol_nm2);
    assert!(
        matches!(result, Err(UcmdError::CoincidentAtoms)),
        "Expected Err(CoincidentAtoms), got {:?}",
        result
    );
}

#[test]
fn test_bond_nonfinite_error() {
    let r_i: [f64; 3] = [f64::NAN, 0.0, 0.0];
    let r_j: [f64; 3] = [0.0, 0.0, 0.0];
    let r0_nm: f64 = 0.1;
    let k_kjmol_nm2: f64 = 100.0;

    let result = bond_energy(r_i, r_j, r0_nm, k_kjmol_nm2);
    assert!(
        matches!(result, Err(UcmdError::NonFiniteCoordinate)),
        "Expected Err(NonFiniteCoordinate), got {:?}",
        result
    );
}
