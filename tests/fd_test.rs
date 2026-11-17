//! Tests of the `fd_test` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::angle::{angle_energy, angle_forces, angle_virial};
use ucmd_ff::bond::{bond_energy, bond_forces, bond_virial};
use ucmd_ff::cell::{Cell, cell_volume};
use ucmd_ff::fd_test::*;
use ucmd_ff::vec3;

fn ok<T>(r: Result<T, ucmd_ff::errors::UcmdError>) -> T {
    match r {
        Ok(v) => v,
        Err(_) => panic!("unexpected error"),
    }
}

#[test]
fn test_numerical_force_bond_matches_analytic() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [0.105, 0.0, 0.0];
    let r0_nm: f64 = 0.1012;
    let k_kjmol_nm2: f64 = 443153.3808;
    let atom_idx: usize = 0;
    let coord_idx: usize = 0;
    let h: f64 = 1e-5;

    let analytic_forces = ok(bond_forces(r_i, r_j, r0_nm, k_kjmol_nm2));
    let analytic_fx = analytic_forces.0[coord_idx];

    let numerical_fx = ok(numerical_force_bond(
        bond_energy,
        r_i,
        r_j,
        r0_nm,
        k_kjmol_nm2,
        atom_idx,
        coord_idx,
        h,
    ));

    let tol: f64 = 1e-3;
    assert!(
        (numerical_fx - analytic_fx).abs() <= tol,
        "got {numerical_fx}, expected {analytic_fx}"
    );
}

#[test]
fn test_numerical_virial_bond_matches_analytic() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [0.105, 0.0, 0.0];
    let r0_nm: f64 = 0.1012;
    let k_kjmol_nm2: f64 = 443153.3808;
    let cell = Cell {
        ax: 1.0,
        bx: 0.0,
        by: 1.0,
        cx: 0.0,
        cy: 0.0,
        cz: 1.0,
    };
    let strain_row: usize = 0;
    let strain_col: usize = 0;
    let h: f64 = 1e-5;

    let analytic_forces = ok(bond_forces(r_i, r_j, r0_nm, k_kjmol_nm2));
    let f_i = analytic_forces.0;
    let r_ij = vec3::sub(r_i, r_j);
    let analytic_w00 = r_ij[0] * f_i[0];

    let numerical_w00 = ok(numerical_virial_bond(
        bond_energy,
        r_i,
        r_j,
        r0_nm,
        k_kjmol_nm2,
        &cell,
        strain_row,
        strain_col,
        h,
    ));

    let tol: f64 = 1e-3;
    assert!(
        (numerical_w00 - analytic_w00).abs() <= tol,
        "got {numerical_w00}, expected {analytic_w00}"
    );
}

#[test]
fn test_compare_force_within_tolerance() {
    let analytic: [f64; 3] = [1.0, 2.0, 3.0];
    let numerical: [f64; 3] = [1.0001, 2.0001, 3.0001];
    let tol: f64 = 0.001;
    assert!(
        compare_force(analytic, numerical, tol),
        "compare_force should return true for vectors within tolerance"
    );
}

#[test]
fn test_compare_virial_within_tolerance() {
    let analytic: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
    let numerical: [[f64; 3]; 3] = [[1.0001, 0.0, 0.0], [0.0, 2.0001, 0.0], [0.0, 0.0, 3.0001]];
    let tol: f64 = 0.001;
    assert!(
        compare_virial(analytic, numerical, tol),
        "compare_virial should return true for matrices within tolerance"
    );
}

#[test]
fn test_numerical_force_angle_matches_analytic() {
    let r1: [f64; 3] = [0.1012, 0.0, 0.0];
    let r2: [f64; 3] = [0.0, 0.0, 0.0];
    let r3: [f64; 3] = [0.0, 0.1012, 0.0];
    let theta0: f64 = 1.9762;
    let k_kjmol_rad2: f64 = 317.5656;
    let atom_idx: usize = 0;
    let coord_idx: usize = 0;
    let h: f64 = 1e-5;

    let analytic_forces = ok(angle_forces(r1, r2, r3, theta0, k_kjmol_rad2));
    let analytic_f1x = analytic_forces.0[coord_idx];

    let numerical_f1x = ok(numerical_force_angle(
        angle_energy,
        r1,
        r2,
        r3,
        theta0,
        k_kjmol_rad2,
        atom_idx,
        coord_idx,
        h,
    ));

    let tol: f64 = 1e-3;
    assert!(
        (numerical_f1x - analytic_f1x).abs() <= tol,
        "got {numerical_f1x}, expected {analytic_f1x}"
    );
}

#[test]
fn test_numerical_virial_angle_matches_analytic() {
    let r1: [f64; 3] = [0.1012, 0.0, 0.0];
    let r2: [f64; 3] = [0.0, 0.0, 0.0];
    let r3: [f64; 3] = [0.0, 0.1012, 0.0];
    let theta0: f64 = 1.9762;
    let k_kjmol_rad2: f64 = 317.5656;
    let cell = Cell {
        ax: 1.0,
        bx: 0.0,
        by: 1.0,
        cx: 0.0,
        cy: 0.0,
        cz: 1.0,
    };
    let strain_row: usize = 0;
    let strain_col: usize = 1;
    let h: f64 = 1e-5;

    let analytic_forces = ok(angle_forces(r1, r2, r3, theta0, k_kjmol_rad2));
    let f1 = analytic_forces.0;
    let f3 = analytic_forces.2;
    let r_12 = vec3::sub(r1, r2);
    let r_23 = vec3::sub(r2, r3);
    let analytic_w01 = r_12[0] * f1[1] + r_23[0] * (-f3[1]);

    let numerical_w01 = ok(numerical_virial_angle(
        angle_energy,
        r1,
        r2,
        r3,
        theta0,
        k_kjmol_rad2,
        &cell,
        strain_row,
        strain_col,
        h,
    ));

    let tol: f64 = 1e-3;
    assert!(
        (numerical_w01 - analytic_w01).abs() <= tol,
        "got {numerical_w01}, expected {analytic_w01}"
    );
}

#[test]
fn test_generic_numerical_forces_bond() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [0.105, 0.0, 0.0];
    let r0_nm: f64 = 0.1012;
    let k_kjmol_nm2: f64 = 443153.3808;
    let h_nm: f64 = 1e-5;

    let positions: [[f64; 3]; 2] = [r_i, r_j];
    let analytic_tuple = ok(bond_forces(r_i, r_j, r0_nm, k_kjmol_nm2));
    let analytic_forces: [[f64; 3]; 2] = [analytic_tuple.0, analytic_tuple.1];

    let energy_closure = |pos: &[[f64; 3]]| -> Result<f64, ucmd_ff::errors::UcmdError> {
        bond_energy(pos[0], pos[1], r0_nm, k_kjmol_nm2)
    };

    let numerical_forces = ok(numerical_forces(&energy_closure, &positions, h_nm));

    let tol: f64 = 1e-3;
    for i in 0..2 {
        let max_diff = max_abs_difference_forces(&analytic_forces[i], &numerical_forces[i]);
        assert!(
            max_diff <= tol,
            "atom {i}: max diff {max_diff} exceeds tolerance {tol}"
        );
    }
}

#[test]
fn test_generic_numerical_forces_angle() {
    let r1: [f64; 3] = [0.1012, 0.0, 0.0];
    let r2: [f64; 3] = [0.0, 0.0, 0.0];
    let r3: [f64; 3] = [0.0, 0.1012, 0.0];
    let theta0: f64 = 1.9762;
    let k_kjmol_rad2: f64 = 317.5656;
    let h_nm: f64 = 1e-5;

    let positions: [[f64; 3]; 3] = [r1, r2, r3];
    let analytic_tuple = ok(angle_forces(r1, r2, r3, theta0, k_kjmol_rad2));
    let analytic_forces: [[f64; 3]; 3] = [analytic_tuple.0, analytic_tuple.1, analytic_tuple.2];

    let energy_closure = |pos: &[[f64; 3]]| -> Result<f64, ucmd_ff::errors::UcmdError> {
        angle_energy(pos[0], pos[1], pos[2], theta0, k_kjmol_rad2)
    };

    let numerical_forces = ok(numerical_forces(&energy_closure, &positions, h_nm));

    let tol: f64 = 1e-3;
    for i in 0..3 {
        let max_diff = max_abs_difference_forces(&analytic_forces[i], &numerical_forces[i]);
        assert!(
            max_diff <= tol,
            "atom {i}: max diff {max_diff} exceeds tolerance {tol}"
        );
    }
}

#[test]
fn test_generic_numerical_virial_bond() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [0.105, 0.0, 0.0];
    let r0_nm: f64 = 0.1012;
    let k_kjmol_nm2: f64 = 443153.3808;
    let eps: f64 = 1e-5;

    let positions: [[f64; 3]; 2] = [r_i, r_j];
    let analytic_virial = ok(bond_virial(r_i, r_j, r0_nm, k_kjmol_nm2));

    let energy_closure = |pos: &[[f64; 3]]| -> Result<f64, ucmd_ff::errors::UcmdError> {
        bond_energy(pos[0], pos[1], r0_nm, k_kjmol_nm2)
    };

    let numerical_virial = ok(numerical_virial(&energy_closure, &positions, eps));

    let max_diff = max_abs_difference_virial(&analytic_virial, &numerical_virial);
    let tol: f64 = 1e-3;
    assert!(
        max_diff <= tol,
        "max diff {max_diff} exceeds tolerance {tol}"
    );
}

#[test]
fn test_generic_numerical_virial_angle() {
    let r1: [f64; 3] = [0.1012, 0.0, 0.0];
    let r2: [f64; 3] = [0.0, 0.0, 0.0];
    let r3: [f64; 3] = [0.0, 0.1012, 0.0];
    let theta0: f64 = 1.9762;
    let k_kjmol_rad2: f64 = 317.5656;
    let eps: f64 = 1e-5;

    let positions: [[f64; 3]; 3] = [r1, r2, r3];
    let analytic_virial = ok(angle_virial(r1, r2, r3, theta0, k_kjmol_rad2));

    let energy_closure = |pos: &[[f64; 3]]| -> Result<f64, ucmd_ff::errors::UcmdError> {
        angle_energy(pos[0], pos[1], pos[2], theta0, k_kjmol_rad2)
    };

    let numerical_virial = ok(numerical_virial(&energy_closure, &positions, eps));

    let max_diff = max_abs_difference_virial(&analytic_virial, &numerical_virial);
    let tol: f64 = 1e-3;
    assert!(
        max_diff <= tol,
        "max diff {max_diff} exceeds tolerance {tol}"
    );
}

#[test]
fn test_numerical_virial_periodic_volume_energy() {
    let cell = Cell {
        ax: 1.0,
        bx: 0.3,
        by: 1.0,
        cx: 0.2,
        cy: 0.4,
        cz: 1.0,
    };
    let positions: [[f64; 3]; 1] = [[0.0, 0.0, 0.0]];
    let p: f64 = 2.5;
    let eps: f64 = 1e-5;

    let volume = ok(cell_volume(&cell));
    let expected_w00: f64 = -p * volume;

    let energy_closure = |_pos: &[[f64; 3]], c: &Cell| -> Result<f64, ucmd_ff::errors::UcmdError> {
        let v = ok(cell_volume(c));
        Ok(p * v)
    };

    let w00 = ok(numerical_virial_periodic(
        &energy_closure,
        &positions,
        &cell,
        0,
        0,
        eps,
    ));
    let tol: f64 = 1e-3;
    assert!(
        (w00 - expected_w00).abs() <= tol,
        "W_00: got {w00}, expected {expected_w00}"
    );

    let w10 = ok(numerical_virial_periodic(
        &energy_closure,
        &positions,
        &cell,
        1,
        0,
        eps,
    ));
    assert!((w10 - 0.0).abs() <= tol, "W_10: got {w10}, expected 0.0");
}

#[test]
fn test_numerical_virial_periodic_cell_independent() {
    let cell = Cell {
        ax: 1.0,
        bx: 0.3,
        by: 1.0,
        cx: 0.2,
        cy: 0.4,
        cz: 1.0,
    };
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [0.105, 0.0, 0.0];
    let r0_nm: f64 = 0.1012;
    let k_kjmol_nm2: f64 = 443153.3808;
    let positions: [[f64; 3]; 2] = [r_i, r_j];
    let eps: f64 = 1e-5;

    let energy_closure = |pos: &[[f64; 3]], _c: &Cell| -> Result<f64, ucmd_ff::errors::UcmdError> {
        bond_energy(pos[0], pos[1], r0_nm, k_kjmol_nm2)
    };

    let w00 = ok(numerical_virial_periodic(
        &energy_closure,
        &positions,
        &cell,
        0,
        0,
        eps,
    ));
    let w10 = ok(numerical_virial_periodic(
        &energy_closure,
        &positions,
        &cell,
        1,
        0,
        eps,
    ));

    let energy_closure_np = |pos: &[[f64; 3]]| -> Result<f64, ucmd_ff::errors::UcmdError> {
        bond_energy(pos[0], pos[1], r0_nm, k_kjmol_nm2)
    };
    let w00_np = ok(numerical_virial(&energy_closure_np, &positions, eps))[0][0];
    let w10_np = ok(numerical_virial(&energy_closure_np, &positions, eps))[1][0];

    let tol: f64 = 1e-3;
    assert!(
        (w00 - w00_np).abs() <= tol,
        "W_00: got {w00}, expected {w00_np}"
    );
    assert!(
        (w10 - w10_np).abs() <= tol,
        "W_10: got {w10}, expected {w10_np}"
    );
}

#[test]
fn test_numerical_virial_periodic_min_image_pair() {
    let cell = Cell {
        ax: 1.0,
        bx: 0.0,
        by: 1.0,
        cx: 0.0,
        cy: 0.0,
        cz: 1.0,
    };
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [0.9, 0.0, 0.0];
    let r0_nm: f64 = 0.1012;
    let k_kjmol_nm2: f64 = 443153.3808;
    let positions: [[f64; 3]; 2] = [r_i, r_j];
    let eps: f64 = 1e-5;

    let energy_closure = |pos: &[[f64; 3]], c: &Cell| -> Result<f64, ucmd_ff::errors::UcmdError> {
        let d = ok(ucmd_ff::cell::minimum_image(pos[0], pos[1], c));
        let r = vec3::norm(d);
        let dr = r - r0_nm;
        Ok(0.5 * k_kjmol_nm2 * dr * dr)
    };

    let w00 = ok(numerical_virial_periodic(
        &energy_closure,
        &positions,
        &cell,
        0,
        0,
        eps,
    ));
    let expected_w00: f64 = 53.178405696;
    let tol: f64 = 1e-3;
    assert!(
        (w00 - expected_w00).abs() <= tol,
        "W_00: got {w00}, expected {expected_w00}"
    );
}

#[test]
fn test_numerical_virial_periodic_invalid_component() {
    let cell = Cell {
        ax: 1.0,
        bx: 0.0,
        by: 1.0,
        cx: 0.0,
        cy: 0.0,
        cz: 1.0,
    };
    let positions: [[f64; 3]; 1] = [[0.0, 0.0, 0.0]];
    let eps: f64 = 1e-5;

    let energy_closure =
        |_pos: &[[f64; 3]], _c: &Cell| -> Result<f64, ucmd_ff::errors::UcmdError> { Ok(0.0) };

    // (0, 1) is not a valid lower-triangular component
    let result = numerical_virial_periodic(&energy_closure, &positions, &cell, 0, 1, eps);
    assert!(
        result.is_err(),
        "expected error for invalid component (0, 1), got {:?}",
        result
    );
}
