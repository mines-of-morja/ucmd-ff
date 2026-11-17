//! Tests of the `min_image_disp` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::cell::Cell;
use ucmd_ff::errors::UcmdError;
use ucmd_ff::min_image_disp::*;

#[test]
fn cubic_cell_wrap_across_a_face() {
    let cell = Cell {
        ax: 2.0,
        bx: 0.0,
        by: 2.0,
        cx: 0.0,
        cy: 0.0,
        cz: 2.0,
    };
    let r_i: [f64; 3] = [1.9, 0.0, 0.0];
    let r_j: [f64; 3] = [0.1, 0.0, 0.0];

    let (d_min, dist) = minimum_image_displacement(r_i, r_j, &cell).unwrap();

    // Raw displacement = (1.8, 0, 0). Fractional s_x = 1.8/2.0 = 0.9.
    // round(0.9) = 1. s - round = -0.1. d_min_x = 2.0 * (-0.1) = -0.2.
    let expected_dx: f64 = 2.0 * (-0.1);
    let expected_dy: f64 = 0.0;
    let expected_dz: f64 = 0.0;
    let expected_dist: f64 = expected_dx.abs();

    let tol: f64 = 1e-14;
    assert!(
        (d_min[0] - expected_dx).abs() <= tol,
        "got {d_min:?}, expected dx {expected_dx}"
    );
    assert!(
        (d_min[1] - expected_dy).abs() <= tol,
        "got {d_min:?}, expected dy {expected_dy}"
    );
    assert!(
        (d_min[2] - expected_dz).abs() <= tol,
        "got {d_min:?}, expected dz {expected_dz}"
    );
    assert!(
        (dist - expected_dist).abs() <= tol,
        "got {dist}, expected {expected_dist}"
    );
}

// corrected after review: d_min = d - n·C, not C·(s-n); for triclinic cell the x,y components are nonzero
#[test]
fn triclinic_cell_wrap_in_c_direction() {
    let cell = Cell {
        ax: 2.0,
        bx: 0.5,
        by: 2.0,
        cx: 0.5,
        cy: 0.5,
        cz: 2.0,
    };
    let r_i: [f64; 3] = [0.0, 0.0, 1.9];
    let r_j: [f64; 3] = [0.0, 0.0, 0.1];

    let (d_min, dist) = minimum_image_displacement(r_i, r_j, &cell).unwrap();

    // d = (0, 0, 1.8). C⁻¹·d gives s = (0, 0, 0.9).
    // round = (0, 0, 1). d_min = d - n·C = (0,0,1.8) - (0.5, 0.5, 2.0) = (-0.5, -0.5, -0.2).
    let expected_dx: f64 = -0.5;
    let expected_dy: f64 = -0.5;
    let expected_dz: f64 = -0.2;
    let expected_dist: f64 =
        (expected_dx * expected_dx + expected_dy * expected_dy + expected_dz * expected_dz).sqrt();

    let tol: f64 = 1e-14;
    assert!(
        (d_min[0] - expected_dx).abs() <= tol,
        "got {d_min:?}, expected dx {expected_dx}"
    );
    assert!(
        (d_min[1] - expected_dy).abs() <= tol,
        "got {d_min:?}, expected dy {expected_dy}"
    );
    assert!(
        (d_min[2] - expected_dz).abs() <= tol,
        "got {d_min:?}, expected dz {expected_dz}"
    );
    assert!(
        (dist - expected_dist).abs() <= tol,
        "got {dist}, expected {expected_dist}"
    );
}

#[test]
fn coincident_atoms_returns_error() {
    let cell = Cell {
        ax: 2.0,
        bx: 0.0,
        by: 2.0,
        cx: 0.0,
        cy: 0.0,
        cz: 2.0,
    };
    let r_i: [f64; 3] = [1.0, 1.0, 1.0];
    let r_j: [f64; 3] = [1.0, 1.0, 1.0];

    let result = minimum_image_displacement(r_i, r_j, &cell);
    assert!(
        matches!(result, Err(UcmdError::CoincidentAtoms)),
        "expected Err(CoincidentAtoms), got {result:?}"
    );
}

#[test]
fn coincident_via_lattice_translation() {
    let cell = Cell {
        ax: 2.0,
        bx: 0.0,
        by: 2.0,
        cx: 0.0,
        cy: 0.0,
        cz: 2.0,
    };
    let r_i: [f64; 3] = [3.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];

    let result = minimum_image_displacement(r_i, r_j, &cell);
    assert!(
        matches!(result, Err(UcmdError::CoincidentAtoms)),
        "expected Err(CoincidentAtoms), got {result:?}"
    );
}

#[test]
fn non_finite_coordinate_returns_error() {
    let cell = Cell {
        ax: 2.0,
        bx: 0.0,
        by: 2.0,
        cx: 0.0,
        cy: 0.0,
        cz: 2.0,
    };
    let r_i: [f64; 3] = [f64::NAN, 0.0, 0.0];
    let r_j: [f64; 3] = [0.0, 0.0, 0.0];

    let result = minimum_image_displacement(r_i, r_j, &cell);
    assert!(
        matches!(result, Err(UcmdError::NonFiniteCoordinate)),
        "expected Err(NonFiniteCoordinate), got {result:?}"
    );
}

#[test]
fn unwrapped_positions_many_cells_apart() {
    let cell = Cell {
        ax: 2.0,
        bx: 0.0,
        by: 2.0,
        cx: 0.0,
        cy: 0.0,
        cz: 2.0,
    };
    let r_i: [f64; 3] = [11.9, 0.0, 0.0];
    let r_j: [f64; 3] = [0.1, 0.0, 0.0];

    let (d_min, dist) = minimum_image_displacement(r_i, r_j, &cell).unwrap();

    // Raw d = (11.8, 0, 0). Fractional s_x = 11.8/2.0 = 5.9.
    // round(5.9) = 6. s - round = -0.1. d_min_x = 2.0 * (-0.1) = -0.2.
    let expected_dx: f64 = 2.0 * (-0.1);
    let expected_dy: f64 = 0.0;
    let expected_dz: f64 = 0.0;
    let expected_dist: f64 = expected_dx.abs();

    let tol: f64 = 1e-14;
    assert!(
        (d_min[0] - expected_dx).abs() <= tol,
        "got {d_min:?}, expected dx {expected_dx}"
    );
    assert!(
        (d_min[1] - expected_dy).abs() <= tol,
        "got {d_min:?}, expected dy {expected_dy}"
    );
    assert!(
        (d_min[2] - expected_dz).abs() <= tol,
        "got {d_min:?}, expected dz {expected_dz}"
    );
    assert!(
        (dist - expected_dist).abs() <= tol,
        "got {dist}, expected {expected_dist}"
    );
}

#[test]
fn invalid_cell_returns_error() {
    let cell = Cell {
        ax: -1.0,
        bx: 0.0,
        by: 2.0,
        cx: 0.0,
        cy: 0.0,
        cz: 2.0,
    };
    let r_i: [f64; 3] = [0.5, 0.0, 0.0];
    let r_j: [f64; 3] = [0.0, 0.0, 0.0];

    let result = minimum_image_displacement(r_i, r_j, &cell);
    assert!(
        matches!(result, Err(UcmdError::InvalidCell)),
        "expected Err(InvalidCell), got {result:?}"
    );
}
