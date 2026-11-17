//! Tests of the `minimum_image` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::cell::Cell;
use ucmd_ff::errors::UcmdError;
use ucmd_ff::minimum_image::*;

#[test]
fn test_cubic_cutoff_pass() {
    let cell = Cell {
        ax: 2.0,
        bx: 0.0,
        by: 2.0,
        cx: 0.0,
        cy: 0.0,
        cz: 2.0,
    };
    let cutoff_nm: f64 = 0.999999;
    let result = validate_cutoff(&cell, cutoff_nm);
    assert!(result.is_ok(), "expected Ok(()), got {:?}", result);
}

#[test]
fn test_cubic_cutoff_fail() {
    let cell = Cell {
        ax: 2.0,
        bx: 0.0,
        by: 2.0,
        cx: 0.0,
        cy: 0.0,
        cz: 2.0,
    };
    let cutoff_nm: f64 = 1.000001;
    let result = validate_cutoff(&cell, cutoff_nm);
    assert!(
        matches!(result, Err(UcmdError::InvalidCell)),
        "expected Err(InvalidCell), got {:?}",
        result
    );
}

#[test]
fn test_sheared_cutoff_fail() {
    let cell = Cell {
        ax: 1.0,
        bx: 0.9,
        by: 1.0,
        cx: 0.9,
        cy: 0.9,
        cz: 1.0,
    };
    // Compute h_min from the specification's arithmetic:
    // V = ax * by * cz = 1.0
    // |b×c| = sqrt((by*cz)^2 + (bx*cz)^2 + (bx*cy - by*cx)^2)
    //       = sqrt(1.0^2 + 0.9^2 + (0.9*0.9 - 1.0*0.9)^2)
    //       = sqrt(1.0 + 0.81 + 0.0081) = sqrt(1.8181)
    // h_a = V / |b×c| = 1.0 / sqrt(1.8181)
    let v: f64 = 1.0 * 1.0 * 1.0;
    let b_cross_c_norm: f64 = (1.0_f64 * 1.0_f64).powi(2)
        + (0.9_f64 * 1.0_f64).powi(2)
        + (0.9_f64 * 0.9_f64 - 1.0_f64 * 0.9_f64).powi(2);
    let h_a: f64 = v / b_cross_c_norm.sqrt();
    let cutoff_nm: f64 = h_a / 2.0 + 0.01;
    let result = validate_cutoff(&cell, cutoff_nm);
    assert!(
        matches!(result, Err(UcmdError::InvalidCell)),
        "expected Err(InvalidCell), got {:?}",
        result
    );
}

#[test]
fn test_invalid_cell_error() {
    let cell = Cell {
        ax: 0.0,
        bx: 0.0,
        by: 1.0,
        cx: 0.0,
        cy: 0.0,
        cz: 1.0,
    };
    let cutoff_nm: f64 = 0.5;
    let result = validate_cutoff(&cell, cutoff_nm);
    assert!(
        matches!(result, Err(UcmdError::InvalidCell)),
        "expected Err(InvalidCell), got {:?}",
        result
    );
}

#[test]
fn test_heights_cubic() {
    let cell = Cell {
        ax: 3.0,
        bx: 0.0,
        by: 3.0,
        cx: 0.0,
        cy: 0.0,
        cz: 3.0,
    };
    let result = perpendicular_heights(&cell);
    assert!(result.is_ok(), "expected Ok, got {:?}", result);
    let heights = result.unwrap();
    let expected: f64 = 3.0;
    let tol: f64 = 1e-12;
    assert!(
        (heights[0] - expected).abs() <= tol,
        "got {heights:?}, expected [{expected}, {expected}, {expected}]"
    );
    assert!(
        (heights[1] - expected).abs() <= tol,
        "got {heights:?}, expected [{expected}, {expected}, {expected}]"
    );
    assert!(
        (heights[2] - expected).abs() <= tol,
        "got {heights:?}, expected [{expected}, {expected}, {expected}]"
    );
}

#[test]
fn test_heights_sheared() {
    let cell = Cell {
        ax: 1.0,
        bx: 0.9,
        by: 1.0,
        cx: 0.9,
        cy: 0.9,
        cz: 1.0,
    };
    let result = perpendicular_heights(&cell);
    assert!(result.is_ok(), "expected Ok, got {:?}", result);
    let heights = result.unwrap();
    // V = ax * by * cz = 1.0
    let v: f64 = 1.0 * 1.0 * 1.0;
    // |b×c| = sqrt((by*cz)^2 + (bx*cz)^2 + (bx*cy - by*cx)^2)
    let bcc: f64 = (1.0_f64 * 1.0_f64).powi(2)
        + (0.9_f64 * 1.0_f64).powi(2)
        + (0.9_f64 * 0.9_f64 - 1.0_f64 * 0.9_f64).powi(2);
    let expected_h_a: f64 = v / bcc.sqrt();
    // |a×c| = sqrt((ax*cz)^2 + (ax*cy)^2)
    let acc: f64 = (1.0_f64 * 1.0_f64).powi(2) + (1.0_f64 * 0.9_f64).powi(2);
    let expected_h_b: f64 = v / acc.sqrt();
    // |a×b| = ax * by
    let expected_h_c: f64 = v / (1.0 * 1.0);
    let tol: f64 = 1e-12;
    assert!(
        (heights[0] - expected_h_a).abs() <= tol,
        "got {heights:?}, expected h_a={expected_h_a}"
    );
    assert!(
        (heights[1] - expected_h_b).abs() <= tol,
        "got {heights:?}, expected h_b={expected_h_b}"
    );
    assert!(
        (heights[2] - expected_h_c).abs() <= tol,
        "got {heights:?}, expected h_c={expected_h_c}"
    );
}

#[test]
fn test_cubic_cutoff_exact_boundary() {
    let cell = Cell {
        ax: 2.0,
        bx: 0.0,
        by: 2.0,
        cx: 0.0,
        cy: 0.0,
        cz: 2.0,
    };
    // Arithmetic: V = 2.0 * 2.0 * 2.0 = 8.0
    // |b×c| = |cross([0,2,0],[0,0,2])| = |[4,0,0]| = 4.0
    // h_a = 8.0 / 4.0 = 2.0
    // |a×c| = |cross([2,0,0],[0,0,2])| = |[0,-4,0]| = 4.0
    // h_b = 8.0 / 4.0 = 2.0
    // |a×b| = |cross([2,0,0],[0,2,0])| = |[0,0,4]| = 4.0
    // h_c = 8.0 / 4.0 = 2.0
    // h_min = 2.0, cutoff = h_min / 2.0 = 1.0
    let cutoff_nm: f64 = 2.0 / 2.0;
    let result = validate_cutoff(&cell, cutoff_nm);
    assert!(
        matches!(result, Err(UcmdError::InvalidCell)),
        "expected Err(InvalidCell) for cutoff = h_min/2 = {cutoff_nm}, got {:?}",
        result
    );
}

#[test]
fn test_triclinic_cutoff_exact_boundary() {
    let cell = Cell {
        ax: 1.0,
        bx: 0.9,
        by: 1.0,
        cx: 0.9,
        cy: 0.9,
        cz: 1.0,
    };
    // Arithmetic:
    // V = ax * by * cz = 1.0 * 1.0 * 1.0 = 1.0
    // a = [1.0, 0.0, 0.0], b = [0.9, 1.0, 0.0], c = [0.9, 0.9, 1.0]
    // b×c = [1.0*1.0 - 0.0*0.9, 0.0*0.9 - 0.9*1.0, 0.9*0.9 - 1.0*0.9]
    //      = [1.0, -0.9, -0.09]
    // |b×c| = sqrt(1.0^2 + 0.9^2 + 0.09^2) = sqrt(1.0 + 0.81 + 0.0081) = sqrt(1.8181)
    // h_a = 1.0 / sqrt(1.8181)
    // a×c = [0.0*1.0 - 0.0*0.9, 0.0*0.9 - 1.0*1.0, 1.0*0.9 - 0.0*0.9]
    //      = [0.0, -1.0, 0.9]
    // |a×c| = sqrt(0.0^2 + 1.0^2 + 0.9^2) = sqrt(1.81)
    // h_b = 1.0 / sqrt(1.81)
    // a×b = [0.0*0.0 - 0.0*1.0, 0.0*0.9 - 1.0*0.0, 1.0*1.0 - 0.0*0.9]
    //      = [0.0, 0.0, 1.0]
    // |a×b| = 1.0
    // h_c = 1.0 / 1.0 = 1.0
    // h_min = h_a = 1.0 / sqrt(1.8181)
    // cutoff = h_min / 2.0 = 1.0 / (2.0 * sqrt(1.8181))
    let v: f64 = 1.0 * 1.0 * 1.0;
    let bcc_sq: f64 = (1.0_f64 * 1.0_f64).powi(2)
        + (0.9_f64 * 1.0_f64).powi(2)
        + (0.9_f64 * 0.9_f64 - 1.0_f64 * 0.9_f64).powi(2);
    let h_a: f64 = v / bcc_sq.sqrt();
    let acc_sq: f64 = (1.0_f64 * 1.0_f64).powi(2) + (1.0_f64 * 0.9_f64).powi(2);
    let h_b: f64 = v / acc_sq.sqrt();
    let h_c: f64 = v / (1.0 * 1.0);
    let h_min: f64 = h_a.min(h_b).min(h_c);
    let cutoff_nm: f64 = h_min / 2.0;
    let result = validate_cutoff(&cell, cutoff_nm);
    assert!(
        matches!(result, Err(UcmdError::InvalidCell)),
        "expected Err(InvalidCell) for cutoff = h_min/2 = {cutoff_nm}, got {:?}",
        result
    );
}
