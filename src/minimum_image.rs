//! Computes perpendicular cell heights and validates that a cutoff radius
//! does not exceed half the smallest perpendicular cell height (R3.5).
//! Units: all length quantities in nm.

use crate::cell::{Cell, validate_cell};
use crate::errors::UcmdError;

/// Computes the perpendicular heights of the cell in the a, b, c directions.
/// Returns `Ok([h_a, h_b, h_c])` if the cell is valid, `Err(UcmdError::InvalidCell)` if invalid.
/// Units: nm.
pub fn perpendicular_heights(cell: &Cell) -> Result<[f64; 3], UcmdError> {
    validate_cell(cell)?;

    let ax: f64 = cell.ax;
    let bx: f64 = cell.bx;
    let by: f64 = cell.by;
    let cx: f64 = cell.cx;
    let cy: f64 = cell.cy;
    let cz: f64 = cell.cz;

    let v: f64 = ax * by * cz;

    // |b × c|
    let bcc_x: f64 = by * cz;
    let bcc_y: f64 = -bx * cz;
    let bcc_z: f64 = bx * cy - by * cx;
    let bcc_norm_sq: f64 = bcc_x * bcc_x + bcc_y * bcc_y + bcc_z * bcc_z;
    let bcc_norm: f64 = bcc_norm_sq.sqrt();

    // |a × c|
    let acc_y: f64 = -ax * cz;
    let acc_z: f64 = ax * cy;
    let acc_norm_sq: f64 = acc_y * acc_y + acc_z * acc_z;
    let acc_norm: f64 = acc_norm_sq.sqrt();

    // |a × b| = ax * by
    let ab_norm: f64 = ax * by;

    let h_a: f64 = v / bcc_norm;
    let h_b: f64 = v / acc_norm;
    let h_c: f64 = v / ab_norm;

    Ok([h_a, h_b, h_c])
}

/// Validates that a cutoff radius does not exceed half the smallest perpendicular cell height.
/// Returns `Ok(())` if the cutoff is within the bound, `Err(UcmdError::InvalidCell)` if it exceeds the bound or the cell is invalid.
/// Units: cutoff_nm in nm.
pub fn validate_cutoff(cell: &Cell, cutoff_nm: f64) -> Result<(), UcmdError> {
    let heights = perpendicular_heights(cell)?;
    let h_min: f64 = heights[0].min(heights[1]).min(heights[2]);
    if cutoff_nm >= h_min / 2.0 {
        Err(UcmdError::InvalidCell)
    } else {
        Ok(())
    }
}
