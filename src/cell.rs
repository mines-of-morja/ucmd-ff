//! Periodic cell validation, volume calculation, and minimum-image displacement.
//! Units: lengths in nm, volume in nm^3, displacements in nm.

use crate::errors::UcmdError;

/// A periodic cell in reduced lower-triangular form.
/// Units: all fields are length in nm.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cell {
    /// Lattice vector a x-component. Units: nm.
    pub ax: f64,
    /// Lattice vector b x-component. Units: nm.
    pub bx: f64,
    /// Lattice vector b y-component. Units: nm.
    pub by: f64,
    /// Lattice vector c x-component. Units: nm.
    pub cx: f64,
    /// Lattice vector c y-component. Units: nm.
    pub cy: f64,
    /// Lattice vector c z-component. Units: nm.
    pub cz: f64,
}

/// Validates a periodic cell.
/// Returns `Ok(())` if valid, `Err(UcmdError::InvalidCell)` if invalid.
/// Units: none.
pub fn validate_cell(cell: &Cell) -> Result<(), UcmdError> {
    let coords = [cell.ax, cell.bx, cell.by, cell.cx, cell.cy, cell.cz];
    for c in coords {
        if !c.is_finite() {
            return Err(UcmdError::NonFiniteCoordinate);
        }
    }
    if cell.ax > 0.0 && cell.by > 0.0 && cell.cz > 0.0 {
        Ok(())
    } else {
        Err(UcmdError::InvalidCell)
    }
}

/// Computes the volume of a periodic cell.
/// Returns `Ok(volume)` if the cell is valid, `Err(UcmdError::InvalidCell)` or `Err(UcmdError::NonFiniteCoordinate)` if invalid.
/// Units: nm^3.
pub fn cell_volume(cell: &Cell) -> Result<f64, UcmdError> {
    validate_cell(cell)?;
    Ok(cell.ax * cell.by * cell.cz)
}

/// Computes the minimum-image displacement from point `r_j` to point `r_i`.
/// Returns the shortest vector `d` such that `r_i - r_j - d` is a lattice vector.
/// Returns `Err(UcmdError::InvalidCell)` or `Err(UcmdError::NonFiniteCoordinate)` if the cell or coordinates are invalid.
/// Units: nm.
pub fn minimum_image(r_i: [f64; 3], r_j: [f64; 3], cell: &Cell) -> Result<[f64; 3], UcmdError> {
    for k in 0..3 {
        if !r_i[k].is_finite() || !r_j[k].is_finite() {
            return Err(UcmdError::NonFiniteCoordinate);
        }
    }
    validate_cell(cell)?;

    let v = [r_i[0] - r_j[0], r_i[1] - r_j[1], r_i[2] - r_j[2]];
    let a = [cell.ax, 0.0, 0.0];
    let b = [cell.bx, cell.by, 0.0];
    let c = [cell.cx, cell.cy, cell.cz];

    // Reduce displacement into the central cell
    let n_c = (v[2] / cell.cz).round();
    if !n_c.is_finite() {
        return Err(UcmdError::NonFiniteCoordinate);
    }
    let v1 = [v[0] - n_c * c[0], v[1] - n_c * c[1], v[2] - n_c * c[2]];

    let n_b = (v1[1] / cell.by).round();
    if !n_b.is_finite() {
        return Err(UcmdError::NonFiniteCoordinate);
    }
    let v2 = [v1[0] - n_b * b[0], v1[1] - n_b * b[1], v1[2] - n_b * b[2]];

    let n_a = (v2[0] / cell.ax).round();
    if !n_a.is_finite() {
        return Err(UcmdError::NonFiniteCoordinate);
    }
    let v3 = [v2[0] - n_a * a[0], v2[1] - n_a * a[1], v2[2] - n_a * a[2]];

    if !v3[0].is_finite() || !v3[1].is_finite() || !v3[2].is_finite() {
        return Err(UcmdError::NonFiniteCoordinate);
    }

    let mut best = [f64::NAN; 3];
    let mut best_norm_sq = f64::INFINITY;

    for m_a in -2..=2 {
        for m_b in -2..=2 {
            for m_c in -2..=2 {
                let ma = m_a as f64;
                let mb = m_b as f64;
                let mc = m_c as f64;

                let l = [
                    ma * a[0] + mb * b[0] + mc * c[0],
                    ma * a[1] + mb * b[1] + mc * c[1],
                    ma * a[2] + mb * b[2] + mc * c[2],
                ];
                let d = [v3[0] - l[0], v3[1] - l[1], v3[2] - l[2]];

                if d[0].is_finite() && d[1].is_finite() && d[2].is_finite() {
                    let norm_sq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                    if norm_sq.is_finite() && norm_sq < best_norm_sq {
                        best_norm_sq = norm_sq;
                        best = d;
                    }
                }
            }
        }
    }

    if !best_norm_sq.is_finite() {
        return Err(UcmdError::NonFiniteCoordinate);
    }

    Ok(best)
}
