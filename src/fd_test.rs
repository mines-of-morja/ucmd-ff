//! Finite-difference test utilities for validating analytic forces and virials
//! against numerical derivatives of energy.
//!
//! Units: positions in nm, energy in kJ/mol, forces in kJ mol^-1 nm^-1,
//! virial in kJ/mol, strain dimensionless.

use crate::cell::Cell;
use crate::errors::UcmdError;

/// Energy of a harmonic bond as a function of (r_i, r_j, r0, k), as `bond::bond_energy`.
/// Units: nm, nm, nm, kJ mol⁻¹ nm⁻²; returns kJ/mol.
pub type BondEnergyFn = fn([f64; 3], [f64; 3], f64, f64) -> Result<f64, UcmdError>;

/// Energy of a harmonic angle as a function of (r1, r2, r3, θ0, k), as `angle::angle_energy`.
/// Units: nm, nm, nm, rad, kJ mol⁻¹ rad⁻²; returns kJ/mol.
pub type AngleEnergyFn = fn([f64; 3], [f64; 3], [f64; 3], f64, f64) -> Result<f64, UcmdError>;

/// Energy of all positions. Units: positions in nm; returns kJ/mol.
pub type EnergyOfPositions<'a> = dyn Fn(&[[f64; 3]]) -> Result<f64, UcmdError> + 'a;

/// Energy of all positions and the periodic cell. Units: nm; returns kJ/mol.
pub type EnergyOfPositionsAndCell<'a> = dyn Fn(&[[f64; 3]], &Cell) -> Result<f64, UcmdError> + 'a;

/// Returns true if the (a, b) strain component keeps the cell lower-triangular.
/// Units: none (categorical).
fn is_valid_strain(a: usize, b: usize) -> bool {
    a >= b
}

/// Computes the numerical force component for a bond using central differences.
/// Units: r_i, r_j in nm; r0_nm in nm; k_kjmol_nm2 in kJ mol^-1 nm^-2; h in nm; return in kJ mol^-1 nm^-1.
#[allow(clippy::too_many_arguments)] // the term's positions and parameters, atom and step
pub fn numerical_force_bond(
    energy_fn: BondEnergyFn,
    r_i: [f64; 3],
    r_j: [f64; 3],
    r0_nm: f64,
    k_kjmol_nm2: f64,
    atom_idx: usize,
    coord_idx: usize,
    h: f64,
) -> Result<f64, UcmdError> {
    if atom_idx > 1 {
        return Err(UcmdError::IndexOutOfRange);
    }
    if coord_idx > 2 {
        return Err(UcmdError::IndexOutOfRange);
    }
    if !h.is_finite() || h <= 0.0 {
        return Err(UcmdError::NonFiniteCoordinate);
    }

    let (ri_p, rj_p) = if atom_idx == 0 {
        (perturb_coord(r_i, coord_idx, h), r_j)
    } else {
        (r_i, perturb_coord(r_j, coord_idx, h))
    };
    let (ri_m, rj_m) = if atom_idx == 0 {
        (perturb_coord(r_i, coord_idx, -h), r_j)
    } else {
        (r_i, perturb_coord(r_j, coord_idx, -h))
    };

    let e_plus = energy_fn(ri_p, rj_p, r0_nm, k_kjmol_nm2)?;
    let e_minus = energy_fn(ri_m, rj_m, r0_nm, k_kjmol_nm2)?;

    Ok(-(e_plus - e_minus) / (2.0 * h))
}

/// Computes the numerical force component for an angle using central differences.
/// Units: r1, r2, r3 in nm; theta0 in rad; k_kjmol_rad2 in kJ mol^-1 rad^-2; h in nm; return in kJ mol^-1 nm^-1.
#[allow(clippy::too_many_arguments)] // the term's positions and parameters, atom and step
pub fn numerical_force_angle(
    energy_fn: AngleEnergyFn,
    r1: [f64; 3],
    r2: [f64; 3],
    r3: [f64; 3],
    theta0: f64,
    k_kjmol_rad2: f64,
    atom_idx: usize,
    coord_idx: usize,
    h: f64,
) -> Result<f64, UcmdError> {
    if atom_idx > 2 {
        return Err(UcmdError::IndexOutOfRange);
    }
    if coord_idx > 2 {
        return Err(UcmdError::IndexOutOfRange);
    }
    if !h.is_finite() || h <= 0.0 {
        return Err(UcmdError::NonFiniteCoordinate);
    }

    let (a1_p, a2_p, a3_p) = match atom_idx {
        0 => (perturb_coord(r1, coord_idx, h), r2, r3),
        1 => (r1, perturb_coord(r2, coord_idx, h), r3),
        2 => (r1, r2, perturb_coord(r3, coord_idx, h)),
        _ => return Err(UcmdError::IndexOutOfRange),
    };
    let (a1_m, a2_m, a3_m) = match atom_idx {
        0 => (perturb_coord(r1, coord_idx, -h), r2, r3),
        1 => (r1, perturb_coord(r2, coord_idx, -h), r3),
        2 => (r1, r2, perturb_coord(r3, coord_idx, -h)),
        _ => return Err(UcmdError::IndexOutOfRange),
    };

    let e_plus = energy_fn(a1_p, a2_p, a3_p, theta0, k_kjmol_rad2)?;
    let e_minus = energy_fn(a1_m, a2_m, a3_m, theta0, k_kjmol_rad2)?;

    Ok(-(e_plus - e_minus) / (2.0 * h))
}

/// Computes the numerical virial component for a bond using central differences under strain.
/// Units: r_i, r_j in nm; r0_nm in nm; k_kjmol_nm2 in kJ mol^-1 nm^-2; cell in nm; h in strain (dimensionless); return in kJ/mol.
#[allow(clippy::too_many_arguments)] // the term's positions and parameters, atom and step
pub fn numerical_virial_bond(
    energy_fn: BondEnergyFn,
    r_i: [f64; 3],
    r_j: [f64; 3],
    r0_nm: f64,
    k_kjmol_nm2: f64,
    cell: &Cell,
    strain_row: usize,
    strain_col: usize,
    h: f64,
) -> Result<f64, UcmdError> {
    if strain_row > 2 {
        return Err(UcmdError::IndexOutOfRange);
    }
    if strain_col > 2 {
        return Err(UcmdError::IndexOutOfRange);
    }
    if !h.is_finite() || h <= 0.0 {
        return Err(UcmdError::NonFiniteCoordinate);
    }

    let _ = cell;

    let ri_p = strain_deform(r_i, strain_row, strain_col, h);
    let rj_p = strain_deform(r_j, strain_row, strain_col, h);
    let ri_m = strain_deform(r_i, strain_row, strain_col, -h);
    let rj_m = strain_deform(r_j, strain_row, strain_col, -h);

    let e_plus = energy_fn(ri_p, rj_p, r0_nm, k_kjmol_nm2)?;
    let e_minus = energy_fn(ri_m, rj_m, r0_nm, k_kjmol_nm2)?;

    Ok(-(e_plus - e_minus) / (2.0 * h))
}

/// Computes the numerical virial component for an angle using central differences under strain.
/// Units: r1, r2, r3 in nm; theta0 in rad; k_kjmol_rad2 in kJ mol^-1 rad^-2; cell in nm; h in strain (dimensionless); return in kJ/mol.
#[allow(clippy::too_many_arguments)] // the term's positions and parameters, atom and step
pub fn numerical_virial_angle(
    energy_fn: AngleEnergyFn,
    r1: [f64; 3],
    r2: [f64; 3],
    r3: [f64; 3],
    theta0: f64,
    k_kjmol_rad2: f64,
    cell: &Cell,
    strain_row: usize,
    strain_col: usize,
    h: f64,
) -> Result<f64, UcmdError> {
    if strain_row > 2 {
        return Err(UcmdError::IndexOutOfRange);
    }
    if strain_col > 2 {
        return Err(UcmdError::IndexOutOfRange);
    }
    if !h.is_finite() || h <= 0.0 {
        return Err(UcmdError::NonFiniteCoordinate);
    }

    let _ = cell;

    let a1_p = strain_deform(r1, strain_row, strain_col, h);
    let a2_p = strain_deform(r2, strain_row, strain_col, h);
    let a3_p = strain_deform(r3, strain_row, strain_col, h);
    let a1_m = strain_deform(r1, strain_row, strain_col, -h);
    let a2_m = strain_deform(r2, strain_row, strain_col, -h);
    let a3_m = strain_deform(r3, strain_row, strain_col, -h);

    let e_plus = energy_fn(a1_p, a2_p, a3_p, theta0, k_kjmol_rad2)?;
    let e_minus = energy_fn(a1_m, a2_m, a3_m, theta0, k_kjmol_rad2)?;

    Ok(-(e_plus - e_minus) / (2.0 * h))
}

/// Compares analytic and numerical force vectors within a tolerance.
/// Units: analytic, numerical in kJ mol^-1 nm^-1; tol in kJ mol^-1 nm^-1.
pub fn compare_force(analytic: [f64; 3], numerical: [f64; 3], tol: f64) -> bool {
    for i in 0..3 {
        let a = analytic[i];
        let n = numerical[i];
        if !a.is_finite() || !n.is_finite() {
            return false;
        }
        if (a - n).abs() > tol {
            return false;
        }
    }
    true
}

/// Compares analytic and numerical virial matrices within a tolerance.
/// Units: analytic, numerical in kJ/mol; tol in kJ/mol.
pub fn compare_virial(analytic: [[f64; 3]; 3], numerical: [[f64; 3]; 3], tol: f64) -> bool {
    for i in 0..3 {
        for j in 0..3 {
            let a = analytic[i][j];
            let n = numerical[i][j];
            if !a.is_finite() || !n.is_finite() {
                return false;
            }
            if (a - n).abs() > tol {
                return false;
            }
        }
    }
    true
}

/// Computes numerical forces for all atoms and coordinates using central differences.
/// Units: positions in nm; h_nm in nm; return forces in kJ mol^-1 nm^-1.
pub fn numerical_forces(
    energy: &EnergyOfPositions<'_>,
    positions: &[[f64; 3]],
    h_nm: f64,
) -> Result<Vec<[f64; 3]>, UcmdError> {
    if !h_nm.is_finite() || h_nm <= 0.0 {
        return Err(UcmdError::NonFiniteCoordinate);
    }

    let n = positions.len();
    let mut forces = Vec::with_capacity(n);
    for _ in 0..n {
        forces.push([0.0_f64; 3]);
    }

    for i in 0..n {
        for c in 0..3 {
            let mut plus = positions.to_vec();
            plus[i][c] += h_nm;
            let mut minus = positions.to_vec();
            minus[i][c] -= h_nm;

            let e_plus = energy(&plus)?;
            let e_minus = energy(&minus)?;

            forces[i][c] = -(e_plus - e_minus) / (2.0_f64 * h_nm);
        }
    }

    Ok(forces)
}

/// Computes the numerical virial matrix using central differences under strain.
/// Units: positions in nm; eps dimensionless; return virial in kJ/mol.
pub fn numerical_virial(
    energy: &EnergyOfPositions<'_>,
    positions: &[[f64; 3]],
    eps: f64,
) -> Result<[[f64; 3]; 3], UcmdError> {
    if !eps.is_finite() || eps <= 0.0 {
        return Err(UcmdError::NonFiniteCoordinate);
    }

    let mut virial = [[0.0_f64; 3]; 3];

    for a in 0..3 {
        for b in 0..3 {
            let mut plus = positions.to_vec();
            let mut minus = positions.to_vec();

            for i in 0..plus.len() {
                plus[i][b] += eps * plus[i][a];
                minus[i][b] -= eps * minus[i][a];
            }

            let e_plus = energy(&plus)?;
            let e_minus = energy(&minus)?;

            virial[a][b] = -(e_plus - e_minus) / (2.0_f64 * eps);
        }
    }

    Ok(virial)
}

/// Computes the maximum absolute difference between analytic and numerical force vectors.
/// Units: analytic, numerical in kJ mol^-1 nm^-1; return in kJ mol^-1 nm^-1.
pub fn max_abs_difference_forces(analytic: &[f64; 3], numerical: &[f64; 3]) -> f64 {
    let mut max_diff = 0.0_f64;
    for c in 0..3 {
        let a = analytic[c];
        let n = numerical[c];
        if !a.is_finite() || !n.is_finite() {
            return f64::INFINITY;
        }
        let d = (a - n).abs();
        if d > max_diff {
            max_diff = d;
        }
    }
    max_diff
}

/// Computes the maximum absolute difference between analytic and numerical virial matrices.
/// Units: analytic, numerical in kJ/mol; return in kJ/mol.
pub fn max_abs_difference_virial(analytic: &[[f64; 3]; 3], numerical: &[[f64; 3]; 3]) -> f64 {
    let mut max_diff = 0.0_f64;
    for i in 0..3 {
        for j in 0..3 {
            let a = analytic[i][j];
            let n = numerical[i][j];
            if !a.is_finite() || !n.is_finite() {
                return f64::INFINITY;
            }
            let d = (a - n).abs();
            if d > max_diff {
                max_diff = d;
            }
        }
    }
    max_diff
}

/// Perturbs a coordinate of a 3D vector by a delta.
/// Units: v in nm; delta in nm; return in nm.
fn perturb_coord(v: [f64; 3], coord_idx: usize, delta: f64) -> [f64; 3] {
    let mut result = v;
    result[coord_idx] += delta;
    result
}

/// Applies a strain deformation to a position vector.
/// Units: r in nm; eps dimensionless; return in nm.
fn strain_deform(r: [f64; 3], strain_row: usize, strain_col: usize, eps: f64) -> [f64; 3] {
    let mut result = r;
    result[strain_col] += eps * r[strain_row];
    result
}

/// Applies a homogeneous strain to a periodic cell, deforming each lattice vector
/// so that the b-th component is increased by eps times the a-th component.
/// Units: cell in nm; eps dimensionless; return in nm.
fn strain_cell(cell: &Cell, a: usize, b: usize, eps: f64) -> Cell {
    let mut c = *cell;
    let mut rows: [[f64; 3]; 3] = [[c.ax, 0.0, 0.0], [c.bx, c.by, 0.0], [c.cx, c.cy, c.cz]];
    for row in rows.iter_mut() {
        row[b] += eps * row[a];
    }
    c.ax = rows[0][0];
    c.bx = rows[1][0];
    c.by = rows[1][1];
    c.cx = rows[2][0];
    c.cy = rows[2][1];
    c.cz = rows[2][2];
    c
}

/// Computes the numerical virial component W_ab for a periodic energy term using
/// central differences under homogeneous strain of both positions and cell.
/// Units: positions in nm; cell in nm; eps dimensionless; return in kJ/mol.
///
/// # Errors
/// * `IndexOutOfRange` — a or b is not in 0..=2.
/// * `InvalidCell` — the (a, b) component would break the lower-triangular cell form.
/// * `NonFiniteCoordinate` — eps is not finite or is not positive.
/// * Propagates errors from the energy function.
pub fn numerical_virial_periodic(
    energy: &EnergyOfPositionsAndCell<'_>,
    positions: &[[f64; 3]],
    cell: &Cell,
    a: usize,
    b: usize,
    eps: f64,
) -> Result<f64, UcmdError> {
    if a > 2 || b > 2 {
        return Err(UcmdError::IndexOutOfRange);
    }
    if !is_valid_strain(a, b) {
        return Err(UcmdError::InvalidCell);
    }
    if !eps.is_finite() || eps <= 0.0 {
        return Err(UcmdError::NonFiniteCoordinate);
    }

    let cell_plus = strain_cell(cell, a, b, eps);
    let cell_minus = strain_cell(cell, a, b, -eps);

    let mut pos_plus: Vec<[f64; 3]> = Vec::with_capacity(positions.len());
    let mut pos_minus: Vec<[f64; 3]> = Vec::with_capacity(positions.len());
    for r in positions.iter() {
        pos_plus.push(strain_deform(*r, a, b, eps));
        pos_minus.push(strain_deform(*r, a, b, -eps));
    }

    let e_plus = energy(&pos_plus, &cell_plus)?;
    let e_minus = energy(&pos_minus, &cell_minus)?;

    Ok(-(e_plus - e_minus) / (2.0_f64 * eps))
}
