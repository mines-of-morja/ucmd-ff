//! Harmonic angle energy, forces, and virial contribution.
//!
//! The central atom is `r2`; terminals are `r1` and `r3`.
//! The angle $\theta$ is the angle between vectors $\mathbf{u} = \mathbf{r}_1 - \mathbf{r}_2$
//! and $\mathbf{v} = \mathbf{r}_3 - \mathbf{r}_2$.
//!
//! Units:
//! - Positions: nm
//! - Angle: rad
//! - Force constant: kJ mol^-1 rad^-2
//! - Energy: kJ/mol
//! - Forces: kJ mol^-1 nm^-1
//! - Virial: kJ/mol

use crate::errors::UcmdError;
use crate::vec3;
use crate::vec3::PerAtom3;

/// Checks if a vector contains any NaN or infinite coordinates.
/// Returns `Err(UcmdError::NonFiniteCoordinate)` if so.
fn check_finite(v: [f64; 3]) -> Result<(), UcmdError> {
    for &x in v.iter() {
        if !x.is_finite() {
            return Err(UcmdError::NonFiniteCoordinate);
        }
    }
    Ok(())
}

/// Computes the energy of a harmonic angle.
/// Units: r1, r2, r3 in nm; theta0 in rad; k_kjmol_rad2 in kJ mol^-1 rad^-2; return in kJ/mol.
pub fn angle_energy(
    r1: [f64; 3],
    r2: [f64; 3],
    r3: [f64; 3],
    theta0: f64,
    k_kjmol_rad2: f64,
) -> Result<f64, UcmdError> {
    check_finite(r1)?;
    check_finite(r2)?;
    check_finite(r3)?;

    let u = vec3::sub(r1, r2);
    let v = vec3::sub(r3, r2);

    let a = vec3::norm(u);
    let b = vec3::norm(v);

    if !a.is_finite() || !b.is_finite() {
        return Err(UcmdError::NonFiniteCoordinate);
    }

    if a == 0.0 || b == 0.0 {
        return Err(UcmdError::CoincidentAtoms);
    }

    let dot_uv = vec3::dot(u, v);
    let c = dot_uv / (a * b);

    // Clamp c to [-1, 1] to handle floating point errors
    let c_clamped = c.clamp(-1.0, 1.0);
    let theta = c_clamped.acos();

    let s_sq: f64 = 1.0 - c_clamped * c_clamped;
    let s = s_sq.sqrt();
    if s == 0.0 {
        return Err(UcmdError::CoincidentAtoms);
    }

    let delta_theta = theta - theta0;
    let energy = 0.5 * k_kjmol_rad2 * delta_theta * delta_theta;

    Ok(energy)
}

/// Computes the forces on atoms 1, 2, and 3 for a harmonic angle.
/// Units: r1, r2, r3 in nm; theta0 in rad; k_kjmol_rad2 in kJ mol^-1 rad^-2; return forces in kJ mol^-1 nm^-1.
pub fn angle_forces(
    r1: [f64; 3],
    r2: [f64; 3],
    r3: [f64; 3],
    theta0: f64,
    k_kjmol_rad2: f64,
) -> Result<PerAtom3, UcmdError> {
    check_finite(r1)?;
    check_finite(r2)?;
    check_finite(r3)?;

    let u = vec3::sub(r1, r2);
    let v = vec3::sub(r3, r2);

    let a = vec3::norm(u);
    let b = vec3::norm(v);

    if !a.is_finite() || !b.is_finite() {
        return Err(UcmdError::NonFiniteCoordinate);
    }

    if a == 0.0 || b == 0.0 {
        return Err(UcmdError::CoincidentAtoms);
    }

    let dot_uv = vec3::dot(u, v);
    let c = dot_uv / (a * b);

    // Clamp c to [-1, 1] to handle floating point errors
    let c_clamped = c.clamp(-1.0, 1.0);
    let theta = c_clamped.acos();

    let s_sq: f64 = 1.0 - c_clamped * c_clamped;
    let s = s_sq.sqrt();
    if s == 0.0 {
        return Err(UcmdError::CoincidentAtoms);
    }

    let delta_theta = theta - theta0;
    let factor = k_kjmol_rad2 * delta_theta / (a * b * s);

    // F1 = factor * (v - c * (b/a) * u)
    let u_hat = vec3::scale(u, 1.0 / a);
    let term1 = vec3::sub(v, vec3::scale(u_hat, c_clamped * b));
    let f1 = vec3::scale(term1, factor);

    // F3 = factor * (u - c * (a/b) * v)
    let v_hat = vec3::scale(v, 1.0 / b);
    let term3 = vec3::sub(u, vec3::scale(v_hat, c_clamped * a));
    let f3 = vec3::scale(term3, factor);

    // F2 = -(F1 + F3)
    let f1_f3_sum = vec3::add(f1, f3);
    let f2 = vec3::scale(f1_f3_sum, -1.0);

    Ok((f1, f2, f3))
}

/// Computes the virial contribution of a harmonic angle.
/// Units: r1, r2, r3 in nm; theta0 in rad; k_kjmol_rad2 in kJ mol^-1 rad^-2; return in kJ/mol.
pub fn angle_virial(
    r1: [f64; 3],
    r2: [f64; 3],
    r3: [f64; 3],
    theta0: f64,
    k_kjmol_rad2: f64,
) -> Result<[[f64; 3]; 3], UcmdError> {
    check_finite(r1)?;
    check_finite(r2)?;
    check_finite(r3)?;

    let u = vec3::sub(r1, r2);
    let v = vec3::sub(r3, r2);

    let a = vec3::norm(u);
    let b = vec3::norm(v);

    if !a.is_finite() || !b.is_finite() {
        return Err(UcmdError::NonFiniteCoordinate);
    }

    if a == 0.0 || b == 0.0 {
        return Err(UcmdError::CoincidentAtoms);
    }

    let dot_uv = vec3::dot(u, v);
    let c = dot_uv / (a * b);

    // Clamp c to [-1, 1] to handle floating point errors
    let c_clamped = c.clamp(-1.0, 1.0);
    let theta = c_clamped.acos();

    let s_sq: f64 = 1.0 - c_clamped * c_clamped;
    let s = s_sq.sqrt();
    if s == 0.0 {
        return Err(UcmdError::CoincidentAtoms);
    }

    let delta_theta = theta - theta0;
    let factor = k_kjmol_rad2 * delta_theta / (a * b * s);

    // F1 = factor * (v - c * (b/a) * u)
    let u_hat = vec3::scale(u, 1.0 / a);
    let term1 = vec3::sub(v, vec3::scale(u_hat, c_clamped * b));
    let f1 = vec3::scale(term1, factor);

    // F3 = factor * (u - c * (a/b) * v)
    let v_hat = vec3::scale(v, 1.0 / b);
    let term3 = vec3::sub(u, vec3::scale(v_hat, c_clamped * a));
    let f3 = vec3::scale(term3, factor);

    // Virial W_ab = F1_a * u_b + F3_a * v_b
    // This matches the finite difference definition W_ab = sum_i F_i_a * r_i_b
    // where r12 = u and r32 = v, and F2 = -F1 - F3.
    let mut w: [[f64; 3]; 3] = [[0.0; 3]; 3];

    for a_idx in 0..3 {
        for b_idx in 0..3 {
            w[a_idx][b_idx] = f1[a_idx] * u[b_idx] + f3[a_idx] * v[b_idx];
        }
    }

    Ok(w)
}
