//! Pure scalar math for the Ryckaert-Bellemans (RB) torsion potential.
//!
//! Given the IUPAC dihedral angle `phi` (rad) and six coefficients `c` (kJ/mol),
//! evaluates the energy (kJ/mol) and its derivative with respect to `phi`
//! (kJ mol^-1 rad^-1). No geometry, no positions, no error type — the functions
//! are polynomials in `cos(phi)` and are well-defined for every real `phi`.
//! The caller is responsible for computing `phi` from positions and for
//! checking geometric validity.

/// Computes the Ryckaert-Bellemans torsion energy.
/// Units: phi in rad; c in kJ/mol; return in kJ/mol.
/// E = sum_{n=0..5} c[n] * (-cos(phi))^n
pub fn rb_energy(phi: f64, c: [f64; 6]) -> f64 {
    let u = -phi.cos();
    // Horner's method for c[0] + c[1]*u + c[2]*u^2 + c[3]*u^3 + c[4]*u^4 + c[5]*u^5
    let result = c[5] * u + c[4];
    let result = result * u + c[3];
    let result = result * u + c[2];
    let result = result * u + c[1];
    result * u + c[0]
}

/// Computes dE/dphi for the Ryckaert-Bellemans torsion energy.
/// Units: phi in rad; c in kJ/mol; return in kJ/mol (per rad, i.e. kJ mol^-1 rad^-1).
/// dE/dphi = sin(phi) * sum_{n=1..5} n * c[n] * (-cos(phi))^(n-1)
pub fn rb_denergy_dphi(phi: f64, c: [f64; 6]) -> f64 {
    let u = -phi.cos();
    // dE/du = c[1] + 2*c[2]*u + 3*c[3]*u^2 + 4*c[4]*u^3 + 5*c[5]*u^4
    // Horner's method:
    let mut result = 5.0_f64 * c[5];
    result = result * u + 4.0_f64 * c[4];
    result = result * u + 3.0_f64 * c[3];
    result = result * u + 2.0_f64 * c[2];
    result = result * u + c[1];
    // dE/dphi = sin(phi) * dE/du
    phi.sin() * result
}
