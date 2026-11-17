//! Converts OPLS-Fourier torsion coefficients to Ryckaert-Bellemans coefficients.
//! Units: input v and output c are both in kJ/mol.

/// Converts OPLS-Fourier torsion coefficients to Ryckaert-Bellemans coefficients.
/// Units: input v and output c are both in kJ/mol.
///
/// # Arguments
/// * `v` — OPLS-Fourier coefficients [V1, V2, V3, V4] in kJ/mol.
///
/// # Returns
/// Ryckaert-Bellemans coefficients [C0, C1, C2, C3, C4, C5] in kJ/mol,
/// where C5 is always 0.
pub fn opls_to_rb(v: [f64; 4]) -> [f64; 6] {
    let v1 = v[0];
    let v2 = v[1];
    let v3 = v[2];
    let v4 = v[3];

    [
        0.5 * (v1 + 2.0 * v2 + v3),
        0.5 * (3.0 * v3 - v1),
        4.0 * v4 - v2,
        -2.0 * v3,
        -4.0 * v4,
        0.0,
    ]
}
