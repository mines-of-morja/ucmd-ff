//! Basic vector arithmetic helpers for `[f64; 3]` arrays.
//!
//! This module provides element-wise and standard vector-algebra operations
//! used throughout the engine for positions, displacements, and forces.
//!
//! Units:
//! - Length: nm
//! - Force: kJ mol^-1 nm^-1
//! - Scalar factors: dimensionless unless otherwise specified
//! - Dot and cross products: units are the product of the input units.
//! - Norm: same units as the input vector.

/// Adds two 3D vectors.
/// Units: depends on inputs (e.g., nm for positions, kJ mol^-1 nm^-1 for forces).
pub fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtracts vector `b` from vector `a`.
/// Units: depends on inputs (e.g., nm for displacements).
pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scales a vector by a scalar factor.
/// Units: output units = input units * scalar units.
pub fn scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Computes the dot product of two vectors.
/// Units: output units = (units of a) * (units of b).
pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Computes the cross product of two vectors.
/// Units: output units = (units of a) * (units of b).
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Computes the Euclidean norm (magnitude) of a vector.
/// Units: same as input units (e.g., nm for displacement).
pub fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Three 3-vectors, one per atom of a three-atom term (forces or gradients).
/// Units: those of the quantity (forces in kJ mol⁻¹ nm⁻¹).
pub type PerAtom3 = ([f64; 3], [f64; 3], [f64; 3]);

/// Four 3-vectors, one per atom of a four-atom term (forces or gradients).
/// Units: those of the quantity (forces in kJ mol⁻¹ nm⁻¹).
pub type PerAtom4 = ([f64; 3], [f64; 3], [f64; 3], [f64; 3]);
