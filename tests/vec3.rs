//! Tests of the `vec3` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::vec3::*;

#[test]
fn test_norm_known_value() {
    let v = [3.0, 4.0, 0.0];
    let got = norm(v);
    let expected: f64 = f64::sqrt(3.0 * 3.0 + 4.0 * 4.0 + 0.0 * 0.0);
    let tol = 1e-12;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

#[test]
fn test_dot_orthogonal() {
    let a = [1.0, 0.0, 0.0];
    let b = [0.0, 1.0, 0.0];
    let got = dot(a, b);
    let expected = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let tol = 1e-12;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

#[test]
fn test_cross_magnitude() {
    let a = [1.0, 0.0, 0.0];
    let b = [0.0, 1.0, 0.0];
    let got = norm(cross(a, b));

    let cx = a[1] * b[2] - a[2] * b[1];
    let cy = a[2] * b[0] - a[0] * b[2];
    let cz = a[0] * b[1] - a[1] * b[0];
    let expected = (cx * cx + cy * cy + cz * cz).sqrt();

    let tol = 1e-12;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

#[test]
fn test_add_sub_inverse() {
    let a = [1.5, -2.0, 3.0];
    let b = [0.5, 1.0, -1.0];
    let got = sub(add(a, b), b);
    let expected = a;
    let tol = 1e-12;

    assert!(
        (got[0] - expected[0]).abs() <= tol,
        "got {got:?}, expected {expected:?}"
    );
    assert!(
        (got[1] - expected[1]).abs() <= tol,
        "got {got:?}, expected {expected:?}"
    );
    assert!(
        (got[2] - expected[2]).abs() <= tol,
        "got {got:?}, expected {expected:?}"
    );
}

#[test]
fn test_scale_distributive() {
    let a = [1.0, 2.0, 3.0];
    let b = [4.0, 5.0, 6.0];
    let s = 2.0;

    let got_lhs = scale(add(a, b), s);
    let got_rhs = add(scale(a, s), scale(b, s));

    let expected = [(a[0] + b[0]) * s, (a[1] + b[1]) * s, (a[2] + b[2]) * s];
    let tol = 1e-12;

    assert!(
        (got_lhs[0] - expected[0]).abs() <= tol,
        "got {got_lhs:?}, expected {expected:?}"
    );
    assert!(
        (got_lhs[1] - expected[1]).abs() <= tol,
        "got {got_lhs:?}, expected {expected:?}"
    );
    assert!(
        (got_lhs[2] - expected[2]).abs() <= tol,
        "got {got_lhs:?}, expected {expected:?}"
    );

    assert!(
        (got_rhs[0] - expected[0]).abs() <= tol,
        "got {got_rhs:?}, expected {expected:?}"
    );
    assert!(
        (got_rhs[1] - expected[1]).abs() <= tol,
        "got {got_rhs:?}, expected {expected:?}"
    );
    assert!(
        (got_rhs[2] - expected[2]).abs() <= tol,
        "got {got_rhs:?}, expected {expected:?}"
    );
}

#[test]
fn test_cross_orthogonality() {
    let a = [1.0, 2.0, 3.0];
    let b = [4.0, 5.0, 6.0];
    let c = cross(a, b);

    let got_a = dot(c, a);
    let got_b = dot(c, b);
    let expected = 0.0;
    let tol = 1e-12;

    assert!(
        (got_a - expected).abs() <= tol,
        "got {got_a}, expected {expected}"
    );
    assert!(
        (got_b - expected).abs() <= tol,
        "got {got_b}, expected {expected}"
    );
}

#[test]
fn test_norm_non_negative() {
    let v = [-1.0, 2.0, -3.0];
    let got = norm(v);
    let expected = 0.0;
    assert!(got >= expected, "got {got}, expected {expected}");
}
