//! Tests of erf and erfc (R4.1: relative error below 1e-13 wherever the value is above 1e-300).
//!
//! Reference values: a 400-digit evaluation of the positive series
//! erf(x) = (2/√π) e^(−x²) Σ 2ⁿ x^(2n+1)/(2n+1)!! in Python `decimal` (π by Machin's formula),
//! erfc = 1 − erf at that precision, rounded to the nearest f64. They were cross-checked against
//! glibc's erfc (Python `math.erfc`) to 2.4e-16 relative.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::erf::{erf, erfc};

/// (x, erf(x), erfc(x)).
const REFERENCE: [(f64, f64, f64); 19] = [
    (0.0, 0.0, 1.0),
    (1e-300, 1.1283791670955126e-300, 1.0),
    (1e-10, 1.1283791670955126e-10, 0.999999999887162),
    (0.001, 0.0011283787909692365, 0.9988716212090307),
    (0.1, 0.1124629160182849, 0.887537083981715),
    (0.5, 0.5204998778130465, 0.4795001221869535),
    (0.9999999999999999, 0.8427007929497148, 0.1572992070502852),
    (1.0, 0.8427007929497149, 0.15729920705028513),
    (1.5, 0.9661051464753108, 0.033894853524689274),
    (2.0, 0.9953222650189527, 0.004677734981047266),
    (3.0, 0.9999779095030014, 2.209049699858544e-05),
    (4.5, 0.9999999998033839, 1.9661604415428876e-10),
    (6.0, 1.0, 2.1519736712498913e-17),
    (10.0, 1.0, 2.088487583762545e-45),
    (15.0, 1.0, 7.212994172451207e-100),
    (20.0, 1.0, 5.395865611607901e-176),
    (26.0, 1.0, 5.663192408856143e-296),
    (-0.5, -0.5204998778130465, 1.5204998778130465),
    (-2.0, -0.9953222650189527, 1.9953222650189528),
];

fn rel_err(got: f64, want: f64) -> f64 {
    if want == 0.0 {
        got.abs()
    } else {
        ((got - want) / want).abs()
    }
}

#[test]
fn erfc_matches_reference_to_1e_13_relative() {
    for (x, _, want) in REFERENCE {
        let got = erfc(x);
        let err = rel_err(got, want);
        assert!(
            err < 1e-13,
            "erfc({x}) = {got:e}, reference {want:e}, relative error {err:e}"
        );
    }
}

#[test]
fn erf_matches_reference_to_1e_13_relative() {
    for (x, want, _) in REFERENCE {
        let got = erf(x);
        let err = rel_err(got, want);
        assert!(
            err < 1e-13,
            "erf({x}) = {got:e}, reference {want:e}, relative error {err:e}"
        );
    }
}

#[test]
fn symmetry_for_negative_arguments() {
    for k in 0..200 {
        let x = 0.037 * k as f64; // arbitrary sample points in [0, 7.4)
        assert_eq!(erf(-x), -erf(x), "erf(-{x})");
        let lhs = erfc(-x);
        let rhs = 2.0 - erfc(x);
        assert!(
            (lhs - rhs).abs() <= 4.0 * f64::EPSILON,
            "erfc(-{x}) = {lhs}, 2 - erfc({x}) = {rhs}"
        );
    }
}

#[test]
fn erf_plus_erfc_is_one() {
    for k in 0..400 {
        let x = -5.0 + 0.025 * k as f64; // arbitrary sample points in [-5, 5)
        let s = erf(x) + erfc(x);
        assert!(
            (s - 1.0).abs() <= 4.0 * f64::EPSILON,
            "erf + erfc at {x} = {s}"
        );
    }
}

#[test]
fn continuous_across_the_method_switch_at_one() {
    // Below 1: series; from 1 on: continued fraction. Neighbouring arguments must agree to the
    // slope times the step (erfc'(1) = -2/sqrt(pi) e^-1 ~ -0.415) plus rounding.
    let below = 1.0_f64 - f64::EPSILON / 2.0;
    let jump = erfc(below) - erfc(1.0);
    let slope = std::f64::consts::FRAC_2_SQRT_PI * (-1.0_f64).exp();
    let slope_step = slope * (1.0 - below);
    assert!(
        (jump - slope_step).abs() < 1e-15,
        "erfc jump at 1: {jump:e}"
    );
}

#[test]
fn limits_and_nan() {
    assert_eq!(erfc(f64::INFINITY), 0.0);
    assert_eq!(erfc(f64::NEG_INFINITY), 2.0);
    assert_eq!(erf(f64::INFINITY), 1.0);
    assert_eq!(erf(f64::NEG_INFINITY), -1.0);
    assert_eq!(erfc(30.0), 0.0); // true value ~ 2.6e-393, below the smallest subnormal
    assert!(erf(f64::NAN).is_nan());
    assert!(erfc(f64::NAN).is_nan());
}

#[test]
fn erfc_is_monotonic_decreasing() {
    let mut prev = erfc(-6.0);
    for k in 1..=3200 {
        let x = -6.0 + 0.01 * k as f64; // arbitrary grid on (-6, 26]
        let v = erfc(x);
        assert!(v <= prev, "erfc not decreasing at {x}: {v:e} > {prev:e}");
        prev = v;
    }
}
