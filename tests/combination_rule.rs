//! Tests of the `combination_rule` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::combination_rule::*;

#[test]
fn identical_types_geometric() {
    let sigma_i: f64 = 0.34;
    let sigma_j: f64 = 0.34;
    let got = combine_sigma(CombinationRule::Geometric, sigma_i, sigma_j);
    let expected: f64 = 0.34;
    let tol: f64 = 1e-15;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

#[test]
fn identical_types_lorentz_berthelot() {
    let sigma_i: f64 = 0.34;
    let sigma_j: f64 = 0.34;
    let got = combine_sigma(CombinationRule::LorentzBerthelot, sigma_i, sigma_j);
    let expected: f64 = 0.34;
    let tol: f64 = 1e-15;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

#[test]
fn geometric_arbitrary_values() {
    let sigma_i: f64 = 0.3;
    let sigma_j: f64 = 0.5;
    let got = combine_sigma(CombinationRule::Geometric, sigma_i, sigma_j);
    let expected = (sigma_i * sigma_j).sqrt();
    let tol: f64 = 1e-15;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

#[test]
fn lorentz_berthelot_arbitrary_values() {
    let sigma_i: f64 = 0.3;
    let sigma_j: f64 = 0.5;
    let got = combine_sigma(CombinationRule::LorentzBerthelot, sigma_i, sigma_j);
    let expected = (sigma_i + sigma_j) / 2.0;
    let tol: f64 = 1e-15;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

#[test]
fn eps_geometric_mean() {
    let eps_i: f64 = 0.6502994552;
    let eps_j: f64 = 0.43932;
    let got = combine_eps(eps_i, eps_j);
    let expected = (eps_i * eps_j).sqrt();
    let tol: f64 = 1e-12;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

#[test]
fn zero_sigma_geometric() {
    let sigma_i: f64 = 0.0;
    let sigma_j: f64 = 0.3165492;
    let got = combine_sigma(CombinationRule::Geometric, sigma_i, sigma_j);
    let expected: f64 = 0.0;
    let tol: f64 = 1e-15;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

#[test]
fn zero_sigma_lorentz_berthelot() {
    let sigma_i: f64 = 0.0;
    let sigma_j: f64 = 0.3165492;
    let got = combine_sigma(CombinationRule::LorentzBerthelot, sigma_i, sigma_j);
    let expected = (sigma_i + sigma_j) / 2.0;
    let tol: f64 = 1e-15;
    assert!(
        (got - expected).abs() <= tol,
        "got {got}, expected {expected}"
    );
}

#[test]
fn symmetry_sigma() {
    let a: f64 = 0.3165492;
    let b: f64 = 0.0;

    // Geometric rule
    let g1 = combine_sigma(CombinationRule::Geometric, a, b);
    let g2 = combine_sigma(CombinationRule::Geometric, b, a);
    let tol: f64 = 1e-15;
    assert!((g1 - g2).abs() <= tol, "geometric: got {g1}, expected {g2}");

    // Lorentz-Berthelot rule
    let l1 = combine_sigma(CombinationRule::LorentzBerthelot, a, b);
    let l2 = combine_sigma(CombinationRule::LorentzBerthelot, b, a);
    assert!(
        (l1 - l2).abs() <= tol,
        "lorentz-berthelot: got {l1}, expected {l2}"
    );
}
