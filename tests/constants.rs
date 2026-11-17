//! Tests of the `constants` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::constants::*;

#[test]
fn gas_constant_matches_tabulated() {
    let expected = 0.00831446261815324;
    let rel = (R_GAS_KJ_MOL_K - expected).abs() / expected;
    assert!(rel <= 1e-12, "got {R_GAS_KJ_MOL_K}, expected {expected}");
}

#[test]
fn coulomb_constant_matches_tabulated() {
    let expected = 138.93545764438196;
    let rel = (COULOMB_KJ_MOL_NM_E2 - expected).abs() / expected;
    assert!(
        rel <= 1e-12,
        "got {COULOMB_KJ_MOL_NM_E2}, expected {expected}"
    );
}

#[test]
fn gas_constant_formula() {
    let expected = K_B_J_K * AVOGADRO_MOL / 1000.0;
    let rel = (R_GAS_KJ_MOL_K - expected).abs() / expected;
    assert!(rel <= 1e-15, "got {R_GAS_KJ_MOL_K}, expected {expected}");
}

#[test]
fn coulomb_constant_formula() {
    let expected = E_CHARGE_C * E_CHARGE_C * AVOGADRO_MOL / (4.0 * std::f64::consts::PI * EPS0_F_M)
        * 1_000_000.0;
    let rel = (COULOMB_KJ_MOL_NM_E2 - expected).abs() / expected;
    assert!(
        rel <= 1e-15,
        "got {COULOMB_KJ_MOL_NM_E2}, expected {expected}"
    );
}

#[test]
fn base_constants_match_codata() {
    assert_eq!(E_CHARGE_C, 1.602176634e-19);
    assert_eq!(AVOGADRO_MOL, 6.02214076e23);
    assert_eq!(K_B_J_K, 1.380649e-23);
    assert_eq!(EPS0_F_M, 8.8541878128e-12);
}

#[test]
fn coulomb_constant_unit_conversion() {
    let si_value = E_CHARGE_C * E_CHARGE_C * AVOGADRO_MOL / (4.0 * std::f64::consts::PI * EPS0_F_M);
    let rel = (COULOMB_KJ_MOL_NM_E2 - si_value * 1_000_000.0).abs() / si_value;
    assert!(
        rel <= 1e-15,
        "got {COULOMB_KJ_MOL_NM_E2}, expected {}",
        si_value * 1_000_000.0
    );
}
