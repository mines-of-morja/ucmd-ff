//! CODATA 2018 physical constants and derived values for UnitcellMD.
//!
//! Base constants are SI values. Derived constants use project units:
//! energy kJ mol^-1, length nm, charge e, temperature K.

/// Elementary charge, in C.
pub const E_CHARGE_C: f64 = 1.602176634e-19;

/// Avogadro constant, in mol^-1.
pub const AVOGADRO_MOL: f64 = 6.02214076e23;

/// Boltzmann constant, in J K^-1.
pub const K_B_J_K: f64 = 1.380649e-23;

/// Vacuum permittivity, in F m^-1.
pub const EPS0_F_M: f64 = 8.8541878128e-12;

/// Molar gas constant, in kJ mol^-1 K^-1.
pub const R_GAS_KJ_MOL_K: f64 = K_B_J_K * AVOGADRO_MOL / 1000.0;

/// Coulomb constant, in kJ mol^-1 nm e^-2.
pub const COULOMB_KJ_MOL_NM_E2: f64 =
    E_CHARGE_C * E_CHARGE_C * AVOGADRO_MOL / (4.0 * std::f64::consts::PI * EPS0_F_M) * 1_000_000.0;
