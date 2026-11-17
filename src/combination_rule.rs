//! Combination rules for deriving pair Lennard-Jones parameters from per-type values.
//!
//! Provides the geometric (OPLS) and Lorentz-Berthelot rules for combining
//! per-type sigma (nm) and epsilon (kJ/mol) values into pair parameters.

/// The combination rule used to derive pair parameters from per-type values.
/// Units: none (enumeration of arithmetic rules).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CombinationRule {
    /// Geometric (OPLS) rule: sigma_ij = sqrt(sigma_i * sigma_j).
    Geometric,
    /// Lorentz-Berthelot rule: sigma_ij = (sigma_i + sigma_j) / 2.
    LorentzBerthelot,
}

/// Computes the pair sigma from per-type sigmas using the given rule.
/// Units: sigma_i, sigma_j in nm; return sigma_ij in nm.
///
/// Geometric: sigma_ij = sqrt(sigma_i * sigma_j).
/// Lorentz-Berthelot: sigma_ij = (sigma_i + sigma_j) / 2.
pub fn combine_sigma(rule: CombinationRule, sigma_i: f64, sigma_j: f64) -> f64 {
    match rule {
        CombinationRule::Geometric => (sigma_i * sigma_j).sqrt(),
        CombinationRule::LorentzBerthelot => (sigma_i + sigma_j) / 2.0,
    }
}

/// Computes the pair epsilon from per-type epsilons.
/// Units: eps_i, eps_j in kJ/mol; return eps_ij in kJ/mol.
///
/// Both rules use the geometric mean: eps_ij = sqrt(eps_i * eps_j).
pub fn combine_eps(eps_i: f64, eps_j: f64) -> f64 {
    (eps_i * eps_j).sqrt()
}
