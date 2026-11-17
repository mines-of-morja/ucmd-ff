//! ucmd-ff: a classical force-field engine for UnitcellMD.
//!
//! Given a parameterized system (`system::System`), positions and an optional periodic cell,
//! `forcefield::ForceField::evaluate` returns the potential energy, the forces F = −∂E/∂r and
//! the virial W = −∂E/∂ε, with the energy split into named terms. Functional forms and
//! conventions follow the UnitcellMD specification (spec/units.md, spec/system-package.md).
//!
//! Units: nm, ps, amu (g/mol), kJ/mol, e, K; forces in kJ mol⁻¹ nm⁻¹; virial in kJ/mol.

/// Computes harmonic angle energy, forces, and virial contribution with error handling for degenerate geometries (R1.6).
pub mod angle;
/// Computes harmonic bond energy, forces, and virial contribution (R1.5).
pub mod bond;
/// Implements periodic cell validation, volume calculation, and minimum-image displacement in lower-triangular form (R1.4).
pub mod cell;
/// Computes pair parameters sigma_ij and eps_ij from per-type values using the geometric or Lorentz-Berthelot rule (R3.1).
pub mod combination_rule;
/// Defines CODATA 2018 physical constants and derived values for gas constant and Coulomb constant in project units (R1.1).
pub mod constants;
/// Coulomb pair interaction: energy, force, and virial for a single pair using the project Coulomb constant (R3.2).
pub mod coulomb_pair;
/// Computes the IUPAC dihedral angle phi in radians and its gradient with respect to the four atom positions (R2.1). Returns an error for collinear atoms where the angle is undefined.
pub mod dihedral_angle;
/// Error function erf and complementary error function erfc to near f64 precision (relative error below 1e-13), by a positive series and a continued fraction (R4.1).
pub mod erf;
/// Defines the library error type for invalid inputs such as non-finite coordinates, coincident atoms, index out of range, invalid cells, invalid parameters and non-neutral periodic systems (R1.3, R4.1).
pub mod errors;
/// Classic Ewald summation of point charges in a periodic cell: parameters, self term, correction for excluded and 1-4 pairs, and the total with forces and virial of every part (R4.1-R4.3).
pub mod ewald;
/// Real-space part of the Ewald sum over all images within the cutoff, with forces and virial (R4.1, R4.2).
pub mod ewald_real;
/// Reciprocal-space part of the Ewald sum with forces and virial (R4.1, R4.2).
pub mod ewald_reciprocal;
/// Builds the pair-classification map from a bond graph: pairs at shortest path length 1 or 2 are excluded, path length 3 are 1-4 scaled, all others are normal (R3.3).
pub mod exclusions;
/// Provides finite-difference test kit for validating forces and virials using central differences and strain deformations (R1.7). Uses function pointers to avoid trait/generic complexity while allowing
pub mod fd_test;
/// Total evaluation of a parameterized system: energy, forces, virial and named energy terms that sum to the energy, non-periodic or periodic with Ewald (R5.2).
pub mod forcefield;
/// Implements the harmonic improper term E = (1/2) k (xi - xi0)^2 with xi - xi0 wrapped to (-pi, pi], including energy, forces, and virial (R2.5).
pub mod harmonic_improper;
/// LJ 12-6 pair interaction: energy, force, and virial for a single pair given sigma_ij and eps_ij (R3.1).
pub mod lj_pair;
/// Compute the minimum-image displacement vector and its magnitude for a pair of atoms in a periodic cell. The displacement is r_i - r_j corrected by integer lattice translations so that each component l
pub mod min_image_disp;
/// Computes the minimum-image displacement between two positions in a periodic cell, and validates that a cutoff does not exceed half the smallest perpendicular cell height (R3.5).
pub mod minimum_image;
/// Evaluates all non-bonded pairs in a non-periodic system with exclusions and 1-4 scaling; returns energy split into LJ and Coulomb, forces, and virial (R3.4).
pub mod nonbonded_nonperiodic;
/// Evaluate all non-bonded pairs in a periodic system using the minimum image and a plain cutoff for LJ; Coulomb uses the same minimum image without cutoff. Applies combination rules to derive pair param
pub mod nonbonded_periodic;
/// Implements the OPLS-Fourier dihedral term E = (1/2)[V1(1+cos phi) + V2(1-cos 2phi) + V3(1+cos 3phi) + V4(1-cos 4phi)] with energy, forces, and virial (R2.3).
pub mod opls_fourier_dihedral;
/// Algebraic conversion from OPLS-Fourier coefficients V[4] to Ryckaert-Bellemans coefficients C[6]. Uses the identities cos(2phi) = 2cos^2(phi)-1, cos(3phi) = 4cos^3(phi)-3cos(phi), cos(4phi) = 8cos^4(p
pub mod opls_to_rb;
/// Given N atoms and a list of bonded pairs, compute the PairClass (Normal, Excluded, OneFour) for every pair (i, j) with i < j by BFS shortest-path on the bond graph. Isolates the pair-classification lo
pub mod pair_class_map;
/// Implements the periodic dihedral term E = k [1 + cos(n*phi - phi0)] with energy, forces, and virial (R2.2).
pub mod periodic_dihedral;
/// Pure trigonometric math for the Ryckaert-Bellemans potential: given the IUPAC dihedral angle phi (rad) and coefficients C[6], compute the energy E = sum_{n=0..5} C_n cos^n(psi) with psi = phi - pi, an
pub mod rb_dihedral_math;
/// The public-facing module for the Ryckaert-Bellemans dihedral term. Combines the existing dihedral_angle/dihedral_gradient (geometry) with rb_dihedral_math (energy and derivative) to produce energy, fo
pub mod ryckaert_bellemans_dihedral;
/// Data model of a parameterized system (molecule types, components, atom types, combination rule, 1-4 scale factors) and its expansion in canonical atom order (R5.1).
pub mod system;
/// Provides vector helpers for [f64; 3] including add, subtract, scale, dot, cross, and norm (R1.2).
pub mod vec3;
