# M4 — Ewald summation for periodic electrostatics: design

Written by Claude (overseer) on 2026-11-17, after the author lifted the rule that only the local
model writes this repository. It replaces the team's module plan of 2026-11-17 (git history:
`design(M4): module plan with 7 modules`). That plan was built on research answers that the
machine check had found wrong, and its errors are not repeated here:

- **Exclusion sign:** the plan gave the excluded-pair correction as +k_e q_i q_j erf(αr)/r. The
  correction must remove that term. The plan's design referee then asked to replace R4.3's
  exclusion identity with a formula derived from the wrong sign.
- **Madelung geometry:** the test put 2 ions in a cubic cell of side a, with Na⁺ at the origin
  and Cl⁻ at (a/2, 0, 0). That is not the rock-salt structure: the ions alternate along x only.
- **Trace check:** the plan checked tr(W_recip) = E_recip/V. That is wrong in dimension and in
  value: tr W_recip = Σ_k e_k (1 − k²/2α²).

## Modules

| Module | Content |
|---|---|
| `erf` | `erf`, `erfc` to relative error below 1e-13 (R4.1): a positive series for \|x\| < 1, a continued fraction (modified Lentz) for x ≥ 1, symmetry for x < 0, and e^(−x²) split exactly with a fused multiply-add. No fitted coefficients. |
| `ewald` | `EwaldParameters` (α, real-space cutoff, reciprocal cutoff; `from_tolerance`), `ScaledPair`, `EwaldPart`, `EwaldResult`, `validate_system` (lengths, finiteness, cell, neutrality within 1e-6 e), `scaled_pairs_from_bonds`, `ewald_self`, `ewald_correction`, `ewald` (all parts and the total). |
| `ewald_real` | Real-space sum over all images within the cutoff. The search is complete for any cutoff, including a cutoff above half the cell height, and it includes every atom's own images. A scaled pair leaves out its nearest image. |
| `ewald_reciprocal` | Reciprocal-space sum over 0 < \|k\| < k_c (half space), with structure factors, forces and the analytic virial. |

## Equations

All energies in kJ/mol, with k_e = `COULOMB_KJ_MOL_NM_E2` and tinfoil boundary conditions. The
pair displacement is r_ij = r_i − r_j, and F_ij is the force on i due to j.

- Real space: k_e Σ_{i<j} Σ_n q_i q_j erfc(α r)/r over r = \|r_ij + n\| < r_c, plus
  ½ k_e Σ_i q_i² Σ_{n≠0} erfc(α\|n\|)/\|n\|. The pair force is F_i = g d with
  g = k_e q_i q_j [erfc(αr)/r + (2α/√π) e^(−α²r²)]/r². The virial is Σ d ⊗ F_i, and
  ½ Σ_n n ⊗ F(n) for an atom's own images.
- Reciprocal space: (2π k_e/V) Σ_{k≠0} f(k)\|S(k)\|², with f = e^(−k²/4α²)/k².
  - Force: F_j = (4π k_e q_j/V) Σ_k f k [C sin(k·r_j) − S cos(k·r_j)].
  - Virial: W_ab = (2π k_e/V) Σ_k f \|S\|² [δ_ab − 2 k_a k_b (1/k² + 1/4α²)]. This follows from
    k·r being invariant under the strain, ∂k_c/∂ε_ab = −δ_ca k_b and ∂V/∂ε_ab = V δ_ab.
- Self term: −k_e (α/√π) Σ q_i². Its forces and virial are zero.
- Correction for a scaled pair (i, j, s), with r the minimum-image distance:
  k_e q_i q_j [s − erf(αr)]/r.
  - s = 0 for an excluded pair: it then has no Coulomb interaction at all.
  - s = the 1-4 Coulomb scale factor for a 1-4 pair.
  - The converged total therefore changes by (s − 1) k_e q_i q_j / r, for any α (R4.3).

## Verification (tests/erf.rs, tests/ewald.rs)

- **erf/erfc:** reference values from a 400-digit evaluation; relative error below 1e-13 at 19
  points from 0 to 26 and at negative arguments. A dense scan of 20,000 points against glibc
  (outside the suite) gave a worst relative error of 3.3e-15.
- **Ewald parts:** the energies of all four parts match an independent Python Ewald sum (the
  harness's oracle) to 1e-9 kJ/mol. The total forces and the six lower virial components match
  its finite differences to 1e-5.
- **Madelung constant:** the rock-salt value from the 8-ion cubic cell is correct to 1e-8
  relative, at two splittings.
- **α independence:** the total energy, forces and virial agree to 1e-9 relative for α from 3.5
  to 9.4 nm⁻¹, with cutoffs from 0.6 nm to 1.6 nm.
- **Exclusions:** the exclusion identity −k_e q_i q_j / r (minimum image) and the 1-4 identity
  (s − 1) k_e q_i q_j / r hold to 1e-10 relative at two values of α.
- **Finite differences:** the forces and the six lower virial components of the real-space,
  reciprocal, correction and total parts match finite differences from the test kit (R1.7).
  The real-space part is also checked with a cutoff above all lattice vectors. Every virial is
  symmetric.
- **Physical consistency:** momentum is conserved, and moving atoms by lattice vectors leaves
  the result unchanged.
- **Error handling:** every error path is tested.
- **Mutation check:** six deliberate bugs were planted, one at a time: the reciprocal force sign,
  the reciprocal virial factor, the self-image energy and virial halves, the correction sign, and
  the minimum image in the correction. Each one makes at least one test fail.
