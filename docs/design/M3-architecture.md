# M3 — Non-bonded pair interactions: module plan

Written by the architect role; checked by the design referee (verdict: revise).

Referee's remarks on the first draft:

- minimum_image must return Result<[f64;3], UcmdError> instead of a bare [f64;3]. The Cell struct does not enforce ax,by,cz>0, so an invalid cell would cause a silent division-by-zero (NaN) in the fractional-coordinate wrapping step, violating project rule 6 (errors propagate). Change the signature to pub fn minimum_image(cell: &Cell, r_i: [f64;3], r_j: [f64;3]) -> Result<[f64;3], UcmdError> and return Err(UcmdError::InvalidCell) when any of ax, by
- cz is non-positive or the cell matrix is singular.

## `lj_pair`

LJ 12-6 pair interaction: energy, force, and virial for a single pair given sigma_ij and eps_ij (R3.1).

- depends on: vec3, errors
- estimated lines: 120

```rust
pub fn lj_pair_energy(sigma_nm: f64, eps_kjmol: f64, r_i: [f64; 3], r_j: [f64; 3]) -> Result<f64, UcmdError>
pub fn lj_pair_forces(sigma_nm: f64, eps_kjmol: f64, r_i: [f64; 3], r_j: [f64; 3]) -> Result<([f64; 3], [f64; 3]), UcmdError>
pub fn lj_pair_virial(sigma_nm: f64, eps_kjmol: f64, r_i: [f64; 3], r_j: [f64; 3]) -> Result<[[f64; 3]; 3], UcmdError>
```

Validation:

- Analytic: at r = 2^(1/6)*sigma, energy = -eps and force = zero (Rahman-1964 argon: sigma=0.34 nm, eps=120*K_B converted to kJ/mol).
- Finite-difference: perturb r_i along each axis by h=1e-8 nm; compare (E(r+dh)-E(r-dh))/(2h) against analytic force component, tolerance 1e-6 kJ mol^-1 nm^-1.
- Virial: for a single pair along x-axis, W_xx = r*F_x; verify by finite-differencing E with respect to a diagonal strain.
- Sum of forces vanishes: F_i + F_j = 0 for a single pair.

## `coulomb_pair`

Coulomb pair interaction: energy, force, and virial for a single pair using the project Coulomb constant (R3.2).

- depends on: vec3, constants, errors
- estimated lines: 100

```rust
pub fn coulomb_pair_energy(q_i: f64, q_j: f64, r_i: [f64; 3], r_j: [f64; 3]) -> Result<f64, UcmdError>
pub fn coulomb_pair_forces(q_i: f64, q_j: f64, r_i: [f64; 3], r_j: [f64; 3]) -> Result<([f64; 3], [f64; 3]), UcmdError>
pub fn coulomb_pair_virial(q_i: f64, q_j: f64, r_i: [f64; 3], r_j: [f64; 3]) -> Result<[[f64; 3]; 3], UcmdError>
```

Validation:

- Analytic: two unit charges at r=1 nm give E = COULOMB_KJ_MOL_NM_E2 = 138.93545764438196 kJ/mol.
- Finite-difference: perturb r_i along each axis by h=1e-8 nm; compare against analytic F = k_e*q_i*q_j/r^2 * r_hat, tolerance 1e-6.
- Virial: W = r_ij outer F_ij; verify by finite-differencing E with respect to diagonal strain.
- Sum of forces vanishes: F_i + F_j = 0.

## `combination_rule`

Computes pair parameters sigma_ij and eps_ij from per-type values using the geometric or Lorentz-Berthelot rule (R3.1).

- depends on: nothing
- estimated lines: 50

```rust
#[derive(Debug, Clone, Copy, PartialEq)] pub enum CombinationRule { Geometric, LorentzBerthelot }
pub fn combine_sigma(rule: CombinationRule, sigma_i: f64, sigma_j: f64) -> f64
pub fn combine_eps(eps_i: f64, eps_j: f64) -> f64
```

Validation:

- Identical types: for sigma_i = sigma_j and eps_i = eps_j, both rules give sigma_ij = sigma_i, eps_ij = eps_i.
- Geometric: sigma_ij = sqrt(sigma_i * sigma_j); verify with arbitrary values (e.g., 0.3 and 0.5 → sqrt(0.15)).
- Lorentz-Berthelot: sigma_ij = (sigma_i + sigma_j)/2; verify with arbitrary values (e.g., 0.3 and 0.5 → 0.4).
- eps is geometric in both rules: eps_ij = sqrt(eps_i * eps_j).

## `exclusions`

Builds the pair-classification map from a bond graph: pairs at shortest path length 1 or 2 are excluded, path length 3 are 1-4 scaled, all others are normal (R3.3).

- depends on: errors
- estimated lines: 130

```rust
#[derive(Debug, Clone, Copy, PartialEq)] pub enum PairClass { Normal, Excluded, OneFour }
pub fn build_exclusion_map(n_atoms: usize, bonds: &[(usize, usize)]) -> Result<Vec<Vec<PairClass>>, UcmdError>
```

Validation:

- Linear chain of 4 atoms (bonds 0-1, 1-2, 2-3): pair (0,3) is OneFour; pairs (0,1),(1,2),(2,3) are Excluded; pair (0,2) is Excluded.
- Ring of 4 atoms (bonds 0-1, 1-2, 2-3, 3-0): pair (0,2) is Excluded (shortest path = 2 via either 0-1-2 or 0-3-2); all adjacent pairs are Excluded.
- Linear chain of 5 atoms: pair (0,3) is OneFour, pair (0,4) is Normal.
- Two disconnected molecules: all inter-molecular pairs are Normal.
- Invalid bond index (>= n_atoms) returns Err(UcmdError::IndexOutOfRange).

## `minimum_image`

Computes the minimum-image displacement between two positions in a periodic cell, and validates that a cutoff does not exceed half the smallest perpendicular cell height (R3.5).

- depends on: cell, vec3, errors
- estimated lines: 140

```rust
pub fn minimum_image(cell: &Cell, r_i: [f64; 3], r_j: [f64; 3]) -> Result<[f64; 3], UcmdError>
pub fn validate_cutoff(cell: &Cell, cutoff_nm: f64) -> Result<(), UcmdError>
```

Validation:

- Orthogonal cubic cell (ax=by=cz=1 nm): r_i=(0.9,0,0), r_j=(0.1,0,0) → d_min=(-0.2,0,0) nm.
- Brute-force: for random point pairs in a random valid triclinic cell, compare minimum_image output against the shortest of all 27 image displacements (±1 in each fractional direction).
- Cutoff bound: for a sheared cell (ax=1, bx=0.9, by=1, cx=0.9, cy=0.9, cz=1), compute h_a, h_b, h_c; verify that cutoff = min(h)/2 + 0.01 is rejected with Err(InvalidCell).
- Orthogonal cell: cutoff = ax/2 - epsilon passes; cutoff = ax/2 + epsilon fails.
- Invalid cell (ax=0 or by=-1): minimum_image returns Err(InvalidCell) rather than producing NaN.

## `nonbonded_nonperiodic`

Evaluates all non-bonded pairs in a non-periodic system with exclusions and 1-4 scaling; returns energy split into LJ and Coulomb, forces, and virial (R3.4).

- depends on: lj_pair, coulomb_pair, combination_rule, exclusions, vec3, errors
- estimated lines: 200

```rust
#[derive(Debug, Clone, PartialEq)] pub struct NonbondedInput { pub positions: Vec<[f64; 3]>, pub types: Vec<usize>, pub charges: Vec<f64>, pub lj_sigma: Vec<f64>, pub lj_eps: Vec<f64>, pub rule: CombinationRule, pub bonds: Vec<(usize, usize)>, pub lj_scale_14: f64, pub coul_scale_14: f64 }
#[derive(Debug, Clone, PartialEq)] pub struct NonbondedResult { pub energy_lj: f64, pub energy_coulomb: f64, pub forces: Vec<[f64; 3]>, pub virial: [[f64; 3]; 3] }
pub fn evaluate_nonperiodic(input: &NonbondedInput) -> Result<NonbondedResult, UcmdError>
```

Validation:

- Two-atom system with known LJ and Coulomb parameters: verify energy equals sum of lj_pair_energy + coulomb_pair_energy; forces and virial match pair-level results.
- Three-atom linear chain (bonds 0-1, 1-2): pair (0,2) is 1-4; verify its LJ energy is scaled by lj_scale_14 and Coulomb by coul_scale_14.
- Sum of all forces vanishes: sum of forces over all atoms = [0,0,0] (Newton's third law).
- Finite-difference: perturb each atom position along each axis by h=1e-8 nm; compare (E(r+dh)-E(r-dh))/(2h) against the returned force component, tolerance 1e-5 kJ mol^-1 nm^-1.
- Virial: finite-difference E with respect to a uniform diagonal strain eps; compare against the returned virial tensor, tolerance 1e-5 kJ/mol.

## `nonbonded_periodic`

Evaluates all non-bonded pairs in a periodic system using the minimum image and a plain cutoff for LJ; Coulomb uses the same minimum image without cutoff (R3.5).

- depends on: lj_pair, coulomb_pair, combination_rule, exclusions, minimum_image, cell, vec3, errors
- estimated lines: 180

```rust
pub fn evaluate_periodic(input: &NonbondedInput, cell: &Cell, lj_cutoff_nm: f64) -> Result<NonbondedResult, UcmdError>
```

Validation:

- Orthogonal cubic cell, two atoms near a boundary: verify that the minimum-image displacement is used (energy matches the non-periodic result computed with the wrapped distance).
- Cutoff: two atoms separated by more than lj_cutoff_nm in minimum-image distance do not contribute to LJ energy; verify energy_lj = 0 for that pair while Coulomb still contributes.
- Cutoff validation: a cutoff exceeding half the smallest perpendicular height returns Err(InvalidCell).
- Finite-difference: perturb positions and verify forces match, tolerance 1e-5.
- Sum of all forces vanishes.
- Translation invariance: shifting all positions by the same lattice vector gives the same energy, forces, and virial.
