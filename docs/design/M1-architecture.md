# M1 — Foundations and harmonic bonded terms: module plan

Written by the architect role; checked by the design referee (verdict: revise).

Referee's remarks on the first draft:

- In module 7 (`fd_test`), the API for `numerical_virial` takes a `cell: &Cell` parameter, but the module's declared dependencies are only `vec3`. It must depend on `cell` to use the `Cell` type.
- In module 7 (`fd_test`), the signature for `numerical_virial` uses `energy_fn: fn(&[[f64; 3]]) -> f64`. This is inconsistent with the bonded modules (5 and 6) which take explicit atom positions (e.g., `r: [f64; 3]`) and parameters (e.g., `k`, `r0_nm`) rather than a full position array. The test kit must be able to call the specific analytic functions (like `bond_energy`) to compute the energy for finite differences. The signature should accept a closure or function pointer that matches the specific energy function's argument list, or the plan should specify how the test kit wraps the specific analytic functions into a generic callable interface consistent with the module APIs.
- In module 7 (`fd_test`), `numerical_force` takes `positions: &[[f64; 3]]` and `atom_idx`. However, the bonded modules (5 and 6) do not take a full position array; they take specific atom coordinates (e.g., `r`, `r0` for bond; `r1`, `r2`, `r3` for angle). The generic `numerical_force` signature cannot directly call `bond_energy` or `angle_energy` without a wrapper that extracts the relevant atoms from the full position array. The plan does not define how this mapping is handled, making the validation 'Verify numerical_force matches analytic force for harmonic bond' impossible to implement as described without additional glue code not specified in the API.
- In module 7 (`fd_test`), `numerical_virial` takes `strain_idx: (usize, usize)`. The milestone requirement R1.7 specifies 'numerical virial from the energy under a small homogeneous strain'. Homogeneous strain is typically a 3x3 matrix (9 components). Using a single `(usize, usize)` pair implies computing the derivative with respect to a single strain component. While this is a valid way to compute individual virial components, the API should clarify if it returns the derivative for that specific component or if it expects the caller to loop over all 9 components. More critically, the signature `numerical_virial(...) -> f64` returns a single float, which is consistent with computing one component. However, the validation 'Verify numerical_virial matches analytic virial for harmonic bond' implies comparing the full 3x3 matrix. The API should either return the full matrix (requiring 9 calls internally or a different signature) or the validation description should be adjusted to reflect that the test calls the function 9 times. Given the requirement to compare analytic and numerical results, the current API is awkward for direct comparison of the full virial tensor. A better API might be `numerical_virial_matrix` or the test kit should provide a helper that assembles the full matrix from 9 calls to `numerical_virial`.

## `constants`

Defines CODATA 2018 physical constants and derived values for gas constant and Coulomb constant in project units (R1.1).

- depends on: nothing
- estimated lines: 40

```rust
pub const E_CHARGE_C: f64
pub const AVOGADRO_MOL: f64
pub const K_B_J_K: f64
pub const EPS0_F_M: f64
pub const R_GAS_KJ_MOL_K: f64
pub const COULOMB_KJ_MOL_NM_E2: f64
```

Validation:

- Compare R_GAS_KJ_MOL_K with 0.00831446261815324 to 1e-12 relative tolerance.
- Compare COULOMB_KJ_MOL_NM_E2 with 138.93545764438196 to 1e-12 relative tolerance.

## `vec3`

Provides vector helpers for [f64; 3] including add, subtract, scale, dot, cross, and norm (R1.2).

- depends on: nothing
- estimated lines: 60

```rust
pub fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3]
pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3]
pub fn scale(v: [f64; 3], s: f64) -> [f64; 3]
pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3]
pub fn norm(v: [f64; 3]) -> f64
```

Validation:

- Verify dot product with known orthogonal vectors yields 0.
- Verify cross product magnitude equals area of parallelogram for known vectors.
- Verify norm of [3, 4, 0] is 5.0.

## `errors`

Defines the library error type for invalid inputs such as non-finite coordinates, coincident atoms, index out of range, and invalid cells (R1.3).

- depends on: nothing
- estimated lines: 30

```rust
pub enum UcmdError { NonFiniteCoordinate, CoincidentAtoms, IndexOutOfRange, InvalidCell }
```

Validation:

- Verify error variants are constructible and matchable in downstream modules.

## `cell`

Implements periodic cell validation, volume calculation, and minimum-image displacement in lower-triangular form (R1.4).

- depends on: vec3, errors
- estimated lines: 120

```rust
pub struct Cell { pub ax: f64, pub bx: f64, pub by: f64, pub cx: f64, pub cy: f64, pub cz: f64 }
pub fn validate_cell(cell: &Cell) -> Result<(), UcmdError>
pub fn cell_volume(cell: &Cell) -> f64
pub fn minimum_image(r_i: [f64; 3], r_j: [f64; 3], cell: &Cell) -> [f64; 3]
```

Validation:

- Compare minimum_image result with brute-force search over 27 neighboring images for random points.
- Verify cell_volume equals ax*by*cz for valid cells.
- Verify validate_cell rejects non-positive ax, by, or cz.

## `bond`

Computes harmonic bond energy, forces, and virial contribution (R1.5).

- depends on: vec3, errors
- estimated lines: 100

```rust
pub fn bond_energy(r_i: [f64; 3], r_j: [f64; 3], r0_nm: f64, k_kjmol_nm2: f64) -> Result<f64, UcmdError>
pub fn bond_forces(r_i: [f64; 3], r_j: [f64; 3], r0_nm: f64, k_kjmol_nm2: f64) -> Result<([f64; 3], [f64; 3]), UcmdError>
pub fn bond_virial(r_i: [f64; 3], r_j: [f64; 3], r0_nm: f64, k_kjmol_nm2: f64) -> Result<[f64; 3; 3], UcmdError>
```

Validation:

- Verify force equals -k(r - r0) * (r - r0)/|r - r0| for arbitrary non-degenerate geometries.
- Verify virial equals outer product of displacement and force on atom i.
- Verify error is returned when atoms coincide (distance = 0).
- Verify momentum conservation: sum of forces on both atoms is zero.

## `angle`

Computes harmonic angle energy, forces, and virial contribution with error handling for degenerate geometries (R1.6).

- depends on: vec3, errors
- estimated lines: 150

```rust
pub fn angle_energy(r1: [f64; 3], r2: [f64; 3], r3: [f64; 3], theta0: f64, k_kjmol_rad2: f64) -> Result<f64, UcmdError>
pub fn angle_forces(r1: [f64; 3], r2: [f64; 3], r3: [f64; 3], theta0: f64, k_kjmol_rad2: f64) -> Result<([f64; 3], [f64; 3], [f64; 3]), UcmdError>
pub fn angle_virial(r1: [f64; 3], r2: [f64; 3], r3: [f64; 3], theta0: f64, k_kjmol_rad2: f64) -> Result<[f64; 3; 3], UcmdError>
```

Validation:

- Verify energy equals 0.5*k*(theta - theta0)^2 for known angles.
- Verify forces via finite differences of energy.
- Verify error is returned when geometry is straight or degenerate (undefined force direction).
- Verify virial via finite differences of energy under small strain.

## `fd_test`

Provides finite-difference test kit for validating forces and virials using central differences and strain deformations (R1.7). Uses function pointers to avoid trait/generic complexity while allowing specific analytic functions to be passed in.

- depends on: vec3, cell
- estimated lines: 200

```rust
pub fn numerical_force_bond(energy_fn: fn([f64; 3], [f64; 3], f64, f64) -> Result<f64, UcmdError>, r_i: [f64; 3], r_j: [f64; 3], r0_nm: f64, k_kjmol_nm2: f64, atom_idx: usize, coord_idx: usize, h: f64) -> f64
pub fn numerical_force_angle(energy_fn: fn([f64; 3], [f64; 3], [f64; 3], f64, f64) -> Result<f64, UcmdError>, r1: [f64; 3], r2: [f64; 3], r3: [f64; 3], theta0: f64, k_kjmol_rad2: f64, atom_idx: usize, coord_idx: usize, h: f64) -> f64
pub fn numerical_virial_bond(energy_fn: fn([f64; 3], [f64; 3], f64, f64) -> Result<f64, UcmdError>, r_i: [f64; 3], r_j: [f64; 3], r0_nm: f64, k_kjmol_nm2: f64, cell: &Cell, strain_row: usize, strain_col: usize, h: f64) -> f64
pub fn numerical_virial_angle(energy_fn: fn([f64; 3], [f64; 3], [f64; 3], f64, f64) -> Result<f64, UcmdError>, r1: [f64; 3], r2: [f64; 3], r3: [f64; 3], theta0: f64, k_kjmol_rad2: f64, cell: &Cell, strain_row: usize, strain_col: usize, h: f64) -> f64
pub fn compare_force(analytic: [f64; 3], numerical: [f64; 3], tol: f64) -> bool
pub fn compare_virial(analytic: [f64; 3; 3], numerical: [f64; 3; 3], tol: f64) -> bool
```

Validation:

- Verify numerical_force_bond matches analytic force for harmonic bond with h=1e-5.
- Verify numerical_virial_bond matches analytic virial component for harmonic bond with h=1e-5.
- Verify compare_force returns true when analytic and numerical are within tolerance.
- Verify compare_virial returns true when analytic and numerical matrices are within tolerance.
