# M2 — Dihedral and improper terms: module plan

Written by the architect role; checked by the design referee (verdict: ok).

## `dihedral_angle`

Computes the IUPAC dihedral angle phi in radians and its gradient with respect to the four atom positions (R2.1). Returns an error for collinear atoms where the angle is undefined.

- depends on: vec3, errors
- estimated lines: 180

```rust
pub fn dihedral_angle(r_i: [f64; 3], r_j: [f64; 3], r_k: [f64; 3], r_l: [f64; 3]) -> Result<f64, UcmdError>
pub fn dihedral_gradient(r_i: [f64; 3], r_j: [f64; 3], r_k: [f64; 3], r_l: [f64; 3]) -> Result<([f64; 3], [f64; 3], [f64; 3], [f64; 3]), UcmdError>
```

Validation:

- Analytic check: trans configuration (phi = pi) and cis configuration (phi = 0) for known planar geometries.
- Finite difference: central difference of dihedral_angle with respect to each position component matches dihedral_gradient within 1e-8.
- Error check: collinear atoms (e.g., all four on a line) return Err(UcmdError::CoincidentAtoms) or similar invalid geometry error.

## `periodic_dihedral`

Implements the periodic dihedral term E = k [1 + cos(n*phi - phi0)] with energy, forces, and virial (R2.2).

- depends on: dihedral_angle, vec3, errors
- estimated lines: 150

```rust
pub fn periodic_dihedral_energy(r_i: [f64; 3], r_j: [f64; 3], r_k: [f64; 3], r_l: [f64; 3], n: i32, phi0: f64, k_kjmol: f64) -> Result<f64, UcmdError>
pub fn periodic_dihedral_forces(r_i: [f64; 3], r_j: [f64; 3], r_k: [f64; 3], r_l: [f64; 3], n: i32, phi0: f64, k_kjmol: f64) -> Result<([f64; 3], [f64; 3], [f64; 3], [f64; 3]), UcmdError>
pub fn periodic_dihedral_virial(r_i: [f64; 3], r_j: [f64; 3], r_k: [f64; 3], r_l: [f64; 3], n: i32, phi0: f64, k_kjmol: f64) -> Result<[[f64; 3]; 3], UcmdError>
```

Validation:

- Analytic check: energy at phi = phi0 is k, at phi = phi0 + pi is 0 (for n=1).
- Finite difference: forces match central difference of energy with respect to positions.
- Finite difference: virial matches central difference of energy under strain deformation.

## `opls_fourier_dihedral`

Implements the OPLS-Fourier dihedral term E = (1/2)[V1(1+cos phi) + V2(1-cos 2phi) + V3(1+cos 3phi) + V4(1-cos 4phi)] with energy, forces, and virial (R2.3).

- depends on: dihedral_angle, vec3, errors
- estimated lines: 160

```rust
pub fn opls_fourier_dihedral_energy(r_i: [f64; 3], r_j: [f64; 3], r_k: [f64; 3], r_l: [f64; 3], v1: f64, v2: f64, v3: f64, v4: f64) -> Result<f64, UcmdError>
pub fn opls_fourier_dihedral_forces(r_i: [f64; 3], r_j: [f64; 3], r_k: [f64; 3], r_l: [f64; 3], v1: f64, v2: f64, v3: f64, v4: f64) -> Result<([f64; 3], [f64; 3], [f64; 3], [f64; 3]), UcmdError>
pub fn opls_fourier_dihedral_virial(r_i: [f64; 3], r_j: [f64; 3], r_k: [f64; 3], r_l: [f64; 3], v1: f64, v2: f64, v3: f64, v4: f64) -> Result<[[f64; 3]; 3], UcmdError>
```

Validation:

- Analytic check: energy at phi = 0 for V1=1, V2=V3=V4=0 is 1.0.
- Finite difference: forces match central difference of energy with respect to positions.
- Finite difference: virial matches central difference of energy under strain deformation.

## `ryckaert_bellemans_dihedral`

Implements the Ryckaert-Bellemans dihedral term E = sum_{n=0..5} C_n cos^n(psi) with psi = phi - pi, including energy, forces, virial, and conversion from OPLS-Fourier coefficients (R2.4).

- depends on: dihedral_angle, vec3, errors
- estimated lines: 200

```rust
pub fn ryckaert_bellemans_dihedral_energy(r_i: [f64; 3], r_j: [f64; 3], r_k: [f64; 3], r_l: [f64; 3], c: [f64; 6]) -> Result<f64, UcmdError>
pub fn ryckaert_bellemans_dihedral_forces(r_i: [f64; 3], r_j: [f64; 3], r_k: [f64; 3], r_l: [f64; 3], c: [f64; 6]) -> Result<([f64; 3], [f64; 3], [f64; 3], [f64; 3]), UcmdError>
pub fn ryckaert_bellemans_dihedral_virial(r_i: [f64; 3], r_j: [f64; 3], r_k: [f64; 3], r_l: [f64; 3], c: [f64; 6]) -> Result<[[f64; 3]; 3], UcmdError>
pub fn opls_to_rb(v: [f64; 4]) -> [f64; 6]
```

Validation:

- Conversion check: for a scan of phi over [0, 2pi], energy from OPLS-Fourier equals energy from Ryckaert-Bellemans using converted coefficients within 1e-10.
- Finite difference: forces match central difference of energy with respect to positions.
- Finite difference: virial matches central difference of energy under strain deformation.

## `harmonic_improper`

Implements the harmonic improper term E = (1/2) k (xi - xi0)^2 with xi - xi0 wrapped to (-pi, pi], including energy, forces, and virial (R2.5).

- depends on: dihedral_angle, vec3, errors
- estimated lines: 150

```rust
pub fn harmonic_improper_energy(r_i: [f64; 3], r_j: [f64; 3], r_k: [f64; 3], r_l: [f64; 3], xi0: f64, k_kjmol_rad2: f64) -> Result<f64, UcmdError>
pub fn harmonic_improper_forces(r_i: [f64; 3], r_j: [f64; 3], r_k: [f64; 3], r_l: [f64; 3], xi0: f64, k_kjmol_rad2: f64) -> Result<([f64; 3], [f64; 3], [f64; 3], [f64; 3]), UcmdError>
pub fn harmonic_improper_virial(r_i: [f64; 3], r_j: [f64; 3], r_k: [f64; 3], r_l: [f64; 3], xi0: f64, k_kjmol_rad2: f64) -> Result<[[f64; 3]; 3], UcmdError>
```

Validation:

- Analytic check: energy is zero when xi = xi0, and equals 0.5*k*(delta)^2 for known deviations.
- Finite difference: forces match central difference of energy with respect to positions.
- Finite difference: virial matches central difference of energy under strain deformation.
- Wrapping check: verify that xi - xi0 is wrapped to (-pi, pi] before squaring.
