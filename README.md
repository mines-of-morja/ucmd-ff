# ucmd-ff

A classical force-field engine for [UnitcellMD](https://github.com/mines-of-morja/UnitcellMD).
Given a parameterized system, positions and an optional periodic cell, it returns:

- the potential energy (kJ/mol);
- the forces F = −∂E/∂r (kJ mol⁻¹ nm⁻¹);
- the virial W = −∂E/∂ε (kJ/mol);
- the energy split into named terms.

It is a plain Rust library with no dependencies. ucmd-engine uses it as one of its
interchangeable in-process force providers. In `run.toml`, `kind = "ucmd-ff"` selects it, and
`kind = "force-field"` selects the engine's own provider.

## What it computes

The functional forms and conventions follow the UnitcellMD specification
(`spec/units.md`, `spec/system-package.md` §5):

- **Bonded terms:**
  - harmonic bonds and angles;
  - periodic, OPLS-Fourier and Ryckaert-Bellemans dihedrals (IUPAC convention);
  - periodic and harmonic impropers.
- **Lennard-Jones 12-6**, with the geometric or Lorentz-Berthelot combination rule.
- **Coulomb** with the CODATA 2018 constant.
- **Exclusions** from the bond graph (path lengths 1 and 2), and scaled 1-4 pairs (path length 3).
- **Non-periodic systems:** all pairs, as a direct sum.
- **Periodic systems** (triclinic cells in reduced lower-triangular form):
  - Lennard-Jones with the minimum image and a plain cutoff below half the smallest cell height;
  - Coulomb by classic Ewald summation (real space, reciprocal space, self term, and the
    correction for excluded and 1-4 pairs), each part with forces and virial;
  - `erf`/`erfc` to a relative error below 1e-13.

Energy terms:
- always: `bond`, `angle`, `dihedral`, `improper`, `lj`;
- non-periodic: `coulomb`;
- periodic: `coulomb_real`, `coulomb_reciprocal`, `coulomb_self`, `coulomb_correction`.

The energy is the sum of the terms.

```rust
use ucmd_ff::forcefield::{ForceField, PeriodicSettings};
use ucmd_ff::ewald::EwaldParameters;

let ff = ForceField::new(&system)?;                  // system: ucmd_ff::system::System
let open = ff.evaluate(&positions, None)?;            // non-periodic
let settings = PeriodicSettings {
    lj_cutoff_nm: 1.0,
    ewald: EwaldParameters::from_tolerance(1.0, 1e-6)?,
};
let periodic = ff.evaluate(&positions, Some((&cell, &settings)))?;
println!("{} kJ/mol, terms {:?}", periodic.energy_kjmol, periodic.energy_terms);
```

## Limits

ucmd-ff is built to be correct before fast:

- the pair loops scale as N²;
- the Ewald sum is the classic one (no PME);
- there are no neighbour lists, no LJ shift and no dispersion correction.

For large production systems, use ucmd-engine's `force-field` provider (PME, neighbour lists,
threads). ucmd-ff is useful as a second, independently written implementation that can check it.
Constraints are not part of the energy; the integrator handles them.

## Verification

`cargo test` runs 227 tests:

- finite-difference forces and virial for every term, using the crate's own test kit;
- the conformance tests of `spec/force-provider.md` §6 on the UnitcellMD fixture water-2 and on a
  molecule with every bonded form, both non-periodic and periodic;
- the Madelung constant of rock salt to 1e-8;
- independence of the total energy from the Ewald splitting parameter;
- the exclusion identity;
- `erf`/`erfc` against 400-digit reference values;
- term-by-term comparisons with independent Python evaluations.

The evidence for each milestone is in `docs/gates/`, and the designs are in `docs/design/`.

## How it was built

ucmd-ff was the test case of an experiment: can a small local language model build a numerical
library when a harness breaks the work into checked steps? The harness was *ucmd-agent*.

- **Milestones M1–M3** (foundations, bonded terms, non-bonded pairs) were written by a local
  27B-parameter model.
  - Every step was checked by machine: compiler, tests, recomputed expected values,
    deliberately planted bugs, two referee rulings before a test could change.
  - An overseer, Claude, steered by changing only the harness. The gate reports of these
    milestones are the harness's records.
- **M4 (Ewald) and M5 (data model, total evaluation)** were written by Claude, at the author's
  request.

The goal and the requirements are in `docs/mission.md`.

## Development

```sh
cargo fmt --all --check
cargo clippy --all-targets -- -D warnings
cargo test
```

Units: nm, ps, amu, kJ/mol, e, K. Licensed under the Apache License, Version 2.0 (`LICENSE`).
