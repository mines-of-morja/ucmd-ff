# M5 — Force-field description and total evaluation: design

Written by Claude (overseer) on 2026-11-17, at the author's request. The team did not plan M5.

## Modules

| Module | Content |
|---|---|
| `system` | Data model of a parameterized system as in spec/system-package.md §3–§5 (R5.1). It holds `AtomType`, `Atom`, `Bond`, `Angle`, `Dihedral` (forms: periodic, OPLS-Fourier, Ryckaert-Bellemans), `Improper` (forms: periodic, harmonic), `MoleculeType`, `Component`, `ForceFieldParameters` (combination rule, 1-4 scale factors, atom types) and `System`. `validate` checks the parameter ranges of §5.3, the indices (E_INDEX) and the references. `topology` expands the system into the canonical atom order of §4. |
| `forcefield` | `ForceField::new(&System)` validates the system and builds the topology and the pair lists. `evaluate(positions, periodic)` returns the energy, forces, virial and named energy terms (R5.2). |

## Evaluation

- **Bonded terms:** `bond`, `angle`, `dihedral`, `improper`, computed by the M1/M2 modules. Each
  term's virial is Σ_k (r_k − r_ref) ⊗ F_k.
- **Non-periodic** (`periodic = None`): `lj` and `coulomb` come from `nonbonded_nonperiodic`
  (M3): all pairs, with exclusions and 1-4 scaling.
- **Periodic** (`Some((cell, settings))`):
  - `lj` comes from `nonbonded_periodic` (M3), using the minimum image and a plain cutoff below
    h_min/2, with the charges switched off.
  - Coulomb comes from the Ewald sum (M4): `coulomb_real`, `coulomb_reciprocal`, `coulomb_self`
    and `coulomb_correction`.
  - Excluded pairs get scale 0 and 1-4 pairs the 1-4 Coulomb scale, both from the bond graph.
  - A bond at least h_min/2 long is an error (`MoleculeSplit`, E_MOLECULE_SPLIT).
- **Energy:** the sum of the terms in the order listed, so the terms add up to it exactly.

## Not included

- **Package reading:** reading `system.json` and `positions.npy` is the engine's task (ucmd-io).
  This crate has no external dependencies and no JSON reader. Angles in the data model are
  radians; the package gives degrees.
- **Constraints:** they belong to the integrator, not to the energy.
- **Performance work:** the pair loops are O(N²), and the bonded terms are evaluated one by one.

## Verification (tests/forcefield.rs)

- **Conformance tests of spec/force-provider.md §6 (R5.3),** on four systems: water-2
  non-periodic (R5.4), water-2 periodic with Ewald, and a test molecule with every bonded form
  and 1-4 pairs, non-periodic and periodic in a triclinic cell. The tests are:
  - finite-difference forces, scanning h from 1e-3 to 1e-7, with relative error ≤ 1e-5 at the
    best step;
  - finite-difference virial: all nine components when non-periodic, the six lower components
    with the cell strained when periodic;
  - energy terms that sum to the energy;
  - |ΣF| ≤ 1e-8 Σ|F|;
  - finite values and the right shapes;
  - bitwise determinism.
- **water-2 against an independent Python evaluation:** the Python script reads the fixture
  files directly. Every energy term matches to 1e-10 relative, and the forces and virial match
  its finite differences.
- **Test molecule against an independent Python evaluation:**
  - all six terms match to 1e-10, with IUPAC dihedrals, wrapped impropers and 8 1-4 pairs;
  - the periodic LJ matches, and so does the periodic Coulomb, computed as the oracle's Ewald sum
    plus the scaled 1-4 pairs.
- **NaCl crystal through the whole evaluator:** the fixture nacl-crystal (64 ions) gives the
  Madelung constant to 1e-8.
- **Canonical order and error paths:** both are tested, including invalid systems, the wrong
  number of positions, non-finite positions, a split molecule, an LJ cutoff that is too long,
  and a non-neutral periodic system.
- **Mutation check:** four deliberate bugs were planted: the 1-4 Coulomb scale taken from LJ, the
  angle vertex swapped, the 1-4 LJ left unscaled, and Coulomb counted twice in the periodic case.
  Each one makes at least one test fail.
