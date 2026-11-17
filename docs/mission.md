# Mission: `ucmd-ff`, the UnitcellMD force-field engine

UnitcellMD is a modular molecular-dynamics project for large periodic systems; its first target
is Li-ion battery electrolytes (LiPF6 in EC:EMC). Its one rule: the engine owns the dynamics,
force providers own the physics.

`ucmd-ff` is the project's own force-field engine: a self-contained Rust library. Given atom
positions, an optional periodic cell and an explicit classical force-field description, it
returns the potential energy (kJ/mol), the forces F = -dE/dr (kJ mol^-1 nm^-1) and the virial
W = -dE/d(eps) (kJ/mol), with the energy split into named terms.

Functional forms (UnitcellMD spec, system-package.md §5): harmonic bonds and angles; periodic,
OPLS-Fourier and Ryckaert-Bellemans dihedrals; harmonic impropers; Lennard-Jones 12-6 with
geometric or Lorentz-Berthelot combination; Coulomb interaction with exclusions and scaled 1-4
pairs; Ewald summation for periodic systems.

Correct before fast. Every term is proven by tests: central finite differences of the energy for
forces and virial, analytic cases, invariances, and reference values from the literature.

Out of scope: integrators, thermostats, file formats, parallelism, performance tuning, and
parameter values for real molecules (those come from cited sources, never from memory).
The library stands alone and uses only the Rust standard library.

---

PROJECT RULES (from the UnitcellMD specification; they override anything found in papers)
1. Units, everywhere: length nm, time ps, mass amu (g/mol), energy kJ/mol, charge e, temperature K,
   force kJ mol^-1 nm^-1. Angles are radians inside the code. [spec/units.md §1]
2. Signs: the force is F_i = -dE/dr_i. The virial is W_ab = -dE/d(eps_ab) for the deformation
   r -> r(1 + eps); for pair forces W = sum over pairs i<j of r_ij (outer product) F_ij, with
   r_ij = r_i - r_j and F_ij the force on i due to j. [spec/units.md §4]
3. Geometry: positions are unwrapped Cartesian [N][3]. A periodic cell is a 3x3 matrix whose ROWS
   are the lattice vectors a, b, c in lower-triangular form: a = (ax, 0, 0), b = (bx, by, 0),
   c = (cx, cy, cz), with ax, by, cz > 0. [spec/units.md §3]
4. Physical constants are defined in ONE module from the CODATA 2018 values; no rounded constant
   is typed anywhere else. [spec/units.md §2]
5. Never invent a physical parameter. A number in a test is either arbitrary test geometry, said
   to be arbitrary, or a parameter cited from a source. [UnitcellMD CLAUDE.md]
6. Errors propagate: never replace an impossible or failed evaluation with zero, a default or a
   stale value; return an error. [UnitcellMD CLAUDE.md, invariant 7]
7. Tests are never weakened: no loosened tolerance, no skipped or deleted test to get a pass.
   [UnitcellMD CLAUDE.md]
8. Code: Rust, f64 only, no unsafe, no external crates, no unwrap/expect/panic in library code;
   every public item has a doc comment that states its units. [UnitcellMD CLAUDE.md]
