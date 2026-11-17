//! Tests of the Ewald sum (M4: R4.1 parts and exclusions, R4.2 forces and virial, R4.3 validation).
//!
//! Test system: a triclinic cell in reduced form, eight point charges with a neutral total,
//! positions partly outside the cell (unwrapped). Cell, charges and positions are arbitrary test
//! geometry. Atoms 0-2 form one "molecule" (bonds 0-1, 0-2), atoms 3-6 a chain (bonds 3-4, 4-5,
//! 5-6, so 3-6 is a 1-4 pair); their pairs at bond distance 1 and 2 are excluded.
//!
//! Independent reference: a classic Ewald sum written separately in Python (its own
//! real-space image loop, a box of reciprocal vectors |h|,|k|,|l| <= k_max instead of a sphere,
//! math.erfc, CODATA 2018 k_e), at settings where every part is converged (raising the cutoffs
//! changes no part by more than 3e-12 kJ/mol); forces and virial by central differences of
//! that sum (h = 1e-5 nm, strain 1e-6).
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code
// Index loops over atoms and axes mirror the formulas they test.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

use ucmd_ff::cell::{Cell, minimum_image};
use ucmd_ff::constants::COULOMB_KJ_MOL_NM_E2;
use ucmd_ff::erf::erfc;
use ucmd_ff::errors::UcmdError;
use ucmd_ff::ewald::{
    EwaldParameters, EwaldPart, ScaledPair, ewald, ewald_correction, ewald_self,
    scaled_pairs_from_bonds,
};
use ucmd_ff::ewald_real::ewald_real;
use ucmd_ff::ewald_reciprocal::ewald_reciprocal;
use ucmd_ff::fd_test::{numerical_forces, numerical_virial_periodic};

fn cell() -> Cell {
    Cell {
        ax: 1.1,
        bx: 0.3,
        by: 1.0,
        cx: -0.2,
        cy: 0.25,
        cz: 0.9,
    }
}

const Q: [f64; 8] = [0.8, -0.4, -0.4, 0.5, -0.5, 1.0, -0.6, -0.4];

fn positions() -> Vec<[f64; 3]> {
    vec![
        [0.10, 0.20, 0.15],
        [0.19, 0.24, 0.17],
        [0.05, 0.28, 0.12],
        [0.70, 0.55, 0.40],
        [0.78, 0.62, 0.47],
        [0.88, 0.58, 0.55],
        [0.97, 0.66, 0.60],
        [1.35, -0.15, 0.95],
    ]
}

fn excluded() -> Vec<ScaledPair> {
    [
        (0, 1),
        (0, 2),
        (1, 2),
        (3, 4),
        (3, 5),
        (4, 5),
        (4, 6),
        (5, 6),
    ]
    .iter()
    .map(|&(i, j)| ScaledPair { i, j, scale: 0.0 })
    .collect()
}

/// Converged settings (alpha = 6 nm^-1, r_c = 0.9 nm: erfc(5.4) = 1.4e-14; k_c = 72 nm^-1:
/// exp(-k_c^2/(4 alpha^2)) = 2.3e-16), the same alpha as the reference.
fn reference_params() -> EwaldParameters {
    EwaldParameters {
        alpha_per_nm: 6.0,
        real_cutoff_nm: 0.9,
        recip_cutoff_per_nm: 72.0,
    }
}

fn max_abs_force(f: &[[f64; 3]]) -> f64 {
    f.iter().flatten().fold(0.0_f64, |m, x| m.max(x.abs()))
}

fn sum_abs_force(f: &[[f64; 3]]) -> f64 {
    f.iter()
        .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
        .sum()
}

/// Strain components that keep a lower-triangular cell lower-triangular (spec/units.md §4).
const LOWER: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (1, 0), (2, 0), (2, 1)];

// ---------------------------------------------------------------- independent reference

#[test]
fn parts_match_independent_ewald_sum() {
    let r = ewald(&positions(), &Q, &cell(), &reference_params(), &excluded()).unwrap();
    // Reference energies (kJ/mol) of the Python Ewald sum, alpha = 6, r_c = 0.9, k_max = 13.
    let cases = [
        ("real", r.real.energy_kjmol, -0.4792310381060577),
        ("reciprocal", r.reciprocal.energy_kjmol, 193.18334058445598),
        ("self", r.self_term.energy_kjmol, -1401.5405712326065),
        (
            "correction (excluded pairs)",
            r.correction.energy_kjmol,
            1019.5804007234583,
        ),
        ("total", r.total.energy_kjmol, -189.25606096279824),
    ];
    for (name, got, want) in cases {
        assert!(
            (got - want).abs() < 1e-9,
            "{name}: {got} vs reference {want}"
        );
    }
}

#[test]
fn total_forces_and_virial_match_independent_reference() {
    // Reference: central differences of the Python sum's converged Ewald total (alpha = 4.5,
    // r_c = 1.25, k_max = 11); the converged total does not depend on these settings.
    let want_f: [[f64; 3]; 8] = [
        [44.74227341404457, -23.46796078285251, -15.222878734277854],
        [-13.247942115413023, 10.236807473518184, 23.73869437946041],
        [-19.818914230995688, -20.43385649699303, 19.067728696597896],
        [265.31395241988776, 120.8189029711093, 25.354787000253506],
        [36.366541679910824, 18.368190781359317, 94.90477002884744],
        [-115.86969130945589, 16.741265535813454, -98.79581998433194],
        [-172.5859252559303, -81.08733791658551, -178.16984570941938],
        [-24.900294629048854, -41.17601159663309, 129.12256436123926],
    ];
    let want_w = [
        ((0, 0), -66.78038214147364),
        ((1, 1), -24.92354828120824),
        ((2, 2), -97.55213058326717),
        ((1, 0), -26.499758661202577),
        ((2, 0), -41.559121967793544),
        ((2, 1), -12.992502476549816),
    ];
    let r = ewald(&positions(), &Q, &cell(), &reference_params(), &excluded()).unwrap();
    for i in 0..8 {
        for c in 0..3 {
            let (got, want) = (r.total.forces[i][c], want_f[i][c]);
            assert!(
                (got - want).abs() < 1e-5,
                "F[{i}][{c}]: {got} vs reference {want}"
            );
        }
    }
    for ((a, b), want) in want_w {
        let got = r.total.virial[a][b];
        assert!(
            (got - want).abs() < 1e-5,
            "W[{a}][{b}]: {got} vs reference {want}"
        );
    }
}

// ---------------------------------------------------------------- R4.3 validation

/// Rock salt in its cubic conventional cell (4 Na+ and 4 Cl-, nearest-neighbour distance a/2).
fn rock_salt(a: f64) -> (Vec<[f64; 3]>, Vec<f64>, Cell) {
    let na = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ];
    let cl = [
        [0.5, 0.0, 0.0],
        [0.0, 0.5, 0.0],
        [0.0, 0.0, 0.5],
        [0.5, 0.5, 0.5],
    ];
    let mut pos = Vec::new();
    let mut q = Vec::new();
    for s in na {
        pos.push([a * s[0], a * s[1], a * s[2]]);
        q.push(1.0);
    }
    for s in cl {
        pos.push([a * s[0], a * s[1], a * s[2]]);
        q.push(-1.0);
    }
    (
        pos,
        q,
        Cell {
            ax: a,
            bx: 0.0,
            by: a,
            cx: 0.0,
            cy: 0.0,
            cz: a,
        },
    )
}

#[test]
fn madelung_constant_of_rock_salt() {
    // R4.3: M = 1.7475645946 to 1e-8 relative. E = -N_pairs M k_e / r0 with 4 ion pairs and
    // r0 = a/2, so M = -E a / (8 k_e). a = 0.564 nm is arbitrary (M does not depend on it).
    let a = 0.564;
    let (pos, q, cell) = rock_salt(a);
    // Two different splittings; the first cutoff exceeds half the cell (complete image search).
    for (r_c, tol) in [(0.8, 1e-15), (0.27, 1e-15)] {
        let params = EwaldParameters::from_tolerance(r_c, tol).unwrap();
        let r = ewald(&pos, &q, &cell, &params, &[]).unwrap();
        let m = -r.total.energy_kjmol * a / (8.0 * COULOMB_KJ_MOL_NM_E2);
        let rel = (m - 1.7475645946) / 1.7475645946;
        assert!(
            rel.abs() < 1e-8,
            "r_c = {r_c}: Madelung constant {m}, relative error {rel:e}"
        );
        // A perfect crystal has zero force on every ion.
        assert!(
            max_abs_force(&r.total.forces) < 1e-9,
            "forces {:?}",
            r.total.forces
        );
    }
}

#[test]
fn total_is_independent_of_alpha() {
    let pos = positions();
    let base = ewald(&pos, &Q, &cell(), &reference_params(), &excluded())
        .unwrap()
        .total;
    // Cutoffs from 0.6 nm (below h_min/2) to 1.6 nm (well above), alpha from about 9.4 to 3.5.
    for r_c in [0.6, 0.75, 1.2, 1.6] {
        let params = EwaldParameters::from_tolerance(r_c, 1e-15).unwrap();
        let t = ewald(&pos, &Q, &cell(), &params, &excluded())
            .unwrap()
            .total;
        let de = t.energy_kjmol - base.energy_kjmol;
        assert!(
            de.abs() < 1e-9 * base.energy_kjmol.abs(),
            "alpha = {}: dE = {de:e}",
            params.alpha_per_nm
        );
        let scale = max_abs_force(&base.forces);
        for i in 0..pos.len() {
            for c in 0..3 {
                let d = t.forces[i][c] - base.forces[i][c];
                assert!(
                    d.abs() < 1e-9 * scale,
                    "alpha = {}: dF[{i}][{c}] = {d:e}",
                    params.alpha_per_nm
                );
            }
        }
        for a in 0..3 {
            for b in 0..3 {
                let d = t.virial[a][b] - base.virial[a][b];
                assert!(
                    d.abs() < 1e-9 * scale,
                    "alpha = {}: dW[{a}][{b}] = {d:e}",
                    params.alpha_per_nm
                );
            }
        }
    }
}

/// Direct Coulomb energy k_e q_i q_j / r of the minimum image of a pair.
fn direct(pos: &[[f64; 3]], i: usize, j: usize) -> f64 {
    let d = minimum_image(pos[i], pos[j], &cell()).unwrap();
    COULOMB_KJ_MOL_NM_E2 * Q[i] * Q[j] / (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

#[test]
fn excluding_a_pair_removes_its_direct_coulomb_energy() {
    // R4.3: excluding a pair changes the converged total by exactly -k_e q_i q_j / r_ij (minimum
    // image), for any alpha. Pair 0-7: its unwrapped distance (1.52 nm) is not the minimum image.
    let pos = positions();
    for r_c in [0.7, 1.3] {
        let params = EwaldParameters::from_tolerance(r_c, 1e-15).unwrap();
        let without = ewald(&pos, &Q, &cell(), &params, &excluded())
            .unwrap()
            .total
            .energy_kjmol;
        for (i, j) in [(0, 7), (3, 6), (2, 5)] {
            let mut more = excluded();
            more.push(ScaledPair { i, j, scale: 0.0 });
            let with = ewald(&pos, &Q, &cell(), &params, &more)
                .unwrap()
                .total
                .energy_kjmol;
            let expected = -direct(&pos, i, j);
            let got = with - without;
            assert!(
                (got - expected).abs() < 1e-10 * expected.abs(),
                "alpha = {}, pair {i}-{j}: dE = {got}, expected {expected}",
                params.alpha_per_nm
            );
        }
    }
}

#[test]
fn scaled_pair_changes_energy_by_scale_minus_one_times_direct() {
    let pos = positions();
    let params = EwaldParameters::from_tolerance(1.0, 1e-15).unwrap();
    let without = ewald(&pos, &Q, &cell(), &params, &excluded())
        .unwrap()
        .total
        .energy_kjmol;
    for s in [0.5, 0.8333333333333334] {
        let mut more = excluded();
        more.push(ScaledPair {
            i: 3,
            j: 6,
            scale: s,
        });
        let with = ewald(&pos, &Q, &cell(), &params, &more)
            .unwrap()
            .total
            .energy_kjmol;
        let expected = (s - 1.0) * direct(&pos, 3, 6);
        let got = with - without;
        assert!(
            (got - expected).abs() < 1e-10 * expected.abs(),
            "s = {s}: dE = {got}, expected {expected}"
        );
    }
    // A scale of exactly 1 is the same as no entry.
    let mut one = excluded();
    one.push(ScaledPair {
        i: 3,
        j: 6,
        scale: 1.0,
    });
    let with = ewald(&pos, &Q, &cell(), &params, &one)
        .unwrap()
        .total
        .energy_kjmol;
    assert!((with - without).abs() < 1e-11 * without.abs());
}

// ---------------------------------------------------------------- R4.2 forces and virial

type Part = fn(&[[f64; 3]], &Cell) -> Result<EwaldPart, UcmdError>;

fn fd_params() -> EwaldParameters {
    EwaldParameters::from_tolerance(0.95, 1e-15).unwrap()
}

fn part_real(p: &[[f64; 3]], c: &Cell) -> Result<EwaldPart, UcmdError> {
    ewald_real(p, &Q, c, &fd_params(), &excluded())
}
/// Real space with a cutoff above every lattice vector (|a| = 1.1, |b| = 1.04, |c| = 0.96 nm), so
/// that the interaction of each atom with its own images enters.
fn part_real_long(p: &[[f64; 3]], c: &Cell) -> Result<EwaldPart, UcmdError> {
    let params = EwaldParameters::from_tolerance(1.5, 1e-15)?;
    ewald_real(p, &Q, c, &params, &excluded())
}
fn part_reciprocal(p: &[[f64; 3]], c: &Cell) -> Result<EwaldPart, UcmdError> {
    ewald_reciprocal(p, &Q, c, &fd_params())
}
fn part_correction(p: &[[f64; 3]], c: &Cell) -> Result<EwaldPart, UcmdError> {
    let mut s = excluded();
    s.push(ScaledPair {
        i: 3,
        j: 6,
        scale: 0.5,
    });
    ewald_correction(p, &Q, c, &fd_params(), &s)
}
fn part_total(p: &[[f64; 3]], c: &Cell) -> Result<EwaldPart, UcmdError> {
    let mut s = excluded();
    s.push(ScaledPair {
        i: 3,
        j: 6,
        scale: 0.5,
    });
    ewald(p, &Q, c, &fd_params(), &s).map(|r| r.total)
}

const PARTS: [(&str, Part); 5] = [
    ("real", part_real),
    ("real, cutoff above the lattice vectors", part_real_long),
    ("reciprocal", part_reciprocal),
    ("correction", part_correction),
    ("total", part_total),
];

#[test]
fn forces_of_every_part_match_finite_differences() {
    let pos = positions();
    let cell = cell();
    for (name, part) in PARTS {
        let analytic = part(&pos, &cell).unwrap().forces;
        let energy = |p: &[[f64; 3]]| part(p, &cell).map(|r| r.energy_kjmol);
        let numeric = numerical_forces(&energy, &pos, 1e-5).unwrap();
        let scale = max_abs_force(&analytic);
        assert!(scale > 1.0, "{name}: forces too small to test ({scale})");
        for i in 0..pos.len() {
            for c in 0..3 {
                let d = analytic[i][c] - numeric[i][c];
                assert!(
                    d.abs() < 1e-7 * scale,
                    "{name}: F[{i}][{c}] analytic {} numerical {}",
                    analytic[i][c],
                    numeric[i][c]
                );
            }
        }
    }
}

#[test]
fn virial_of_every_part_matches_finite_differences() {
    let pos = positions();
    let cell = cell();
    for (name, part) in PARTS {
        let analytic = part(&pos, &cell).unwrap().virial;
        let energy = |p: &[[f64; 3]], c: &Cell| part(p, c).map(|r| r.energy_kjmol);
        let scale = analytic
            .iter()
            .flatten()
            .fold(0.0_f64, |m, x| m.max(x.abs()));
        assert!(scale > 1.0, "{name}: virial too small to test ({scale})");
        for (a, b) in LOWER {
            let numeric = numerical_virial_periodic(&energy, &pos, &cell, a, b, 1e-6).unwrap();
            let d = analytic[a][b] - numeric;
            assert!(
                d.abs() < 1e-7 * scale,
                "{name}: W[{a}][{b}] analytic {} numerical {numeric}",
                analytic[a][b]
            );
        }
        // Rotation invariance: the virial is symmetric.
        for a in 0..3 {
            for b in 0..a {
                let d = analytic[a][b] - analytic[b][a];
                assert!(
                    d.abs() < 1e-11 * scale,
                    "{name}: W not symmetric at ({a},{b}): {d:e}"
                );
            }
        }
    }
}

#[test]
fn self_term_has_no_forces_or_virial_and_the_closed_form_energy() {
    let params = reference_params();
    let s = ewald_self(&positions(), &Q, &cell(), &params).unwrap();
    let sum_q2: f64 = Q.iter().map(|q| q * q).sum();
    let expected =
        -COULOMB_KJ_MOL_NM_E2 * params.alpha_per_nm / std::f64::consts::PI.sqrt() * sum_q2;
    assert!((s.energy_kjmol - expected).abs() < 1e-12 * expected.abs());
    assert!(s.forces.iter().flatten().all(|&f| f == 0.0));
    assert!(s.virial.iter().flatten().all(|&w| w == 0.0));
}

#[test]
fn momentum_is_conserved() {
    let r = part_total(&positions(), &cell()).unwrap();
    let mut sum = [0.0_f64; 3];
    for f in &r.forces {
        for c in 0..3 {
            sum[c] += f[c];
        }
    }
    let total = sum_abs_force(&r.forces);
    for c in 0..3 {
        assert!(
            sum[c].abs() <= 1e-10 * total,
            "sum F[{c}] = {:e} (sum |F| = {total})",
            sum[c]
        );
    }
}

#[test]
fn unwrapped_positions_give_the_same_result() {
    // Moving atoms by lattice vectors changes nothing (positions are unwrapped, spec/units.md §3).
    let pos = positions();
    let c = cell();
    let mut moved = pos.clone();
    let shift = |n: [f64; 3]| {
        [
            n[0] * c.ax + n[1] * c.bx + n[2] * c.cx,
            n[1] * c.by + n[2] * c.cy,
            n[2] * c.cz,
        ]
    };
    for (i, n) in [
        (3, [1.0, 2.0, -1.0]),
        (7, [-3.0, 0.0, 4.0]),
        (1, [0.0, -5.0, 2.0]),
    ] {
        let s = shift(n);
        for k in 0..3 {
            moved[i][k] += s[k];
        }
    }
    let a = part_total(&pos, &c).unwrap();
    let b = part_total(&moved, &c).unwrap();
    assert!((a.energy_kjmol - b.energy_kjmol).abs() < 1e-10 * a.energy_kjmol.abs());
    let scale = max_abs_force(&a.forces);
    for i in 0..pos.len() {
        for k in 0..3 {
            assert!(
                (a.forces[i][k] - b.forces[i][k]).abs() < 1e-9 * scale,
                "F[{i}][{k}]"
            );
        }
    }
    for x in 0..3 {
        for y in 0..3 {
            assert!(
                (a.virial[x][y] - b.virial[x][y]).abs() < 1e-9 * scale,
                "W[{x}][{y}]"
            );
        }
    }
}

// ---------------------------------------------------------------- parameters and errors

#[test]
fn parameters_from_tolerance() {
    let p = EwaldParameters::from_tolerance(0.9, 1e-12).unwrap();
    let e = erfc(p.alpha_per_nm * p.real_cutoff_nm);
    assert!((e - 1e-12).abs() < 1e-12 * 1e-9, "erfc(alpha r_c) = {e:e}");
    let k = 2.0 * p.alpha_per_nm * (-(1e-12_f64).ln()).sqrt();
    assert!((p.recip_cutoff_per_nm - k).abs() < 1e-12 * k);
    for (r_c, tol) in [
        (0.0, 1e-10),
        (-1.0, 1e-10),
        (f64::NAN, 1e-10),
        (1.0, 0.0),
        (1.0, 1.0),
        (1.0, -0.5),
    ] {
        assert_eq!(
            EwaldParameters::from_tolerance(r_c, tol),
            Err(UcmdError::InvalidParameter),
            "{r_c} {tol}"
        );
    }
}

#[test]
fn net_charge_beyond_1e_6_is_an_error() {
    let pos = positions();
    let params = reference_params();
    let mut q = Q.to_vec();
    q[7] += 2e-6;
    assert_eq!(
        ewald(&pos, &q, &cell(), &params, &[]).err(),
        Some(UcmdError::NonNeutralSystem)
    );
    assert_eq!(
        ewald_real(&pos, &q, &cell(), &params, &[]).err(),
        Some(UcmdError::NonNeutralSystem)
    );
    assert_eq!(
        ewald_reciprocal(&pos, &q, &cell(), &params).err(),
        Some(UcmdError::NonNeutralSystem)
    );
    let mut q = Q.to_vec();
    q[7] += 5e-7;
    assert!(ewald(&pos, &q, &cell(), &params, &[]).is_ok());
}

#[test]
fn invalid_parameters_are_errors() {
    let pos = positions();
    let good = reference_params();
    let bad = [
        EwaldParameters {
            alpha_per_nm: 0.0,
            ..good
        },
        EwaldParameters {
            alpha_per_nm: -1.0,
            ..good
        },
        EwaldParameters {
            real_cutoff_nm: 0.0,
            ..good
        },
        EwaldParameters {
            recip_cutoff_per_nm: f64::NAN,
            ..good
        },
        EwaldParameters {
            real_cutoff_nm: f64::INFINITY,
            ..good
        },
    ];
    for p in bad {
        assert_eq!(
            ewald(&pos, &Q, &cell(), &p, &[]).err(),
            Some(UcmdError::InvalidParameter),
            "{p:?}"
        );
    }
}

#[test]
fn invalid_inputs_are_errors() {
    let pos = positions();
    let p = reference_params();
    let c = cell();
    // Lengths differ.
    assert_eq!(
        ewald(&pos[..7], &Q, &c, &p, &[]).err(),
        Some(UcmdError::IndexOutOfRange)
    );
    // Non-finite position and charge.
    let mut bad = pos.clone();
    bad[2][1] = f64::NAN;
    assert_eq!(
        ewald(&bad, &Q, &c, &p, &[]).err(),
        Some(UcmdError::NonFiniteCoordinate)
    );
    let mut q = Q.to_vec();
    q[0] = f64::INFINITY;
    assert_eq!(
        ewald(&pos, &q, &c, &p, &[]).err(),
        Some(UcmdError::InvalidParameter)
    );
    // Invalid cell.
    let flat = Cell { by: 0.0, ..c };
    assert_eq!(
        ewald(&pos, &Q, &flat, &p, &[]).err(),
        Some(UcmdError::InvalidCell)
    );
    // Two charged atoms at the same place (here modulo a lattice vector).
    let mut same = pos.clone();
    same[1] = [pos[0][0] + c.ax, pos[0][1], pos[0][2]];
    assert_eq!(
        ewald(&same, &Q, &c, &p, &[]).err(),
        Some(UcmdError::CoincidentAtoms)
    );
    // Bad scaled pairs.
    let cases = [
        (
            vec![ScaledPair {
                i: 0,
                j: 8,
                scale: 0.0,
            }],
            UcmdError::IndexOutOfRange,
        ),
        (
            vec![ScaledPair {
                i: 2,
                j: 2,
                scale: 0.0,
            }],
            UcmdError::InvalidParameter,
        ),
        (
            vec![ScaledPair {
                i: 0,
                j: 1,
                scale: 1.5,
            }],
            UcmdError::InvalidParameter,
        ),
        (
            vec![ScaledPair {
                i: 0,
                j: 1,
                scale: f64::NAN,
            }],
            UcmdError::InvalidParameter,
        ),
        (
            vec![
                ScaledPair {
                    i: 0,
                    j: 1,
                    scale: 0.0,
                },
                ScaledPair {
                    i: 1,
                    j: 0,
                    scale: 0.0,
                },
            ],
            UcmdError::InvalidParameter,
        ),
    ];
    for (pairs, err) in cases {
        assert_eq!(
            ewald(&pos, &Q, &c, &p, &pairs).err(),
            Some(err),
            "{pairs:?}"
        );
    }
}

#[test]
fn scaled_pairs_from_a_bond_graph() {
    // Chain 0-1-2-3-4 plus an isolated atom 5.
    let bonds = [(0, 1), (1, 2), (2, 3), (3, 4)];
    let pairs = scaled_pairs_from_bonds(6, &bonds, 0.5).unwrap();
    let got: Vec<(usize, usize, f64)> = pairs.iter().map(|p| (p.i, p.j, p.scale)).collect();
    let expected = vec![
        (0, 1, 0.0),
        (0, 2, 0.0),
        (0, 3, 0.5),
        (1, 2, 0.0),
        (1, 3, 0.0),
        (1, 4, 0.5),
        (2, 3, 0.0),
        (2, 4, 0.0),
        (3, 4, 0.0),
    ];
    assert_eq!(got, expected);
    // A 1-4 scale of 1 adds no 1-4 pair.
    let full = scaled_pairs_from_bonds(6, &bonds, 1.0).unwrap();
    assert_eq!(full.len(), 7);
    assert!(full.iter().all(|p| p.scale == 0.0));
    assert_eq!(
        scaled_pairs_from_bonds(6, &[(0, 6)], 0.5).err(),
        Some(UcmdError::IndexOutOfRange)
    );
    assert_eq!(
        scaled_pairs_from_bonds(6, &bonds, 1.2).err(),
        Some(UcmdError::InvalidParameter)
    );
}
