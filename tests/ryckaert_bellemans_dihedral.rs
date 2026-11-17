//! Tests of the `ryckaert_bellemans_dihedral` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code
// Index loops over atoms and axes mirror the formulas they test.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

use ucmd_ff::dihedral_angle::dihedral_angle;
use ucmd_ff::errors::UcmdError;
use ucmd_ff::rb_dihedral_math::rb_energy;
use ucmd_ff::ryckaert_bellemans_dihedral::*;
use ucmd_ff::vec3::{add, norm};

#[test]
fn fd_forces_rb_dihedral() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];
    let r_k: [f64; 3] = [1.0, 1.0, 0.0];
    let r_l: [f64; 3] = [1.5, 1.3, 0.4];
    let c: [f64; 6] = [1.0, -0.5, 0.3, 0.1, 0.0, 0.0];
    let h: f64 = 1e-6;

    let (f_i, f_j, f_k, f_l) = ryckaert_bellemans_dihedral_forces(r_i, r_j, r_k, r_l, c).unwrap();

    let positions: [[f64; 3]; 4] = [r_i, r_j, r_k, r_l];
    let forces: [[f64; 3]; 4] = [f_i, f_j, f_k, f_l];

    for atom in 0..4 {
        for d in 0..3 {
            let mut pos_plus = positions;
            let mut pos_minus = positions;
            pos_plus[atom][d] += h;
            pos_minus[atom][d] -= h;

            let e_plus = ryckaert_bellemans_dihedral_energy(
                pos_plus[0],
                pos_plus[1],
                pos_plus[2],
                pos_plus[3],
                c,
            )
            .unwrap();
            let e_minus = ryckaert_bellemans_dihedral_energy(
                pos_minus[0],
                pos_minus[1],
                pos_minus[2],
                pos_minus[3],
                c,
            )
            .unwrap();

            let fd_force = -(e_plus - e_minus) / (2.0 * h);
            let got = forces[atom][d];

            assert!(
                (got - fd_force).abs() <= 1e-8,
                "atom {atom}, dir {d}: got {got}, expected {fd_force}"
            );
        }
    }
}

#[test]
fn fd_virial_rb_dihedral() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];
    let r_k: [f64; 3] = [1.0, 1.0, 0.0];
    let r_l: [f64; 3] = [1.5, 1.3, 0.4];
    let c: [f64; 6] = [1.0, -0.5, 0.3, 0.1, 0.0, 0.0];
    let eps: f64 = 1e-6;

    let w = ryckaert_bellemans_dihedral_virial(r_i, r_j, r_k, r_l, c).unwrap();

    let positions: [[f64; 3]; 4] = [r_i, r_j, r_k, r_l];

    for a in 0..3 {
        for b in 0..3 {
            // Deform: for each atom, add eps * r_atom[a] to coordinate b
            let mut pos_plus = positions;
            let mut pos_minus = positions;
            for atom in 0..4 {
                pos_plus[atom][b] += eps * positions[atom][a];
                pos_minus[atom][b] -= eps * positions[atom][a];
            }

            let e_plus = ryckaert_bellemans_dihedral_energy(
                pos_plus[0],
                pos_plus[1],
                pos_plus[2],
                pos_plus[3],
                c,
            )
            .unwrap();
            let e_minus = ryckaert_bellemans_dihedral_energy(
                pos_minus[0],
                pos_minus[1],
                pos_minus[2],
                pos_minus[3],
                c,
            )
            .unwrap();

            let fd_virial = -(e_plus - e_minus) / (2.0 * eps);
            let got = w[a][b];

            assert!(
                (got - fd_virial).abs() <= 1e-8,
                "W[{a}][{b}]: got {got}, expected {fd_virial}"
            );
        }
    }
}

#[test]
fn translation_invariance_energy() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];
    let r_k: [f64; 3] = [1.0, 1.0, 0.0];
    let r_l: [f64; 3] = [1.5, 1.3, 0.4];
    let c: [f64; 6] = [2.0, -1.0, 0.5, 0.2, 0.0, 0.0];
    let t: [f64; 3] = [10.0, -3.0, 7.0];

    let e_orig = ryckaert_bellemans_dihedral_energy(r_i, r_j, r_k, r_l, c).unwrap();

    let r_i_t = add(r_i, t);
    let r_j_t = add(r_j, t);
    let r_k_t = add(r_k, t);
    let r_l_t = add(r_l, t);

    let e_trans = ryckaert_bellemans_dihedral_energy(r_i_t, r_j_t, r_k_t, r_l_t, c).unwrap();

    assert!(
        (e_orig - e_trans).abs() <= 1e-12,
        "got {e_trans}, expected {e_orig}"
    );
}

#[test]
fn rotation_invariance_forces() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];
    let r_k: [f64; 3] = [1.0, 1.0, 0.0];
    let r_l: [f64; 3] = [1.5, 1.3, 0.4];
    let c: [f64; 6] = [1.0, -0.5, 0.3, 0.1, 0.0, 0.0];

    // 30 degrees about z-axis
    let angle: f64 = 30.0_f64.to_radians();
    let ca = angle.cos();
    let sa = angle.sin();

    // Rotation matrix R (rows are the new basis vectors)
    // R = [[ca, -sa, 0], [sa, ca, 0], [0, 0, 1]]
    let r: [[f64; 3]; 3] = [[ca, -sa, 0.0], [sa, ca, 0.0], [0.0, 0.0, 1.0]];

    // Apply rotation to all positions
    let positions: [[f64; 3]; 4] = [r_i, r_j, r_k, r_l];

    let mut rotated_positions: [[f64; 3]; 4] = [[0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3]];
    for atom in 0..4 {
        for d in 0..3 {
            rotated_positions[atom][d] = r[0][d] * positions[atom][0]
                + r[1][d] * positions[atom][1]
                + r[2][d] * positions[atom][2];
        }
    }

    let (f_i, f_j, f_k, f_l) = ryckaert_bellemans_dihedral_forces(r_i, r_j, r_k, r_l, c).unwrap();
    let (f_i_r, f_j_r, f_k_r, f_l_r) = ryckaert_bellemans_dihedral_forces(
        rotated_positions[0],
        rotated_positions[1],
        rotated_positions[2],
        rotated_positions[3],
        c,
    )
    .unwrap();

    let forces_orig: [[f64; 3]; 4] = [f_i, f_j, f_k, f_l];
    let forces_rot: [[f64; 3]; 4] = [f_i_r, f_j_r, f_k_r, f_l_r];

    for atom in 0..4 {
        for d in 0..3 {
            // Expected: R * F_original
            let expected = r[0][d] * forces_orig[atom][0]
                + r[1][d] * forces_orig[atom][1]
                + r[2][d] * forces_orig[atom][2];
            let got = forces_rot[atom][d];

            assert!(
                (got - expected).abs() <= 1e-10,
                "atom {atom}, dir {d}: got {got}, expected {expected}"
            );
        }
    }
}

#[test]
fn momentum_conservation() {
    let r_i: [f64; 3] = [0.3, 0.1, 0.2];
    let r_j: [f64; 3] = [1.3, 0.1, 0.2];
    let r_k: [f64; 3] = [1.3, 1.1, 0.2];
    let r_l: [f64; 3] = [1.8, 1.4, 0.6];
    let c: [f64; 6] = [1.5, -0.7, 0.4, 0.1, 0.0, 0.0];

    let (f_i, f_j, f_k, f_l) = ryckaert_bellemans_dihedral_forces(r_i, r_j, r_k, r_l, c).unwrap();

    let forces: [[f64; 3]; 4] = [f_i, f_j, f_k, f_l];

    // Sum of forces
    let mut sum: [f64; 3] = [0.0; 3];
    let mut total_norm: f64 = 0.0;
    for atom in 0..4 {
        for d in 0..3 {
            sum[d] += forces[atom][d];
        }
        total_norm += norm(forces[atom]);
    }

    let sum_norm = norm(sum);

    assert!(
        sum_norm <= 1e-10 * total_norm,
        "sum of forces {:?} has norm {sum_norm}, total norm {total_norm}",
        sum
    );
}

#[test]
fn opls_to_rb_energy_agreement() {
    let v: [f64; 4] = [1.0, 2.0, 3.0, 4.0];
    let c: [f64; 6] = opls_to_rb(v);

    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];
    let r_k: [f64; 3] = [1.0, 1.0, 0.0];

    let t_values: [f64; 11] = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 5.5];

    for &t in t_values.iter() {
        let r_l: [f64; 3] = [1.0 + t.cos(), 1.0 + t.sin(), 0.5 * t.sin()];

        let phi = dihedral_angle(r_i, r_j, r_k, r_l).unwrap();

        // OPLS-Fourier energy
        let e_opls = 0.5
            * (v[0] * (1.0 + phi.cos())
                + v[1] * (1.0 - (2.0 * phi).cos())
                + v[2] * (1.0 + (3.0 * phi).cos())
                + v[3] * (1.0 - (4.0 * phi).cos()));

        // RB energy
        let e_rb = rb_energy(phi, c);

        assert!(
            (e_opls - e_rb).abs() < 1e-10,
            "t={t}: E_opls={e_opls}, E_rb={e_rb}, diff={}",
            (e_opls - e_rb).abs()
        );
    }
}

#[test]
fn collinear_returns_error() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];
    let r_k: [f64; 3] = [2.0, 0.0, 0.0];
    let r_l: [f64; 3] = [3.0, 0.0, 0.0];
    let c: [f64; 6] = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0];

    let e_result = ryckaert_bellemans_dihedral_energy(r_i, r_j, r_k, r_l, c);
    assert!(
        matches!(e_result, Err(UcmdError::CoincidentAtoms)),
        "energy: expected Err(CoincidentAtoms), got {e_result:?}"
    );

    let f_result = ryckaert_bellemans_dihedral_forces(r_i, r_j, r_k, r_l, c);
    assert!(
        matches!(f_result, Err(UcmdError::CoincidentAtoms)),
        "forces: expected Err(CoincidentAtoms), got {f_result:?}"
    );

    let w_result = ryckaert_bellemans_dihedral_virial(r_i, r_j, r_k, r_l, c);
    assert!(
        matches!(w_result, Err(UcmdError::CoincidentAtoms)),
        "virial: expected Err(CoincidentAtoms), got {w_result:?}"
    );
}

#[test]
fn zero_coefficients_give_zero() {
    let r_i: [f64; 3] = [0.0, 0.0, 0.0];
    let r_j: [f64; 3] = [1.0, 0.0, 0.0];
    let r_k: [f64; 3] = [1.0, 1.0, 0.0];
    let r_l: [f64; 3] = [1.5, 1.3, 0.4];
    let c: [f64; 6] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

    let e = ryckaert_bellemans_dihedral_energy(r_i, r_j, r_k, r_l, c).unwrap();
    assert!(e.abs() <= 1e-15, "energy: got {e}, expected 0.0");

    let (f_i, f_j, f_k, f_l) = ryckaert_bellemans_dihedral_forces(r_i, r_j, r_k, r_l, c).unwrap();

    let forces: [[f64; 3]; 4] = [f_i, f_j, f_k, f_l];
    for atom in 0..4 {
        for d in 0..3 {
            assert!(
                forces[atom][d].abs() <= 1e-15,
                "force atom {atom} dir {d}: got {}, expected 0.0",
                forces[atom][d]
            );
        }
    }

    let w = ryckaert_bellemans_dihedral_virial(r_i, r_j, r_k, r_l, c).unwrap();
    for a in 0..3 {
        for b in 0..3 {
            assert!(
                w[a][b].abs() <= 1e-15,
                "virial W[{a}][{b}]: got {}, expected 0.0",
                w[a][b]
            );
        }
    }
}
