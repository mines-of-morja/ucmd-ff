//! Tests of the `nonbonded_periodic` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code
// Index loops over atoms and axes mirror the formulas they test.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

use ucmd_ff::cell::*;
use ucmd_ff::combination_rule::*;
use ucmd_ff::constants::*;
use ucmd_ff::errors::*;
use ucmd_ff::fd_test::*;
use ucmd_ff::lj_pair::*;
use ucmd_ff::min_image_disp::*;
use ucmd_ff::nonbonded_nonperiodic::*;
use ucmd_ff::nonbonded_periodic::*;

fn make_cell_cubic(size: f64) -> Cell {
    Cell {
        ax: size,
        bx: 0.0,
        by: size,
        cx: 0.0,
        cy: 0.0,
        cz: size,
    }
}

fn make_input(
    positions: Vec<[f64; 3]>,
    charges: Vec<f64>,
    lj_sigma: Vec<f64>,
    lj_eps: Vec<f64>,
    bonds: Vec<(usize, usize)>,
    rule: CombinationRule,
    lj_scale_14: f64,
    coul_scale_14: f64,
) -> NonbondedInput {
    let n = positions.len();
    NonbondedInput {
        positions,
        charges,
        lj_sigma,
        lj_eps,
        lj_scale_14,
        coul_scale_14,
        bonds,
        rule,
        types: vec![0; n],
    }
}

#[test]
fn test_lj_minimum_energy() {
    let sigma: f64 = 0.34;
    let eps: f64 = 1.0;
    let r_min: f64 = (2.0_f64).powf(1.0 / 6.0) * sigma;
    let cell = make_cell_cubic(2.0);
    let input = make_input(
        vec![[0.0, 0.0, 0.0], [r_min, 0.0, 0.0]],
        vec![0.0, 0.0],
        vec![sigma, sigma],
        vec![eps, eps],
        vec![],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let result = evaluate_periodic(&input, &cell, 0.9).unwrap();
    let expected: f64 = -eps;
    assert!(
        (result.energy_lj_kjmol - expected).abs() <= 1e-12,
        "got {}, expected {}",
        result.energy_lj_kjmol,
        expected
    );
}

#[test]
fn test_cutoff_excludes_lj_keeps_coulomb() {
    let sigma: f64 = 0.34;
    let eps: f64 = 1.0;
    let r: f64 = 1.5;
    let cell = make_cell_cubic(3.0);
    let input = make_input(
        vec![[0.0, 0.0, 0.0], [r, 0.0, 0.0]],
        vec![1.0, -1.0],
        vec![sigma, sigma],
        vec![eps, eps],
        vec![],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let result = evaluate_periodic(&input, &cell, 1.0).unwrap();
    assert!(
        result.energy_lj_kjmol.abs() <= 1e-12,
        "got {}, expected 0.0",
        result.energy_lj_kjmol
    );
    let expected_coulomb: f64 = -COULOMB_KJ_MOL_NM_E2 * 1.0 * 1.0 / r;
    assert!(
        (result.energy_coulomb_kjmol - expected_coulomb).abs() <= 1e-12,
        "got {}, expected {}",
        result.energy_coulomb_kjmol,
        expected_coulomb
    );
}

#[test]
fn test_cutoff_validation_error() {
    let cell = make_cell_cubic(2.0);
    let input = make_input(
        vec![[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]],
        vec![0.0, 0.0],
        vec![0.34, 0.34],
        vec![1.0, 1.0],
        vec![],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let result = evaluate_periodic(&input, &cell, 1.1);
    assert!(
        matches!(result, Err(UcmdError::InvalidCell)),
        "expected Err(InvalidCell), got {:?}",
        result
    );
}

#[test]
fn test_minimum_image_wrapping() {
    let sigma: f64 = 0.34;
    let eps: f64 = 1.0;
    let cell = make_cell_cubic(2.0);
    let input = make_input(
        vec![[1.9, 0.0, 0.0], [0.1, 0.0, 0.0]],
        vec![0.0, 0.0],
        vec![sigma, sigma],
        vec![eps, eps],
        vec![],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let result = evaluate_periodic(&input, &cell, 0.99).unwrap();
    let expected = lj_pair_energy(sigma, eps, [0.2, 0.0, 0.0], [0.0, 0.0, 0.0]).unwrap();
    assert!(
        (result.energy_lj_kjmol - expected).abs() <= 1e-10,
        "got {}, expected {}",
        result.energy_lj_kjmol,
        expected
    );
}

#[test]
fn test_finite_difference_forces() {
    let sigma: f64 = 0.34;
    let eps: f64 = 1.0;
    let positions = vec![[0.1, 0.2, 0.3], [1.5, 0.4, 0.6], [0.8, 1.1, 0.1]];
    let charges = vec![1.0, -1.0, 0.5];
    let cell = make_cell_cubic(2.0);
    let input = make_input(
        positions.clone(),
        charges.clone(),
        vec![sigma, sigma, sigma],
        vec![eps, eps, eps],
        vec![(0, 1)],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let result = evaluate_periodic(&input, &cell, 0.49).unwrap();
    let energy_fn = |pos: &[[f64; 3]]| -> Result<f64, UcmdError> {
        let mut inp = input.clone();
        inp.positions = pos.to_vec();
        evaluate_periodic(&inp, &cell, 0.49).map(|r| r.energy_lj_kjmol + r.energy_coulomb_kjmol)
    };
    let numerical = numerical_forces(&energy_fn, &positions, 1e-6).unwrap();
    let max_diff = result
        .forces
        .iter()
        .zip(numerical.iter())
        .map(|(a, n)| {
            [
                (a[0] - n[0]).abs(),
                (a[1] - n[1]).abs(),
                (a[2] - n[2]).abs(),
            ]
            .into_iter()
            .fold(0.0_f64, |m, d| m.max(d))
        })
        .fold(0.0_f64, |m, d| m.max(d));
    assert!(
        max_diff <= 1e-5,
        "max abs diff {} exceeds tolerance 1e-5",
        max_diff
    );
}

#[test]
fn test_finite_difference_virial() {
    let sigma: f64 = 0.34;
    let eps: f64 = 1.0;
    let positions = vec![[0.1, 0.2, 0.3], [1.5, 0.4, 0.6], [0.8, 1.1, 0.1]];
    let charges = vec![1.0, -1.0, 0.5];
    let cell = make_cell_cubic(2.0);
    let input = make_input(
        positions.clone(),
        charges.clone(),
        vec![sigma, sigma, sigma],
        vec![eps, eps, eps],
        vec![(0, 1)],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let result = evaluate_periodic(&input, &cell, 0.49).unwrap();
    let energy_fn = |pos: &[[f64; 3]], c: &Cell| -> Result<f64, UcmdError> {
        let mut inp = input.clone();
        inp.positions = pos.to_vec();
        evaluate_periodic(&inp, c, 0.49).map(|r| r.energy_lj_kjmol + r.energy_coulomb_kjmol)
    };
    for (a, b) in [(0, 0), (1, 0), (2, 0), (2, 1), (1, 1), (2, 2)] {
        let numerical =
            numerical_virial_periodic(&energy_fn, &positions, &cell, a, b, 1e-6).unwrap();
        let analytic = result.virial[a][b];
        assert!(
            (analytic - numerical).abs() <= 1e-5,
            "virial[{}][{}]: got {}, expected {}",
            a,
            b,
            analytic,
            numerical
        );
    }
}

#[test]
fn test_translation_invariance() {
    let sigma: f64 = 0.34;
    let eps: f64 = 1.0;
    let positions = vec![[0.1, 0.2, 0.3], [1.5, 0.4, 0.6], [0.8, 1.2, 0.1]];
    let charges = vec![1.0, -1.0, 0.5];
    let cell = make_cell_cubic(2.0);
    let input1 = make_input(
        positions.clone(),
        charges.clone(),
        vec![sigma, sigma, sigma],
        vec![eps, eps, eps],
        vec![(0, 1)],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let shifted: Vec<[f64; 3]> = positions.iter().map(|p| [p[0] + 2.0, p[1], p[2]]).collect();
    let input2 = make_input(
        shifted,
        charges.clone(),
        vec![sigma, sigma, sigma],
        vec![eps, eps, eps],
        vec![(0, 1)],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let r1 = evaluate_periodic(&input1, &cell, 0.49).unwrap();
    let r2 = evaluate_periodic(&input2, &cell, 0.49).unwrap();
    let tol: f64 = 1e-12;
    assert!(
        (r1.energy_lj_kjmol - r2.energy_lj_kjmol).abs() <= tol,
        "LJ energy: got {}, expected {}",
        r1.energy_lj_kjmol,
        r2.energy_lj_kjmol
    );
    assert!(
        (r1.energy_coulomb_kjmol - r2.energy_coulomb_kjmol).abs() <= tol,
        "Coulomb energy: got {}, expected {}",
        r1.energy_coulomb_kjmol,
        r2.energy_coulomb_kjmol
    );
    for k in 0..3 {
        for c in 0..3 {
            assert!(
                (r1.forces[k][c] - r2.forces[k][c]).abs() <= tol,
                "force[{}][{}]: got {}, expected {}",
                k,
                c,
                r1.forces[k][c],
                r2.forces[k][c]
            );
        }
    }
    for a in 0..3 {
        for b in 0..3 {
            assert!(
                (r1.virial[a][b] - r2.virial[a][b]).abs() <= tol,
                "virial[{}][{}]: got {}, expected {}",
                a,
                b,
                r1.virial[a][b],
                r2.virial[a][b]
            );
        }
    }
}

#[test]
fn test_sum_of_forces_vanishes() {
    let sigma: f64 = 0.34;
    let eps: f64 = 1.0;
    let positions = vec![
        [0.1, 0.2, 0.3],
        [1.5, 0.4, 0.6],
        [0.8, 1.2, 0.1],
        [1.1, 0.9, 1.4],
    ];
    let charges = vec![1.0, -1.0, 0.5, -0.5];
    let cell = make_cell_cubic(2.0);
    let input = make_input(
        positions,
        charges,
        vec![sigma, sigma, sigma, sigma],
        vec![eps, eps, eps, eps],
        vec![(0, 1)],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let result = evaluate_periodic(&input, &cell, 0.49).unwrap();
    let tol: f64 = 1e-10;
    for c in 0..3 {
        let sum: f64 = result.forces.iter().map(|f| f[c]).sum();
        assert!(
            sum.abs() <= tol,
            "sum of forces in component {} is {}, exceeds tolerance {}",
            c,
            sum,
            tol
        );
    }
}

// corrected after review: cutoff 1.0 equals h_min/2 and is rejected by validate_cutoff; 0.99 satisfies the strict inequality and excludes the same pairs
#[test]
fn test_type_indexed_parameters() {
    let cell = make_cell_cubic(2.0);
    let positions = vec![
        [0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, 0.0, 1.5],
        [0.0, 1.5, 0.0],
        [1.5, 0.0, 0.0],
    ];
    let types = vec![0, 1, 0, 1, 0];
    let lj_sigma = vec![0.34, 0.40];
    let lj_eps = vec![1.0, 2.0];
    let lj_cutoff: f64 = 0.99;
    let input = NonbondedInput {
        positions: positions.clone(),
        charges: vec![0.0; 5],
        lj_sigma: lj_sigma.clone(),
        lj_eps: lj_eps.clone(),
        lj_scale_14: 0.5,
        coul_scale_14: 0.5,
        bonds: vec![],
        rule: CombinationRule::LorentzBerthelot,
        types: types.clone(),
    };
    let result = evaluate_periodic(&input, &cell, lj_cutoff).unwrap();
    let n = positions.len();
    let mut expected: f64 = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let (d_min, r) = minimum_image_displacement(positions[i], positions[j], &cell).unwrap();
            if r < lj_cutoff {
                let sigma_ij = (lj_sigma[types[i]] + lj_sigma[types[j]]) / 2.0;
                let eps_ij = (lj_eps[types[i]] * lj_eps[types[j]]).sqrt();
                let r_j_mi = [
                    positions[i][0] - d_min[0],
                    positions[i][1] - d_min[1],
                    positions[i][2] - d_min[2],
                ];
                expected += lj_pair_energy(sigma_ij, eps_ij, positions[i], r_j_mi).unwrap();
            }
        }
    }
    assert!(
        (result.energy_lj_kjmol - expected).abs() <= 1e-12,
        "LJ energy: got {}, expected {}",
        result.energy_lj_kjmol,
        expected
    );
}

#[test]
fn test_type_indexed_parameters_geometric() {
    let cell = make_cell_cubic(2.0);
    let positions = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let types = vec![0, 1, 0];
    let lj_sigma = vec![0.34, 0.40];
    let lj_eps = vec![1.0, 2.0];
    let lj_cutoff: f64 = 0.99;
    let input = NonbondedInput {
        positions: positions.clone(),
        charges: vec![0.0; 3],
        lj_sigma: lj_sigma.clone(),
        lj_eps: lj_eps.clone(),
        lj_scale_14: 0.5,
        coul_scale_14: 0.5,
        bonds: vec![],
        rule: CombinationRule::Geometric,
        types: types.clone(),
    };
    let result = evaluate_periodic(&input, &cell, lj_cutoff).unwrap();
    let n = positions.len();
    let mut expected: f64 = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let (d_min, r) = minimum_image_displacement(positions[i], positions[j], &cell).unwrap();
            if r < lj_cutoff {
                let sigma_ij = (lj_sigma[types[i]] * lj_sigma[types[j]]).sqrt();
                let eps_ij = (lj_eps[types[i]] * lj_eps[types[j]]).sqrt();
                let r_j_mi = [
                    positions[i][0] - d_min[0],
                    positions[i][1] - d_min[1],
                    positions[i][2] - d_min[2],
                ];
                expected += lj_pair_energy(sigma_ij, eps_ij, positions[i], r_j_mi).unwrap();
            }
        }
    }
    assert!(
        (result.energy_lj_kjmol - expected).abs() <= 1e-12,
        "LJ energy: got {}, expected {}",
        result.energy_lj_kjmol,
        expected
    );
}

#[test]
fn test_cutoff_equal_boundary_cubic() {
    let cell = make_cell_cubic(2.0);
    let input = make_input(
        vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
        vec![0.0, 0.0],
        vec![0.34, 0.34],
        vec![1.0, 1.0],
        vec![],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let result = evaluate_periodic(&input, &cell, 1.0);
    assert!(
        matches!(result, Err(UcmdError::InvalidCell)),
        "expected Err(InvalidCell) for cutoff == h_min/2, got {:?}",
        result
    );
}

#[test]
fn test_cutoff_below_boundary_cubic() {
    let cell = make_cell_cubic(2.0);
    let input = make_input(
        vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
        vec![0.0, 0.0],
        vec![0.34, 0.34],
        vec![1.0, 1.0],
        vec![],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let result = evaluate_periodic(&input, &cell, 0.99);
    assert!(
        result.is_ok(),
        "expected Ok for cutoff < h_min/2, got {:?}",
        result
    );
}

#[test]
fn test_cutoff_above_boundary_cubic() {
    let cell = make_cell_cubic(2.0);
    let input = make_input(
        vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
        vec![0.0, 0.0],
        vec![0.34, 0.34],
        vec![1.0, 1.0],
        vec![],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let result = evaluate_periodic(&input, &cell, 1.01);
    assert!(
        matches!(result, Err(UcmdError::InvalidCell)),
        "expected Err(InvalidCell) for cutoff > h_min/2, got {:?}",
        result
    );
}

#[test]
fn test_cutoff_equal_boundary_triclinic() {
    let cell = Cell {
        ax: 2.0,
        bx: 0.5,
        by: 2.0,
        cx: 0.3,
        cy: 0.4,
        cz: 2.0,
    };
    let heights = ucmd_ff::minimum_image::perpendicular_heights(&cell).unwrap();
    let h_min: f64 = heights[0].min(heights[1]).min(heights[2]);
    let input = make_input(
        vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
        vec![0.0, 0.0],
        vec![0.34, 0.34],
        vec![1.0, 1.0],
        vec![],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let result = evaluate_periodic(&input, &cell, h_min / 2.0);
    assert!(
        matches!(result, Err(UcmdError::InvalidCell)),
        "expected Err(InvalidCell) for cutoff == h_min/2, got {:?}",
        result
    );
}

#[test]
fn test_cutoff_below_boundary_triclinic() {
    let cell = Cell {
        ax: 2.0,
        bx: 0.5,
        by: 2.0,
        cx: 0.3,
        cy: 0.4,
        cz: 2.0,
    };
    let heights = ucmd_ff::minimum_image::perpendicular_heights(&cell).unwrap();
    let h_min: f64 = heights[0].min(heights[1]).min(heights[2]);
    let input = make_input(
        vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
        vec![0.0, 0.0],
        vec![0.34, 0.34],
        vec![1.0, 1.0],
        vec![],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let result = evaluate_periodic(&input, &cell, h_min / 2.0 - 0.01);
    assert!(
        result.is_ok(),
        "expected Ok for cutoff < h_min/2, got {:?}",
        result
    );
}

#[test]
fn test_cutoff_above_boundary_triclinic() {
    let cell = Cell {
        ax: 2.0,
        bx: 0.5,
        by: 2.0,
        cx: 0.3,
        cy: 0.4,
        cz: 2.0,
    };
    let heights = ucmd_ff::minimum_image::perpendicular_heights(&cell).unwrap();
    let h_min: f64 = heights[0].min(heights[1]).min(heights[2]);
    let input = make_input(
        vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
        vec![0.0, 0.0],
        vec![0.34, 0.34],
        vec![1.0, 1.0],
        vec![],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let result = evaluate_periodic(&input, &cell, h_min / 2.0 + 0.01);
    assert!(
        matches!(result, Err(UcmdError::InvalidCell)),
        "expected Err(InvalidCell) for cutoff > h_min/2, got {:?}",
        result
    );
}

#[test]
fn test_cutoff_validation_uses_all_three_heights() {
    let cell = Cell {
        ax: 3.0,
        bx: 0.0,
        by: 3.0,
        cx: 0.0,
        cy: 0.0,
        cz: 1.0,
    };
    let input = make_input(
        vec![[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]],
        vec![0.0, 0.0],
        vec![0.34, 0.34],
        vec![1.0, 1.0],
        vec![],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let result = evaluate_periodic(&input, &cell, 1.5);
    assert!(
        matches!(result, Err(UcmdError::InvalidCell)),
        "expected Err(InvalidCell) because cutoff 1.5 > min_height/2 = 0.5, got {:?}",
        result
    );
}

#[test]
fn test_cutoff_validation_uses_height_a() {
    let cell = Cell {
        ax: 0.5,
        bx: 0.0,
        by: 4.0,
        cx: 0.0,
        cy: 0.0,
        cz: 4.0,
    };
    let input = make_input(
        vec![[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]],
        vec![0.0, 0.0],
        vec![0.34],
        vec![1.0],
        vec![],
        CombinationRule::LorentzBerthelot,
        0.5,
        0.5,
    );
    let result = evaluate_periodic(&input, &cell, 0.75);
    assert!(
        matches!(result, Err(UcmdError::InvalidCell)),
        "expected Err(InvalidCell) because cutoff 0.75 > h_a/2 = 0.25, got {:?}",
        result
    );
}
