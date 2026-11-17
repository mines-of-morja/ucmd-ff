//! Tests of the `cell` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::cell::*;
use ucmd_ff::errors::*;

fn assert_close_vec3(got: [f64; 3], expected: [f64; 3], tol: f64) {
    for k in 0..3 {
        assert!(
            (got[k] - expected[k]).abs() <= tol,
            "component {k}: got {got:?}, expected {expected:?}"
        );
    }
}

fn reduce_displacement(v: [f64; 3], cell: &Cell) -> [f64; 3] {
    let a = [cell.ax, 0.0, 0.0];
    let b = [cell.bx, cell.by, 0.0];
    let c = [cell.cx, cell.cy, cell.cz];

    let n_c = (v[2] / cell.cz).round() as i32;
    let v1 = [
        v[0] - n_c as f64 * c[0],
        v[1] - n_c as f64 * c[1],
        v[2] - n_c as f64 * c[2],
    ];

    let n_b = (v1[1] / cell.by).round() as i32;
    let v2 = [
        v1[0] - n_b as f64 * b[0],
        v1[1] - n_b as f64 * b[1],
        v1[2] - n_b as f64 * b[2],
    ];

    let n_a = (v2[0] / cell.ax).round() as i32;
    [
        v2[0] - n_a as f64 * a[0],
        v2[1] - n_a as f64 * a[1],
        v2[2] - n_a as f64 * a[2],
    ]
}

fn brute_force_minimum_image(r_i: [f64; 3], r_j: [f64; 3], cell: &Cell) -> [f64; 3] {
    let v = [r_i[0] - r_j[0], r_i[1] - r_j[1], r_i[2] - r_j[2]];
    let v_red = reduce_displacement(v, cell);

    let a = [cell.ax, 0.0, 0.0];
    let b = [cell.bx, cell.by, 0.0];
    let c = [cell.cx, cell.cy, cell.cz];

    let mut best = [f64::INFINITY; 3];
    let mut best_norm = f64::INFINITY;

    for n_a in -2..=2 {
        for n_b in -2..=2 {
            for n_c in -2..=2 {
                let la = [n_a as f64 * a[0], n_a as f64 * a[1], n_a as f64 * a[2]];
                let lb = [n_b as f64 * b[0], n_b as f64 * b[1], n_b as f64 * b[2]];
                let lc = [n_c as f64 * c[0], n_c as f64 * c[1], n_c as f64 * c[2]];
                let l = [
                    la[0] + lb[0] + lc[0],
                    la[1] + lb[1] + lc[1],
                    la[2] + lb[2] + lc[2],
                ];
                let d = [v_red[0] - l[0], v_red[1] - l[1], v_red[2] - l[2]];
                let norm_sq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                if norm_sq < best_norm {
                    best_norm = norm_sq;
                    best = d;
                }
            }
        }
    }
    best
}

#[test]
fn test_validate_cell_rejects_non_positive() {
    let cell = Cell {
        ax: -1.0,
        bx: 0.0,
        by: 1.0,
        cx: 0.0,
        cy: 0.0,
        cz: 1.0,
    };
    assert!(matches!(validate_cell(&cell), Err(UcmdError::InvalidCell)));
}

#[test]
fn test_validate_cell_accepts_valid() {
    let cell = Cell {
        ax: 2.0,
        bx: 0.5,
        by: 3.0,
        cx: 0.2,
        cy: 0.1,
        cz: 4.0,
    };
    assert!(validate_cell(&cell).is_ok());
}

#[test]
fn test_cell_volume_correct() {
    let cell = Cell {
        ax: 2.0,
        bx: 0.5,
        by: 3.0,
        cx: 0.2,
        cy: 0.1,
        cz: 4.0,
    };
    let expected = 2.0 * 3.0 * 4.0;
    let result = cell_volume(&cell);
    assert!(result.is_ok());
    let got = result.ok().unwrap_or(f64::NAN);
    assert!(
        (got - expected).abs() <= 1e-12,
        "got {got}, expected {expected}"
    );
}

// revised: brute-force helper now implements the 125-image search over the reduced displacement
#[test]
fn test_minimum_image_brute_force() {
    let cell = Cell {
        ax: 5.0,
        bx: 1.0,
        by: 5.0,
        cx: 0.5,
        cy: 0.5,
        cz: 5.0,
    };
    let r_i = [1.2, 2.3, 3.4];
    let r_j = [4.1, 1.8, 2.9];
    let expected = brute_force_minimum_image(r_i, r_j, &cell);
    let result = minimum_image(r_i, r_j, &cell);
    assert!(result.is_ok());
    let got = result.ok().unwrap_or([f64::NAN; 3]);
    assert_close_vec3(got, expected, 1e-12);
}

#[test]
fn test_minimum_image_translation_invariant() {
    let cell = Cell {
        ax: 5.0,
        bx: 1.0,
        by: 5.0,
        cx: 0.5,
        cy: 0.5,
        cz: 5.0,
    };
    let r_i = [1.2, 2.3, 3.4];
    let r_j = [4.1, 1.8, 2.9];
    let a = [cell.ax, 0.0, 0.0];
    let r_i_t = [r_i[0] + a[0], r_i[1] + a[1], r_i[2] + a[2]];
    let r_j_t = [r_j[0] + a[0], r_j[1] + a[1], r_j[2] + a[2]];

    let result1 = minimum_image(r_i, r_j, &cell);
    assert!(result1.is_ok());
    let d1 = result1.ok().unwrap_or([f64::NAN; 3]);

    let result2 = minimum_image(r_i_t, r_j_t, &cell);
    assert!(result2.is_ok());
    let d2 = result2.ok().unwrap_or([f64::NAN; 3]);

    assert_close_vec3(d1, d2, 1e-12);
}

#[test]
fn test_minimum_image_symmetry() {
    let cell = Cell {
        ax: 5.0,
        bx: 1.0,
        by: 5.0,
        cx: 0.5,
        cy: 0.5,
        cz: 5.0,
    };
    let r_i = [1.2, 2.3, 3.4];
    let r_j = [4.1, 1.8, 2.9];

    let result1 = minimum_image(r_i, r_j, &cell);
    assert!(result1.is_ok());
    let d1 = result1.ok().unwrap_or([f64::NAN; 3]);

    let result2 = minimum_image(r_j, r_i, &cell);
    assert!(result2.is_ok());
    let d2 = result2.ok().unwrap_or([f64::NAN; 3]);

    let neg_d2 = [-d2[0], -d2[1], -d2[2]];
    assert_close_vec3(d1, neg_d2, 1e-12);
}

#[test]
fn test_minimum_image_zero_displacement() {
    let cell = Cell {
        ax: 5.0,
        bx: 1.0,
        by: 5.0,
        cx: 0.5,
        cy: 0.5,
        cz: 5.0,
    };
    let r = [1.0, 2.0, 3.0];
    let expected = [0.0, 0.0, 0.0];
    let result = minimum_image(r, r, &cell);
    assert!(result.is_ok());
    let got = result.ok().unwrap_or([f64::NAN; 3]);
    assert_close_vec3(got, expected, 1e-15);
}

#[test]
fn test_minimum_image_far_apart() {
    let cell = Cell {
        ax: 5.0,
        bx: 1.0,
        by: 5.0,
        cx: 0.5,
        cy: 0.5,
        cz: 5.0,
    };
    let r_i = [1.2, 2.3, 3.4];
    let r_j = [4.1, 1.8, 2.9];

    let a = [cell.ax, 0.0, 0.0];
    let b = [cell.bx, cell.by, 0.0];
    let c = [cell.cx, cell.cy, cell.cz];

    let r_i_far = [
        r_i[0] + 3.0 * a[0] + 3.0 * b[0] + 3.0 * c[0],
        r_i[1] + 3.0 * a[1] + 3.0 * b[1] + 3.0 * c[1],
        r_i[2] + 3.0 * a[2] + 3.0 * b[2] + 3.0 * c[2],
    ];
    let r_j_far = [
        r_j[0] + 7.0 * a[0] + 7.0 * b[0] + 7.0 * c[0],
        r_j[1] + 7.0 * a[1] + 7.0 * b[1] + 7.0 * c[1],
        r_j[2] + 7.0 * a[2] + 7.0 * b[2] + 7.0 * c[2],
    ];

    let result1 = minimum_image(r_i, r_j, &cell);
    assert!(result1.is_ok());
    let d1 = result1.ok().unwrap_or([f64::NAN; 3]);

    let result2 = minimum_image(r_i_far, r_j_far, &cell);
    assert!(result2.is_ok());
    let d2 = result2.ok().unwrap_or([f64::NAN; 3]);

    assert_close_vec3(d1, d2, 1e-12);
}

#[test]
fn test_minimum_image_lattice_shift() {
    let cell = Cell {
        ax: 5.0,
        bx: 1.0,
        by: 5.0,
        cx: 0.5,
        cy: 0.5,
        cz: 5.0,
    };
    let r_i = [1.2, 2.3, 3.4];
    let r_j = [4.1, 1.8, 2.9];

    let a = [cell.ax, 0.0, 0.0];
    let b = [cell.bx, cell.by, 0.0];
    let c = [cell.cx, cell.cy, cell.cz];

    let result_base = minimum_image(r_i, r_j, &cell);
    assert!(result_base.is_ok());
    let d_base = result_base.ok().unwrap_or([f64::NAN; 3]);

    for n_a in -5..=5 {
        for n_b in -5..=5 {
            for n_c in -5..=5 {
                let r_i_shifted = [
                    r_i[0] + n_a as f64 * a[0] + n_b as f64 * b[0] + n_c as f64 * c[0],
                    r_i[1] + n_a as f64 * a[1] + n_b as f64 * b[1] + n_c as f64 * c[1],
                    r_i[2] + n_a as f64 * a[2] + n_b as f64 * b[2] + n_c as f64 * c[2],
                ];
                let result = minimum_image(r_i_shifted, r_j, &cell);
                assert!(result.is_ok());
                let d = result.ok().unwrap_or([f64::NAN; 3]);
                assert_close_vec3(d, d_base, 1e-12);
            }
        }
    }
}

#[test]
fn test_minimum_image_125_images() {
    let cell = Cell {
        ax: 5.0,
        bx: 1.0,
        by: 5.0,
        cx: 0.5,
        cy: 0.5,
        cz: 5.0,
    };
    let r_i = [1.2, 2.3, 3.4];
    let r_j = [4.1, 1.8, 2.9];

    let result = minimum_image(r_i, r_j, &cell);
    assert!(result.is_ok());
    let d = result.ok().unwrap_or([f64::NAN; 3]);
    let d_norm_sq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];

    let v = [r_i[0] - r_j[0], r_i[1] - r_j[1], r_i[2] - r_j[2]];
    let v_red = reduce_displacement(v, &cell);

    let a = [cell.ax, 0.0, 0.0];
    let b = [cell.bx, cell.by, 0.0];
    let c = [cell.cx, cell.cy, cell.cz];

    for n_a in -2..=2 {
        for n_b in -2..=2 {
            for n_c in -2..=2 {
                let la = [n_a as f64 * a[0], n_a as f64 * a[1], n_a as f64 * a[2]];
                let lb = [n_b as f64 * b[0], n_b as f64 * b[1], n_b as f64 * b[2]];
                let lc = [n_c as f64 * c[0], n_c as f64 * c[1], n_c as f64 * c[2]];
                let l = [
                    la[0] + lb[0] + lc[0],
                    la[1] + lb[1] + lc[1],
                    la[2] + lb[2] + lc[2],
                ];
                let candidate = [v_red[0] - l[0], v_red[1] - l[1], v_red[2] - l[2]];
                let candidate_norm_sq = candidate[0] * candidate[0]
                    + candidate[1] * candidate[1]
                    + candidate[2] * candidate[2];
                assert!(
                    d_norm_sq <= candidate_norm_sq + 1e-12,
                    "result norm_sq {d_norm_sq} > candidate norm_sq {candidate_norm_sq} for n_a={n_a}, n_b={n_b}, n_c={n_c}"
                );
            }
        }
    }
}

#[test]
fn test_minimum_image_highly_skewed() {
    let cell = Cell {
        ax: 1.0,
        bx: 100.0,
        by: 1.0,
        cx: 100.0,
        cy: 1.0,
        cz: 1.0,
    };
    let r_i = [50.0, 0.5, 0.5];
    let r_j = [0.0, 0.0, 0.0];

    let expected = brute_force_minimum_image(r_i, r_j, &cell);
    let result = minimum_image(r_i, r_j, &cell);
    assert!(result.is_ok());
    let got = result.ok().unwrap_or([f64::NAN; 3]);
    assert_close_vec3(got, expected, 1e-12);
}

#[test]
fn test_minimum_image_skewed_cell_x_component() {
    let cell = Cell {
        ax: 10.0,
        bx: 1.0,
        by: 10.0,
        cx: 0.0,
        cy: 0.0,
        cz: 10.0,
    };
    let r_i = [5.2, 5.0, 0.0];
    let r_j = [0.0, 0.0, 0.0];
    let expected = [4.2, -5.0, 0.0];
    let result = minimum_image(r_i, r_j, &cell);
    assert!(result.is_ok());
    let got = result.ok().unwrap_or([f64::NAN; 3]);
    assert_close_vec3(got, expected, 1e-12);
}
