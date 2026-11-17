//! Computes the IUPAC dihedral angle and its Cartesian gradient for four atoms.
//!
//! Units: positions in nm; angle in rad; gradients in rad nm^-1.

use crate::errors::UcmdError;
use crate::vec3::PerAtom4;
use crate::vec3::{cross, dot, norm, scale, sub};

type Mat3 = [[f64; 3]; 3];

/// Computes the IUPAC dihedral angle phi for atoms i, j, k, l.
/// Units: positions in nm; return in rad.
pub fn dihedral_angle(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
) -> Result<f64, UcmdError> {
    let state = compute_state(r_i, r_j, r_k, r_l)?;
    let mut phi = state.y.atan2(state.x);
    if phi <= -std::f64::consts::PI {
        phi = std::f64::consts::PI;
    }
    Ok(phi)
}

/// Computes the gradient of the IUPAC dihedral angle phi with respect to the four positions.
/// Units: positions in nm; return gradients in rad nm^-1.
pub fn dihedral_gradient(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
) -> Result<PerAtom4, UcmdError> {
    let state = compute_state(r_i, r_j, r_k, r_l)?;
    let d = state.x * state.x + state.y * state.y;
    if d <= 0.0 {
        return Err(UcmdError::CoincidentAtoms);
    }

    let (g_n1_i, g_n1_j, g_n1_k, g_n1_l) = compute_n1_gradients(&state);
    let (g_n2_i, g_n2_j, g_n2_k, g_n2_l) = compute_n2_gradients(&state);
    let (g_v_i, g_v_j, g_v_k, g_v_l) = compute_v_gradients(&state);

    let n1_x_n2 = cross(state.n1, state.n2);

    let grad_x_i = vec_add(dot_vec_mat(state.n2, g_n1_i), dot_vec_mat(state.n1, g_n2_i));
    let grad_x_j = vec_add(dot_vec_mat(state.n2, g_n1_j), dot_vec_mat(state.n1, g_n2_j));
    let grad_x_k = vec_add(dot_vec_mat(state.n2, g_n1_k), dot_vec_mat(state.n1, g_n2_k));
    let grad_x_l = vec_add(dot_vec_mat(state.n2, g_n1_l), dot_vec_mat(state.n1, g_n2_l));

    let grad_y_i = vec_add(
        dot_vec_mat(
            state.v,
            add_mat(
                cross_mat_vec(g_n1_i, state.n2),
                vec_cross_mat(state.n1, g_n2_i),
            ),
        ),
        dot_vec_mat(n1_x_n2, g_v_i),
    );
    let grad_y_j = vec_add(
        dot_vec_mat(
            state.v,
            add_mat(
                cross_mat_vec(g_n1_j, state.n2),
                vec_cross_mat(state.n1, g_n2_j),
            ),
        ),
        dot_vec_mat(n1_x_n2, g_v_j),
    );
    let grad_y_k = vec_add(
        dot_vec_mat(
            state.v,
            add_mat(
                cross_mat_vec(g_n1_k, state.n2),
                vec_cross_mat(state.n1, g_n2_k),
            ),
        ),
        dot_vec_mat(n1_x_n2, g_v_k),
    );
    let grad_y_l = vec_add(
        dot_vec_mat(
            state.v,
            add_mat(
                cross_mat_vec(g_n1_l, state.n2),
                vec_cross_mat(state.n1, g_n2_l),
            ),
        ),
        dot_vec_mat(n1_x_n2, g_v_l),
    );

    let g_i = [
        (state.x * grad_y_i[0] - state.y * grad_x_i[0]) / d,
        (state.x * grad_y_i[1] - state.y * grad_x_i[1]) / d,
        (state.x * grad_y_i[2] - state.y * grad_x_i[2]) / d,
    ];
    let g_j = [
        (state.x * grad_y_j[0] - state.y * grad_x_j[0]) / d,
        (state.x * grad_y_j[1] - state.y * grad_x_j[1]) / d,
        (state.x * grad_y_j[2] - state.y * grad_x_j[2]) / d,
    ];
    let g_k = [
        (state.x * grad_y_k[0] - state.y * grad_x_k[0]) / d,
        (state.x * grad_y_k[1] - state.y * grad_x_k[1]) / d,
        (state.x * grad_y_k[2] - state.y * grad_x_k[2]) / d,
    ];
    let g_l = [
        (state.x * grad_y_l[0] - state.y * grad_x_l[0]) / d,
        (state.x * grad_y_l[1] - state.y * grad_x_l[1]) / d,
        (state.x * grad_y_l[2] - state.y * grad_x_l[2]) / d,
    ];

    Ok((g_i, g_j, g_k, g_l))
}

struct State {
    x: f64,
    y: f64,
    n1: [f64; 3],
    n2: [f64; 3],
    v: [f64; 3],
    u: [f64; 3],
    w: [f64; 3],
    l1: f64,
    l2: f64,
    l3: f64,
}

fn check_finite(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
) -> Result<(), UcmdError> {
    for v in [r_i, r_j, r_k, r_l] {
        for c in v {
            if !c.is_finite() {
                return Err(UcmdError::NonFiniteCoordinate);
            }
        }
    }
    Ok(())
}

fn compute_state(
    r_i: [f64; 3],
    r_j: [f64; 3],
    r_k: [f64; 3],
    r_l: [f64; 3],
) -> Result<State, UcmdError> {
    check_finite(r_i, r_j, r_k, r_l)?;

    let b1 = sub(r_j, r_i);
    let b2 = sub(r_k, r_j);
    let b3 = sub(r_l, r_k);

    let l1 = norm(b1);
    let l2 = norm(b2);
    let l3 = norm(b3);

    if l1 <= 1e-12 || l2 <= 1e-12 || l3 <= 1e-12 {
        return Err(UcmdError::CoincidentAtoms);
    }

    let u = scale(b1, 1.0 / l1);
    let v = scale(b2, 1.0 / l2);
    let w = scale(b3, 1.0 / l3);

    let n1 = cross(u, v);
    let n2 = cross(v, w);

    let n1_norm = norm(n1);
    let n2_norm = norm(n2);

    if n1_norm <= 1e-12 || n2_norm <= 1e-12 {
        return Err(UcmdError::CoincidentAtoms);
    }

    let x = dot(n1, n2);
    let y = dot(cross(n1, n2), v);

    Ok(State {
        x,
        y,
        n1,
        n2,
        v,
        u,
        w,
        l1,
        l2,
        l3,
    })
}

fn unit_vec_gradient_matrix(u: [f64; 3], l: f64, sign: f64) -> Mat3 {
    let mut g = [[0.0; 3]; 3];
    for a in 0..3 {
        for c in 0..3 {
            let delta = if c == a { 1.0 } else { 0.0 };
            g[c][a] = sign * (delta - u[c] * u[a]) / l;
        }
    }
    g
}

fn zero_mat() -> Mat3 {
    [[0.0; 3]; 3]
}

fn add_mat(a: Mat3, b: Mat3) -> Mat3 {
    let mut out = [[0.0; 3]; 3];
    for c in 0..3 {
        for a_idx in 0..3 {
            out[c][a_idx] = a[c][a_idx] + b[c][a_idx];
        }
    }
    out
}

fn cross_mat_vec(mat: Mat3, vec: [f64; 3]) -> Mat3 {
    let mut out = [[0.0; 3]; 3];
    for a in 0..3 {
        let col = [mat[0][a], mat[1][a], mat[2][a]];
        let cv = cross(col, vec);
        for c in 0..3 {
            out[c][a] = cv[c];
        }
    }
    out
}

fn vec_cross_mat(vec: [f64; 3], mat: Mat3) -> Mat3 {
    let mut out = [[0.0; 3]; 3];
    for a in 0..3 {
        let col = [mat[0][a], mat[1][a], mat[2][a]];
        let cv = cross(vec, col);
        for c in 0..3 {
            out[c][a] = cv[c];
        }
    }
    out
}

fn dot_vec_mat(vec: [f64; 3], mat: Mat3) -> [f64; 3] {
    let mut out = [0.0; 3];
    for a in 0..3 {
        let col = [mat[0][a], mat[1][a], mat[2][a]];
        out[a] = dot(vec, col);
    }
    out
}

fn vec_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn compute_n1_gradients(state: &State) -> (Mat3, Mat3, Mat3, Mat3) {
    let g_u_i = unit_vec_gradient_matrix(state.u, state.l1, -1.0);
    let g_u_j = unit_vec_gradient_matrix(state.u, state.l1, 1.0);
    let g_v_j = unit_vec_gradient_matrix(state.v, state.l2, -1.0);
    let g_v_k = unit_vec_gradient_matrix(state.v, state.l2, 1.0);

    let g_n1_i = cross_mat_vec(g_u_i, state.v);
    let g_n1_j = add_mat(cross_mat_vec(g_u_j, state.v), vec_cross_mat(state.u, g_v_j));
    let g_n1_k = vec_cross_mat(state.u, g_v_k);
    let g_n1_l = zero_mat();

    (g_n1_i, g_n1_j, g_n1_k, g_n1_l)
}

fn compute_n2_gradients(state: &State) -> (Mat3, Mat3, Mat3, Mat3) {
    let g_v_j = unit_vec_gradient_matrix(state.v, state.l2, -1.0);
    let g_v_k = unit_vec_gradient_matrix(state.v, state.l2, 1.0);
    let g_w_k = unit_vec_gradient_matrix(state.w, state.l3, -1.0);
    let g_w_l = unit_vec_gradient_matrix(state.w, state.l3, 1.0);

    let g_n2_i = zero_mat();
    let g_n2_j = cross_mat_vec(g_v_j, state.w);
    let g_n2_k = add_mat(cross_mat_vec(g_v_k, state.w), vec_cross_mat(state.v, g_w_k));
    let g_n2_l = vec_cross_mat(state.v, g_w_l);

    (g_n2_i, g_n2_j, g_n2_k, g_n2_l)
}

fn compute_v_gradients(state: &State) -> (Mat3, Mat3, Mat3, Mat3) {
    let g_v_j = unit_vec_gradient_matrix(state.v, state.l2, -1.0);
    let g_v_k = unit_vec_gradient_matrix(state.v, state.l2, 1.0);

    let g_v_i = zero_mat();
    let g_v_l = zero_mat();

    (g_v_i, g_v_j, g_v_k, g_v_l)
}
