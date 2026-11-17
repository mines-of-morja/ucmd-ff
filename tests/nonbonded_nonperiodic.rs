//! Tests of the `nonbonded_nonperiodic` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code
// Index loops over atoms and axes mirror the formulas they test.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

use ucmd_ff::combination_rule::*;
use ucmd_ff::coulomb_pair::*;
use ucmd_ff::lj_pair::*;
use ucmd_ff::nonbonded_nonperiodic::*;
use ucmd_ff::vec3::*;

fn rel_close(got: f64, exp: f64, tol: f64) -> bool {
    (got - exp).abs() <= tol * exp.abs().max(1e-30)
}

fn vec_close(got: [f64; 3], exp: [f64; 3], tol: f64) -> bool {
    for a in 0..3 {
        if !rel_close(got[a], exp[a], tol) {
            return false;
        }
    }
    true
}

#[test]
fn two_atom_energy_matches_pair_sum() {
    let p = [[0.0, 0.0, 0.0], [0.4, 0.0, 0.0]];
    let input = NonbondedInput {
        positions: p.to_vec(),
        types: vec![0, 0],
        charges: vec![1.0, -1.0],
        lj_sigma: vec![0.34],
        lj_eps: vec![1.9998],
        rule: CombinationRule::Geometric,
        bonds: vec![],
        lj_scale_14: 1.0,
        coul_scale_14: 1.0,
    };
    let r = evaluate_nonperiodic(&input).unwrap();
    let el = lj_pair_energy(0.34, 1.9998, p[0], p[1]).unwrap();
    let ec = coulomb_pair_energy(1.0, -1.0, p[0], p[1]).unwrap();
    let (fil, fjl) = lj_pair_forces(0.34, 1.9998, p[0], p[1]).unwrap();
    let (fic, fjc) = coulomb_pair_forces(1.0, -1.0, p[0], p[1]).unwrap();
    let fi = add(fil, fic);
    let fj = add(fjl, fjc);
    let rij = sub(p[0], p[1]);
    let mut w = [[0.0; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            w[a][b] = rij[a] * fi[b];
        }
    }
    assert!(
        rel_close(r.energy_lj, el, 1e-12),
        "lj: got {}, exp {}",
        r.energy_lj,
        el
    );
    assert!(
        rel_close(r.energy_coulomb, ec, 1e-12),
        "coul: got {}, exp {}",
        r.energy_coulomb,
        ec
    );
    assert!(
        vec_close(r.forces[0], fi, 1e-12),
        "f0: got {:?}, exp {:?}",
        r.forces[0],
        fi
    );
    assert!(
        vec_close(r.forces[1], fj, 1e-12),
        "f1: got {:?}, exp {:?}",
        r.forces[1],
        fj
    );
    for a in 0..3 {
        for b in 0..3 {
            assert!(
                rel_close(r.virial[a][b], w[a][b], 1e-12),
                "w[{}][{}]: got {}, exp {}",
                a,
                b,
                r.virial[a][b],
                w[a][b]
            );
        }
    }
}

// corrected after review: 3-atom chain has no 1-4 pair; use 4-atom chain so pair (0,3) is 1-4
#[test]
fn three_atom_chain_14_scaling() {
    let p = [
        [0.0, 0.0, 0.0],
        [0.3, 0.0, 0.0],
        [0.6, 0.0, 0.0],
        [0.9, 0.0, 0.0],
    ];
    let input = NonbondedInput {
        positions: p.to_vec(),
        types: vec![0, 0, 0, 0],
        charges: vec![1.0, 0.0, 0.0, -1.0],
        lj_sigma: vec![0.34],
        lj_eps: vec![1.9998],
        rule: CombinationRule::Geometric,
        bonds: vec![(0, 1), (1, 2), (2, 3)],
        lj_scale_14: 0.5,
        coul_scale_14: 0.8333,
    };
    let r = evaluate_nonperiodic(&input).unwrap();
    let el03 = lj_pair_energy(0.34, 1.9998, p[0], p[3]).unwrap();
    let el_exp = 0.5 * el03;
    let ec03 = coulomb_pair_energy(1.0, -1.0, p[0], p[3]).unwrap();
    let ec_exp = 0.8333 * ec03;
    assert!(
        rel_close(r.energy_lj, el_exp, 1e-12),
        "lj: got {}, exp {}",
        r.energy_lj,
        el_exp
    );
    assert!(
        rel_close(r.energy_coulomb, ec_exp, 1e-12),
        "coul: got {}, exp {}",
        r.energy_coulomb,
        ec_exp
    );
}

#[test]
fn newtons_third_law() {
    let p = [
        [0.0, 0.0, 0.0],
        [0.5, 0.2, 0.1],
        [1.0, 0.3, 0.2],
        [0.3, 0.7, 0.4],
        [1.2, 0.1, 0.5],
    ];
    let input = NonbondedInput {
        positions: p.to_vec(),
        types: vec![0, 0, 1, 1, 0],
        charges: vec![0.5, -0.5, 0.3, -0.3, 0.1],
        lj_sigma: vec![0.34, 0.28],
        lj_eps: vec![1.9998, 0.5],
        rule: CombinationRule::LorentzBerthelot,
        bonds: vec![(0, 1), (1, 2), (2, 3), (3, 4)],
        lj_scale_14: 0.5,
        coul_scale_14: 0.8333,
    };
    let r = evaluate_nonperiodic(&input).unwrap();
    let mut sum = [0.0; 3];
    for i in 0..5 {
        for a in 0..3 {
            sum[a] += r.forces[i][a];
        }
    }
    for a in 0..3 {
        assert!(sum[a].abs() <= 1e-12, "sum[{}]: got {}", a, sum[a]);
    }
}

#[test]
fn finite_difference_force() {
    let p = [
        [0.0, 0.0, 0.0],
        [0.5, 0.2, 0.1],
        [1.0, 0.3, 0.2],
        [0.3, 0.7, 0.4],
        [1.2, 0.1, 0.5],
    ];
    let types = vec![0, 0, 1, 1, 0];
    let charges = vec![0.5, -0.5, 0.3, -0.3, 0.1];
    let lj_sigma = vec![0.34, 0.28];
    let lj_eps = vec![1.9998, 0.5];
    let rule = CombinationRule::LorentzBerthelot;
    let bonds = vec![(0, 1), (1, 2), (2, 3), (3, 4)];
    let lj_scale_14 = 0.5;
    let coul_scale_14 = 0.8333;
    let input = NonbondedInput {
        positions: p.to_vec(),
        types: types.clone(),
        charges: charges.clone(),
        lj_sigma: lj_sigma.clone(),
        lj_eps: lj_eps.clone(),
        rule,
        bonds: bonds.clone(),
        lj_scale_14,
        coul_scale_14,
    };
    let r = evaluate_nonperiodic(&input).unwrap();
    let h: f64 = 1e-8;
    for k in 0..5 {
        for a in 0..3 {
            let mut pp = p.to_vec();
            let mut pm = p.to_vec();
            pp[k][a] += h;
            pm[k][a] -= h;
            let ip = NonbondedInput {
                positions: pp,
                types: types.clone(),
                charges: charges.clone(),
                lj_sigma: lj_sigma.clone(),
                lj_eps: lj_eps.clone(),
                rule,
                bonds: bonds.clone(),
                lj_scale_14,
                coul_scale_14,
            };
            let im = NonbondedInput {
                positions: pm,
                types: types.clone(),
                charges: charges.clone(),
                lj_sigma: lj_sigma.clone(),
                lj_eps: lj_eps.clone(),
                rule,
                bonds: bonds.clone(),
                lj_scale_14,
                coul_scale_14,
            };
            let ep = evaluate_nonperiodic(&ip).unwrap();
            let em = evaluate_nonperiodic(&im).unwrap();
            let fd =
                -(ep.energy_lj + ep.energy_coulomb - em.energy_lj - em.energy_coulomb) / (2.0 * h);
            assert!(
                (r.forces[k][a] - fd).abs() <= 1e-5,
                "F[{}][{}]: got {}, fd {}",
                k,
                a,
                r.forces[k][a],
                fd
            );
        }
    }
}

#[test]
fn finite_difference_virial() {
    let p = [
        [0.0, 0.0, 0.0],
        [0.5, 0.2, 0.1],
        [1.0, 0.3, 0.2],
        [0.3, 0.7, 0.4],
        [1.2, 0.1, 0.5],
    ];
    let types = vec![0, 0, 1, 1, 0];
    let charges = vec![0.5, -0.5, 0.3, -0.3, 0.1];
    let lj_sigma = vec![0.34, 0.28];
    let lj_eps = vec![1.9998, 0.5];
    let rule = CombinationRule::LorentzBerthelot;
    let bonds = vec![(0, 1), (1, 2), (2, 3), (3, 4)];
    let lj_scale_14 = 0.5;
    let coul_scale_14 = 0.8333;
    let input = NonbondedInput {
        positions: p.to_vec(),
        types: types.clone(),
        charges: charges.clone(),
        lj_sigma: lj_sigma.clone(),
        lj_eps: lj_eps.clone(),
        rule,
        bonds: bonds.clone(),
        lj_scale_14,
        coul_scale_14,
    };
    let r = evaluate_nonperiodic(&input).unwrap();
    let eps: f64 = 1e-8;
    for a in 0..3 {
        for b in 0..3 {
            let mut pp = p.to_vec();
            let mut pm = p.to_vec();
            for i in 0..5 {
                pp[i][b] += eps * p[i][a];
                pm[i][b] -= eps * p[i][a];
            }
            let ip = NonbondedInput {
                positions: pp,
                types: types.clone(),
                charges: charges.clone(),
                lj_sigma: lj_sigma.clone(),
                lj_eps: lj_eps.clone(),
                rule,
                bonds: bonds.clone(),
                lj_scale_14,
                coul_scale_14,
            };
            let im = NonbondedInput {
                positions: pm,
                types: types.clone(),
                charges: charges.clone(),
                lj_sigma: lj_sigma.clone(),
                lj_eps: lj_eps.clone(),
                rule,
                bonds: bonds.clone(),
                lj_scale_14,
                coul_scale_14,
            };
            let ep = evaluate_nonperiodic(&ip).unwrap();
            let em = evaluate_nonperiodic(&im).unwrap();
            let fd = -(ep.energy_lj + ep.energy_coulomb - em.energy_lj - em.energy_coulomb)
                / (2.0 * eps);
            assert!(
                (r.virial[a][b] - fd).abs() <= 1e-5,
                "W[{}][{}]: got {}, fd {}",
                a,
                b,
                r.virial[a][b],
                fd
            );
        }
    }
}

#[test]
fn zero_sigma_zero_charge() {
    let input = NonbondedInput {
        positions: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        types: vec![0, 0],
        charges: vec![0.0, 0.0],
        lj_sigma: vec![0.0],
        lj_eps: vec![1.0],
        rule: CombinationRule::Geometric,
        bonds: vec![],
        lj_scale_14: 0.5,
        coul_scale_14: 0.5,
    };
    let r = evaluate_nonperiodic(&input).unwrap();
    assert!(r.energy_lj == 0.0, "lj: got {}", r.energy_lj);
    assert!(r.energy_coulomb == 0.0, "coul: got {}", r.energy_coulomb);
    for i in 0..2 {
        for a in 0..3 {
            assert!(
                r.forces[i][a] == 0.0,
                "f[{}][{}]: got {}",
                i,
                a,
                r.forces[i][a]
            );
        }
    }
    for a in 0..3 {
        for b in 0..3 {
            assert!(
                r.virial[a][b] == 0.0,
                "w[{}][{}]: got {}",
                a,
                b,
                r.virial[a][b]
            );
        }
    }
}
