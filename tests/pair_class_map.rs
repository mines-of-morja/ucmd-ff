//! Tests of the `pair_class_map` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::errors::UcmdError;
use ucmd_ff::exclusions::PairClass;
use ucmd_ff::pair_class_map::*;

#[test]
fn linear_chain_4_atom() {
    let bonds = [(0, 1), (1, 2), (2, 3)];
    let map = build(4, &bonds).unwrap();

    assert_eq!(map.class_at(0, 1).unwrap(), PairClass::Excluded);
    assert_eq!(map.class_at(0, 2).unwrap(), PairClass::Excluded);
    assert_eq!(map.class_at(0, 3).unwrap(), PairClass::OneFour);
    assert_eq!(map.class_at(1, 2).unwrap(), PairClass::Excluded);
    assert_eq!(map.class_at(1, 3).unwrap(), PairClass::Excluded);
    assert_eq!(map.class_at(2, 3).unwrap(), PairClass::Excluded);
}

#[test]
fn ring_4_atom_all_excluded() {
    let bonds = [(0, 1), (1, 2), (2, 3), (3, 0)];
    let map = build(4, &bonds).unwrap();

    assert_eq!(map.class_at(0, 1).unwrap(), PairClass::Excluded);
    assert_eq!(map.class_at(0, 2).unwrap(), PairClass::Excluded);
    assert_eq!(map.class_at(0, 3).unwrap(), PairClass::Excluded);
    assert_eq!(map.class_at(1, 2).unwrap(), PairClass::Excluded);
    assert_eq!(map.class_at(1, 3).unwrap(), PairClass::Excluded);
    assert_eq!(map.class_at(2, 3).unwrap(), PairClass::Excluded);
}

#[test]
fn two_disconnected_atoms_normal() {
    let bonds: [(usize, usize); 0] = [];
    let map = build(2, &bonds).unwrap();

    assert_eq!(map.class_at(0, 1).unwrap(), PairClass::Normal);
}

#[test]
fn index_out_of_range_returns_error() {
    let bonds = [(0, 5)];
    let result = build(3, &bonds);
    assert!(matches!(result, Err(UcmdError::IndexOutOfRange)));
}

#[test]
fn self_bond_ignored() {
    let bonds_with_self = [(0, 1), (1, 1), (1, 2)];
    let bonds_without_self = [(0, 1), (1, 2)];

    let map_with = build(3, &bonds_with_self).unwrap();
    let map_without = build(3, &bonds_without_self).unwrap();

    assert_eq!(map_with.class_at(0, 1).unwrap(), PairClass::Excluded);
    assert_eq!(map_with.class_at(0, 2).unwrap(), PairClass::Excluded);
    assert_eq!(map_with.class_at(1, 2).unwrap(), PairClass::Excluded);

    assert_eq!(map_with.classes, map_without.classes);
}

#[test]
fn class_at_symmetric_swap() {
    let bonds = [(0, 1), (1, 2), (2, 3), (3, 4)];
    let map = build(5, &bonds).unwrap();

    for i in 0..5 {
        for j in (i + 1)..5 {
            let c1 = map.class_at(i, j).unwrap();
            let c2 = map.class_at(j, i).unwrap();
            assert_eq!(c1, c2, "swap invariance failed for pair ({i}, {j})");
        }
    }
}

#[test]
fn disconnected_molecules_all_normal() {
    let bonds = [(0, 1), (1, 2), (3, 4), (4, 5)];
    let map = build(6, &bonds).unwrap();

    // Cross-molecule pairs are all Normal
    let cross_pairs = [
        (0, 3),
        (0, 4),
        (0, 5),
        (1, 3),
        (1, 4),
        (1, 5),
        (2, 3),
        (2, 4),
        (2, 5),
    ];
    for &(i, j) in &cross_pairs {
        assert_eq!(
            map.class_at(i, j).unwrap(),
            PairClass::Normal,
            "cross-molecule pair ({i}, {j}) should be Normal"
        );
    }

    // Within-molecule pairs follow the chain rule
    assert_eq!(map.class_at(0, 1).unwrap(), PairClass::Excluded);
    assert_eq!(map.class_at(0, 2).unwrap(), PairClass::Excluded);
    assert_eq!(map.class_at(1, 2).unwrap(), PairClass::Excluded);
    assert_eq!(map.class_at(3, 4).unwrap(), PairClass::Excluded);
    assert_eq!(map.class_at(3, 5).unwrap(), PairClass::Excluded);
    assert_eq!(map.class_at(4, 5).unwrap(), PairClass::Excluded);
}

#[test]
fn build_n1_empty_classes_no_panic() {
    let bonds: [(usize, usize); 0] = [];
    let result = build(1, &bonds);
    let map = result.expect("build(1, &[]) should succeed, not panic");
    assert_eq!(map.n, 1);
    assert!(map.classes.is_empty());
    assert_eq!(map.class_at(0, 0).unwrap(), PairClass::Excluded);
}
