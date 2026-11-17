//! Tests of the `exclusions` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::errors::*;
use ucmd_ff::exclusions::*;

#[test]
fn test_linear_chain_4() {
    let bonds = [(0, 1), (1, 2), (2, 3)];
    let map = build_exclusion_map(4, &bonds).unwrap();

    assert_eq!(map[0][1], PairClass::Excluded, "d(0,1)=1 → Excluded");
    assert_eq!(map[0][2], PairClass::Excluded, "d(0,2)=2 → Excluded");
    assert_eq!(map[0][3], PairClass::OneFour, "d(0,3)=3 → OneFour");
    assert_eq!(map[1][2], PairClass::Excluded, "d(1,2)=1 → Excluded");
    assert_eq!(map[1][3], PairClass::Excluded, "d(1,3)=2 → Excluded");
    assert_eq!(map[2][3], PairClass::Excluded, "d(2,3)=1 → Excluded");
}

#[test]
fn test_ring_4_shortest_path() {
    let bonds = [(0, 1), (1, 2), (2, 3), (3, 0)];
    let map = build_exclusion_map(4, &bonds).unwrap();

    assert_eq!(
        map[0][2],
        PairClass::Excluded,
        "d(0,2)=2 via 0-1-2 or 0-3-2 → Excluded"
    );
    assert_eq!(map[0][1], PairClass::Excluded, "d(0,1)=1 → Excluded");
    assert_eq!(map[0][3], PairClass::Excluded, "d(0,3)=1 → Excluded");
    assert_eq!(
        map[1][3],
        PairClass::Excluded,
        "d(1,3)=2 via 1-2-3 or 1-0-3 → Excluded"
    );
}

#[test]
fn test_linear_chain_5_path_beyond_3() {
    let bonds = [(0, 1), (1, 2), (2, 3), (3, 4)];
    let map = build_exclusion_map(5, &bonds).unwrap();

    assert_eq!(map[0][3], PairClass::OneFour, "d(0,3)=3 → OneFour");
    assert_eq!(map[0][4], PairClass::Normal, "d(0,4)=4 → Normal");
    assert_eq!(map[1][4], PairClass::OneFour, "d(1,4)=3 → OneFour");
    assert_eq!(map[0][2], PairClass::Excluded, "d(0,2)=2 → Excluded");
}

#[test]
fn test_disconnected_molecules() {
    let bonds = [(0, 1), (2, 3)];
    let map = build_exclusion_map(4, &bonds).unwrap();

    assert_eq!(map[0][1], PairClass::Excluded, "d(0,1)=1 → Excluded");
    assert_eq!(map[2][3], PairClass::Excluded, "d(2,3)=1 → Excluded");
    assert_eq!(map[0][2], PairClass::Normal, "no path → Normal");
    assert_eq!(map[0][3], PairClass::Normal, "no path → Normal");
    assert_eq!(map[1][2], PairClass::Normal, "no path → Normal");
    assert_eq!(map[1][3], PairClass::Normal, "no path → Normal");
}

#[test]
fn test_invalid_index_returns_error() {
    let bonds = [(0, 5)];
    let result = build_exclusion_map(3, &bonds);
    assert!(
        matches!(result, Err(UcmdError::IndexOutOfRange)),
        "expected Err(IndexOutOfRange), got {result:?}"
    );
}

#[test]
fn test_symmetry_and_diagonal() {
    let bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (0, 3)];
    let map = build_exclusion_map(5, &bonds).unwrap();

    for i in 0..5 {
        assert_eq!(
            map[i][i],
            PairClass::Normal,
            "diagonal [{}][{}] should be Normal",
            i,
            i
        );
        for j in 0..5 {
            assert_eq!(
                map[i][j], map[j][i],
                "symmetry violated at [{}][{}] vs [{}][{}]",
                i, j, j, i
            );
        }
    }
}

#[test]
fn test_zero_atoms() {
    let bonds: &[(usize, usize)] = &[];
    let map = build_exclusion_map(0, bonds).unwrap();
    assert!(
        map.is_empty(),
        "expected empty Vec, got {:?} entries",
        map.len()
    );
}
