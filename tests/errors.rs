//! Tests of the `errors` module.
#![allow(clippy::unwrap_used, clippy::expect_used)] // test code

use ucmd_ff::errors::*;

#[test]
fn non_finite_coordinate_variant_is_matchable() {
    let e = UcmdError::NonFiniteCoordinate;
    let got = matches!(&e, UcmdError::NonFiniteCoordinate) as i64;
    let expected = 1i64;
    assert!(
        (got - expected).abs() <= 0,
        "got {got}, expected {expected}"
    );
}

#[test]
fn coincident_atoms_variant_is_matchable() {
    let e = UcmdError::CoincidentAtoms;
    let got = matches!(&e, UcmdError::CoincidentAtoms) as i64;
    let expected = 1i64;
    assert!(
        (got - expected).abs() <= 0,
        "got {got}, expected {expected}"
    );
}

#[test]
fn index_out_of_range_variant_is_matchable() {
    let e = UcmdError::IndexOutOfRange;
    let got = matches!(&e, UcmdError::IndexOutOfRange) as i64;
    let expected = 1i64;
    assert!(
        (got - expected).abs() <= 0,
        "got {got}, expected {expected}"
    );
}

#[test]
fn invalid_cell_variant_is_matchable() {
    let e = UcmdError::InvalidCell;
    let got = matches!(&e, UcmdError::InvalidCell) as i64;
    let expected = 1i64;
    assert!(
        (got - expected).abs() <= 0,
        "got {got}, expected {expected}"
    );
}

#[test]
fn variants_are_distinct() {
    let a = UcmdError::NonFiniteCoordinate;
    let b = UcmdError::CoincidentAtoms;
    let c = UcmdError::IndexOutOfRange;
    let d = UcmdError::InvalidCell;
    let pairs = [(&a, &b), (&a, &c), (&a, &d), (&b, &c), (&b, &d), (&c, &d)];
    for (x, y) in pairs {
        let got = (std::mem::discriminant(x) != std::mem::discriminant(y)) as i64;
        let expected = 1i64;
        assert!(
            (got - expected).abs() <= 0,
            "got {got}, expected {expected}"
        );
    }
}

#[test]
fn error_propagates_through_result() {
    let r: Result<(), UcmdError> = Err(UcmdError::InvalidCell);
    let got = (r.is_err() && matches!(&r, Err(UcmdError::InvalidCell))) as i64;
    let expected = 1i64;
    assert!(
        (got - expected).abs() <= 0,
        "got {got}, expected {expected}"
    );
}

#[test]
fn every_variant_has_a_message_and_is_an_error() {
    let all = [
        UcmdError::NonFiniteCoordinate,
        UcmdError::CoincidentAtoms,
        UcmdError::IndexOutOfRange,
        UcmdError::InvalidCell,
        UcmdError::InvalidParameter,
        UcmdError::NonNeutralSystem,
        UcmdError::MoleculeSplit,
    ];
    let mut seen = std::collections::HashSet::new();
    for e in all {
        let text = e.to_string();
        assert!(!text.is_empty());
        assert!(seen.insert(text), "duplicate message for {e:?}");
        let _: &dyn std::error::Error = &e;
    }
}
