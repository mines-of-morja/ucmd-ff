//! Compact pair-classification map for non-bonded pair evaluation.
//! Units: none (categorical).
//!
//! Wraps `exclusions::build_exclusion_map` into a flat vector containing only
//! the N(N−1)/2 entries for pairs (i, j) with i < j, ordered lexicographically.

use crate::errors::UcmdError;
use crate::exclusions::{PairClass, build_exclusion_map};

/// Compact pair-classification map for N atoms.
/// Units: none (categorical).
///
/// Stores the `PairClass` for every pair (i, j) with i < j in a flat vector
/// of length N*(N-1)/2. Pairs are ordered lexicographically: (0,1), (0,2), …, (0,N-1), (1,2), …
#[derive(Debug, Clone, PartialEq)]
pub struct PairClassMap {
    /// Total number of atoms. Units: none (count).
    pub n: usize,
    /// Flat array of pair classes, one entry per pair (i, j) with i < j.
    /// Units: none (categorical).
    pub classes: Vec<PairClass>,
}

/// Builds a `PairClassMap` for `n` atoms given a list of bonded pairs.
/// Units: none (atom indices and categorical output).
///
/// # Arguments
/// * `n` — total number of atoms (≥ 0).
/// * `bonds` — bonded pairs as (i, j) index tuples; self-bonds are ignored.
///
/// # Returns
/// A `PairClassMap` with `n*(n-1)/2` entries.
///
/// # Errors
/// * `IndexOutOfRange` — any bond index is ≥ `n`.
pub fn build(n: usize, bonds: &[(usize, usize)]) -> Result<PairClassMap, UcmdError> {
    let matrix = build_exclusion_map(n, bonds)?;

    let mut classes = Vec::with_capacity(n * n.saturating_sub(1) / 2);
    for (i, row) in matrix.iter().enumerate() {
        classes.extend_from_slice(&row[i + 1..]);
    }

    Ok(PairClassMap { n, classes })
}

impl PairClassMap {
    /// Returns the `PairClass` for the pair (i, j).
    /// Units: none (categorical).
    ///
    /// If i == j, returns `Ok(PairClass::Excluded)` (self-pairs never interact).
    /// If i > j, the arguments are swapped before lookup.
    /// If i or j ≥ n, returns `Err(UcmdError::IndexOutOfRange)`.
    pub fn class_at(&self, i: usize, j: usize) -> Result<PairClass, UcmdError> {
        if i >= self.n || j >= self.n {
            return Err(UcmdError::IndexOutOfRange);
        }
        if i == j {
            return Ok(PairClass::Excluded);
        }
        let (i, j) = if i > j { (j, i) } else { (i, j) };
        let n = self.n;
        let idx = i * (2 * n - 1 - i) / 2 + (j - i - 1);
        Ok(self.classes[idx])
    }
}
