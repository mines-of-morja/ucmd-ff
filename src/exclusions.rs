//! Builds a symmetric pair-classification matrix from a bond graph.
//! Units: none (atom indices and categorical output).
//!
//! For every pair of distinct atoms (i, j), the class is determined by the
//! shortest path length d in the undirected bond graph:
//! - d = 1 or d = 2 → `Excluded` (no non-bonded interaction)
//! - d = 3 → `OneFour` (scaled interaction)
//! - d > 3 or no path exists → `Normal` (full interaction)
//!
//! The diagonal (i = i) is always `Normal`. Self-bonds are silently ignored.

use crate::errors::UcmdError;

/// Classification of a non-bonded pair in force-field evaluation.
/// Units: none (categorical).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PairClass {
    /// Full non-bonded interaction (LJ and Coulomb at full strength).
    Normal,
    /// No non-bonded interaction (shortest path length 1 or 2).
    Excluded,
    /// Scaled non-bonded interaction (shortest path length 3, i.e. 1-4 pair).
    OneFour,
}

/// Builds the pair-classification matrix from a bond graph.
/// Units: none (atom indices and categorical output).
///
/// # Arguments
/// * `n_atoms` — total number of atoms (≥ 0).
/// * `bonds` — bonded pairs as (i, j) index tuples; self-bonds are ignored.
///
/// # Returns
/// An `n_atoms × n_atoms` matrix; entry `[i][j]` is the `PairClass` for pair (i, j).
///
/// # Errors
/// * `IndexOutOfRange` — any bond index is ≥ `n_atoms`.
pub fn build_exclusion_map(
    n_atoms: usize,
    bonds: &[(usize, usize)],
) -> Result<Vec<Vec<PairClass>>, UcmdError> {
    // Validate all bond indices
    for &(i, j) in bonds {
        if i >= n_atoms || j >= n_atoms {
            return Err(UcmdError::IndexOutOfRange);
        }
    }

    // Build undirected adjacency list; self-bonds are ignored
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n_atoms];
    for &(i, j) in bonds {
        if i != j {
            adj[i].push(j);
            adj[j].push(i);
        }
    }

    // Initialize the matrix with Normal
    let mut map: Vec<Vec<PairClass>> = vec![vec![PairClass::Normal; n_atoms]; n_atoms];

    // For each source atom, BFS up to depth 3 to find shortest path lengths
    for s in 0..n_atoms {
        // dist[t] = shortest path from s to t; usize::MAX means unreachable
        let mut dist: Vec<usize> = vec![usize::MAX; n_atoms];
        dist[s] = 0;

        // BFS using a Vec as a queue with a head pointer
        let mut queue: Vec<usize> = vec![s];
        let mut head = 0;
        while head < queue.len() {
            let u = queue[head];
            head += 1;
            if dist[u] >= 3 {
                continue;
            }
            for &v in &adj[u] {
                if dist[v] == usize::MAX {
                    dist[v] = dist[u] + 1;
                    queue.push(v);
                }
            }
        }

        // Classify each pair (s, t)
        for t in 0..n_atoms {
            if t == s {
                continue;
            }
            let d = dist[t];
            if d != usize::MAX {
                match d {
                    1 | 2 => map[s][t] = PairClass::Excluded,
                    3 => map[s][t] = PairClass::OneFour,
                    _ => {}
                }
            }
        }
    }

    Ok(map)
}
