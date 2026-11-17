# ucmd-ff — agent notes

The shared rules are in `../CLAUDE.md`, in the UnitcellMD umbrella repo one level up in the
workspace. ucmd-ff is a force-field library that ucmd-engine uses as an in-process provider
(adapter crate `ucmd-providers-ff` in ucmd-engine, `kind = "ucmd-ff"`).

## Commands

```sh
cargo fmt --all --check
cargo clippy --all-targets -- -D warnings
cargo test
```

When chaining checks in a shell, use `set -o pipefail`.

## Repo rules

- **No dependencies.** The crate builds anywhere, without network access.
- **One constants module.** Physical constants live only in `constants.rs`, from CODATA 2018.
- **Errors propagate.** Library code has no `unwrap`/`expect`/`panic`; clippy enforces this.
  Every public item documents its units (`missing_docs`). Integration-test files start with
  `#![allow(clippy::unwrap_used, clippy::expect_used)]`.
- **Float range checks reject NaN.** Check finiteness first, or use positive forms
  (`x > 0.0`).
- **Every energy term is tested by finite differences:** forces, and the strain virial with the
  cell strained too in periodic systems (`fd_test`). Every reference value in a test names its
  source.
- **Never weaken a test.** Never loosen a tolerance, and never skip or delete a test to get a
  pass.
- **Changes to the public API also change the adapter** in ucmd-engine
  (`crates/ucmd-providers-ff`). The engine pins a ucmd-ff tag.
