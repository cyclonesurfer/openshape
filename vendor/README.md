# Patched copies of dependencies

Each directory here is a published crates.io package, copied unmodified
except for the patches listed below, and used in place of the published
crate through `[patch.crates-io]` in the workspace's `Cargo.toml` (ADR
0022). The copies keep their own licences (the files in each directory);
cargo-deny checks them like any dependency. They are not workspace
members: our lints, formatting and tests do not apply to them.

On each upgrade of a crate, replace its copy with the new release and
reapply the patches, or drop the copy if upstream has fixed the defect.

## spade 2.15.1 (MIT OR Apache-2.0)

Source: the crates.io package `spade` 2.15.1
(https://github.com/Stoeoef/spade). Copied: `src`, `images` (included by
the documentation), `README.md`, the licences, and `Cargo.toml` without
the example and benchmark targets and the development dependencies (none
of them are part of the copy).

Patch 1, `src/delaunay_core/bulk_load.rs` (`try_bulk_load_cdt`): the
vertices and constraint edges deferred during the bulk load are kept in
`BTreeSet`s instead of `std::collections::HashSet`s, so they are inserted
in the same order on every run. The standard library seeds hash-set
iteration order at random per process, and the order changes the
triangulation of points on a common circle, so the same polygon was
triangulated differently from run to run. Test:
`meshes_and_booleans_are_the_same_every_time` in
`crates/kernel/tests/golden.rs` fails with the published crate.
