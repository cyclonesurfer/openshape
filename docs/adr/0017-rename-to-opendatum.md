# ADR 0017 — The product is renamed OpenDatum

Status: accepted, 2026-10-12.

## Context

The project started as OpenShape. That name sounds too close to the name of
the commercial product whose feature set we benchmark against, in the same
field, for software styled after its workflow; a clean-room project should
not invite that confusion (docs/CLEANROOM.md). Pat had candidate names
checked against the US trademark register, crates.io, GitHub and the web
(the project's naming notes record each round) and chose **OpenDatum** on
2026-10-12. A manual search of the EU and WIPO registers, which could not
be run from the build environment, is still owed before any public launch;
it does not block the code.

The rename was planned to cover the crate and binary names, the window
title, the documentation and the document file extension. Some choices
were left open; this ADR settles them.

## Decisions

1. **Name in text:** `OpenDatum` in prose and the UI, `opendatum` in
   package, binary and file names, `OPENDATUM_` for environment variables.
   Crate packages become `opendatum-oplog`, `opendatum-kernel`, and so on;
   the binary becomes `opendatum`, and the application type
   `OpenShapeApp` becomes `OpenDatumApp`. Crate directories keep their
   role names (`crates/kernel`, …).
2. **The kernel type** `OsKernel` (the abbreviation of the old name) becomes
   `BrepKernel`. Options were `OdKernel` (the same abbreviation trick,
   cryptic) or keeping `OsKernel` (a stale name in every test). A name that
   says what it is (the B-rep kernel behind the `Kernel` trait) survives any
   later rename. Recommended and taken.
3. **Document extension:** `.odatum`. Options were `.oddoc` (the old
   pattern, but reads as "odd doc"), `.opendatum` (long) and `.odatum`
   (short, distinctive, plainly ours). Taken: `.odatum`. Opening does not
   look at the extension, so documents saved as `.osdoc` still open.
4. **Container magic:** new documents are written with
   `"OpenDatum-oplog"`; the reader also accepts the old `"OpenShape-oplog"`,
   so documents saved before the rename still open (the container format is
   otherwise unchanged, version 1). This is the one place the old name
   stays in code, as a documented legacy constant.
5. **Environment variables:** `OPENSHAPE_WGPU_BACKEND` and
   `OPENSHAPE_RENDER_TEST_OUT` become `OPENDATUM_…`. No compatibility
   alias: they are developer switches, and the hardware scripts are updated
   in the same commit.
6. **Repository:** the GitHub repository stays `cyclonesurfer/openshape`
   until Pat renames it (that needs admin rights; GitHub redirects old URLs,
   so existing clones keep pulling). `Cargo.toml`'s `repository` field keeps
   the current URL until then. Local clone directories (`~/OpenShape` on
   odin, `~/openshape` on huginn) are Pat's and are not renamed by us; the
   docs that give their paths keep them.
7. **History stays as it was:** commit messages are not rewritten. The
   documentation (STATUS.md with its dated entries, ADRs, design notes)
   describes the code as it is and uses the new name throughout: a package
   or command named in an old entry is the same one under its new name.

## Consequences

Every crate's package name and every `use openshape_…` path changes in one
commit, made on its own so it is easy to review and to bisect around. The
clean-room name guard (`crates/ui/tests/cleanroom.rs`) is unchanged and
keeps blocking the reference product's name. A repository-wide search
after the rename finds the old name only where decisions 4 and 6 keep it
(the legacy magic and its tests and docs, the repository URL, the clone
directories), in CLAUDE.md's one-line note of the rename and in this ADR.
