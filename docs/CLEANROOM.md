# Clean-room policy

OpenShape is an independent implementation of a parametric 3D CAD system.
This document states how it is built so that the origin of every part of it
can be shown. It binds every contributor, human or AI agent.

## Functional reference

OpenShape's feature goals and its general workflow are benchmarked against
commercial browser-based parametric CAD systems, in particular Onshape
(a trademark of PTC Inc.; OpenShape is not affiliated with, endorsed by or
derived from it or its owner). "Functional reference" means only the
publicly known feature set and the publicly observable way such a system is
used: which modelling operations exist, and how a user generally moves
through them. It never means its implementation.

The name of that product appears only in this document, in `NOTICE` and in
the contributor instructions (`CLAUDE.md`), to state this boundary. It never
appears in the product: not in the code, the user interface, the
documentation shipped with the program, file formats or metadata. The test
`crates/ui/tests/cleanroom.rs` fails the build if it does.

## What is never used

- Source code of any proprietary CAD system, in any form.
- Decompiled, disassembled or otherwise reverse-engineered binaries,
  network protocols or file formats.
- Leaked, confidential or access-restricted material of any kind, including
  material obtained under a licence or terms of use that forbid this use.
- Proprietary file formats without a public specification. Such formats are
  listed in [FORMATS.md](FORMATS.md) as not implemented, with the reason.
- Icons, logos, colour palettes, screenshots, artwork or other visual
  assets of any other product, and trademarks or product names in the
  product.

## What is used

- **Our own reasoning and code.** All source code in this repository is
  written for OpenShape. Third-party code enters only as declared Cargo
  dependencies, recorded in `Cargo.lock`. `cargo deny check` (configuration
  in `deny.toml`, rationale in [ADR 0009](adr/0009-dependency-licensing.md))
  fails the build unless every dependency's licence is compatible with
  GPL-3.0-or-later and it comes from crates.io.
- **Published literature and standards.** Hard algorithms get a design note
  in [design/](design/) before they are implemented, and each note cites the
  textbooks, papers and standards it follows (for example Piegl & Tiller,
  *The NURBS Book*, for NURBS; Madsen, Nielsen & Tingleff and Nocedal &
  Wright for the constraint solver). File formats are implemented only from
  public specifications, cited in [FORMATS.md](FORMATS.md).
- **Open-source prior art.** A design note may name an open-source project
  whose published approach it compares with. That names prior art. It does
  not mean code was taken: no code from such projects is in OpenShape.
- **Common user-interface conventions.** The application uses the layout
  conventions shared across parametric CAD software: a feature list on the
  left with a rollback marker, a command toolbar, a graphics area with an
  orientation cube, element tabs, floating feature dialogs, context menus
  and a command search. These are functional conventions of the field, not
  any one product's design ([ADR 0012](adr/0012-ui-look.md)).

## Visual assets

- Every icon is an original drawing made for OpenShape, drawn in code
  (`crates/ui-kit/src/icons.rs`); no image files are used. Each icon's
  documentation says what it depicts.
- The colour palette is OpenShape's own, defined in one place
  (`crates/ui-kit/src/theme.rs`).
- Text and a few symbols are drawn with the fonts bundled with the egui
  toolkit (Ubuntu Light, Hack, Noto Emoji and emoji-icon-font), under their
  own open licences (Ubuntu Font Licence 1.0, MIT with the Bitstream Vera
  licence, SIL OFL 1.1 and MIT respectively). See `NOTICE`.

## Records

The record of how OpenShape was built is the git history, the architecture
decision records in [adr/](adr/), the design notes in [design/](design/),
[STATUS.md](../STATUS.md) and this policy. Contributors describe their own
decisions and sources in commit messages and documents. Anyone who finds
material in the repository that breaks this policy must report it so it can
be removed and rewritten from clean sources.
