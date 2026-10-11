# ADR 0012 — Application look: familiar CAD layout, OpenShape's own visual language

Status: accepted (Phase 2)

## Context

The project owner asked for the application to look and feel as close as
possible to the browser-based parametric CAD systems its users already know,
so they can move to OpenShape without relearning the workflow. The project is
clean-room ([CLEANROOM.md](../CLEANROOM.md)): no other product's code, icons,
logos, colours, artwork or name may be used, and the record of how the look
was made must show that.

Two different things are easy to conflate:

* **Layout and workflow**: where the feature list, toolbar, graphics area
  and tabs are, how a feature is created and edited, what a right-click
  offers. These are functional conventions shared across parametric CAD
  software, and they are what makes the program familiar.
* **Visual identity (trade dress)**: a product's particular icons, logo,
  brand colours and overall artwork. These belong to their owners.

Options:

| | familiar to users | clean-room | effort |
|---|---|---|---|
| A. egui's default look and our own ad-hoc layout | no | yes | none |
| B. the field's common layout and workflow, drawn in OpenShape's own visual language | yes | yes | moderate |
| C. a close visual imitation of one product | yes | no | moderate |

## Decision

Option B. The application shell follows the common parametric-CAD layout
and workflow, with icons, colours and name that are OpenShape's own.

Layout:

* **Top bar**: application menu, document name, undo and redo, a command
  search box (also opened with Ctrl+K), display and view controls.
* **Feature toolbar** under it: icon buttons for the commands of the active
  element (sketch, extrude, revolve, templates …). In sketch mode it shows
  the sketch tools, constraints, and Finish. Every button has a tooltip with
  its name and shortcut. Only commands that are implemented appear.
* **Feature list** on the left: the default geometry (origin and the three
  standard planes), then the features in order with their status, then the
  **rollback bar**, a draggable marker after the last feature that is
  regenerated; new features are inserted at the bar. Below it, the parts.
* **Graphics area** in the middle: light gradient background, the
  **view cube** in the top-right corner, and **feature dialogs** floating at
  its top-left with accept and cancel buttons in their header.
* **Element tabs** along the bottom, with a button to add an element (only
  element kinds that are implemented are offered).
* **Context menus** on features, parts, tabs and the graphics area.
* Light, compact, crisp: thin borders, small radii, one blue accent for
  selection and primary actions, a lighter blue for preselection.

Visual language:

* **Icons** are original drawings made for OpenShape, in code
  (`crates/ui-kit/src/icons.rs`), designed on a 24-unit grid as thin line
  art. Each icon's documentation says what it depicts. No icon is copied,
  traced or imitated from any product or icon set, and no image files are
  used.
* **Colours, spacing and text sizes** are tokens in one place
  (`crates/ui-kit/src/theme.rs`, `Theme::light()`); every part of the UI,
  including sketch and selection colours, takes its colours from the theme,
  so OpenShape can be re-themed (a dark theme later) by adding a `Theme`.
* **Name and logo** are OpenShape's. No other product's name appears in the
  product; `crates/ui/tests/cleanroom.rs` enforces it.

The look lives in a separate crate, `openshape-ui-kit` (theme, icons, view
cube, offscreen egui rendering for screenshots and tests), so the
application crate only composes it.

## Consequences

* The feature list's rollback bar needs model support: a rollback position
  in the part studio, stored in the op log like any edit, which
  regeneration and feature insertion respect.
* Headless screenshots of the whole application (rendered offscreen through
  egui-wgpu) become part of the visual checks, since the container has no
  display.
* Changing the look later means editing tokens and icons in one crate, not
  hunting colours through the application.
