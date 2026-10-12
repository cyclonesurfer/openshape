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
* **Mouse and keys** follow the field's common convention, the same in a
  Part Studio and in sketch mode: right-drag orbits, middle-drag (or
  Ctrl+right-drag) pans, the wheel zooms, a click selects, a double-click
  edits, a right-click without dragging opens the context menu. In sketch
  mode left-drag moves geometry or box-selects; in a Part Studio left-drag
  also orbits until parts can be box-selected. Enter accepts a dialog,
  Escape cancels it; every command's shortcut is in its tooltip and in the
  command search.
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
* **Palette: the feel, not the values.** The project owner asked for the
  colours to feel like the browser-based CAD tools users know (2026-10-11).
  That feel is a functional convention of the field: light neutral
  surfaces, dark text, thin grey borders, one blue accent for selection and
  primary actions, a lighter blue for preselection. OpenShape follows it
  with values chosen here, from that description and these criteria, never
  by sampling, eyedropping or copying another product's colours, screenshots
  or style sheets:
  * neutrals are cool greys stepped for hierarchy (chrome `#EEF0F3` behind
    white panels, headers `#F6F7F9`, borders `#D5D9DF`/`#AEB4BD`);
  * the accent `#2B6CD4` is a mid blue dark enough for white text on it and
    for blue text on white (both 5.0:1, WCAG 2.x AA), with a hover shade,
    a soft tint for selected rows and a lighter preselection blue;
  * body text is 15.8:1 on panels, secondary text 6.0:1, status colours
    for error text at least 4.5:1 and for icons at least 3:1;
  * sketch colours keep the field's meaning (constrained geometry near
    black, free geometry blue, conflicts red) and selection is shown with
    a halo so it stays distinct from free geometry.
  `crates/ui-kit/tests/theme.rs` pins these contrast ratios, so a later
  change to the palette has to keep them.
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

## Addendum (2026-10-12): the window's title bar

GNOME on Wayland leaves window decorations to each application (it offers
no server-side decorations); KDE and X11 window managers draw them
themselves. With eframe's default features off, winit fell back to a bare
dark frame with no title on GNOME (seen on the first desktop run). Options:
winit's GNOME-style (Adwaita) title bar with the title; the same without the
title text; window buttons drawn in our own top bar; or nothing. Pat chose
the GNOME-style bar with the title. `crates/ui` enables winit's
`wayland-csd-adwaita` feature (it uses no winit API directly). That bar
comes from the `sctk-adwaita` crate, which renders the title with
`ab_glyph`, whose font parser is `ttf-parser` (declared unmaintained,
RUSTSEC-2026-0192; allowed as a transitive dependency by `deny.toml`, and
used here only on the desktop's own UI font), and which asks `fc-match`,
`gsettings` and `dbus-send` for the font and the light or dark preference.
Checked under a headless Weston, which also leaves decorations to the
application.
