# ADR 0002 — UI toolkit: egui (via eframe + egui-wgpu)

Status: accepted (Phase 0)

## Context

We need a native Rust GUI on Linux (Wayland and X11) that is GPL-compatible,
hosts a wgpu 3D viewport efficiently, and can grow to a dense CAD UI (feature
tree, property panels, dialogs auto-generated from feature parameters, code
editor, command palette, docking, multi-window, HiDPI, themes).

| | egui | iced | slint |
|---|---|---|---|
| License | MIT/Apache-2.0 | MIT | GPL-3.0 / royalty-free / commercial (GPL option is compatible) |
| Model | immediate mode | Elm architecture (retained) | declarative DSL (`.slint`), retained |
| wgpu integration | first-class: `egui-wgpu` shares the device; native textures and paint callbacks | wgpu renderer with a `shader` widget | wgpu interop exists but the renderer abstraction is less direct |
| Widgets for CAD UIs | rich ecosystem: docking (`egui_dock`), tables, code editor, plots | fewer complex widgets; API churn between releases | good standard widgets; custom widgets in DSL |
| Dynamic UI (dialogs generated from scripts) | trivial in immediate mode | possible, verbose | awkward: DSL is compiled ahead of time |
| Accessibility | AccessKit | AccessKit (partial) | own accessibility layer |

## Decision

egui with eframe's wgpu backend (wayland + x11 features). The viewport is
rendered by `openshape-render` into an offscreen texture registered with
egui-wgpu, so the renderer is shared verbatim with headless tests (ADR 0007).

## Consequences

* Auto-generated feature dialogs (needed for the custom-feature language)
  fall out naturally from immediate mode.
* Very large trees/tables must be virtualized (`ScrollArea::show_rows`).
* egui pins a wgpu major version; `openshape-render` must track it
  (currently wgpu 30 with egui 0.36).
