//! OpenDatum's visual language in one place.
//!
//! - [`theme`]: the colour, spacing and type tokens ([`Theme`]) and how they
//!   map onto egui's style ([`Theme::apply`], [`font_definitions`]).
//! - [`icons`]: OpenDatum's own line-art icon set ([`Icon`], [`paint_icon`])
//!   and the square toolbar button that shows them ([`icon_button`]).
//! - [`view_cube`]: the view-orientation cube for the graphics area, with
//!   pure picking ([`view_cube::pick`]) and target angles
//!   ([`view_cube::CubeTarget::yaw_pitch`]).
//! - [`palette`]: the command palette ([`CommandPalette`] listing
//!   [`PaletteEntry`] commands with fuzzy matching,
//!   [`palette::match_score`]) and the top-bar [`search_field`] that opens
//!   it.
//! - [`snapshot`]: offscreen egui rendering ([`UiSnapshot`]) for screenshots
//!   of the app and for pixel tests, on the same backend selection as the
//!   3D renderer.
//!
//! The crate depends on egui and the render crate's camera and headless
//! device only; it knows nothing about documents or the kernel.

pub mod icons;
pub mod palette;
pub mod snapshot;
pub mod theme;
pub mod view_cube;

pub use icons::{icon_button, paint_icon, Icon};
pub use palette::{search_field, CommandPalette, PaletteEntry, PaletteOutcome};
pub use snapshot::UiSnapshot;
pub use theme::{font_definitions, Theme};
