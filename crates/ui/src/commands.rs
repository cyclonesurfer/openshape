//! The application's commands: one list that the toolbar, menus, keyboard
//! shortcuts and the command palette all draw from, so a command has one
//! name, icon, shortcut and enabled rule wherever it appears.

use egui::{Key, KeyboardShortcut, Modifiers};
use openshape_render::DisplayMode;
use openshape_ui_kit::view_cube::{CubeFace, CubeTarget};
use openshape_ui_kit::Icon;

use crate::dialogs::Template;
use crate::planes::StdPlane;
use crate::sketch_mode::{ConstraintTool, Tool};

/// A command.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Command {
    /// Starts a new, empty document.
    NewDocument,
    /// Opens a document file.
    OpenDocument,
    /// Saves the document.
    SaveDocument,
    /// Exports every part as STL.
    ExportStl,
    /// Undo.
    Undo,
    /// Redo.
    Redo,
    /// Adds a Part Studio tab.
    NewPartStudio,
    /// Opens the command palette.
    CommandSearch,
    /// A new sketch on the selected plane (asks for one when none is).
    Sketch,
    /// A new sketch on a default plane.
    SketchOn(StdPlane),
    /// Extrude from sketch regions.
    Extrude,
    /// Revolve sketch regions about a sketch line.
    Revolve,
    /// A feature from a template dialog.
    Template(Template),
    /// Moves the rollback bar to the end of the feature list.
    RollToEnd,
    /// Fits the view to the model.
    FitView,
    /// Looks at the model from a standard direction.
    View(CubeTarget),
    /// Sets how surfaces and edges are drawn.
    Display(DisplayMode),
    /// Perspective projection.
    Perspective,
    /// Orthographic projection.
    Orthographic,
    /// Shows or hides sketches outside sketch mode.
    ToggleSketches,
    /// Shows or hides the default planes.
    TogglePlanes,
    /// Picks a sketch tool.
    Tool(Tool),
    /// Toggles construction geometry.
    Construction,
    /// Shows or hides the curvature combs of the selected curves.
    CurvatureComb,
    /// Adds a constraint to the sketch selection.
    Constrain(ConstraintTool),
    /// Dimensions the sketch selection.
    Dimension,
    /// Returns dimension labels to automatic placement: the selected
    /// dimensions', or every placed one when no dimension is selected.
    ResetLabels,
    /// Deletes the sketch selection.
    DeleteSelection,
    /// Accepts the sketch and leaves sketch mode.
    FinishSketch,
    /// Reverts everything done in sketch mode and leaves it.
    CancelSketch,
}

/// Which commands make sense right now.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    /// Working in a Part Studio's feature list.
    PartStudio,
    /// Editing a sketch.
    Sketch,
}

/// The view directions offered as commands.
const VIEWS: [(CubeTarget, &str); 7] = [
    (CubeTarget::HOME, "Isometric view"),
    (CubeTarget::Face(CubeFace::Front), "Front view"),
    (CubeTarget::Face(CubeFace::Back), "Back view"),
    (CubeTarget::Face(CubeFace::Top), "Top view"),
    (CubeTarget::Face(CubeFace::Bottom), "Bottom view"),
    (CubeTarget::Face(CubeFace::Left), "Left view"),
    (CubeTarget::Face(CubeFace::Right), "Right view"),
];

impl Command {
    /// Every command, in palette order for an empty query.
    pub fn all() -> Vec<Command> {
        let mut v = vec![
            Command::Sketch,
            Command::SketchOn(StdPlane::Top),
            Command::SketchOn(StdPlane::Front),
            Command::SketchOn(StdPlane::Right),
            Command::Extrude,
            Command::Revolve,
        ];
        v.extend(Template::ALL.map(Command::Template));
        v.extend([Command::FinishSketch, Command::CancelSketch]);
        v.extend(Tool::ALL.map(Command::Tool));
        v.push(Command::Construction);
        v.push(Command::CurvatureComb);
        v.push(Command::Dimension);
        v.push(Command::ResetLabels);
        v.extend(ConstraintTool::ALL.map(Command::Constrain));
        v.push(Command::DeleteSelection);
        v.extend([
            Command::Undo,
            Command::Redo,
            Command::RollToEnd,
            Command::NewPartStudio,
            Command::FitView,
        ]);
        v.extend(VIEWS.map(|(t, _)| Command::View(t)));
        v.extend([
            Command::Display(DisplayMode::Shaded),
            Command::Display(DisplayMode::ShadedWithEdges),
            Command::Display(DisplayMode::Wireframe),
            Command::Perspective,
            Command::Orthographic,
            Command::ToggleSketches,
            Command::TogglePlanes,
            Command::SaveDocument,
            Command::OpenDocument,
            Command::ExportStl,
            Command::NewDocument,
        ]);
        v
    }

    /// Display name.
    pub fn label(self) -> String {
        match self {
            Command::NewDocument => "New document".into(),
            Command::OpenDocument => "Open document…".into(),
            Command::SaveDocument => "Save document…".into(),
            Command::ExportStl => "Export STL…".into(),
            Command::Undo => "Undo".into(),
            Command::Redo => "Redo".into(),
            Command::NewPartStudio => "New Part Studio".into(),
            Command::CommandSearch => "Search commands".into(),
            Command::Sketch => "Sketch".into(),
            Command::SketchOn(p) => format!("Sketch on {} plane", p.label()),
            Command::Extrude => "Extrude".into(),
            Command::Revolve => "Revolve".into(),
            Command::Template(t) => t.label().into(),
            Command::RollToEnd => "Roll history to end".into(),
            Command::FitView => "Fit view".into(),
            Command::View(t) => VIEWS
                .iter()
                .find(|(v, _)| *v == t)
                .map_or("View", |(_, l)| l)
                .into(),
            Command::Display(DisplayMode::Shaded) => "Shaded".into(),
            Command::Display(DisplayMode::ShadedWithEdges) => "Shaded with edges".into(),
            Command::Display(DisplayMode::Wireframe) => "Wireframe".into(),
            Command::Perspective => "Perspective".into(),
            Command::Orthographic => "Orthographic".into(),
            Command::ToggleSketches => "Show or hide sketches".into(),
            Command::TogglePlanes => "Show or hide planes".into(),
            Command::Tool(t) => t.label().into(),
            Command::Construction => "Construction".into(),
            Command::CurvatureComb => "Curvature comb".into(),
            Command::Constrain(c) => c.label().into(),
            Command::Dimension => "Dimension".into(),
            Command::ResetLabels => "Reset dimension label positions".into(),
            Command::DeleteSelection => "Delete selection".into(),
            Command::FinishSketch => "Finish sketch".into(),
            Command::CancelSketch => "Cancel sketch".into(),
        }
    }

    /// Muted group name shown after the label in the palette.
    pub fn group(self) -> &'static str {
        match self {
            Command::NewDocument
            | Command::OpenDocument
            | Command::SaveDocument
            | Command::ExportStl
            | Command::NewPartStudio => "Document",
            Command::Undo | Command::Redo | Command::CommandSearch => "Edit",
            Command::Sketch
            | Command::SketchOn(_)
            | Command::Extrude
            | Command::Revolve
            | Command::Template(_)
            | Command::RollToEnd => "Feature",
            Command::FitView
            | Command::View(_)
            | Command::Display(_)
            | Command::Perspective
            | Command::Orthographic
            | Command::ToggleSketches
            | Command::TogglePlanes => "View",
            Command::Tool(_)
            | Command::Construction
            | Command::CurvatureComb
            | Command::Dimension
            | Command::ResetLabels
            | Command::DeleteSelection
            | Command::FinishSketch
            | Command::CancelSketch => "Sketch",
            Command::Constrain(_) => "Constraint",
        }
    }

    /// Extra words the palette matches.
    pub fn keywords(self) -> &'static [&'static str] {
        match self {
            Command::Extrude => &["pad", "boss", "pull", "prism"],
            Command::Revolve => &["turn", "lathe", "spin"],
            Command::Sketch | Command::SketchOn(_) => &["draw", "profile", "2d"],
            Command::Template(_) => &["primitive", "template"],
            Command::FitView => &["zoom", "all", "extents"],
            Command::View(_) => &["orientation", "camera", "look"],
            Command::ExportStl => &["mesh", "print", "3d print"],
            Command::SaveDocument => &["write", "file"],
            Command::OpenDocument => &["load", "file"],
            Command::Dimension => &["measure", "distance", "length", "angle", "radius"],
            Command::ResetLabels => &["automatic", "place", "move", "annotation"],
            Command::DeleteSelection => &["remove", "erase"],
            Command::FinishSketch => &["accept", "done", "exit", "ok"],
            Command::CancelSketch => &["discard", "abort", "exit"],
            Command::RollToEnd => &["rollback", "history"],
            Command::CurvatureComb => &[
                "show curvature",
                "hide curvature",
                "porcupine",
                "analysis",
                "smoothness",
                "continuity",
            ],
            Command::NewPartStudio => &["tab", "element"],
            _ => &[],
        }
    }

    /// Icon, when the command has one.
    pub fn icon(self) -> Option<Icon> {
        Some(match self {
            Command::NewDocument => Icon::New,
            Command::OpenDocument => Icon::Open,
            Command::SaveDocument => Icon::Save,
            Command::ExportStl => Icon::Export,
            Command::Undo => Icon::Undo,
            Command::Redo => Icon::Redo,
            Command::NewPartStudio => Icon::PartStudio,
            Command::CommandSearch => Icon::Search,
            Command::Sketch | Command::SketchOn(_) => Icon::Sketch,
            Command::Extrude => Icon::Extrude,
            Command::Revolve => Icon::Revolve,
            Command::Template(Template::Box) => Icon::Box,
            Command::Template(Template::Cylinder) => Icon::Cylinder,
            Command::Template(Template::Plate) => Icon::Extrude,
            Command::Template(Template::Ring | Template::Vase) => Icon::Revolve,
            Command::RollToEnd => Icon::Rollback,
            Command::FitView => Icon::ZoomFit,
            Command::View(_) => Icon::Home,
            Command::Display(DisplayMode::Shaded) => Icon::Shaded,
            Command::Display(DisplayMode::ShadedWithEdges) => Icon::ShadedEdges,
            Command::Display(DisplayMode::Wireframe) => Icon::Wireframe,
            Command::Perspective => Icon::Perspective,
            Command::Orthographic => Icon::Orthographic,
            Command::ToggleSketches => Icon::Sketch,
            Command::TogglePlanes => Icon::Plane,
            Command::Tool(t) => match t {
                Tool::Select => Icon::Select,
                Tool::Line => Icon::Line,
                Tool::Rectangle => Icon::Rectangle,
                Tool::CenterRectangle => Icon::CenterRectangle,
                Tool::AlignedRectangle => Icon::AlignedRectangle,
                Tool::Circle => Icon::Circle,
                Tool::Arc => Icon::Arc3Point,
                Tool::CenterArc => Icon::CenterArc,
                Tool::TangentArc => Icon::TangentArc,
                Tool::Slot => Icon::Slot,
                Tool::Polygon => Icon::Polygon,
                Tool::Point => Icon::Point,
                Tool::Spline => Icon::Spline,
                Tool::ControlSpline => Icon::ControlSpline,
                Tool::Trim => Icon::Trim,
                Tool::Extend => Icon::Extend,
                Tool::Split => Icon::Split,
                Tool::SketchFillet => Icon::SketchFillet,
                Tool::Mirror => Icon::Mirror,
                Tool::Offset => Icon::Offset,
                Tool::LinearPattern => Icon::LinearPattern,
                Tool::CircularPattern => Icon::CircularPattern,
                Tool::Ellipse => Icon::Ellipse,
                Tool::Text => Icon::Text,
                Tool::EllipticalArc => Icon::EllipticalArc,
                Tool::Conic => Icon::Conic,
                Tool::Move => Icon::Move,
                Tool::Rotate => Icon::Rotate,
                Tool::Scale => Icon::Scale,
            },
            Command::Construction => Icon::Construction,
            Command::CurvatureComb => Icon::CurvatureComb,
            Command::Constrain(c) => match c {
                ConstraintTool::Coincident => Icon::Coincident,
                ConstraintTool::Horizontal => Icon::Horizontal,
                ConstraintTool::Vertical => Icon::Vertical,
                ConstraintTool::Parallel => Icon::Parallel,
                ConstraintTool::Perpendicular => Icon::Perpendicular,
                ConstraintTool::Tangent => Icon::Tangent,
                ConstraintTool::Curvature => Icon::Curvature,
                ConstraintTool::Equal => Icon::Equal,
                ConstraintTool::Midpoint => Icon::Midpoint,
                ConstraintTool::Concentric => Icon::Concentric,
                ConstraintTool::Symmetric => Icon::Symmetric,
                ConstraintTool::Normal => Icon::Normal,
                ConstraintTool::Fix => Icon::Fix,
            },
            Command::Dimension | Command::ResetLabels => Icon::Dimension,
            Command::DeleteSelection => Icon::Close,
            Command::FinishSketch => Icon::FinishSketch,
            Command::CancelSketch => Icon::Close,
        })
    }

    /// Keyboard shortcut.
    pub fn shortcut(self) -> Option<KeyboardShortcut> {
        let key = |k| KeyboardShortcut::new(Modifiers::NONE, k);
        let cmd = |k| KeyboardShortcut::new(Modifiers::COMMAND, k);
        Some(match self {
            Command::OpenDocument => cmd(Key::O),
            Command::SaveDocument => cmd(Key::S),
            Command::Undo => cmd(Key::Z),
            Command::Redo => KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::Z),
            Command::CommandSearch => cmd(Key::K),
            Command::Sketch => key(Key::S),
            Command::Extrude => key(Key::E),
            Command::Revolve => key(Key::V),
            Command::FitView => key(Key::F),
            Command::Tool(t) => key(t.key()?),
            Command::Construction => key(Key::Q),
            Command::Dimension => key(Key::D),
            Command::DeleteSelection => key(Key::Delete),
            _ => return None,
        })
    }

    /// The shortcut as text, e.g. "Ctrl+Shift+Z".
    pub fn shortcut_text(self) -> Option<String> {
        let s = self.shortcut()?;
        let mut parts = Vec::new();
        if s.modifiers.command || s.modifiers.ctrl {
            parts.push("Ctrl".to_string());
        }
        if s.modifiers.alt {
            parts.push("Alt".to_string());
        }
        if s.modifiers.shift {
            parts.push("Shift".to_string());
        }
        parts.push(s.logical_key.name().to_string());
        Some(parts.join("+"))
    }

    /// The tooltip: name, then shortcut.
    pub fn tooltip(self) -> String {
        match self.shortcut_text() {
            Some(k) => format!("{}  ({k})", self.label()),
            None => self.label(),
        }
    }

    /// True if the command belongs to `mode` (others are not offered).
    pub fn available_in(self, mode: Mode) -> bool {
        match self.group() {
            "Sketch" | "Constraint" => mode == Mode::Sketch,
            "Feature" => mode == Mode::PartStudio,
            _ => true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn commands_are_unique_and_named() {
        let all = Command::all();
        let set: HashSet<Command> = all.iter().copied().collect();
        assert_eq!(set.len(), all.len());
        let labels: HashSet<String> = all.iter().map(|c| c.label()).collect();
        assert_eq!(labels.len(), all.len(), "labels are distinct");
        for c in all {
            assert!(c.icon().is_some(), "{c:?} has an icon");
        }
    }

    #[test]
    fn shortcuts_do_not_clash_within_a_mode() {
        for mode in [Mode::PartStudio, Mode::Sketch] {
            let mut seen = HashSet::new();
            for c in Command::all().into_iter().filter(|c| c.available_in(mode)) {
                if let Some(s) = c.shortcut() {
                    assert!(seen.insert(s), "{c:?} reuses {s:?} in {mode:?}");
                }
            }
        }
    }

    #[test]
    fn shortcut_text_reads_like_a_menu() {
        assert_eq!(Command::Redo.shortcut_text().unwrap(), "Ctrl+Shift+Z");
        assert_eq!(Command::Extrude.shortcut_text().unwrap(), "E");
        assert_eq!(Command::Tool(Tool::Line).tooltip(), "Line  (L)");
    }
}
