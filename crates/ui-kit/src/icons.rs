//! OpenShape's icon set.
//!
//! Every icon is an original drawing made for OpenShape in code: line art
//! designed on a 24-unit grid from generic metaphors (a page for "new", a
//! pencil over a plane for "sketch", a profile with an arrow for
//! "extrude", and so on), drawn with egui shapes only. There are no image
//! files and no third-party artwork. Each [`Icon`] variant's documentation
//! says what its drawing depicts.
//!
//! Drawing rules: the geometry stays inside a 2-unit padding (units 2..22),
//! strokes are 1.5 units wide and scale with the icon, line ends and sharp
//! corners are rounded, and filled areas are translucent tints of the same
//! colour. [`paint_icon`] fits the 24-unit square into any rectangle; icons
//! are designed to stay legible from 16 to 32 px. Use opaque colours: a
//! translucent colour shows where rounded joins overlap the lines.

use crate::theme::Theme;
use egui::{
    Color32, CornerRadius, Painter, Pos2, Rect, Response, Sense, Shape, Stroke, Ui, Vec2,
    WidgetInfo, WidgetType,
};

/// One icon of OpenShape's set. Draw it with [`paint_icon`] or show it in a
/// toolbar with [`icon_button`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Icon {
    // ---- App ----
    /// A page with a folded top-right corner and a plus sign on it.
    New,
    /// An open folder: the back of the folder with its tab, and the front
    /// flap tilted open.
    Open,
    /// A downward arrow dropping into an open tray.
    Save,
    /// An arrow leaving an open box towards the top right.
    Export,
    /// A U-turn arrow: a stroke curving round from the bottom and ending in
    /// an arrowhead pointing left.
    Undo,
    /// The mirror image of [`Icon::Undo`]: the arrowhead points right.
    Redo,
    /// A magnifying glass: a ring with a handle to the bottom right.
    Search,
    /// Three horizontal bars.
    Menu,
    /// A plus sign.
    Plus,
    /// A diagonal cross.
    Close,
    /// A check mark (accept).
    Check,
    /// A chevron pointing down (expanded).
    ChevronDown,
    /// A chevron pointing right (collapsed).
    ChevronRight,
    /// An almond-shaped eye outline with a round pupil (visible).
    Eye,
    /// The eye of [`Icon::Eye`] struck through by a diagonal line (hidden).
    EyeOff,
    /// A triangle with an exclamation mark.
    Warning,
    /// An octagon with an exclamation mark.
    Error,
    /// A circle with a lower-case "i" made of a dot and a bar.
    Info,
    /// A cog wheel with eight teeth around a round hub (settings).
    Gear,

    // ---- View ----
    /// Four corner brackets framing a small square (fit the model to the
    /// view).
    ZoomFit,
    /// A house: pitched roof, walls and a door (home view).
    Home,
    /// A cube in oblique projection drawn only with three tinted faces, no
    /// edges (shaded display).
    Shaded,
    /// The tinted cube of [`Icon::Shaded`] with its visible edges outlined
    /// (shaded with edges).
    ShadedEdges,
    /// All twelve edges of a cube in oblique projection, no faces
    /// (wireframe display).
    Wireframe,
    /// A trapezoid floor split into three strips whose lines converge
    /// towards the top (perspective projection).
    Perspective,
    /// A parallelogram floor split into three strips whose lines stay
    /// parallel (orthographic projection).
    Orthographic,

    // ---- Tree and tabs ----
    /// A small isometric cube inside a rounded square frame (a part
    /// studio tab).
    PartStudio,
    /// A parallelogram: a flat plane seen at an angle.
    Plane,
    /// Three axes from a common origin dot: one up, one right, one
    /// receding up and to the right.
    Origin,
    /// An isometric cube with tinted faces and outlined edges (a solid
    /// part).
    Part,
    /// A closed folder with a tab.
    Folder,
    /// A pencil whose tip rests on a flat plane.
    Sketch,
    /// A flat profile with an arrow rising from it (extrude).
    Extrude,
    /// A dashed vertical axis with a tinted rectangular profile beside its
    /// upper part and an elliptical arrow sweeping round it below (revolve).
    Revolve,
    /// An isometric cube outline: hexagon silhouette and three inner edges.
    Box,
    /// A cylinder: elliptical top, two straight sides and the front half of
    /// the bottom ellipse.
    Cylinder,
    /// A clock face whose rim is an arrow turning anticlockwise (roll the
    /// history back).
    Rollback,

    // ---- Sketch tools ----
    /// A pointer arrow (select).
    Select,
    /// A diagonal line segment with a dot at each end.
    Line,
    /// A rectangle with dots on two opposite corners.
    Rectangle,
    /// A circle with a dot at its centre.
    Circle,
    /// An arc through three dotted points.
    Arc3Point,
    /// A quarter arc with a dot at its centre and a dashed radius.
    CenterArc,
    /// A dot inside a gapped cross hair.
    Point,
    /// A smooth S-shaped curve through four dotted points.
    Spline,
    /// A smooth curve pulled towards four dotted control points joined by a
    /// dashed control polygon (a control-point spline).
    ControlSpline,
    /// A wide ellipse with a dot at its centre and a dashed semi-major
    /// axis.
    Ellipse,
    /// An arc over the top of a wide ellipse with a dot at its centre and
    /// each end, the centre joined to the right end by a dashed semi-axis.
    EllipticalArc,
    /// A full rounded arch inside a dashed triangle: dots at its two ends
    /// and a small ring at the apex where its end tangents meet (a conic
    /// arc and its control point).
    Conic,
    /// A capital A drawn in strokes inside a dashed box, with a dot at the
    /// box's bottom-left corner (sketch text in its box).
    Text,
    /// A rectangle with a dot at its centre and one corner, joined by a
    /// dashed half-diagonal.
    CenterRectangle,
    /// A rectangle tilted up to the right with dots on the two ends of its
    /// lower-left side.
    AlignedRectangle,
    /// A horizontal line running into a quarter arc that turns up to the
    /// right, with a dot where they meet.
    TangentArc,
    /// A stadium (two half-circles joined by straight sides) with a dot at
    /// each half-circle's centre.
    Slot,
    /// A regular hexagon with a dot at its centre.
    Polygon,
    /// Two crossing lines; one arm past the crossing is dashed (cut away).
    Trim,
    /// A vertical boundary line; a horizontal line reaches it, its new part
    /// dashed and ending in a solid arrowhead.
    Extend,
    /// A diagonal line cut in two by a gap, a dot at each side of the gap.
    Split,
    /// A corner rounded by an arc, with the original sharp corner dashed.
    SketchFillet,
    /// An L-shaped corner and a smaller copy nested inside it.
    Offset,
    /// A dashed vertical axis with an open angle on the left and its mirror
    /// image on the right, under an arrow pointing right.
    Mirror,
    /// Three small squares in a row above an arrow pointing right.
    LinearPattern,
    /// Five small circles spaced round a centre dot.
    CircularPattern,
    /// A cross of two lines with a solid arrowhead at each of its four
    /// ends.
    Move,
    /// Three quarters of a circle about a centre dot, ending in an open
    /// arrowhead.
    Rotate,
    /// A small square in a dashed large one sharing its lower left corner,
    /// with an arrow from the small square's corner towards the large one's.
    Scale,
    /// A dashed diagonal line (construction geometry).
    Construction,
    /// Two extension lines with a dimension line between them ending in
    /// solid arrowheads.
    Dimension,
    /// Two segments meeting at one large dot (coincident points).
    Coincident,
    /// A horizontal segment with a dot at each end.
    Horizontal,
    /// A vertical segment with a dot at each end.
    Vertical,
    /// Two parallel slanted lines.
    Parallel,
    /// A vertical line standing on a horizontal one, with a right-angle
    /// mark.
    Perpendicular,
    /// A circle with a horizontal line touching its top at a dot.
    Tangent,
    /// A line running into a curve at a dot, the curve's curvature comb
    /// rising from zero at the joint.
    Curvature,
    /// Two equal horizontal bars (an equals sign).
    Equal,
    /// A segment with end ticks, a dot at its middle and a solid triangle
    /// pointing at the dot.
    Midpoint,
    /// Two concentric circles around a centre dot.
    Concentric,
    /// Two mirror-image triangles on either side of a dashed axis.
    Symmetric,
    /// An arrow rising perpendicular from a gently curved line, with a
    /// right-angle mark.
    Normal,
    /// A point held by a triangular support on a hatched ground line
    /// (fixed).
    Fix,
    /// A check mark above a flat plane (finish the sketch).
    FinishSketch,
}

impl Icon {
    /// Every icon, in declaration order.
    pub const ALL: &'static [Icon] = &[
        Icon::New,
        Icon::Open,
        Icon::Save,
        Icon::Export,
        Icon::Undo,
        Icon::Redo,
        Icon::Search,
        Icon::Menu,
        Icon::Plus,
        Icon::Close,
        Icon::Check,
        Icon::ChevronDown,
        Icon::ChevronRight,
        Icon::Eye,
        Icon::EyeOff,
        Icon::Warning,
        Icon::Error,
        Icon::Info,
        Icon::Gear,
        Icon::ZoomFit,
        Icon::Home,
        Icon::Shaded,
        Icon::ShadedEdges,
        Icon::Wireframe,
        Icon::Perspective,
        Icon::Orthographic,
        Icon::PartStudio,
        Icon::Plane,
        Icon::Origin,
        Icon::Part,
        Icon::Folder,
        Icon::Sketch,
        Icon::Extrude,
        Icon::Revolve,
        Icon::Box,
        Icon::Cylinder,
        Icon::Rollback,
        Icon::Select,
        Icon::Line,
        Icon::Rectangle,
        Icon::Circle,
        Icon::Arc3Point,
        Icon::CenterArc,
        Icon::Point,
        Icon::Spline,
        Icon::ControlSpline,
        Icon::Ellipse,
        Icon::EllipticalArc,
        Icon::Conic,
        Icon::Text,
        Icon::CenterRectangle,
        Icon::AlignedRectangle,
        Icon::TangentArc,
        Icon::Slot,
        Icon::Polygon,
        Icon::Trim,
        Icon::Extend,
        Icon::Split,
        Icon::SketchFillet,
        Icon::Offset,
        Icon::Mirror,
        Icon::LinearPattern,
        Icon::CircularPattern,
        Icon::Move,
        Icon::Rotate,
        Icon::Scale,
        Icon::Construction,
        Icon::Dimension,
        Icon::Coincident,
        Icon::Horizontal,
        Icon::Vertical,
        Icon::Parallel,
        Icon::Perpendicular,
        Icon::Tangent,
        Icon::Curvature,
        Icon::Equal,
        Icon::Midpoint,
        Icon::Concentric,
        Icon::Symmetric,
        Icon::Normal,
        Icon::Fix,
        Icon::FinishSketch,
    ];

    /// Human-readable name, e.g. `"Center arc"` (for tooltips and lists).
    pub fn name(self) -> &'static str {
        match self {
            Icon::New => "New",
            Icon::Open => "Open",
            Icon::Save => "Save",
            Icon::Export => "Export",
            Icon::Undo => "Undo",
            Icon::Redo => "Redo",
            Icon::Search => "Search",
            Icon::Menu => "Menu",
            Icon::Plus => "Add",
            Icon::Close => "Close",
            Icon::Check => "Accept",
            Icon::ChevronDown => "Expanded",
            Icon::ChevronRight => "Collapsed",
            Icon::Eye => "Show",
            Icon::EyeOff => "Hide",
            Icon::Warning => "Warning",
            Icon::Error => "Error",
            Icon::Info => "Information",
            Icon::Gear => "Settings",
            Icon::ZoomFit => "Zoom to fit",
            Icon::Home => "Home view",
            Icon::Shaded => "Shaded",
            Icon::ShadedEdges => "Shaded with edges",
            Icon::Wireframe => "Wireframe",
            Icon::Perspective => "Perspective",
            Icon::Orthographic => "Orthographic",
            Icon::PartStudio => "Part studio",
            Icon::Plane => "Plane",
            Icon::Origin => "Origin",
            Icon::Part => "Part",
            Icon::Folder => "Folder",
            Icon::Sketch => "Sketch",
            Icon::Extrude => "Extrude",
            Icon::Revolve => "Revolve",
            Icon::Box => "Box",
            Icon::Cylinder => "Cylinder",
            Icon::Rollback => "Rollback",
            Icon::Select => "Select",
            Icon::Line => "Line",
            Icon::Rectangle => "Rectangle",
            Icon::Circle => "Circle",
            Icon::Arc3Point => "Three-point arc",
            Icon::CenterArc => "Center arc",
            Icon::Point => "Point",
            Icon::Spline => "Spline",
            Icon::ControlSpline => "Control-point spline",
            Icon::Ellipse => "Ellipse",
            Icon::EllipticalArc => "Elliptical arc",
            Icon::Conic => "Conic",
            Icon::Text => "Text",
            Icon::CenterRectangle => "Center rectangle",
            Icon::AlignedRectangle => "Aligned rectangle",
            Icon::TangentArc => "Tangent arc",
            Icon::Slot => "Slot",
            Icon::Polygon => "Polygon",
            Icon::Trim => "Trim",
            Icon::Extend => "Extend",
            Icon::Split => "Split",
            Icon::SketchFillet => "Sketch fillet",
            Icon::Offset => "Offset",
            Icon::Mirror => "Mirror",
            Icon::LinearPattern => "Linear pattern",
            Icon::CircularPattern => "Circular pattern",
            Icon::Move => "Move",
            Icon::Rotate => "Rotate",
            Icon::Scale => "Scale",
            Icon::Construction => "Construction",
            Icon::Dimension => "Dimension",
            Icon::Coincident => "Coincident",
            Icon::Horizontal => "Horizontal",
            Icon::Vertical => "Vertical",
            Icon::Parallel => "Parallel",
            Icon::Perpendicular => "Perpendicular",
            Icon::Tangent => "Tangent",
            Icon::Curvature => "Curvature",
            Icon::Equal => "Equal",
            Icon::Midpoint => "Midpoint",
            Icon::Concentric => "Concentric",
            Icon::Symmetric => "Symmetric",
            Icon::Normal => "Normal",
            Icon::Fix => "Fix",
            Icon::FinishSketch => "Finish sketch",
        }
    }
}

/// Draws `icon` in `color`, fitted into `rect`: the icon's 24-unit square
/// is centred in `rect` with the side of `rect`'s shorter dimension. Nothing
/// is drawn outside that square.
pub fn paint_icon(painter: &Painter, rect: Rect, icon: Icon, color: Color32) {
    let side = rect.width().min(rect.height());
    if !(side.is_finite() && side > 0.0) {
        return;
    }
    let scale = side / 24.0;
    let pen = Pen {
        painter,
        origin: rect.center() - Vec2::splat(side * 0.5),
        scale,
        color,
        width: 1.5 * scale,
        mirror: false,
    };
    draw(&pen, icon);
}

/// A square toolbar button ([`Theme::toolbar_button`] points) showing
/// `icon` at [`Theme::toolbar_icon`] points, with `tooltip` on hover.
///
/// States: at rest the icon is in [`Theme::text`] on no background;
/// hovered adds a [`Theme::hover_soft`] background; pressed or `selected`
/// shows [`Theme::accent_soft`] behind the icon (in [`Theme::accent`] when
/// selected); disabled draws the icon in [`Theme::text_disabled`] and never
/// reports clicks.
pub fn icon_button(
    ui: &mut Ui,
    theme: &Theme,
    icon: Icon,
    tooltip: &str,
    selected: bool,
    enabled: bool,
) -> Response {
    let sense = if enabled {
        Sense::click()
    } else {
        Sense::hover()
    };
    let (rect, response) = ui.allocate_exact_size(Vec2::splat(theme.toolbar_button), sense);
    response.widget_info(|| WidgetInfo::selected(WidgetType::Button, enabled, selected, tooltip));
    if ui.is_rect_visible(rect) {
        let radius = CornerRadius::same(theme.rounding.round().clamp(0.0, 255.0) as u8);
        let hovered = enabled && response.hovered();
        let pressed = enabled && response.is_pointer_button_down_on();
        let background = if selected || pressed {
            Some(theme.accent_soft)
        } else if hovered {
            Some(theme.hover_soft)
        } else {
            None
        };
        if let Some(fill) = background {
            ui.painter().rect_filled(rect, radius, fill);
        }
        let ink = if !enabled {
            theme.text_disabled
        } else if selected {
            theme.accent
        } else {
            theme.text
        };
        let icon_rect = Rect::from_center_size(rect.center(), Vec2::splat(theme.toolbar_icon));
        paint_icon(ui.painter(), icon_rect, icon, ink);
    }
    response.on_hover_text(tooltip)
}

// ---------------------------------------------------------------------------
// Drawing helpers. Coordinates are in icon units: (0, 0) is the top-left of
// the 24-unit square, x to the right, y down. Angles are in degrees,
// counter-clockwise as seen on screen, 0 pointing right.

type P = (f32, f32);

/// Whether an icon coordinate lies in the drawing area (units 2..22).
fn in_bounds(v: f32) -> bool {
    (2.0 - 1e-3..=22.0 + 1e-3).contains(&v)
}

/// Whether the disc or circle at `c` with radius `r` lies in the drawing
/// area.
fn disc_in_bounds(c: P, r: f32) -> bool {
    in_bounds(c.0 - r) && in_bounds(c.0 + r) && in_bounds(c.1 - r) && in_bounds(c.1 + r)
}

struct Pen<'a> {
    painter: &'a Painter,
    origin: Pos2,
    scale: f32,
    color: Color32,
    width: f32,
    mirror: bool,
}

/// Point on the ellipse centred at `c` with radii `rx`, `ry` at angle `deg`.
fn on_ellipse(c: P, rx: f32, ry: f32, deg: f32) -> P {
    let (s, co) = deg.to_radians().sin_cos();
    (c.0 + rx * co, c.1 - ry * s)
}

/// Points of an elliptical arc from `a0` to `a1` degrees (either direction).
fn ellipse_arc(c: P, rx: f32, ry: f32, a0: f32, a1: f32) -> Vec<P> {
    let n = ((a1 - a0).abs() / 8.0).ceil().max(2.0) as usize;
    (0..=n)
        .map(|i| on_ellipse(c, rx, ry, a0 + (a1 - a0) * i as f32 / n as f32))
        .collect()
}

/// Points of a circular arc from `a0` to `a1` degrees.
fn arc_pts(c: P, r: f32, a0: f32, a1: f32) -> Vec<P> {
    ellipse_arc(c, r, r, a0, a1)
}

/// Samples of the quadratic Bezier curve `p0`, `ctrl`, `p1`.
fn quad(p0: P, ctrl: P, p1: P) -> Vec<P> {
    let n = 14;
    (0..=n)
        .map(|i| {
            let t = i as f32 / n as f32;
            let (a, b, c) = ((1.0 - t) * (1.0 - t), 2.0 * t * (1.0 - t), t * t);
            (
                a * p0.0 + b * ctrl.0 + c * p1.0,
                a * p0.1 + b * ctrl.1 + c * p1.1,
            )
        })
        .collect()
}

/// A smooth curve through `pts` (uniform Catmull-Rom, end tangents from the
/// neighbouring segment).
fn catmull_rom(pts: &[P]) -> Vec<P> {
    let n = pts.len();
    let at = |i: isize| pts[i.clamp(0, n as isize - 1) as usize];
    let mut out = vec![pts[0]];
    for i in 0..n as isize - 1 {
        let (p0, p1, p2, p3) = (at(i - 1), at(i), at(i + 1), at(i + 2));
        let steps = 10;
        for s in 1..=steps {
            let t = s as f32 / steps as f32;
            let (t2, t3) = (t * t, t * t * t);
            let f = |a: f32, b: f32, c: f32, d: f32| {
                0.5 * (2.0 * b
                    + (-a + c) * t
                    + (2.0 * a - 5.0 * b + 4.0 * c - d) * t2
                    + (-a + 3.0 * b - 3.0 * c + d) * t3)
            };
            out.push((f(p0.0, p1.0, p2.0, p3.0), f(p0.1, p1.1, p2.1, p3.1)));
        }
    }
    out
}

impl Pen<'_> {
    fn p(&self, (x, y): P) -> Pos2 {
        debug_assert!(
            in_bounds(x) && in_bounds(y),
            "icon point ({x}, {y}) outside the 2..22 drawing area"
        );
        let x = if self.mirror { 24.0 - x } else { x };
        self.origin + Vec2::new(x, y) * self.scale
    }

    fn mirrored(&self) -> Pen<'_> {
        Pen {
            painter: self.painter,
            origin: self.origin,
            scale: self.scale,
            color: self.color,
            width: self.width,
            mirror: !self.mirror,
        }
    }

    /// Round dot of the stroke's width (caps and round joins).
    fn cap(&self, p: Pos2) {
        self.painter
            .add(Shape::circle_filled(p, self.width * 0.5, self.color));
    }

    /// Strokes the screen-space polyline `pts`; rounds its sharp corners,
    /// and its ends when `caps` is set.
    fn stroke_screen(&self, pts: &[Pos2], caps: bool) {
        if pts.len() < 2 {
            return;
        }
        self.painter.add(Shape::line(
            pts.to_vec(),
            Stroke::new(self.width, self.color),
        ));
        if caps {
            self.cap(pts[0]);
            self.cap(pts[pts.len() - 1]);
        }
        // Egui bevels or mitres corners; a disc on each sharp corner turns
        // that into a round join.
        let cos_limit = 40f32.to_radians().cos();
        for w in pts.windows(3) {
            let (d0, d1) = (w[1] - w[0], w[2] - w[1]);
            let (l0, l1) = (d0.length(), d1.length());
            if l0 > 0.0 && l1 > 0.0 && d0.dot(d1) < cos_limit * l0 * l1 {
                self.cap(w[1]);
            }
        }
    }

    /// Open polyline with round ends and joins.
    fn poly(&self, pts: &[P]) {
        let pts: Vec<Pos2> = pts.iter().map(|&q| self.p(q)).collect();
        self.stroke_screen(&pts, true);
    }

    /// Line segment with round ends.
    fn line(&self, a: P, b: P) {
        self.poly(&[a, b]);
    }

    /// Closed polygon outline with round joins.
    fn closed(&self, pts: &[P]) {
        if pts.len() < 3 {
            self.poly(pts);
            return;
        }
        // Stroked as an open path that starts and ends in the middle of the
        // first edge, so every corner gets a round join instead of a mitre
        // spike.
        let q: Vec<Pos2> = pts.iter().map(|&p| self.p(p)).collect();
        let mid = q[0].lerp(q[1], 0.5);
        let mut path = Vec::with_capacity(q.len() + 2);
        path.push(mid);
        path.extend_from_slice(&q[1..]);
        path.push(q[0]);
        path.push(mid);
        self.stroke_screen(&path, false);
    }

    /// Axis-aligned rectangle outline.
    fn rect(&self, x0: f32, y0: f32, x1: f32, y1: f32) {
        self.closed(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]);
    }

    /// Circle outline.
    fn circle(&self, c: P, r: f32) {
        debug_assert!(disc_in_bounds(c, r), "circle {c:?} r {r} too large");
        self.painter.add(Shape::circle_stroke(
            self.p(c),
            r * self.scale,
            Stroke::new(self.width, self.color),
        ));
    }

    /// Filled disc.
    fn dot(&self, c: P, r: f32) {
        debug_assert!(disc_in_bounds(c, r), "dot {c:?} r {r} too large");
        self.painter
            .add(Shape::circle_filled(self.p(c), r * self.scale, self.color));
    }

    /// Circular arc from `a0` to `a1` degrees.
    fn arc(&self, c: P, r: f32, a0: f32, a1: f32) {
        self.poly(&arc_pts(c, r, a0, a1));
    }

    /// Convex polygon filled with the pen colour at `alpha` (0..=1).
    fn fill(&self, pts: &[P], alpha: f32) {
        let pts: Vec<Pos2> = pts.iter().map(|&q| self.p(q)).collect();
        self.painter.add(Shape::convex_polygon(
            pts,
            self.color.gamma_multiply(alpha),
            Stroke::NONE,
        ));
    }

    /// Dashed segment: dashes of roughly `dash` units (caps included) with
    /// gaps of roughly `gap`, evenly spread so both ends are dash ends.
    fn dashed(&self, a: P, b: P, dash: f32, gap: f32) {
        let (dx, dy) = (b.0 - a.0, b.1 - a.1);
        let len = (dx * dx + dy * dy).sqrt();
        if len <= 0.0 {
            return;
        }
        // Caps add half a stroke at each end of every dash.
        let w = 1.5;
        let (dash, gap) = ((dash - w).max(0.1), gap + w);
        let n = ((len + gap) / (dash + gap)).round().max(1.0);
        let d = ((len - (n - 1.0) * gap) / n).max(0.1);
        let (ux, uy) = (dx / len, dy / len);
        for i in 0..n as usize {
            let s = i as f32 * (d + gap);
            self.line(
                (a.0 + ux * s, a.1 + uy * s),
                (a.0 + ux * (s + d), a.1 + uy * (s + d)),
            );
        }
    }

    /// Open arrowhead (a chevron) with its tip at `tip`, pointing along
    /// `deg`, with arms `len` units long.
    fn chevron(&self, tip: P, deg: f32, len: f32) {
        let back = deg + 180.0;
        let a = on_ellipse(tip, len, len, back + 42.0);
        let b = on_ellipse(tip, len, len, back - 42.0);
        self.poly(&[a, tip, b]);
    }

    /// Solid triangular arrowhead with its tip at `tip`, pointing along
    /// `deg`, `len` long and `half_width` wide at the base on each side.
    fn solid_head(&self, tip: P, deg: f32, len: f32, half_width: f32) {
        let base = on_ellipse(tip, len, len, deg + 180.0);
        let l = on_ellipse(base, half_width, half_width, deg + 90.0);
        let r = on_ellipse(base, half_width, half_width, deg - 90.0);
        self.fill(&[tip, l, r], 1.0);
    }

    /// Isometric cube centred at `c`, `r` units from centre to each corner of
    /// the hexagon; `tints` fills the top, left and right faces.
    fn iso_cube(&self, c: P, r: f32, tints: Option<[f32; 3]>) {
        let h = r * 0.866_025_4;
        let top = (c.0, c.1 - r);
        let ur = (c.0 + h, c.1 - r * 0.5);
        let lr = (c.0 + h, c.1 + r * 0.5);
        let bottom = (c.0, c.1 + r);
        let ll = (c.0 - h, c.1 + r * 0.5);
        let ul = (c.0 - h, c.1 - r * 0.5);
        if let Some([t, l, rt]) = tints {
            self.fill(&[top, ur, c, ul], t);
            self.fill(&[ul, c, bottom, ll], l);
            self.fill(&[c, ur, lr, bottom], rt);
        }
        self.closed(&[top, ur, lr, bottom, ll, ul]);
        self.line(c, ul);
        self.line(c, ur);
        self.line(c, bottom);
    }

    /// Faces of the oblique cube used by the display-mode icons: front,
    /// top, right.
    fn oblique_faces() -> [[P; 4]; 3] {
        let (x0, x1, y0, y1, ox, oy) = (3.5, 15.5, 8.5, 20.5, 5.0, -5.0);
        [
            [(x0, y0), (x1, y0), (x1, y1), (x0, y1)],
            [(x0, y0), (x0 + ox, y0 + oy), (x1 + ox, y0 + oy), (x1, y0)],
            [(x1, y0), (x1 + ox, y0 + oy), (x1 + ox, y1 + oy), (x1, y1)],
        ]
    }

    /// A flat slab-like plane along the bottom of the icon.
    fn ground_plane(&self) {
        self.closed(&[(2.5, 21.5), (5.5, 17.0), (21.5, 17.0), (18.5, 21.5)]);
    }
}

fn draw(pen: &Pen<'_>, icon: Icon) {
    match icon {
        // ---------------- App ----------------
        Icon::New => {
            pen.closed(&[
                (5.5, 2.5),
                (13.5, 2.5),
                (18.5, 7.5),
                (18.5, 21.5),
                (5.5, 21.5),
            ]);
            pen.poly(&[(13.5, 2.5), (13.5, 7.5), (18.5, 7.5)]);
            pen.line((12.0, 11.0), (12.0, 18.0));
            pen.line((8.5, 14.5), (15.5, 14.5));
        }
        Icon::Open => {
            pen.poly(&[
                (3.0, 19.5),
                (3.0, 4.5),
                (9.0, 4.5),
                (11.0, 6.5),
                (18.0, 6.5),
                (18.0, 10.0),
            ]);
            pen.closed(&[(3.0, 19.5), (6.5, 10.0), (21.5, 10.0), (18.0, 19.5)]);
        }
        Icon::Save => {
            pen.line((12.0, 3.0), (12.0, 14.5));
            pen.chevron((12.0, 14.5), -90.0, 6.0);
            pen.poly(&[(3.5, 13.0), (3.5, 20.0), (20.5, 20.0), (20.5, 13.0)]);
        }
        Icon::Export => {
            pen.poly(&[
                (11.0, 4.0),
                (4.0, 4.0),
                (4.0, 20.0),
                (20.0, 20.0),
                (20.0, 13.0),
            ]);
            pen.line((10.5, 13.5), (20.0, 4.0));
            pen.poly(&[(13.5, 4.0), (20.0, 4.0), (20.0, 10.5)]);
        }
        Icon::Undo => draw_undo(pen),
        Icon::Redo => draw_undo(&pen.mirrored()),
        Icon::Search => {
            pen.circle((10.0, 10.0), 6.5);
            pen.line((14.8, 14.8), (20.5, 20.5));
        }
        Icon::Menu => {
            for y in [6.0, 12.0, 18.0] {
                pen.line((4.0, y), (20.0, y));
            }
        }
        Icon::Plus => {
            pen.line((12.0, 5.0), (12.0, 19.0));
            pen.line((5.0, 12.0), (19.0, 12.0));
        }
        Icon::Close => {
            pen.line((6.0, 6.0), (18.0, 18.0));
            pen.line((18.0, 6.0), (6.0, 18.0));
        }
        Icon::Check => {
            pen.poly(&[(4.5, 12.5), (9.5, 17.5), (19.5, 6.5)]);
        }
        Icon::ChevronDown => {
            pen.poly(&[(6.5, 9.0), (12.0, 14.5), (17.5, 9.0)]);
        }
        Icon::ChevronRight => {
            pen.poly(&[(9.0, 6.5), (14.5, 12.0), (9.0, 17.5)]);
        }
        Icon::Eye => draw_eye(pen),
        Icon::EyeOff => {
            draw_eye(pen);
            pen.line((4.5, 19.5), (19.5, 4.5));
        }
        Icon::Warning => {
            pen.closed(&[(12.0, 3.5), (21.0, 19.5), (3.0, 19.5)]);
            pen.line((12.0, 9.0), (12.0, 13.5));
            pen.dot((12.0, 16.5), 1.1);
        }
        Icon::Error => {
            let oct: Vec<P> = (0..8)
                .map(|k| on_ellipse((12.0, 12.0), 9.4, 9.4, 22.5 + 45.0 * k as f32))
                .collect();
            pen.closed(&oct);
            pen.line((12.0, 7.5), (12.0, 12.5));
            pen.dot((12.0, 16.0), 1.1);
        }
        Icon::Info => {
            pen.circle((12.0, 12.0), 9.0);
            pen.dot((12.0, 7.8), 1.1);
            pen.line((12.0, 11.0), (12.0, 16.5));
        }
        Icon::Gear => {
            let mut teeth = Vec::with_capacity(32);
            for k in 0..8 {
                let c = 45.0 * k as f32;
                teeth.push(on_ellipse((12.0, 12.0), 7.0, 7.0, c - 19.0));
                teeth.push(on_ellipse((12.0, 12.0), 9.4, 9.4, c - 10.0));
                teeth.push(on_ellipse((12.0, 12.0), 9.4, 9.4, c + 10.0));
                teeth.push(on_ellipse((12.0, 12.0), 7.0, 7.0, c + 19.0));
            }
            pen.closed(&teeth);
            pen.circle((12.0, 12.0), 3.0);
        }

        // ---------------- View ----------------
        Icon::ZoomFit => {
            pen.poly(&[(3.0, 8.0), (3.0, 3.0), (8.0, 3.0)]);
            pen.poly(&[(16.0, 3.0), (21.0, 3.0), (21.0, 8.0)]);
            pen.poly(&[(21.0, 16.0), (21.0, 21.0), (16.0, 21.0)]);
            pen.poly(&[(8.0, 21.0), (3.0, 21.0), (3.0, 16.0)]);
            pen.rect(8.5, 8.5, 15.5, 15.5);
        }
        Icon::Home => {
            pen.poly(&[(3.0, 11.5), (12.0, 3.5), (21.0, 11.5)]);
            pen.poly(&[(5.5, 9.5), (5.5, 20.5), (18.5, 20.5), (18.5, 9.5)]);
            pen.poly(&[(10.0, 20.5), (10.0, 14.5), (14.0, 14.5), (14.0, 20.5)]);
        }
        Icon::Shaded => {
            let [front, top, right] = Pen::oblique_faces();
            pen.fill(&front, 0.55);
            pen.fill(&top, 0.25);
            pen.fill(&right, 0.85);
        }
        Icon::ShadedEdges => {
            let [front, top, right] = Pen::oblique_faces();
            pen.fill(&front, 0.3);
            pen.fill(&top, 0.12);
            pen.fill(&right, 0.5);
            pen.closed(&front);
            pen.closed(&top);
            pen.closed(&right);
        }
        Icon::Wireframe => {
            let [front, top, right] = Pen::oblique_faces();
            pen.closed(&front);
            pen.closed(&top);
            pen.closed(&right);
            // The three hidden edges meet at the back-bottom-left corner.
            let hidden = (front[3].0 + 5.0, front[3].1 - 5.0);
            pen.line(hidden, front[3]);
            pen.line(hidden, top[1]);
            pen.line(hidden, right[2]);
        }
        Icon::Perspective => {
            // Floor strips converging towards a vanishing point above.
            pen.closed(&[(2.5, 20.0), (8.5, 6.0), (15.5, 6.0), (21.5, 20.0)]);
            pen.line((8.83, 20.0), (10.83, 6.0));
            pen.line((15.17, 20.0), (13.17, 6.0));
        }
        Icon::Orthographic => {
            // The same floor strips, kept parallel.
            pen.closed(&[(2.5, 20.0), (6.5, 6.0), (21.5, 6.0), (17.5, 20.0)]);
            pen.line((7.5, 20.0), (11.5, 6.0));
            pen.line((12.5, 20.0), (16.5, 6.0));
        }

        // ---------------- Tree and tabs ----------------
        Icon::PartStudio => {
            pen.closed(&rounded_square(3.0, 3.0, 21.0, 21.0, 3.0));
            pen.iso_cube((12.0, 12.0), 5.5, None);
        }
        Icon::Plane => {
            pen.closed(&[(2.5, 18.0), (8.0, 6.0), (21.5, 6.0), (16.0, 18.0)]);
        }
        Icon::Origin => {
            let o = (8.0, 16.0);
            pen.line(o, (8.0, 3.0));
            pen.line(o, (21.0, 16.0));
            pen.line(o, (16.5, 7.5));
            pen.dot(o, 2.4);
        }
        Icon::Part => {
            pen.iso_cube((12.0, 12.0), 9.0, Some([0.15, 0.4, 0.65]));
        }
        Icon::Folder => {
            pen.closed(&[
                (3.0, 4.5),
                (9.0, 4.5),
                (11.0, 6.5),
                (21.0, 6.5),
                (21.0, 19.5),
                (3.0, 19.5),
            ]);
            pen.line((3.0, 10.0), (21.0, 10.0));
        }
        Icon::Sketch => {
            pen.ground_plane();
            draw_pencil(pen, (8.5, 19.0));
        }
        Icon::Extrude => {
            pen.ground_plane();
            pen.line((12.0, 19.25), (12.0, 3.5));
            pen.chevron((12.0, 3.5), 90.0, 5.5);
        }
        Icon::Revolve => {
            pen.dashed((12.0, 2.5), (12.0, 21.5), 2.6, 1.8);
            let profile = [(14.5, 3.5), (19.5, 3.5), (19.5, 10.0), (14.5, 10.0)];
            pen.fill(&profile, 0.3);
            pen.closed(&profile);
            // Sweep from the back-left round the front to the back-right,
            // leaving a gap at the axis.
            let pts = ellipse_arc((12.0, 16.0), 8.5, 3.5, 102.0, 420.0);
            let end = pts[pts.len() - 1];
            pen.poly(&pts);
            pen.chevron(end, 167.0, 3.5);
        }
        Icon::Box => {
            pen.iso_cube((12.0, 12.0), 9.0, None);
        }
        Icon::Cylinder => {
            let mut top = ellipse_arc((12.0, 6.5), 7.5, 3.0, 0.0, 360.0);
            top.pop();
            pen.closed(&top);
            pen.line((4.5, 6.5), (4.5, 17.5));
            pen.line((19.5, 6.5), (19.5, 17.5));
            pen.poly(&ellipse_arc((12.0, 17.5), 7.5, 3.0, 180.0, 360.0));
        }
        Icon::Rollback => {
            let pts = arc_pts((12.0, 12.0), 8.0, -50.0, 215.0);
            let end = pts[pts.len() - 1];
            pen.poly(&pts);
            // Tangent of an anticlockwise arc at 215 degrees points along
            // 305 degrees.
            pen.chevron(end, 305.0, 3.5);
            pen.poly(&[(12.0, 7.5), (12.0, 12.0), (15.0, 14.0)]);
        }

        // ---------------- Sketch tools ----------------
        Icon::Select => {
            pen.closed(&[
                (6.5, 3.0),
                (6.5, 19.0),
                (10.3, 15.3),
                (13.0, 21.0),
                (15.8, 19.8),
                (13.1, 14.2),
                (18.3, 14.2),
            ]);
        }
        Icon::Line => {
            pen.line((5.0, 19.0), (19.0, 5.0));
            pen.dot((5.0, 19.0), 2.0);
            pen.dot((19.0, 5.0), 2.0);
        }
        Icon::Rectangle => {
            pen.rect(4.0, 6.5, 20.0, 17.5);
            pen.dot((4.0, 17.5), 2.0);
            pen.dot((20.0, 6.5), 2.0);
        }
        Icon::Circle => {
            pen.circle((12.0, 12.0), 8.5);
            pen.dot((12.0, 12.0), 1.6);
        }
        Icon::Arc3Point => {
            // Circle through (4, 17), (12, 8.5) and (20, 17).
            let (c, r) = ((12.0, 16.51), 8.015);
            pen.arc(c, r, 183.5, -3.5);
            for q in [(4.0, 17.0), (12.0, 8.5), (20.0, 17.0)] {
                pen.dot(q, 1.9);
            }
        }
        Icon::CenterArc => {
            let c = (5.0, 19.0);
            pen.arc(c, 14.0, 0.0, 90.0);
            pen.dashed(c, (19.0, 19.0), 2.8, 1.6);
            pen.dot(c, 2.0);
        }
        Icon::Point => {
            pen.dot((12.0, 12.0), 2.6);
            pen.line((12.0, 3.5), (12.0, 7.5));
            pen.line((12.0, 16.5), (12.0, 20.5));
            pen.line((3.5, 12.0), (7.5, 12.0));
            pen.line((16.5, 12.0), (20.5, 12.0));
        }
        Icon::Spline => {
            let knots = [(4.0, 16.0), (8.5, 7.5), (15.5, 16.5), (20.0, 8.0)];
            pen.poly(&catmull_rom(&knots));
            for q in knots {
                pen.dot(q, 1.7);
            }
        }
        Icon::ControlSpline => {
            let ctrl = [(4.0, 18.0), (6.5, 5.0), (17.5, 19.0), (20.0, 6.0)];
            for w in ctrl.windows(2) {
                pen.dashed(w[0], w[1], 2.4, 1.8);
            }
            // The cubic Bézier the four points control.
            let curve: Vec<P> = (0..=24)
                .map(|i| {
                    let t = i as f32 / 24.0;
                    let s = 1.0 - t;
                    let b = [s * s * s, 3.0 * s * s * t, 3.0 * s * t * t, t * t * t];
                    let x = (0..4).map(|k| b[k] * ctrl[k].0).sum();
                    let y = (0..4).map(|k| b[k] * ctrl[k].1).sum();
                    (x, y)
                })
                .collect();
            pen.poly(&curve);
            for q in ctrl {
                pen.dot(q, 1.7);
            }
        }
        Icon::Ellipse => {
            let c = (12.0, 12.0);
            let mut outline = ellipse_arc(c, 9.5, 6.0, 0.0, 360.0);
            outline.pop();
            pen.closed(&outline);
            pen.dashed(c, (21.5, 12.0), 2.6, 1.6);
            pen.dot(c, 1.6);
        }
        Icon::EllipticalArc => {
            let (c, rx, ry) = ((12.0, 15.0), 8.0, 6.5);
            pen.poly(&ellipse_arc(c, rx, ry, 0.0, 200.0));
            pen.dashed(c, (20.0, 15.0), 2.6, 1.6);
            pen.dot(c, 1.6);
            for a in [0.0, 200.0] {
                pen.dot(on_ellipse(c, rx, ry, a), 1.9);
            }
        }
        Icon::Conic => {
            let (a, k, b) = ((4.5, 19.0), (12.0, 4.5), (19.5, 19.0));
            pen.dashed(a, k, 2.4, 1.8);
            pen.dashed(k, b, 2.4, 1.8);
            // The rational quadratic with weight 2 (ρ = ⅔): a full arch.
            let w = 2.0;
            let curve: Vec<P> = (0..=24)
                .map(|i| {
                    let t = i as f32 / 24.0;
                    let s = 1.0 - t;
                    let r = [s * s, 2.0 * s * t * w, t * t];
                    let d = r[0] + r[1] + r[2];
                    (
                        (r[0] * a.0 + r[1] * k.0 + r[2] * b.0) / d,
                        (r[0] * a.1 + r[1] * k.1 + r[2] * b.1) / d,
                    )
                })
                .collect();
            pen.poly(&curve);
            pen.dot(a, 1.9);
            pen.dot(b, 1.9);
            pen.circle(k, 2.0);
        }
        Icon::Text => {
            let (l, t, r, b) = (4.0, 4.5, 20.0, 19.0);
            pen.dashed((l, b), (r, b), 2.4, 1.6);
            pen.dashed((r, b), (r, t), 2.4, 1.6);
            pen.dashed((r, t), (l, t), 2.4, 1.6);
            pen.dashed((l, t), (l, b), 2.4, 1.6);
            pen.poly(&[(7.5, 16.5), (12.0, 7.5), (16.5, 16.5)]);
            pen.line((9.25, 13.0), (14.75, 13.0));
            pen.dot((l, b), 2.0);
        }
        Icon::CenterRectangle => {
            pen.rect(4.0, 6.5, 20.0, 17.5);
            pen.dashed((12.0, 12.0), (20.0, 6.5), 2.6, 1.6);
            pen.dot((12.0, 12.0), 2.0);
            pen.dot((20.0, 6.5), 2.0);
        }
        Icon::AlignedRectangle => {
            // Side (12, −6) and its perpendicular (3.6, 7.2).
            let (a, b) = ((4.0, 13.5), (16.0, 7.5));
            pen.closed(&[a, b, (19.6, 14.7), (7.6, 20.7)]);
            pen.dot(a, 2.0);
            pen.dot(b, 2.0);
        }
        Icon::TangentArc => {
            pen.line((3.0, 17.0), (11.0, 17.0));
            pen.arc((11.0, 9.0), 8.0, -90.0, 30.0);
            pen.dot((11.0, 17.0), 2.0);
        }
        Icon::Slot => {
            let mut pts = arc_pts((16.0, 12.0), 5.0, -90.0, 90.0);
            pts.extend(arc_pts((8.0, 12.0), 5.0, 90.0, 270.0));
            pts.push((16.0, 17.0));
            pen.poly(&pts);
            pen.dot((8.0, 12.0), 1.6);
            pen.dot((16.0, 12.0), 1.6);
        }
        Icon::Polygon => {
            let pts: Vec<P> = (0..6)
                .map(|k| on_ellipse((12.0, 12.0), 9.0, 9.0, 60.0 * k as f32))
                .collect();
            pen.closed(&pts);
            pen.dot((12.0, 12.0), 1.6);
        }
        Icon::Trim => {
            pen.line((11.0, 3.0), (11.0, 21.0));
            pen.line((3.0, 12.0), (11.0, 12.0));
            pen.dashed((14.0, 12.0), (21.0, 12.0), 2.6, 1.8);
        }
        Icon::Extend => {
            pen.line((19.5, 3.0), (19.5, 21.0));
            pen.line((3.0, 12.0), (10.0, 12.0));
            pen.dashed((10.0, 12.0), (13.5, 12.0), 2.0, 1.5);
            pen.solid_head((18.5, 12.0), 0.0, 4.5, 2.6);
        }
        Icon::Split => {
            pen.line((3.0, 19.0), (9.5, 14.5));
            pen.line((14.5, 9.5), (21.0, 5.0));
            pen.dot((9.5, 14.5), 2.0);
            pen.dot((14.5, 9.5), 2.0);
        }
        Icon::SketchFillet => {
            pen.line((4.5, 21.0), (4.5, 12.0));
            pen.arc((12.5, 12.0), 8.0, 180.0, 90.0);
            pen.line((12.5, 4.0), (21.0, 4.0));
            pen.dashed((4.5, 9.0), (4.5, 4.0), 2.0, 1.5);
            pen.dashed((4.5, 4.0), (9.5, 4.0), 2.0, 1.5);
        }
        Icon::Offset => {
            pen.poly(&[(3.5, 21.0), (3.5, 3.5), (21.0, 3.5)]);
            pen.poly(&[(9.0, 21.0), (9.0, 9.0), (21.0, 9.0)]);
        }
        Icon::Mirror => {
            pen.dashed((12.0, 7.5), (12.0, 21.5), 3.0, 1.8);
            pen.poly(&[(9.5, 9.0), (4.0, 14.5), (9.5, 20.0)]);
            pen.poly(&[(14.5, 9.0), (20.0, 14.5), (14.5, 20.0)]);
            pen.line((5.5, 4.5), (18.0, 4.5));
            pen.chevron((18.5, 4.5), 0.0, 3.0);
        }
        Icon::LinearPattern => {
            for x in [3.0, 9.5, 16.0] {
                pen.rect(x, 5.0, x + 5.0, 10.0);
            }
            pen.line((3.0, 17.0), (17.0, 17.0));
            pen.solid_head((21.0, 17.0), 0.0, 4.5, 2.6);
        }
        Icon::CircularPattern => {
            for k in 0..5 {
                let c = on_ellipse((12.0, 12.0), 7.5, 7.5, 90.0 + 72.0 * k as f32);
                pen.circle(c, 2.0);
            }
            pen.dot((12.0, 12.0), 1.6);
        }
        Icon::Move => {
            pen.line((12.0, 6.5), (12.0, 17.5));
            pen.line((6.5, 12.0), (17.5, 12.0));
            pen.solid_head((12.0, 2.5), 90.0, 4.5, 2.8);
            pen.solid_head((12.0, 21.5), 270.0, 4.5, 2.8);
            pen.solid_head((2.5, 12.0), 180.0, 4.5, 2.8);
            pen.solid_head((21.5, 12.0), 0.0, 4.5, 2.8);
        }
        Icon::Rotate => {
            // Counter-clockwise from the top round to just below the
            // right, where the arrowhead points along the circle.
            let (c, r) = ((12.0, 12.0), 7.5);
            pen.arc(c, r, 100.0, 350.0);
            pen.chevron(on_ellipse(c, r, r, 350.0), 80.0, 4.0);
            pen.dot(c, 1.6);
        }
        Icon::Scale => {
            pen.rect(3.5, 13.5, 10.5, 20.5);
            pen.dashed((3.5, 13.5), (3.5, 3.5), 2.6, 1.8);
            pen.dashed((3.5, 3.5), (20.5, 3.5), 2.6, 1.8);
            pen.dashed((20.5, 3.5), (20.5, 20.5), 2.6, 1.8);
            pen.dashed((10.5, 20.5), (20.5, 20.5), 2.6, 1.8);
            pen.line((10.5, 13.5), (14.5, 9.5));
            pen.solid_head((17.0, 7.0), 45.0, 4.5, 2.6);
        }
        Icon::Construction => {
            pen.dashed((4.0, 20.0), (20.0, 4.0), 4.0, 3.2);
        }
        Icon::Dimension => {
            pen.line((4.5, 20.5), (4.5, 5.0));
            pen.line((19.5, 20.5), (19.5, 5.0));
            pen.line((9.0, 11.0), (15.0, 11.0));
            pen.solid_head((5.5, 11.0), 180.0, 4.5, 2.6);
            pen.solid_head((18.5, 11.0), 0.0, 4.5, 2.6);
        }
        Icon::Coincident => {
            let p = (12.0, 9.0);
            pen.line((3.5, 20.0), p);
            pen.line((20.5, 20.0), p);
            pen.dot(p, 3.2);
        }
        Icon::Horizontal => {
            pen.line((4.5, 12.0), (19.5, 12.0));
            pen.dot((4.5, 12.0), 2.3);
            pen.dot((19.5, 12.0), 2.3);
        }
        Icon::Vertical => {
            pen.line((12.0, 4.5), (12.0, 19.5));
            pen.dot((12.0, 4.5), 2.3);
            pen.dot((12.0, 19.5), 2.3);
        }
        Icon::Parallel => {
            pen.line((4.0, 17.0), (11.5, 4.0));
            pen.line((12.5, 20.0), (20.0, 7.0));
        }
        Icon::Perpendicular => {
            pen.line((11.0, 3.5), (11.0, 20.0));
            pen.line((3.5, 20.0), (20.5, 20.0));
            pen.poly(&[(11.0, 15.0), (16.0, 15.0), (16.0, 20.0)]);
        }
        Icon::Tangent => {
            pen.circle((12.0, 14.5), 6.5);
            pen.line((3.0, 8.0), (21.0, 8.0));
            pen.dot((12.0, 8.0), 2.0);
        }
        Icon::Curvature => {
            // A cubic leaving a line curvature continuously (its first
            // three control points on the line), with its curvature comb:
            // normals to the outside, as long as the curvature, their tips
            // joined. The comb starts at zero at the joint, as the line's.
            let ctrl = [(7.0, 9.0), (10.5, 9.0), (14.0, 9.0), (20.0, 20.0)];
            let bez = |t: f32| {
                let s = 1.0 - t;
                let b = [s * s * s, 3.0 * s * s * t, 3.0 * s * t * t, t * t * t];
                let d1 = [
                    -3.0 * s * s,
                    3.0 * s * (s - 2.0 * t),
                    3.0 * t * (2.0 * s - t),
                    3.0 * t * t,
                ];
                let d2 = [
                    6.0 * s,
                    6.0 * (3.0 * t - 2.0),
                    6.0 * (1.0 - 3.0 * t),
                    6.0 * t,
                ];
                let at = |w: [f32; 4]| {
                    let x: f32 = (0..4).map(|k| w[k] * ctrl[k].0).sum();
                    let y: f32 = (0..4).map(|k| w[k] * ctrl[k].1).sum();
                    (x, y)
                };
                (at(b), at(d1), at(d2))
            };
            pen.line((2.5, 9.0), ctrl[0]);
            let curve: Vec<P> = (0..=24).map(|i| bez(i as f32 / 24.0).0).collect();
            pen.poly(&curve);
            let mut tips = vec![ctrl[0]];
            for i in 1..=5 {
                let (p, d1, d2) = bez(i as f32 / 5.0);
                let n = d1.0.hypot(d1.1);
                let k = (d1.0 * d2.1 - d1.1 * d2.0) / (n * n * n);
                // The left normal turned to the side away from the centre
                // of curvature.
                let (nx, ny) = (-d1.1 / n, d1.0 / n);
                let len = (k.abs() * 45.0).min(5.0) * -k.signum();
                let tip = (p.0 + nx * len, p.1 + ny * len);
                pen.line(p, tip);
                tips.push(tip);
            }
            pen.poly(&tips);
            pen.dot(ctrl[0], 1.8);
        }
        Icon::Equal => {
            pen.line((5.0, 9.0), (19.0, 9.0));
            pen.line((5.0, 15.0), (19.0, 15.0));
        }
        Icon::Midpoint => {
            pen.line((3.5, 16.0), (20.5, 16.0));
            pen.line((3.5, 13.0), (3.5, 19.0));
            pen.line((20.5, 13.0), (20.5, 19.0));
            pen.dot((12.0, 16.0), 2.5);
            pen.fill(&[(8.5, 6.0), (15.5, 6.0), (12.0, 11.0)], 1.0);
        }
        Icon::Concentric => {
            pen.circle((12.0, 12.0), 8.5);
            pen.circle((12.0, 12.0), 4.5);
            pen.dot((12.0, 12.0), 1.4);
        }
        Icon::Symmetric => {
            pen.dashed((12.0, 2.5), (12.0, 21.5), 3.0, 1.8);
            pen.closed(&[(9.5, 6.5), (9.5, 17.5), (3.0, 17.5)]);
            pen.closed(&[(14.5, 6.5), (14.5, 17.5), (21.0, 17.5)]);
        }
        Icon::Normal => {
            // An arc whose top is at (12, 18).
            pen.arc((12.0, 28.0), 10.0, 37.0, 143.0);
            pen.line((12.0, 18.0), (12.0, 4.0));
            pen.chevron((12.0, 4.0), 90.0, 5.0);
            pen.poly(&[(12.0, 14.5), (15.5, 14.5), (15.5, 18.6)]);
        }
        Icon::Fix => {
            pen.closed(&[(12.0, 6.0), (17.5, 14.0), (6.5, 14.0)]);
            pen.line((3.5, 14.0), (20.5, 14.0));
            for x in [7.0, 11.0, 15.0, 19.0] {
                pen.line((x, 14.0), (x - 3.5, 19.5));
            }
            pen.dot((12.0, 6.0), 2.4);
        }
        Icon::FinishSketch => {
            pen.ground_plane();
            pen.poly(&[(5.5, 9.0), (10.0, 13.5), (18.5, 3.5)]);
        }
    }
}

fn draw_undo(pen: &Pen<'_>) {
    let mut pts = vec![(5.5, 9.0), (14.0, 9.0)];
    pts.extend(arc_pts((14.0, 14.0), 5.0, 90.0, -90.0).into_iter().skip(1));
    pts.push((9.0, 19.0));
    pen.poly(&pts);
    pen.chevron((5.0, 9.0), 180.0, 5.0);
}

fn draw_eye(pen: &Pen<'_>) {
    let mut outline = quad((2.5, 12.0), (12.0, 2.5), (21.5, 12.0));
    let lower = quad((21.5, 12.0), (12.0, 21.5), (2.5, 12.0));
    outline.extend(&lower[1..lower.len() - 1]);
    pen.closed(&outline);
    pen.circle((12.0, 12.0), 3.0);
}

/// A pencil pointing down-left with its tip at `tip`.
fn draw_pencil(pen: &Pen<'_>, tip: P) {
    use std::f32::consts::FRAC_1_SQRT_2 as H;
    let d = (H, -H);
    let n = (H * 2.3, H * 2.3);
    let along = |t: f32| (tip.0 + d.0 * t, tip.1 + d.1 * t);
    let off = |p: P, s: f32| (p.0 + n.0 * s, p.1 + n.1 * s);
    let (cone, end) = (along(4.0), along(15.0));
    pen.closed(&[
        tip,
        off(cone, 1.0),
        off(end, 1.0),
        off(end, -1.0),
        off(cone, -1.0),
    ]);
    pen.line(off(cone, 1.0), off(cone, -1.0));
    let band = along(12.0);
    pen.line(off(band, 1.0), off(band, -1.0));
}

/// Outline points of a square with rounded corners.
fn rounded_square(x0: f32, y0: f32, x1: f32, y1: f32, r: f32) -> Vec<P> {
    let mut pts = Vec::new();
    pts.extend(arc_pts((x1 - r, y0 + r), r, 90.0, 0.0));
    pts.extend(arc_pts((x1 - r, y1 - r), r, 0.0, -90.0));
    pts.extend(arc_pts((x0 + r, y1 - r), r, -90.0, -180.0));
    pts.extend(arc_pts((x0 + r, y0 + r), r, 180.0, 90.0));
    pts
}
