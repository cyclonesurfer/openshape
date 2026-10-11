//! The view cube: a small cube in the graphics area's corner that turns with
//! the camera and snaps the view to a face, edge or corner direction when
//! clicked.
//!
//! The cube is drawn with egui's painter in an orthographic projection that
//! uses only the camera's rotation ([`Camera::right_dir`],
//! [`Camera::up_dir`], [`Camera::back_dir`]), so it never shows perspective
//! distortion. Its local axes are the world axes: the front face looks
//! towards -Y, the right face towards +X and the top face towards +Z,
//! matching [`openshape_render::ViewPreset`].
//!
//! Every visible face is split into a 3x3 grid: the centre cell picks the
//! face, the side cells pick the edge shared with the neighbouring face and
//! the corner cells pick the corner, so 26 directions are reachable
//! ([`CubeTarget`]).

use crate::icons::{paint_icon, Icon};
use crate::theme::Theme;
use egui::{Color32, CornerRadius, FontId, Pos2, Rect, Sense, Shape, Stroke, Ui, Vec2};
use glam::DVec3;
use openshape_render::{Camera, ViewPreset};

/// Half-width of a face's centre cell, in face coordinates where the face
/// spans `-1..=1`; the bands outside it (20 % of the face width on each
/// side) are the edge and corner cells.
pub const CENTER_HALF: f64 = 0.6;

/// Fraction of the cube's maximal on-screen radius used for its size, so it
/// never touches the corners of its rectangle (where the home button sits).
const FILL: f64 = 0.86;

/// Faces with `normal . back` below this are treated as edge-on: not drawn
/// and not pickable.
const MIN_FACING: f64 = 1e-3;

/// Slack on the face boundary when picking, so a point exactly on a cube
/// edge (such as the screen centre when looking along an edge) still hits.
const PICK_SLACK: f64 = 1e-9;

/// Labels are hidden on faces turned further than this from the viewer
/// (`normal . back` below it).
const LABEL_MIN_FACING: f64 = 0.3;

/// A face of the view cube.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum CubeFace {
    /// Faces -Y (seen from the front view).
    Front,
    /// Faces +Y.
    Back,
    /// Faces +Z.
    Top,
    /// Faces -Z.
    Bottom,
    /// Faces -X.
    Left,
    /// Faces +X.
    Right,
}

impl CubeFace {
    /// All six faces.
    pub const ALL: [CubeFace; 6] = [
        CubeFace::Front,
        CubeFace::Back,
        CubeFace::Top,
        CubeFace::Bottom,
        CubeFace::Left,
        CubeFace::Right,
    ];

    /// Outward normal as a sign triple.
    pub fn normal(self) -> [i8; 3] {
        match self {
            CubeFace::Front => [0, -1, 0],
            CubeFace::Back => [0, 1, 0],
            CubeFace::Top => [0, 0, 1],
            CubeFace::Bottom => [0, 0, -1],
            CubeFace::Left => [-1, 0, 0],
            CubeFace::Right => [1, 0, 0],
        }
    }

    /// The face's label as drawn on the cube.
    pub fn label(self) -> &'static str {
        match self {
            CubeFace::Front => "Front",
            CubeFace::Back => "Back",
            CubeFace::Top => "Top",
            CubeFace::Bottom => "Bottom",
            CubeFace::Left => "Left",
            CubeFace::Right => "Right",
        }
    }

    /// The matching standard view.
    pub fn view_preset(self) -> ViewPreset {
        match self {
            CubeFace::Front => ViewPreset::Front,
            CubeFace::Back => ViewPreset::Back,
            CubeFace::Top => ViewPreset::Top,
            CubeFace::Bottom => ViewPreset::Bottom,
            CubeFace::Left => ViewPreset::Left,
            CubeFace::Right => ViewPreset::Right,
        }
    }

    fn normal_vec(self) -> DVec3 {
        to_vec(self.normal())
    }

    /// World direction that is "up" in the face's label: +Z on the four
    /// side faces, +Y on the top face (so it reads upright from the top
    /// view) and -Y on the bottom face (likewise from the bottom view).
    fn label_up(self) -> DVec3 {
        match self {
            CubeFace::Top => DVec3::Y,
            CubeFace::Bottom => -DVec3::Y,
            _ => DVec3::Z,
        }
    }

    /// In-face axes `(u, v)`: `u = label_up x normal` (the label's right),
    /// `v = label_up`.
    fn axes(self) -> (DVec3, DVec3) {
        let up = self.label_up();
        (up.cross(self.normal_vec()), up)
    }
}

/// One of the twelve cube edges, named by its two faces.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum CubeEdge {
    /// Between front and top.
    FrontTop,
    /// Between front and bottom.
    FrontBottom,
    /// Between front and left.
    FrontLeft,
    /// Between front and right.
    FrontRight,
    /// Between back and top.
    BackTop,
    /// Between back and bottom.
    BackBottom,
    /// Between back and left.
    BackLeft,
    /// Between back and right.
    BackRight,
    /// Between top and left.
    TopLeft,
    /// Between top and right.
    TopRight,
    /// Between bottom and left.
    BottomLeft,
    /// Between bottom and right.
    BottomRight,
}

impl CubeEdge {
    /// All twelve edges.
    pub const ALL: [CubeEdge; 12] = [
        CubeEdge::FrontTop,
        CubeEdge::FrontBottom,
        CubeEdge::FrontLeft,
        CubeEdge::FrontRight,
        CubeEdge::BackTop,
        CubeEdge::BackBottom,
        CubeEdge::BackLeft,
        CubeEdge::BackRight,
        CubeEdge::TopLeft,
        CubeEdge::TopRight,
        CubeEdge::BottomLeft,
        CubeEdge::BottomRight,
    ];

    /// The two faces meeting at this edge.
    pub fn faces(self) -> [CubeFace; 2] {
        use CubeFace::*;
        match self {
            CubeEdge::FrontTop => [Front, Top],
            CubeEdge::FrontBottom => [Front, Bottom],
            CubeEdge::FrontLeft => [Front, Left],
            CubeEdge::FrontRight => [Front, Right],
            CubeEdge::BackTop => [Back, Top],
            CubeEdge::BackBottom => [Back, Bottom],
            CubeEdge::BackLeft => [Back, Left],
            CubeEdge::BackRight => [Back, Right],
            CubeEdge::TopLeft => [Top, Left],
            CubeEdge::TopRight => [Top, Right],
            CubeEdge::BottomLeft => [Bottom, Left],
            CubeEdge::BottomRight => [Bottom, Right],
        }
    }
}

/// One of the eight cube corners, named by its three faces.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum CubeCorner {
    /// Front, top and left.
    FrontTopLeft,
    /// Front, top and right: the isometric view direction.
    FrontTopRight,
    /// Front, bottom and left.
    FrontBottomLeft,
    /// Front, bottom and right.
    FrontBottomRight,
    /// Back, top and left.
    BackTopLeft,
    /// Back, top and right.
    BackTopRight,
    /// Back, bottom and left.
    BackBottomLeft,
    /// Back, bottom and right.
    BackBottomRight,
}

impl CubeCorner {
    /// All eight corners.
    pub const ALL: [CubeCorner; 8] = [
        CubeCorner::FrontTopLeft,
        CubeCorner::FrontTopRight,
        CubeCorner::FrontBottomLeft,
        CubeCorner::FrontBottomRight,
        CubeCorner::BackTopLeft,
        CubeCorner::BackTopRight,
        CubeCorner::BackBottomLeft,
        CubeCorner::BackBottomRight,
    ];

    /// The three faces meeting at this corner.
    pub fn faces(self) -> [CubeFace; 3] {
        use CubeFace::*;
        match self {
            CubeCorner::FrontTopLeft => [Front, Top, Left],
            CubeCorner::FrontTopRight => [Front, Top, Right],
            CubeCorner::FrontBottomLeft => [Front, Bottom, Left],
            CubeCorner::FrontBottomRight => [Front, Bottom, Right],
            CubeCorner::BackTopLeft => [Back, Top, Left],
            CubeCorner::BackTopRight => [Back, Top, Right],
            CubeCorner::BackBottomLeft => [Back, Bottom, Left],
            CubeCorner::BackBottomRight => [Back, Bottom, Right],
        }
    }
}

/// A view direction the cube can snap to: one of 6 faces, 12 edges or
/// 8 corners.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum CubeTarget {
    /// Look straight at a face.
    Face(CubeFace),
    /// Look at an edge, halfway (45 degrees) between its two faces.
    Edge(CubeEdge),
    /// Look at a corner along the cube's diagonal.
    Corner(CubeCorner),
}

impl CubeTarget {
    /// The isometric corner (front, top, right), as chosen by the home
    /// button; equal to [`ViewPreset::Isometric`].
    pub const HOME: CubeTarget = CubeTarget::Corner(CubeCorner::FrontTopRight);

    /// All 26 targets: faces, then edges, then corners.
    pub fn all() -> Vec<CubeTarget> {
        CubeFace::ALL
            .into_iter()
            .map(CubeTarget::Face)
            .chain(CubeEdge::ALL.into_iter().map(CubeTarget::Edge))
            .chain(CubeCorner::ALL.into_iter().map(CubeTarget::Corner))
            .collect()
    }

    /// Direction from the cube's centre towards the eye as a sign triple:
    /// the sum of the normals of the target's faces.
    pub fn direction(self) -> [i8; 3] {
        let faces: &[CubeFace] = match self {
            CubeTarget::Face(f) => &[f],
            CubeTarget::Edge(e) => &e.faces(),
            CubeTarget::Corner(c) => &c.faces(),
        };
        faces.iter().fold([0; 3], |acc, f| {
            let n = f.normal();
            [acc[0] + n[0], acc[1] + n[1], acc[2] + n[2]]
        })
    }

    /// The target with the given sign triple (each component -1, 0 or 1,
    /// not all zero).
    pub fn from_direction(d: [i8; 3]) -> Option<CubeTarget> {
        if d.iter().any(|c| !(-1..=1).contains(c)) || d == [0; 3] {
            return None;
        }
        CubeTarget::all().into_iter().find(|t| t.direction() == d)
    }

    /// Camera `(yaw, pitch)` in radians that looks along this target
    /// (see [`Camera`] for the conventions). Faces give exactly the
    /// [`ViewPreset`] angles, edges look 45 degrees between their faces, and
    /// [`CubeTarget::HOME`] gives [`ViewPreset::Isometric`].
    pub fn yaw_pitch(self) -> (f64, f64) {
        if let CubeTarget::Face(f) = self {
            return f.view_preset().yaw_pitch();
        }
        if self == CubeTarget::HOME {
            return ViewPreset::Isometric.yaw_pitch();
        }
        // back_dir = (sin yaw cos pitch, -cos yaw cos pitch, sin pitch).
        let d = to_vec(self.direction());
        let horizontal = (d.x * d.x + d.y * d.y).sqrt();
        if horizontal == 0.0 {
            // Straight up or down; same yaw as the top and bottom presets.
            return (0.0, d.z.signum() * std::f64::consts::FRAC_PI_2);
        }
        ((d.x).atan2(-d.y), (d.z / horizontal).atan())
    }
}

fn to_vec(v: [i8; 3]) -> DVec3 {
    DVec3::new(v[0] as f64, v[1] as f64, v[2] as f64)
}

/// The cube's half-size (centre to face) in points for `rect`: the largest
/// that keeps every orientation inside the circle inscribed in `rect`, less
/// a margin.
pub fn cube_half_size(rect: Rect) -> f64 {
    let half = 0.5 * rect.width().min(rect.height()) as f64;
    half * FILL / 3f64.sqrt()
}

/// Screen-space basis of the cube for `camera`: projections of the world
/// axes X, Y and Z, in points, for a cube of half-size `s`.
fn screen_axes(camera: &Camera, s: f64) -> [Vec2; 3] {
    let (r, u) = (camera.right_dir(), camera.up_dir());
    [DVec3::X, DVec3::Y, DVec3::Z]
        .map(|a| Vec2::new((a.dot(r) * s) as f32, (-(a.dot(u)) * s) as f32))
}

/// Projects the cube-local point `p` (the cube spans `-1..=1` on each axis)
/// to the screen for `camera`, with the cube centred in `rect`.
pub fn project(camera: &Camera, rect: Rect, p: [f64; 3]) -> Pos2 {
    let [ax, ay, az] = screen_axes(camera, cube_half_size(rect));
    rect.center() + ax * p[0] as f32 + ay * p[1] as f32 + az * p[2] as f32
}

/// A pick on one face: face, in-face coordinates and depth towards the eye.
struct FaceHit {
    face: CubeFace,
    a: f64,
    b: f64,
    depth: f64,
}

/// In-face coordinates `(a, b)` (each `-1..=1` on the face) of screen point
/// `pos` on the plane of `face`, or `None` if the face is edge-on.
fn face_coords(camera: &Camera, rect: Rect, face: CubeFace, pos: Pos2) -> Option<(f64, f64)> {
    let s = cube_half_size(rect);
    let (r, u) = (camera.right_dir(), camera.up_dir());
    let to_screen = |v: DVec3| (v.dot(r) * s, -v.dot(u) * s);
    let (fu, fv) = face.axes();
    let n = face.normal_vec();
    let (cx, cy) = to_screen(n);
    let (ux, uy) = to_screen(fu);
    let (vx, vy) = to_screen(fv);
    let det = ux * vy - vx * uy;
    if !(s > 0.0 && det.abs() >= 1e-9 * s * s) {
        return None;
    }
    let px = (pos.x - rect.center().x) as f64 - cx;
    let py = (pos.y - rect.center().y) as f64 - cy;
    let a = (px * vy - vx * py) / det;
    let b = (ux * py - px * uy) / det;
    Some((a, b))
}

/// Which band of a face coordinate: -1, 0 (centre) or 1.
fn band(t: f64) -> i8 {
    if t < -CENTER_HALF {
        -1
    } else if t > CENTER_HALF {
        1
    } else {
        0
    }
}

/// The cube target under `pos` for the cube centred in `rect` (see
/// [`show`]), or `None` if `pos` is not on a visible face.
pub fn pick(camera: &Camera, rect: Rect, pos: Pos2) -> Option<CubeTarget> {
    let back = camera.back_dir();
    let mut best: Option<FaceHit> = None;
    for face in CubeFace::ALL {
        let n = face.normal_vec();
        if n.dot(back) < MIN_FACING {
            continue;
        }
        let Some((a, b)) = face_coords(camera, rect, face, pos) else {
            continue;
        };
        // Written so NaN (from a non-finite pointer position) misses too.
        if !(a.abs() <= 1.0 + PICK_SLACK && b.abs() <= 1.0 + PICK_SLACK) {
            continue;
        }
        let (fu, fv) = face.axes();
        let depth = (n + fu * a + fv * b).dot(back);
        if best.as_ref().is_none_or(|h| depth > h.depth) {
            best = Some(FaceHit { face, a, b, depth });
        }
    }
    let hit = best?;
    let (fu, fv) = hit.face.axes();
    let d = hit.face.normal_vec() + fu * band(hit.a) as f64 + fv * band(hit.b) as f64;
    CubeTarget::from_direction([d.x as i8, d.y as i8, d.z as i8])
}

/// The home button's square in the top-left corner of `rect`.
fn home_rect(rect: Rect) -> Rect {
    let side = (rect.width().min(rect.height()) * 0.18).clamp(14.0, 24.0);
    Rect::from_min_size(rect.min, Vec2::splat(side))
}

/// The cell of `face` that belongs to `target` (as band indices along the
/// face's axes), if `face` is one of the target's faces.
fn target_cell(face: CubeFace, target: CubeTarget) -> Option<(i8, i8)> {
    let d = to_vec(target.direction());
    if d.dot(face.normal_vec()) < 0.5 {
        return None;
    }
    let (fu, fv) = face.axes();
    Some((d.dot(fu).round() as i8, d.dot(fv).round() as i8))
}

/// Draws the view cube centred in `rect` for `camera`, with a small home
/// button in `rect`'s top-left corner, and returns the target the user
/// clicked (the home button returns [`CubeTarget::HOME`]).
///
/// Faces turned away from the viewer are culled; visible faces are shaded
/// by a fixed light from the theme's panel colours, outlined with
/// [`Theme::border_strong`] and labelled when they face the viewer enough.
/// The cell under the pointer, on every visible face it touches, is
/// highlighted with [`Theme::preselect`].
pub fn show(ui: &mut Ui, rect: Rect, camera: &Camera, theme: &Theme) -> Option<CubeTarget> {
    // An automatic id (unique per widget position in `ui`) without moving
    // the layout cursor: the cube is usually an overlay on the viewport.
    let id = ui.auto_id_with("openshape_view_cube");
    ui.skip_ahead_auto_ids(1);
    let response = ui.interact(rect, id, Sense::click());
    let home = home_rect(rect);
    let pointer = response.hover_pos();
    let over_home = pointer.is_some_and(|p| home.contains(p));
    let hovered = if over_home {
        None
    } else {
        pointer.and_then(|p| pick(camera, rect, p))
    };
    if hovered.is_some() || over_home {
        ui.ctx().set_cursor_icon(egui::CursorIcon::PointingHand);
    }

    paint(ui.painter(), rect, camera, theme, hovered);
    paint_home(ui.painter(), home, theme, over_home);

    if response.clicked() {
        let pos = response.interact_pointer_pos()?;
        if home.contains(pos) {
            return Some(CubeTarget::HOME);
        }
        return pick(camera, rect, pos);
    }
    None
}

/// Paints the cube (without the home button) with `hovered` highlighted.
pub fn paint(
    painter: &egui::Painter,
    rect: Rect,
    camera: &Camera,
    theme: &Theme,
    hovered: Option<CubeTarget>,
) {
    let back = camera.back_dir();
    let s = cube_half_size(rect);
    // A light from the upper left, slightly towards the viewer, fixed in
    // view space so the top face reads brightest in the default views.
    let light = (camera.up_dir() * 0.75 - camera.right_dir() * 0.3 + back * 0.6).normalize();
    let edge = Stroke::new(1.0, theme.border_strong);
    let to_screen = |p: DVec3| project(camera, rect, [p.x, p.y, p.z]);

    for face in CubeFace::ALL {
        let n = face.normal_vec();
        let facing = n.dot(back);
        if facing < MIN_FACING {
            continue;
        }
        let (fu, fv) = face.axes();
        let corner = |a: f64, b: f64| to_screen(n + fu * a + fv * b);
        let quad = vec![
            corner(-1.0, -1.0),
            corner(1.0, -1.0),
            corner(1.0, 1.0),
            corner(-1.0, 1.0),
        ];
        let lit = n.dot(light).clamp(0.0, 1.0) as f32;
        let fill = lerp_color(theme.border, theme.panel, 0.35 + 0.65 * lit);
        painter.add(Shape::convex_polygon(quad.clone(), fill, Stroke::NONE));

        if let Some((i, j)) = hovered.and_then(|t| target_cell(face, t)) {
            let range = |k: i8| match k {
                -1 => (-1.0, -CENTER_HALF),
                0 => (-CENTER_HALF, CENTER_HALF),
                _ => (CENTER_HALF, 1.0),
            };
            let ((a0, a1), (b0, b1)) = (range(i), range(j));
            let cell = vec![
                corner(a0, b0),
                corner(a1, b0),
                corner(a1, b1),
                corner(a0, b1),
            ];
            painter.add(Shape::convex_polygon(
                cell,
                theme.preselect.gamma_multiply(0.6),
                Stroke::NONE,
            ));
        }

        if facing >= LABEL_MIN_FACING {
            paint_label(painter, theme, face, s, &corner);
        }
        painter.add(Shape::closed_line(quad, edge));
    }
}

/// Draws `face`'s label centred on the face, rotated with the face's
/// projected right direction and shrunk to fit the face's projected size.
fn paint_label(
    painter: &egui::Painter,
    theme: &Theme,
    face: CubeFace,
    s: f64,
    corner: &dyn Fn(f64, f64) -> Pos2,
) {
    let center = corner(0.0, 0.0);
    let right = corner(1.0, 0.0) - center;
    let up = center - corner(0.0, 1.0);
    let (right_len, up_len) = (right.length(), up.length());
    if right_len < 1.0 || up_len < 1.0 {
        return;
    }
    // Width available along the label's baseline and height across it: the
    // centre cell's projected size, less a small margin.
    let avail_w = 2.0 * CENTER_HALF as f32 * right_len * 1.25;
    let avail_h = up_len * 0.9;
    let reference = (s as f32 * 0.42).clamp(6.0, 16.0);
    let font = FontId::proportional(reference);
    let galley = painter.layout_no_wrap(face.label().to_owned(), font, theme.text_muted);
    let size = galley.size();
    let scale = (avail_w / size.x).min(avail_h / size.y).min(1.0);
    let font_size = reference * scale;
    if font_size < 6.0 {
        return;
    }
    let galley = painter.layout_no_wrap(
        face.label().to_owned(),
        FontId::proportional(font_size),
        theme.text_muted,
    );
    let size = galley.size();
    let angle = right.y.atan2(right.x);
    let rot = egui::emath::Rot2::from_angle(angle);
    let pos = center - rot * (size * 0.5);
    painter.add(egui::epaint::TextShape::new(pos, galley, theme.text_muted).with_angle(angle));
}

fn paint_home(painter: &egui::Painter, rect: Rect, theme: &Theme, hovered: bool) {
    let radius = CornerRadius::same(theme.rounding.round().clamp(0.0, 255.0) as u8);
    if hovered {
        painter.rect_filled(rect, radius, theme.hover_soft);
    }
    let ink = if hovered {
        theme.accent
    } else {
        theme.text_muted
    };
    paint_icon(painter, rect.shrink(rect.width() * 0.12), Icon::Home, ink);
}

fn lerp_color(a: Color32, b: Color32, t: f32) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    let mix = |x: u8, y: u8| (x as f32 + (y as f32 - x as f32) * t).round() as u8;
    Color32::from_rgb(mix(a.r(), b.r()), mix(a.g(), b.g()), mix(a.b(), b.b()))
}
