//! Sketch offset (design note 0007 §10), driven through the editor: chains
//! of lines and arcs and circles offset to the pointer's side, the distance
//! typed afterwards as a dimension, the images following their sources with
//! no degrees of freedom or redundancy added, refusals of offsets that
//! would turn back or cross, and a fuzz checked against the kernel's own
//! curve intersections.

mod common;

use common::*;
use opendatum_app::sketch_mode::{Snap, SnapTarget, Tool};
use opendatum_app::sketch_ops::offset;
use opendatum_app::sketch_ops::trim::Preview;
use opendatum_app::sketch_ops::{shapes, Batch};
use opendatum_app::Model;
use opendatum_core_model::sketch::{Constraint, Geometry};
use opendatum_core_model::{DVec2, FeatureId, PlaneRef, SketchData, SketchId};
use opendatum_kernel::{Kernel, ProfileCurve};

use opendatum_sketch_solver::{ResolvedCurve, SolveStatus};
use proptest::prelude::*;

fn redundant(m: &Model, f: FeatureId) -> usize {
    m.sketch_result(f).unwrap().solution.redundant.len()
}

fn solved(m: &Model, f: FeatureId) -> SketchData {
    m.sketch_result(f).unwrap().solution.sketch.clone()
}

/// An offset constraint: its key, distance, pieces and smooth points.
type OffsetEntry = (
    SketchId,
    f64,
    Vec<(SketchId, SketchId, bool)>,
    Vec<SketchId>,
);

/// The offset constraints of a sketch.
fn offsets(m: &Model, f: FeatureId) -> Vec<OffsetEntry> {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .constraints
        .iter()
        .filter_map(|c| match &c.constraint {
            Constraint::Offset {
                distance,
                pieces,
                smooth,
            } => Some((c.id, *distance, pieces.clone(), smooth.clone())),
            _ => None,
        })
        .collect()
}

/// How far the offset's images are from their sources' offset curves, from
/// the solved geometry with formulas of its own: each image end at the
/// distance from its source's line or circle, on the given side.
fn off_by(s: &SketchData, d: f64, pieces: &[(SketchId, SketchId, bool)]) -> f64 {
    let p = |k: SketchId| s.point(k).unwrap();
    let mut worst = 0.0f64;
    for &(src, img, left) in pieces {
        let side = if left { 1.0 } else { -1.0 };
        match (
            &s.entity(src).unwrap().geometry,
            &s.entity(img).unwrap().geometry,
        ) {
            (Geometry::Line { start, end }, Geometry::Line { start: a, end: b }) => {
                let u = (p(*end) - p(*start)).normalize();
                for x in [*a, *b] {
                    let dist = u.perp_dot(p(x) - p(*start));
                    worst = worst.max((dist - side * d).abs());
                }
                // Same direction.
                assert!((p(*b) - p(*a)).dot(u) > 0.0, "an image line turned back");
            }
            (
                Geometry::Arc { center, start, .. },
                Geometry::Arc {
                    center: c2,
                    start: a,
                    end: b,
                },
            ) => {
                assert_eq!(center, c2, "an image arc shares its source's centre");
                let r = (p(*start) - p(*center)).length();
                for x in [*a, *b] {
                    worst = worst.max(((p(x) - p(*center)).length() - (r - side * d)).abs());
                }
            }
            (
                Geometry::Circle { center, radius },
                Geometry::Circle {
                    center: c2,
                    radius: r2,
                },
            ) => {
                assert_eq!(center, c2);
                worst = worst.max((r2 - (radius - side * d)).abs());
            }
            other => panic!("offset pair of different kinds: {other:?}"),
        }
    }
    worst
}

/// The kernel's profile curve of a solved sketch curve.
fn profile(c: &ResolvedCurve) -> ProfileCurve {
    opendatum_features::kernel_curve(c).expect("a valid curve")
}

/// True if two of the images cross or touch anywhere but at an end point
/// they share, by the kernel's curve intersection (independent of the
/// solver's own check).
fn images_cross(m: &Model, f: FeatureId, pieces: &[(SketchId, SketchId, bool)]) -> bool {
    let s = solved(m, f);
    let all = curves(m, f, true);
    let shape = |k: SketchId| all.iter().find(|(id, _)| *id == k).unwrap().1.clone();
    let ends = |k: SketchId| match s.entity(k).unwrap().geometry {
        Geometry::Line { start, end } | Geometry::Arc { start, end, .. } => vec![start, end],
        _ => Vec::new(),
    };
    for (i, a) in pieces.iter().enumerate() {
        for b in &pieces[i + 1..] {
            let shared: Vec<DVec2> = ends(a.1)
                .into_iter()
                .filter(|p| ends(b.1).contains(p))
                .map(|p| s.point(p).unwrap())
                .collect();
            let contacts = m
                .kernel()
                .curve_contacts(&profile(&shape(a.1)), &profile(&shape(b.1)))
                .unwrap();
            if contacts
                .iter()
                .any(|c| shared.iter().all(|x| (c.point - *x).length() > 1e-6))
            {
                return true;
            }
        }
    }
    false
}

#[test]
fn a_rectangle_is_offset_to_the_pointer_and_the_distance_typed() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Rectangle,
        &[(-20.0, -10.0), (20.0, 10.0)],
    );
    let (before, lines) = (
        dof(&m, f),
        count(&m, f, |g| matches!(g, Geometry::Line { .. })),
    );
    assert_eq!(lines, 4);
    ed.set_tool(Tool::Offset);
    assert_eq!(
        ed.hint(),
        "Click a line, arc or circle to offset it with the curves joined to it"
    );
    // Hovering a side shows its chain: the whole rectangle.
    hover(&m, &mut ed, &map, 0.0, -10.0);
    let plan = ed
        .curve_edit
        .clone()
        .expect("a curve")
        .plan
        .expect("a chain");
    let Preview::Shapes { curves: shown, .. } = plan.preview else {
        panic!()
    };
    assert_eq!(shown.len(), 4);
    click(&mut m, &mut ed, &map, 0.0, -10.0);
    assert_eq!(ed.message.as_deref(), Some("Chain of 4 curves chosen"));
    // 5 mm below the bottom side: outside.
    hover(&m, &mut ed, &map, 3.0, -15.0);
    let plan = ed.curve_edit.clone().unwrap().plan.expect("an offset");
    let Preview::Shapes { curves: shown, .. } = plan.preview else {
        panic!()
    };
    assert!(shown
        .iter()
        .any(|c| matches!(c, ResolvedCurve::Line { start, end }
        if (start.y + 15.0 * MM).abs() < 1e-9 && (end.y + 15.0 * MM).abs() < 1e-9
            && ((start.x - end.x).abs() - 50.0 * MM).abs() < 1e-9)));
    click(&mut m, &mut ed, &map, 3.0, -15.0);
    assert_eq!(ed.message.as_deref(), Some("Offset 4 curves"));
    // The distance is open for typing.
    let prompt = ed.prompt.as_mut().expect("the distance prompt");
    assert_eq!(prompt.text, "5");
    prompt.text = "3".into();
    ed.commit_prompt(&mut m);
    let o = offsets(&m, f);
    assert_eq!(o.len(), 1);
    assert!((o[0].1 - 3.0 * MM).abs() < 1e-15);
    assert!(o[0].3.is_empty(), "a rectangle has no smooth joints");
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), before);
    assert!(off_by(&solved(&m, f), o[0].1, &o[0].2) < 1e-12);
    // A 46 × 26 mm outline around the 40 × 20 mm rectangle.
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 2);
    assert!((areas[0] - 800.0 * MM * MM).abs() < 1e-12, "{areas:?}");
    assert!((areas[1] - 396.0 * MM * MM).abs() < 1e-12, "{areas:?}");
    let x = extrude(&mut m, &map, &[(0.0, -11.5)], "10 mm");
    assert!((volume(&m, x) - 396.0 * MM * MM * 10.0 * MM).abs() < 1e-15);
    // The tool is ready for the next chain.
    assert_eq!(ed.tool, Tool::Offset);
    assert!(ed.offset_chain.is_none());
}

#[test]
fn the_images_follow_their_sources() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Rectangle,
        &[(-20.0, -10.0), (20.0, 10.0)],
    );
    ed.set_tool(Tool::Offset);
    click(&mut m, &mut ed, &map, -20.0, 0.0);
    click(&mut m, &mut ed, &map, -16.0, 0.0);
    ed.prompt = None;
    let o = offsets(&m, f);
    assert!((o[0].1 - 4.0 * MM).abs() < 1e-15, "inside, 4 mm");
    // Dragging a corner of the rectangle: the inner outline follows.
    ed.set_tool(Tool::Select);
    let corner = at(&map, 20.0, 10.0);
    ed.drag_start(&mut m, &map, corner);
    ed.drag_to(&m, &map, at(&map, 30.0, 20.0));
    ed.drag_end(&mut m, &map, at(&map, 30.0, 20.0));
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let s = solved(&m, f);
    assert!(off_by(&s, o[0].1, &o[0].2) < 1e-12);
    let areas = region_areas(&m, f);
    // 50 × 30 outside, 42 × 22 inside.
    assert!(
        areas.iter().any(|a| (a - 924.0 * MM * MM).abs() < 1e-9),
        "{areas:?}"
    );
}

#[test]
fn a_slot_offsets_through_its_smooth_joints() {
    for (x, y, r) in [(0.0, 8.0, 8.0), (0.0, 2.0, 2.0)] {
        let (mut m, mut ed, map) = setup();
        let f = ed.feature;
        draw(
            &mut m,
            &mut ed,
            &map,
            Tool::Slot,
            &[(-20.0, 0.0), (20.0, 0.0), (0.0, 5.0)],
        );
        let before = dof(&m, f);
        ed.set_tool(Tool::Offset);
        click(&mut m, &mut ed, &map, 0.0, 5.0);
        assert_eq!(ed.message.as_deref(), Some("Chain of 4 curves chosen"));
        click(&mut m, &mut ed, &map, x, y);
        ed.prompt = None;
        let o = offsets(&m, f);
        assert_eq!(o.len(), 1);
        assert_eq!(o[0].3.len(), 4, "four smooth joints");
        assert_eq!(status(&m, f), SolveStatus::Solved);
        assert_eq!(redundant(&m, f), 0);
        assert_eq!(dof(&m, f), before);
        assert!(off_by(&solved(&m, f), o[0].1, &o[0].2) < 1e-12);
        // The offset slot: 40 mm between centres, radius r, around or
        // inside the 5 mm slot; the regions are the inner slot and the ring
        // between the two.
        let slot = |r: f64| (40.0 * 2.0 * r + std::f64::consts::PI * r * r) * MM * MM;
        let (inner, outer) = if r > 5.0 {
            (slot(5.0), slot(r))
        } else {
            (slot(r), slot(5.0))
        };
        let mut areas = region_areas(&m, f);
        areas.sort_by(f64::total_cmp);
        let mut want = [inner, outer - inner];
        want.sort_by(f64::total_cmp);
        assert_eq!(areas.len(), 2);
        for (a, w) in areas.iter().zip(want) {
            assert!((a - w).abs() < 1e-12, "{areas:?} {want:?}");
        }
    }
}

#[test]
fn an_open_chain_with_a_tangent_arc_and_a_corner() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    // (−30, 0) → (0, 0), a quarter arc about (0, 10) to (10, 10), then up to
    // (10, 30): smooth at (0, 0), a corner at (10, 10)… which is smooth
    // too for a vertical line: make the last line lean instead.
    let mut b = Batch::default();
    let l0 = b.line(DVec2::new(-30.0, 0.0) * MM, DVec2::ZERO, false);
    let arc = b.arc(
        DVec2::new(0.0, 10.0) * MM,
        DVec2::ZERO,
        DVec2::new(10.0, 10.0) * MM,
        false,
    );
    let l1 = b.line(
        DVec2::new(10.0, 10.0) * MM,
        DVec2::new(20.0, 30.0) * MM,
        false,
    );
    b.coincident(l0.end, arc.start);
    b.coincident(arc.end, l1.start);
    b.tangent_at(l0.curve, arc.curve, arc.start);
    m.edit_sketch(f, b.edits).unwrap();
    let before = dof(&m, f);
    ed.set_tool(Tool::Offset);
    click(&mut m, &mut ed, &map, -15.0, 0.0);
    assert_eq!(ed.message.as_deref(), Some("Chain of 3 curves chosen"));
    // 2 mm below the first line: outside the bend.
    click(&mut m, &mut ed, &map, -15.0, -2.0);
    ed.prompt = None;
    let o = offsets(&m, f);
    assert_eq!(o[0].2.len(), 3);
    assert_eq!(o[0].3.len(), 1, "one smooth joint");
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), before);
    let s = solved(&m, f);
    assert!(off_by(&s, o[0].1, &o[0].2) < 1e-12);
    // The free ends lie on their sources' end normals.
    let img0 = o[0].2[0].1;
    let Geometry::Line { start, .. } = s.entity(img0).unwrap().geometry else {
        panic!()
    };
    assert!((s.point(start).unwrap() - DVec2::new(-30.0, -2.0) * MM).length() < 1e-12);
    assert!(!images_cross(&m, f, &o[0].2));
}

#[test]
fn a_circle_offsets_inwards_and_outwards() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Circle,
        &[(0.0, 0.0), (10.0, 0.0)],
    );
    ed.set_tool(Tool::Offset);
    click(&mut m, &mut ed, &map, 10.0, 0.0);
    assert_eq!(ed.message.as_deref(), Some("Curve chosen"));
    click(&mut m, &mut ed, &map, 6.0, 0.0);
    assert_eq!(ed.message.as_deref(), Some("Offset 1 curve"));
    ed.prompt = None;
    click(&mut m, &mut ed, &map, 0.0, 10.0);
    click(&mut m, &mut ed, &map, 0.0, 12.0);
    ed.prompt = None;
    let radii: Vec<f64> = curves(&m, f, true)
        .iter()
        .filter_map(|(_, c)| match c {
            ResolvedCurve::Circle { radius, .. } => Some(radius / MM),
            _ => None,
        })
        .collect();
    let mut radii = radii;
    radii.sort_by(f64::total_cmp);
    assert_eq!(radii.len(), 3);
    for (r, want) in radii.iter().zip([6.0, 10.0, 12.0]) {
        assert!((r - want).abs() < 1e-9, "{radii:?}");
    }
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // The centre snapped to the origin: only the first circle's radius is
    // free, and the offsets add nothing.
    assert_eq!(dof(&m, f), 1);
}

#[test]
fn offsets_that_cannot_be_made_are_refused() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Rectangle,
        &[(-20.0, -10.0), (20.0, 10.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Spline,
        &[(30.0, 0.0), (40.0, 10.0), (50.0, 0.0)],
    );
    ed.set_tool(Tool::Offset);
    // A spline is refused.
    hover(&m, &mut ed, &map, 30.0, 0.0);
    let e = ed.curve_edit.clone();
    if let Some(e) = e {
        if let Err(why) = e.plan {
            assert!(why.contains("spline"), "{why}");
        }
    }
    // Inwards by more than half the height: the short sides turn back.
    click(&mut m, &mut ed, &map, 0.0, -10.0);
    hover(&m, &mut ed, &map, 0.0, 2.0);
    let why = ed.curve_edit.clone().unwrap().plan.unwrap_err();
    assert!(why.contains("too large"), "{why}");
    let before = m.sketch_def(f).unwrap().sketch.clone();
    click(&mut m, &mut ed, &map, 0.0, 2.0);
    assert!(ed.message.as_deref().unwrap().starts_with("Not done"));
    assert_eq!(m.sketch_def(f).unwrap().sketch, before);
    // Escape chooses another curve, then leaves the tool.
    ed.escape();
    assert!(ed.offset_chain.is_none());
    assert_eq!(ed.tool, Tool::Offset);
    ed.escape();
    assert_eq!(ed.tool, Tool::Select);
}

#[test]
fn deleting_an_image_keeps_the_rest_offset() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Rectangle,
        &[(-20.0, -10.0), (20.0, 10.0)],
    );
    ed.set_tool(Tool::Offset);
    click(&mut m, &mut ed, &map, 0.0, -10.0);
    click(&mut m, &mut ed, &map, 0.0, -15.0);
    ed.prompt = None;
    let o = offsets(&m, f);
    let image = o[0].2[1].1;
    ed.set_tool(Tool::Select);
    ed.selection = vec![image];
    ed.delete_selection(&mut m);
    let o = offsets(&m, f);
    assert_eq!(o[0].2.len(), 3);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert!(off_by(&solved(&m, f), o[0].1, &o[0].2) < 1e-12);
}

// ---------------------------------------------------------------------------
// Fuzz
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
enum Shape {
    Rectangle(DVec2, DVec2),
    Polygon(DVec2, DVec2, u32),
    Slot(DVec2, DVec2, DVec2),
    Circle(DVec2, f64),
    Chain(Vec<DVec2>),
}

fn shape() -> impl Strategy<Value = Shape> {
    let p = || (-50.0..50.0f64, -50.0..50.0f64).prop_map(|(x, y)| DVec2::new(x, y));
    prop_oneof![
        (p(), p()).prop_map(|(c, k)| Shape::Rectangle(c, k)),
        (p(), p(), 3u32..8).prop_map(|(c, v, n)| Shape::Polygon(c, v, n)),
        (p(), p(), p()).prop_map(|(a, b, s)| Shape::Slot(a, b, s)),
        (p(), 1.0..30.0f64).prop_map(|(c, r)| Shape::Circle(c, r)),
        proptest::collection::vec(p(), 2..6).prop_map(Shape::Chain),
    ]
}

fn add_shape(m: &mut Model, f: FeatureId, s: &Shape) {
    let free = |p: DVec2| Snap {
        at: p * MM,
        target: SnapTarget::Free,
    };
    let mut b = Batch::default();
    let made = match s {
        Shape::Rectangle(c, k) => shapes::centre_rectangle(&mut b, &free(*c), &free(*k), false),
        Shape::Polygon(c, v, n) => shapes::polygon(&mut b, &free(*c), &free(*v), *n, false),
        Shape::Slot(a, e, side) => shapes::slot(&mut b, &free(*a), &free(*e), *side * MM, false),
        Shape::Circle(c, r) => {
            let center = b.point(*c * MM);
            b.entity(
                Geometry::Circle {
                    center,
                    radius: r * MM,
                },
                false,
            );
            Ok(())
        }
        Shape::Chain(pts) => {
            let mut prev: Option<SketchId> = None;
            for w in pts.windows(2) {
                let l = b.line(w[0] * MM, w[1] * MM, false);
                if let Some(p) = prev {
                    b.coincident(p, l.start);
                }
                prev = Some(l.end);
            }
            Ok(())
        }
    };
    if made.is_ok() {
        m.edit_sketch(f, b.edits).unwrap();
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// Random shapes and chains offset by random distances to random sides
    /// through the editor's planner: every offset either commits, leaving
    /// the sketch solved without redundancy or added freedom, the sources
    /// where they were and every image at the distance with no two images
    /// crossing (by the kernel's intersections), or is refused with a
    /// message naming the offset; nothing else changes the sketch.
    #[test]
    fn random_offsets_are_exact_or_refused(
        shapes_in in proptest::collection::vec(shape(), 1..4),
        pick in 0usize..64,
        frac in 0.005..0.4f64,
        left in any::<bool>(),
    ) {
        let mut m = Model::new();
        let f = m.add_sketch(PlaneRef::Top);
        for s in &shapes_in {
            add_shape(&mut m, f, s);
        }
        prop_assume!(status(&m, f) == SolveStatus::Solved);
        prop_assume!(redundant(&m, f) == 0);
        let res = m.sketch_result(f).unwrap();
        let shown: Vec<(SketchId, ResolvedCurve)> =
            res.curves.iter().map(|(id, _, c)| (*id, c.clone())).collect();
        let sketch = res.solution.sketch.clone();
        let candidates: Vec<SketchId> = shown
            .iter()
            .filter(|(_, c)| !matches!(c, ResolvedCurve::Spline { .. }))
            .map(|(id, _)| *id)
            .collect();
        prop_assume!(!candidates.is_empty());
        let target = candidates[pick % candidates.len()];
        let chain = offset::chain(&sketch, target).unwrap();
        let d = frac * 30.0 * MM;
        let before = (dof(&m, f), m.sketch_def(f).unwrap().sketch.clone());
        let planned = offset::offset(&sketch, &shown, &chain, d, left);
        let Ok(o) = planned else {
            // The planner refuses only offsets it cannot construct.
            let why = planned.unwrap_err();
            prop_assert!(why.contains("do not meet") || why.contains("larger than"), "{}", why);
            return Ok(());
        };
        match m.edit_sketch_checked(f, o.plan.edits.clone()) {
            Ok(_) => {
                prop_assert_eq!(status(&m, f), SolveStatus::Solved);
                prop_assert_eq!(redundant(&m, f), 0);
                prop_assert_eq!(dof(&m, f), before.0);
                let after = solved(&m, f);
                for e in &sketch.entities {
                    let Geometry::Point { at } = e.geometry else { continue };
                    prop_assert!((after.point(e.id).unwrap() - at).length() < 1e-9);
                }
                let all = offsets(&m, f);
                let (_, dist, pieces, _) = all.last().unwrap();
                prop_assert!(off_by(&after, *dist, pieces) < 1e-9);
                prop_assert!(!images_cross(&m, f, pieces), "the images cross");
            }
            Err(why) => {
                prop_assert!(why.contains("offset"), "{}", why);
                prop_assert_eq!(&m.sketch_def(f).unwrap().sketch, &before.1);
            }
        }
    }
}
