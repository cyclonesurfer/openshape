//! References to the faces of earlier parts (ADR 0019, design note 0013
//! §10 step 3): sketches on planar faces get the frame the ADR gives them
//! and follow their faces through edits; extrudes end on parallel faces;
//! references to what is gone, later, curved or not parallel fail with a
//! message; documents saved before references existed regenerate as they
//! did.

use std::f64::consts::PI;

use opendatum_core_model::naming::{FaceName, Hint, Name, Reference, Role, Source, Use};
use opendatum_core_model::sketch::{Entity, Geometry};
use opendatum_core_model::{
    CurveDef, DVec2, DVec3, Feature, FeatureId, InsertAt, ModelOp, OpId, PlaneRef, ProfileDef,
    ProfileSource, ReplicaId, Session, SketchData, SketchDef, SketchId,
};
use opendatum_features::{face_frame, regenerate, FeatureStatus, PartBody, StudioResult};
use opendatum_kernel::{Body, BrepKernel, FacePlane, Frame, Kernel, TessellationOptions};

const MM: f64 = 1e-3;

fn k() -> BrepKernel {
    BrepKernel::default()
}

fn insert(s: &mut Session, name: &str, f: Feature) -> FeatureId {
    let element = s.first_studio().unwrap();
    FeatureId(s.apply(ModelOp::InsertFeature {
        element,
        at: InsertAt::End,
        name: name.into(),
        feature: f,
    }))
}

fn edit(s: &mut Session, feature: FeatureId, definition: Feature) {
    let element = s.first_studio().unwrap();
    let op = s.apply(ModelOp::EditFeature {
        element,
        feature,
        definition,
    });
    assert!(s.document().conflicts.iter().all(|c| c.op != op));
}

fn regen(s: &Session) -> StudioResult<Body> {
    let el = s.document().element(s.first_studio().unwrap()).unwrap();
    regenerate(&k(), &el.studio)
}

fn ok(r: &StudioResult<Body>, f: FeatureId) {
    assert_eq!(r.status(f), Some(&FeatureStatus::Ok), "{:?}", r.status(f));
}

fn failed(r: &StudioResult<Body>, f: FeatureId, why: &str) {
    match r.status(f) {
        Some(FeatureStatus::Failed(m)) => assert!(m.contains(why), "{m:?} lacks {why:?}"),
        other => panic!("expected a failure with {why:?}, got {other:?}"),
    }
}

fn part(r: &StudioResult<Body>, f: FeatureId) -> &PartBody<Body> {
    r.bodies
        .iter()
        .find(|b| b.feature == f)
        .unwrap_or_else(|| panic!("no part; status {:?}", r.status(f)))
}

fn volume(r: &StudioResult<Body>, f: FeatureId) -> f64 {
    k().mass_properties(&part(r, f).body, 1e-10).unwrap().volume
}

fn bounds(r: &StudioResult<Body>, f: FeatureId) -> (DVec3, DVec3) {
    k().tessellate(&part(r, f).body, &TessellationOptions::default())
        .unwrap()
        .bounds()
        .unwrap()
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * (1.0 + b.abs())
}

fn assert_close(a: f64, b: f64) {
    assert!(close(a, b), "{a} vs {b}");
}

fn assert_vec(a: DVec3, b: DVec3) {
    assert!((a - b).length() <= 1e-9 * (1.0 + b.length()), "{a} vs {b}");
}

/// The extrude's mesh lies between `lo` and `hi` (within the default
/// tessellation's chord tolerance, which is exact for straight sides).
fn assert_bounds(r: &StudioResult<Body>, f: FeatureId, lo: DVec3, hi: DVec3) {
    let (a, b) = bounds(r, f);
    assert!(
        (a - lo).abs().max_element() < 1e-9 && (b - hi).abs().max_element() < 1e-9,
        "bounds {a}..{b}, expected {lo}..{hi}"
    );
}

/// A carrier reference to the one face of `f`'s part whose name passes
/// `pick`, with a hint inside it.
fn face_ref(r: &StudioResult<Body>, f: FeatureId, pick: impl Fn(&FaceName) -> bool) -> Reference {
    let p = part(r, f);
    let names = p.names.as_ref().expect("named");
    let found: Vec<_> = names.faces.iter().filter(|l| pick(&l.name)).collect();
    assert_eq!(found.len(), 1, "faces matching: {}", found.len());
    let mesh = k()
        .tessellate(&p.body, &TessellationOptions::default())
        .unwrap();
    let fm = &mesh.faces[found[0].pieces[0].0 as usize];
    let t = fm.triangles[0].map(|i| fm.positions[i as usize]);
    let point = (t[0] + t[1] + t[2]) / 3.0;
    let direction = (t[1] - t[0]).cross(t[2] - t[0]).normalize();
    Reference {
        name: Name::Face(found[0].name.clone()),
        hint: Hint { point, direction },
        usage: Use::Carrier,
    }
}

fn box_side(axis: u8, positive: bool) -> impl Fn(&FaceName) -> bool {
    move |n| n.role == Role::BoxSide { axis, positive }
}

fn cap(end: bool) -> impl Fn(&FaceName) -> bool {
    move |n| n.role == Role::Cap { end }
}

/// The side face swept from the sketch entity `entity`.
fn side_of(entity: SketchId) -> impl Fn(&FaceName) -> bool {
    move |n| {
        n.role == Role::Side
            && n.sources
                .iter()
                .any(|s| matches!(s, Source::SketchEdge(e) if e.entity == entity))
    }
}

fn add_box(s: &mut Session, corner: DVec3, size: DVec3) -> FeatureId {
    insert(s, "Box", Feature::Box { corner, size })
}

/// Sketch entity ids: the sketches here are built whole, with ids made up
/// in an op of their own, so a redefinition keeps them.
fn sid(sketch: u64, index: u32) -> SketchId {
    SketchId {
        op: OpId {
            lamport: 1000 + sketch,
            replica: ReplicaId(77),
        },
        index,
    }
}

/// A sketch of an axis-aligned rectangle (`x0, y0` to `x1, y1`, mm, in
/// plane coordinates) whose sides are entities 4 (bottom), 5 (right),
/// 6 (top) and 7 (left), plus circles `(cx, cy, r)` (mm).
fn rect_sketch(n: u64, plane: PlaneRef, r: [f64; 4], circles: &[(f64, f64, f64)]) -> Feature {
    let [x0, y0, x1, y1] = r;
    let mut entities = Vec::new();
    let corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)];
    for (i, (x, y)) in corners.iter().enumerate() {
        entities.push(Entity {
            id: sid(n, i as u32),
            geometry: Geometry::Point {
                at: DVec2::new(*x, *y) * MM,
            },
            construction: false,
        });
    }
    for i in 0..4u32 {
        entities.push(Entity {
            id: sid(n, 4 + i),
            geometry: Geometry::Line {
                start: sid(n, i),
                end: sid(n, (i + 1) % 4),
            },
            construction: false,
        });
    }
    for (j, &(cx, cy, rad)) in circles.iter().enumerate() {
        let c = sid(n, 8 + 2 * j as u32);
        entities.push(Entity {
            id: c,
            geometry: Geometry::Point {
                at: DVec2::new(cx, cy) * MM,
            },
            construction: false,
        });
        entities.push(Entity {
            id: sid(n, 9 + 2 * j as u32),
            geometry: Geometry::Circle {
                center: c,
                radius: rad * MM,
            },
            construction: false,
        });
    }
    let data = SketchData {
        entities,
        constraints: Vec::new(),
    };
    Feature::Sketch(SketchDef::new(plane, data))
}

/// A sketch of one circle, entity 1 about point 0.
fn circle_sketch(n: u64, plane: PlaneRef, cx: f64, cy: f64, r: f64) -> Feature {
    let data = SketchData {
        entities: vec![
            Entity {
                id: sid(n, 0),
                geometry: Geometry::Point {
                    at: DVec2::new(cx, cy) * MM,
                },
                construction: false,
            },
            Entity {
                id: sid(n, 1),
                geometry: Geometry::Circle {
                    center: sid(n, 0),
                    radius: r * MM,
                },
                construction: false,
            },
        ],
        constraints: Vec::new(),
    };
    Feature::Sketch(SketchDef::new(plane, data))
}

/// Extrudes every region of a sketch.
fn extrude_def(
    s: &Session,
    sketch: FeatureId,
    distance: f64,
    symmetric: bool,
    up_to: Option<Box<Reference>>,
) -> Feature {
    let r = regen(s);
    let regions = r
        .sketch(sketch)
        .unwrap()
        .regions
        .as_ref()
        .expect("regions")
        .keys
        .clone();
    Feature::Extrude {
        profile: ProfileSource::Sketch { sketch, regions },
        distance: distance * MM,
        symmetric,
        up_to,
    }
}

fn extrude(s: &mut Session, sketch: FeatureId, distance: f64) -> FeatureId {
    let f = extrude_def(s, sketch, distance, false, None);
    insert(s, "Extrude", f)
}

fn extrude_up_to(s: &mut Session, sketch: FeatureId, face: Reference) -> FeatureId {
    let f = extrude_def(s, sketch, 1.0, false, Some(Box::new(face)));
    insert(s, "Extrude", f)
}

/// An inline cylinder on the Top plane, as the application's template
/// makes one.
fn cylinder(s: &mut Session, cx: f64, cy: f64, r: f64, h: f64) -> FeatureId {
    insert(
        s,
        "Cylinder",
        Feature::Extrude {
            profile: ProfileSource::Inline(ProfileDef {
                plane: PlaneRef::Top,
                curves: vec![CurveDef::Circle {
                    center: DVec2::new(cx, cy) * MM,
                    radius: r * MM,
                }],
            }),
            distance: h * MM,
            symmetric: false,
            up_to: None,
        },
    )
}

#[test]
fn documents_saved_before_references_regenerate_unchanged() {
    // Saved by the code before ADR 0019: a box, a rectangle sketched on the
    // Front plane and extruded, and a symmetric inline cylinder.
    let s = Session::open(include_bytes!("fixtures/pre-0019.oplog"), ReplicaId(9)).unwrap();
    assert!(s.document().conflicts.is_empty());
    let studio = &s
        .document()
        .element(s.first_studio().unwrap())
        .unwrap()
        .studio;
    assert_eq!(studio.features.len(), 4);
    let r = regenerate(&k(), &studio.clone());
    let mut volumes = Vec::new();
    for f in &studio.features {
        ok(&r, f.id);
        match &f.feature {
            Feature::Extrude { up_to, .. } => assert!(up_to.is_none()),
            Feature::Sketch(def) => assert_eq!(def.plane, PlaneRef::Front),
            _ => {}
        }
        if !matches!(f.feature, Feature::Sketch(_)) {
            assert!(part(&r, f.id).names.is_ok());
            volumes.push(volume(&r, f.id));
        }
    }
    let want = [
        0.04 * 0.03 * 0.02,
        20.0 * 10.0 * 15.0 * MM.powi(3),
        PI * 0.01f64.powi(2) * 0.025,
    ];
    for (v, w) in volumes.iter().zip(want) {
        assert_close(*v, w);
    }
}

#[test]
fn face_frames_follow_the_rule_of_adr_0019() {
    let at = |n: DVec3| {
        face_frame(&FacePlane {
            point: n * 0.3 + DVec3::new(0.7, -0.2, 0.1).cross(n),
            normal: n,
        })
    };
    // Faces parallel to the standard planes get their frames, offset.
    let cases = [
        (DVec3::Z, DVec3::X, DVec3::Y),
        (-DVec3::Z, DVec3::X, -DVec3::Y),
        (-DVec3::Y, DVec3::X, DVec3::Z),
        (DVec3::Y, -DVec3::X, DVec3::Z),
        (DVec3::X, DVec3::Y, DVec3::Z),
        (-DVec3::X, -DVec3::Y, DVec3::Z),
    ];
    for (n, x, y) in cases {
        let f = at(n);
        assert_vec(f.z, n);
        assert_vec(f.x, x);
        assert_vec(f.y, y);
        assert_vec(f.origin, n * 0.3);
    }
    assert_eq!(
        at(DVec3::Z),
        Frame {
            origin: DVec3::Z * 0.3,
            ..Frame::top()
        }
    );
    // Tilted faces: right-handed and orthonormal, origin the projection of
    // the world origin, upright when facing sideways, x along the world X
    // projection when facing up or down.
    for i in 0..200 {
        let a = i as f64 * 0.731;
        let b = (i as f64 * 0.377).sin() * PI;
        let n = DVec3::new(b.cos() * a.cos(), b.cos() * a.sin(), b.sin());
        let f = at(n);
        assert_vec(f.z, n);
        assert_vec(f.x.cross(f.y), f.z);
        assert!(f.x.dot(f.y).abs() < 1e-12 && (f.x.length() - 1.0).abs() < 1e-12);
        assert!(
            f.origin.cross(n).length() < 1e-12,
            "the origin is on the normal"
        );
        assert_close(f.origin.dot(n), 0.3);
        if n.z.abs() >= n.x.abs().max(n.y.abs()) {
            let projected = (DVec3::X - n * n.x).normalize();
            assert_vec(f.x, projected);
        } else {
            assert!(f.x.z.abs() < 1e-12, "x is horizontal: {}", f.x);
            assert!(f.y.z > 0.0, "y points up: {}", f.y);
        }
    }
}

#[test]
fn sketches_on_box_faces_sit_on_them_seen_from_outside() {
    let corner = DVec3::new(0.01, 0.02, 0.03);
    let size = DVec3::new(0.04, 0.03, 0.02);
    for axis in 0..3u8 {
        for positive in [false, true] {
            let mut s = Session::new_document(ReplicaId(1));
            let b = add_box(&mut s, corner, size);
            let r = regen(&s);
            let face = face_ref(&r, b, box_side(axis, positive));
            let sk = insert(
                &mut s,
                "Sketch",
                rect_sketch(
                    1,
                    PlaneRef::Face(Box::new(face)),
                    [-5.0, -2.5, 5.0, 2.5],
                    &[],
                ),
            );
            let e = extrude(&mut s, sk, 4.0);
            let r = regen(&s);
            ok(&r, sk);
            ok(&r, e);
            let mut n = DVec3::ZERO;
            n[axis as usize] = if positive { 1.0 } else { -1.0 };
            let level = if positive {
                corner[axis as usize] + size[axis as usize]
            } else {
                corner[axis as usize]
            };
            let frame = r.sketch(sk).unwrap().frame;
            let want = at_level(n, level);
            for (a, b) in [
                (frame.origin, want.origin),
                (frame.x, want.x),
                (frame.y, want.y),
                (frame.z, want.z),
            ] {
                assert_vec(a, b);
            }
            // The extrude sits on the face, outside the box.
            assert_close(volume(&r, e), 10.0 * 5.0 * 4.0 * MM.powi(3));
            let (lo, hi) = bounds(&r, e);
            let (near, far) = if positive {
                (lo[axis as usize], hi[axis as usize])
            } else {
                (-hi[axis as usize], -lo[axis as usize])
            };
            let level = if positive { level } else { -level };
            assert_close(near, level);
            assert_close(far, level + 0.004);
        }
    }
}

/// The frame [`face_frame`] gives a face with outward normal `n` (an axis)
/// at signed level `level` along its axis.
fn at_level(n: DVec3, level: f64) -> Frame {
    let axis_dir = n.abs();
    face_frame(&FacePlane {
        point: axis_dir * level,
        normal: n,
    })
}

#[test]
fn a_sketch_follows_its_face_through_edits() {
    let mut s = Session::new_document(ReplicaId(1));
    let size = DVec3::new(0.04, 0.03, 0.02);
    let b = add_box(&mut s, DVec3::ZERO, size);
    let face = face_ref(&regen(&s), b, box_side(2, true));
    let sk = insert(
        &mut s,
        "Sketch",
        rect_sketch(
            1,
            PlaneRef::Face(Box::new(face)),
            [5.0, 5.0, 15.0, 10.0],
            &[],
        ),
    );
    let e = extrude(&mut s, sk, 5.0);
    let r = regen(&s);
    assert_bounds(
        &r,
        e,
        DVec3::new(0.005, 0.005, 0.02),
        DVec3::new(0.015, 0.01, 0.025),
    );
    edit(
        &mut s,
        b,
        Feature::Box {
            corner: DVec3::new(0.0, 0.0, -0.01),
            size: DVec3::new(0.04, 0.03, 0.045),
        },
    );
    let r = regen(&s);
    ok(&r, sk);
    ok(&r, e);
    assert_bounds(
        &r,
        e,
        DVec3::new(0.005, 0.005, 0.035),
        DVec3::new(0.015, 0.01, 0.04),
    );
    // A box built downwards still has its +z side on top.
    edit(
        &mut s,
        b,
        Feature::Box {
            corner: DVec3::new(0.0, 0.0, 0.05),
            size: DVec3::new(0.04, 0.03, -0.05),
        },
    );
    let r = regen(&s);
    ok(&r, e);
    assert_bounds(
        &r,
        e,
        DVec3::new(0.005, 0.005, 0.05),
        DVec3::new(0.015, 0.01, 0.055),
    );
}

#[test]
fn sketches_on_sketched_parts_follow_their_caps_and_sides() {
    let mut s = Session::new_document(ReplicaId(1));
    let base = insert(
        &mut s,
        "Base sketch",
        rect_sketch(1, PlaneRef::Top, [0.0, 0.0, 40.0, 20.0], &[]),
    );
    let block = extrude(&mut s, base, 10.0);
    let r = regen(&s);
    let top = face_ref(&r, block, cap(true));
    let right = face_ref(&r, block, side_of(sid(1, 5)));
    let boss_sketch = insert(
        &mut s,
        "Boss sketch",
        circle_sketch(2, PlaneRef::Face(Box::new(top)), 20.0, 10.0, 5.0),
    );
    let boss = extrude(&mut s, boss_sketch, 5.0);
    let tab_sketch = insert(
        &mut s,
        "Tab sketch",
        rect_sketch(
            3,
            PlaneRef::Face(Box::new(right)),
            [5.0, 2.0, 10.0, 6.0],
            &[],
        ),
    );
    let tab = extrude(&mut s, tab_sketch, 3.0);
    let r = regen(&s);
    for f in [boss_sketch, boss, tab_sketch, tab] {
        ok(&r, f);
    }
    assert_close(volume(&r, boss), PI * 25.0 * 5.0 * MM.powi(3));
    let (lo, hi) = bounds(&r, boss);
    assert_close(lo.z, 0.01);
    assert_close(hi.z, 0.015);
    // The right side faces +x: its sketch is upright, x along +y.
    assert_bounds(
        &r,
        tab,
        DVec3::new(0.04, 0.005, 0.002),
        DVec3::new(0.043, 0.01, 0.006),
    );
    // A taller, longer block: the boss rises with the top cap and the tab
    // moves out with the right side.
    let longer = rect_sketch(1, PlaneRef::Top, [0.0, 0.0, 55.0, 20.0], &[]);
    edit(&mut s, base, longer);
    let taller = extrude_def(&s, base, 15.0, false, None);
    edit(&mut s, block, taller);
    let r = regen(&s);
    for f in [block, boss_sketch, boss, tab_sketch, tab] {
        ok(&r, f);
    }
    let (lo, hi) = bounds(&r, boss);
    assert_close(lo.z, 0.015);
    assert_close(hi.z, 0.02);
    assert_bounds(
        &r,
        tab,
        DVec3::new(0.055, 0.005, 0.002),
        DVec3::new(0.058, 0.01, 0.006),
    );
    // A hole in the base region splits nothing the references use: both
    // still resolve, the cap through its region's stored outer loop.
    let holed = rect_sketch(
        1,
        PlaneRef::Top,
        [0.0, 0.0, 55.0, 20.0],
        &[(45.0, 10.0, 3.0)],
    );
    edit(&mut s, base, holed);
    let r = regen(&s);
    for f in [block, boss_sketch, boss, tab_sketch, tab] {
        ok(&r, f);
    }
}

#[test]
fn sketches_on_inline_profiles_use_their_caps() {
    let mut s = Session::new_document(ReplicaId(1));
    let c = cylinder(&mut s, 0.0, 0.0, 10.0, 25.0);
    let r = regen(&s);
    let top = face_ref(&r, c, cap(true));
    let bottom = face_ref(&r, c, cap(false));
    let pin = insert(
        &mut s,
        "Pin sketch",
        circle_sketch(1, PlaneRef::Face(Box::new(top)), 0.0, 0.0, 2.0),
    );
    let pin_e = extrude(&mut s, pin, 6.0);
    let foot = insert(
        &mut s,
        "Foot sketch",
        circle_sketch(2, PlaneRef::Face(Box::new(bottom)), 0.0, 0.0, 4.0),
    );
    let foot_e = extrude(&mut s, foot, 3.0);
    let side = face_ref(&r, c, |n| n.role == Role::Side);
    let on_side = insert(
        &mut s,
        "On the side",
        circle_sketch(3, PlaneRef::Face(Box::new(side)), 0.0, 0.0, 1.0),
    );
    let r = regen(&s);
    ok(&r, pin_e);
    ok(&r, foot_e);
    let (lo, hi) = bounds(&r, pin_e);
    assert_close(lo.z, 0.025);
    assert_close(hi.z, 0.031);
    // The bottom cap faces down: its sketch extrudes downwards.
    let (lo, hi) = bounds(&r, foot_e);
    assert_close(lo.z, -0.003);
    assert_close(hi.z, 0.0);
    failed(&r, on_side, "sketch face is not planar");
}

#[test]
fn extrudes_end_on_parallel_faces() {
    let mut s = Session::new_document(ReplicaId(1));
    let b = add_box(&mut s, DVec3::ZERO, DVec3::new(0.04, 0.03, 0.02));
    let low = add_box(
        &mut s,
        DVec3::new(0.0, 0.05, -0.05),
        DVec3::new(0.01, 0.01, 0.02),
    );
    let c = cylinder(&mut s, -30.0, 0.0, 5.0, 10.0);
    let sk = insert(
        &mut s,
        "Sketch",
        rect_sketch(1, PlaneRef::Top, [50.0, 0.0, 60.0, 10.0], &[]),
    );
    let r = regen(&s);
    let top = face_ref(&r, b, box_side(2, true));
    let up = extrude_up_to(&mut s, sk, top.clone());
    let down = extrude_up_to(&mut s, sk, face_ref(&r, low, box_side(2, true)));
    let r = regen(&s);
    ok(&r, up);
    ok(&r, down);
    assert_close(volume(&r, up), 10.0 * 10.0 * 20.0 * MM.powi(3));
    assert_bounds(
        &r,
        up,
        DVec3::new(0.05, 0.0, 0.0),
        DVec3::new(0.06, 0.01, 0.02),
    );
    assert_bounds(
        &r,
        down,
        DVec3::new(0.05, 0.0, -0.03),
        DVec3::new(0.06, 0.01, 0.0),
    );
    // The extrude follows its face.
    edit(
        &mut s,
        b,
        Feature::Box {
            corner: DVec3::ZERO,
            size: DVec3::new(0.04, 0.03, 0.035),
        },
    );
    let r = regen(&s);
    assert_bounds(
        &r,
        up,
        DVec3::new(0.05, 0.0, 0.0),
        DVec3::new(0.06, 0.01, 0.035),
    );

    // Faces it cannot end on.
    let cases = [
        (face_ref(&r, b, box_side(2, false)), "in the sketch plane"),
        (face_ref(&r, b, box_side(0, true)), "not parallel"),
        (face_ref(&r, c, |n| n.role == Role::Side), "not planar"),
    ];
    for (face, why) in cases {
        let f = extrude_up_to(&mut s, sk, face);
        failed(&regen(&s), f, why);
    }
    let sym = extrude_def(&s, sk, 5.0, true, Some(Box::new(top)));
    let sym = insert(&mut s, "Symmetric", sym);
    failed(&regen(&s), sym, "symmetric");
}

#[test]
fn references_to_parts_that_are_gone_or_later_fail() {
    let mut s = Session::new_document(ReplicaId(1));
    let b = add_box(&mut s, DVec3::ZERO, DVec3::new(0.04, 0.03, 0.02));
    let face = face_ref(&regen(&s), b, box_side(2, true));
    let sk = insert(
        &mut s,
        "Sketch",
        rect_sketch(
            1,
            PlaneRef::Face(Box::new(face.clone())),
            [5.0, 5.0, 15.0, 10.0],
            &[],
        ),
    );
    let e = extrude(&mut s, sk, 5.0);
    let flat = insert(
        &mut s,
        "Flat sketch",
        rect_sketch(2, PlaneRef::Top, [50.0, 0.0, 60.0, 10.0], &[]),
    );
    let up = extrude_up_to(&mut s, flat, face);
    let r = regen(&s);
    for f in [sk, e, up] {
        ok(&r, f);
    }
    let element = s.first_studio().unwrap();
    let gone = "is on a part that is missing, suppressed, failed or later in the list";

    // Suppressed.
    s.apply(ModelOp::SetSuppressed {
        element,
        feature: b,
        suppressed: true,
    });
    let r = regen(&s);
    failed(&r, sk, &format!("the sketch face {gone}"));
    failed(&r, e, "sketch regions");
    failed(&r, up, &format!("the face to extrude up to {gone}"));
    s.undo();
    ok(&regen(&s), e);

    // Moved after the features that use it.
    s.apply(ModelOp::MoveFeature {
        element,
        feature: b,
        to: InsertAt::End,
    });
    let r = regen(&s);
    failed(&r, sk, gone);
    failed(&r, up, gone);
    s.undo();

    // Deleted.
    s.apply(ModelOp::DeleteFeature {
        element,
        feature: b,
    });
    let r = regen(&s);
    failed(&r, sk, gone);
    failed(&r, up, gone);
    s.undo();
    let r = regen(&s);
    for f in [sk, e, up] {
        ok(&r, f);
    }
}

#[test]
fn references_to_faces_that_no_longer_exist_fail() {
    let mut s = Session::new_document(ReplicaId(1));
    let base = insert(
        &mut s,
        "Base sketch",
        rect_sketch(1, PlaneRef::Top, [0.0, 0.0, 40.0, 20.0], &[]),
    );
    let block = extrude(&mut s, base, 10.0);
    let r = regen(&s);
    let right = face_ref(&r, block, side_of(sid(1, 5)));
    let on_right = insert(
        &mut s,
        "On the right",
        circle_sketch(2, PlaneRef::Face(Box::new(right)), 10.0, 5.0, 2.0),
    );
    ok(&regen(&s), on_right);
    // Redraw the base without its right side: the face is gone, and the
    // reference fails rather than land on another side.
    let data = match rect_sketch(1, PlaneRef::Top, [0.0, 0.0, 40.0, 20.0], &[]) {
        Feature::Sketch(mut def) => {
            def.sketch.entities.retain(|e| e.id != sid(1, 5));
            def.sketch.entities.push(Entity {
                id: sid(1, 20),
                geometry: Geometry::Line {
                    start: sid(1, 1),
                    end: sid(1, 2),
                },
                construction: false,
            });
            Feature::Sketch(def)
        }
        _ => unreachable!(),
    };
    edit(&mut s, base, data);
    let r = regen(&s);
    ok(&r, block);
    failed(&r, on_right, "the sketch face no longer exists");
}

#[test]
fn a_face_in_pieces_is_one_plane() {
    // A line hanging from the middle of the base's bottom side splits that
    // side in two: the block's front face is two kernel faces joined
    // across a seam, one face to the user.
    let split = |on: bool| match rect_sketch(1, PlaneRef::Top, [0.0, 0.0, 40.0, 20.0], &[]) {
        Feature::Sketch(mut def) => {
            if on {
                for (i, y) in [(31, 0.0), (32, -10.0)] {
                    def.sketch.entities.push(Entity {
                        id: sid(1, i),
                        geometry: Geometry::Point {
                            at: DVec2::new(20.0, y) * MM,
                        },
                        construction: false,
                    });
                }
                def.sketch.entities.push(Entity {
                    id: sid(1, 30),
                    geometry: Geometry::Line {
                        start: sid(1, 31),
                        end: sid(1, 32),
                    },
                    construction: false,
                });
            }
            Feature::Sketch(def)
        }
        _ => unreachable!(),
    };
    let mut s = Session::new_document(ReplicaId(1));
    let base = insert(&mut s, "Base sketch", split(true));
    let block = extrude(&mut s, base, 10.0);
    let r = regen(&s);
    let front = face_ref(&r, block, side_of(sid(1, 4)));
    let names = part(&r, block).names.as_ref().unwrap();
    let Name::Face(front_name) = &front.name else {
        unreachable!()
    };
    let logical = names.faces.iter().find(|l| l.name == *front_name).unwrap();
    assert_eq!(logical.pieces.len(), 2, "the front face is in two pieces");
    let sk = insert(
        &mut s,
        "On the front",
        rect_sketch(
            2,
            PlaneRef::Face(Box::new(front)),
            [10.0, 2.0, 30.0, 8.0],
            &[],
        ),
    );
    let e = extrude(&mut s, sk, 4.0);
    let r = regen(&s);
    ok(&r, sk);
    ok(&r, e);
    // The front faces −y: its sketch has the Front plane's frame.
    assert_bounds(
        &r,
        e,
        DVec3::new(0.01, -0.004, 0.002),
        DVec3::new(0.03, 0.0, 0.008),
    );
    // Without the hanging line the front is one piece; the sketch stays.
    edit(&mut s, base, split(false));
    let r = regen(&s);
    ok(&r, sk);
    assert_bounds(
        &r,
        e,
        DVec3::new(0.01, -0.004, 0.002),
        DVec3::new(0.03, 0.0, 0.008),
    );
}
