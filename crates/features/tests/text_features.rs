//! Sketch text through regeneration (ADR 0015): a word's letters are
//! regions that extrude to their exact area, and an extrude keeps its
//! letters when the text is resized, turned or extended.

use std::collections::BTreeSet;

use opendatum_core_model::sketch::{Constraint, ConstraintEntry, Entity, Geometry, TextPiece};
use opendatum_core_model::{
    DVec2, Feature, FeatureId, InsertAt, ModelOp, PlaneRef, ProfileSource, RegionKey, ReplicaId,
    Session, SketchData, SketchDef, SketchEdit, SketchId,
};
use opendatum_features::{regenerate, FeatureStatus, StudioResult};
use opendatum_kernel::{Body, BrepKernel, Kernel, SketchArrangement};
use opendatum_sketch_solver::SolveStatus;
use opendatum_text::{lay_out, Font, Layout};

const MM: f64 = 1e-3;

fn insert(s: &mut Session, name: &str, f: Feature) -> FeatureId {
    let element = s.first_studio().unwrap();
    FeatureId(s.apply(ModelOp::InsertFeature {
        element,
        at: InsertAt::End,
        name: name.into(),
        feature: f,
    }))
}

fn edit(s: &mut Session, sketch: FeatureId, edits: Vec<SketchEdit>) {
    let element = s.first_studio().unwrap();
    let op = s.apply(ModelOp::EditSketch {
        element,
        feature: sketch,
        edits,
    });
    assert!(
        s.document().conflicts.iter().all(|c| c.op != op),
        "{:?}",
        s.document().conflicts
    );
}

fn regen(s: &Session) -> StudioResult<Body> {
    let el = s.document().element(s.first_studio().unwrap()).unwrap();
    regenerate(&BrepKernel::default(), &el.studio)
}

fn outline(layout: &Layout) -> Vec<Vec<TextPiece>> {
    layout
        .contours
        .iter()
        .map(|c| {
            c.iter()
                .map(|p| TextPiece {
                    degree: p.degree,
                    points: p.points.clone(),
                })
                .collect()
        })
        .collect()
}

/// Ids of a placed text.
struct Placed {
    text: SketchId,
    height: SketchId,
    angle: SketchId,
}

/// A text in `font` with its bottom-left corner fixed at the origin, its
/// box `h` mm high and its bottom edge at `angle` to a fixed horizontal
/// reference line.
fn place(s: &mut Session, sketch: FeatureId, content: &str, font: Font, h: f64) -> Placed {
    let layout = lay_out(content, font).unwrap();
    let mut edits = Vec::new();
    let mut next = 0u32;
    let mut id = || {
        next += 1;
        SketchId::new_in_batch(next - 1)
    };
    let add = |edits: &mut Vec<SketchEdit>, id: SketchId, g: Geometry<SketchId>, c: bool| {
        edits.push(SketchEdit::AddEntity(Entity {
            id,
            geometry: g,
            construction: c,
        }));
    };
    let w = h * layout.aspect;
    let at = [(0.0, 0.0), (w, 0.0), (w, h), (0.0, h)];
    let corners = at.map(|_| id());
    for (k, (x, y)) in corners.iter().zip(at) {
        add(
            &mut edits,
            *k,
            Geometry::Point {
                at: DVec2::new(x, y) * MM,
            },
            false,
        );
    }
    let text = id();
    add(
        &mut edits,
        text,
        Geometry::Text {
            corners,
            content: content.into(),
            font: font.name().into(),
            aspect: layout.aspect,
            mirrored: false,
            outline: outline(&layout),
        },
        false,
    );
    let lines = [0, 1, 2, 3].map(|i| {
        let l = id();
        add(
            &mut edits,
            l,
            Geometry::Line {
                start: corners[i],
                end: corners[(i + 1) % 4],
            },
            true,
        );
        l
    });
    // A fixed horizontal reference line below the text.
    let (r0, r1, reference) = (id(), id(), id());
    add(
        &mut edits,
        r0,
        Geometry::Point {
            at: DVec2::new(0.0, -5.0) * MM,
        },
        false,
    );
    add(
        &mut edits,
        r1,
        Geometry::Point {
            at: DVec2::new(10.0, -5.0) * MM,
        },
        false,
    );
    add(
        &mut edits,
        reference,
        Geometry::Line { start: r0, end: r1 },
        true,
    );
    let mut c = |edits: &mut Vec<SketchEdit>, constraint: Constraint<SketchId>| {
        let k = id();
        edits.push(SketchEdit::AddConstraint(ConstraintEntry {
            id: k,
            constraint,
        }));
        k
    };
    for p in [r0, r1] {
        let at = match p {
            p if p == r0 => DVec2::new(0.0, -5.0) * MM,
            _ => DVec2::new(10.0, -5.0) * MM,
        };
        c(&mut edits, Constraint::Fix { point: p, at });
    }
    c(
        &mut edits,
        Constraint::Fix {
            point: corners[0],
            at: DVec2::ZERO,
        },
    );
    let height = c(
        &mut edits,
        Constraint::Length {
            line: lines[3],
            value: h * MM,
        },
    );
    let angle = c(
        &mut edits,
        Constraint::Angle {
            a: reference,
            b: lines[0],
            value: 0.0,
        },
    );
    let element = s.first_studio().unwrap();
    let op = s.apply(ModelOp::EditSketch {
        element,
        feature: sketch,
        edits,
    });
    assert!(s.document().conflicts.iter().all(|c| c.op != op));
    let real = |k: SketchId| SketchId { op, index: k.index };
    Placed {
        text: real(text),
        height: real(height),
        angle: real(angle),
    }
}

/// The regions holding ink: every region but the counters (a counter's
/// outer boundary is another region's hole).
fn ink(a: &SketchArrangement) -> Vec<usize> {
    let set = |l: &opendatum_kernel::RegionLoop| -> BTreeSet<usize> {
        l.edges.iter().map(|(e, _)| *e).collect()
    };
    let holes: Vec<BTreeSet<usize>> = a
        .regions
        .iter()
        .flat_map(|r| r.holes.iter().map(set))
        .collect();
    (0..a.regions.len())
        .filter(|&i| !holes.contains(&set(&a.regions[i].outer)))
        .collect()
}

/// The ink area of a layout, by Green's theorem over its pieces (box
/// units², 1 high), independent of the kernel.
fn ink_area(layout: &Layout) -> f64 {
    fn de_boor(degree: usize, ctrl: &[DVec2], u: f64) -> DVec2 {
        let n = ctrl.len() - 1;
        let spans = n + 1 - degree;
        let mut knots = vec![0.0; degree + 1];
        for i in 1..spans {
            knots.push(i as f64 / spans as f64);
        }
        knots.extend(std::iter::repeat_n(1.0, degree + 1));
        let mut k = degree;
        while k < n && u >= knots[k + 1] {
            k += 1;
        }
        let mut d: Vec<DVec2> = (0..=degree).map(|j| ctrl[j + k - degree]).collect();
        for r in 1..=degree {
            for j in (r..=degree).rev() {
                let i = j + k - degree;
                let den = knots[i + degree + 1 - r] - knots[i];
                let a = if den == 0.0 {
                    0.0
                } else {
                    (u - knots[i]) / den
                };
                d[j] = d[j - 1] * (1.0 - a) + d[j] * a;
            }
        }
        d[degree]
    }
    let mut twice = 0.0;
    for p in layout.contours.iter().flatten() {
        let n = 4000;
        let mut prev = p.points[0];
        for k in 1..=n {
            let q = de_boor(p.degree as usize, &p.points, k as f64 / n as f64);
            twice += prev.perp_dot(q);
            prev = q;
        }
    }
    // TrueType outer contours run clockwise: ink has negative signed area.
    -twice * 0.5
}

fn keys(r: &StudioResult<Body>, sketch: FeatureId) -> (Vec<RegionKey>, Vec<usize>) {
    let regions = r.sketch(sketch).unwrap().regions.as_ref().unwrap();
    (regions.keys.clone(), ink(&regions.arrangement))
}

fn volume(r: &StudioResult<Body>, f: FeatureId) -> f64 {
    let b = r
        .bodies
        .iter()
        .find(|b| b.feature == f)
        .unwrap_or_else(|| panic!("no body; status {:?}", r.status(f)));
    BrepKernel::default()
        .mass_properties(&b.body, 1e-9)
        .unwrap()
        .volume
}

fn setup(content: &str, font: Font, h: f64) -> (Session, FeatureId, Placed, FeatureId, f64) {
    let mut s = Session::new_document(ReplicaId(1));
    let sk = insert(
        &mut s,
        "Sketch",
        Feature::Sketch(SketchDef::new(PlaneRef::Top, SketchData::default())),
    );
    let placed = place(&mut s, sk, content, font, h);
    let r = regen(&s);
    assert_eq!(r.status(sk), Some(&FeatureStatus::Ok), "{:?}", r.status(sk));
    let sr = r.sketch(sk).unwrap();
    assert_eq!(sr.solution.status, SolveStatus::Solved);
    assert_eq!(sr.solution.dof, 0);
    let (all, ink) = keys(&r, sk);
    let picked: Vec<RegionKey> = ink.iter().map(|&i| all[i].clone()).collect();
    let e = insert(
        &mut s,
        "Extrude",
        Feature::Extrude {
            profile: ProfileSource::Sketch {
                sketch: sk,
                regions: picked,
            },
            distance: 2.0 * MM,
            symmetric: false,
        },
    );
    let area = ink_area(&lay_out(content, font).unwrap()) * (h * MM).powi(2);
    (s, sk, placed, e, area)
}

fn assert_close(actual: f64, expected: f64, rel: f64) {
    assert!(
        (actual - expected).abs() <= rel * expected.abs(),
        "{actual} vs {expected}"
    );
}

#[test]
fn a_word_extrudes_to_its_letters() {
    for font in [Font::Sans, Font::SerifBold, Font::Mono] {
        let (s, sk, _, e, area) = setup("OpenDatum", font, 10.0);
        let r = regen(&s);
        assert_eq!(
            r.status(e),
            Some(&FeatureStatus::Ok),
            "{font:?}: {:?}",
            r.status(e)
        );
        // 9 letters of ink and 5 counters (O, p, e, D, a).
        let (all, ink) = keys(&r, sk);
        assert_eq!((all.len(), ink.len()), (14, 9), "{font:?}");
        assert_close(volume(&r, e), area * 2.0 * MM, 1e-6);
        // One side face per outline piece, two caps per letter.
        let layout = lay_out("OpenDatum", font).unwrap();
        let pieces: usize = layout.contours.iter().map(Vec::len).sum();
        let body = &r.bodies.iter().find(|b| b.feature == e).unwrap().body;
        let t = BrepKernel::default().topology(body);
        assert_eq!(
            (t.faces, t.edges, t.vertices),
            (pieces + 18, 3 * pieces, 2 * pieces),
            "{font:?}"
        );
    }
}

#[test]
fn the_letters_follow_a_resized_turned_and_extended_text() {
    let (mut s, sk, placed, e, area) = setup("Gear 42", Font::Sans, 8.0);
    let before = keys(&regen(&s), sk).0;
    // Taller: the area grows with the square of the height.
    edit(
        &mut s,
        sk,
        vec![SketchEdit::SetValue {
            constraint: placed.height,
            value: 12.0 * MM,
        }],
    );
    let r = regen(&s);
    assert_eq!(r.status(e), Some(&FeatureStatus::Ok), "{:?}", r.status(e));
    assert_close(volume(&r, e), area * 2.25 * 2.0 * MM, 1e-6);
    // Every key still matches exactly: pieces are named by their place in
    // the outline, not by where they are.
    assert_eq!(
        keys(&r, sk).0.iter().map(|k| &k.loops).collect::<Vec<_>>(),
        before.iter().map(|k| &k.loops).collect::<Vec<_>>()
    );
    // Turned by 40°: same letters, same area.
    edit(
        &mut s,
        sk,
        vec![SketchEdit::SetValue {
            constraint: placed.angle,
            value: 40f64.to_radians(),
        }],
    );
    let r = regen(&s);
    assert_eq!(r.status(e), Some(&FeatureStatus::Ok), "{:?}", r.status(e));
    assert_close(volume(&r, e), area * 2.25 * 2.0 * MM, 1e-6);
    // Extended with more characters: the extrude keeps exactly the letters
    // it had, and the new ones are not added to it.
    let longer = lay_out("Gear 42 teeth", Font::Sans).unwrap();
    edit(
        &mut s,
        sk,
        vec![SketchEdit::SetText {
            entity: placed.text,
            content: "Gear 42 teeth".into(),
            font: Font::Sans.name().into(),
            aspect: longer.aspect,
            outline: outline(&longer),
        }],
    );
    let r = regen(&s);
    assert_eq!(r.status(e), Some(&FeatureStatus::Ok), "{:?}", r.status(e));
    assert_close(volume(&r, e), area * 2.25 * 2.0 * MM, 1e-6);
    assert!(keys(&r, sk).0.len() > before.len());
}

#[test]
fn malformed_text_is_refused_by_the_document() {
    let mut s = Session::new_document(ReplicaId(1));
    let sk = insert(
        &mut s,
        "Sketch",
        Feature::Sketch(SketchDef::new(PlaneRef::Top, SketchData::default())),
    );
    let placed = place(&mut s, sk, "ok", Font::Sans, 5.0);
    let element = s.first_studio().unwrap();
    let op = s.apply(ModelOp::EditSketch {
        element,
        feature: sk,
        edits: vec![SketchEdit::SetText {
            entity: placed.text,
            content: "bad".into(),
            font: "Liberation Sans".into(),
            aspect: -1.0,
            outline: Vec::new(),
        }],
    });
    let conflict = s
        .document()
        .conflicts
        .iter()
        .find(|c| c.op == op)
        .expect("refused");
    assert!(
        format!("{conflict:?}").contains("aspect ratio"),
        "{conflict:?}"
    );
}
