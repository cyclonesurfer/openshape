//! Sketch regions: planar arrangements of curves (design note 0006, §6).

mod arrangement_util;

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use arrangement_util::*;
use openshape_kernel::arrangement::{arrange, region_containing, region_profile_with_tolerance};
use openshape_kernel::{DVec2, Kernel, KernelError, ProfileCurve, SketchArrangement, Tolerance};

#[test]
fn empty_input_has_no_regions() {
    let a = k().sketch_regions(&[]).unwrap();
    assert_eq!(a, SketchArrangement::default());
}

#[test]
fn rectangle_with_diagonals_has_four_regions() {
    let (w, h) = (4.0, 3.0);
    let mut c = polygon(&[v(0.0, 0.0), v(w, 0.0), v(w, h), v(0.0, h)]);
    c.push(line(v(0.0, 0.0), v(w, h)));
    c.push(line(v(w, 0.0), v(0.0, h)));
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 4);
    for r in &a.regions {
        assert_rel(r.area, w * h / 4.0, 1e-12);
        assert!(r.holes.is_empty());
        assert_eq!(r.outer.edges.len(), 3);
    }
    assert_eq!(a.vertices.len(), 5);
    assert_eq!(a.edges.len(), 8);
    let centre = a
        .vertices
        .iter()
        .find(|x| (x.point - v(w / 2.0, h / 2.0)).length() < 1e-12)
        .expect("diagonals cross at the centre");
    assert_eq!(centre.curves, vec![4, 5]);
    assert_rel(centre.params[0], 0.5, 1e-12);
    assert_rel(centre.params[1], 0.5, 1e-12);
    extrude_all(&c, &a);
}

#[test]
fn two_crossing_circles_have_three_regions() {
    let c = vec![circle(v(0.0, 0.0), 1.0), circle(v(1.0, 0.0), 1.0)];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 3);
    let lens = 2.0 * (0.5f64).acos() - 0.5 * 3.0f64.sqrt();
    let mut areas = sorted_areas(&a);
    areas.sort_by(f64::total_cmp);
    assert_rel(areas[0], lens, 1e-12);
    assert_rel(areas[1], PI - lens, 1e-12);
    assert_rel(areas[2], PI - lens, 1e-12);
    assert_eq!(a.vertices.len(), 2);
    assert_eq!(a.edges.len(), 4);
    for x in &a.vertices {
        assert_eq!(x.curves, vec![0, 1]);
        assert!((x.point.x - 0.5).abs() < 1e-14);
        for (ci, &t) in x.curves.iter().zip(&x.params) {
            assert!((0.0..TAU).contains(&t));
            let ProfileCurve::Circle { center, radius } = &c[*ci] else {
                unreachable!()
            };
            assert!((*center + *radius * DVec2::from_angle(t) - x.point).length() < 1e-14);
        }
    }
    extrude_all(&c, &a);
    // The union of all three regions is the outline of both discs.
    let all: Vec<usize> = (0..3).collect();
    let prof = k().region_profile(&a, &c, &all).unwrap();
    assert_eq!(prof.len(), 2);
    let body = extrude_profile(&prof, 1.0);
    assert_rel(volume(&body), 2.0 * PI - lens, 1e-9);
}

#[test]
fn internally_tangent_circles() {
    let c = vec![circle(v(0.0, 0.0), 2.0), circle(v(1.0, 0.0), 1.0)];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 2);
    let areas = sorted_areas(&a);
    assert_rel(areas[0], PI, 1e-12);
    assert_rel(areas[1], 3.0 * PI, 1e-12);
    assert_eq!(a.vertices.len(), 1);
    assert!((a.vertices[0].point - v(2.0, 0.0)).length() < 1e-12);
    // The inner disc extrudes; the crescent touches itself at the tangency
    // (a non-manifold profile), which extrude rejects.
    let inner = a
        .regions
        .iter()
        .position(|r| (r.area - PI).abs() < 1e-9)
        .unwrap();
    extrude_region(&c, &a, &[inner]);
    // Both together are the big disc.
    let prof = k().region_profile(&a, &c, &[0, 1]).unwrap();
    assert_rel(volume(&extrude_profile(&prof, 1.0)), 4.0 * PI, 1e-9);
}

#[test]
fn externally_tangent_circles() {
    let c = vec![circle(v(0.0, 0.0), 1.0), circle(v(2.0, 0.0), 1.0)];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 2);
    for r in &a.regions {
        assert_rel(r.area, PI, 1e-12);
    }
    assert_eq!(a.vertices.len(), 1);
    extrude_all(&c, &a);
}

#[test]
fn line_tangent_to_circle_dangles() {
    let c = vec![circle(v(0.0, 0.0), 1.0), line(v(-2.0, 1.0), v(2.0, 1.0))];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 1);
    assert_rel(a.regions[0].area, PI, 1e-12);
    // The tangency is a vertex; the dangling line is pruned.
    assert!(a.edges.iter().all(|e| e.curve == 0));
    assert!(a
        .vertices
        .iter()
        .any(|x| x.curves == vec![0, 1] && (x.point - v(0.0, 1.0)).length() < 1e-12));
    extrude_all(&c, &a);
}

#[test]
fn circle_inscribed_in_square_makes_no_slivers() {
    let mut c = polygon(&[v(-1.0, -1.0), v(1.0, -1.0), v(1.0, 1.0), v(-1.0, 1.0)]);
    c.push(circle(v(0.0, 0.0), 1.0));
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 5);
    let areas = sorted_areas(&a);
    for x in &areas[..4] {
        assert_rel(*x, (4.0 - PI) / 4.0, 1e-12);
    }
    assert_rel(areas[4], PI, 1e-12);
    extrude_all(&c, &a);
}

#[test]
fn overlapping_collinear_segments_have_no_duplicate_edges() {
    // Rectangle 3 x 1 whose bottom side is drawn as two overlapping lines.
    let c = vec![
        line(v(0.0, 0.0), v(2.0, 0.0)),
        line(v(1.0, 0.0), v(3.0, 0.0)),
        line(v(3.0, 0.0), v(3.0, 1.0)),
        line(v(3.0, 1.0), v(0.0, 1.0)),
        line(v(0.0, 1.0), v(0.0, 0.0)),
    ];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 1);
    assert_rel(a.regions[0].area, 3.0, 1e-12);
    assert_no_duplicate_edges(&c, &a);
    // Bottom: [0,1] and [1,2] of line 0, [2,3] of line 1.
    let bottom: Vec<_> = a.edges.iter().filter(|e| e.curve <= 1).collect();
    assert_eq!(bottom.len(), 3);
    assert_eq!(bottom.iter().filter(|e| e.curve == 1).count(), 1);
    extrude_all(&c, &a);
}

#[test]
fn coincident_arcs_and_circles_have_no_duplicate_edges() {
    let c = vec![
        circle(v(0.0, 0.0), 1.0),
        arc(v(0.0, 0.0), 1.0, 0.0, FRAC_PI_2),
        circle(v(0.0, 0.0), 1.0),
        arc(v(0.0, 0.0), 1.0, PI, -FRAC_PI_2),
    ];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 1);
    assert_rel(a.regions[0].area, PI, 1e-12);
    assert!(a.edges.iter().all(|e| e.curve == 0));
    // The arcs' three distinct end points split the circle.
    assert_eq!(a.edges.len(), 3);
    assert_no_duplicate_edges(&c, &a);
    extrude_all(&c, &a);
}

#[test]
fn dangling_lines_produce_no_regions() {
    let a = arrange_checked(&[line(v(0.0, 0.0), v(1.0, 0.0))]);
    assert!(a.regions.is_empty() && a.edges.is_empty());
    assert_eq!(a.vertices.len(), 2);
    let open = vec![
        line(v(0.0, 0.0), v(1.0, 0.0)),
        line(v(1.0, 0.0), v(1.0, 1.0)),
        line(v(1.0, 1.0), v(0.0, 1.0)),
    ];
    assert!(arrange_checked(&open).regions.is_empty());
    // A square with a tail and a spike inside.
    let mut c = polygon(&[v(0.0, 0.0), v(1.0, 0.0), v(1.0, 1.0), v(0.0, 1.0)]);
    c.push(line(v(1.0, 1.0), v(2.0, 2.0)));
    c.push(line(v(0.0, 0.0), v(0.5, 0.25)));
    c.push(line(v(0.25, 0.75), v(0.75, 0.5)));
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 1);
    assert_rel(a.regions[0].area, 1.0, 1e-12);
    assert!(a.edges.iter().all(|e| e.curve < 4));
    extrude_all(&c, &a);
}

#[test]
fn nested_shapes_get_holes() {
    let mut c = vec![circle(v(0.0, 0.0), 1.0)];
    c.extend(polygon(&[
        v(-2.0, -2.0),
        v(2.0, -2.0),
        v(2.0, 2.0),
        v(-2.0, 2.0),
    ]));
    c.push(circle(v(0.0, 0.0), 4.0));
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 3);
    let by_area = |x: f64| {
        a.regions
            .iter()
            .find(|r| (r.area - x).abs() < 1e-9 * x)
            .unwrap_or_else(|| panic!("no region of area {x}"))
    };
    assert!(by_area(PI).holes.is_empty());
    assert_eq!(by_area(16.0 - PI).holes.len(), 1);
    assert_eq!(by_area(16.0 * PI - 16.0).holes.len(), 1);
    extrude_all(&c, &a);
    // Disjoint islands.
    let c = vec![circle(v(0.0, 0.0), 1.0), circle(v(5.0, 0.0), 2.0)];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 2);
    assert!(a.regions.iter().all(|r| r.holes.is_empty()));
    extrude_all(&c, &a);
}

#[test]
fn spline_crossing_a_line_several_times() {
    let pts: Vec<DVec2> = (0..=60)
        .map(|i| {
            let x = 3.0 * PI * i as f64 / 60.0;
            v(x, x.sin())
        })
        .collect();
    let s = ProfileCurve::interpolate(&pts).unwrap();
    let c = vec![s.clone(), line(v(-1.0, 0.0), v(3.0 * PI + 1.0, 0.0))];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 3, "{:?}", sorted_areas(&a));
    // Independent estimate: |∫ y dx| between consecutive zeros of the spline.
    let samples = spline_samples(&s, 200_000);
    let total: f64 = samples
        .windows(2)
        .map(|w| 0.5 * (w[0].y.abs() + w[1].y.abs()) * (w[1].x - w[0].x).abs())
        .sum();
    let sum: f64 = a.regions.iter().map(|r| r.area).sum();
    assert!((sum - total).abs() < 1e-6, "{sum} vs {total}");
    for r in &a.regions {
        assert!((r.area - 2.0).abs() < 1e-2, "lobe area {}", r.area);
    }
    extrude_all(&c, &a);
}

#[test]
fn self_intersecting_figure_eight_spline() {
    // Lemniscate of Gerono: x = cos t, y = sin t cos t; each lobe has area 2/3.
    let pts: Vec<DVec2> = (0..=64)
        .map(|i| {
            let t = TAU * i as f64 / 64.0;
            v(t.cos(), t.sin() * t.cos())
        })
        .collect();
    let s = ProfileCurve::interpolate(&pts).unwrap();
    let c = vec![s];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 2, "{:?}", sorted_areas(&a));
    for r in &a.regions {
        assert!((r.area - 2.0 / 3.0).abs() < 1e-3, "lobe area {}", r.area);
    }
    let x = a
        .vertices
        .iter()
        .find(|x| x.point.length() < 1e-3)
        .expect("crossing vertex");
    assert_eq!(x.curves, vec![0]);
    extrude_all(&c, &a);
}

fn index_of_area(a: &SketchArrangement, x: f64) -> usize {
    a.regions
        .iter()
        .position(|r| (r.area - x).abs() <= 1e-9 * x)
        .unwrap_or_else(|| panic!("no region of area {x}: {:?}", sorted_areas(a)))
}

#[test]
fn region_containing_locates_points() {
    let t = tol();
    let c = vec![circle(v(0.0, 0.0), 3.0), circle(v(0.0, 0.0), 1.0)];
    let a = arrange_checked(&c);
    let ring = index_of_area(&a, 8.0 * PI);
    let disc = index_of_area(&a, PI);
    assert_eq!(a.regions[ring].holes.len(), 1);
    let at = |p: DVec2| region_containing(&a, &c, p, &Tolerance::default());
    assert_eq!(at(v(2.0, 0.0)), Some(ring));
    assert_eq!(at(v(0.0, -2.5)), Some(ring));
    // Inside the ring's hole: the disc filling it, not the ring.
    assert_eq!(at(v(0.0, 0.0)), Some(disc));
    assert_eq!(at(v(0.5, 0.3)), Some(disc));
    // Outside everything.
    assert_eq!(at(v(4.0, 0.0)), None);
    assert_eq!(at(v(-2.2, 2.2)), None);
    assert_eq!(at(v(1e9, -1e9)), None);
    // On and near the boundaries.
    for p in [v(1.0, 0.0), v(0.0, -1.0), v(3.0, 0.0), v(-3.0, 0.0)] {
        assert_eq!(at(p), None, "{p}");
        assert_eq!(at(p * (1.0 + 0.5 * t / p.length())), None, "{p}");
        assert_eq!(at(p * (1.0 - 0.5 * t / p.length())), None, "{p}");
    }
    assert_eq!(at(v(1.0 + 3.0 * t, 0.0)), Some(ring));
    assert_eq!(at(v(1.0 - 3.0 * t, 0.0)), Some(disc));
    assert_eq!(at(v(0.0, 3.0 - 3.0 * t)), Some(ring));
    assert_eq!(at(v(0.0, 3.0 + 3.0 * t)), None);
    // Unusable input.
    assert_eq!(at(v(f64::NAN, 0.0)), None);
    assert_eq!(
        region_containing(&a, &[], v(2.0, 0.0), &Tolerance::default()),
        None
    );
    let other = vec![
        line(v(0.0, 0.0), v(1.0, 0.0)),
        line(v(0.0, 0.0), v(0.0, 1.0)),
    ];
    assert_eq!(
        region_containing(&a, &other, v(2.0, 0.0), &Tolerance::default()),
        None
    );
    let zero = Tolerance {
        linear: 0.0,
        angular: 1e-9,
    };
    assert_eq!(region_containing(&a, &c, v(2.0, 0.0), &zero), None);

    // Four triangles meeting at a vertex.
    let (w, h) = (4.0, 3.0);
    let mut c = polygon(&[v(0.0, 0.0), v(w, 0.0), v(w, h), v(0.0, h)]);
    c.push(line(v(0.0, 0.0), v(w, h)));
    c.push(line(v(w, 0.0), v(0.0, h)));
    let a = arrange_checked(&c);
    let at = |p: DVec2| region_containing(&a, &c, p, &Tolerance::default());
    let found: Vec<Option<usize>> = [v(2.0, 0.5), v(3.5, 1.5), v(2.0, 2.5), v(0.5, 1.5)]
        .into_iter()
        .map(at)
        .collect();
    let mut distinct: Vec<usize> = found
        .iter()
        .map(|x| x.expect("inside a triangle"))
        .collect();
    distinct.sort_unstable();
    distinct.dedup();
    assert_eq!(distinct.len(), 4);
    assert_eq!(at(v(w / 2.0, h / 2.0)), None, "the crossing vertex");
    assert_eq!(at(v(1.0, 0.75)), None, "on a diagonal");
    assert_eq!(at(v(1.0, 0.75 + 0.5 * t)), None, "near a diagonal");
    assert_eq!(at(v(0.0, 1.0)), None, "on a side");
}

#[test]
fn invalid_curves_are_reported_by_index() {
    let good = || line(v(0.0, 0.0), v(1.0, 0.0));
    let cases: Vec<(usize, ProfileCurve)> = vec![
        (1, line(v(1.0, 1.0), v(1.0, 1.0))),
        (1, line(v(1.0, 1.0), v(1.0 + 0.5 * tol(), 1.0))),
        (1, line(v(f64::NAN, 1.0), v(1.0, 1.0))),
        (1, line(v(f64::INFINITY, 1.0), v(1.0, 1.0))),
        (1, line(v(1e300, 1.0), v(1.0, 1.0))),
        (1, circle(v(0.0, 0.0), 0.0)),
        (1, circle(v(0.0, 0.0), -1.0)),
        (1, circle(v(0.0, f64::NAN), 1.0)),
        (1, arc(v(0.0, 0.0), 1.0, 0.0, 0.0)),
        (1, arc(v(0.0, 0.0), 1.0, 0.0, 7.0)),
        (1, arc(v(0.0, 0.0), 1.0, f64::INFINITY, 1.0)),
        (
            1,
            ProfileCurve::Spline {
                degree: 2,
                knots: vec![0.0, 0.0, 1.0, 0.5, 1.0, 1.0],
                points: vec![v(0.0, 0.0), v(1.0, 1.0), v(2.0, 0.0)],
                weights: None,
            },
        ),
        (
            1,
            ProfileCurve::Spline {
                degree: 2,
                knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                points: vec![v(0.0, 0.0), v(1.0, 1.0), v(2.0, 0.0)],
                weights: Some(vec![1.0, -1.0, 1.0]),
            },
        ),
        (
            1,
            ProfileCurve::Spline {
                degree: 2,
                knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                points: vec![v(0.0, 0.0), v(1.0, 1.0)],
                weights: None,
            },
        ),
        (
            1,
            ProfileCurve::Spline {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                points: vec![v(2.0, 2.0), v(2.0, 2.0)],
                weights: None,
            },
        ),
    ];
    for (n, (_, bad)) in cases.into_iter().enumerate() {
        let i = n % 3;
        let mut c = vec![good(), good(), circle(v(0.0, 0.0), 1.0)];
        c[i] = bad.clone();
        match k().sketch_regions(&c) {
            Err(KernelError::InvalidInput(m)) => {
                assert!(m.contains(&format!("curve {i}")), "{bad:?}: {m}")
            }
            other => panic!("{bad:?}: {other:?}"),
        }
    }
    for linear in [0.0, -1.0, f64::NAN] {
        let t = Tolerance {
            linear,
            angular: 1e-9,
        };
        assert!(matches!(
            arrange(&[good()], &t),
            Err(KernelError::InvalidInput(_))
        ));
    }
}

#[test]
fn region_profile_rejects_bad_selections() {
    let c = polygon(&[v(0.0, 0.0), v(1.0, 0.0), v(1.0, 1.0), v(0.0, 1.0)]);
    let a = arrange_checked(&c);
    let bad = |r: Result<Vec<ProfileCurve>, KernelError>| {
        assert!(matches!(r, Err(KernelError::InvalidInput(_))), "{r:?}")
    };
    bad(k().region_profile(&a, &c, &[]));
    bad(k().region_profile(&a, &c, &[1]));
    bad(k().region_profile(&a, &c, &[0, usize::MAX]));
    bad(k().region_profile(&a, &c[..2], &[0]));
    bad(k().region_profile(&a, &[circle(v(0.0, 0.0), 1.0)], &[0]));
    let mut broken = a.clone();
    broken.regions[0].outer.edges.swap(0, 1);
    bad(k().region_profile(&broken, &c, &[0]));
    let mut broken = a.clone();
    broken.edges[0].t1 = 2.0;
    bad(k().region_profile(&broken, &c, &[0]));
    // Duplicate indices are ignored.
    let prof = k().region_profile(&a, &c, &[0, 0]).unwrap();
    assert_eq!(prof.len(), 4);
}

/// Signed area of a profile loop sequence, from fine samples.
fn profile_area(prof: &[ProfileCurve]) -> f64 {
    let pts: Vec<DVec2> = prof
        .iter()
        .flat_map(|c| {
            let (a, b) = domain(c);
            let mut s = piece_samples(c, a, b, true, 512);
            s.pop();
            s
        })
        .collect();
    shoelace(&pts)
}

/// The curves of a profile join end to end into closed loops, in order.
fn assert_chained(prof: &[ProfileCurve]) {
    let ends: Vec<(DVec2, DVec2)> = prof
        .iter()
        .map(|c| {
            let (a, b) = domain(c);
            (curve_point(c, a), curve_point(c, b))
        })
        .collect();
    let n = ends.len();
    let mut start = 0;
    for i in 0..n {
        if (ends[i].1 - ends[start].0).length() < 1e-9 {
            start = i + 1;
            continue;
        }
        assert!(
            i + 1 < n && (ends[i].1 - ends[i + 1].0).length() < 1e-9,
            "profile curve {i} is not joined to the next"
        );
    }
    assert_eq!(start, n, "profile does not end with a closed loop");
}

#[test]
fn unions_of_regions_drop_shared_edges() {
    let (w, h) = (4.0, 3.0);
    let mut c = polygon(&[v(0.0, 0.0), v(w, 0.0), v(w, h), v(0.0, h)]);
    c.push(line(v(0.0, 0.0), v(w, h)));
    c.push(line(v(w, 0.0), v(0.0, h)));
    let a = arrange_checked(&c);
    let at = |p: DVec2| region_containing(&a, &c, p, &Tolerance::default()).unwrap();
    let (bottom, right, top) = (at(v(2.0, 0.5)), at(v(3.5, 1.5)), at(v(2.0, 2.5)));
    let prof = k().region_profile(&a, &c, &[bottom, right]).unwrap();
    assert_eq!(prof.len(), 4);
    assert_chained(&prof);
    assert_rel(profile_area(&prof), w * h / 2.0, 1e-12);
    extrude_region(&c, &a, &[bottom, right]);
    let all: Vec<usize> = (0..4).collect();
    let prof = k().region_profile(&a, &c, &all).unwrap();
    assert_eq!(prof.len(), 4);
    assert!(prof.iter().all(|p| {
        let ProfileCurve::Line { start, end } = p else {
            return false;
        };
        start.x == end.x || start.y == end.y
    }));
    extrude_region(&c, &a, &all);
    // Opposite triangles touch only at the centre: a pinched union.
    let prof = k().region_profile(&a, &c, &[bottom, top]).unwrap();
    assert_eq!(prof.len(), 6);
    assert_rel(profile_area(&prof), w * h / 2.0, 1e-12);
    if let Err(e) = try_extrude(&prof, 1.0) {
        assert!(matches!(e, KernelError::InvalidInput(_)), "{e}");
    }

    // A ring and the disc filling its hole: one disc, no hole.
    let c = vec![circle(v(0.0, 0.0), 3.0), circle(v(1.0, 0.0), 1.0)];
    let a = arrange_checked(&c);
    let prof = k().region_profile(&a, &c, &[0, 1]).unwrap();
    assert_eq!(prof.len(), 2, "the closed circle edge leaves as two halves");
    assert_chained(&prof);
    assert_rel(profile_area(&prof), 9.0 * PI, 1e-5);
    extrude_region(&c, &a, &[0, 1]);
    // The ring alone: outer loop counter-clockwise, hole clockwise.
    let ring = index_of_area(&a, 8.0 * PI);
    let prof = k().region_profile(&a, &c, &[ring]).unwrap();
    assert_eq!(prof.len(), 4);
    assert_rel(profile_area(&prof[..2]), 9.0 * PI, 1e-5);
    assert_rel(profile_area(&prof[2..]), -PI, 1e-5);
    extrude_region(&c, &a, &[ring]);
}

#[test]
fn bridges_bound_nothing() {
    // Two squares joined by a segment.
    let mut c = polygon(&[v(0.0, 0.0), v(1.0, 0.0), v(1.0, 1.0), v(0.0, 1.0)]);
    c.extend(polygon(&[
        v(3.0, 0.0),
        v(4.0, 0.0),
        v(4.0, 1.0),
        v(3.0, 1.0),
    ]));
    c.push(line(v(1.0, 0.5), v(3.0, 0.5)));
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 2);
    let bridge: Vec<_> = a.edges.iter().filter(|e| e.curve == 8).collect();
    assert_eq!(bridge.len(), 1);
    assert_eq!((bridge[0].left, bridge[0].right), (None, None));
    extrude_all(&c, &a);
    let both = k().region_profile(&a, &c, &[0, 1]).unwrap();
    assert!(both.iter().all(
        |p| !matches!(p, ProfileCurve::Line { start, end } if start.y == 0.5 && end.y == 0.5)
    ));
    assert_rel(volume(&extrude_profile(&both, 1.0)), 2.0, 1e-9);

    // A square inside a square, joined by a segment: the ring's boundary
    // runs along the bridge both ways, and its profile leaves it out.
    let mut c = polygon(&[v(0.0, 0.0), v(4.0, 0.0), v(4.0, 4.0), v(0.0, 4.0)]);
    c.extend(polygon(&[
        v(1.0, 1.0),
        v(3.0, 1.0),
        v(3.0, 3.0),
        v(1.0, 3.0),
    ]));
    c.push(line(v(0.0, 2.0), v(1.0, 2.0)));
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 2);
    let ring = index_of_area(&a, 12.0);
    let b = a.edges.iter().position(|e| e.curve == 8).unwrap();
    assert_eq!(
        (a.edges[b].left, a.edges[b].right),
        (Some(ring), Some(ring))
    );
    let prof = k().region_profile(&a, &c, &[ring]).unwrap();
    assert_eq!(prof.len(), 10);
    assert_chained(&prof);
    assert_rel(profile_area(&prof), 12.0, 1e-12);
    assert_eq!(
        region_containing(&a, &c, v(0.5, 2.0), &Tolerance::default()),
        None
    );
    assert_eq!(
        region_containing(&a, &c, v(0.5, 2.1), &Tolerance::default()),
        Some(ring)
    );
    extrude_all(&c, &a);
}

#[test]
fn vertex_params_follow_the_conventions() {
    // Upper half circle drawn clockwise from (-1, 0), its diameter, a
    // vertical line through both, and a circle crossing at angle 3π/2.
    let c = vec![
        arc(v(0.0, 0.0), 1.0, PI, -PI),
        line(v(1.0, 0.0), v(-1.0, 0.0)),
        line(v(0.0, -1.0), v(0.0, 2.0)),
        circle(v(0.0, 0.5), 1.5),
    ];
    let a = arrange_checked(&c);
    let at = |p: DVec2| {
        a.vertices
            .iter()
            .find(|x| (x.point - p).length() < 1e-12)
            .unwrap_or_else(|| panic!("no vertex at {p}"))
    };
    let top = at(v(0.0, 1.0));
    assert_eq!(top.curves, vec![0, 2]);
    assert_rel(top.params[0], FRAC_PI_2, 1e-14);
    assert_rel(top.params[1], 2.0 / 3.0, 1e-14);
    let bottom = at(v(0.0, -1.0));
    assert_eq!(bottom.curves, vec![2, 3]);
    assert_eq!(bottom.params[0], 0.0);
    assert_rel(bottom.params[1], 3.0 * FRAC_PI_2, 1e-14);
    let mid = at(v(0.0, 0.0));
    assert_eq!(mid.curves, vec![1, 2]);
    assert_rel(mid.params[0], 0.5, 1e-14);
    assert_rel(mid.params[1], 1.0 / 3.0, 1e-14);
    let end = at(v(1.0, 0.0));
    assert_eq!(end.curves, vec![0, 1]);
    assert_rel(end.params[0], PI, 1e-14);
    assert_eq!(end.params[1], 0.0);
    extrude_all(&c, &a);
}

#[test]
fn circles_tangent_at_one_point() {
    // Three circles tangent at the origin on the same side, and a line
    // tangent to all of them there: ordered by curvature.
    let c = vec![
        circle(v(1.0, 0.0), 1.0),
        circle(v(2.0, 0.0), 2.0),
        circle(v(3.0, 0.0), 3.0),
        line(v(0.0, -4.0), v(0.0, 4.0)),
    ];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 3);
    let areas = sorted_areas(&a);
    assert_rel(areas[0], PI, 1e-12);
    assert_rel(areas[1], 3.0 * PI, 1e-12);
    assert_rel(areas[2], 5.0 * PI, 1e-12);
    // The origin, and the dangling line's ends.
    assert_eq!(a.vertices.len(), 3);
    let o = a
        .vertices
        .iter()
        .find(|x| x.point.length() < 1e-12)
        .unwrap();
    assert_eq!(o.curves, vec![0, 1, 2, 3]);
    assert!(a.edges.iter().all(|e| e.curve < 3));
    extrude_all(&c, &a);
    // Circles tangent from opposite sides, plus one through the point.
    let c = vec![
        circle(v(-1.0, 0.0), 1.0),
        circle(v(2.0, 0.0), 2.0),
        circle(v(0.0, 1.0), 1.0),
    ];
    let a = arrange_checked(&c);
    let lens = |r1: f64, r2: f64, d: f64| {
        let a1 = ((d * d + r1 * r1 - r2 * r2) / (2.0 * d * r1)).acos();
        let a2 = ((d * d + r2 * r2 - r1 * r1) / (2.0 * d * r2)).acos();
        r1 * r1 * (a1 - a1.sin() * a1.cos()) + r2 * r2 * (a2 - a2.sin() * a2.cos())
    };
    let l1 = lens(1.0, 1.0, 2f64.sqrt());
    let l2 = lens(2.0, 1.0, 5f64.sqrt());
    let total = PI + 4.0 * PI + PI - l1 - l2;
    let sum: f64 = a.regions.iter().map(|r| r.area).sum();
    assert_rel(sum, total, 1e-12);
    assert_eq!(a.regions.len(), 5);
    extrude_all(&c, &a);
}

#[test]
fn rational_circle_coincides_with_circle() {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let (r, o) = (1.5, v(0.5, -0.25));
    let pts: Vec<DVec2> = [
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 0.0),
        (-1.0, -1.0),
        (0.0, -1.0),
        (1.0, -1.0),
        (1.0, 0.0),
    ]
    .iter()
    .map(|&(x, y)| o + r * v(x, y))
    .collect();
    let nurbs = ProfileCurve::Spline {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        points: pts,
        weights: Some(vec![1.0, s, 1.0, s, 1.0, s, 1.0, s, 1.0]),
    };
    // On its own: the exact disc.
    let c = vec![nurbs.clone()];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 1);
    assert_rel(a.regions[0].area, PI * r * r, 1e-12);
    extrude_all(&c, &a);
    // Together with the same circle: one region, one copy of each edge.
    let c = vec![circle(o, r), nurbs, line(o - v(r, 0.0), o + v(r, 0.0))];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 2, "{:?}", sorted_areas(&a));
    for x in sorted_areas(&a) {
        assert_rel(x, 0.5 * PI * r * r, 1e-12);
    }
    assert_no_duplicate_edges(&c, &a);
    assert!(
        a.edges.iter().all(|e| e.curve != 1),
        "the circle (lower index) is kept"
    );
    assert_eq!(a.edges.len(), 3);
    assert_eq!(a.vertices.len(), 2);
    extrude_all(&c, &a);
}

#[test]
fn identical_splines_and_spline_tangencies() {
    let pts: Vec<DVec2> = (0..=12)
        .map(|i| {
            let x = -1.0 + i as f64 / 6.0;
            v(x, x * x + 0.3 * (3.0 * x).sin())
        })
        .collect();
    let s = ProfileCurve::interpolate(&pts).unwrap();
    let (p0, p1) = (pts[0], pts[12]);
    // The same spline twice, closed by a line: one region.
    let c = vec![s.clone(), s.clone(), line(p1, p0)];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 1, "{:?}", sorted_areas(&a));
    assert_no_duplicate_edges(&c, &a);
    assert!(a.edges.iter().all(|e| e.curve != 1));
    // No leftover vertices inside the overlap: the spline is one edge.
    assert_eq!(a.edges.len(), 2);
    assert_eq!(a.vertices.len(), 2);
    extrude_all(&c, &a);
    // Reversed copy.
    let ProfileCurve::Spline {
        degree,
        knots,
        points,
        weights,
    } = s.clone()
    else {
        unreachable!()
    };
    let rev = ProfileCurve::Spline {
        degree,
        knots: knots.iter().rev().map(|u| 1.0 - u).collect(),
        points: points.into_iter().rev().collect(),
        weights,
    };
    let c = vec![rev, s.clone(), line(p1, p0)];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 1, "{:?}", sorted_areas(&a));
    assert_no_duplicate_edges(&c, &a);
    assert_eq!(a.edges.len(), 2);
    assert!(a.edges.iter().all(|e| e.curve != 1));
    extrude_all(&c, &a);

    // A parabola y = x² (exact quadratic Bézier), the chord y = 1 and the
    // tangent y = 0 at its apex: area 4/3, the tangent line dangles.
    let par = ProfileCurve::Spline {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        points: vec![v(-1.0, 1.0), v(0.0, -1.0), v(1.0, 1.0)],
        weights: None,
    };
    let c = vec![
        par,
        line(v(1.0, 1.0), v(-1.0, 1.0)),
        line(v(-2.0, 0.0), v(2.0, 0.0)),
    ];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 1);
    assert_rel(a.regions[0].area, 4.0 / 3.0, 1e-12);
    // A tangency is located to within the tolerance (here ~1e-8).
    let apex = a
        .vertices
        .iter()
        .find(|x| x.point.length() < tol())
        .expect("tangency vertex");
    assert_eq!(apex.curves, vec![0, 2]);
    assert!((apex.params[0] - 0.5).abs() < 1e-6);
    assert!((apex.params[1] - 0.5).abs() < 1e-6);
    extrude_all(&c, &a);
}

#[test]
fn arcs_and_lines_close_regions() {
    let c = vec![
        arc(v(0.0, 0.0), 1.0, 0.0, FRAC_PI_2),
        arc(v(0.0, 0.0), 1.0, PI, -FRAC_PI_2),
        arc(v(0.0, 0.0), 1.0, PI, FRAC_PI_2),
        arc(v(0.0, 0.0), 1.0, 0.0, -FRAC_PI_2),
    ];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 1);
    assert_rel(a.regions[0].area, PI, 1e-12);
    extrude_all(&c, &a);
    // A D shape and a slot (two half circles and two lines).
    let c = vec![
        arc(v(0.0, 0.0), 2.0, -FRAC_PI_2, PI),
        line(v(0.0, 2.0), v(0.0, -2.0)),
        arc(v(5.0, 0.0), 0.5, FRAC_PI_2, PI),
        line(v(5.0, -0.5), v(7.0, -0.5)),
        arc(v(7.0, 0.0), 0.5, -FRAC_PI_2, PI),
        line(v(7.0, 0.5), v(5.0, 0.5)),
    ];
    let a = arrange_checked(&c);
    let areas = sorted_areas(&a);
    assert_eq!(areas.len(), 2);
    assert_rel(areas[0], 0.25 * PI + 2.0, 1e-12);
    assert_rel(areas[1], 2.0 * PI, 1e-12);
    extrude_all(&c, &a);
}

#[test]
fn edges_below_the_profile_resolution_are_absorbed() {
    // A full-turn arc and a line ending 1.2 tolerances from where it touches
    // the arc: two vertices 1.2 tolerances apart, joined by a tiny arc
    // piece, which the profile builder (merging ends within ten
    // tolerances) cannot take. The profile leaves it out.
    let t = tol();
    let c = vec![
        arc(v(0.0, 3.0), 1.0, 0.0, -TAU),
        line(v(0.9 * t, 2.0 - 0.3 * t), v(-1.0, 2.0)),
    ];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 1);
    assert_eq!(a.regions[0].outer.edges.len(), 3);
    let tiny = a.edges.iter().filter(|e| e.t1 - e.t0 < 10.0 * t).count();
    assert_eq!(tiny, 1);
    let prof = k().region_profile(&a, &c, &[0]).unwrap();
    assert_eq!(prof.len(), 2);
    extrude_all(&c, &a);
    // A disc smaller than the profile resolution has no profile.
    let c = vec![circle(v(0.0, 0.0), 4.0 * t)];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 1);
    assert!(matches!(
        k().region_profile(&a, &c, &[0]),
        Err(KernelError::InvalidInput(_))
    ));
    // At a finer tolerance it has one.
    let fine = Tolerance {
        linear: 1e-9,
        angular: 1e-9,
    };
    let a = arrange(&c, &fine).unwrap();
    let prof = region_profile_with_tolerance(&a, &c, &[0], &fine).unwrap();
    assert_eq!(prof.len(), 2);
    let bad = Tolerance {
        linear: f64::NAN,
        angular: 1e-9,
    };
    assert!(matches!(
        region_profile_with_tolerance(&a, &c, &[0], &bad),
        Err(KernelError::InvalidInput(_))
    ));
}

#[test]
fn chains_of_edges_below_the_profile_resolution_are_rejected() {
    // The unit square with part of its bottom side drawn as a hundred
    // segments of five tolerances each: every one of them is degenerate to
    // the profile builder and together they span far more than its
    // resolution, so leaving them out would change the region. The profile
    // is refused instead of being wrong.
    let t = tol();
    let x0 = 0.5;
    let x1 = x0 + 100.0 * 5.0 * t;
    let mut c = vec![line(v(0.0, 0.0), v(x0, 0.0))];
    c.extend((0..100).map(|i| {
        let a = x0 + 5.0 * t * i as f64;
        line(v(a, 0.0), v(if i == 99 { x1 } else { a + 5.0 * t }, 0.0))
    }));
    c.extend(polygon(&[
        v(x1, 0.0),
        v(1.0, 0.0),
        v(1.0, 1.0),
        v(0.0, 1.0),
    ]));
    c.pop(); // polygon() closes back to (x1, 0); close at the origin instead.
    c.push(line(v(0.0, 1.0), v(0.0, 0.0)));
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 1);
    assert_eq!(a.edges.len(), 105);
    assert_rel(a.regions[0].area, 1.0, 1e-12);
    match k().region_profile(&a, &c, &[0]) {
        Err(KernelError::InvalidInput(m)) => {
            assert!(m.contains("shorter than the profile resolution"), "{m}")
        }
        other => panic!("expected InvalidInput, got {other:?}"),
    }
    // Drawn as one segment, the side is fine.
    let mut c2 = vec![line(v(0.0, 0.0), v(x1, 0.0))];
    c2.extend_from_slice(&c[101..]);
    let a2 = arrange_checked(&c2);
    assert_eq!(a2.regions.len(), 1);
    extrude_region(&c2, &a2, &[0]);
}

#[test]
fn duplicate_curves_do_not_move_vertices() {
    // A closed spline drawn twice, crossing a circle. The two copies touch
    // everywhere, which adds overlap contacts placed only to about the
    // tolerance along the curves; the vertices must still be the exact
    // crossings, and the region profiles must join exactly (a joint gap of
    // a fraction of the tolerance makes extrude bend the curves, changing
    // the volume).
    let c0 = v(0.5705899629558933, 0.19297138976103784);
    let star = [
        (0.1835120869843455, 0.7202327873321246),
        (0.19314654623395242, 0.6636225695002943),
        (0.29478914146621615, 1.194142843604957),
        (0.34026620678644726, 1.2866370508855403),
        (0.2067259846484217, 0.9804325940698583),
    ];
    let mut pts: Vec<DVec2> = star
        .iter()
        .enumerate()
        .map(|(k, &(f, r))| c0 + r * DVec2::from_angle(TAU * (k as f64 + 0.8 * f) / 5.0))
        .collect();
    pts.push(pts[0]);
    let sp = ProfileCurve::interpolate(&pts).unwrap();
    let disc = circle(
        v(0.3703018374793569, -1.0380711104343028),
        1.0763643642661234,
    );
    let single = arrange_checked(&[disc.clone(), sp.clone()]);
    let c = vec![disc, sp.clone(), sp];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 3);
    assert_eq!(a.vertices.len(), single.vertices.len());
    for (x, y) in a.vertices.iter().zip(&single.vertices) {
        assert!((x.point - y.point).length() < 1e-12, "{x:?} vs {y:?}");
    }
    for r in 0..a.regions.len() {
        let prof = k().region_profile(&a, &c, &[r]).unwrap();
        for (i, cur) in prof.iter().enumerate() {
            let next = &prof[(i + 1) % prof.len()];
            let gap = curve_point(cur, domain(cur).1) - curve_point(next, domain(next).0);
            assert!(gap.length() < 1e-12, "region {r} joint {i}: gap {gap}");
        }
    }
    assert_eq!(extrude_all(&c, &a), 3);
}

#[test]
fn spline_overlapping_part_of_another_spline() {
    use openshape_nurbs::NurbsCurve;
    let pts: Vec<DVec2> = (0..=12)
        .map(|i| {
            let x = -1.0 + i as f64 / 6.0;
            v(x, x * x + 0.3 * (3.0 * x).sin())
        })
        .collect();
    let s = ProfileCurve::interpolate(&pts).unwrap();
    let ProfileCurve::Spline {
        degree,
        knots,
        points,
        ..
    } = &s
    else {
        unreachable!()
    };
    // The middle third of the spline as a spline of its own, with its own
    // knots (split by knot insertion).
    let n = NurbsCurve::new(
        *degree,
        knots.clone(),
        points.iter().map(|p| p.extend(0.0)).collect(),
        None,
    )
    .unwrap();
    let (lo, hi) = n.domain();
    let rest = n.try_split(lo + (hi - lo) / 3.0).unwrap().1;
    let mid = rest.try_split(lo + 2.0 * (hi - lo) / 3.0).unwrap().0;
    let part = ProfileCurve::Spline {
        degree: mid.degree(),
        knots: mid.knots().as_slice().to_vec(),
        points: mid.control_points().iter().map(|p| p.truncate()).collect(),
        weights: None,
    };
    let chord = line(pts[12], pts[0]);
    let alone = arrange_checked(&[s.clone(), chord.clone()]);
    assert_eq!(alone.regions.len(), 1);
    let c = vec![part, s, chord];
    let a = arrange_checked(&c);
    assert_eq!(a.regions.len(), 1, "{:?}", sorted_areas(&a));
    assert_rel(a.regions[0].area, alone.regions[0].area, 1e-12);
    assert_no_duplicate_edges(&c, &a);
    // The overlap is kept on the part (lower index), the rest on the spline.
    assert_eq!(a.edges.iter().filter(|e| e.curve == 0).count(), 1);
    assert_eq!(a.edges.iter().filter(|e| e.curve == 1).count(), 2);
    extrude_all(&c, &a);
}
