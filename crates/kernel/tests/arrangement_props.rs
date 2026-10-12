//! Property tests for sketch regions (design note 0006, §6) against
//! independent oracles.

mod arrangement_util;

use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_2, TAU};
use std::sync::atomic::{AtomicUsize, Ordering};

use arrangement_util::oracle::{in_polygon, polyline_distance, CellGrid, Shape, Shapes};
use arrangement_util::raster::Raster;
use arrangement_util::*;
use openshape_kernel::arrangement::region_containing;
use openshape_kernel::{DVec2, Kernel, KernelError, ProfileCurve, SketchArrangement, Tolerance};
use proptest::prelude::*;
use proptest::test_runner::TestCaseError;

fn vec2(lo: f64, hi: f64) -> impl Strategy<Value = DVec2> {
    (lo..hi, lo..hi).prop_map(|(x, y)| v(x, y))
}

/// A polygon with vertices in angular order about `c`: star-shaped about
/// `c` and counter-clockwise when every angular gap is under π, which
/// jittered triangles need not be (they may leave `c` outside and run
/// clockwise; the oracle takes either direction).
fn star(c: DVec2, pts: &[(f64, f64)]) -> Vec<DVec2> {
    let n = pts.len();
    pts.iter()
        .enumerate()
        .map(|(k, &(f, r))| {
            let a = TAU * (k as f64 + 0.8 * f) / n as f64;
            c + r * DVec2::from_angle(a)
        })
        .collect()
}

fn shape() -> impl Strategy<Value = Shape> {
    prop_oneof![
        (vec2(-2.0, 2.0), 0.2f64..1.5).prop_map(|(c, r)| Shape::Circle(c, r)),
        (
            vec2(-2.0, 2.0),
            prop::collection::vec((0.0f64..1.0, 0.3f64..1.5), 3..8)
        )
            .prop_map(|(c, p)| Shape::Poly(star(c, &p))),
    ]
}

/// Chords of the circles among `shapes`: they lie inside the union, so they
/// split regions without changing the enclosed area.
fn chords(shapes: &[Shape], spec: &[(usize, f64, f64)]) -> Vec<ProfileCurve> {
    spec.iter()
        .filter_map(|&(i, a, b)| match &shapes[i % shapes.len()] {
            Shape::Circle(c, r) if (a - b).abs() > 0.05 => Some(line(
                *c + *r * DVec2::from_angle(a),
                *c + *r * DVec2::from_angle(b),
            )),
            _ => None,
        })
        .collect()
}

/// Polygons of every region's loops (outer first).
fn region_polys(a: &SketchArrangement, curves: &[ProfileCurve]) -> Vec<Vec<Vec<DVec2>>> {
    a.regions
        .iter()
        .map(|r| {
            std::iter::once(&r.outer)
                .chain(&r.holes)
                .map(|l| loop_polygon(a, curves, l, 512))
                .collect()
        })
        .collect()
}

fn in_region(p: DVec2, loops: &[Vec<DVec2>]) -> bool {
    in_polygon(p, &loops[0]) && !loops[1..].iter().any(|h| in_polygon(p, h))
}

/// Every sample point away from the curves lies in exactly the regions the
/// oracle predicts (one if enclosed, none otherwise), and
/// `region_containing` agrees.
fn check_partition(
    a: &SketchArrangement,
    curves: &[ProfileCurve],
    samples: &[DVec2],
    margin: f64,
    enclosed: impl Fn(DVec2) -> bool,
) {
    let polys = region_polys(a, curves);
    for &p in samples {
        let near = polys
            .iter()
            .flatten()
            .any(|l| polyline_distance(p, l) < margin);
        if near {
            continue;
        }
        let inside: Vec<usize> = (0..polys.len())
            .filter(|&r| in_region(p, &polys[r]))
            .collect();
        assert!(inside.len() <= 1, "{p:?} is in regions {inside:?}");
        assert_eq!(
            !inside.is_empty(),
            enclosed(p),
            "{p:?}: in regions {inside:?}"
        );
        assert_eq!(
            region_containing(a, curves, p, &Tolerance::default()),
            inside.first().copied(),
            "region_containing at {p:?}"
        );
    }
}

/// Extrudes every region, allowing a rejection only for a non-manifold
/// (pinched) profile.
fn extrude_every_region(curves: &[ProfileCurve], a: &SketchArrangement) {
    extrude_all(curves, a);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// Random polygons and circles with chords: the regions tile exactly the
    /// area the curves enclose, every point is in the right region, and
    /// every region extrudes to a valid solid of volume area × height.
    #[test]
    fn polygons_and_circles(
        shapes in prop::collection::vec(shape(), 1..6),
        chord_spec in prop::collection::vec((0usize..8, 0.0f64..TAU, 0.0f64..TAU), 0..4),
        samples in prop::collection::vec(vec2(-3.6, 3.6), 48),
    ) {
        let oracle = Shapes::new(shapes.clone());
        let mut curves = oracle.curves();
        curves.extend(chords(&shapes, &chord_spec));
        let a = arrange_checked(&curves);
        let sum: f64 = a.regions.iter().map(|r| r.area).sum();
        let want = oracle.filled_area();
        prop_assert!((sum - want).abs() <= 1e-9 * want, "region areas {sum} vs enclosed area {want}");
        check_partition(&a, &curves, &samples, 1e-4, |p| oracle.filled_contains(p));
        extrude_every_region(&curves, &a);
    }

    /// Axis-aligned rectangles and segments on an integer grid, drawn with
    /// overlapping and duplicated pieces in random directions: exact
    /// degeneracies everywhere (collinear overlaps, T-junctions, shared
    /// corners, pinches, bridges). The faces must match the cell oracle
    /// exactly, at scales 1, 1e-4 and 1e3.
    #[test]
    fn grid_rectangles_and_segments(
        rects in prop::collection::vec((0i64..7, 0i64..7, 1i64..5, 1i64..5), 1..5),
        segs in prop::collection::vec((0i64..9, 0i64..9, 1i64..5, any::<bool>()), 0..5),
        splits in prop::collection::vec((any::<bool>(), any::<bool>(), 0i64..4), 32),
        scale in prop::sample::select(vec![1.0, 1e-4, 1e3]),
    ) {
        let n = 8;
        let mut lines: Vec<(i64, i64, i64, i64)> = Vec::new();
        for &(x, y, w, h) in &rects {
            let (x1, y1) = ((x + w).min(n), (y + h).min(n));
            lines.extend([(x, y, x1, y), (x1, y, x1, y1), (x1, y1, x, y1), (x, y1, x, y)]);
        }
        for &(x, y, l, horiz) in &segs {
            let (x, y) = (x.min(n), y.min(n));
            if horiz {
                lines.push((x, y, (x + l).min(n), y));
            } else {
                lines.push((x, y, x, (y + l).min(n)));
            }
        }
        lines.retain(|&(x0, y0, x1, y1)| (x0, y0) != (x1, y1));
        let grid = CellGrid { n, segs: lines.clone() };
        let mut curves = Vec::new();
        for (k, &(x0, y0, x1, y1)) in lines.iter().enumerate() {
            let (rev, split, at) = splits[k % splits.len()];
            let (a, b) = (v(x0 as f64, y0 as f64), v(x1 as f64, y1 as f64));
            let (a, b) = if rev { (b, a) } else { (a, b) };
            let len = (x1 - x0).abs().max((y1 - y0).abs());
            if split && len >= 2 {
                // Two overlapping pieces covering the segment.
                let m = (at % (len - 1)) + 1;
                let d = (b - a) / len as f64;
                curves.push(line(a * scale, (a + d * (m + 1).min(len) as f64) * scale));
                curves.push(line((a + d * (m - 1) as f64) * scale, b * scale));
            } else {
                curves.push(line(a * scale, b * scale));
            }
        }
        let a = arrange_checked(&curves);
        let faces = grid.faces();
        prop_assert_eq!(a.regions.len(), faces.len(), "region count");
        let mut want: Vec<f64> = faces.iter().map(|f| f.len() as f64 * scale * scale).collect();
        want.sort_by(f64::total_cmp);
        let got = sorted_areas(&a);
        for (g, w) in got.iter().zip(&want) {
            prop_assert!((g - w).abs() <= 1e-9 * w, "areas {got:?} vs {want:?}");
        }
        // Each region's interior point is in a face of the same area, and
        // no face is matched twice.
        let mut matched = vec![false; faces.len()];
        for r in &a.regions {
            let p = r.interior_point / scale;
            let cell = (p.x.floor() as i64, p.y.floor() as i64);
            let f = faces.iter().position(|f| f.contains(&cell));
            prop_assert!(f.is_some(), "interior point {:?} is outside every face", p);
            let f = f.unwrap();
            prop_assert!(!matched[f], "two regions in one face");
            matched[f] = true;
            prop_assert!((faces[f].len() as f64 * scale * scale - r.area).abs() <= 1e-9 * r.area);
        }
        extrude_every_region(&curves, &a);
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    /// Scaling the sketch by 1e-4 or 1e3 scales areas by the square and
    /// keeps the topology, for inputs whose features are well above the
    /// tolerance at every scale.
    #[test]
    fn scale_invariance(
        shapes in prop::collection::vec(shape(), 1..5),
        chord_spec in prop::collection::vec((0usize..8, 0.0f64..TAU, 0.0f64..TAU), 0..3),
    ) {
        let oracle = Shapes::new(shapes.clone());
        let mut curves = oracle.curves();
        curves.extend(chords(&shapes, &chord_spec));
        let a = arrange_checked(&curves);
        prop_assume!(feature_size(&a, &curves) > 0.02);
        for s in [1e-4, 1e3] {
            let scaled: Vec<ProfileCurve> = curves.iter().map(|c| scale_curve(c, s)).collect();
            let b = arrange_checked(&scaled);
            prop_assert_eq!(a.vertices.len(), b.vertices.len());
            prop_assert_eq!(a.edges.len(), b.edges.len());
            prop_assert_eq!(a.regions.len(), b.regions.len());
            let (x, y) = (sorted_areas(&a), sorted_areas(&b));
            for (p, q) in x.iter().zip(&y) {
                prop_assert!((q / (s * s) - p).abs() <= 1e-9 * p, "{x:?} vs {y:?} at scale {s}");
            }
            for (e, f) in a.edges.iter().zip(&b.edges) {
                prop_assert_eq!((e.curve, e.start, e.end, e.left, e.right), (f.curve, f.start, f.end, f.left, f.right));
            }
            extrude_every_region(&scaled, &b);
        }
    }
}

fn scale_curve(c: &ProfileCurve, s: f64) -> ProfileCurve {
    match c {
        ProfileCurve::Line { start, end } => line(*start * s, *end * s),
        ProfileCurve::Circle { center, radius } => circle(*center * s, radius * s),
        ProfileCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => arc(*center * s, radius * s, *start_angle, *sweep_angle),
        ProfileCurve::Spline {
            degree,
            knots,
            points,
            weights,
        } => ProfileCurve::Spline {
            degree: *degree,
            knots: knots.clone(),
            points: points.iter().map(|p| *p * s).collect(),
            weights: weights.clone(),
        },
    }
}

/// Smallest distance between two vertices, between a vertex and an edge
/// not incident to it, between two edges sharing no vertex, between edges
/// sharing a vertex at a near-contact away from it, and the width of every
/// two-edge sliver (edges as fine polylines).
fn feature_size(a: &SketchArrangement, curves: &[ProfileCurve]) -> f64 {
    let polys: Vec<Vec<DVec2>> = a
        .edges
        .iter()
        .map(|e| {
            (0..=256)
                .map(|i| curve_point(&curves[e.curve], e.t0 + (e.t1 - e.t0) * i as f64 / 256.0))
                .collect()
        })
        .collect();
    let to_poly = |p: DVec2, poly: &[DVec2]| {
        poly.windows(2)
            .map(|w| {
                let ab = w[1] - w[0];
                let s = ((p - w[0]).dot(ab) / ab.length_squared().max(1e-300)).clamp(0.0, 1.0);
                (w[0] + ab * s - p).length()
            })
            .fold(f64::INFINITY, f64::min)
    };
    let mut d = f64::INFINITY;
    for i in 0..a.vertices.len() {
        for j in i + 1..a.vertices.len() {
            d = d.min((a.vertices[i].point - a.vertices[j].point).length());
        }
        for (k, e) in a.edges.iter().enumerate() {
            if e.start != i && e.end != i {
                d = d.min(to_poly(a.vertices[i].point, &polys[k]));
            }
        }
    }
    for i in 0..a.edges.len() {
        for j in 0..a.edges.len() {
            let (e, f) = (&a.edges[i], &a.edges[j]);
            let shared: Vec<DVec2> = [e.start, e.end]
                .into_iter()
                .filter(|x| *x == f.start || *x == f.end)
                .map(|x| a.vertices[x].point)
                .collect();
            if i == j {
                continue;
            }
            if shared.is_empty() {
                if i < j {
                    for p in &polys[i] {
                        d = d.min(to_poly(*p, &polys[j]));
                    }
                }
                continue;
            }
            // Edges meeting at a vertex: the gap to the other edge, sampled
            // along this one, grows away from the shared vertices. A local
            // minimum away from them is a near-contact elsewhere, and two
            // edges joining the same two vertices bound a sliver whose width
            // is the largest gap.
            let gaps: Vec<f64> = polys[i].iter().map(|p| to_poly(*p, &polys[j])).collect();
            for k in 1..gaps.len() - 1 {
                let far = shared.iter().all(|v| (polys[i][k] - *v).length() > 0.02);
                if far && gaps[k] < gaps[k - 1] && gaps[k] <= gaps[k + 1] {
                    d = d.min(gaps[k]);
                }
            }
            let both = e.start != e.end && shared.len() == 2;
            if both {
                d = d.min(gaps.iter().copied().fold(0.0, f64::max));
            }
        }
    }
    d
}

/// A line, circle or arc on the integer lattice: lines run along the axes
/// or diagonals between lattice points, circles and arcs have lattice
/// centres and integer radii, and arcs start and end at quarter turns.
#[derive(Clone, Debug)]
enum Lat {
    Line(i64, i64, usize, i64),
    Circle(i64, i64, i64),
    Arc(i64, i64, i64, i64, i64),
}

const DIRS: [(i64, i64); 8] = [
    (1, 0),
    (1, 1),
    (0, 1),
    (-1, 1),
    (-1, 0),
    (-1, -1),
    (0, -1),
    (1, -1),
];

fn lat() -> impl Strategy<Value = Lat> {
    prop_oneof![
        3 => (0i64..5, 0i64..5, 0usize..8, 1i64..5).prop_map(|(x, y, d, l)| Lat::Line(x, y, d, l)),
        2 => (0i64..5, 0i64..5, 1i64..4).prop_map(|(x, y, r)| Lat::Circle(x, y, r)),
        2 => (0i64..5, 0i64..5, 1i64..4, 0i64..4, prop_oneof![-4i64..=-1, 1i64..=4])
            .prop_map(|(x, y, r, s, w)| Lat::Arc(x, y, r, s, w)),
    ]
}

/// Perturbations in multiples of the tolerance.
fn jitter() -> impl Strategy<Value = f64> {
    prop_oneof![
        6 => Just(0.0),
        2 => prop::sample::select(vec![-0.3, -0.1, 0.1, 0.3]),
        1 => prop::sample::select(vec![-3.0, -1.1, -0.9, 0.9, 1.1, 3.0]),
    ]
}

/// The curve of `l` at scale `s`, its defining points moved by `j` tolerances.
fn lat_curve(l: &Lat, s: f64, j: [f64; 4]) -> ProfileCurve {
    let t = tol();
    let p = |x: i64, y: i64| v(x as f64, y as f64) * s;
    match *l {
        Lat::Line(x, y, d, len) => {
            let (dx, dy) = DIRS[d];
            line(
                p(x, y) + v(j[0], j[1]) * t,
                p(x + dx * len, y + dy * len) + v(j[2], j[3]) * t,
            )
        }
        Lat::Circle(x, y, r) => circle(p(x, y) + v(j[0], j[1]) * t, r as f64 * s + j[2] * t),
        Lat::Arc(x, y, r, a, w) => {
            let rr = r as f64 * s + j[2] * t;
            arc(
                p(x, y) + v(j[0], j[1]) * t,
                rr,
                a as f64 * FRAC_PI_2 + j[3] * t / rr,
                w as f64 * FRAC_PI_2,
            )
        }
    }
}

/// Lattice inputs rejected as unresolvable at the tolerance (allowed only
/// for perturbations near or above the tolerance).
static LATTICE_REJECTED: AtomicUsize = AtomicUsize::new(0);

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    /// Lines, circles and arcs on an integer lattice: exact tangencies,
    /// overlaps, shared end points and curves through crossings everywhere,
    /// optionally perturbed by fractions or multiples of the tolerance, at
    /// scales 1e-3, 1 and 1e3. Checked against a raster of the exact curves
    /// (no two lattice curves come closer than 0.12 without touching):
    /// * every raster component lies in one region or none, and the
    ///   unbounded component in none;
    /// * unperturbed, or perturbed by at most 0.3 tolerances: the
    ///   arrangement must succeed (faces thinner than the tolerance that
    ///   such small moves open up are absorbed), every large enclosed
    ///   component lies in a region, and no region covers two separate
    ///   components;
    /// * perturbed by about the tolerance or more: an `InvalidInput` error
    ///   is allowed, a wrong region is not.
    #[test]
    fn lattice_lines_circles_and_arcs(
        items in prop::collection::vec((lat(), [jitter(), jitter(), jitter(), jitter()]), 1..8),
        perturb in prop::bool::weighted(0.5),
        scale in prop::sample::select(vec![1.0, 1e-3, 1e3]),
    ) {
        check_lattice(&items, perturb, scale)?;
    }
}

/// The checks of [`lattice_lines_circles_and_arcs`] for one input.
fn check_lattice(
    items: &[(Lat, [f64; 4])],
    perturb: bool,
    scale: f64,
) -> Result<(), TestCaseError> {
    let exact: Vec<ProfileCurve> = items
        .iter()
        .map(|(l, _)| lat_curve(l, 1.0, [0.0; 4]))
        .collect();
    let curves: Vec<ProfileCurve> = items
        .iter()
        .map(|(l, j)| lat_curve(l, scale, if perturb { *j } else { [0.0; 4] }))
        .collect();
    let big = perturb && items.iter().any(|(_, j)| j.iter().any(|x| x.abs() > 0.5));
    let a = match k().sketch_regions(&curves) {
        Ok(a) => a,
        Err(KernelError::InvalidInput(m)) if big => {
            let n = LATTICE_REJECTED.fetch_add(1, Ordering::Relaxed) + 1;
            eprintln!("LATTICE_REJECTED #{n}: {m}");
            return Ok(());
        }
        Err(e) => return Err(TestCaseError::fail(format!("{e}"))),
    };
    check_invariants(&curves, &a);
    let raster = Raster::new(&exact, v(-5.0, -5.0), v(9.0, 9.0), 1.0 / 64.0);
    let mut owner: HashMap<usize, usize> = HashMap::new();
    for (id, pts) in raster.samples(6).iter().enumerate() {
        if pts.is_empty() {
            continue;
        }
        let got: Vec<Option<usize>> = pts
            .iter()
            .map(|p| region_containing(&a, &curves, *p * scale, &Tolerance::default()))
            .collect();
        prop_assert!(
            got.iter().all(|g| *g == got[0]),
            "one raster component in regions {:?} at {:?}",
            got,
            pts
        );
        if raster.outer[id] {
            prop_assert_eq!(got[0], None, "unbounded area at {:?} in a region", pts[0]);
            continue;
        }
        if big || raster.sizes[id] < 256 {
            continue;
        }
        let Some(r) = got[0] else {
            return Err(TestCaseError::fail(format!(
                "enclosed area at {:?} is in no region",
                pts[0]
            )));
        };
        if let Some(other) = owner.insert(r, id) {
            return Err(TestCaseError::fail(format!(
                "region {r} covers separate areas at {:?} and {:?}",
                pts[0],
                raster.samples(1)[other]
            )));
        }
    }
    let jmax = items
        .iter()
        .flat_map(|(_, j)| j.iter())
        .fold(0.0f64, |m, x| m.max(x.abs()));
    extrude_all_within(&curves, &a, if perturb { 3.0 * jmax * tol() } else { 0.0 });
    Ok(())
}

/// Two circles and the line through both their crossings, each moved by
/// 0.3 tolerances (found by the lattice test): the three crossings land
/// 1.3 to 1.4 tolerances apart and bound a triangle 0.4 tolerances across,
/// once refused as a degenerate face. It merges into one vertex.
#[test]
fn lattice_near_concurrent_crossings() {
    let items = [
        (Lat::Circle(3, 2, 3), [0.0, 0.3, 0.0, 0.0]),
        (Lat::Arc(1, 0, 3, 0, -2), [0.3, -0.3, 0.3, 0.0]),
        (Lat::Line(2, 1, 7, 2), [0.0, 0.0, 0.0, -0.3]),
    ];
    for scale in [1.0, 1e-3, 1e3] {
        check_lattice(&items, true, scale).unwrap();
    }
}

/// Four arcs tangent at one point, one of them twice (a copy moved by 0.3
/// tolerances), at millimetre scale (found by the lattice test). Nearly
/// tangent curves stay within the tolerance of each other for microns, and
/// contacts chained along that stretch split off a sliver between two of
/// them 6 µm long and up to 1.1 tolerances wide, once refused as a region
/// thinner than the tolerance. It is absorbed as an overlap.
#[test]
fn lattice_tangent_arcs_leave_no_sliver() {
    let items = [
        (Lat::Arc(1, 2, 2, 0, 2), [0.0, -0.1, 0.0, 0.0]),
        (Lat::Arc(3, 3, 1, 2, 2), [-0.1, 0.0, 0.0, 0.0]),
        (Lat::Arc(0, 3, 2, 3, 1), [0.0; 4]),
        (Lat::Arc(4, 3, 2, 2, 1), [0.0; 4]),
        (Lat::Line(3, 1, 3, 2), [0.0, 0.1, 0.0, 0.0]),
        (Lat::Arc(4, 3, 2, 2, -4), [0.0, -0.3, -0.3, 0.0]),
    ];
    for scale in [1.0, 1e-3, 1e3] {
        check_lattice(&items, true, scale).unwrap();
    }
}

/// Three circles tangent at one point (two inside each other, one outside)
/// and two lines starting there, moved by up to three tolerances (found by
/// the lattice test at scale 10³ and, shrinking, at scale 1). The contacts
/// make two vertices three tolerances apart joined by pieces of the
/// tangent curves that coincide to rounding: their area's sign was noise,
/// so one read as a second outer boundary, and pieces seen at their
/// middles inside the vertices were ordered by noise (Euler's formula
/// failed). Points inside a vertex now count as at it when overlaps are
/// found, and faces with an area of rounding noise are absorbed before
/// the checks: the first sketch gets its four regions at every scale; the
/// second, whose line starts three tolerances above the tangency, leaves
/// a triangle 0.9 tolerances thick there and is refused, naming where.
#[test]
fn lattice_tangent_circles_with_lines_at_the_tangency() {
    let cases = [
        [
            (Lat::Arc(0, 1, 1, 2, -2), [0.0; 4]),
            (Lat::Line(0, 2, 7, 1), [0.0, 3.0, 0.0, 0.0]),
            (Lat::Circle(0, 0, 2), [-0.3, 0.0, 0.0, 0.0]),
            (Lat::Arc(0, 3, 1, 3, -4), [0.0, 0.0, 0.0, 1.1]),
            (Lat::Line(0, 2, 1, 1), [-0.3, 0.0, 0.0, 0.0]),
        ],
        [
            (Lat::Arc(0, 1, 1, 2, -2), [0.0; 4]),
            (Lat::Line(0, 2, 7, 1), [0.0, 3.0, 0.0, 0.0]),
            (Lat::Circle(0, 0, 2), [-0.3, 0.0, 0.0, 0.0]),
            (Lat::Arc(0, 3, 1, 3, 1), [0.0, 0.0, 0.0, 1.1]),
            (Lat::Line(0, 4, 6, 2), [0.0, 0.0, -0.3, -0.3]),
        ],
    ];
    for items in &cases {
        for scale in [1.0, 1e-3, 1e3] {
            check_lattice(items, true, scale).unwrap();
        }
    }
}

/// A chord of a circle on a longer line that diverges from it by about a
/// tolerance, with a line tangent to the circle at its end.
#[test]
fn lattice_chord_along_a_nearly_collinear_line() {
    let cases = [[
        (Lat::Circle(1, 1, 1), [0.0; 4]),
        (Lat::Line(1, 2, 7, 1), [0.0, 0.0, 0.0, 1.1]),
        (Lat::Line(3, 0, 3, 3), [0.0, 0.0, 3.0, 0.0]),
        (Lat::Line(0, 0, 0, 1), [0.0; 4]),
    ]];
    for items in &cases {
        for scale in [1.0, 1e-3, 1e3] {
            check_lattice(items, true, scale).unwrap_or_else(|e| panic!("scale {scale}: {e}"));
        }
    }
}

/// A closed star-shaped spline about `c` through the star polygon's
/// corners, optionally with random weights (then it passes through the first
/// and last points only).
fn star_spline(c: DVec2, pts: &[(f64, f64)], weights: Option<Vec<f64>>) -> ProfileCurve {
    let mut p = star(c, pts);
    p.push(p[0]);
    let ProfileCurve::Spline {
        degree,
        knots,
        points,
        ..
    } = ProfileCurve::interpolate(&p).unwrap()
    else {
        unreachable!()
    };
    let weights = weights.map(|w| (0..points.len()).map(|i| w[i % w.len()]).collect());
    ProfileCurve::Spline {
        degree,
        knots,
        points,
        weights,
    }
}

/// True if the closed polygon has no two crossing non-adjacent sides.
fn simple(p: &[DVec2]) -> bool {
    let n = p.len();
    let cross = |a: DVec2, b: DVec2, c: DVec2, d: DVec2| {
        let o = |p: DVec2, q: DVec2, r: DVec2| (q - p).perp_dot(r - p);
        (o(a, b, c) > 0.0) != (o(a, b, d) > 0.0) && (o(c, d, a) > 0.0) != (o(c, d, b) > 0.0)
    };
    for i in 0..n {
        for j in i + 2..n {
            if i == 0 && j == n - 1 {
                continue;
            }
            if cross(p[i], p[(i + 1) % n], p[j], p[(j + 1) % n]) {
                return false;
            }
        }
    }
    true
}

#[derive(Clone, Debug)]
enum SplineShape {
    /// A star spline, with an optional overlapping copy: 0 none, 1 the
    /// same, 2 reversed, 3 a third of it.
    Star(DVec2, Vec<(f64, f64)>, Option<Vec<f64>>, u8),
    Circle(DVec2, f64),
}

fn spline_shape() -> impl Strategy<Value = SplineShape> {
    prop_oneof![
        (vec2(-1.5, 1.5), 0.2f64..1.2).prop_map(|(c, r)| SplineShape::Circle(c, r)),
        (
            vec2(-1.5, 1.5),
            prop::collection::vec((0.0f64..0.6, 0.5f64..1.3), 4..9),
            prop::option::weighted(0.4, prop::collection::vec(0.5f64..2.0, 1..6)),
            prop_oneof![5 => Just(0u8), 1 => 1u8..4],
        )
            .prop_map(|(c, p, w, d)| SplineShape::Star(c, p, w, d)),
    ]
}

/// An overlapping copy of a spline: itself (`kind` 1), reversed (2), or
/// its piece from 40% to 75% of its domain, split off by knot insertion (3).
fn spline_copy(sp: &ProfileCurve, kind: u8) -> ProfileCurve {
    use openshape_nurbs::NurbsCurve;
    let ProfileCurve::Spline {
        degree,
        knots,
        points,
        weights,
    } = sp
    else {
        unreachable!()
    };
    let n = NurbsCurve::new(
        *degree,
        knots.clone(),
        points.iter().map(|p| p.extend(0.0)).collect(),
        weights.clone(),
    )
    .unwrap();
    let n = match kind {
        2 => n.reversed(),
        3 => {
            let (lo, hi) = n.domain();
            let rest = n.try_split(lo + 0.4 * (hi - lo)).unwrap().1;
            rest.try_split(lo + 0.75 * (hi - lo)).unwrap().0
        }
        _ => n,
    };
    ProfileCurve::Spline {
        degree: n.degree(),
        knots: n.knots().as_slice().to_vec(),
        points: n.control_points().iter().map(|p| p.truncate()).collect(),
        weights: if n.is_rational() {
            Some(n.weights())
        } else {
            None
        },
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    /// Closed polynomial and rational splines with circles, some drawn
    /// twice (the same, reversed, or in part): the regions tile the
    /// enclosed area (against a union of fine polygons), every point is in
    /// the right region, and every region extrudes.
    #[test]
    fn splines_and_circles(
        shapes in prop::collection::vec(spline_shape(), 1..5),
        samples in prop::collection::vec(vec2(-3.0, 3.0), 48),
    ) {
        let mut curves = Vec::new();
        let mut polys = Vec::new();
        for s in &shapes {
            match s {
                SplineShape::Circle(c, r) => {
                    curves.push(circle(*c, *r));
                    polys.push(Shape::Circle(*c, *r));
                }
                SplineShape::Star(c, p, w, d) => {
                    let sp = star_spline(*c, p, w.clone());
                    let mut poly = spline_samples(&sp, 2048);
                    poly.pop();
                    prop_assume!(simple(&poly) && shoelace(&poly) > 0.0);
                    if *d > 0 {
                        curves.push(spline_copy(&sp, *d));
                    }
                    curves.push(sp);
                    polys.push(Shape::Poly(poly));
                }
            }
        }
        let oracle = Shapes::new(polys);
        let a = arrange_checked(&curves);
        let sum: f64 = a.regions.iter().map(|r| r.area).sum();
        let want = oracle.filled_area();
        prop_assert!((sum - want).abs() <= 1e-4 * want, "region areas {sum} vs enclosed area {want}");
        check_partition(&a, &curves, &samples, 5e-3, |p| oracle.filled_contains(p));
        extrude_every_region(&curves, &a);
    }
}

/// Values that break naive code: non-finite, huge, tiny, and lattice values.
fn wild() -> impl Strategy<Value = f64> {
    prop_oneof![
        4 => (-4i32..5).prop_map(f64::from),
        3 => -3.0f64..3.0,
        2 => prop::sample::select(vec![
            0.0, -0.0, 1e-300, 5e-8, 1e-7, 2e-7, 1e-3, 0.5, 1e6, 1e12, 1e300, -1e300,
            f64::MIN_POSITIVE, f64::MAX, f64::MIN, f64::EPSILON,
            f64::NAN, f64::INFINITY, f64::NEG_INFINITY,
        ]),
        1 => any::<f64>(),
    ]
}

fn wild_curve() -> impl Strategy<Value = ProfileCurve> {
    let p = || (wild(), wild()).prop_map(|(x, y)| v(x, y));
    prop_oneof![
        (p(), p()).prop_map(|(a, b)| line(a, b)),
        (p(), wild()).prop_map(|(c, r)| circle(c, r)),
        (p(), wild(), wild(), wild()).prop_map(|(c, r, s, w)| arc(c, r, s, w)),
        (
            prop_oneof![4 => 0usize..5, 1 => 0usize..40],
            prop::collection::vec(p(), 0..9),
            prop_oneof![
                prop::collection::vec(wild(), 0..16),
                prop::collection::vec(-3i32..4, 0..16).prop_map(|mut k| {
                    k.sort_unstable();
                    k.into_iter().map(f64::from).collect()
                }),
            ],
            prop::option::of(prop::collection::vec(wild(), 0..9)),
        )
            .prop_map(|(degree, points, knots, weights)| ProfileCurve::Spline {
                degree,
                knots,
                points,
                weights
            }),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// Garbage and borderline curves never panic; whatever succeeds is
    /// self-consistent and its regions' profiles and points work (a profile
    /// may be refused only for geometry below the profile builder's
    /// resolution of ten tolerances).
    #[test]
    fn garbage_never_panics(curves in prop::collection::vec(wild_curve(), 0..6), p in (wild(), wild())) {
        let a = match k().sketch_regions(&curves) {
            Ok(a) => a,
            Err(KernelError::InvalidInput(_)) => return Ok(()),
            Err(e) => return Err(TestCaseError::fail(format!("{e}"))),
        };
        let profile_ok = |sel: &[usize]| match k().region_profile(&a, &curves, sel) {
            Ok(prof) => !prof.is_empty(),
            Err(KernelError::InvalidInput(m)) => m.contains("profile resolution"),
            Err(_) => false,
        };
        for (r, reg) in a.regions.iter().enumerate() {
            prop_assert!(reg.area.is_finite() && reg.area > 0.0);
            prop_assert_eq!(region_containing(&a, &curves, reg.interior_point, &Tolerance::default()), Some(r));
            prop_assert!(profile_ok(&[r]), "region {}", r);
        }
        let _ = region_containing(&a, &curves, v(p.0, p.1), &Tolerance::default());
        if !a.regions.is_empty() {
            let all: Vec<usize> = (0..a.regions.len()).collect();
            prop_assert!(profile_ok(&all));
        }
    }
}
