//! Persistent names (design note 0013) of boxes, extrudes and revolves
//! made by the kernel, with sketch keys made up for their profiles: right
//! by an independent geometric check, the same across dimension changes,
//! revolve angles and curve order, and resolved as §6 of the note says.

use std::collections::BTreeSet;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

use opendatum_core_model::naming::{
    Along, EdgeName, FaceName, Hint, Name, ProfilePoint, Reference, Role, Source, Use,
};
use opendatum_core_model::sketch::{EdgeKey, SketchId, VertexKey};
use opendatum_core_model::{FeatureId, RegionKey};
use opendatum_kernel::{
    Adjacency, Body, BodyMesh, BrepKernel, DVec2, DVec3, EdgeIx, ExtrudeSpec, FaceIx, Frame,
    Kernel, Profile, ProfileCurve, RevolveSpec, TessellationOptions, VertexIx,
};
use opendatum_naming::{
    BodyNames, BodySource, Made, NamingError, Pick, PickedRegion, ResolveError, Resolved,
};
use opendatum_oplog::{OpId, ReplicaId};

fn k() -> BrepKernel {
    BrepKernel::default()
}

/// Mesh distances below this count as "on" in the geometric checks: the
/// tessellation's chord tolerance with room to spare, far below the
/// distances between the test models' entities.
const CLOSE: f64 = 2e-3;

fn opts() -> TessellationOptions {
    TessellationOptions {
        chord_tolerance: 1e-3,
        max_angle: 0.3,
    }
}

// ---------------------------------------------------------------------------
// Made-up sketches.

fn sid(i: u32) -> SketchId {
    SketchId {
        op: OpId {
            lamport: 7,
            replica: ReplicaId(3),
        },
        index: i,
    }
}

fn feature(n: u64) -> FeatureId {
    FeatureId(OpId {
        lamport: n,
        replica: ReplicaId(1),
    })
}

/// The sketch vertex where the sketch curves `curves` meet.
fn v(curves: &[u32]) -> VertexKey {
    let mut entities: Vec<SketchId> = curves.iter().map(|&c| sid(c)).collect();
    entities.sort();
    VertexKey {
        entities,
        ordinal: 0,
        parts: Vec::new(),
    }
}

fn edge(curve: u32, forward: bool, from: VertexKey, to: VertexKey) -> EdgeKey {
    EdgeKey {
        entity: sid(curve),
        forward,
        from,
        to,
        part: 0,
    }
}

/// A profile as the features crate makes one from sketch regions: curves
/// drawn along the boundary of the region they bound, the sketch edge of
/// each, and the keys of the regions selected.
#[derive(Clone, Debug)]
struct Sketch {
    plane: Frame,
    curves: Vec<ProfileCurve>,
    keys: Vec<EdgeKey>,
    regions: Vec<RegionKey>,
}

impl Sketch {
    fn new(plane: Frame) -> Sketch {
        Sketch {
            plane,
            curves: Vec::new(),
            keys: Vec::new(),
            regions: Vec::new(),
        }
    }

    fn curve(&mut self, c: ProfileCurve, key: &EdgeKey) {
        self.curves.push(c);
        self.keys.push(key.clone());
    }

    fn region(&mut self, loops: &[&[&EdgeKey]]) {
        self.regions.push(region_key(loops));
    }

    fn profile(&self) -> Profile {
        Profile {
            plane: self.plane,
            curves: self.curves.clone(),
        }
    }

    /// The curves with sketch edge `key`.
    fn curves_of<'a>(&'a self, key: &'a EdgeKey) -> impl Iterator<Item = &'a ProfileCurve> + 'a {
        self.curves
            .iter()
            .zip(&self.keys)
            .filter(move |(_, k)| *k == key)
            .map(|(c, _)| c)
    }

    /// Moves curve `i` (and its key) to the end, keeping the rest in order.
    fn rotated(&self, i: usize) -> Sketch {
        let mut s = self.clone();
        let c = s.curves.remove(i);
        let key = s.keys.remove(i);
        s.curves.push(c);
        s.keys.push(key);
        s
    }
}

fn region_key(loops: &[&[&EdgeKey]]) -> RegionKey {
    RegionKey {
        loops: loops
            .iter()
            .map(|l| l.iter().map(|&k| k.clone()).collect())
            .collect(),
        interior: DVec2::ZERO,
    }
}

fn line(a: DVec2, b: DVec2) -> ProfileCurve {
    ProfileCurve::Line { start: a, end: b }
}

fn arc(center: DVec2, radius: f64, start_angle: f64, sweep_angle: f64) -> ProfileCurve {
    ProfileCurve::Arc {
        center,
        radius,
        start_angle,
        sweep_angle,
    }
}

/// The keys of the slot-and-square sketch (see [`slot_and_square`]).
struct SlotKeys {
    bottom: EdgeKey,
    right: EdgeKey,
    top: EdgeKey,
    left: EdgeKey,
    hole: EdgeKey,
    square: [EdgeKey; 4],
}

fn slot_keys() -> SlotKeys {
    let (a, b, c, d) = (v(&[0, 3]), v(&[0, 1]), v(&[1, 2]), v(&[2, 3]));
    let s = [v(&[5, 8]), v(&[5, 6]), v(&[6, 7]), v(&[7, 8])];
    SlotKeys {
        bottom: edge(0, true, a.clone(), b.clone()),
        right: edge(1, true, b, c.clone()),
        // Drawn from left to right; the boundary runs against it.
        top: edge(2, false, c, d.clone()),
        left: edge(3, true, d, a),
        // A circle: one sketch edge from its one vertex round to it,
        // clockwise as the boundary of a hole.
        hole: edge(4, false, v(&[4]), v(&[4])),
        square: [0, 1, 2, 3].map(|i| {
            edge(
                5 + i,
                true,
                s[i as usize].clone(),
                s[(i as usize + 1) % 4].clone(),
            )
        }),
    }
}

/// A slot of length `l` and radius `r` (two lines and two half circles)
/// with a round hole, and a separate square: two regions, a hole, arcs of
/// exactly π and a closed sketch edge cut in two, as the region profile
/// cuts it. Sketch curves 0 to 3 are the slot's bottom, right, top and
/// left, 4 the hole, 5 to 8 the square's sides.
fn slot_and_square(plane: Frame, l: f64, r: f64, regions: [bool; 2]) -> Sketch {
    let ks = slot_keys();
    let mut s = Sketch::new(plane);
    if regions[0] {
        s.curve(line(DVec2::new(0.0, -r), DVec2::new(l, -r)), &ks.bottom);
        s.curve(arc(DVec2::new(l, 0.0), r, -FRAC_PI_2, PI), &ks.right);
        s.curve(line(DVec2::new(l, r), DVec2::new(0.0, r)), &ks.top);
        s.curve(arc(DVec2::ZERO, r, FRAC_PI_2, PI), &ks.left);
        let c = DVec2::new(0.5 * l, 0.0);
        s.curve(arc(c, 0.5 * r, 0.0, -PI), &ks.hole);
        s.curve(arc(c, 0.5 * r, -PI, -PI), &ks.hole);
        s.region(&[&[&ks.bottom, &ks.right, &ks.top, &ks.left], &[&ks.hole]]);
    }
    if regions[1] {
        let x0 = 8.0;
        let p = [
            DVec2::new(x0, -1.0),
            DVec2::new(x0 + 2.0, -1.0),
            DVec2::new(x0 + 2.0, 1.0),
            DVec2::new(x0, 1.0),
        ];
        for i in 0..4 {
            s.curve(line(p[i], p[(i + 1) % 4]), &ks.square[i]);
        }
        s.region(&[&[&ks.square[0], &ks.square[1], &ks.square[2], &ks.square[3]]]);
    }
    s
}

fn tilted() -> Frame {
    let z = DVec3::new(0.3, -0.4, 0.866).normalize();
    let x = DVec3::Y.cross(z).normalize();
    Frame {
        origin: DVec3::new(0.5, -1.0, 2.0),
        x,
        y: z.cross(x),
        z,
    }
}

// ---------------------------------------------------------------------------
// Bodies and their names.

#[derive(Clone, Copy, Debug)]
enum Sweep {
    Extrude(ExtrudeSpec),
    Revolve(RevolveSpec),
}

fn extrude(distance: f64, symmetric: bool) -> Sweep {
    Sweep::Extrude(ExtrudeSpec {
        distance,
        symmetric,
    })
}

fn revolve(axis_origin: DVec2, axis_dir: DVec2, angle: f64) -> Sweep {
    Sweep::Revolve(RevolveSpec {
        axis_origin,
        axis_dir,
        angle,
    })
}

/// A body with its names, mesh and adjacency.
struct Named {
    body: Body,
    names: BodyNames,
    mesh: BodyMesh,
    adj: Adjacency,
}

fn make(s: &Sketch, sweep: Sweep) -> Body {
    let p = s.profile();
    match sweep {
        Sweep::Extrude(e) => k().extrude(&p, &e),
        Sweep::Revolve(r) => k().revolve(&p, &r),
    }
    .unwrap_or_else(|e| panic!("{sweep:?}: {e}"))
}

fn name_body(
    s: &Sketch,
    body: &Body,
    f: FeatureId,
    stored: &[RegionKey],
) -> Result<BodyNames, NamingError> {
    let picked: Vec<PickedRegion> = stored
        .iter()
        .zip(&s.regions)
        .map(|(stored, current)| PickedRegion { stored, current })
        .collect();
    BodyNames::new(
        k().origins(body),
        &k().adjacency(body),
        &BodySource {
            feature: f,
            made: Made::Sweep {
                curves: &s.keys,
                regions: &picked,
            },
        },
    )
}

/// Makes, names and checks a sweep with the regions' current keys stored.
fn named(s: &Sketch, sweep: Sweep) -> Named {
    named_as(s, sweep, &s.regions)
}

fn named_as(s: &Sketch, sweep: Sweep, stored: &[RegionKey]) -> Named {
    let body = make(s, sweep);
    let names = name_body(s, &body, feature(1), stored).unwrap();
    let n = Named {
        mesh: k().tessellate(&body, &opts()).unwrap(),
        adj: k().adjacency(&body),
        body,
        names,
    };
    check(s, sweep, stored, &n);
    n
}

// ---------------------------------------------------------------------------
// The geometric oracle.

fn tri_distance(p: DVec3, a: DVec3, b: DVec3, c: DVec3) -> f64 {
    // Closest point on a triangle (Ericson, Real-Time Collision Detection
    // 5.1.5), by its barycentric regions.
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let (d1, d2) = (ab.dot(ap), ac.dot(ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return (p - a).length();
    }
    let bp = p - b;
    let (d3, d4) = (ab.dot(bp), ac.dot(bp));
    if d3 >= 0.0 && d4 <= d3 {
        return (p - b).length();
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return (p - (a + ab * (d1 / (d1 - d3)))).length();
    }
    let cp = p - c;
    let (d5, d6) = (ab.dot(cp), ac.dot(cp));
    if d6 >= 0.0 && d5 <= d6 {
        return (p - c).length();
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return (p - (a + ac * (d2 / (d2 - d6)))).length();
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (p - (b + (c - b) * w)).length();
    }
    let denom = 1.0 / (va + vb + vc);
    (p - (a + ab * (vb * denom) + ac * (vc * denom))).length()
}

fn segment_distance(p: DVec3, a: DVec3, b: DVec3) -> f64 {
    let ab = b - a;
    let t = ((p - a).dot(ab) / ab.length_squared()).clamp(0.0, 1.0);
    (p - (a + ab * t)).length()
}

fn face_distance(mesh: &BodyMesh, f: FaceIx, p: DVec3) -> f64 {
    let m = &mesh.faces[f.0 as usize];
    m.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| m.positions[i as usize]);
            tri_distance(p, a, b, c)
        })
        .fold(f64::INFINITY, f64::min)
}

fn edge_distance(mesh: &BodyMesh, e: EdgeIx, p: DVec3) -> f64 {
    mesh.edges[e.0 as usize]
        .windows(2)
        .map(|w| segment_distance(p, w[0], w[1]))
        .fold(f64::INFINITY, f64::min)
}

/// Where a vertex is.
fn vertex_point(mesh: &BodyMesh, v: VertexIx) -> DVec3 {
    mesh.vertices[v.0 as usize]
}

/// A point inside an edge, away from its ends.
fn edge_point(mesh: &BodyMesh, e: EdgeIx) -> DVec3 {
    let l = &mesh.edges[e.0 as usize];
    if l.len() > 2 {
        l[l.len() / 2]
    } else {
        (l[0] + l[1]) * 0.5
    }
}

/// Pieces an arc of this sweep is split into (at most π each).
fn pieces(sweep: f64) -> usize {
    ((sweep.abs() / PI) - 1e-9).ceil().max(1.0) as usize
}

/// The point at fraction `f` along piece `piece` of a curve, in the
/// curve's own direction.
fn piece_point(c: &ProfileCurve, piece: u32, f: f64) -> DVec2 {
    let on_arc = |center: DVec2, radius: f64, start: f64, sweep: f64| {
        let a = start + sweep * (piece as f64 + f) / pieces(sweep) as f64;
        center + DVec2::from_angle(a) * radius
    };
    match *c {
        ProfileCurve::Line { start, end } => start + (end - start) * f,
        ProfileCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => on_arc(center, radius, start_angle, sweep_angle),
        ProfileCurve::Circle { center, radius } => on_arc(center, radius, 0.0, TAU),
        ProfileCurve::Spline { .. } => unreachable!("test profiles have no splines"),
    }
}

fn curve_distance(c: &ProfileCurve, q: DVec2) -> f64 {
    match *c {
        ProfileCurve::Line { start, end } => {
            let d = end - start;
            let t = ((q - start).dot(d) / d.length_squared()).clamp(0.0, 1.0);
            (q - (start + d * t)).length()
        }
        ProfileCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => {
            let rel = q - center;
            // The angle from the arc's start in its direction, in [0, 2π).
            let a = (rel.y.atan2(rel.x) - start_angle) * sweep_angle.signum();
            if a.rem_euclid(TAU) <= sweep_angle.abs() {
                (rel.length() - radius).abs()
            } else {
                let p0 = piece_point(c, 0, 0.0);
                let p1 = piece_point(c, pieces(sweep_angle) as u32 - 1, 1.0);
                (q - p0).length().min((q - p1).length())
            }
        }
        ProfileCurve::Circle { center, radius } => ((q - center).length() - radius).abs(),
        ProfileCurve::Spline { .. } => unreachable!("test profiles have no splines"),
    }
}

fn wrap(a: f64) -> f64 {
    let a = a.rem_euclid(TAU);
    if a > PI {
        a - TAU
    } else {
        a
    }
}

/// A point of a swept body taken back to the profile: its profile point
/// and how far along the sweep it is (the height for an extrude, the
/// angle for a revolve, NaN on a revolve's axis).
fn unsweep(s: &Sketch, sweep: Sweep, x: DVec3) -> (DVec2, f64) {
    let plane = &s.plane;
    match sweep {
        Sweep::Extrude(_) => {
            let rel = x - plane.origin;
            (
                DVec2::new(rel.dot(plane.x), rel.dot(plane.y)),
                rel.dot(plane.z),
            )
        }
        Sweep::Revolve(r) => {
            let d2 = r.axis_dir.normalize();
            let p0 = piece_point(&s.curves[0], 0, 0.5) - r.axis_origin;
            let side = (d2.x * p0.y - d2.y * p0.x).signum();
            let n2 = DVec2::new(-d2.y, d2.x) * side;
            let (o, a, n3) = (
                plane.point2(r.axis_origin),
                plane.vector2(d2),
                plane.vector2(n2),
            );
            let rel = x - o;
            let h = rel.dot(a);
            let radial = rel - a * h;
            let rad = radial.length();
            let phi = if rad < 1e-9 {
                f64::NAN
            } else {
                radial.dot(a.cross(n3)).atan2(radial.dot(n3))
            };
            (r.axis_origin + d2 * h + n2 * rad, phi)
        }
    }
}

/// The sweep coordinate of its start (`end` false) or end.
fn sweep_end(sweep: Sweep, end: bool) -> f64 {
    match sweep {
        Sweep::Extrude(e) => {
            let base = if e.symmetric { -0.5 } else { 0.0 };
            e.distance * (base + if end { 1.0 } else { 0.0 })
        }
        Sweep::Revolve(r) => {
            if end {
                r.angle
            } else {
                0.0
            }
        }
    }
}

/// Whether a sweep coordinate is the start's or end's.
fn at_end(sweep: Sweep, w: f64, end: bool) -> bool {
    let target = sweep_end(sweep, end);
    match sweep {
        Sweep::Extrude(_) => (w - target).abs() < 1e-9,
        // On the axis every angle is.
        Sweep::Revolve(_) => w.is_nan() || wrap(w - target).abs() < 1e-9,
    }
}

/// The profile points a name's point may be (a closed sketch edge cut in
/// two has two curves with its key).
fn profile_points(s: &Sketch, p: &ProfilePoint) -> Vec<DVec2> {
    let f = if p.end { 1.0 } else { 0.0 };
    s.curves_of(&p.edge)
        .map(|c| piece_point(c, p.piece, f))
        .collect()
}

/// The 2D bounds of the curves of a region's outer loop.
fn loop_bounds(s: &Sketch, keys: &[EdgeKey]) -> (DVec2, DVec2) {
    let mut lo = DVec2::splat(f64::INFINITY);
    let mut hi = DVec2::splat(f64::NEG_INFINITY);
    for key in keys {
        for c in s.curves_of(key) {
            let n = match c {
                ProfileCurve::Arc { sweep_angle, .. } => pieces(*sweep_angle),
                ProfileCurve::Circle { .. } => 2,
                _ => 1,
            };
            for piece in 0..n as u32 {
                for i in 0..=64 {
                    let q = piece_point(c, piece, i as f64 / 64.0);
                    lo = lo.min(q);
                    hi = hi.max(q);
                }
            }
        }
    }
    (lo, hi)
}

/// The logical faces within [`CLOSE`] of a point.
fn faces_at(n: &Named, p: DVec3) -> BTreeSet<FaceName> {
    n.names
        .faces
        .iter()
        .filter(|f| {
            f.pieces
                .iter()
                .any(|&q| face_distance(&n.mesh, q, p) < CLOSE)
        })
        .map(|f| f.name.clone())
        .collect()
}

/// Every check of a swept body's names against its geometry.
fn check(s: &Sketch, sweep: Sweep, stored: &[RegionKey], n: &Named) {
    let names = &n.names;
    let o = k().origins(&n.body);
    // Every kernel entity is in at most one logical entity, every face in
    // one; the kernel's seams are exactly the edges in no logical edge.
    let mut seen = BTreeSet::new();
    for f in &names.faces {
        for p in &f.pieces {
            assert!(seen.insert(*p), "face {p:?} twice");
        }
    }
    assert_eq!(seen.len(), o.faces.len(), "every face named");
    let mut seen = BTreeSet::new();
    for e in &names.edges {
        for p in &e.pieces {
            assert!(seen.insert(*p), "edge {p:?} twice");
        }
    }
    for (e, &seam) in n.mesh.seams.iter().enumerate() {
        let e = EdgeIx(e as u32);
        assert_eq!(names.edge(e).is_none(), seam, "edge {e:?}: seam {seam}");
        assert_eq!(seen.contains(&e), !seam);
    }
    let all: BTreeSet<_> = names.faces.iter().map(|f| &f.name).collect();
    assert_eq!(all.len(), names.faces.len(), "face names unique");
    let all: BTreeSet<_> = names.edges.iter().map(|e| &e.name).collect();
    assert_eq!(all.len(), names.edges.len(), "edge names unique");
    let all: BTreeSet<_> = names.vertices.iter().map(|v| &v.name).collect();
    assert_eq!(all.len(), names.vertices.len(), "vertex names unique");

    // Faces: a side lies on its sketch edge swept, a cap at its end of the
    // sweep and within the regions it names.
    for f in &names.faces {
        assert_eq!(f.name.feature, feature(1));
        let points: Vec<DVec3> = f
            .pieces
            .iter()
            .flat_map(|p| n.mesh.faces[p.0 as usize].positions.iter().copied())
            .collect();
        assert!(!points.is_empty());
        match (&f.name.role, &f.name.sources[..]) {
            (Role::Side, sources) => {
                // Pieces of one surface look like one face, named by all
                // their sketch edges.
                let keys: Vec<&EdgeKey> = sources
                    .iter()
                    .map(|src| match src {
                        Source::SketchEdge(k) => k,
                        other => panic!("{:?}: a side made from {other:?}", f.name),
                    })
                    .collect();
                assert!(!keys.is_empty());
                for &x in &points {
                    let (q, _) = unsweep(s, sweep, x);
                    let d = keys
                        .iter()
                        .flat_map(|k| s.curves_of(k))
                        .map(|c| curve_distance(c, q))
                        .fold(f64::INFINITY, f64::min);
                    assert!(d < CLOSE, "{:?}: {x} is {d} off its curves", f.name);
                }
            }
            (Role::Cap { end }, sources) => {
                assert!(!sources.is_empty());
                let mut lo = DVec2::splat(f64::INFINITY);
                let mut hi = DVec2::splat(f64::NEG_INFINITY);
                for src in sources {
                    let Source::Region(outer) = src else {
                        panic!("{:?}: a cap source that is no region", f.name)
                    };
                    assert!(
                        stored.iter().any(|r| r.loops[0] == *outer),
                        "{:?}: a region not selected",
                        f.name
                    );
                    // The current keys of the stored region.
                    let i = stored.iter().position(|r| r.loops[0] == *outer).unwrap();
                    let (a, b) = loop_bounds(s, &s.regions[i].loops[0]);
                    lo = lo.min(a);
                    hi = hi.max(b);
                }
                for &x in &points {
                    let (q, w) = unsweep(s, sweep, x);
                    assert!(at_end(sweep, w, *end), "{:?}: {x} not at its end", f.name);
                    assert!(
                        q.cmpge(lo - CLOSE).all() && q.cmple(hi + CLOSE).all(),
                        "{:?}: {x} outside its regions",
                        f.name
                    );
                }
            }
            (role, sources) => panic!("unexpected face name {role:?} {sources:?}"),
        }
    }

    // Edges: inside, they touch exactly the faces they name; the
    // tie-break says where they are.
    for e in &names.edges {
        assert!(matches!(e.name.faces.len(), 1 | 2), "{:?}", e.name);
        for p in &e.pieces {
            let x = edge_point(&n.mesh, *p);
            assert_eq!(faces_at(n, x), e.name.faces, "edge {p:?} at {x}");
            let (q, w) = unsweep(s, sweep, x);
            match &e.name.along {
                None => {}
                Some(Along::Profile { edge, end }) => {
                    assert!(at_end(sweep, w, *end), "{:?}", e.name);
                    let d = s
                        .curves_of(edge)
                        .map(|c| curve_distance(c, q))
                        .fold(f64::INFINITY, f64::min);
                    assert!(d < CLOSE, "{:?}", e.name);
                }
                Some(Along::Swept(point)) => {
                    let ok = profile_points(s, point)
                        .iter()
                        .any(|&c| (c - q).length() < 1e-9);
                    assert!(ok, "{:?}: at {q}", e.name);
                }
                Some(other) => panic!("{:?}: unexpected tie-break {other:?}", e.name),
            }
        }
    }

    // Vertices: they touch exactly the faces they name; the tie-break
    // says where they are; no vertex is inside a logical edge.
    for vx in &names.vertices {
        assert_eq!(vx.pieces.len(), 1, "a vertex is one kernel vertex");
        let x = vertex_point(&n.mesh, vx.pieces[0]);
        assert_eq!(
            faces_at(n, x),
            vx.name.faces,
            "vertex {:?} at {x}",
            vx.pieces
        );
        match &vx.name.along {
            None => {}
            Some(Along::Point { point, end }) => {
                let (q, w) = unsweep(s, sweep, x);
                assert!(at_end(sweep, w, *end), "{:?}", vx.name);
                let ok = profile_points(s, point)
                    .iter()
                    .any(|&c| (c - q).length() < 1e-9);
                assert!(ok, "{:?}: at {q}", vx.name);
            }
            Some(other) => panic!("{:?}: unexpected tie-break {other:?}", vx.name),
        }
    }
    // The mesh's vertices are where their edges end.
    for (v, edges) in n.adj.vertex_edges.iter().enumerate() {
        let x = n.mesh.vertices[v];
        for e in edges {
            let l = &n.mesh.edges[e.0 as usize];
            let d = (l[0] - x).length().min((*l.last().unwrap() - x).length());
            assert!(d < 1e-9, "vertex {v} is {d} from an end of edge {e:?}");
        }
    }
    // A kernel edge is in a logical edge where it parts two logical faces,
    // a kernel vertex is a logical one where three meet.
    for e in 0..n.mesh.edges.len() {
        let e = EdgeIx(e as u32);
        let faces = faces_at(n, edge_point(&n.mesh, e));
        assert_eq!(names.edge(e).is_some(), faces.len() >= 2, "edge {e:?}");
    }
    for v in 0..n.adj.vertex_edges.len() {
        let v = VertexIx(v as u32);
        let faces = faces_at(n, vertex_point(&n.mesh, v));
        assert_eq!(names.vertex(v).is_some(), faces.len() >= 3, "vertex {v:?}");
    }
}

// ---------------------------------------------------------------------------
// Names of each kind of body.

fn side(key: &EdgeKey) -> FaceName {
    FaceName {
        feature: feature(1),
        role: Role::Side,
        sources: vec![Source::SketchEdge(key.clone())],
    }
}

fn cap(end: bool, regions: &[&RegionKey]) -> FaceName {
    let mut sources: Vec<Source> = regions
        .iter()
        .map(|r| Source::Region(r.loops[0].clone()))
        .collect();
    sources.sort();
    FaceName {
        feature: feature(1),
        role: Role::Cap { end },
        sources,
    }
}

fn counts(n: &Named) -> (usize, usize, usize) {
    (
        n.names.faces.len(),
        n.names.edges.len(),
        n.names.vertices.len(),
    )
}

#[test]
fn box_sides_are_named_by_their_axes() {
    for (frame, size) in [
        (Frame::WORLD, DVec3::new(1.0, 2.0, 3.0)),
        (tilted(), DVec3::new(-1.0, 0.5, -2.0)),
        (tilted(), DVec3::new(2.0, -0.5, 1.0)),
    ] {
        let body = k().make_box(&frame, size).unwrap();
        let names = BodyNames::new(
            k().origins(&body),
            &k().adjacency(&body),
            &BodySource {
                feature: feature(1),
                made: Made::Box {
                    z_positive: size.z > 0.0,
                },
            },
        )
        .unwrap();
        assert_eq!(
            (names.faces.len(), names.edges.len(), names.vertices.len()),
            (6, 12, 8)
        );
        let mesh = k().tessellate(&body, &opts()).unwrap();
        let local = |x: DVec3| {
            let r = x - frame.origin;
            DVec3::new(r.dot(frame.x), r.dot(frame.y), r.dot(frame.z))
        };
        let (lo, hi) = (DVec3::ZERO.min(size), DVec3::ZERO.max(size));
        let mut roles = BTreeSet::new();
        for f in &names.faces {
            let Role::BoxSide { axis, positive } = f.name.role else {
                panic!("{:?}", f.name)
            };
            assert!(f.name.sources.is_empty());
            assert_eq!(f.pieces.len(), 1);
            let want = if positive { hi } else { lo }[axis as usize];
            for &x in &mesh.faces[f.pieces[0].0 as usize].positions {
                assert!(
                    (local(x)[axis as usize] - want).abs() < 1e-9,
                    "{:?} at {x}",
                    f.name
                );
            }
            roles.insert((axis, positive));
        }
        assert_eq!(roles.len(), 6);
        // Edges between sides of two axes, vertices among sides of three.
        for e in &names.edges {
            assert_eq!(e.name.along, None);
            let axes: BTreeSet<u8> = e
                .name
                .faces
                .iter()
                .map(|f| match f.role {
                    Role::BoxSide { axis, .. } => axis,
                    _ => unreachable!(),
                })
                .collect();
            assert_eq!(axes.len(), 2);
        }
        for v in &names.vertices {
            assert_eq!(v.name.along, None);
            assert_eq!(v.name.faces.len(), 3);
        }
    }
}

#[test]
fn extrude_names() {
    let ks = slot_keys();
    for plane in [Frame::WORLD, tilted()] {
        for sweep in [extrude(0.75, false), extrude(-0.5, true)] {
            let s = slot_and_square(plane, 3.0, 1.0, [true, true]);
            let n = named(&s, sweep);
            // The hole's two curves are one face, its two cap edges one
            // edge each, and no vertex.
            assert_eq!(counts(&n), (13, 26, 16));
            let hole = n
                .names
                .faces
                .iter()
                .find(|f| f.name == side(&ks.hole))
                .unwrap();
            assert_eq!(hole.pieces.len(), 2);
            let (slot, square) = (&s.regions[0], &s.regions[1]);
            for end in [false, true] {
                for r in [slot, square] {
                    assert!(n.names.faces.iter().any(|f| f.name == cap(end, &[r])));
                }
            }
            // Edges and vertices are told apart by their faces alone.
            assert!(n.names.edges.iter().all(|e| e.name.along.is_none()));
            assert!(n.names.vertices.iter().all(|v| v.name.along.is_none()));
            let bottom_top = EdgeName {
                faces: [side(&ks.bottom), cap(true, &[slot])].into(),
                along: None,
            };
            assert!(n.names.edges.iter().any(|e| e.name == bottom_top));
        }
    }
}

/// A lens of two arcs (each less than π): the two edges between its side
/// faces, and the two vertices at each end, share their faces.
fn lens(plane: Frame) -> Sketch {
    let (p, q) = (
        v(&[0, 1]),
        VertexKey {
            ordinal: 1,
            ..v(&[0, 1])
        },
    );
    let ka = edge(0, true, p.clone(), q.clone());
    let kb = edge(1, true, q, p);
    // Arc 0 from (0, -1) to (0, 1) bulging right, arc 1 back bulging
    // left, each of 120°.
    let c = 1.0 / 3.0f64.sqrt();
    let r = 2.0 * c;
    let half = PI / 3.0;
    let mut s = Sketch::new(plane);
    s.curve(arc(DVec2::new(-c, 0.0), r, -half, 2.0 * half), &ka);
    s.curve(arc(DVec2::new(c, 0.0), r, PI - half, 2.0 * half), &kb);
    s.region(&[&[&ka, &kb]]);
    s
}

#[test]
fn edges_and_vertices_sharing_faces_are_told_apart() {
    for plane in [Frame::WORLD, tilted()] {
        let s = lens(plane);
        let n = named(&s, extrude(0.5, false));
        assert_eq!(counts(&n), (4, 6, 4));
        let lateral: Vec<&EdgeName> = n
            .names
            .edges
            .iter()
            .map(|e| &e.name)
            .filter(|e| e.faces.iter().all(|f| f.role == Role::Side))
            .collect();
        assert_eq!(lateral.len(), 2);
        assert!(lateral
            .iter()
            .all(|e| matches!(e.along, Some(Along::Swept(_)))));
        assert!(n
            .names
            .vertices
            .iter()
            .all(|v| matches!(v.name.along, Some(Along::Point { .. }))));
        // Edges along the caps are told apart by their faces.
        assert_eq!(
            n.names
                .edges
                .iter()
                .filter(|e| e.name.along.is_none())
                .count(),
            4
        );
    }
}

/// A rectangle beside the y axis, for revolves about it.
fn ring(plane: Frame, x0: f64, x1: f64) -> Sketch {
    let c = [v(&[0, 3]), v(&[0, 1]), v(&[1, 2]), v(&[2, 3])];
    let keys: Vec<EdgeKey> = (0..4)
        .map(|i| {
            edge(
                i,
                true,
                c[i as usize].clone(),
                c[(i as usize + 1) % 4].clone(),
            )
        })
        .collect();
    let p = [
        DVec2::new(x0, 0.0),
        DVec2::new(x1, 0.0),
        DVec2::new(x1, 1.0),
        DVec2::new(x0, 1.0),
    ];
    let mut s = Sketch::new(plane);
    for i in 0..4 {
        s.curve(line(p[i], p[(i + 1) % 4]), &keys[i]);
    }
    s.region(&[&[&keys[0], &keys[1], &keys[2], &keys[3]]]);
    s
}

/// A half disc on the y axis: its arc and its diameter on the axis.
fn half_disc(plane: Frame) -> Sketch {
    let (bottom, top) = (
        v(&[0, 1]),
        VertexKey {
            ordinal: 1,
            ..v(&[0, 1])
        },
    );
    let ka = edge(0, true, bottom.clone(), top.clone());
    let kd = edge(1, false, top, bottom);
    let mut s = Sketch::new(plane);
    s.curve(arc(DVec2::ZERO, 1.0, -FRAC_PI_2, PI), &ka);
    s.curve(line(DVec2::new(0.0, 1.0), DVec2::new(0.0, -1.0)), &kd);
    s.region(&[&[&ka, &kd]]);
    s
}

/// A circle off the y axis, as the region profile gives it: two halves.
fn disc(plane: Frame) -> Sketch {
    let kc = edge(0, true, v(&[0]), v(&[0]));
    let c = DVec2::new(2.0, 0.0);
    let mut s = Sketch::new(plane);
    s.curve(arc(c, 0.5, 0.0, PI), &kc);
    s.curve(arc(c, 0.5, PI, PI), &kc);
    s.region(&[&[&kc]]);
    s
}

#[test]
fn revolve_names() {
    let y = DVec2::Y;
    for plane in [Frame::WORLD, tilted()] {
        // A partial revolve of a rectangle, of one and two pieces, either
        // way round: four sides and two caps whatever the pieces.
        for angle in [0.5 * PI, 0.9 * PI, 1.1 * PI, -1.2 * PI] {
            let n = named(&ring(plane, 1.0, 2.0), revolve(DVec2::ZERO, y, angle));
            assert_eq!(counts(&n), (6, 12, 8), "angle {angle}");
        }
        // A full revolve: a tube of four faces and four circular edges.
        let n = named(&ring(plane, 1.0, 2.0), revolve(DVec2::ZERO, y, TAU));
        assert_eq!(counts(&n), (4, 4, 0));
        // A half disc revolved: the two points on the axis share their
        // faces and are told apart by the profile points they are.
        for angle in [0.75 * PI, 1.5 * PI] {
            let n = named(&half_disc(plane), revolve(DVec2::ZERO, y, angle));
            assert_eq!(counts(&n), (3, 3, 2), "angle {angle}");
            let alongs: BTreeSet<_> = n.names.vertices.iter().map(|v| &v.name.along).collect();
            assert_eq!(alongs.len(), 2);
            assert!(alongs
                .iter()
                .all(|a| matches!(a, Some(Along::Point { end: false, .. }))));
        }
        // A sphere and a torus: one face, nothing else.
        let n = named(&half_disc(plane), revolve(DVec2::ZERO, y, TAU));
        assert_eq!(counts(&n), (1, 0, 0));
        let n = named(&disc(plane), revolve(DVec2::ZERO, y, TAU));
        assert_eq!(counts(&n), (1, 0, 0));
        // A partial torus: the circle's two halves are one side face.
        let n = named(&disc(plane), revolve(DVec2::ZERO, y, 1.5 * PI));
        assert_eq!(counts(&n), (3, 2, 0));
    }
}

/// The names of a body, as a set.
fn name_set(n: &Named) -> BTreeSet<Name> {
    let mut out = BTreeSet::new();
    out.extend(n.names.faces.iter().map(|f| Name::Face(f.name.clone())));
    out.extend(n.names.edges.iter().map(|e| Name::Edge(e.name.clone())));
    out.extend(
        n.names
            .vertices
            .iter()
            .map(|v| Name::Vertex(v.name.clone())),
    );
    out
}

/// Each name's entity's bounds.
fn bounds(n: &Named) -> Vec<(Name, DVec3, DVec3)> {
    let b = |pts: &mut dyn Iterator<Item = DVec3>| {
        pts.fold(
            (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY)),
            |(lo, hi), p| (lo.min(p), hi.max(p)),
        )
    };
    let mut out = Vec::new();
    for f in &n.names.faces {
        let (lo, hi) = b(&mut f
            .pieces
            .iter()
            .flat_map(|p| n.mesh.faces[p.0 as usize].positions.iter().copied()));
        out.push((Name::Face(f.name.clone()), lo, hi));
    }
    for e in &n.names.edges {
        let (lo, hi) = b(&mut e
            .pieces
            .iter()
            .flat_map(|p| n.mesh.edges[p.0 as usize].iter().copied()));
        out.push((Name::Edge(e.name.clone()), lo, hi));
    }
    for v in &n.names.vertices {
        let x = vertex_point(&n.mesh, v.pieces[0]);
        out.push((Name::Vertex(v.name.clone()), x, x));
    }
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

#[test]
fn names_do_not_change_with_dimensions_curve_order_or_pieces() {
    let base = slot_and_square(Frame::WORLD, 3.0, 1.0, [true, true]);
    let n = named(&base, extrude(0.75, false));
    let want = name_set(&n);
    // Dimensions.
    for (l, r, d) in [(4.0, 1.0, 0.75), (3.0, 1.3, 0.75), (2.0, 0.6, 2.0)] {
        let s = slot_and_square(Frame::WORLD, l, r, [true, true]);
        assert_eq!(name_set(&named(&s, extrude(d, false))), want);
    }
    // Curve order: each entity is where it was.
    let at = bounds(&n);
    for i in 0..base.curves.len() {
        let m = named(&base.rotated(i), extrude(0.75, false));
        for ((a, lo, hi), (b, lo2, hi2)) in at.iter().zip(bounds(&m)) {
            assert_eq!(*a, b);
            assert!(
                (*lo - lo2).length() < 1e-9 && (*hi - hi2).length() < 1e-9,
                "{a:?} moved with curve order"
            );
        }
    }
    let lens0 = named(&lens(Frame::WORLD), extrude(0.5, false));
    let lens1 = named(&lens(Frame::WORLD).rotated(0), extrude(0.5, false));
    assert_eq!(bounds(&lens0).len(), bounds(&lens1).len());
    for ((a, lo, hi), (b, lo2, hi2)) in bounds(&lens0).into_iter().zip(bounds(&lens1)) {
        assert_eq!(a, b);
        assert!((lo - lo2).length() < 1e-9 && (hi - hi2).length() < 1e-9);
    }
    // Pieces: a revolve of less and more than π names the same entities.
    let s = ring(Frame::WORLD, 1.0, 2.0);
    let names: Vec<BTreeSet<Name>> = [0.9 * PI, 1.1 * PI, 1.9 * PI]
        .into_iter()
        .map(|a| name_set(&named(&s, revolve(DVec2::ZERO, DVec2::Y, a))))
        .collect();
    assert_eq!(names[0], names[1]);
    assert_eq!(names[0], names[2]);
}

// ---------------------------------------------------------------------------
// Resolution.

/// Whether a hint lies on a kernel entity of a body.
fn lies_on(n: &Named) -> impl Fn(Pick, &Hint) -> bool + '_ {
    move |pick: Pick, hint: &Hint| {
        let d = match pick {
            Pick::Face(f) => face_distance(&n.mesh, f, hint.point),
            Pick::Edge(e) => edge_distance(&n.mesh, e, hint.point),
            Pick::Vertex(v) => (vertex_point(&n.mesh, v) - hint.point).length(),
        };
        d < 1e-9
    }
}

/// A reference to a kernel entity as the user picks it: its logical
/// entity's name and a point on the piece picked.
fn pick(n: &Named, p: Pick, usage: Use) -> Reference {
    let point = match p {
        Pick::Face(f) => {
            let m = &n.mesh.faces[f.0 as usize];
            let [a, b, c] = m.triangles[m.triangles.len() / 2].map(|i| m.positions[i as usize]);
            (a + b + c) / 3.0
        }
        Pick::Edge(e) => edge_point(&n.mesh, e),
        Pick::Vertex(v) => vertex_point(&n.mesh, v),
    };
    Reference {
        name: n.names.name(p).expect("a named entity"),
        hint: Hint {
            point,
            direction: DVec3::ZERO,
        },
        usage,
    }
}

fn face_named(n: &Named, name: &FaceName) -> Vec<Pick> {
    let f = n.names.faces.iter().find(|f| f.name == *name);
    f.unwrap_or_else(|| panic!("no face {name:?}"))
        .pieces
        .iter()
        .map(|&p| Pick::Face(p))
        .collect()
}

fn resolve(n: &Named, r: &Reference) -> Result<Resolved, ResolveError> {
    n.names.resolve(r, &lies_on(n))
}

/// Every pick of a body, logical entity by logical entity.
fn picks(n: &Named) -> Vec<Vec<Pick>> {
    let mut out: Vec<Vec<Pick>> = Vec::new();
    out.extend(
        n.names
            .faces
            .iter()
            .map(|f| f.pieces.iter().map(|&p| Pick::Face(p)).collect()),
    );
    out.extend(
        n.names
            .edges
            .iter()
            .map(|e| e.pieces.iter().map(|&p| Pick::Edge(p)).collect()),
    );
    out.extend(
        n.names
            .vertices
            .iter()
            .map(|v| v.pieces.iter().map(|&p| Pick::Vertex(p)).collect()),
    );
    out
}

#[test]
fn every_entity_resolves_to_itself() {
    let bodies = [
        named(
            &slot_and_square(tilted(), 3.0, 1.0, [true, true]),
            extrude(0.75, false),
        ),
        named(&lens(Frame::WORLD), extrude(0.5, true)),
        named(
            &ring(tilted(), 1.0, 2.0),
            revolve(DVec2::ZERO, DVec2::Y, 1.5 * PI),
        ),
        named(
            &half_disc(Frame::WORLD),
            revolve(DVec2::ZERO, DVec2::Y, 0.75 * PI),
        ),
    ];
    for n in &bodies {
        for entity in picks(n) {
            for &p in &entity {
                for usage in [Use::Carrier, Use::Set, Use::Single] {
                    let r = pick(n, p, usage);
                    let got = resolve(n, &r).unwrap();
                    assert_eq!(got.pieces, entity, "{:?}", r.name);
                    assert!(!got.moved, "{:?}", r.name);
                }
            }
        }
    }
}

#[test]
fn references_follow_dimension_changes_and_revolve_pieces() {
    let before = named(
        &slot_and_square(Frame::WORLD, 3.0, 1.0, [true, true]),
        extrude(0.75, false),
    );
    let after = named(
        &slot_and_square(Frame::WORLD, 4.0, 1.0, [true, true]),
        extrude(0.75, false),
    );
    let ks = slot_keys();
    for entity in picks(&before) {
        let r = pick(&before, entity[0], Use::Single);
        let got = resolve(&after, &r).unwrap();
        let want: Vec<Pick> = picks(&after)
            .into_iter()
            .find(|e| after.names.name(e[0]).as_ref() == Some(&r.name))
            .unwrap();
        assert_eq!(got.pieces, want);
    }
    // The right end moved with the length; the square did not.
    let right = pick(
        &before,
        face_named(&before, &side(&ks.right))[0],
        Use::Carrier,
    );
    assert!(resolve(&after, &right).unwrap().moved);
    let square = pick(
        &before,
        face_named(&before, &side(&ks.square[1]))[0],
        Use::Carrier,
    );
    let got = resolve(&after, &square).unwrap();
    assert!(!got.moved);

    // A revolve from 0.9π to 1.1π: its side faces become two pieces each,
    // references to them take both.
    let s = ring(Frame::WORLD, 1.0, 2.0);
    let a = named(&s, revolve(DVec2::ZERO, DVec2::Y, 0.9 * PI));
    let b = named(&s, revolve(DVec2::ZERO, DVec2::Y, 1.1 * PI));
    for entity in picks(&a) {
        let r = pick(&a, entity[0], Use::Single);
        let got = resolve(&b, &r).unwrap();
        let want = picks(&b)
            .into_iter()
            .find(|e| b.names.name(e[0]).as_ref() == Some(&r.name))
            .unwrap();
        assert_eq!(got.pieces, want);
    }
    let outer = face_named(&b, &side(&s.keys[1]));
    assert_eq!(outer.len(), 2);
}

/// The slot sketch after a new line 9 was drawn to end where the bottom
/// line ends: that vertex's key, and so the keys of the edges meeting
/// there, change; the geometry does not.
fn slot_extended() -> (Sketch, EdgeKey) {
    let mut s = slot_and_square(Frame::WORLD, 3.0, 1.0, [true, false]);
    let ks = slot_keys();
    let b = v(&[0, 1, 9]);
    let bottom = EdgeKey {
        to: b.clone(),
        ..ks.bottom.clone()
    };
    let right = EdgeKey {
        from: b,
        ..ks.right.clone()
    };
    s.keys[0] = bottom.clone();
    s.keys[1] = right.clone();
    s.regions = vec![region_key(&[
        &[&bottom, &right, &ks.top, &ks.left],
        &[&ks.hole],
    ])];
    (s, bottom)
}

/// The slot sketch after a line 9 split its bottom line at x = 1.
fn slot_split() -> (Sketch, [EdgeKey; 2]) {
    let s0 = slot_and_square(Frame::WORLD, 3.0, 1.0, [true, false]);
    let ks = slot_keys();
    let m = v(&[0, 9]);
    let a = EdgeKey {
        to: m.clone(),
        ..ks.bottom.clone()
    };
    let b = EdgeKey {
        from: m,
        ..ks.bottom.clone()
    };
    let mut s = Sketch::new(Frame::WORLD);
    s.curve(line(DVec2::new(0.0, -1.0), DVec2::new(1.0, -1.0)), &a);
    s.curve(line(DVec2::new(1.0, -1.0), DVec2::new(3.0, -1.0)), &b);
    for i in 1..s0.curves.len() {
        s.curve(s0.curves[i].clone(), &s0.keys[i]);
    }
    s.region(&[&[&a, &b, &ks.right, &ks.top, &ks.left], &[&ks.hole]]);
    (s, [a, b])
}

#[test]
fn references_follow_sketch_edges_split_or_extended() {
    let ks = slot_keys();
    let before = named(
        &slot_and_square(Frame::WORLD, 3.0, 1.0, [true, false]),
        extrude(0.75, false),
    );
    let bottom = face_named(&before, &side(&ks.bottom))[0];
    let top_edge = before
        .names
        .edges
        .iter()
        .find(|e| {
            e.name.faces.contains(&side(&ks.bottom))
                && e.name
                    .faces
                    .iter()
                    .any(|f| f.role == Role::Cap { end: true })
        })
        .unwrap()
        .pieces[0];
    let corner = before
        .names
        .vertices
        .iter()
        .find(|v| {
            v.name.faces.contains(&side(&ks.bottom))
                && v.name.faces.contains(&side(&ks.right))
                && v.name
                    .faces
                    .iter()
                    .any(|f| f.role == Role::Cap { end: true })
        })
        .unwrap()
        .pieces[0];

    // Extended: the one descendant, whatever the use.
    let (s, new_bottom) = slot_extended();
    let after = named(&s, extrude(0.75, false));
    for usage in [Use::Carrier, Use::Set, Use::Single] {
        let got = resolve(&after, &pick(&before, bottom, usage)).unwrap();
        assert_eq!(got.pieces, face_named(&after, &side(&new_bottom)));
        let got = resolve(&after, &pick(&before, Pick::Edge(top_edge), usage)).unwrap();
        let name = after.names.name(got.pieces[0]).unwrap();
        let Name::Edge(e) = name else { panic!() };
        assert!(e.faces.contains(&side(&new_bottom)));
        let got = resolve(&after, &pick(&before, Pick::Vertex(corner), usage)).unwrap();
        assert_eq!(got.pieces.len(), 1);
        assert!(!got.moved, "the corner is where it was");
    }

    // Split by a vertex where nothing else meets it: the pieces look like
    // one face, named by both sketch edges; the reference finds it,
    // whatever its use, and the edge along it.
    let (s, [a, b]) = slot_split();
    let after = named(&s, extrude(0.75, false));
    assert_eq!(counts(&after), counts(&before));
    let mut sources = vec![Source::SketchEdge(a.clone()), Source::SketchEdge(b.clone())];
    sources.sort();
    let merged = FaceName {
        sources,
        ..side(&a)
    };
    let want = face_named(&after, &merged);
    assert_eq!(want.len(), 2);
    for usage in [Use::Carrier, Use::Set, Use::Single] {
        let got = resolve(&after, &pick(&before, bottom, usage)).unwrap();
        assert_eq!(got.pieces, want);
        let got = resolve(&after, &pick(&before, Pick::Edge(top_edge), usage)).unwrap();
        assert_eq!(got.pieces.len(), 2, "one edge of two pieces");
    }
}

/// A bar [0, 4] × [0, 1] of sketch lines 0 (bottom), 1 (right), 2 (top,
/// drawn left to right) and 3 (left, drawn upwards); with `cut`, lines 9
/// and 10 at x = 1.5 and x = 2.5 cut it in three and the outer two are
/// selected, so its bottom line bounds two separate faces. The bottom
/// line's sketch edges (twice the one without the cut).
fn bar(cut: bool) -> (Sketch, [EdgeKey; 2]) {
    let (c0, c1, c2, c3) = (v(&[0, 3]), v(&[0, 1]), v(&[1, 2]), v(&[2, 3]));
    let p = |x: f64, y: f64| DVec2::new(x, y);
    let right = edge(1, true, c1.clone(), c2.clone());
    let left = edge(3, false, c3.clone(), c0.clone());
    let mut s = Sketch::new(Frame::WORLD);
    if !cut {
        let bottom = edge(0, true, c0, c1);
        let top = edge(2, false, c2, c3);
        s.curve(line(p(0.0, 0.0), p(4.0, 0.0)), &bottom);
        s.curve(line(p(4.0, 0.0), p(4.0, 1.0)), &right);
        s.curve(line(p(4.0, 1.0), p(0.0, 1.0)), &top);
        s.curve(line(p(0.0, 1.0), p(0.0, 0.0)), &left);
        s.region(&[&[&bottom, &right, &top, &left]]);
        return (s, [bottom.clone(), bottom]);
    }
    let (m0, m1, n0, n1) = (v(&[0, 9]), v(&[0, 10]), v(&[2, 9]), v(&[2, 10]));
    let bl = edge(0, true, c0.clone(), m0.clone());
    let up = edge(9, true, m0, n0.clone());
    let tl = edge(2, false, n0, c3);
    let br = edge(0, true, m1.clone(), c1);
    let tr = edge(2, false, c2, n1.clone());
    let down = edge(10, false, n1, m1);
    s.curve(line(p(0.0, 0.0), p(1.5, 0.0)), &bl);
    s.curve(line(p(1.5, 0.0), p(1.5, 1.0)), &up);
    s.curve(line(p(1.5, 1.0), p(0.0, 1.0)), &tl);
    s.curve(line(p(0.0, 1.0), p(0.0, 0.0)), &left);
    s.curve(line(p(2.5, 0.0), p(4.0, 0.0)), &br);
    s.curve(line(p(4.0, 0.0), p(4.0, 1.0)), &right);
    s.curve(line(p(4.0, 1.0), p(2.5, 1.0)), &tr);
    s.curve(line(p(2.5, 1.0), p(2.5, 0.0)), &down);
    s.region(&[&[&bl, &up, &tl, &left]]);
    s.region(&[&[&br, &right, &tr, &down]]);
    (s, [bl, br])
}

#[test]
fn references_to_a_face_cut_in_two_take_both_or_the_one_hinted() {
    let (s, [bottom, _]) = bar(false);
    let before = named(&s, extrude(0.5, false));
    let (s, [bl, br]) = bar(true);
    let after = named(&s, extrude(0.5, false));
    assert_eq!(counts(&after), (12, 24, 16));
    let face = face_named(&before, &side(&bottom))[0];
    let both: Vec<Pick> = [
        face_named(&after, &side(&bl)),
        face_named(&after, &side(&br)),
    ]
    .concat();
    for usage in [Use::Carrier, Use::Set] {
        let got = resolve(&after, &pick(&before, face, usage)).unwrap();
        assert_eq!(got.pieces, both);
    }
    let mut r = pick(&before, face, Use::Single);
    for (x, want) in [(0.5, &bl), (3.0, &br)] {
        r.hint.point = DVec3::new(x, 0.0, 0.25);
        let got = resolve(&after, &r).unwrap();
        assert_eq!(got.pieces, face_named(&after, &side(want)));
        assert!(!got.moved);
    }
    r.hint.point = DVec3::new(2.0, 0.0, 0.25);
    assert_eq!(resolve(&after, &r), Err(ResolveError::Ambiguous(2)));
    // The edge between it and the left side: one descendant.
    let left = &s.keys[3];
    let corner = |n: &Named, bottom: &EdgeKey| {
        n.names
            .edges
            .iter()
            .find(|e| e.name.faces == [side(bottom), side(left)].into())
            .unwrap()
            .pieces
            .clone()
    };
    let got = resolve(
        &after,
        &pick(
            &before,
            Pick::Edge(corner(&before, &bottom)[0]),
            Use::Single,
        ),
    )
    .unwrap();
    let want: Vec<Pick> = corner(&after, &bl).into_iter().map(Pick::Edge).collect();
    assert_eq!(got.pieces, want);
    // The edge between it and the top cap: the cap was selected anew, a
    // region of other sketch edges, so nothing descends from it and the
    // reference fails rather than guess (a cap keeps its name while the
    // stored region key finds its region; see caps_follow_their_regions).
    let edge = before
        .names
        .edges
        .iter()
        .find(|e| {
            e.name.faces.contains(&side(&bottom))
                && e.name
                    .faces
                    .iter()
                    .any(|f| f.role == Role::Cap { end: true })
        })
        .unwrap()
        .pieces[0];
    let r = pick(&before, Pick::Edge(edge), Use::Set);
    assert_eq!(resolve(&after, &r), Err(ResolveError::Missing));
}

#[test]
fn references_to_what_is_gone_fail() {
    let ks = slot_keys();
    let both = named(
        &slot_and_square(Frame::WORLD, 3.0, 1.0, [true, true]),
        extrude(0.75, false),
    );
    let slot = named(
        &slot_and_square(Frame::WORLD, 3.0, 1.0, [true, false]),
        extrude(0.75, false),
    );
    // The square is no longer extruded: its faces, edges and vertices are
    // gone, and nothing of the slot stands in for them.
    for entity in picks(&both) {
        let name = both.names.name(entity[0]).unwrap();
        let of_square = match &name {
            Name::Face(f) => f.sources.iter().any(|s| match s {
                Source::SketchEdge(k) => ks.square.contains(k),
                Source::Region(l) => l[0] == ks.square[0],
            }),
            Name::Edge(e) => e.faces.iter().any(|f| {
                f.sources
                    .iter()
                    .any(|s| matches!(s, Source::SketchEdge(k) if ks.square.contains(k)))
            }),
            Name::Vertex(v) => v.faces.iter().any(|f| {
                f.sources
                    .iter()
                    .any(|s| matches!(s, Source::SketchEdge(k) if ks.square.contains(k)))
            }),
        };
        let r = pick(&both, entity[0], Use::Single);
        if of_square {
            assert_eq!(resolve(&slot, &r), Err(ResolveError::Missing), "{name:?}");
        } else {
            assert!(resolve(&slot, &r).is_ok(), "{name:?}");
        }
    }
    // A side face whose sketch edge was deleted, picked where another side
    // now is: never substituted.
    let mut r = pick(
        &both,
        face_named(&both, &side(&ks.square[0]))[0],
        Use::Carrier,
    );
    let bottom = face_named(&slot, &side(&ks.bottom))[0];
    r.hint = pick(&slot, bottom, Use::Carrier).hint;
    assert_eq!(resolve(&slot, &r), Err(ResolveError::Missing));
    // Another feature's face.
    let mut r = pick(&both, bottom, Use::Carrier);
    let Name::Face(f) = &mut r.name else { panic!() };
    f.feature = feature(2);
    assert_eq!(resolve(&slot, &r), Err(ResolveError::Missing));
}

/// Two unit squares side by side, sharing sketch curve 3 (x = 1); either
/// or both selected.
fn two_squares(left: bool, right: bool) -> (Sketch, [RegionKey; 2]) {
    // Corners: 0 (0,0), 1 (1,0), 2 (2,0), 3 (2,1), 4 (1,1), 5 (0,1).
    // Curves: 0 bottom (0,0)-(2,0), 1 right x=2, 2 top (2,1)-(0,1),
    // 3 middle (1,0)-(1,1), 4 left x=0.
    let c = [
        v(&[0, 4]),
        v(&[0, 3]),
        v(&[0, 1]),
        v(&[1, 2]),
        v(&[2, 3]),
        v(&[2, 4]),
    ];
    let e = |curve, fwd, a: usize, b: usize| edge(curve, fwd, c[a].clone(), c[b].clone());
    let (b0, b1, r, t1, t0, l) = (
        e(0, true, 0, 1),
        e(0, true, 1, 2),
        e(1, true, 2, 3),
        e(2, true, 3, 4),
        e(2, true, 4, 5),
        e(4, false, 5, 0),
    );
    let (m_up, m_down) = (e(3, true, 1, 4), e(3, false, 4, 1));
    let keys = [
        region_key(&[&[&b0, &m_up, &t0, &l]]),
        region_key(&[&[&b1, &r, &t1, &m_down]]),
    ];
    let p = |i: usize| {
        [
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 0.0),
            DVec2::new(2.0, 0.0),
            DVec2::new(2.0, 1.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(0.0, 1.0),
        ][i]
    };
    let mut s = Sketch::new(Frame::WORLD);
    let mut add = |a: usize, b: usize, key: &EdgeKey| s.curve(line(p(a), p(b)), key);
    match (left, right) {
        (true, false) => {
            add(0, 1, &b0);
            add(1, 4, &m_up);
            add(4, 5, &t0);
            add(5, 0, &l);
        }
        (false, true) => {
            add(1, 2, &b1);
            add(2, 3, &r);
            add(3, 4, &t1);
            add(4, 1, &m_down);
        }
        _ => {
            add(0, 1, &b0);
            add(1, 2, &b1);
            add(2, 3, &r);
            add(3, 4, &t1);
            add(4, 5, &t0);
            add(5, 0, &l);
        }
    }
    s.regions = [left, right]
        .iter()
        .zip(&keys)
        .filter(|(on, _)| **on)
        .map(|(_, k)| k.clone())
        .collect();
    (s, keys)
}

#[test]
fn caps_follow_their_regions() {
    let (s, [left, right]) = two_squares(true, false);
    let one = named(&s, extrude(0.5, false));
    let top = face_named(&one, &cap(true, &[&left]))[0];
    let r = pick(&one, top, Use::Single);

    // Both selected: one cap covering both; the reference finds it.
    let (s, _) = two_squares(true, true);
    let merged = named(&s, extrude(0.5, false));
    // The bottom and top lines' pieces each look like one face.
    assert_eq!(counts(&merged), (6, 12, 8));
    assert_eq!(
        merged
            .names
            .faces
            .iter()
            .filter(|f| f.name.sources.len() == 2 && f.name.role == Role::Side)
            .count(),
        2
    );
    let got = resolve(&merged, &r).unwrap();
    assert_eq!(
        got.pieces,
        face_named(&merged, &cap(true, &[&left, &right]))
    );
    assert!(!got.moved);
    // And back: a reference to the merged cap finds the left one's.
    let r2 = pick(&merged, got.pieces[0], Use::Single);
    assert_eq!(
        resolve(&one, &r2).unwrap().pieces,
        face_named(&one, &cap(true, &[&left]))
    );
    // A cap of a region not selected is never the reference's.
    let (s, _) = two_squares(false, true);
    let other = named(&s, extrude(0.5, false));
    assert_eq!(resolve(&other, &r), Err(ResolveError::Missing));

    // Adding a separate region keeps the slot's cap's name.
    let slot = named(
        &slot_and_square(Frame::WORLD, 3.0, 1.0, [true, false]),
        extrude(0.75, false),
    );
    let both = named(
        &slot_and_square(Frame::WORLD, 3.0, 1.0, [true, true]),
        extrude(0.75, false),
    );
    let r = pick(
        &slot,
        face_named(
            &slot,
            &cap(
                false,
                &[&slot_and_square(Frame::WORLD, 3.0, 1.0, [true, false]).regions[0]],
            ),
        )[0],
        Use::Single,
    );
    let got = resolve(&both, &r).unwrap();
    assert_eq!(both.names.name(got.pieces[0]), Some(r.name.clone()));

    // A cap keeps its name while the stored key finds its region, though
    // the region's edges changed (the bottom line extended).
    let stored = slot_and_square(Frame::WORLD, 3.0, 1.0, [true, false]).regions;
    let (s, _) = slot_extended();
    let n = named_as(&s, extrude(0.75, false), &stored);
    let got = resolve(&n, &r).unwrap();
    assert_eq!(n.names.name(got.pieces[0]), Some(r.name.clone()));
    // Selected anew, the stored key is the new one: found by its curves.
    let n = named(&s, extrude(0.75, false));
    let got = resolve(&n, &r).unwrap();
    assert_ne!(n.names.name(got.pieces[0]), Some(r.name.clone()));
    assert!(matches!(got.pieces[0], Pick::Face(_)));
}

#[test]
fn inconsistent_inputs_are_refused() {
    let s = slot_and_square(Frame::WORLD, 3.0, 1.0, [true, true]);
    let body = make(&s, extrude(0.75, false));
    // The square's region not selected.
    let mut t = s.clone();
    t.regions.pop();
    assert!(matches!(
        name_body(&t, &body, feature(1), &t.regions),
        Err(NamingError::Inconsistent(_))
    ));
    // A curve without a sketch edge.
    let mut t = s.clone();
    t.keys.pop();
    assert_eq!(
        name_body(&t, &body, feature(1), &t.regions).unwrap_err(),
        NamingError::MissingCurve(s.keys.len() as u32 - 1)
    );
    // A selected region the profile does not have.
    let mut t = s.clone();
    t.regions.push(two_squares(true, false).1[0].clone());
    assert!(matches!(
        name_body(&t, &body, feature(1), &t.regions),
        Err(NamingError::Inconsistent(_))
    ));
    // Two curves with one sketch edge that are not one face.
    let mut t = s.clone();
    let ks = slot_keys();
    t.keys[0] = ks.top.clone();
    assert!(name_body(&t, &body, feature(1), &t.regions).is_err());
}
