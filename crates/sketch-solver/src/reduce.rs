//! Elimination of trivial equations and decomposition into components
//! (design note 0005 §5 step 2).
//!
//! Equalities of two unknowns (coincident points, horizontal and vertical
//! constraints, concentric centres, equal circle radii) and unknowns fixed
//! to constants (fixed points, circle radius and diameter dimensions) are
//! eliminated by substitution, as in SolveSpace's solver: the unknowns of an
//! equality class are replaced by one, and constant classes leave the
//! system. This is exact (the eliminated equations are linear and are
//! satisfied identically) and it typically removes most of a sketch's
//! unknowns. Rank and null space transfer exactly: with `M` the matrix that
//! expands class values to slots, `null(J_full) = M · null(J_reduced)` and
//! `rank(J_full) = (independent eliminations) + rank(J_reduced)`.
//!
//! Equations are processed in priority order. An elimination already implied
//! by earlier ones is redundant; one that contradicts them (two different
//! constants for one class) is *demoted*: it stays in the numerical system
//! as an ordinary equation, where it shows up as a conflict.
//!
//! The remaining equations are split into connected components of the
//! equation–unknown incidence graph; each is solved on its own.
//!
//! Then images of derived constraints that only their defining linear
//! equation references (a pattern's copy points) are *peeled* ([`peel`]):
//! such an equation is the identity in the image, so it fixes the image
//! explicitly from the other unknowns and leaves the numerical system,
//! and the image is computed after the solve by back-substitution.
//! Peeling repeats, so copies of copies peel too. It is exact, and rank
//! and null space transfer exactly: each peeled equation adds one to the
//! rank and one unknown, and the image's null-space row is the same linear
//! combination of the other unknowns' rows. Without it, a pattern of a
//! hundred hexagons is one dense component of about 2,400 unknowns.

use crate::build::{Problem, Sub};
use crate::eq::EqKind;

/// What happened to an equation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Status {
    /// An ordinary equation, solved numerically.
    General,
    /// Eliminated by substitution.
    Accepted,
    /// An elimination implied by earlier eliminations.
    Redundant,
    /// An elimination that contradicts earlier ones; solved numerically.
    Demoted,
    /// Peeled: the only equation referencing an image unknown, which it
    /// fixes explicitly ([`peel`]).
    Explicit,
}

/// A peeled equation and the unknown it fixes.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Peeled {
    /// The equation (linear).
    pub eq: usize,
    /// The reduced unknown it fixes.
    pub unknown: usize,
    /// The sum of the equation's coefficients on that unknown's slots.
    pub coef: f64,
}

/// An equality class of slots.
#[derive(Clone, Debug)]
pub(crate) struct Class {
    /// Member slots, ascending.
    pub members: Vec<usize>,
    /// The constant value (scaled, unscaled), if fixed.
    pub constant: Option<(f64, f64)>,
}

/// A connected set of unknowns and the equations over them.
#[derive(Clone, Debug)]
pub(crate) struct Component {
    /// Reduced unknowns, ascending.
    pub unknowns: Vec<usize>,
    /// Equations, in priority order.
    pub eqs: Vec<usize>,
}

/// The reduced system.
#[derive(Clone, Debug)]
pub(crate) struct Reduced {
    /// Per equation.
    pub status: Vec<Status>,
    /// Slot → class.
    pub class_of: Vec<usize>,
    /// Classes, ordered by smallest member slot.
    pub classes: Vec<Class>,
    /// Reduced unknown → class.
    pub unknowns: Vec<usize>,
    /// Class → reduced unknown (`None` for constant classes).
    pub unknown_of_class: Vec<Option<usize>>,
    /// Components with at least one equation.
    pub components: Vec<Component>,
    /// Numerical equations with no unknowns left (all their slots constant).
    pub constant_eqs: Vec<usize>,
    /// Unknowns that no numerical equation references.
    pub free_unknowns: Vec<usize>,
    /// Peeled equations in back-substitution order: each one's unknown
    /// depends only on unknowns solved numerically, free, or peeled before
    /// it.
    pub peeled: Vec<Peeled>,
}

fn find(uf: &mut [usize], mut x: usize) -> usize {
    while uf[x] != x {
        uf[x] = uf[uf[x]];
        x = uf[x];
    }
    x
}

/// Eliminates substitution equations and decomposes the rest. Two constants
/// within `const_tol` (scaled) are equal.
pub(crate) fn reduce(prob: &Problem, const_tol: f64) -> Reduced {
    let n = prob.n_slots();
    let mut uf: Vec<usize> = (0..n).collect();
    // Constant attached to each root.
    let mut constant: Vec<Option<(f64, f64)>> = vec![None; n];
    let mut status = Vec::with_capacity(prob.eqs.len());
    for eq in &prob.eqs {
        let st = match eq.sub {
            None => Status::General,
            Some(Sub::Equal(a, b)) => {
                let (ra, rb) = (find(&mut uf, a), find(&mut uf, b));
                if ra == rb {
                    Status::Redundant
                } else {
                    match (constant[ra], constant[rb]) {
                        (Some(ca), Some(cb)) => {
                            if (ca.0 - cb.0).abs() <= const_tol {
                                Status::Redundant
                            } else {
                                Status::Demoted
                            }
                        }
                        (ca, cb) => {
                            let (lo, hi) = (ra.min(rb), ra.max(rb));
                            uf[hi] = lo;
                            constant[lo] = ca.or(cb);
                            constant[hi] = None;
                            Status::Accepted
                        }
                    }
                }
            }
            Some(Sub::Const { slot, value, abs }) => {
                let r = find(&mut uf, slot);
                match constant[r] {
                    Some(c) => {
                        if (c.0 - value).abs() <= const_tol {
                            Status::Redundant
                        } else {
                            Status::Demoted
                        }
                    }
                    None => {
                        constant[r] = Some((value, abs));
                        Status::Accepted
                    }
                }
            }
        };
        status.push(st);
    }
    // Classes in order of their smallest slot (= the root, by construction).
    let mut class_of = vec![usize::MAX; n];
    let mut classes: Vec<Class> = Vec::new();
    for s in 0..n {
        let r = find(&mut uf, s);
        if class_of[r] == usize::MAX {
            class_of[r] = classes.len();
            classes.push(Class {
                members: Vec::new(),
                constant: constant[r],
            });
        }
        let c = class_of[r];
        class_of[s] = c;
        classes[c].members.push(s);
    }
    let mut unknowns = Vec::new();
    let mut unknown_of_class = vec![None; classes.len()];
    for (c, cl) in classes.iter().enumerate() {
        if cl.constant.is_none() {
            unknown_of_class[c] = Some(unknowns.len());
            unknowns.push(c);
        }
    }
    let (components, constant_eqs, free_unknowns) = decompose(
        prob,
        &status,
        &class_of,
        &unknown_of_class,
        unknowns.len(),
        &[],
    );
    Reduced {
        status,
        class_of,
        classes,
        unknowns,
        unknown_of_class,
        components,
        constant_eqs,
        free_unknowns,
        peeled: Vec::new(),
    }
}

/// The reduced unknowns of each equation that is solved numerically.
fn equation_unknowns(
    prob: &Problem,
    status: &[Status],
    class_of: &[usize],
    unknown_of_class: &[Option<usize>],
) -> Vec<Vec<usize>> {
    prob.eqs
        .iter()
        .enumerate()
        .map(|(i, eq)| {
            let mut us: Vec<usize> = Vec::new();
            if matches!(status[i], Status::General | Status::Demoted) {
                for s in eq.kind.slots(&prob.splines) {
                    if let Some(u) = unknown_of_class[class_of[s]] {
                        us.push(u);
                    }
                }
                us.sort_unstable();
                us.dedup();
            }
            us
        })
        .collect()
}

/// Connected components of the numerical equations, the equations left
/// with no unknowns, and the unknowns no numerical equation references
/// (`peeled` unknowns are neither: they are computed afterwards).
fn decompose(
    prob: &Problem,
    status: &[Status],
    class_of: &[usize],
    unknown_of_class: &[Option<usize>],
    nu: usize,
    peeled: &[usize],
) -> (Vec<Component>, Vec<usize>, Vec<usize>) {
    let mut cuf: Vec<usize> = (0..nu).collect();
    let eq_unknowns = equation_unknowns(prob, status, class_of, unknown_of_class);
    let mut referenced = vec![false; nu];
    for &u in peeled {
        referenced[u] = true;
    }
    let mut constant_eqs = Vec::new();
    for (i, us) in eq_unknowns.iter().enumerate() {
        if matches!(status[i], Status::General | Status::Demoted) && us.is_empty() {
            constant_eqs.push(i);
        }
        for &u in us {
            referenced[u] = true;
        }
        for w in us.windows(2) {
            let (a, b) = (find(&mut cuf, w[0]), find(&mut cuf, w[1]));
            if a != b {
                let (lo, hi) = (a.min(b), a.max(b));
                cuf[hi] = lo;
            }
        }
    }
    let mut is_peeled = vec![false; nu];
    for &u in peeled {
        is_peeled[u] = true;
    }
    let mut comp_of_root = vec![usize::MAX; nu];
    let mut components: Vec<Component> = Vec::new();
    let mut free_unknowns = Vec::new();
    for (u, &is_referenced) in referenced.iter().enumerate().take(nu) {
        if !is_referenced {
            free_unknowns.push(u);
            continue;
        }
        if is_peeled[u] {
            continue;
        }
        let r = find(&mut cuf, u);
        if comp_of_root[r] == usize::MAX {
            comp_of_root[r] = components.len();
            components.push(Component {
                unknowns: Vec::new(),
                eqs: Vec::new(),
            });
        }
        components[comp_of_root[r]].unknowns.push(u);
    }
    for (i, us) in eq_unknowns.iter().enumerate() {
        if let Some(&u) = us.first() {
            let r = find(&mut cuf, u);
            components[comp_of_root[r]].eqs.push(i);
        }
    }
    (components, constant_eqs, free_unknowns)
}

/// Peels the linear equations that alone reference an unknown marked
/// `peelable` (an image of a derived constraint, not dragged): the
/// equation fixes it explicitly, so both leave the numerical system
/// (module documentation). Repeats until nothing more peels; the order is
/// recorded reversed, for back-substitution.
pub(crate) fn peel(prob: &Problem, red: &mut Reduced, peelable: &[bool]) {
    let nu = red.unknowns.len();
    let eq_unknowns = equation_unknowns(prob, &red.status, &red.class_of, &red.unknown_of_class);
    let mut uses = vec![0usize; nu];
    let mut eqs_of: Vec<Vec<usize>> = vec![Vec::new(); nu];
    for (i, us) in eq_unknowns.iter().enumerate() {
        for &u in us {
            uses[u] += 1;
            eqs_of[u].push(i);
        }
    }
    let mut order: Vec<Peeled> = Vec::new();
    let mut gone = vec![false; prob.eqs.len()];
    let mut queue: Vec<usize> = (0..nu).filter(|&u| peelable[u] && uses[u] == 1).collect();
    while let Some(u) = queue.pop() {
        if uses[u] != 1 {
            continue;
        }
        let Some(&i) = eqs_of[u].iter().find(|&&i| !gone[i]) else {
            continue;
        };
        if red.status[i] != Status::General {
            continue;
        }
        let EqKind::Lin { terms, .. } = &prob.eqs[i].kind else {
            continue;
        };
        let coef: f64 = terms
            .iter()
            .filter(|&&(s, _)| red.unknown_of_class[red.class_of[s]] == Some(u))
            .map(|&(_, c)| c)
            .sum();
        let largest = terms.iter().map(|&(_, c)| c.abs()).fold(0.0, f64::max);
        // The images' own coefficient is 1 (patterns) or ±1 (offset
        // radii); anything much smaller is a merged class and not peeled.
        if !coef.is_finite() || coef.abs() < 1e-3 * largest {
            continue;
        }
        gone[i] = true;
        order.push(Peeled {
            eq: i,
            unknown: u,
            coef,
        });
        for &w in &eq_unknowns[i] {
            uses[w] -= 1;
            if w != u && peelable[w] && uses[w] == 1 {
                queue.push(w);
            }
        }
    }
    if order.is_empty() {
        return;
    }
    for p in &order {
        red.status[p.eq] = Status::Explicit;
    }
    order.reverse();
    let peeled: Vec<usize> = order.iter().map(|p| p.unknown).collect();
    let (components, constant_eqs, free_unknowns) = decompose(
        prob,
        &red.status,
        &red.class_of,
        &red.unknown_of_class,
        nu,
        &peeled,
    );
    red.components = components;
    red.constant_eqs = constant_eqs;
    red.free_unknowns = free_unknowns;
    red.peeled = order;
}
