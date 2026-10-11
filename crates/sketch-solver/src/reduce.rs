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

use crate::build::{Problem, Sub};

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
    // Components over the numerical equations.
    let nu = unknowns.len();
    let mut cuf: Vec<usize> = (0..nu).collect();
    let mut eq_unknowns: Vec<Vec<usize>> = Vec::with_capacity(prob.eqs.len());
    let mut referenced = vec![false; nu];
    let mut constant_eqs = Vec::new();
    for (i, eq) in prob.eqs.iter().enumerate() {
        let mut us: Vec<usize> = Vec::new();
        if matches!(status[i], Status::General | Status::Demoted) {
            for s in eq.kind.slots(&prob.splines) {
                if let Some(u) = unknown_of_class[class_of[s]] {
                    us.push(u);
                }
            }
            us.sort_unstable();
            us.dedup();
            if us.is_empty() {
                constant_eqs.push(i);
            }
            for &u in &us {
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
        eq_unknowns.push(us);
    }
    let mut comp_of_root = vec![usize::MAX; nu];
    let mut components: Vec<Component> = Vec::new();
    let mut free_unknowns = Vec::new();
    for (u, &is_referenced) in referenced.iter().enumerate().take(nu) {
        if !is_referenced {
            free_unknowns.push(u);
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
    Reduced {
        status,
        class_of,
        classes,
        unknowns,
        unknown_of_class,
        components,
        constant_eqs,
        free_unknowns,
    }
}
