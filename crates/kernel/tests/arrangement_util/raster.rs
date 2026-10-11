//! A raster oracle for arrangements of arbitrary curves.
//!
//! The curves are drawn as walls about three pixels thick on a fine grid,
//! and the free pixels are flooded into 4-connected components. A free
//! pixel is more than a pixel from every curve, so each component lies
//! inside a single face of the arrangement, and the component touching the
//! frame lies in the unbounded face. When no two curves come closer than a
//! few pixels without touching (true for the lattice configurations of the
//! property tests), every face of more than a few pixels is one large
//! component, so large components and bounded faces correspond one to one.

use std::collections::VecDeque;

use openshape_kernel::{DVec2, ProfileCurve};

use super::{curve_point, domain};

/// Flooded raster of a set of curves.
pub struct Raster {
    lo: DVec2,
    px: f64,
    nx: usize,
    ny: usize,
    /// Component of each pixel, `u32::MAX` on a wall.
    comp: Vec<u32>,
    /// Pixel count of each component.
    pub sizes: Vec<usize>,
    /// True for components touching the frame (in the unbounded face).
    pub outer: Vec<bool>,
}

impl Raster {
    /// Rasterizes `curves` over `[lo, hi]` with square pixels of side `px`.
    pub fn new(curves: &[ProfileCurve], lo: DVec2, hi: DVec2, px: f64) -> Raster {
        let nx = ((hi.x - lo.x) / px).ceil() as usize;
        let ny = ((hi.y - lo.y) / px).ceil() as usize;
        let mut wall = vec![false; nx * ny];
        for c in curves {
            let (a, b) = domain(c);
            let coarse: Vec<DVec2> = (0..=256)
                .map(|i| curve_point(c, a + (b - a) * i as f64 / 256.0))
                .collect();
            let len: f64 = coarse.windows(2).map(|w| (w[1] - w[0]).length()).sum();
            // Steps of at most a quarter pixel (with a margin for uneven
            // parameterization).
            let n = ((len / px) * 16.0).ceil().max(64.0) as usize;
            for i in 0..=n {
                let p = curve_point(c, a + (b - a) * i as f64 / n as f64);
                let ix = ((p.x - lo.x) / px).floor() as i64;
                let iy = ((p.y - lo.y) / px).floor() as i64;
                for dx in -1..=1 {
                    for dy in -1..=1 {
                        let (x, y) = (ix + dx, iy + dy);
                        if x >= 0 && y >= 0 && (x as usize) < nx && (y as usize) < ny {
                            wall[y as usize * nx + x as usize] = true;
                        }
                    }
                }
            }
        }
        let mut comp = vec![u32::MAX; nx * ny];
        let mut sizes = Vec::new();
        let mut outer = Vec::new();
        let mut q = VecDeque::new();
        for s in 0..nx * ny {
            if wall[s] || comp[s] != u32::MAX {
                continue;
            }
            let id = sizes.len() as u32;
            comp[s] = id;
            q.push_back(s);
            let (mut size, mut touches) = (0, false);
            while let Some(k) = q.pop_front() {
                size += 1;
                let (x, y) = (k % nx, k / nx);
                if x == 0 || y == 0 || x == nx - 1 || y == ny - 1 {
                    touches = true;
                }
                let mut push = |n: usize| {
                    if !wall[n] && comp[n] == u32::MAX {
                        comp[n] = id;
                        q.push_back(n);
                    }
                };
                if x > 0 {
                    push(k - 1);
                }
                if x + 1 < nx {
                    push(k + 1);
                }
                if y > 0 {
                    push(k - nx);
                }
                if y + 1 < ny {
                    push(k + nx);
                }
            }
            sizes.push(size);
            outer.push(touches);
        }
        Raster {
            lo,
            px,
            nx,
            ny,
            comp,
            sizes,
            outer,
        }
    }

    /// Up to `per` sample points of each component, at pixel centres whose
    /// 5 × 5 neighbourhood is free (so at least two pixels from any curve),
    /// spread over the component.
    pub fn samples(&self, per: usize) -> Vec<Vec<DVec2>> {
        let mut deep: Vec<Vec<DVec2>> = vec![Vec::new(); self.sizes.len()];
        for y in 2..self.ny.saturating_sub(2) {
            for x in 2..self.nx.saturating_sub(2) {
                let id = self.comp[y * self.nx + x];
                if id == u32::MAX {
                    continue;
                }
                let clear = (y - 2..=y + 2)
                    .all(|yy| (x - 2..=x + 2).all(|xx| self.comp[yy * self.nx + xx] == id));
                if clear {
                    deep[id as usize]
                        .push(self.lo + DVec2::new(x as f64 + 0.5, y as f64 + 0.5) * self.px);
                }
            }
        }
        deep.into_iter()
            .map(|d| {
                if d.len() <= per {
                    d
                } else {
                    (0..per)
                        .map(|i| d[i * d.len() / per + d.len() / (2 * per)])
                        .collect()
                }
            })
            .collect()
    }
}
