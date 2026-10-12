//! Headless end-to-end demo: builds parts through the op log, renders them
//! offscreen with the shipping renderer, orbits the camera, and checks that
//! undo/redo, save/open and STL export behave. Used by CI (on Mesa's software
//! Vulkan) and runnable with `opendatum --headless-demo <dir>`.

use std::path::{Path, PathBuf};

use opendatum_core_model::{
    CurveDef, DVec2, DVec3, Feature, PlaneRef, ProfileDef, ProfileSource, RevolveAxis,
};
use opendatum_kernel::Kernel;
use opendatum_render::{write_png, Camera, OffscreenRenderer, RgbaImage};

use crate::dialogs::{FeatureDialog, Template};
use crate::Model;

/// Result of the demo.
#[derive(Debug, Default)]
pub struct DemoReport {
    /// Adapter description (name, backend, driver).
    pub adapter: String,
    /// Images written.
    pub images: Vec<PathBuf>,
    /// Named checks and whether they passed.
    pub checks: Vec<(String, bool)>,
}

impl DemoReport {
    /// True if every check passed.
    pub fn ok(&self) -> bool {
        !self.checks.is_empty() && self.checks.iter().all(|(_, ok)| *ok)
    }

    fn check(&mut self, name: &str, ok: bool) {
        self.checks.push((name.to_string(), ok));
    }
}

const W: u32 = 800;
const H: u32 = 600;

fn differs(a: &RgbaImage, b: &RgbaImage) -> usize {
    a.pixels
        .chunks(4)
        .zip(b.pixels.chunks(4))
        .filter(|(p, q)| p != q)
        .count()
}

/// A spline vase revolved about a vertical axis at `x = x0` (front plane).
pub fn vase_at(x0: f64) -> Feature {
    let (h, r) = (0.05, 0.012);
    let p = |x: f64, y: f64| DVec2::new(x0 + x, y);
    Feature::Revolve {
        profile: ProfileSource::Inline(ProfileDef {
            plane: PlaneRef::Front,
            curves: vec![
                CurveDef::Spline {
                    through: vec![
                        p(r, 0.0),
                        p(1.6 * r, 0.3 * h),
                        p(0.8 * r, 0.7 * h),
                        p(1.1 * r, h),
                    ],
                },
                CurveDef::Line {
                    start: p(1.1 * r, h),
                    end: p(0.0, h),
                },
                CurveDef::Line {
                    start: p(0.0, h),
                    end: p(0.0, 0.0),
                },
                CurveDef::Line {
                    start: p(0.0, 0.0),
                    end: p(r, 0.0),
                },
            ],
        }),
        axis: RevolveAxis::Inline {
            origin: DVec2::new(x0, 0.0),
            dir: DVec2::Y,
        },
        angle: std::f64::consts::TAU,
    }
}

/// The demo's parts, in insertion order: box, cylinder, ring, vase.
fn demo_parts() -> Result<Vec<(&'static str, Feature)>, String> {
    let mut box_dialog = FeatureDialog::new(Template::Box);
    box_dialog.fields[0].text = "40 mm".into();
    box_dialog.fields[1].text = "30 mm".into();
    box_dialog.fields[2].text = "20 mm".into();
    let mut box_feature = box_dialog.build()?;
    if let Feature::Box { corner, .. } = &mut box_feature {
        *corner = DVec3::new(-0.02, -0.015, 0.0);
    }
    let mut cyl = FeatureDialog::new(Template::Cylinder);
    cyl.fields[2].text = "60 mm".into();
    cyl.fields[3].text = "0 mm".into();
    let mut ring = FeatureDialog::new(Template::Ring);
    ring.fields[0].text = "35 mm".into();
    ring.fields[1].text = "4 mm".into();
    ring.fields[2].text = "270 deg".into();
    Ok(vec![
        ("Box 1", box_feature),
        ("Cylinder 1", cyl.build()?),
        ("Ring 1", ring.build()?),
        ("Vase 1", vase_at(-0.07)),
    ])
}

/// Runs the demo, writing PNGs, an STL and a native document into `out`.
pub fn run(out: &Path) -> Result<DemoReport, String> {
    std::fs::create_dir_all(out).map_err(|e| e.to_string())?;
    let mut report = DemoReport::default();
    let mut gpu = OffscreenRenderer::new().map_err(|e| e.to_string())?;
    let info = gpu.adapter_info();
    report.adapter = format!("{} ({:?}, {})", info.name, info.backend, info.driver_info);

    let save = |name: &str, img: &RgbaImage, report: &mut DemoReport| -> Result<(), String> {
        let p = out.join(name);
        write_png(&p, img).map_err(|e| e.to_string())?;
        report.images.push(p);
        Ok(())
    };
    let mut model = Model::new();
    let mut render = |model: &Model, camera: &Camera| -> Result<RgbaImage, String> {
        gpu.render_meshes_rgba(model.meshes(), camera, W, H)
            .map_err(|e| e.to_string())
    };

    // 1. A kernel-made box, created through the op log.
    let (n, f) = demo_parts()?.into_iter().next().ok_or("parts")?;
    model.add_feature(n, f);
    report.check("box regenerates", model.bodies().bodies.len() == 1);
    let volume = model
        .kernel()
        .mass_properties(&model.bodies().bodies[0].body, 1e-12)
        .map(|m| m.volume)
        .unwrap_or(0.0);
    report.check("box volume is 24 cm³", (volume - 24e-6).abs() < 1e-12);
    // Frame the region all demo parts will occupy so the camera stays fixed
    // and undo/redo images can be compared pixel for pixel.
    let mut scratch = Model::new();
    for (name, f) in demo_parts()? {
        scratch.add_feature(name, f);
    }
    let (lo, hi) = scratch.bounds().ok_or("no geometry")?;
    let mut camera = Camera::isometric();
    camera.fit(lo, hi);
    let a = render(&model, &camera)?;
    save("01-box.png", &a, &mut report)?;
    let bg = a.pixel(2, 2);
    let center = a.pixel(W / 2, H / 2);
    report.check(
        "box is drawn (centre differs from background)",
        center != bg,
    );

    // 2. Orbit.
    camera.orbit(140.0, -60.0);
    let b = render(&model, &camera)?;
    save("02-orbit.png", &b, &mut report)?;
    report.check("orbit changes the view", differs(&a, &b) > 1000);

    // 3. More parts: an extruded cylinder, a revolved ring, a spline vase.
    let mut parts = demo_parts()?.into_iter().skip(1);
    let (n, f) = parts.next().ok_or("parts")?;
    model.add_feature(n, f);
    let c1 = render(&model, &camera)?;
    let (n, f) = parts.next().ok_or("parts")?;
    model.add_feature(n, f);
    let c2 = render(&model, &camera)?;
    let (n, f) = parts.next().ok_or("parts")?;
    model.add_feature(n, f);
    let c3 = render(&model, &camera)?;
    save("03-four-parts.png", &c3, &mut report)?;
    report.check(
        "all four features regenerate",
        model.bodies().bodies.len() == 4
            && model
                .bodies()
                .statuses
                .iter()
                .all(|(_, s)| *s == opendatum_features::FeatureStatus::Ok),
    );

    // 4. Undo/redo through the op log; images must match the earlier states
    // exactly.
    model.undo();
    let u1 = render(&model, &camera)?;
    model.undo();
    let u2 = render(&model, &camera)?;
    model.undo();
    let u3 = render(&model, &camera)?;
    save("04-undo-to-box.png", &u3, &mut report)?;
    report.check("undo 1 matches the 3-part image", differs(&u1, &c2) == 0);
    report.check("undo 2 matches the 2-part image", differs(&u2, &c1) == 0);
    report.check(
        "undo 3 matches the orbited box image",
        differs(&u3, &b) == 0,
    );
    model.redo();
    model.redo();
    model.redo();
    let r3 = render(&model, &camera)?;
    save("05-redo.png", &r3, &mut report)?;
    report.check("redo restores the 4-part image", differs(&r3, &c3) == 0);
    report.check(
        "undo/redo are ops in the log (nothing removed)",
        model.session().log().len() == 1 + 4 + 6,
    );

    // 5. Save the op log, reopen it, render the same image.
    let doc = out.join("demo.odatum");
    model.save(&doc)?;
    let reopened = Model::open(&doc)?;
    let ro = render(&reopened, &camera)?;
    report.check("save/open reproduces the model", differs(&ro, &c3) == 0);

    // 6. STL export.
    let stl = out.join("parts.stl");
    model.export_stl(&stl)?;
    let len = std::fs::metadata(&stl).map(|m| m.len()).unwrap_or(0);
    report.check("STL export written", len > 84);

    // 7. A fitted overview from the top-right.
    let mut overview = Camera::isometric();
    if let Some((lo, hi)) = model.bounds() {
        overview.fit(lo, hi);
    }
    overview.orbit(-60.0, 20.0);
    let ov = render(&model, &overview)?;
    save("06-overview.png", &ov, &mut report)?;
    Ok(report)
}

/// Frame-time measurements from [`bench()`].
#[derive(Debug, Default)]
pub struct BenchReport {
    /// Adapter description (name, backend, driver).
    pub adapter: String,
    /// Render target size in pixels.
    pub size: (u32, u32),
    /// Triangles drawn per frame.
    pub triangles: usize,
    /// Edge line segments drawn per frame.
    pub edge_segments: usize,
    /// Wall-clock time of each measured frame, in milliseconds, in order.
    pub frame_ms: Vec<f64>,
}

impl BenchReport {
    /// The `q` quantile (0..=1) of the frame times, in milliseconds.
    pub fn quantile(&self, q: f64) -> f64 {
        let mut v = self.frame_ms.clone();
        v.sort_by(f64::total_cmp);
        if v.is_empty() {
            return f64::NAN;
        }
        v[((v.len() - 1) as f64 * q.clamp(0.0, 1.0)).round() as usize]
    }

    /// Mean frame time in milliseconds.
    pub fn mean(&self) -> f64 {
        self.frame_ms.iter().sum::<f64>() / self.frame_ms.len() as f64
    }

    /// A plain-text summary.
    pub fn summary(&self) -> String {
        format!(
            "adapter: {}\n\
             scene: demo parts (4 bodies), {} triangles, {} edge segments, \
             {}x{}, 4x MSAA\n\
             frames: {} (orbiting a full turn), each timed to GPU completion, no readback\n\
             frame time ms: mean {:.3}, median {:.3}, p95 {:.3}, p99 {:.3}, max {:.3}\n\
             frames per second (from mean): {:.1}\n",
            self.adapter,
            self.triangles,
            self.edge_segments,
            self.size.0,
            self.size.1,
            self.frame_ms.len(),
            self.mean(),
            self.quantile(0.5),
            self.quantile(0.95),
            self.quantile(0.99),
            self.quantile(1.0),
            1000.0 / self.mean(),
        )
    }
}

/// Measures offscreen frame times for the demo scene: `frames` frames at
/// `width` x `height`, orbiting a full turn, each waited on until the GPU is
/// done. Twenty untimed warm-up frames come first.
pub fn bench(frames: u32, width: u32, height: u32) -> Result<BenchReport, String> {
    let mut gpu = OffscreenRenderer::new().map_err(|e| e.to_string())?;
    let info = gpu.adapter_info();
    let mut model = Model::new();
    for (name, f) in demo_parts()? {
        model.add_feature(name, f);
    }
    gpu.set_meshes(model.meshes()).map_err(|e| e.to_string())?;
    let stats = gpu.scene().stats();
    let (lo, hi) = model.bounds().ok_or("no geometry")?;
    let mut camera = Camera::isometric();
    camera.fit(lo, hi);
    for _ in 0..20 {
        gpu.render_frame(&camera, width, height)
            .map_err(|e| e.to_string())?;
    }
    let step = std::f64::consts::TAU / frames.max(1) as f64;
    let mut frame_ms = Vec::with_capacity(frames as usize);
    for _ in 0..frames {
        camera.orbit_radians(step, 0.0);
        let t = std::time::Instant::now();
        gpu.render_frame(&camera, width, height)
            .map_err(|e| e.to_string())?;
        frame_ms.push(t.elapsed().as_secs_f64() * 1e3);
    }
    Ok(BenchReport {
        adapter: format!("{} ({:?}, {})", info.name, info.backend, info.driver_info),
        size: (width, height),
        triangles: stats.triangles,
        edge_segments: stats.edge_segments,
        frame_ms,
    })
}
