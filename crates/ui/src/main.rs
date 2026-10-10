//! OpenShape desktop application entry point.
//!
//! ```text
//! openshape                         # windowed app (Wayland or X11)
//! openshape <file.osdoc>            # open a document
//! openshape --headless-demo <dir>   # offscreen end-to-end demo, writes PNGs
//! openshape --bench <dir> [frames]  # offscreen frame-time benchmark
//! ```

use std::path::PathBuf;

use openshape_app::{app::OpenShapeApp, demo, Model};

fn main() -> std::process::ExitCode {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().map(String::as_str) == Some("--headless-demo") {
        let out = PathBuf::from(args.get(1).map(String::as_str).unwrap_or("out"));
        return match demo::run(&out) {
            Ok(report) => {
                println!("adapter: {}", report.adapter);
                for (name, ok) in &report.checks {
                    println!("[{}] {name}", if *ok { "PASS" } else { "FAIL" });
                }
                for p in &report.images {
                    println!("wrote {}", p.display());
                }
                if report.ok() {
                    std::process::ExitCode::SUCCESS
                } else {
                    std::process::ExitCode::FAILURE
                }
            }
            Err(e) => {
                eprintln!("headless demo failed: {e}");
                std::process::ExitCode::FAILURE
            }
        };
    }
    if args.first().map(String::as_str) == Some("--bench") {
        let out = PathBuf::from(args.get(1).map(String::as_str).unwrap_or("bench-out"));
        let frames = match args.get(2).map(|s| s.parse::<u32>()) {
            None => 600,
            Some(Ok(n)) if n > 0 => n,
            Some(_) => {
                eprintln!("--bench: frame count must be a positive integer");
                return std::process::ExitCode::FAILURE;
            }
        };
        return match demo::bench(frames, 1920, 1080) {
            Ok(report) => {
                let text = report.summary();
                print!("{text}");
                let csv: String = report
                    .frame_ms
                    .iter()
                    .enumerate()
                    .map(|(i, t)| format!("{i},{t:.4}\n"))
                    .collect();
                let written = std::fs::create_dir_all(&out)
                    .and_then(|_| std::fs::write(out.join("bench.txt"), &text))
                    .and_then(|_| {
                        std::fs::write(out.join("frame-times.csv"), format!("frame,ms\n{csv}"))
                    });
                match written {
                    Ok(()) => {
                        println!("wrote {}", out.join("bench.txt").display());
                        std::process::ExitCode::SUCCESS
                    }
                    Err(e) => {
                        eprintln!("bench: cannot write {}: {e}", out.display());
                        std::process::ExitCode::FAILURE
                    }
                }
            }
            Err(e) => {
                eprintln!("bench failed: {e}");
                std::process::ExitCode::FAILURE
            }
        };
    }
    if args.first().map(String::as_str) == Some("--version") {
        println!("openshape {}", env!("CARGO_PKG_VERSION"));
        return std::process::ExitCode::SUCCESS;
    }
    let model = match args.first() {
        Some(path) => match Model::open(std::path::Path::new(path)) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("cannot open {path}: {e}");
                return std::process::ExitCode::FAILURE;
            }
        },
        None => Model::new(),
    };
    let options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_title("OpenShape")
            .with_inner_size([1280.0, 800.0]),
        ..Default::default()
    };
    match eframe::run_native(
        "OpenShape",
        options,
        Box::new(|cc| Ok(Box::new(OpenShapeApp::new(cc, model)))),
    ) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("openshape: {e}");
            std::process::ExitCode::FAILURE
        }
    }
}
