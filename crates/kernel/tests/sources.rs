//! Checks on the kernel's own sources.

use std::path::{Path, PathBuf};

/// Every `.rs` file under `dir`, in a fixed order.
fn sources(dir: &Path, out: &mut Vec<PathBuf>) {
    let mut entries: Vec<PathBuf> = std::fs::read_dir(dir)
        .expect("a source directory")
        .map(|e| e.expect("an entry").path())
        .collect();
    entries.sort();
    for path in entries {
        if path.is_dir() {
            sources(&path, out);
        } else if path.extension().is_some_and(|x| x == "rs") {
            out.push(path);
        }
    }
}

/// The kernel's results depend on its inputs alone: no source of the
/// kernel or of its NURBS crate reads the process's environment. A
/// switch left from debugging once turned off the meshes' error estimate
/// when a variable was set (found by review, 2026-10-15).
#[test]
fn the_kernel_reads_no_environment() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut files = Vec::new();
    sources(&root.join("src"), &mut files);
    sources(&root.join("../nurbs/src"), &mut files);
    assert!(files.len() > 20, "found only {} sources", files.len());
    let found: Vec<String> = files
        .iter()
        .flat_map(|path| {
            let text = std::fs::read_to_string(path).expect("a source file");
            text.lines()
                .enumerate()
                .filter(|(_, line)| line.contains("env::var") || line.contains("var_os("))
                .map(|(i, line)| format!("{}:{}: {}", path.display(), i + 1, line.trim()))
                .collect::<Vec<_>>()
        })
        .collect();
    assert!(found.is_empty(), "{found:#?}");
}
