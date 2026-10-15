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

/// The places the modules that `file` declares out of line under
/// `#[cfg(test)]` may live: a file, or a directory all of whose sources
/// belong to the module. Built only for tests, they cannot change the
/// kernel's results.
///
/// The reading is line by line and errs towards scanning: an attribute
/// counts only written whole on its own line as `#[cfg(test)]`, and any
/// item between it and the declaration ends it.
fn test_only_modules(file: &Path, text: &str) -> Vec<PathBuf> {
    let dir = file.parent().expect("a source's directory");
    let stem = file.file_stem().and_then(|s| s.to_str()).unwrap_or("");
    // Modules declared in `mod.rs`, `lib.rs` and `main.rs` live beside
    // them; those declared in `name.rs` live under `name/`.
    let base = if matches!(stem, "mod" | "lib" | "main") {
        dir.to_path_buf()
    } else {
        dir.join(stem)
    };
    let mut found = Vec::new();
    let mut attributes: Vec<String> = Vec::new();
    for line in text.lines() {
        let line: String = line.split_whitespace().collect::<Vec<_>>().join(" ");
        if line.is_empty() || line.starts_with("//") {
            continue;
        }
        if line.starts_with("#[") {
            attributes.push(line);
            continue;
        }
        let declared = line
            .strip_prefix("pub(crate) ")
            .or_else(|| line.strip_prefix("pub(super) "))
            .or_else(|| line.strip_prefix("pub "))
            .unwrap_or(&line)
            .strip_prefix("mod ")
            .and_then(|rest| rest.strip_suffix(';'));
        if let Some(name) = declared {
            if attributes.iter().any(|a| a == "#[cfg(test)]") {
                // `#[path]` on an out-of-line module names a file relative
                // to the declaring file's directory.
                let path = attributes.iter().find_map(|a| {
                    a.strip_prefix("#[path = \"")
                        .and_then(|rest| rest.strip_suffix("\"]"))
                });
                match path {
                    Some(path) => found.push(dir.join(path)),
                    None => {
                        found.push(base.join(format!("{name}.rs")));
                        found.push(base.join(name));
                    }
                }
            }
        }
        attributes.clear();
    }
    found
}

/// The kernel's results depend on its inputs alone: no source of the
/// kernel or of its NURBS crate reads the process's environment, when it
/// runs or when it is built (`env!`, `option_env!`). A switch left from
/// debugging once turned off the meshes' error estimate when a variable
/// was set (found by review, 2026-10-15).
///
/// Modules built only for tests are left out: the rounding tests read
/// their sample sizes from the environment, as the integration tests'
/// fuzzes do.
#[test]
fn the_kernel_reads_no_environment() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut files = Vec::new();
    sources(&root.join("src"), &mut files);
    sources(&root.join("../nurbs/src"), &mut files);
    assert!(files.len() > 20, "found only {} sources", files.len());
    let texts: Vec<String> = files
        .iter()
        .map(|path| std::fs::read_to_string(path).expect("a source file"))
        .collect();
    let excluded: Vec<PathBuf> = files
        .iter()
        .zip(&texts)
        .flat_map(|(path, text)| test_only_modules(path, text))
        .collect();
    let left_out = |path: &Path| excluded.iter().any(|e| path.starts_with(e));
    // The parser finds the test-only modules there are.
    for name in [
        "boolean/dd.rs",
        "boolean/testing.rs",
        "boolean/chart_rounding.rs",
    ] {
        assert!(
            left_out(&root.join("src").join(name)),
            "{name} is not left out"
        );
    }
    let found: Vec<String> = files
        .iter()
        .zip(&texts)
        .filter(|(path, _)| !left_out(path))
        .flat_map(|(path, text)| {
            text.lines()
                .enumerate()
                .filter(|(_, line)| {
                    ["std::env", "env::", "env!(", "var_os("]
                        .iter()
                        .any(|p| line.contains(p))
                })
                .map(|(i, line)| format!("{}:{}: {}", path.display(), i + 1, line.trim()))
                .collect::<Vec<_>>()
        })
        .collect();
    assert!(found.is_empty(), "{found:#?}");
}

/// The reading of `#[cfg(test)]` declarations takes test-only modules,
/// with and without `#[path]`, in `mod.rs` and in `name.rs` files, and
/// nothing else.
#[test]
fn test_only_modules_are_read_from_their_declarations() {
    let text = "\
//! A module.
mod plain;
#[cfg(test)]
mod dd;
/// Helpers.
#[cfg(test)]
pub(crate) mod testing;
#[cfg(not(test))]
mod production;
#[cfg(test)]
fn helper() {}
mod after_an_item;
#[cfg(test)]
#[path = \"chart_rounding.rs\"]
mod rounding;
";
    let at = |file: &str| test_only_modules(Path::new(file), text);
    let dir = |names: &[&str]| -> Vec<PathBuf> { names.iter().map(PathBuf::from).collect() };
    assert_eq!(
        at("src/boolean/mod.rs"),
        dir(&[
            "src/boolean/dd.rs",
            "src/boolean/dd",
            "src/boolean/testing.rs",
            "src/boolean/testing",
            "src/boolean/chart_rounding.rs",
        ])
    );
    assert_eq!(
        at("src/boolean/chart.rs"),
        dir(&[
            "src/boolean/chart/dd.rs",
            "src/boolean/chart/dd",
            "src/boolean/chart/testing.rs",
            "src/boolean/chart/testing",
            "src/boolean/chart_rounding.rs",
        ])
    );
}
