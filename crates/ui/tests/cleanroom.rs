//! Clean-room guard (docs/CLEANROOM.md): the name of the commercial product
//! OpenShape's feature goals are benchmarked against appears nowhere in the
//! repository except the three files that state the policy.

use std::path::{Path, PathBuf};

/// The name, assembled so this file does not contain it.
fn name() -> String {
    ["on", "sha", "pe"].concat()
}

/// Files allowed to name it, to state the boundary (relative to the root).
const ALLOWED: [&str; 3] = ["CLAUDE.md", "NOTICE", "docs/CLEANROOM.md"];

/// Build output and other generated, untracked folders.
const SKIPPED_DIRS: [&str; 4] = ["target", ".git", "out", "gpu-check-out"];

fn files(dir: &Path, out: &mut Vec<PathBuf>) {
    for entry in std::fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_string_lossy().to_string();
        if path.is_dir() {
            if !SKIPPED_DIRS.contains(&name.as_str()) {
                files(&path, out);
            }
        } else {
            out.push(path);
        }
    }
}

/// Byte offsets of whole-word, case-insensitive occurrences of `word` (so
/// `RegionShape` does not count).
fn occurrences(text: &str, word: &str) -> Vec<usize> {
    let lower = text.to_ascii_lowercase();
    let bytes = lower.as_bytes();
    let is_word = |b: u8| b.is_ascii_alphanumeric() || b == b'_';
    lower
        .match_indices(word)
        .map(|(i, _)| i)
        .filter(|&i| {
            let before = i.checked_sub(1).map(|j| bytes[j]);
            let after = bytes.get(i + word.len()).copied();
            !before.is_some_and(is_word) && !after.is_some_and(is_word)
        })
        .collect()
}

#[test]
fn the_reference_product_is_not_named_outside_the_policy_files() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = root.canonicalize().unwrap();
    let mut all = Vec::new();
    files(&root, &mut all);
    assert!(all.len() > 50, "found only {} files", all.len());
    let word = name();
    let mut found = Vec::new();
    for path in all {
        let rel = path
            .strip_prefix(&root)
            .unwrap()
            .to_string_lossy()
            .replace('\\', "/");
        if ALLOWED.contains(&rel.as_str()) {
            continue;
        }
        let bytes = std::fs::read(&path).unwrap();
        let text = String::from_utf8_lossy(&bytes);
        for i in occurrences(&text, &word) {
            let line = text[..i].matches('\n').count() + 1;
            found.push(format!("{rel}:{line}"));
        }
    }
    assert!(
        found.is_empty(),
        "docs/CLEANROOM.md: the reference product's name appears in:\n{}",
        found.join("\n")
    );
    // The policy files themselves do name it, so the scan can see it.
    let policy = std::fs::read_to_string(root.join("docs/CLEANROOM.md")).unwrap();
    assert!(!occurrences(&policy, &word).is_empty());
}

#[test]
fn whole_word_matching() {
    let w = name();
    assert_eq!(occurrences("RegionShape", &w), Vec::<usize>::new());
    assert_eq!(occurrences("OpenShape", &w), Vec::<usize>::new());
    assert_eq!(
        occurrences(&format!("an {} tab", w.to_uppercase()), &w),
        vec![3]
    );
    assert_eq!(occurrences(&format!("({w})"), &w), vec![1]);
}
