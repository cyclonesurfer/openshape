# Task prompts for Claude Code on a GPU machine

These prompts are for Claude Code running directly on a machine with a real
GPU. Today that machine is **huginn** (Fedora 44, Radeon 8060S, RADV), with
the clone at `~/openshape`. They are the agent-side counterpart of
[RUNNING.md](RUNNING.md). Start Claude Code in the clone (`cd ~/openshape &&
claude`) and paste one prompt. Replace anything in `<angle brackets>`.

Each prompt makes the agent read the ground rules below first. That keeps the
prompts short, and the rules live in one place.

## Ground rules

Every prompt below sends the agent here first. Follow these on every task.

1. **The project's rules are in CLAUDE.md.** The ones that matter most on
   the box:
   - Never make a failing check pass by weakening it. Do not skip, ignore,
     delete or loosen a test, add `#[ignore]`, or raise a tolerance.
   - Sometimes a tolerance really is the right answer, because a test
     assumes something the graphics APIs don't guarantee. That is Pat's
     decision. Name it in the report as a judgment call, with the evidence
     (spec text, measured numbers), and don't apply it yourself unless asked.
   - Report exactly what happened: the commands run, their exit status, and
     numbers copied from the output. A failure is never summarized as a
     pass. Say what you did not run.
   - Stop and escalate on a correctness problem you can't fix properly. Tell
     Pat the specifics and the options instead of working around it.
   - Clean room: OpenShape is an independent implementation. Don't consult
     or copy code, shaders, assets or file formats from any commercial CAD
     product.
   - Before any commit, all four checks must pass: `cargo fmt --all`,
     `cargo clippy --workspace --all-targets -- -D warnings`,
     `cargo test --workspace` and `cargo deny check`.
2. **Stay inside the clone.**
   - Don't change drivers, system settings or files outside `~/openshape`
     (scratch output under `/tmp` is fine).
   - If a tool is missing, install only what RUNNING.md §1 lists for Fedora
     (`sudo dnf install …`). List every package you installed in the report.
3. **Git.**
   - Never push to `main` or rewrite published history.
   - Validation tasks change no code and make no commits.
   - When a task allows fixes, work on a branch named `gpu/<topic>` and push
     only that branch. The cloud project reviews and merges it.
4. **Use the real GPU.**
   - Export `OPENSHAPE_WGPU_BACKEND=vulkan` for cargo commands
     (`scripts/gpu-check.sh` sets it itself).
   - The `adapter:` lines must name the AMD card (RADV). If one says
     llvmpipe or lavapipe, the run used the CPU renderer: stop and report
     that first.
5. **Report back.** Pat pastes your report into the OpenShape project, so
   finish with one fenced block of at most 60 lines. It contains:
   - the commit (`git rev-parse --short HEAD`) and the adapter line;
   - PASS/FAIL per step, with test totals;
   - each failure's assertion message and `file:line`;
   - the numbers the task asks for;
   - any packages installed;
   - your reading of each failure (cause, evidence, proposed fix), labelled
     as your assessment.

## 1. Standard GPU check

Use this after anything lands on `main` that touches rendering, or whenever a
fresh real-hardware result is wanted.

```text
In ~/openshape, run `git switch main && git pull --ff-only`, then read
CLAUDE.md and the "Ground rules" section of docs/huginn-tasks.md and follow
them. Run scripts/gpu-check.sh. It builds, runs every test, the headless demo,
the 1080p benchmark and the 1,000-instance timing on the real GPU, and writes
gpu-check-out/report.txt. Read the report and the log of every failed step in
gpu-check-out/. Change no code.
Report: the commit, the adapter, the RESULT line, PASS/FAIL per step, test
totals, the demo checks, the benchmark mean/p95/p99, and for each failure the
assertion, file:line and your assessment of the cause.
```

## 2. Regression check of a branch or pull request

Use this before merging rendering or GPU-facing changes.

```text
In ~/openshape, run `git switch main && git pull --ff-only`, read CLAUDE.md
and the "Ground rules" section of docs/huginn-tasks.md, and follow them. Then
check out <branch name, or PR number via `gh pr checkout <n>`>. With
OPENSHAPE_WGPU_BACKEND=vulkan, run:
  cargo build --release --workspace --all-targets
  cargo test --release -p openshape-render -- --nocapture
  cargo test --release --workspace --no-fail-fast
  cargo run --release -p openshape -- --headless-demo /tmp/os-demo
For every test that fails on the branch, run the same test on main and say
whether the failure is new. Change no code.
Report per command: exit status, totals, failures with assertion text and
file:line, the adapter line, and new versus already failing on main.
```

## 3. Performance run

Use this to track frame times against the roadmap goal of 60 fps (16.7 ms) at
1,000+ parts on RADV.

```text
In ~/openshape, run `git switch main && git pull --ff-only`, read CLAUDE.md
and the "Ground rules" section of docs/huginn-tasks.md, and follow them. Close
other GPU-heavy programs. With OPENSHAPE_WGPU_BACKEND=vulkan, run each of
these three times:
  cargo run --release -p openshape -- --bench /tmp/os-bench 600
  cargo test --release -p openshape-render --test offscreen thousand_instances_timing -- --nocapture
Change no code.
Report for each run: the adapter, the scene (triangles, edge segments,
resolution), the frame time mean/median/p95/p99/max and fps from bench.txt,
and the 1,000-instance upload, first-frame and steady times. Then give the
median of the three and compare it with the previous numbers in STATUS.md and
with the 16.7 ms target. The benchmark scene is the 4-part demo. Only the
1,000-instance timing speaks to the 1,000-part goal, so don't claim more than
it shows. Note anything that could skew the numbers: power profile, other
load, thermal throttling.
```

## 4. Investigate a failure that only happens on the GPU

Use this only when Pat asks for a fix on the box. It is the one prompt that
allows code changes, on a `gpu/<topic>` branch.

```text
In ~/openshape, run `git switch main && git pull --ff-only`, read CLAUDE.md
and the "Ground rules" section of docs/huginn-tasks.md, and follow them.
<test name> fails on this GPU but passes on CI's software renderer
(lavapipe). Reproduce it with OPENSHAPE_WGPU_BACKEND=vulkan. Set
OPENSHAPE_RENDER_TEST_OUT=/tmp/os-render to save the test images, then find
the root cause. Decide which case it is:
(a) a bug in OpenShape: fix the code;
(b) a driver bug: write a minimal reproduction and report it, without
    working around it silently;
(c) a test that assumes something the graphics APIs don't guarantee: stop
    and give Pat the evidence and a proposed documented bound to decide on.
    Don't apply the bound.
Work on a branch gpu/<topic>. Before committing, run fmt, clippy with
-D warnings, the full test suite and cargo deny; push only that branch.
Report: the cause and its evidence, the case (a/b/c), what you changed and
the branch name, the test results on this GPU, and any open questions.
```
