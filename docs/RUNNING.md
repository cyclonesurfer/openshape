# Running OpenDatum

How to build and run OpenDatum on Linux: the desktop app on an AMD/RADV
machine, the developer checks, the offscreen demo and benchmark, and a
one-command check for a GPU machine with no display.

## 1. Prerequisites (once per machine)

Pick your distribution. Run as root, or prefix with `sudo`.

Debian / Ubuntu:

```sh
apt-get update
apt-get install -y build-essential pkg-config git curl \
  mesa-vulkan-drivers libvulkan1 vulkan-tools \
  libxkbcommon-dev libwayland-dev libx11-dev libxcursor-dev libxrandr-dev libxi-dev
```

Fedora:

```sh
dnf install -y gcc gcc-c++ make pkgconf-pkg-config git curl \
  mesa-vulkan-drivers vulkan-loader vulkan-tools \
  libxkbcommon-devel wayland-devel libX11-devel libXcursor-devel libXrandr-devel libXi-devel
```

Arch:

```sh
pacman -S --needed base-devel git curl vulkan-radeon vulkan-icd-loader vulkan-tools \
  libxkbcommon wayland libx11 libxcursor libxrandr libxi
```

Rust, as the user who will build (the repository pins the stable channel):

```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
. "$HOME/.cargo/env"
```

A distribution's packaged Rust (Fedora's `rust` and `cargo`, say) also
builds OpenDatum, but it ignores `rust-toolchain.toml` and may be older
than CI's stable; `scripts/gpu-check.sh` prints the `rustc` it used. The
suite passes on Rust 1.95 (Fedora 44's) as well as on current stable.

Confirm Vulkan sees the GPU. The output should name your AMD card with driver
`radv`, not `llvmpipe`:

```sh
vulkaninfo --summary | grep -E 'deviceName|driverName'
```

If it lists only `llvmpipe`, RADV isn't installed or your user can't open the
GPU's render node: add the user to the `render` group
(`usermod -aG render $USER`, then log in again) or run as root.

## 2. Get the code

```sh
git clone https://github.com/cyclonesurfer/openshape.git
cd openshape
```

Later: `git pull` inside the clone.

## 3. Desktop app (AMD/RADV workstation)

```sh
cargo run --release -p opendatum                 # empty document
cargo run --release -p opendatum -- part.odatum   # open a saved document
```

The first build compiles everything in release mode and takes a while; later
runs start immediately. The app runs natively on Wayland or X11 and uses
Vulkan by default. `OPENDATUM_WGPU_BACKEND=gl` forces OpenGL if Vulkan
misbehaves (`vulkan` and `any` are the other values; wgpu's own
`WGPU_BACKEND` is honoured when it is unset). The status bar at the bottom
right names the adapter and backend in use, and `RUST_LOG=opendatum_app=info`
prints them with the driver in the terminal.

Controls: the toolbar adds sketches (S), extrudes (E), revolves (V), a box
and templates (cylinder, plate with hole, ring, spline vase); Ctrl+K opens
the command search. Right-drag orbits (left-drag too, outside sketch mode),
middle-drag or Ctrl+right-drag pans, the scroll wheel zooms, F fits the
view, and the view cube turns it. Ctrl+Z undoes, Ctrl+Shift+Z or Ctrl+Y
redoes. Right-click a feature in the list to edit, rename, suppress, move or
delete it; drag the rollback bar to roll the history back. The menu at the
top left saves and opens `.odatum` documents and exports STL.
`cargo run -p opendatum -- --screenshots out/` writes screenshots of the
whole app rendered offscreen.

The windowed app has so far been tested only offscreen in CI; please report
anything that looks wrong on real hardware.

## 4. Developer checks

These four must pass before every commit (see CLAUDE.md):

```sh
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo deny check        # first: cargo install cargo-deny --locked
```

## 5. Offscreen demo and benchmark

Both render offscreen and need no display.

```sh
cargo run --release -p opendatum -- --headless-demo out/demo
cargo run --release -p opendatum -- --bench out/bench 600
```

The demo prints the adapter, then one `[PASS]` or `[FAIL]` line per check, and
exits non-zero if any check fails. It builds a box through the op log, orbits,
adds a cylinder, a 270° ring and a spline vase, undoes and redoes (each image
must match the earlier state pixel for pixel), saves and reopens the
document, and exports STL. Output in `out/demo/`: `01-box.png`,
`02-orbit.png`, `03-four-parts.png`, `04-undo-to-box.png`, `05-redo.png`,
`06-overview.png`, `demo.odatum`, `parts.stl`.

The benchmark renders the demo scene (four parts, 36,656 triangles) at
1920x1080 with 4x MSAA for the given number of frames while orbiting, timing
each frame until the GPU finishes. It prints mean, median, p95, p99 and max
frame times and writes `out/bench/bench.txt` and `out/bench/frame-times.csv`.
For reference, Mesa's software renderer in a cloud container averages about
67 ms per frame (15 fps); a discrete AMD card should be well under 16.7 ms
(60 fps).

## 6. Checking a GPU machine in one step (works on a headless server)

After steps 1 and 2:

```sh
scripts/gpu-check.sh
```

It builds in release mode, runs the test suite, the headless demo, the
benchmark and the renderer's 1,000-instance timing test, all on Vulkan, and
prints one line per step:

| step | passes when |
|---|---|
| `build` | the workspace compiles in release mode |
| `tests` | every unit, golden and render test passes |
| `demo` | every headless-demo check prints `[PASS]` |
| `bench` | the frame-time benchmark completes |
| `render-timing` | the 1,000-instance render test passes |

It ends with `RESULT: all steps passed` or lists what failed, with the last
30 lines of each failing log. If it prints `WARNING: rendering ran on Mesa's
software driver`, the run did not use the GPU; see the note in step 1.

Everything lands in `gpu-check-out/`: `report.txt` (the summary to send
back), the demo images in `demo/`, the benchmark in `bench/` and a log per
step. `gpu-check-out.tar.gz`, next to that folder, holds all of it, the
finished report included. To copy the images from a server to
your own computer, run this on your computer (adjust user, host and path):

```sh
scp -r root@10.0.20.10:opendatum/gpu-check-out/demo ./opendatum-gpu-demo
```

## 7. With Claude Code on the machine

To work hands-on on a GPU machine, run Claude Code inside the clone
(`cd openshape && claude`). CLAUDE.md gives it the project's rules and
commands, and it works against the same repository as the cloud project.

[huginn-tasks.md](huginn-tasks.md) has ready-to-paste task prompts for it:
the standard GPU check, a regression check of a branch or pull request, a
performance run, and investigating a failure that only happens on the GPU.
Each prompt holds the agent to the project's rules and ends with a report to
paste back into the project.
