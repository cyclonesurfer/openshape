# Running OpenShape's checks on a real GPU

This builds OpenShape on a Linux machine with an AMD GPU (RADV), runs the test
suite, renders the headless demo on the GPU and measures frame times. It works
on a headless server: nothing needs a display.

Everything below runs on the GPU machine itself, over your own SSH session.

## 1. Prerequisites (once)

Pick your distribution. As root, or prefix with `sudo`.

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

The X11/Wayland libraries are only needed to build the desktop app; the checks
themselves render offscreen.

Rust (as the user who will build; the repository pins the stable channel):

```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
. "$HOME/.cargo/env"
```

Confirm the GPU is visible to Vulkan. The output should name your AMD card with
driver `radv`, not `llvmpipe`:

```sh
vulkaninfo --summary | grep -E 'deviceName|driverName'
```

If it only lists `llvmpipe`, RADV isn't installed or the user can't open the
GPU's render node: run as root, or add the user to the `render` group
(`usermod -aG render $USER`, then log in again).

## 2. Run the checks

```sh
git clone https://github.com/cyclonesurfer/openshape.git
cd openshape
scripts/gpu-check.sh
```

To re-run later: `cd openshape && git pull && scripts/gpu-check.sh`.

The first run compiles everything in release mode, which takes a while.

## 3. Reading the result

The script prints one line per step and ends with `RESULT:`.

| step | passes when |
|---|---|
| `build` | the workspace compiles in release mode |
| `tests` | every unit, golden and render test passes |
| `demo` | all headless-demo checks print `[PASS]` (box built through the op log, orbit, four parts, undo/redo pixel-exact, save/open, STL) |
| `bench` | the frame-time benchmark completes |

The `adapter:` line must name the AMD card (RADV). If the script prints
`WARNING: rendering ran on Mesa's software driver`, the run did not use the
GPU; see the note in step 1.

The benchmark renders the demo scene (four parts, about 37k triangles) at
1920x1080 with 4x MSAA for 600 frames while orbiting, timing each frame until
the GPU finishes. For reference, Mesa's software renderer in a cloud
container gives a mean of about 63 ms (16 fps); a discrete AMD card should be
well under 16.7 ms (60 fps).

## 4. What to send back

- `gpu-check-out/report.txt`: paste its contents.
- `gpu-check-out/demo/*.png`: the six rendered images. To copy them to your
  own computer, run this there:

  ```sh
  scp -r root@10.0.20.10:openshape/gpu-check-out/demo ./openshape-gpu-demo
  ```

  (Adjust the user and path if you cloned somewhere other than the home
  directory.) `gpu-check-out.tar.gz` contains the report, logs, images and
  per-frame times (`bench/frame-times.csv`) in one file.

If a step fails, its last 30 log lines are in the report; the full log is
`gpu-check-out/<step>.log`.
