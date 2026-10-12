# ADR 0007 — Rendering: wgpu, viewport-as-texture, offscreen-testable

Status: accepted (Phase 0)

* `opendatum-render` uses wgpu (version pinned to egui-wgpu's, now 30).
  Primary target is Vulkan on Mesa RADV (AMD); GL is the fallback backend;
  software Vulkan (lavapipe) is used for CI and headless tests.
* The renderer has no dependency on the kernel or egui: it consumes neutral
  `MeshData` and a `Camera`.
* The UI renders the viewport into an offscreen texture each frame
  (own depth buffer, 4× MSAA, resolve), registered with egui-wgpu as a native
  texture. The same `SceneRenderer` path is used by `OffscreenRenderer`, so
  tests and the headless demo exercise exactly the shipping code.
* Large-assembly performance (60 fps at 1,000+ parts) will come from
  instancing, LOD and GPU-side culling (Phase 6); the offscreen 1,000-instance
  timing test tracks it.
