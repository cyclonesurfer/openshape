// OpenShape viewport shader: gradient background, shaded meshes, edge lines.
//
// Lighting is computed in linear space. All colours in the uniform block are
// already linear; mesh vertex colours arrive sRGB-encoded (Unorm8x4) and are
// decoded in the vertex stage. If the colour target is not an sRGB format the
// fragment stages encode to sRGB themselves (params.z > 0.5).
//
// Depth is reverse-Z (near = 1, far = 0, compare Greater).

struct Globals {
    view: mat4x4<f32>,
    proj: mat4x4<f32>,
    // xyz: unit direction towards the key light in view space; w: intensity.
    key_light: vec4<f32>,
    // Hemisphere ambient: colour for normals pointing to world +Z / -Z.
    sky: vec4<f32>,
    ground: vec4<f32>,
    edge_color: vec4<f32>,
    bg_top: vec4<f32>,
    bg_bottom: vec4<f32>,
    // x: edge view-space scale, y: edge view-space z offset,
    // z: 1 if the shader must sRGB-encode its output, w: 1 if perspective.
    params: vec4<f32>,
    // x: specular strength, y: shininess.
    material: vec4<f32>,
    // xy: target size in pixels, z: half the edge line width in pixels.
    viewport: vec4<f32>,
};

@group(0) @binding(0) var<uniform> g: Globals;

fn srgb_to_linear(c: vec3<f32>) -> vec3<f32> {
    let lo = c / 12.92;
    let hi = pow((c + vec3<f32>(0.055)) / 1.055, vec3<f32>(2.4));
    return select(hi, lo, c <= vec3<f32>(0.04045));
}

fn linear_to_srgb(c: vec3<f32>) -> vec3<f32> {
    let x = clamp(c, vec3<f32>(0.0), vec3<f32>(1.0));
    let lo = x * 12.92;
    let hi = 1.055 * pow(x, vec3<f32>(1.0 / 2.4)) - vec3<f32>(0.055);
    return select(hi, lo, x <= vec3<f32>(0.0031308));
}

fn encode_output(c: vec3<f32>) -> vec3<f32> {
    if g.params.z > 0.5 {
        return linear_to_srgb(c);
    }
    return c;
}

// ---------------------------------------------------------------------------
// Background: one full-screen triangle with a vertical gradient.

struct BgOut {
    @builtin(position) clip: vec4<f32>,
    // 1 at the top edge of the viewport, 0 at the bottom.
    @location(0) t: f32,
};

@vertex
fn vs_background(@builtin(vertex_index) i: u32) -> BgOut {
    let x = f32(i32(i & 1u) * 4 - 1);
    let y = f32(i32(i >> 1u) * 4 - 1);
    var out: BgOut;
    out.clip = vec4<f32>(x, y, 0.0, 1.0);
    out.t = y * 0.5 + 0.5;
    return out;
}

@fragment
fn fs_background(in: BgOut) -> @location(0) vec4<f32> {
    let c = mix(g.bg_bottom.rgb, g.bg_top.rgb, clamp(in.t, 0.0, 1.0));
    return vec4<f32>(encode_output(c), 1.0);
}

// ---------------------------------------------------------------------------
// Shaded meshes.

struct MeshIn {
    @location(0) position: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) color: vec4<f32>,
};

struct MeshOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) normal_world: vec3<f32>,
    @location(1) pos_view: vec3<f32>,
    @location(2) color: vec4<f32>,
};

@vertex
fn vs_mesh(v: MeshIn) -> MeshOut {
    let pv = g.view * vec4<f32>(v.position, 1.0);
    var out: MeshOut;
    out.clip = g.proj * pv;
    out.normal_world = v.normal;
    out.pos_view = pv.xyz;
    out.color = vec4<f32>(srgb_to_linear(v.color.rgb), v.color.a);
    return out;
}

@fragment
fn fs_mesh(in: MeshOut) -> @location(0) vec4<f32> {
    // Direction towards the eye in view space.
    var v = vec3<f32>(0.0, 0.0, 1.0);
    if g.params.w > 0.5 {
        v = normalize(-in.pos_view);
    }
    var nw = in.normal_world;
    let len = length(nw);
    var n = v;
    if len > 1e-20 {
        nw = nw / len;
        n = normalize((g.view * vec4<f32>(nw, 0.0)).xyz);
    } else {
        nw = vec3<f32>(0.0, 0.0, 1.0);
    }
    // Two-sided lighting: back faces of open shells are lit like front faces.
    if dot(n, v) < 0.0 {
        n = -n;
        nw = -nw;
    }
    let hemi = 0.5 + 0.5 * nw.z;
    let ambient = mix(g.ground.rgb, g.sky.rgb, hemi);
    let l = g.key_light.xyz;
    let ndl = max(dot(n, l), 0.0);
    let h = normalize(l + v);
    let spec = select(0.0, pow(max(dot(n, h), 0.0), g.material.y), ndl > 0.0) * g.material.x;
    let base = in.color.rgb;
    let rgb = base * (ambient + vec3<f32>(g.key_light.w * ndl)) + vec3<f32>(spec);
    return vec4<f32>(encode_output(rgb), in.color.a);
}

// ---------------------------------------------------------------------------
// Edge lines, drawn as screen-space quads of a fixed pixel width (one
// instance per segment, six vertices) so they look the same on every driver.
// End points are pulled slightly towards the eye so lines win the depth test
// against the faces they lie on. In perspective the view-space position is
// scaled towards the eye (a relative depth bias that keeps the screen
// position); in orthographic it is shifted along +Z by a constant.

fn edge_clip(position: vec3<f32>) -> vec4<f32> {
    let pv = g.view * vec4<f32>(position, 1.0);
    let biased = vec3<f32>(pv.xyz * g.params.x) + vec3<f32>(0.0, 0.0, g.params.y);
    return g.proj * vec4<f32>(biased, 1.0);
}

@vertex
fn vs_edge(
    @builtin(vertex_index) vi: u32,
    @location(0) a: vec3<f32>,
    @location(1) b: vec3<f32>,
) -> @builtin(position) vec4<f32> {
    var ca = edge_clip(a);
    var cb = edge_clip(b);
    // Keep both ends in front of the eye before dividing by w; the hardware
    // near-plane clip handles the rest.
    let min_w = 1e-6;
    if ca.w < min_w && cb.w < min_w {
        return vec4<f32>(2.0, 2.0, 2.0, 1.0);
    }
    if ca.w < min_w {
        ca = mix(ca, cb, (min_w - ca.w) / (cb.w - ca.w));
    } else if cb.w < min_w {
        cb = mix(cb, ca, (min_w - cb.w) / (ca.w - cb.w));
    }
    // Corners: (end, side) = (a,-) (a,+) (b,-) | (b,-) (a,+) (b,+).
    let at_b = vi == 2u || vi == 3u || vi == 5u;
    let side = select(-1.0, 1.0, vi == 1u || vi == 4u || vi == 5u);
    let size = g.viewport.xy;
    let pa = ca.xy / ca.w * 0.5 * size;
    let pb = cb.xy / cb.w * 0.5 * size;
    var dir = vec2<f32>(1.0, 0.0);
    let len = length(pb - pa);
    if len > 1e-6 {
        dir = (pb - pa) / len;
    }
    let normal = vec2<f32>(-dir.y, dir.x);
    let half_width = g.viewport.z;
    // Square caps (extend by half the width) so polyline joins have no gaps.
    let along = select(-1.0, 1.0, at_b);
    let offset_px = normal * (side * half_width) + dir * (along * half_width);
    let c = select(ca, cb, at_b);
    return vec4<f32>(c.xy + offset_px * 2.0 / size * c.w, c.z, c.w);
}

@fragment
fn fs_edge() -> @location(0) vec4<f32> {
    return vec4<f32>(encode_output(g.edge_color.rgb), g.edge_color.a);
}
