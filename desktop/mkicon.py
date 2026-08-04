#!/usr/bin/env python3
"""Render the FluLens mark to a PNG with no third-party libraries.

There is no ImageMagick, PIL, or node on this machine, and the browser pane's
downloads are sandboxed away from the filesystem. A PNG is a zlib-compressed
pixel buffer with four chunks, so writing one directly is less work than
installing a toolchain for it.

Geometry matches the 24-unit viewBox the header SVG uses, scaled up, so the app
icon and the in-app mark are the same drawing.
"""
import math, struct, zlib, sys

S = 1024                      # macOS wants 1024 for the source icon
K = S / 24.0                  # the mark is authored on a 24-unit grid

BG   = (0x14, 0x17, 0x1d, 255)
LENS = (0x4d, 0xab, 0xf7, 255)
CORE = (0x8b, 0x5c, 0xf6, 255)
KNOB = (0xa7, 0x8b, 0xfa, 255)
# The glint is a pale tint of the lens blue, not white: it has to read as light
# ON this glass rather than as a scratch through it, and it sits directly beside
# the rim, where pure white would be a second colour in a two-colour mark.
GLINT = (0xd0, 0xeb, 0xff, 255)

def blend(dst, src, a):
    """src over dst at coverage a — antialiasing, so edges are not stepped."""
    return tuple(int(round(d + (s - d) * a)) for d, s in zip(dst[:3], src[:3])) + (255,)

def dist_seg(px, py, x1, y1, x2, y2):
    dx, dy = x2 - x1, y2 - y1
    L2 = dx * dx + dy * dy
    t = 0.0 if L2 == 0 else max(0.0, min(1.0, ((px - x1) * dx + (py - y1) * dy) / L2))
    return math.hypot(px - (x1 + t * dx), py - (y1 + t * dy))

def dist_arc(px, py, cx, cy, R, lo, hi):
    """Distance to a circular arc — the SVG `A` command with round caps.

    Inside the angular sweep the nearest point is radial, so the distance is
    |r - R|; outside it, the nearest point is an endpoint, which is exactly what
    stroke-linecap="round" draws. Thresholding this at a constant gives a stroke
    of uniform width, unlike a radial falloff.
    """
    dx, dy = px - cx, py - cy
    th = math.atan2(dy, dx)                            # y grows downward
    while th < lo - math.pi: th += 2 * math.pi
    while th > lo + math.pi: th -= 2 * math.pi
    if lo <= th <= hi:
        return abs(math.hypot(dx, dy) - R)
    return min(math.hypot(px - (cx + math.cos(a) * R), py - (cy + math.sin(a) * R))
               for a in (lo, hi))

# Rounded-square tile: a superellipse is what macOS icons read as, and a plain
# square looks wrong beside every other icon in the Dock.
def tile_alpha(px, py):
    r = S * 0.22
    x = min(px, S - px); y = min(py, S - py)
    if x >= r or y >= r:
        return 1.0
    d = math.hypot(r - x, r - y)
    return max(0.0, min(1.0, r - d + 0.5))

knobs = [(10.5 + math.cos(i * 2 * math.pi / 9) * 4.9,
          10.5 + math.sin(i * 2 * math.pi / 9) * 4.9) for i in range(9)]

rows = []
for py in range(S):
    row = bytearray()
    uy = py / K
    for px in range(S):
        ux = px / K
        ta = tile_alpha(px, py)
        if ta <= 0:
            row += bytes((0, 0, 0, 0)); continue
        c = BG

        # lens ring, stroke width 2 units centred on r=7.5
        d = abs(math.hypot(ux - 10.5, uy - 10.5) - 7.5)
        a = max(0.0, min(1.0, (1.0 - d) * K))          # half-width 1.0 unit
        if a > 0: c = blend(c, LENS, a)

        # handle, stroke width 2.6, round caps
        d = dist_seg(ux, uy, 16, 16, 21.5, 21.5)
        a = max(0.0, min(1.0, (1.3 - d) * K))
        if a > 0: c = blend(c, LENS, a)

        # Glass reflection — the F7 "glass highlight" from the logo sheet, drawn
        # at its own geometry: an arc through the upper left, where a light
        # source above and left catches a curved lens.
        #
        # It was first drawn as a radial gradient fading out at both ends, on
        # the theory that a soft glint is more lenslike than a painted stripe.
        # It was INVISIBLE at Dock size — 204/255 at the single brightest pixel
        # and 71 by the ends, which averages to grey mush the moment the icon is
        # resampled to 32 px. A highlight has to be a constant-width,
        # constant-alpha stroke to survive downscaling, for the same reason the
        # ring and the handle are strokes. Match the SVG exactly: uniform fill,
        # round caps, antialiased only at its edge.
        #
        # Drawn after the rim but BEFORE the virion, so the specimen sits on top
        # of the glass rather than under the highlight — the way round a real
        # lens works.
        d = dist_arc(ux, uy, 10.5, 10.5, 6.2, -2.53, -1.73)
        a = max(0.0, min(1.0, (0.75 - d) * K))         # stroke width 1.5 units
        if a > 0: c = blend(c, GLINT, a)

        # virion core
        d = math.hypot(ux - 10.5, uy - 10.5)
        a = max(0.0, min(1.0, (3.4 - d) * K))
        if a > 0: c = blend(c, CORE, a)

        # nine knobs
        for kx, ky in knobs:
            d = math.hypot(ux - kx, uy - ky)
            a = max(0.0, min(1.0, (1.05 - d) * K))
            if a > 0:
                c = blend(c, KNOB, a); break

        row += bytes((c[0], c[1], c[2], int(round(255 * ta))))
    rows.append(row)

raw = b"".join(b"\x00" + bytes(r) for r in rows)       # filter byte 0 per scanline

def chunk(tag, data):
    return (struct.pack(">I", len(data)) + tag + data
            + struct.pack(">I", zlib.crc32(tag + data) & 0xffffffff))

png = (b"\x89PNG\r\n\x1a\n"
       + chunk(b"IHDR", struct.pack(">IIBBBBB", S, S, 8, 6, 0, 0, 0))
       + chunk(b"IDAT", zlib.compress(raw, 9))
       + chunk(b"IEND", b""))

open(sys.argv[1], "wb").write(png)
print(f"wrote {sys.argv[1]}  {len(png)} bytes  {S}x{S}")
