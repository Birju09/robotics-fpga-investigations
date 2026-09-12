#!/usr/bin/env python3
"""
Render the manipulator's DH frame assignment as SVG.

Writes docs/figures/dh_frames.svg and dh_frames_dark.svg, which README.md
embeds through a <picture> element so the figure follows the reader's theme.

Every dimension and every frame in the output is computed from ik_model's DH
table by calling ik_model.fk_all() - nothing here restates a link length or a
frame orientation as a literal.  Editing DH_A/DH_ALPHA/DH_D and re-running is
the whole update procedure, and a figure that disagrees with the kinematics is
therefore not a state this script can reach.

SVG rather than a raster format, and hand-emitted rather than matplotlib: the
output is resolution-independent, and the repository's only Python dependency
stays numpy.  The drawing primitives needed here are a projected line, an
arrowhead and a text anchor, which is not enough to justify a plotting stack.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ik_model as M  # noqa: E402

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                       "docs", "figures")

# Configuration drawn in panel A.  Chosen only so that no two frame origins
# project on top of each other at the viewing angle below; nothing about the
# figure's content depends on it.
Q_VIEW = np.array([0.55, -1.20, 0.55, 0.60, 0.70, 0.40])

# Panel B is the same arm with theta1 = 0 and the wrist joints zeroed.  On an
# arm with a lateral offset this does NOT put the chain in one plane: frames 0
# to 2 sit at y = 0 and frames 3 to 6 at y = -d3.  The two planes are parallel,
# though, so a projection onto x-z is still undistorted for every link - only
# the d3 step itself is edge-on, and it is annotated rather than drawn.
Q_PLANAR = np.array([0.0, -1.20, 0.55, 0.0, 0.0, 0.0])

# Camera azimuth is negative so that +x projects to the right of the frame: at
# a positive azimuth a forward reach runs down-left across the page and reads
# as an arm hanging rather than extending.
AZIM_DEG, ELEV_DEG = -58.0, 22.0
AXIS_LEN = 0.050          # metres, length of a drawn frame axis
GRID_HALF = 0.14          # metres, half-extent of the decorative floor grid
W, H = 1560, 820
PANEL_A = (0, 0, 900, H)
PANEL_B = (900, 0, 660, H)

LIGHT = dict(
    ink="#161b22", muted="#5b6673", faint="#c9d1d9", grid="#e6eaef",
    link="#3d4653", joint_fill="#ffffff", halo="#ffffff",
    x="#c92a3d", y="#1a8f3c", z="#1f5fd0",
    dim="#a5591b", l3="#6d3bd1", panel="#f6f8fa", rule="#d8dee4",
)
DARK = dict(
    ink="#e6edf3", muted="#9aa4b2", faint="#3a434e", grid="#22272e",
    link="#b3bdc9", joint_fill="#0d1117", halo="#0d1117",
    x="#ff8592", y="#57d364", z="#7aa7ff",
    dim="#e3a53f", l3="#bd9bff", panel="#161b22", rule="#30363d",
)

SANS = "system-ui,-apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif"
SERIF = "Georgia,'Times New Roman',serif"


# ---------------------------------------------------------------------------
# SVG primitives
# ---------------------------------------------------------------------------
def esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


class Canvas:
    def __init__(self):
        self.parts = []

    def add(self, s):
        self.parts.append(s)

    def line(self, p0, p1, stroke, width=2.0, dash=None, opacity=1.0,
             cap="round"):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(f'<line x1="{p0[0]:.2f}" y1="{p0[1]:.2f}" '
                 f'x2="{p1[0]:.2f}" y2="{p1[1]:.2f}" stroke="{stroke}" '
                 f'stroke-width="{width}" stroke-linecap="{cap}" '
                 f'opacity="{opacity}"{d}/>')

    def poly(self, pts, fill, opacity=1.0, stroke="none", width=1.0):
        s = " ".join(f"{x:.2f},{y:.2f}" for x, y in pts)
        self.add(f'<polygon points="{s}" fill="{fill}" stroke="{stroke}" '
                 f'stroke-width="{width}" opacity="{opacity}"/>')

    def path(self, d, stroke, width=2.0, fill="none", dash=None, opacity=1.0):
        da = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(f'<path d="{d}" fill="{fill}" stroke="{stroke}" '
                 f'stroke-width="{width}" stroke-linecap="round" '
                 f'stroke-linejoin="round" opacity="{opacity}"{da}/>')

    def circle(self, p, r, fill, stroke="none", width=1.5, opacity=1.0):
        self.add(f'<circle cx="{p[0]:.2f}" cy="{p[1]:.2f}" r="{r}" '
                 f'fill="{fill}" stroke="{stroke}" stroke-width="{width}" '
                 f'opacity="{opacity}"/>')

    def text(self, p, s, size=15, fill="#000", anchor="start", weight="400",
             family=SANS, style="normal", opacity=1.0, halo=None,
             baseline="middle"):
        common = (f'x="{p[0]:.2f}" y="{p[1]:.2f}" font-size="{size}" '
                  f'font-family="{family}" font-style="{style}" '
                  f'font-weight="{weight}" text-anchor="{anchor}" '
                  f'dominant-baseline="{baseline}"')
        if halo:
            # Painted twice: a fat stroke of the background colour under the
            # glyphs keeps a label legible where it crosses a link.
            self.add(f'<text {common} fill="none" stroke="{halo}" '
                     f'stroke-width="4.5" stroke-linejoin="round" '
                     f'opacity="{opacity}">{esc(s)}</text>')
        self.add(f'<text {common} fill="{fill}" opacity="{opacity}">'
                 f'{esc(s)}</text>')

    def arrow(self, p0, p1, colour, width=2.4, head=9.0, opacity=1.0):
        v = np.array(p1, float) - np.array(p0, float)
        n = float(np.hypot(*v))
        if n < 1e-9:
            return
        u = v / n
        base = np.array(p1, float) - u * head
        self.line(p0, base, colour, width, opacity=opacity, cap="butt")
        perp = np.array([-u[1], u[0]])
        self.poly([tuple(p1), tuple(base + perp * head * 0.38),
                   tuple(base - perp * head * 0.38)], colour, opacity=opacity)

    def svg(self, w, h):
        body = "\n".join(self.parts)
        return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" '
                f'height="{h}" viewBox="0 0 {w} {h}" '
                f'role="img" aria-label="Denavit-Hartenberg frame assignment '
                f'for the 6R manipulator">\n{body}\n</svg>\n')


# ---------------------------------------------------------------------------
# Projection
# ---------------------------------------------------------------------------
class View:
    """Orthographic projection onto a fitted panel rectangle."""

    def __init__(self, azim_deg, elev_deg):
        a, e = np.radians(azim_deg), np.radians(elev_deg)
        self.right = np.array([-np.sin(a), np.cos(a), 0.0])
        self.up = np.array([-np.cos(a) * np.sin(e), -np.sin(a) * np.sin(e),
                            np.cos(e)])
        self.scale = 1.0
        self.ox = self.oy = 0.0

    def flat(self, p):
        p = np.asarray(p, float)
        return np.array([p @ self.right, p @ self.up])

    def fit(self, pts, rect, margin):
        x0, y0, w, h = rect
        f = np.array([self.flat(p) for p in pts])
        lo, hi = f.min(axis=0), f.max(axis=0)
        span = np.maximum(hi - lo, 1e-6)
        self.scale = min((w - 2 * margin) / span[0], (h - 2 * margin) / span[1])
        mid = (lo + hi) / 2
        self.ox = x0 + w / 2 - mid[0] * self.scale
        self.oy = y0 + h / 2 + mid[1] * self.scale

    def __call__(self, p):
        f = self.flat(p)
        return (self.ox + f[0] * self.scale, self.oy - f[1] * self.scale)


class PlanarView:
    """x-z of the world onto the panel, for the dimensioned side elevation."""

    def __init__(self):
        self.scale = 1.0
        self.ox = self.oy = 0.0

    def flat(self, p):
        p = np.asarray(p, float)
        return np.array([p[0], p[2]])

    def fit(self, pts, rect, margin):
        x0, y0, w, h = rect
        f = np.array([self.flat(p) for p in pts])
        lo, hi = f.min(axis=0), f.max(axis=0)
        span = np.maximum(hi - lo, 1e-6)
        self.scale = min((w - 2 * margin) / span[0], (h - 2 * margin) / span[1])
        mid = (lo + hi) / 2
        self.ox = x0 + w / 2 - mid[0] * self.scale
        self.oy = y0 + h / 2 + mid[1] * self.scale

    def __call__(self, p):
        f = self.flat(p)
        return (self.ox + f[0] * self.scale, self.oy - f[1] * self.scale)


# ---------------------------------------------------------------------------
# Panel A: 3D frame assignment
# ---------------------------------------------------------------------------
# Where each frame's {i} label sits relative to its origin, in screen pixels.
# Tuned against the rendered output and connected by a leader line, because
# two pairs of frames are close enough that no automatic rule places them
# well: {2} and {3} are a3 = 2 cm apart, and {4} and {5} coincide exactly.
FRAME_LABEL_DXY = {
    0: (-30, 26), 1: (-42, 16), 2: (-46, -12), 3: (38, -14),
    4: (-56, -22), 5: (-58, 12), 6: (36, 6),
}

# Only the z axes are labelled: they are the joint axes, and labelling all 21
# arrows made the elbow and wrist unreadable.  x and y are identified by the
# colour key instead.
AXIS_LABEL_DXY = {0: (0, 0)}


def draw_panel_a(c, C):
    Ts = M.fk_all(Q_VIEW)
    origins = [T[:3, 3] for T in Ts]

    pts = []
    for T in Ts:
        o = T[:3, 3]
        pts.append(o)
        for k in range(3):
            pts.append(o + AXIS_LEN * T[:3, k])
    pts += [np.array([0, 0, 0.0]), np.array([GRID_HALF, GRID_HALF, 0.0]),
            np.array([-GRID_HALF, -GRID_HALF, 0.0])]

    v = View(AZIM_DEG, ELEV_DEG)
    v.fit(pts, (PANEL_A[0] + 26, PANEL_A[1] + 92, PANEL_A[2] - 52,
                PANEL_A[3] - 300), 34)

    # --- ground plane grid, z = 0 ---
    g = GRID_HALF
    step = 0.05
    t = -g
    while t <= g + 1e-9:
        c.line(v([t, -g, 0]), v([t, g, 0]), C["grid"], 1.0)
        c.line(v([-g, t, 0]), v([g, t, 0]), C["grid"], 1.0)
        t += step

    # --- base ---
    br = 0.055
    ring = [v([br * np.cos(a), br * np.sin(a), 0])
            for a in np.linspace(0, 2 * np.pi, 48)]
    c.poly(ring, C["faint"], opacity=0.30)
    c.path("M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in ring) + " Z",
           C["muted"], 1.4)

    # --- links ---
    seg = [v(o) for o in origins]
    for i in range(6):
        if np.linalg.norm(origins[i + 1] - origins[i]) < 1e-9:
            continue
        c.line(seg[i], seg[i + 1], C["link"], 7.0)

    # --- frames ---
    for i, T in enumerate(Ts):
        o = T[:3, 3]
        po = v(o)
        for k, key in enumerate(("x", "y", "z")):
            tip = v(o + AXIS_LEN * T[:3, k])
            c.arrow(po, tip, C[key], 2.4, 9.0)
        c.circle(po, 5.0, C["joint_fill"], C["ink"], 2.0)

    # Labels last, so that no arrow is drawn over one.
    for i, T in enumerate(Ts):
        po = np.array(v(T[:3, 3]))
        dx, dy = FRAME_LABEL_DXY[i]
        lp = po + np.array([dx, dy], float)
        n = np.hypot(dx, dy) or 1.0
        c.line(po + np.array([dx, dy]) / n * 9,
               lp - np.array([dx, dy]) / n * 13, C["muted"], 1.0, opacity=0.55)
        c.text(lp, "{%d}" % i, 17, C["ink"], "middle", "700", SANS,
               halo=C["halo"])

    # --- joint-variable callouts: joint i turns about z of frame i-1 ---
    for i in range(6):
        T = Ts[i]
        o = T[:3, 3]
        po = np.array(v(o))
        tip = np.array(v(o + AXIS_LEN * T[:3, 2]))
        d = tip - po
        n = np.hypot(*d) or 1.0
        perp = np.array([-d[1], d[0]]) / n
        # Put the label to whichever side of the joint axis points away from
        # the arm's centroid, so it never lands on top of a link.
        cen = np.mean([v(p) for p in origins], axis=0)
        if np.dot(perp, po - cen) < 0:
            perp = -perp
        lp = po + d * 0.62 + perp * 26
        c.text(lp, f"θ{i + 1}", 17, C["z"], "middle", "700", SERIF, "italic",
               halo=C["halo"])

    # --- tool approach axis ---
    o6 = Ts[6][:3, 3]
    tip = v(o6 + AXIS_LEN * 1.7 * Ts[6][:3, 2])
    c.line(v(o6), tip, C["z"], 1.6, dash="5 4", opacity=0.75)
    c.text((tip[0] + 26, tip[1] + 34), "approach axis", 13, C["muted"],
           "middle", "500", SANS, halo=C["halo"])

    # --- titles ---
    c.text((PANEL_A[0] + 34, 44), "A.  Frame assignment", 22, C["ink"],
           "start", "700")
    c.text((PANEL_A[0] + 34, 72),
           "standard (distal) DH; frame {i} attached to link i, "
           "joint i turns about z of frame {i−1}",
           14, C["muted"], "start", "400")

    # --- legend ---
    lx, ly = PANEL_A[0] + 34, H - 176
    c.text((lx, ly - 26), "axis colour key", 13, C["ink"], "start", "700")
    for k, (key, lab) in enumerate((("x", "x — common normal"),
                                    ("y", "y — completes the triad"),
                                    ("z", "z — joint axis"))):
        y = ly + k * 23
        c.arrow((lx, y), (lx + 28, y), C[key], 2.6, 8.5)
        c.text((lx + 38, y), lab, 13, C["muted"], "start", "500")

    nx = lx + 250
    c.text((nx, ly - 26), "the two structural features that matter", 13,
           C["ink"], "start", "700")
    for k, s in enumerate((
            "{4} and {5} share an origin exactly, since a₅ = d₅ = 0.  That is",
            "what makes the wrist spherical, and it is why the analytic solver",
            "can decouple position from orientation.  {2} and {3} are set apart",
            "by the lateral offset d₃, which keeps the wrist centre off the",
            "joint-1 axis where θ₁ would be undefined.")):
        c.text((nx, ly + k * 20), s, 13, C["muted"], "start", "400")

    c.text((lx, H - 66),
           f"drawn at q = ({', '.join('%.2f' % x for x in Q_VIEW)}) rad; "
           "the configuration is arbitrary and affects nothing but legibility",
           12, C["faint"], "start", "400", SERIF, "italic")


# ---------------------------------------------------------------------------
# Panel B: dimensioned side elevation
# ---------------------------------------------------------------------------
def dim_line(c, C, p0, p1, label, off=(0, 0), tick=6, size=14, colour=None):
    """A dimension line with end ticks and a label at its midpoint."""
    colour = colour or C["dim"]
    p0 = np.array(p0, float)
    p1 = np.array(p1, float)
    d = p1 - p0
    n = np.hypot(*d) or 1.0
    u = d / n
    perp = np.array([-u[1], u[0]])
    c.line(p0, p1, colour, 1.5)
    for p in (p0, p1):
        c.line(p - perp * tick, p + perp * tick, colour, 1.5)
    mid = (p0 + p1) / 2 + np.array(off, float)
    c.text(mid, label, size, colour, "middle", "600", SERIF, "italic",
           halo=C["halo"])


def draw_panel_b(c, C):
    Ts = M.fk_all(Q_PLANAR)
    o = [T[:3, 3] for T in Ts]
    # Every frame must sit at y = 0 or y = -d3, i.e. in one of two planes
    # parallel to x-z; anything else would make the projection below lie.
    for i, p in enumerate(o):
        assert min(abs(p[1]), abs(abs(p[1]) - M.D3)) < 1e-9, \
            f"frame {i} is at y = {p[1]:.4f}, neither 0 nor -d3; the side " \
            "elevation would be foreshortened"

    # The drawing keeps to the upper two thirds; the lower third is the note
    # block, and the left inset is where the d1 dimension line lives.
    v = PlanarView()
    v.fit(o, (PANEL_B[0] + 112, PANEL_B[1] + 112, PANEL_B[2] - 236, 452), 24)

    P = [np.array(v(p)) for p in o]

    # Elbow decomposition: o2 -> o3 is a3, o3 -> o4 is d4, and the closing
    # side is the effective forearm L3 = hypot(a3, d4) at PHI from a3.  Both
    # are read from ik_model, not recomputed here.
    c.line(P[0], P[1], C["link"], 6.0)
    c.line(P[1], P[2], C["link"], 6.0)
    c.line(P[2], P[3], C["link"], 6.0)
    c.line(P[3], P[4], C["link"], 6.0)
    c.line(P[4], P[6], C["link"], 6.0)

    # L3 closing side of the forearm triangle
    c.line(P[2], P[4], C["l3"], 2.2, dash="7 5")

    # base ground hatch
    gy = P[0][1]
    c.line((P[0][0] - 70, gy), (P[0][0] + 70, gy), C["muted"], 2.0)
    for k in range(-6, 7):
        x = P[0][0] + k * 11
        c.line((x, gy), (x - 9, gy + 9), C["faint"], 1.4)

    for i, p in enumerate(P):
        if i == 5:
            continue
        c.circle(p, 5.0, C["joint_fill"], C["ink"], 2.0)

    # --- dimensions, all values from ik_model ---
    o0, o1, o2, o3, o4, o6 = o[0], o[1], o[2], o[3], o[4], o[6]

    # d1: vertical rise of frame 1 above the base plane
    xoff = -62
    a = np.array(v([o0[0], 0, o0[2]])) + np.array([xoff, 0])
    b = np.array(v([o0[0], 0, o1[2]])) + np.array([xoff, 0])
    dim_line(c, C, a, b, f"d₁ = {M.D1:g}", off=(-26, 0))
    c.line(np.array(v(o1)), b, C["dim"], 1.0, dash="4 4", opacity=0.7)

    # a1: horizontal offset of frame 1.  Zero on this arm, so there is no
    # dimension to draw - the shoulder axis meets the base axis - and a
    # zero-length dimension line would read as a drawing error.
    if abs(M.A1) > 1e-9:
        yoff = -34
        a = np.array(v([o0[0], 0, o1[2]])) + np.array([0, yoff])
        b = np.array(v([o1[0], 0, o1[2]])) + np.array([0, yoff])
        dim_line(c, C, a, b, f"a₁ = {M.A1:g}", off=(0, -16), tick=5)
        c.line(np.array(v(o1)), b, C["dim"], 1.0, dash="4 4", opacity=0.7)
    else:
        c.text(np.array(v(o1)) + np.array([16, -18]),
               f"a₁ = {M.A1:g}: joint axes 1 and 2 intersect", 13, C["dim"],
               "start", "600", SERIF, "italic", halo=C["halo"])

    # a2: the upper arm itself
    mid = (P[1] + P[2]) / 2
    d = P[2] - P[1]
    n = np.hypot(*d) or 1.0
    perp = np.array([-d[1], d[0]]) / n
    c.text(mid + perp * 22, f"a₂ = {M.A2:g}", 15, C["dim"], "middle",
           "600", SERIF, "italic", halo=C["halo"])

    # a3 is 2 cm and its own segment is barely longer than the joint dot, so
    # it gets a leader line out to clear space rather than an inline label.
    mid = (P[2] + P[3]) / 2
    lab = mid + np.array([26, -46])
    c.line(mid, lab + np.array([-6, 8]), C["dim"], 1.0, opacity=0.75)
    c.text(lab, f"a₃ = {M.A3:g}", 14, C["dim"], "start", "600", SERIF,
           "italic", halo=C["halo"])

    mid = (P[3] + P[4]) / 2
    d = P[4] - P[3]
    n = np.hypot(*d) or 1.0
    perp = np.array([-d[1], d[0]]) / n
    c.text(mid - perp * 26, f"d₄ = {M.D4:g}", 15, C["dim"], "middle",
           "600", SERIF, "italic", halo=C["halo"])
    # L3's label rides at a third of the way along rather than the midpoint,
    # which keeps it clear of the labels crowding the wrist centre.
    d = P[4] - P[2]
    n = np.hypot(*d) or 1.0
    perp = np.array([-d[1], d[0]]) / n
    c.text(P[2] + d * 0.55 + perp * 34, f"L₃ = {M.L3:.4f}", 15, C["l3"],
           "middle", "700", SERIF, "italic", halo=C["halo"])

    # d6: the tool segment is 8 cm and its midpoint sits inside the wrist
    # marker, so this label is led out to the right as well.
    mid = (P[4] + P[6]) / 2
    lab = mid + np.array([44, -34])
    c.line(mid, lab + np.array([-6, 8]), C["dim"], 1.0, opacity=0.75)
    c.text(lab, f"d₆ = {M.D6:g}", 14, C["dim"], "start", "600", SERIF,
           "italic", halo=C["halo"])

    # wrist centre callout: the point the analytic solver locates first
    c.circle(P[4], 10.0, "none", C["l3"], 1.6, opacity=0.9)
    c.text(P[4] + np.array([-16, 28]), "wrist centre", 12, C["muted"], "end",
           "500", SANS, halo=C["halo"])
    c.text(P[6] + np.array([12, 14]), "tool", 12, C["muted"], "start",
           "500", SANS, halo=C["halo"])

    # --- titles and notes ---
    c.text((PANEL_B[0] + 30, 44), "B.  Link dimensions", 22, C["ink"],
           "start", "700")
    c.text((PANEL_B[0] + 30, 72),
           "side elevation, θ₁ = 0, wrist zeroed; metres",
           14, C["muted"], "start", "400")

    # d3 is perpendicular to the page in this view, so it gets a marker on the
    # frame-3 origin and a note rather than a dimension line.
    c.circle(P[3], 11.0, "none", C["dim"], 1.6, opacity=0.9)
    c.line(P[3], P[3] + np.array([34, 40]), C["dim"], 1.0, opacity=0.75)
    c.text(P[3] + np.array([38, 44]),
           f"d₃ = {M.D3:g} into the page", 13, C["dim"], "start", "600",
           SERIF, "italic", halo=C["halo"])

    ny = H - 232
    c.text((PANEL_B[0] + 30, ny), "the effective forearm", 13, C["ink"],
           "start", "700")
    c.text((PANEL_B[0] + 30, ny + 28),
           f"L₃ = hypot(a₃, d₄) = {M.L3:.4f}", 15, C["l3"], "start", "700",
           SERIF, "italic")
    c.text((PANEL_B[0] + 30, ny + 52),
           f"φ = atan2(d₄, a₃) = {np.degrees(M.PHI):.2f}°", 15, C["l3"],
           "start", "700", SERIF, "italic")
    sign = "-" if M.PHI > 0 else "+"
    for k, s in enumerate((
            f"a₃ is {M.A3 * 100:.2f} cm against d₄'s {M.D4 * 100:.1f} cm, so the forearm is",
            f"nearly straight and φ sits {90 - np.degrees(M.PHI):.1f}° off a right angle.",
            "ik_analytic reduces the position problem to a planar two-link",
            f"chain in a₂ and L₃, solves for γ, and recovers θ₃ = γ {sign} φ.",
            "Its sensitivity goes as 1/|sin γ| - the elbow conditioning",
            "number - and, because d₃ is non-zero, θ₁ has a second one of",
            "its own: see shoulder_conditioning() in model/ik_model.py.")):
        c.text((PANEL_B[0] + 30, ny + 82 + k * 19), s, 13, C["muted"],
               "start", "400")


# ---------------------------------------------------------------------------
def build(C):
    c = Canvas()
    c.add(f'<rect width="{W}" height="{H}" fill="none"/>')
    c.add(f'<rect x="{PANEL_B[0]}" y="24" width="{PANEL_B[2] - 24}" '
          f'height="{H - 48}" rx="10" fill="{C["panel"]}" opacity="0.6"/>')
    c.line((PANEL_B[0] - 8, 40), (PANEL_B[0] - 8, H - 40), C["rule"], 1.2)
    draw_panel_a(c, C)
    draw_panel_b(c, C)

    c.text((34, H - 26),
           "Generated by model/plot_dh.py from the DH table in "
           "model/ik_model.py.  All lengths in metres.",
           12, C["faint"], "start", "400")
    return c.svg(W, H)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for name, palette in (("dh_frames.svg", LIGHT),
                          ("dh_frames_dark.svg", DARK)):
        p = os.path.join(OUT_DIR, name)
        with open(p, "w") as f:
            f.write(build(palette))
        print(f"  {name:22s} -> {os.path.relpath(p)}")


if __name__ == "__main__":
    print("Rendering DH frame figure:")
    main()
    print("done.")
