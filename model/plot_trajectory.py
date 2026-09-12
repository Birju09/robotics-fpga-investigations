#!/usr/bin/env python3
"""
Render the pentagon trajectory workload (Workload B) as SVG.

Writes docs/figures/trajectory.svg and trajectory_dark.svg, which README.md
embeds through a <picture> element alongside the DH frame figure.

Every coordinate is recomputed from model/gen_vectors.py's TRAJ_* constants
and model/robot.py's ROBOT geometry - nothing here restates the centre,
radius, height or reach as a literal, so the figure cannot silently disagree
with the workload it is describing.

Panel A is a top-down (x-y) view of the base plane: the pentagon path, its
vertices, the sample density along one edge, and the wrist-centre reachable
annulus for scale. Panel B is a schematic of one edge showing the linear yaw
interpolation between vertices, which is the change that removed the
step-at-vertex orientation discontinuity from an earlier version of this
workload.

Drawing primitives (Canvas, esc) are imported from plot_dh rather than
reimplemented.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_vectors as G  # noqa: E402
import robot as R         # noqa: E402
from plot_dh import Canvas, LIGHT, DARK, SANS  # noqa: E402

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                       "docs", "figures")

W, H = 1560, 760
PANEL_A = (0, 0, 900, H)
PANEL_B = (900, 0, 660, H)


class View2D:
    """Orthographic x-y projection fitted to a panel rectangle."""

    def __init__(self):
        self.scale = 1.0
        self.ox = self.oy = 0.0

    def fit(self, pts, rect, margin):
        x0, y0, w, h = rect
        f = np.asarray(pts, float)
        lo, hi = f.min(axis=0), f.max(axis=0)
        span = np.maximum(hi - lo, 1e-6)
        self.scale = min((w - 2 * margin) / span[0], (h - 2 * margin) / span[1])
        mid = (lo + hi) / 2
        self.ox = x0 + w / 2 - mid[0] * self.scale
        self.oy = y0 + h / 2 + mid[1] * self.scale

    def __call__(self, p):
        x, y = p[0], p[1]
        return (self.ox + x * self.scale, self.oy - y * self.scale)


def _vertices():
    cx, cy = G.TRAJ_CENTRE
    return [(cx + G.TRAJ_RADIUS * np.cos(np.pi / 2 + 2 * np.pi * k / 5),
             cy + G.TRAJ_RADIUS * np.sin(np.pi / 2 + 2 * np.pi * k / 5))
            for k in range(5)]


def _path_samples():
    """One lap's worth of (x, y) points at TRAJ_STEPS per edge."""
    vtx = _vertices()
    pts = []
    for k in range(5):
        x0, y0 = vtx[k]
        x1, y1 = vtx[(k + 1) % 5]
        for j in range(G.TRAJ_STEPS):
            t = j / G.TRAJ_STEPS
            pts.append((x0 + (x1 - x0) * t, y0 + (y1 - y0) * t))
    return pts


def draw_panel_a(c, C):
    vtx = _vertices()
    samples = _path_samples()
    cx, cy = G.TRAJ_CENTRE

    # Wrist-centre reachable annulus, for scale: the pentagon sits well inside
    # it, which is the point of the search that sited TRAJ_CENTRE/TRAJ_RADIUS
    # in the first place (see the comment above those constants).
    ro, ri = R.ROBOT.reach_outer, R.ROBOT.reach_inner
    ring_pts = [(ro * np.cos(a), ro * np.sin(a))
                for a in np.linspace(0, 2 * np.pi, 96)]
    ring_pts += [(0, 0), (ri, 0)]

    v = View2D()
    fit_pts = ring_pts + vtx + [(cx, cy)]
    v.fit(fit_pts, (PANEL_A[0] + 30, PANEL_A[1] + 90, PANEL_A[2] - 60,
                    PANEL_A[3] - 260), 20)

    c.text((PANEL_A[0] + 30, 42), "A", 20, C["ink"], "start", "700", SANS)
    c.text((PANEL_A[0] + 58, 42),
           "Workload B: pentagon trajectory, base-plane view", 16, C["muted"],
           "start", "500")

    # --- reachable envelope. reach_inner is 4.8e-4 m on this arm - the whole
    # point of the lateral shoulder offset (README §4.5, §2.1) is that it
    # shrinks the unreachable singular region to a pinhole rather than a
    # cylinder - so it is not drawn: at this scale it would not be visible
    # and a dashed circle at that radius would only misrepresent it as larger
    # than it is.
    outer = [v(p) for p in ((ro * np.cos(a), ro * np.sin(a))
                            for a in np.linspace(0, 2 * np.pi, 96))]
    c.path("M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in outer) + " Z",
           C["faint"], 1.2, dash="2,4")
    c.text(v((ro * 0.55, ro * 0.72)), "wrist-centre reach", 12, C["muted"],
           "middle", "400")
    c.text(v((ro * 0.55, ro * 0.72 - 0.045)), f"(outer radius {ro:.3f} m; "
           f"inner exclusion {ri * 1000:.2f} mm, not to scale)", 11,
           C["faint"], "middle", "400")

    # --- base origin ---
    o = v((0, 0))
    c.circle(o, 4.5, C["ink"])
    c.text((o[0], o[1] + 20), "base origin", 12, C["muted"], "middle")

    # --- the path itself ---
    pverts = [v(p) for p in vtx]
    closed = pverts + [pverts[0]]
    c.path("M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in closed),
           C["z"], 2.6)

    # --- samples along the path, so the figure shows the actual sampling
    # density rather than just the polygon it approximates ---
    for p in samples:
        sp = v(p)
        c.circle(sp, 2.0, C["z"], opacity=0.55)

    # --- vertices, labelled with their commanded yaw ---
    cpt = v((cx, cy))
    for k, (p, sp) in enumerate(zip(vtx, pverts)):
        c.circle(sp, 6.5, C["panel"], stroke=C["dim"], width=2.2)
        c.circle(sp, 2.4, C["dim"])
        dx, dy = (sp[0] - cpt[0]), (sp[1] - cpt[1])
        n = np.hypot(dx, dy) or 1.0
        lx, ly = sp[0] + dx / n * 34, sp[1] + dy / n * 34
        c.text((lx, ly - 9), f"v{k}", 13, C["ink"], "middle", "600")
        c.text((lx, ly + 10), f"ψ={G.TRAJ_YAW_DEG[k]:g}°", 12,
               C["muted"], "middle", "400")

    # --- centre marker. The radius value is in the caption block below
    # rather than dimensioned here - at this panel's scale (the pentagon
    # against the full reach envelope) any leader line to a radius label
    # runs into a vertex label on one side of the pentagon or the other. ---
    c.circle(cpt, 2.5, C["muted"])

    ny = PANEL_A[3] - 150
    perim = 2 * G.TRAJ_RADIUS * np.sin(np.pi / 5) * 5
    step_mm = 1000.0 * perim / (5 * G.TRAJ_STEPS)
    for k, s in enumerate((
            f"centre ({cx:.2f}, {cy:.2f}) m, z = {G.TRAJ_Z:.2f} m, "
            f"tool axis vertical",
            f"circumradius {G.TRAJ_RADIUS:.2f} m, perimeter {perim:.3f} m",
            f"{G.TRAJ_STEPS} samples/edge -> {5 * G.TRAJ_STEPS} samples/lap, "
            f"{step_mm:.2f} mm/sample",
            f"yaw {' - '.join('%g°' % y for y in G.TRAJ_YAW_DEG)}, "
            f"interpolated linearly per edge",
    )):
        c.text((PANEL_A[0] + 30, ny + k * 20), s, 13, C["muted"], "start")


def draw_panel_b(c, C):
    c.add(f'<rect x="{PANEL_B[0]}" y="24" width="{PANEL_B[2] - 24}" '
          f'height="{H - 48}" rx="10" fill="{C["panel"]}" opacity="0.6"/>')

    x0, y0 = PANEL_B[0] + 50, 100
    x1 = PANEL_B[0] + PANEL_B[2] - 90
    y1 = y0

    c.text((PANEL_B[0] + 30, 42), "B", 20, C["ink"], "start", "700")
    c.text((PANEL_B[0] + 58, 42), "Orientation along one edge", 16,
           C["muted"], "start", "500")

    # One edge, drawn straight and to a distorted scale purely as a schematic
    # for the yaw ramp - it is not another projection of panel A.
    c.line((x0, y0), (x1, y1), C["link"], 4.0)
    c.circle((x0, y0), 7, C["panel"], stroke=C["dim"], width=2.4)
    c.circle((x1, y1), 7, C["panel"], stroke=C["dim"], width=2.4)
    c.text((x0, y0 - 22), "vertex k", 13, C["ink"], "middle", "600")
    c.text((x1, y1 - 22), "vertex k+1", 13, C["ink"], "middle", "600")
    c.text((x0, y0 + 24), "ψ = 0°", 12, C["muted"], "middle")
    c.text((x1, y1 + 24), "ψ = 90°", 12, C["muted"], "middle")

    n = G.TRAJ_STEPS if G.TRAJ_STEPS <= 12 else 10
    for j in range(1, n):
        t = j / n
        px = x0 + (x1 - x0) * t
        psi = t * 90.0
        r = 16
        ay0 = y0 + 56
        c.circle((px, ay0), 3.0, C["z"])
        ang = np.radians(90 - psi)
        tip = (px + r * np.cos(ang), ay0 - r * np.sin(ang))
        c.arrow((px, ay0), tip, C["z"], width=1.6, head=5.0)

    c.text((x0, y0 + 96),
           "tool-frame yaw arrow, sampled every 1/%d of the edge" %
           (G.TRAJ_STEPS if G.TRAJ_STEPS <= 12 else 10), 12, C["muted"],
           "start")

    ny = y0 + 150
    for k, s in enumerate((
            "Yaw is interpolated LINEARLY in psi along each edge, not held",
            "constant and stepped at the vertex.  An earlier version of this",
            "workload switched orientation discontinuously at each vertex;",
            "that produced a step input to the DLS solver at every vertex",
            "and was corrected to the ramp shown here, which is what a",
            "tool-frame trajectory generator would actually emit.",
            "",
            "The ramp is why the tool never turns faster than "
            f"{90.0 / G.TRAJ_STEPS:.2f}°/sample -",
            "a servo-scale increment rather than a discontinuity - and it is",
            "part of why the tracking workload (§5.3) converges in so",
            "few iterations: nothing about the commanded pose jumps.",
    )):
        c.text((PANEL_B[0] + 30, ny + k * 19), s, 13, C["muted"], "start")


def build(C):
    c = Canvas()
    c.add(f'<rect width="{W}" height="{H}" fill="none"/>')
    c.line((PANEL_B[0] - 8, 40), (PANEL_B[0] - 8, H - 40), C["rule"], 1.2)
    draw_panel_a(c, C)
    draw_panel_b(c, C)
    c.text((34, H - 20),
           "Generated by model/plot_trajectory.py from the TRAJ_* constants "
           "in model/gen_vectors.py.  All lengths in metres.",
           12, C["faint"], "start", "400")
    return c.svg(W, H, label="Pentagon trajectory workload: base-plane path "
                             "with sample density, and the linear yaw ramp "
                             "between vertices")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for name, palette in (("trajectory.svg", LIGHT),
                          ("trajectory_dark.svg", DARK)):
        p = os.path.join(OUT_DIR, name)
        with open(p, "w") as f:
            f.write(build(palette))
        print(f"  {name:22s} -> {os.path.relpath(p)}")


if __name__ == "__main__":
    print("Rendering trajectory figure:")
    main()
    print("done.")
