#!/usr/bin/env python3
"""
Generate sw/src/ik_regmap.h from the register maps Vitis HLS emits.

Vitis HLS assigns AXI4-Lite offsets by argument order and width, so any change
to a kernel signature silently moves them.  Rather than transcribe the offsets
by hand and let them rot, this pulls them out of the generated driver headers:

    hls/build/<kernel>/**/impl/ip/drivers/*/src/x<kernel>_hw.h

Run after `make -C hls ip`:

    python3 scripts/gen_regmap.py

The header it writes is consumed by sw/src/ik_driver.c.  main.c additionally
verifies the map at runtime, so a stale header is caught before it can turn
into plausible-looking wrong answers.
"""

import os
import re
import sys
import glob

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HLS_BUILD = os.path.join(ROOT, "hls", "build")
OUT = os.path.join(ROOT, "sw", "src", "ik_regmap.h")

KERNELS = [
    "mat_mul_kernel",
    "mat_inv_kernel",
    "ik_analytic_kernel",
    "ik_dls_kernel",
]

DEFINE_RE = re.compile(r"^\s*#define\s+(\S+)\s+(0x[0-9a-fA-F]+|\d+)\s*$")


def find_hw_header(kernel):
    pats = [
        os.path.join(HLS_BUILD, "**", "x%s_hw.h" % kernel),
        os.path.join(HLS_BUILD, "**", "X%s_hw.h" % kernel),
    ]
    hits = []
    for p in pats:
        hits.extend(glob.glob(p, recursive=True))
    if not hits:
        return None
    # Prefer the shallowest path: deeper copies tend to be stale build debris.
    hits.sort(key=lambda p: (len(p.split(os.sep)), p))
    return hits[0]


def parse(path):
    out = {}
    with open(path, "r", errors="replace") as f:
        for line in f:
            m = DEFINE_RE.match(line)
            if m:
                out[m.group(1)] = int(m.group(2), 0)
    return out


def main():
    if not os.path.isdir(HLS_BUILD):
        sys.exit("No hls/build directory. Run 'make -C hls ip' first.")

    blocks = []
    missing = []

    for k in KERNELS:
        hdr = find_hw_header(k)
        if hdr is None:
            missing.append(k)
            continue
        defs = parse(hdr)
        prefix = "X%s_CTRL_" % k.upper()
        entries = {n: v for n, v in defs.items() if n.startswith(prefix)}
        if not entries:
            # Some releases use the bus name rather than CTRL.
            entries = {n: v for n, v in defs.items()
                       if n.startswith("X%s_" % k.upper()) and "ADDR" in n}
        if not entries:
            missing.append(k)
            continue

        rel = os.path.relpath(hdr, ROOT)
        lines = ["/* ---- %s ---- */" % k,
                 "/*   from %s */" % rel]
        for n in sorted(entries, key=lambda x: (entries[x], x)):
            lines.append("#define %-56s 0x%02x" % (n, entries[n]))
        blocks.append("\n".join(lines))
        print("  %-22s %-3d offsets  <- %s" % (k, len(entries), rel))

    if missing:
        print("\nWARNING: no generated header found for: %s" % ", ".join(missing))
        print("         Run 'make -C hls ip' for those kernels and re-run this.")
        if len(missing) == len(KERNELS):
            sys.exit(1)

    body = "\n\n".join(blocks)
    with open(OUT, "w") as f:
        f.write(
            "/*\n"
            " * GENERATED FILE - do not edit.\n"
            " *\n"
            " * Produced by scripts/gen_regmap.py from the register maps Vitis HLS\n"
            " * emits alongside each packaged IP.  Re-run that script after any\n"
            " * change to a kernel's argument list; HLS assigns these offsets by\n"
            " * argument order and width, so adding one scalar shifts everything\n"
            " * after it.\n"
            " */\n"
            "#ifndef IK_REGMAP_H\n"
            "#define IK_REGMAP_H\n\n"
            "#define IK_REGMAP_GENERATED 1\n\n"
            + body +
            "\n\n#endif /* IK_REGMAP_H */\n")
    print("\nwrote %s" % os.path.relpath(OUT, ROOT))


if __name__ == "__main__":
    print("Extracting AXI4-Lite register maps:")
    main()
