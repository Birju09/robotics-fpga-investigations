#!/usr/bin/env python3
"""
Vitis 2025.2 - build the bare-metal timing application.

Run with the Vitis Python interpreter, not a system one:

    vitis -s scripts/build_vitis.py
    vitis -s scripts/build_vitis.py --xsa build/vivado/ik_platform.xsa --clean

Produces build/vitis/ik_app/build/ik_app.elf.

Sources are staged with shutil, not the API's import helper, since the
helper's name/signature has moved between releases. Both components are
deleted and recreated every run: the platform caches the hardware handoff
from its XSA, so reusing one after a Vivado rebuild would silently compile
against the previous bitstream's xparameters.h.
"""

import argparse
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PLATFORM = "ik_platform"
APP = "ik_app"
CPU = "ps7_cortexa9_0"        # Zynq-7000 has Cortex-A9s, not A53s
DOMAIN = "standalone_a9"

# The bare-metal app reuses the HLS kernel sources in their float
# configuration, so the PS reference runs the same algorithm as the PL.
SW_SOURCES = [
    ("sw/src", "main.c"),
    ("sw/src", "ik_driver.c"),
    ("sw/src", "ik_sw_ref.cpp"),
    ("sw/src", "ik_driver.h"),
    ("sw/src", "ik_regmap.h"),
    ("sw/src", "ik_vectors.h"),
    ("hls/src", "kinematics.cpp"),
    ("hls/src", "matmul.cpp"),
    ("hls/src", "matinv.cpp"),
    ("hls/src", "spd.cpp"),
    ("hls/src", "ik_analytic.cpp"),
    ("hls/src", "ik_dls.cpp"),
]

HEADER_DIRS = ["hls/include"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--xsa", default=os.path.join(ROOT, "build", "vivado",
                                                 "ik_platform.xsa"))
    p.add_argument("--workspace", default=os.path.join(ROOT, "build", "vitis"))
    p.add_argument("--clean", action="store_true",
                   help="delete the whole workspace first, BSP included; the "
                        "platform and app components are recreated on every "
                        "run regardless")
    # `vitis -s` passes its own arguments through; ignore anything unknown.
    args, unknown = p.parse_known_args()
    if unknown:
        print("note: ignoring extra arguments %s" % unknown)
    return args


def stage_sources(dest):
    """Copy the application sources and the shared headers into the component."""
    os.makedirs(dest, exist_ok=True)
    n = 0
    for sub, name in SW_SOURCES:
        src = os.path.join(ROOT, sub, name)
        if not os.path.exists(src):
            sys.exit("missing source: %s\n"
                     "  (ik_vectors.h comes from model/gen_vectors.py,\n"
                     "   ik_regmap.h from scripts/gen_regmap.py)" % src)
        shutil.copy2(src, os.path.join(dest, name))
        n += 1

    inc = os.path.join(dest, "include")
    os.makedirs(inc, exist_ok=True)
    for d in HEADER_DIRS:
        for f in os.listdir(os.path.join(ROOT, d)):
            if f.endswith((".hpp", ".h")):
                shutil.copy2(os.path.join(ROOT, d, f), os.path.join(inc, f))
                n += 1
    print("  staged %d files into %s" % (n, dest))


def report_kernel_xpar(workspace):
    """Print the XPAR_* base-address macros the platform generated.

    main.c's #if cascade over known xparameters.h spellings can only say
    what it didn't find if none match; this prints what's actually there.
    """
    hits = []
    for dirpath, _dirs, files in os.walk(workspace):
        if "xparameters.h" in files:
            hits.append(os.path.join(dirpath, "xparameters.h"))
    if not hits:
        print("  note: no xparameters.h found under %s" % workspace)
        return
    names = ("KERNEL", "MAT_MUL", "MAT_INV", "IK_ANALYTIC", "IK_DLS")
    for h in sorted(hits, key=len):
        lines = []
        with open(h, "r", errors="replace") as fh:
            for line in fh:
                if line.startswith("#define XPAR_") and any(n in line for n in names):
                    lines.append(line.rstrip())
        if lines:
            print("  kernel base addresses in %s:" % h)
            for line in lines:
                print("      %s" % line)
            return
    print("  note: no kernel XPAR_* macros in any xparameters.h under %s"
          % workspace)
    print("        the bitstream may not contain the kernels you expect;")
    print("        check 'INFO: building kernels' in the Vivado log.")


def drop_component(client, workspace, name):
    """Remove a component left from an earlier run (else create_*_component
    dies with ALREADY_EXISTS). Rebuilding avoids silently compiling against a
    stale xparameters.h after a Vivado rebuild.
    """
    path = os.path.join(workspace, name)
    if not os.path.isdir(path):
        return
    print("  removing existing component '%s'" % name)
    try:
        client.delete_component(name=name)
        return
    except Exception as e:                             # noqa: BLE001
        # API name has moved between releases; fall back to removing the dir.
        print("    (delete_component failed: %s - removing the directory)" % e)
    shutil.rmtree(path, ignore_errors=True)


def main():
    args = parse_args()

    if not os.path.exists(args.xsa):
        sys.exit("XSA not found: %s\n"
                 "  Build it first:  vivado -mode batch -source scripts/build_vivado.tcl"
                 % args.xsa)

    try:
        import vitis
    except ImportError:
        sys.exit("This script must run under the Vitis interpreter:\n"
                 "    vitis -s scripts/build_vitis.py")

    if args.clean and os.path.isdir(args.workspace):
        print("removing %s" % args.workspace)
        shutil.rmtree(args.workspace)
    os.makedirs(args.workspace, exist_ok=True)

    client = vitis.create_client()
    client.set_workspace(args.workspace)

    # ---------------- platform ----------------
    print("creating platform component '%s' from %s" % (PLATFORM, args.xsa))
    drop_component(client, args.workspace, PLATFORM)
    plat = client.create_platform_component(
        name=PLATFORM,
        hw_design=args.xsa,
        os="standalone",
        cpu=CPU,
        domain_name=DOMAIN,
    )
    plat.build()

    xpfm = os.path.join(args.workspace, PLATFORM, "export", PLATFORM,
                        PLATFORM + ".xpfm")
    if not os.path.exists(xpfm):
        # Layout has moved between releases; find it rather than guess.
        hits = []
        for dirpath, _dirs, files in os.walk(os.path.join(args.workspace, PLATFORM)):
            for f in files:
                if f.endswith(".xpfm"):
                    hits.append(os.path.join(dirpath, f))
        if not hits:
            sys.exit("platform built but no .xpfm found under %s" % args.workspace)
        xpfm = sorted(hits, key=len)[0]
    print("  platform: %s" % xpfm)
    report_kernel_xpar(os.path.join(args.workspace, PLATFORM))

    # ---------------- application ----------------
    print("creating app component '%s'" % APP)
    drop_component(client, args.workspace, APP)
    app = client.create_app_component(
        name=APP,
        platform=xpfm,
        domain=DOMAIN,
        template="empty_application",
    )

    app_src = os.path.join(args.workspace, APP, "src")
    stage_sources(app_src)

    # Absolute include path, not '${_ide_ws}/...': the 2025.2 CMake build
    # doesn't substitute that variable, so it expanded to nothing and every
    # kernel .cpp failed to find its own header (error pointed at the wrong
    # place). The -Wno- flags suppress HLS pragma/label warnings that mean
    # nothing to a host compiler and previously buried a real diagnostic
    # under 500+ lines of noise.
    flags = ('-I%s -DIK_USE_FLOAT -O2 -Wno-unused-label -Wno-unknown-pragmas'
             % os.path.join(app_src, "include"))
    for key, val in (("USER_COMPILE_OTHER_FLAGS", flags),
                     ("USER_COMPILE_DEBUG_LEVEL", "-g")):
        try:
            app.set_app_config(key=key, values=val)
        except Exception as e:                     # noqa: BLE001
            print("  warning: set_app_config(%s) failed: %s" % (key, e))
            print("           set it by hand in the Vitis IDE if the build fails.")

    print("building...")
    app = client.get_component(name=APP)
    app.build()

    elf = None
    for dirpath, _dirs, files in os.walk(os.path.join(args.workspace, APP)):
        for f in files:
            if f.endswith(".elf"):
                elf = os.path.join(dirpath, f)
    if elf:
        print("\nELF: %s" % elf)
        print("\nRun it with:")
        print("    vitis -s -c 'import xsdb; ...'   or use the Vitis IDE debugger")
        print("    (connect, program %s, then load and run the ELF)" % args.xsa)
    else:
        print("\nbuild finished but no .elf was found; check the Vitis log.")

    vitis.dispose()


if __name__ == "__main__":
    main()
