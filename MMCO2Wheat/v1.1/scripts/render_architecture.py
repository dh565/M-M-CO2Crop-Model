#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial
# research use with citation; redistribution of modified versions and commercial use
# require the copyright holders' consent.
"""Regenerate docs/ARCHITECTURE.png from docs/ARCHITECTURE.dot.

docs/ARCHITECTURE.dot is the source of truth for the module call and data-flow diagram.
This script draws that graph and writes docs/ARCHITECTURE.svg, then converts the SVG to
docs/ARCHITECTURE.png. Before drawing anything it checks its own layout table against
the nodes and edges the .dot file declares and fails if the two disagree, so the shipped
image cannot drift from its source. Every edge drawn must exist in the .dot. Two-way
pairs in the .dot (A -> B and B -> A) are drawn as one double-headed arrow whose label
carries a "down:" and an "up:" half; the script reports how many .dot edges that merging
accounts for, and which .dot edges, if any, it does not draw.

Usage
-----
    python3 scripts/render_architecture.py            # write the SVG and the PNG
    python3 scripts/render_architecture.py --no-png   # write the SVG only
    python3 scripts/render_architecture.py --keep-svg # write both, keeping the SVG

Requirements
------------
Python 3 standard library only. The SVG step is portable. The PNG step requires macOS,
because it converts the SVG with the system renderer /usr/bin/qlmanage; on any other
platform use --no-png and convert docs/ARCHITECTURE.svg with a tool of your choice.

The shipped docs/ARCHITECTURE.png was produced by this script through that system
renderer, not by Graphviz and not by matplotlib, neither of which was installed where the
release was built. With Graphviz available the same graph can be rendered directly from
its source instead:

    dot -Tpng docs/ARCHITECTURE.dot -o docs/ARCHITECTURE.png
"""

import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DOT = os.path.join(ROOT, "docs", "ARCHITECTURE.dot")
SVG = os.path.join(ROOT, "docs", "ARCHITECTURE.svg")
PNG = os.path.join(ROOT, "docs", "ARCHITECTURE.png")

W, H = 1860, 1860  # square: the macOS thumbnailer fits a square frame
CONTENT_H = 1320   # the graph itself occupies the top CONTENT_H pixels

# ---------------------------------------------------------------------------- #
# Bands: one row per time step. (title, y of the band label, y top, y bottom)
# ---------------------------------------------------------------------------- #
BANDS = [
    ("NO TIME STEP   declarations, input sizing, CSV helper", 62, 60, 170),
    ("SEASON   the program and its daily loop", 236, 234, 352),
    ("DAILY   one call per simulated day", 430, 428, 560),
    ("HOURLY   24 steps per day, with a nested fixed point", 636, 634, 766),
    ("HOURLY ENERGY-BALANCE SUB-STEP   once per Newton iterate", 842, 840, 972),
    ("ADAPTIVE SUB-STEP WITHIN ONE DAY", 1048, 1046, 1178),
]

# ---------------------------------------------------------------------------- #
# Nodes: id -> (x, y, w, h, band, lines)
# ---------------------------------------------------------------------------- #
NODES = {
    "TYPES":        (40, 86, 210, 70, 0, ["TYPES.f90", "working-precision kinds"]),
    "SHARED":       (280, 86, 330, 70, 0, ["SHARED.f90", "global parameter and state store"]),
    "count_lines":  (940, 86, 230, 70, 0, ["count_lines.f90", "counts input records"]),
    "csv_file":     (1200, 86, 200, 70, 0, ["csv_file.f90", "minimal CSV writer"]),
    "csv_file_1d":  (1430, 86, 240, 70, 0, ["csv_file_1d.f90", "include body of csv_file"]),

    "main_MMCO2Crop": (340, 258, 900, 88, 1,
                       ["main_MMCO2Crop.f90   -   PROGRAM mmco2crop",
                        "reads the decks, refuses retired options, runs the day loop,",
                        "writes the 88 output files and the season summary"]),

    "phenology_module": (180, 452, 260, 100, 2,
                         ["phenology_module.f90", "thermal time, vernalization,",
                          "stage dates, dvs 0-2"]),
    "growth_module":    (470, 452, 310, 100, 2,
                         ["growth_module.f90", "carbon and nitrogen economy,",
                          "organ canopy, yield formation"]),
    "soiln_module":     (1010, 452, 260, 100, 2,
                         ["soiln_module.f90", "residue and microbial", "C-N turnover"]),
    "mineral_module":   (1300, 452, 260, 100, 2,
                         ["mineral_module.f90", "six-element grain", "mineral economy"]),

    "meteo":        (180, 658, 220, 100, 3,
                     ["meteo.f90", "hourly weather from", "a daily record"]),
    "th":           (430, 658, 120, 100, 3, ["th.f90", "solar", "geometry"]),
    "solar":        (575, 658, 180, 100, 3, ["solar.f90", "two-leaf PAR", "partition"]),
    "CANOPY_MODEL": (790, 658, 250, 100, 3,
                     ["CANOPY_MODEL.f90", "canopy integration and",
                      "canopy-temperature solve"]),
    "assim":        (1070, 658, 220, 100, 3,
                     ["assim.f90", "leaf gas exchange,", "one leaf class"]),

    "energy_balance_module": (700, 864, 320, 100, 4,
                              ["energy_balance_module.f90",
                               "three-node canopy / soil /", "airspace energy balance"]),
    "radiation_module":      (1060, 864, 260, 100, 4,
                              ["radiation_module.f90", "component shortwave,",
                               "gray-body longwave"]),
    "soil_heat_module":      (1350, 864, 260, 100, 4,
                              ["soil_heat_module.f90", "1-D soil heat", "conduction"]),

    "hydro": (760, 1070, 380, 100, 5,
              ["hydro.f90   mixed-form Richards soil water,",
               "root uptake, infiltration, drainage,",
               "adaptive implicit sub-steps"]),
}

# ---------------------------------------------------------------------------- #
# Edges. Each entry:
#   (src, src_side, src_frac, dst, dst_side, dst_frac, bowx, bowy,
#    label, lab_dx, lab_dy, style, heads)
# side in l r t b; frac is the position along that side, 0..1.
# bowx / bowy displace the Bezier mid control point.
# style: solid | dashed ; heads: 1 (src -> dst) or 2 (double-headed)
# Every pair drawn with heads=2 stands for the two .dot edges named in DOT_PAIRS.
# ---------------------------------------------------------------------------- #
E = [
    ("TYPES", "r", 0.5, "SHARED", "l", 0.5, 0, 0, "wp, ik", 0, 52, "solid", 1),
    ("SHARED", "b", 0.5, "main_MMCO2Crop", "t", 0.10, 0, 0,
     "all parameters, switches and state\nread and written by 15 of the 20 units",
     -40, 0, "dashed", 2),
    ("count_lines", "b", 0.5, "main_MMCO2Crop", "t", 0.80, 0, 0, "NDays", 0, 0, "solid", 1),
    ("main_MMCO2Crop", "t", 0.95, "csv_file", "b", 0.5, 0, 0,
     "daily energy-\nbalance row", 0, 0, "solid", 1),
    ("csv_file", "r", 0.5, "csv_file_1d", "l", 0.5, 0, 0, "include", 0, 52, "dashed", 1),

    ("main_MMCO2Crop", "b", 0.06, "phenology_module", "t", 0.5, 0, 0,
     "down: Tmax, Tmin, DOY\nup: dvs 0-2, stage dates", 0, 0, "solid", 2),
    ("main_MMCO2Crop", "b", 0.24, "growth_module", "t", 0.5, 0, 0,
     "down: assimilate, dvs, N supply\nup: organ DM, LAI, N pools, grain", 0, 0, "solid", 2),
    ("main_MMCO2Crop", "b", 0.84, "soiln_module", "t", 0.5, 0, 0,
     "down: residue, Tsoil, theta\nup: mineralization / immobilization", 0, 0, "solid", 2),
    ("main_MMCO2Crop", "b", 0.97, "mineral_module", "t", 0.5, 0, 0,
     "down: organ dDM, Tr\nup: grain mg kg-1", 92, 26, "solid", 2),

    ("main_MMCO2Crop", "b", 0.639, "CANOPY_MODEL", "t", 0.5, 0, 0,
     "down: LAI, leaf N, theta, CO2\nup: An, Tr and Esoil demand,\nTcan, gsw, Ci", 0, 0, "solid", 2),
    ("main_MMCO2Crop", "r", 0.5, "hydro", "r", 0.45, 520, 0,
     "down: Tr and Esoil demand,\nrain, irrigation, root depth\nup: theta, Tr and E actual,\ndrainage, runoff",
     140, 330, "solid", 2),
    ("main_MMCO2Crop", "r", 0.88, "soil_heat_module", "r", 0.15, 700, -120,
     "initialize_soil_heat:\nlayer theta, Tdeep", -4, -78, "dashed", 1),

    ("CANOPY_MODEL", "l", 0.12, "th", "r", 0.12, 0, -90,
     "down: DOY, latitude\nup: daylength, declination", 0, -70, "solid", 2),
    ("CANOPY_MODEL", "b", 0.02, "meteo", "b", 0.5, 0, 96,
     "down: daily weather record\nup: Tair, VPD, SW down, sky emissivity", -80, -40, "solid", 2),
    ("meteo", "r", 0.25, "solar", "l", 0.25, 0, -46,
     "solar elevation,\ndirect / diffuse split", 0, -54, "solid", 1),
    ("solar", "r", 0.3, "CANOPY_MODEL", "l", 0.3, 0, -44,
     "Kdr, sunlit and shaded LAI,\nabsorbed PAR per class", 120, 150, "solid", 1),
    ("CANOPY_MODEL", "r", 0.45, "assim", "l", 0.45, 0, 46,
     "down: Tleaf, PAR, Narea, Ca, VPD\nup: An, A, gsw, Ci, N assimilated", 138, 58, "solid", 2),

    ("CANOPY_MODEL", "b", 0.5, "energy_balance_module", "t", 0.5, 0, 0,
     "down: Tcan, rcs, LAI, wind\nup: Tac, e_ac, Tsoil, Hc, Hs,\nLEc, LEs, raa/rca/rsa/rss",
     -230, 0, "solid", 2),
    ("energy_balance_module", "r", 0.35, "radiation_module", "l", 0.35, 0, -40,
     "down: Tcan, Tsoil, Kdr\nup: Rnc, Rns, G, rad_* ledger", 168, -56, "solid", 2),
    ("radiation_module", "r", 0.5, "soil_heat_module", "l", 0.5, 0, 0,
     "", 0, 0, "solid", 0),
    ("energy_balance_module", "b", 0.95, "soil_heat_module", "b", 0.3, 0, 76,
     "down: surface flux G, k_face\nup: Tsoil profile, Cv, lambda", 118, 88, "solid", 2),

    ("hydro", "t", 0.25, "energy_balance_module", "b", 0.5, 0, 0,
     "surface head, soil\nlatent-heat supply limit", -74, 0, "solid", 1),
    ("hydro", "t", 0.72, "CANOPY_MODEL", "b", 0.98, 0, 0,
     "soil water\npotential", 0, 96, "solid", 1),
    ("hydro", "l", 0.55, "growth_module", "b", 0.08, -620, 0,
     "water-stress state", -268, 14, "solid", 1),
    ("growth_module", "b", 0.03, "hydro", "l", 0.85, -500, 0,
     "root depth,\nTr demand", -305, -111, "solid", 1),
    ("growth_module", "r", 0.12, "CANOPY_MODEL", "t", 0.08, 40, -30,
     "LAI, leaf N", -60, 90, "solid", 1),
    ("growth_module", "b", 0.9, "mineral_module", "b", 0.3, 0, 72,
     "organ dDM, senescence DM", -44, 0, "solid", 1),
]

# the .dot edges that each double-headed arrow above stands for, and the one
# unlabelled connector (radiation_module -> soil_heat_module) that is drawn as a plain
# line because the .dot carries no such edge -- it is removed before the subset check.
NO_DOT = {("radiation_module", "soil_heat_module")}

FONT = "Helvetica, Arial, sans-serif"


# ---------------------------------------------------------------------------- #
# .dot parsing (only the forms this project's .dot uses)
# ---------------------------------------------------------------------------- #
def parse_dot(path):
    text = open(path, encoding="utf-8").read()
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    nodes, edges = set(), set()
    for line in text.split("\n"):
        line = line.strip()
        if not line or line.startswith("//"):
            continue
        m = re.match(r"^([A-Za-z_][A-Za-z_0-9]*)\s*->\s*([A-Za-z_][A-Za-z_0-9]*)", line)
        if m:
            edges.add((m.group(1), m.group(2)))
            if "dir=both" in line:
                edges.add((m.group(2), m.group(1)))
            continue
        m = re.match(r"^([A-Za-z_][A-Za-z_0-9]*)\s*\[label=", line)
        if m and m.group(1) not in ("graph", "node", "edge"):
            nodes.add(m.group(1))
    return nodes, edges


def anchor(node, side, frac):
    x, y, w, h, _band, _lines = NODES[node]
    if side == "l":
        return (x, y + h * frac, -1, 0)
    if side == "r":
        return (x + w, y + h * frac, 1, 0)
    if side == "t":
        return (x + w * frac, y, 0, -1)
    if side == "b":
        return (x + w * frac, y + h, 0, 1)
    raise ValueError(side)


def esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def main(argv):
    make_png = "--no-png" not in argv

    dot_nodes, dot_edges = parse_dot(DOT)

    # --- consistency checks against the source of truth ---------------------- #
    missing = dot_nodes - set(NODES)
    extra = set(NODES) - dot_nodes
    if missing or extra:
        print("FAIL: node sets differ.")
        if missing:
            print("  in .dot but not in this layout: %s" % sorted(missing))
        if extra:
            print("  in this layout but not in .dot: %s" % sorted(extra))
        return 2
    print("nodes: %d in .dot, %d in layout, sets identical" % (len(dot_nodes), len(NODES)))

    drawn = set()
    for e in E:
        s, _ss, _sf, d, _ds, _df = e[0], e[1], e[2], e[3], e[4], e[5]
        heads = e[12]
        if (s, d) in NO_DOT:
            continue
        drawn.add((s, d))
        if heads == 2:
            drawn.add((d, s))
    not_in_dot = drawn - dot_edges
    if not_in_dot:
        print("FAIL: edges drawn that the .dot does not declare: %s" % sorted(not_in_dot))
        return 2
    not_drawn = dot_edges - drawn
    print("edges: %d in .dot, %d represented by %d drawn arrows (%d double-headed)"
          % (len(dot_edges), len(drawn), len([x for x in E if x[12] != 0]),
             len([x for x in E if x[12] == 2])))
    if not_drawn:
        print("NOT drawn in the PNG (still in the .dot): %s" % sorted(not_drawn))
    else:
        print("every .dot edge is represented in the PNG")

    # --- SVG ---------------------------------------------------------------- #
    o = []
    o.append('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
             'viewBox="0 0 %d %d">' % (W, H, W, H))
    o.append('<defs>'
             '<marker id="a" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
             'markerHeight="7" orient="auto-start-reverse">'
             '<path d="M 0 0 L 10 5 L 0 10 z" fill="#333"/></marker>'
             '<marker id="ad" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
             'markerHeight="7" orient="auto-start-reverse">'
             '<path d="M 0 0 L 10 5 L 0 10 z" fill="#777"/></marker>'
             '</defs>')
    o.append('<rect width="%d" height="%d" fill="#fbfbf9"/>' % (W, H))
    o.append('<text x="%d" y="34" font-family="%s" font-size="22" font-weight="600" '
             'text-anchor="middle" fill="#111">M&amp;M-CO2Wheat v1.0 - module call and '
             'data-flow graph</text>' % (W // 2, FONT))

    # bands
    for title, ylab, y0, y1 in BANDS:
        o.append('<rect x="14" y="%d" width="%d" height="%d" rx="10" fill="#f0efe9" '
                 'stroke="#d8d6cc"/>' % (y0, W - 28, y1 - y0))
        o.append('<text x="26" y="%d" font-family="%s" font-size="11.5" '
                 'font-weight="700" letter-spacing="0.6" fill="#8a8575">%s</text>'
                 % (ylab + 14, FONT, esc(title)))

    # edges first, so boxes sit on top of the line ends
    for (s, ss, sf, d, ds, df, bx, by, lab, ldx, ldy, style, heads) in E:
        x1, y1, nx1, ny1 = anchor(s, ss, sf)
        x2, y2, nx2, ny2 = anchor(d, ds, df)
        c1 = (x1 + nx1 * 46 + bx * 0.5, y1 + ny1 * 46 + by * 0.5)
        c2 = (x2 + nx2 * 46 + bx * 0.5, y2 + ny2 * 46 + by * 0.5)
        col = "#777" if style == "dashed" else "#444"
        dash = ' stroke-dasharray="7 5"' if style == "dashed" else ""
        mk = "ad" if style == "dashed" else "a"
        start = ' marker-start="url(#%s)"' % mk if heads == 2 else ""
        end = ' marker-end="url(#%s)"' % mk if heads >= 1 else ""
        o.append('<path d="M %.1f %.1f C %.1f %.1f %.1f %.1f %.1f %.1f" fill="none" '
                 'stroke="%s" stroke-width="1.6"%s%s%s/>'
                 % (x1, y1, c1[0], c1[1], c2[0], c2[1], x2, y2, col, dash, start, end))
        if not lab:
            continue
        # label at the Bezier midpoint (t = 0.5), displaced by the given offset
        mx = 0.125 * x1 + 0.375 * c1[0] + 0.375 * c2[0] + 0.125 * x2 + ldx
        my = 0.125 * y1 + 0.375 * c1[1] + 0.375 * c2[1] + 0.125 * y2 + ldy
        rows = lab.split("\n")
        wpx = max(len(r) for r in rows) * 5.35 + 12
        hpx = len(rows) * 13 + 7
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="4" '
                 'fill="#fbfbf9" fill-opacity="0.93" stroke="#e2e0d6"/>'
                 % (mx - wpx / 2, my - hpx / 2, wpx, hpx))
        ty = my - hpx / 2 + 16
        for r in rows:
            o.append('<text x="%.1f" y="%.1f" font-family="%s" font-size="10.5" '
                     'text-anchor="middle" fill="#3a3a3a">%s</text>'
                     % (mx, ty, FONT, esc(r)))
            ty += 13

    # boxes
    for name, (x, y, w, h, _b, lines) in NODES.items():
        fill = "#ffffff"
        stroke = "#555"
        swid = 1.4
        if name in ("SHARED", "TYPES"):
            fill = "#eef2f7"
        if name == "main_MMCO2Crop":
            fill = "#e8eff1"
            stroke = "#2f4f5a"
            swid = 2.0
        o.append('<rect x="%d" y="%d" width="%d" height="%d" rx="8" fill="%s" '
                 'stroke="%s" stroke-width="%.1f"/>' % (x, y, w, h, fill, stroke, swid))
        ty = y + h / 2 - (len(lines) - 1) * 8 + 5
        for i, ln in enumerate(lines):
            bold = ' font-weight="700"' if i == 0 else ""
            size = 12.5 if i == 0 else 11
            o.append('<text x="%d" y="%.1f" font-family="%s" font-size="%.1f"%s '
                     'text-anchor="middle" fill="#111">%s</text>'
                     % (x + w // 2, ty, FONT, size, bold, esc(ln)))
            ty += 16

    # legend panel, in the square frame's spare space below the graph
    ly = CONTENT_H + 24
    o.append('<rect x="14" y="%d" width="%d" height="%d" rx="10" fill="#f6f5f0" '
             'stroke="#d8d6cc"/>' % (ly, W - 28, H - ly - 50))
    o.append('<text x="40" y="%d" font-family="%s" font-size="15" font-weight="700" '
             'fill="#333">How to read this diagram</text>' % (ly + 34, FONT))
    legend = [
        "Each grey band is one TIME STEP. A module sits in the band on which it is "
        "called, and the bands run from the top down in the order the program enters "
        "them: declarations, the season and its day loop,",
        "the daily crop and soil step, the hourly canopy step, the hourly "
        "energy-balance sub-step inside it, and the adaptive Richards sub-step that "
        "the day loop drives.",
        "",
        "Solid arrow: data flow, labelled with the main variables that cross it.",
        "Double-headed arrow: a two-way pair of flows, drawn once; its label is split "
        "into a \u201cdown:\u201d half (toward the lower module) and an \u201cup:\u201d half.",
        "Dashed arrow: the global store SHARED, or a textual include.",
        "Shaded boxes TYPES.f90 and SHARED.f90: declarations only, no process law.",
        "",
        "Source of truth: docs/ARCHITECTURE.dot (39 edges, 19 nodes). This PNG was "
        "rendered from that file by scripts/render_architecture.py, which fails if "
        "the two disagree. Per-module detail: docs/MODULES.md.",
    ]
    ty = ly + 68
    for line in legend:
        if line:
            o.append('<text x="40" y="%d" font-family="%s" font-size="12.5" '
                     'fill="#444">%s</text>' % (ty, FONT, esc(line)))
        ty += 24
    o.append("</svg>")

    with open(SVG, "w", encoding="utf-8") as fh:
        fh.write("\n".join(o))
    print("wrote %s (%d bytes)" % (SVG, os.path.getsize(SVG)))

    if not make_png:
        return 0
    tmpdir = os.path.join(os.path.dirname(SVG), "_qlpng")
    os.makedirs(tmpdir, exist_ok=True)
    subprocess.run(["/usr/bin/qlmanage", "-t", "-s", str(W), "-o", tmpdir, SVG],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    produced = os.path.join(tmpdir, os.path.basename(SVG) + ".png")
    if not os.path.exists(produced):
        print("FAIL: qlmanage produced no PNG")
        return 3
    os.replace(produced, PNG)
    os.rmdir(tmpdir)
    print("wrote %s (%d bytes)" % (PNG, os.path.getsize(PNG)))
    # The SVG is an intermediate on this path: docs/ARCHITECTURE.png is what the archive
    # ships and what MD5SUMS lists. Leaving the SVG behind made a post-regeneration
    # MD5SUMS check report one extra file, so it is removed unless the caller asks to
    # keep it (--keep-svg) or asked for the SVG alone (--no-png).
    if "--keep-svg" not in argv:
        os.remove(SVG)
        print("removed the intermediate %s (pass --keep-svg to keep it)" % SVG)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
