#!/usr/bin/env python3
"""
snakemap.py  –  Snakemake rulegraph → tube-map SVG (+ optional animated GIF)
=============================================================================

Usage
-----
# Pipe directly from snakemake (most common):
    snakemake -s rules/convert2fastq.smk --forceall --rulegraph \\
        | python3 snakemap.py -o convert2fastq_tubemap.svg

# From a saved DOT file:
    python3 snakemap.py -i rulegraph.dot -o pipeline.svg
    python3 snakemap.py -i rulegraph.dot -o pipeline.gif   # animated reveal

# Merge two pipelines side-by-side into one map (stacked vertically):
    python3 snakemap.py \\
        -i convert2fastq.dot -i fastqProcessing.dot \\
        --pipeline-title "BCL → FASTQ" \\
        --pipeline-title "FASTQ processing" \\
        -o combined_tubemap.svg

Conda requirements
------------------
    conda install networkx          # required (graph layout)
    conda install pillow            # only needed for --output *.gif
    conda install -c conda-forge cairosvg   # optional, needed for GIF
    # alternatives for SVG rasterisation: rsvg-convert or inkscape in PATH
"""

from __future__ import annotations

import argparse
import html
import re
import sys
import textwrap
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# Palette  (TfL-inspired tube-line colours)
# ─────────────────────────────────────────────────────────────────────────────
LINE_COLOURS = [
    "#E32017",  # Central (red)
    "#0098D4",  # Victoria (sky blue)
    "#007229",  # District (green)
    "#9B0056",  # Metropolitan (maroon)
    "#EE7C0E",  # Overground (orange)
    "#9B59B6",  # Elizabeth / violet
    "#F3A9BB",  # Hammersmith (pink)
    "#00A4A7",  # Jubilee-teal
    "#B36305",  # Bakerloo (brown)
    "#FFD300",  # Circle (yellow – used last, low contrast)
]

BG_COLOUR     = "#F5F5E8"   # classic tube-map parchment
STATION_FILL  = "#FFFFFF"
STATION_STROKE = "#1A1A1A"
FONT           = "Arial, Helvetica, sans-serif"

# ─── Layout constants ────────────────────────────────────────────────────────
TRACK_W    = 7
STATION_R  = 10
JUNC_R     = 14
COL_GAP    = 200
ROW_GAP    = 90
MARGIN_X   = 70
MARGIN_Y   = 85
LABEL_DX   = 15
LABEL_ANGLE = -40


# ─────────────────────────────────────────────────────────────────────────────
# 1. DOT parser
# ─────────────────────────────────────────────────────────────────────────────

def parse_dot(text: str, id_offset: int = 0) -> tuple[dict[str, str], list[tuple[str, str]]]:
    """
    Parse a Snakemake --rulegraph DOT string.
    Returns (nodes, edges).  node ids are string integers offset by id_offset.
    """
    nodes: dict[str, str] = {}
    edges: list[tuple[str, str]] = []

    for m in re.finditer(
        r"""^\s*(\d+)\s*\[.*?label\s*=\s*["\\]?\s*([^"\\,\]\n]+)""",
        text, re.M | re.I,
    ):
        nid   = str(int(m.group(1)) + id_offset)
        label = m.group(2).strip().strip('"').strip("\\").strip()
        if label:
            nodes[nid] = label

    for m in re.finditer(r"^\s*(\d+)\s*->\s*(\d+)", text, re.M):
        src = str(int(m.group(1)) + id_offset)
        dst = str(int(m.group(2)) + id_offset)
        edges.append((src, dst))

    if not nodes:
        sys.exit(
            "ERROR: no nodes found in the DOT input.\n"
            "Make sure you are passing output from:  snakemake --forceall --rulegraph"
        )
    return nodes, edges


# ─────────────────────────────────────────────────────────────────────────────
# 2. Layer assignment (longest-path)
# ─────────────────────────────────────────────────────────────────────────────

def assign_layers(nodes: dict, edges: list) -> dict[str, int]:
    import networkx as nx
    G = nx.DiGraph()
    G.add_nodes_from(nodes.keys())
    G.add_edges_from(edges)
    if not nx.is_directed_acyclic_graph(G):
        for cycle in list(nx.simple_cycles(G)):
            if len(cycle) >= 2:
                G.remove_edge(cycle[-1], cycle[0])
    layer: dict[str, int] = {}
    for n in nx.topological_sort(G):
        preds = list(G.predecessors(n))
        layer[n] = 0 if not preds else 1 + max(layer.get(p, 0) for p in preds)
    return layer


# ─────────────────────────────────────────────────────────────────────────────
# 3. Position assignment
# ─────────────────────────────────────────────────────────────────────────────

def assign_positions(
    nodes: dict,
    edges: list,
    layer: dict[str, int],
    y_base: float = 0.0,
) -> dict[str, tuple[float, float]]:
    import networkx as nx
    G = nx.DiGraph(edges)
    by_layer: dict[int, list[str]] = {}
    for nid, lyr in layer.items():
        by_layer.setdefault(lyr, []).append(nid)

    # Barycentre sort – forward pass
    for lyr in range(1, max(by_layer, default=0) + 1):
        if lyr not in by_layer:
            continue
        prev = {n: i for i, n in enumerate(by_layer.get(lyr - 1, []))}
        def bary(n: str) -> float:
            preds = [p for p in G.predecessors(n) if p in prev]
            return sum(prev[p] for p in preds) / len(preds) if preds else 1e9
        by_layer[lyr].sort(key=bary)

    max_col = max(len(v) for v in by_layer.values()) if by_layer else 1
    centre_y = y_base + MARGIN_Y + (max_col - 1) * ROW_GAP / 2

    pos: dict[str, tuple[float, float]] = {}
    for lyr, nids in by_layer.items():
        x = MARGIN_X + lyr * COL_GAP
        count = len(nids)
        for i, nid in enumerate(nids):
            y = centre_y + (i - (count - 1) / 2) * ROW_GAP
            pos[nid] = (x, y)
    return pos


# ─────────────────────────────────────────────────────────────────────────────
# 4. Colour assignment
# ─────────────────────────────────────────────────────────────────────────────

def assign_colours(
    nodes: dict, edges: list, palette_offset: int = 0
) -> dict[str, str]:
    import networkx as nx
    G = nx.DiGraph(edges)
    colour_map: dict[str, str] = {}
    comps = sorted(nx.weakly_connected_components(G), key=len, reverse=True)
    for ci, comp in enumerate(comps):
        sub = G.subgraph(comp)
        sources = [n for n in sub.nodes if sub.in_degree(n) == 0] or [next(iter(comp))]
        visited: set[str] = set()
        for si, src in enumerate(sources):
            c = LINE_COLOURS[(ci + si + palette_offset) % len(LINE_COLOURS)]
            for n in nx.bfs_tree(sub, src).nodes:
                if n not in visited:
                    colour_map[n] = c
                    visited.add(n)
        fallback = LINE_COLOURS[(ci + palette_offset) % len(LINE_COLOURS)]
        for n in comp:
            colour_map.setdefault(n, fallback)
    return colour_map


# ─────────────────────────────────────────────────────────────────────────────
# 5. SVG rendering
# ─────────────────────────────────────────────────────────────────────────────

def _wrap_label(label: str, max_chars: int = 14) -> list[str]:
    words = re.split(r"[_\-]", label)
    lines: list[str] = []
    cur = ""
    for w in words:
        candidate = (cur + " " + w).strip() if cur else w
        if cur and len(candidate) > max_chars:
            lines.append(cur)
            cur = w
        else:
            cur = candidate
    if cur:
        lines.append(cur)
    return lines or [label]


def _track_path(x1: float, y1: float, x2: float, y2: float) -> str:
    """
    Tube-map 45° routing:
    Horizontal → 45° diagonal → Horizontal
    """
    if abs(y1 - y2) < 4:
        return f"M {x1:.1f} {y1:.1f} L {x2:.1f} {y2:.1f}"
    dy   = abs(y2 - y1)
    dx   = x2 - x1
    horiz = max(0.0, (dx - dy) / 2)
    sign = 1 if y2 > y1 else -1
    mid_x1 = x1 + horiz
    mid_x2 = mid_x1 + dy
    return (
        f"M {x1:.1f} {y1:.1f} "
        f"L {mid_x1:.1f} {y1:.1f} "
        f"L {mid_x2:.1f} {y2:.1f} "
        f"L {x2:.1f} {y2:.1f}"
    )


def render_svg(
    nodes_all:    dict[str, str],
    edges_all:    list[tuple[str, str]],
    pos_all:      dict[str, tuple[float, float]],
    colours_all:  dict[str, str],
    titles:       list[str],
    pipeline_ids: list[set[str]],
) -> str:
    import networkx as nx
    G   = nx.DiGraph(edges_all)
    deg = {n: G.in_degree(n) + G.out_degree(n) for n in nodes_all}

    xs = [p[0] for p in pos_all.values()]
    ys = [p[1] for p in pos_all.values()]
    W  = int(max(xs) + MARGIN_X + 180)
    H  = int(max(ys) + MARGIN_Y + 70)

    combined_title = (
        " → ".join(titles) if len(titles) > 1
        else (titles[0] if titles else "Snakemake pipeline")
    )

    out: list[str] = []

    # ── SVG header ───────────────────────────────────────────────────────────
    out.append(f"""\
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     width="{W}" height="{H}"
     viewBox="0 0 {W} {H}"
     font-family="{FONT}">
  <title>{html.escape(combined_title)}</title>
  <desc>Tube-map style Snakemake pipeline diagram – snakemap.py</desc>
  <rect width="{W}" height="{H}" fill="{BG_COLOUR}"/>
""")

    # ── Title / pipeline banners ──────────────────────────────────────────────
    if len(titles) == 1:
        out.append(
            f'  <text x="22" y="34" font-size="17" font-weight="bold" fill="#1a1a1a">'
            f'{html.escape(titles[0])}</text>\n'
            f'  <text x="22" y="52" font-size="10" fill="#888">snakemap.py · tube-map view</text>\n'
        )
    else:
        for pi, (pids, ptitle) in enumerate(zip(pipeline_ids, titles)):
            p_ys = [pos_all[n][1] for n in pids if n in pos_all]
            if not p_ys:
                continue
            banner_y = min(p_ys) - 44
            colour = LINE_COLOURS[(pi * 3) % len(LINE_COLOURS)]
            # Coloured strip behind title
            out.append(
                f'  <rect x="16" y="{banner_y - 14:.0f}" width="14" height="20" '
                f'rx="3" fill="{colour}"/>\n'
                f'  <text x="36" y="{banner_y:.0f}" font-size="13" font-weight="bold" fill="#333">'
                f'{html.escape(ptitle)}</text>\n'
            )
        out.append(
            f'  <text x="22" y="24" font-size="9" fill="#999">snakemap.py · tube-map view</text>\n'
        )

    # ── Tracks (draw thicker coloured lines) ────────────────────────────────
    out.append("\n  <!-- tracks -->\n")
    for src, dst in edges_all:
        if src not in pos_all or dst not in pos_all:
            continue
        x1, y1 = pos_all[src]
        x2, y2 = pos_all[dst]
        c = colours_all.get(src, LINE_COLOURS[0])
        d = _track_path(x1, y1, x2, y2)
        out.append(
            f'  <path d="{d}" fill="none" stroke="{c}" '
            f'stroke-width="{TRACK_W}" stroke-linecap="round" stroke-linejoin="round"/>\n'
        )

    # ── Track highlight (white centreline, rail-like effect) ─────────────────
    out.append("\n  <!-- track highlights -->\n")
    for src, dst in edges_all:
        if src not in pos_all or dst not in pos_all:
            continue
        x1, y1 = pos_all[src]
        x2, y2 = pos_all[dst]
        d = _track_path(x1, y1, x2, y2)
        out.append(
            f'  <path d="{d}" fill="none" stroke="white" '
            f'stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round" opacity="0.5"/>\n'
        )

    # ── Stations ────────────────────────────────────────────────────────────
    out.append("\n  <!-- stations -->\n")
    for nid, label in nodes_all.items():
        if nid not in pos_all:
            continue
        x, y  = pos_all[nid]
        c     = colours_all.get(nid, LINE_COLOURS[0])
        d_val = deg.get(nid, 1)
        r     = JUNC_R if d_val > 3 else STATION_R
        is_src = G.in_degree(nid) == 0
        is_dst = G.out_degree(nid) == 0

        if is_src or is_dst:
            # Terminal: filled blob + small white inner dot
            out.append(
                f'  <circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" '
                f'fill="{c}" stroke="{STATION_STROKE}" stroke-width="2.5"/>\n'
                f'  <circle cx="{x:.1f}" cy="{y:.1f}" r="{r * 0.38:.1f}" '
                f'fill="white" stroke="none"/>\n'
            )
        elif d_val > 3:
            # Interchange: white disc + thick coloured ring + small inner circle
            out.append(
                f'  <circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" '
                f'fill="{STATION_FILL}" stroke="{c}" stroke-width="{TRACK_W}"/>\n'
            )
        else:
            # Standard: white disc + thin coloured ring
            out.append(
                f'  <circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" '
                f'fill="{STATION_FILL}" stroke="{c}" stroke-width="3.5"/>\n'
            )

        # ── Angled station label ─────────────────────────────────────────────
        lines  = _wrap_label(label)
        lx     = x + r + LABEL_DX
        ly     = y + 2
        out.append(f'  <g transform="translate({lx:.1f},{ly:.1f}) rotate({LABEL_ANGLE})">\n')
        for li, line_text in enumerate(lines):
            dy_px   = li * 14
            escaped = html.escape(line_text)
            # White halo for legibility over tracks
            out.append(
                f'    <text x="0" y="{dy_px}" font-size="11" font-weight="bold" '
                f'fill="white" stroke="white" stroke-width="4" paint-order="stroke">{escaped}</text>\n'
                f'    <text x="0" y="{dy_px}" font-size="11" font-weight="bold" fill="{c}">{escaped}</text>\n'
            )
        out.append("  </g>\n")

    # ── Legend ───────────────────────────────────────────────────────────────
    seen: dict[str, str] = {}
    for nid, c in colours_all.items():
        if c not in seen and nid in nodes_all:
            seen[c] = nodes_all[nid]

    leg_x0 = 20
    leg_y0 = H - 16 - len(seen) * 18
    out.append("\n  <!-- legend -->\n")
    for idx, (c, rep) in enumerate(seen.items()):
        ly = leg_y0 + idx * 18
        out.append(
            f'  <rect x="{leg_x0}" y="{ly}" width="30" height="9" rx="4.5" fill="{c}"/>\n'
            f'  <text x="{leg_x0 + 38}" y="{ly + 7.5}" font-size="9" fill="#555">'
            f'{html.escape(rep)}</text>\n'
        )

    out.append("</svg>\n")
    return "".join(out)


# ─────────────────────────────────────────────────────────────────────────────
# 6. Animated GIF (optional, requires Pillow + SVG rasteriser)
# ─────────────────────────────────────────────────────────────────────────────

def svg_to_gif(svg_text: str, out_path: Path, n_frames: int = 50) -> None:
    try:
        from PIL import Image
    except ImportError:
        sys.exit("Pillow not found.  Install with:  conda install pillow")

    import io, shutil, subprocess

    png_bytes: bytes | None = None

    # Try cairosvg first (best quality)
    try:
        import cairosvg  # type: ignore
        png_bytes = cairosvg.svg2png(bytestring=svg_text.encode())
    except ImportError:
        pass

    # Fall back to command-line tools
    if png_bytes is None:
        for cmd in (
            ["rsvg-convert", "--format=png"],
            ["inkscape", "--pipe", "--export-type=png"],
        ):
            if shutil.which(cmd[0]):
                r = subprocess.run(cmd, input=svg_text.encode(), capture_output=True)
                if r.returncode == 0 and r.stdout:
                    png_bytes = r.stdout
                    break

    if png_bytes is None:
        print(
            "WARNING: no SVG rasteriser found (cairosvg / rsvg-convert / inkscape).\n"
            "         GIF not created – SVG file is saved instead.",
            file=sys.stderr,
        )
        return

    base = Image.open(io.BytesIO(png_bytes)).convert("RGBA")
    W, H = base.size
    bg   = Image.new("RGBA", (W, H), (245, 245, 232, 255))

    frames: list[Image.Image] = []
    for i in range(n_frames):
        reveal_x = int(W * (i + 1) / n_frames)
        frame    = bg.copy()
        crop     = base.crop((0, 0, reveal_x, H))
        frame.paste(crop, (0, 0), crop)
        frames.append(frame.convert("RGB").quantize(colors=256, dither=Image.FLOYDSTEINBERG))

    for _ in range(25):           # hold final frame
        frames.append(frames[-1])

    frames[0].save(
        out_path,
        save_all=True, append_images=frames[1:],
        loop=0, duration=55, optimize=True,
    )
    print(f"GIF  →  {out_path}  ({len(frames)} frames, {W}×{H})", file=sys.stderr)


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        prog="snakemap.py",
        description="Convert Snakemake rulegraph DOT → tube-map SVG/GIF",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        Examples
        --------
          # Single pipeline piped from snakemake
          snakemake -s rules/convert2fastq.smk --forceall --rulegraph \\
              | python3 snakemap.py -o convert2fastq_tubemap.svg

          # From a saved DOT file
          python3 snakemap.py -i rulegraph.dot -o pipeline.svg

          # Animated GIF (requires Pillow + cairosvg or rsvg-convert)
          python3 snakemap.py -i rulegraph.dot -o pipeline.gif

          # Two pipelines stacked on one map
          python3 snakemap.py \\
              -i convert2fastq.dot -i fastqProcessing.dot \\
              --pipeline-title "BCL → FASTQ" \\
              --pipeline-title "FASTQ processing" \\
              -o combined_tubemap.svg

          # Wider spacing for crowded graphs
          python3 snakemap.py -i rulegraph.dot -o pipeline.svg --col-gap 240 --row-gap 100
        """),
    )
    ap.add_argument("-i", "--input",  action="append", default=[],
                    metavar="DOT",
                    help="Input DOT file (- for stdin). Repeat for multiple pipelines.")
    ap.add_argument("-o", "--output", default="tubemap.svg",
                    help="Output .svg or .gif (default: tubemap.svg)")
    ap.add_argument("--title", default="Snakemake pipeline",
                    help="Map title (single-pipeline mode)")
    ap.add_argument("--pipeline-title", action="append", default=[],
                    dest="pipeline_titles", metavar="TITLE",
                    help="Per-pipeline label (repeat once per -i)")
    ap.add_argument("--gif-frames", type=int, default=50,
                    help="GIF animation frames (default: 50)")
    ap.add_argument("--col-gap", type=int, default=COL_GAP,
                    help=f"Horizontal layer spacing px (default: {COL_GAP})")
    ap.add_argument("--row-gap", type=int, default=ROW_GAP,
                    help=f"Vertical station spacing px (default: {ROW_GAP})")
    args = ap.parse_args()

    # Apply spacing overrides – update the module-level constants used by layout functions
    import sys as _sys
    _this = _sys.modules[__name__]
    _this.COL_GAP = args.col_gap
    _this.ROW_GAP = args.row_gap

    inputs = args.input or ["-"]
    dot_texts: list[str] = []
    for src in inputs:
        dot_texts.append(sys.stdin.read() if src == "-" else Path(src).read_text())

    all_nodes:    dict[str, str]            = {}
    all_edges:    list[tuple[str, str]]     = []
    all_colours:  dict[str, str]            = {}
    all_pos:      dict[str, tuple[float, float]] = {}
    titles:       list[str]                 = []
    pipeline_ids: list[set[str]]            = []

    id_offset = 0
    y_cursor  = 0.0

    for pi, dot_text in enumerate(dot_texts):
        if not dot_text.strip():
            continue
        nodes, edges = parse_dot(dot_text, id_offset=id_offset)
        print(f"Pipeline {pi+1}: {len(nodes)} rules, {len(edges)} deps", file=sys.stderr)
        id_offset += len(nodes) + 1

        # Pipeline title
        if pi < len(args.pipeline_titles):
            ptitle = args.pipeline_titles[pi]
        elif len(dot_texts) == 1:
            ptitle = args.title
        else:
            # Guess from terminal (sink) nodes
            import networkx as nx
            G_tmp = nx.DiGraph(edges)
            sinks = [nodes[n] for n in nodes if G_tmp.out_degree(n) == 0]
            ptitle = f"Pipeline {pi + 1}" + (f" [{sinks[0]}]" if sinks else "")
        titles.append(ptitle)

        layer = assign_layers(nodes, edges)
        pos   = assign_positions(nodes, edges, layer, y_base=y_cursor)
        cols  = assign_colours(nodes, edges, palette_offset=pi * 3)

        # Advance y_cursor past this pipeline's bounding box + gap
        if pos:
            bottom = max(v[1] for v in pos.values())
            y_cursor = bottom + ROW_GAP * 2.2

        pipeline_ids.append(set(nodes.keys()))
        all_nodes.update(nodes)
        all_edges.extend(edges)
        all_colours.update(cols)
        all_pos.update(pos)

    if not all_nodes:
        sys.exit("ERROR: no valid pipeline data found in input(s)")

    svg_text = render_svg(
        all_nodes, all_edges, all_pos, all_colours, titles, pipeline_ids
    )

    out      = Path(args.output)
    svg_path = out if out.suffix.lower() == ".svg" else out.with_suffix(".svg")
    svg_path.write_text(svg_text, encoding="utf-8")
    print(f"SVG  →  {svg_path}", file=sys.stderr)

    if out.suffix.lower() == ".gif":
        svg_to_gif(svg_text, out, n_frames=args.gif_frames)


if __name__ == "__main__":
    main()