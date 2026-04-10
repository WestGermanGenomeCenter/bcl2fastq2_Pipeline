#!/usr/bin/env python3
"""
illumina_run_report.py
======================
Generates a multi-page PDF QC report from an Illumina runfolder.
Supports NextSeq 2000 and NovaSeq X Plus.

Usage:
    python illumina_run_report.py <runfolder> <output.pdf>

Conda packages:
    conda install -c conda-forge matplotlib pandas numpy reportlab lxml
"""

import argparse
import os
import re
import sys
import warnings
import xml.etree.ElementTree as ET
from datetime import datetime
from io import BytesIO
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib.enums import TA_CENTER
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    Image as RLImage, PageBreak, HRFlowable,
)

# ── Colour palette ────────────────────────────────────────────────────────────
C = {
    "blue":   "#1f77b4",
    "orange": "#ff7f0e",
    "green":  "#2ca02c",
    "red":    "#d62728",
    "purple": "#9467bd",
    "grey":   "#7f7f7f",
    "teal":   "#17becf",
    "brown":  "#8c564b",
}
SURFACE_COLORS = {"Top": C["blue"], "Bottom": C["orange"]}
SWATH_COLORS   = [C["blue"], C["orange"], C["green"], C["purple"],
                  C["teal"], C["brown"], C["red"], "#e377c2"]

RL_ACCENT   = colors.HexColor("#1f77b4")
RL_LIGHT_BG = colors.HexColor("#f0f4f8")
RL_DARK     = colors.HexColor("#1a1a2e")
RL_GREY     = colors.HexColor("#cccccc")
RL_PASS     = colors.HexColor("#2ca02c")
RL_FAIL     = colors.HexColor("#d62728")

PAGE_W, PAGE_H = A4
MARGIN = 1.8 * cm


# ═══════════════════════════════════════════════════════════════════════════════
# File discovery
# ═══════════════════════════════════════════════════════════════════════════════

def find_file(runfolder: Path, patterns: list) -> "Path | None":
    """Search runfolder root and one level of subdirectories."""
    candidates = [runfolder] + [d for d in runfolder.iterdir() if d.is_dir()]
    for d in candidates:
        for pat in patterns:
            hits = sorted(d.glob(pat))
            if hits:
                return hits[0]
    return None


def warn_skip(name: str, reason: str = "file not found"):
    warnings.warn(f"[SKIP] {name}: {reason}", UserWarning, stacklevel=3)


# ═══════════════════════════════════════════════════════════════════════════════
# Instrument detection
# ═══════════════════════════════════════════════════════════════════════════════

def detect_instrument(runfolder: Path, recipe: "dict | None",
                      run_info: "dict | None") -> tuple:
    """
    Returns (model_string, serial_string).
    Only RunInfo.xml / RunParameters.xml is used as a reliable model source.
    Serial is extracted from the runfolder name but model is NOT inferred
    from prefix (too unreliable across firmware/naming conventions).
    """
    # Serial from runfolder name (YYMMDD_<SERIAL>_RUN_FLOWCELL)
    serial = "N/A"
    parts  = runfolder.name.split("_")
    if len(parts) >= 2:
        serial = parts[1]

    # RunInfo / RunParameters XML is the only reliable source for model
    if run_info:
        if run_info.get("serial"):
            serial = run_info["serial"]
        if run_info.get("instrument_type"):
            return run_info["instrument_type"], serial

    # No reliable instrument model source available
    return "N/A", serial


def load_run_info(runfolder: Path) -> "dict | None":
    """Parse RunInfo.xml or RunParameters.xml for instrument/run metadata."""
    for fname in ("RunInfo.xml", "RunParameters.xml"):
        path = find_file(runfolder, [fname])
        if path is None:
            continue
        try:
            root = ET.parse(path).getroot()
            info = {}
            for tag in ("Instrument", "InstrumentId", "InstrumentSerialNumber",
                        "SerialNumber"):
                el = root.find(f".//{tag}")
                if el is not None and el.text:
                    info["serial"] = el.text.strip()
            for tag in ("InstrumentType", "ApplicationName", "Application",
                        "InstrumentName"):
                el = root.find(f".//{tag}")
                if el is not None and el.text:
                    info["instrument_type"] = el.text.strip()
            for tag in ("FlowcellId", "FlowCellId", "Flowcell"):
                el = root.find(f".//{tag}")
                if el is not None and el.text:
                    info["flowcell"] = el.text.strip()
            reads = []
            for r in root.findall(".//Read"):
                nc = r.attrib.get("NumCycles") or r.attrib.get("Cycles", "?")
                reads.append({"cycles": nc,
                               "indexed": r.attrib.get("IsIndexedRead","N") == "Y"})
            if reads:
                info["reads"] = reads
            if info:
                return info
        except Exception:
            continue
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# Data loaders
# ═══════════════════════════════════════════════════════════════════════════════

def load_primary_metrics(runfolder: Path) -> "pd.DataFrame | None":
    path = find_file(runfolder, ["PrimaryAnalysisMetrics.csv"])
    if path is None:
        warn_skip("PrimaryAnalysisMetrics"); return None
    try:
        df = pd.read_csv(path)
        df.columns = df.columns.str.strip()
        df["Value"]  = pd.to_numeric(df["Value"], errors="coerce")
        df["Unit"]   = df["Unit"].astype(str).str.strip()
        df["Metric"] = df["Metric"].astype(str).str.strip()
        return df
    except Exception as e:
        warn_skip("PrimaryAnalysisMetrics", str(e)); return None


def load_xy_stage(runfolder: Path) -> "pd.DataFrame | None":
    path = find_file(runfolder, ["XYStagePositionErrorsReport.csv"])
    if path is None:
        warn_skip("XYStagePositionErrorsReport"); return None
    try:
        df = pd.read_csv(path)
        df.columns = df.columns.str.strip()
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        return df.dropna(subset=["Cycle"])
    except Exception as e:
        warn_skip("XYStagePositionErrorsReport", str(e)); return None


def load_autofocus(runfolder: Path) -> "pd.DataFrame | None":
    path = find_file(runfolder, ["AutofocusReport.csv"])
    if path is None:
        warn_skip("AutofocusReport"); return None
    try:
        df = pd.read_csv(path)
        df.columns = df.columns.str.strip()

        # Force-coerce every column except Surface
        for col in df.columns:
            if col != "Surface":
                df[col] = pd.to_numeric(df[col], errors="coerce")

        # Replace sentinel 99999 with NaN
        for col in df.columns:
            if "Stability" in col or "Stallguard" in col:
                df[col] = df[col].where(df[col] < 90000, np.nan)

        # z_clamped is boolean (0/1 or True/False) – not useful as histogram
        if "z_clamped" in df.columns:
            df.drop(columns=["z_clamped"], inplace=True)

        # Seeding z: nm → µm, sentinel -1 → NaN
        if "seeding Target z_nm" in df.columns:
            z_raw = df["seeding Target z_nm"]
            df["seeding_z_um"] = z_raw.where(z_raw > 0, np.nan) / 1000.0

        # X/Y stage position: nm → mm
        if "x_nm" in df.columns:
            df["x_mm"] = df["x_nm"] / 1e6
        if "y_nm" in df.columns:
            df["y_mm"] = df["y_nm"] / 1e6

        return df
    except Exception as e:
        warn_skip("AutofocusReport", str(e)); return None


def load_autocenter(runfolder: Path) -> "pd.DataFrame | None":
    path = find_file(runfolder, ["*AutoCenterData.csv"])
    if path is None:
        warn_skip("AutoCenterData"); return None
    try:
        rows, header = [], None
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line.strip():
                    continue
                if line.startswith("Timestamp") and header is None:
                    header = [c.strip() for c in line.split(",")]
                    continue
                if line.startswith("Timestamp"):
                    continue
                if header and not line.startswith("AutoCenter"):
                    fields = [f.strip() for f in line.split(",")]
                    if len(fields) == len(header):
                        rows.append(fields)
        if not rows or header is None:
            warn_skip("AutoCenterData", "no parseable rows"); return None
        df = pd.DataFrame(rows, columns=header)
        for col in df.columns:
            if col in ("Timestamp", "Surface"):
                continue
            df[col] = pd.to_numeric(
                df[col].astype(str).str.replace(r"\s*mm", "", regex=True),
                errors="coerce")
        return df
    except Exception as e:
        warn_skip("AutoCenterData", str(e)); return None


def load_autotilt(runfolder: Path) -> "dict | None":
    path = find_file(runfolder, ["*AutotiltData.csv"])
    if path is None:
        warn_skip("AutotiltData"); return None
    try:
        with open(path) as fh:
            content = fh.read()
        result = {}

        # ZCam positions
        zcam_rows, in_zcam = [], False
        for line in content.splitlines():
            if line.strip() == "ZCam Calculation Summary":
                in_zcam = True; continue
            if in_zcam:
                if line.startswith("Position"):
                    continue
                m = re.match(r"^(\d+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)", line)
                if m:
                    zcam_rows.append({"Position": int(m.group(1)),
                                      "X_mm": float(m.group(2)),
                                      "Y_mm": float(m.group(3))})
                elif zcam_rows:
                    break
        if zcam_rows:
            result["zcam"] = pd.DataFrame(zcam_rows)

        # Tilt angles
        m2 = re.search(
            r"ThetaX\(Milliradians\),ThetaY\(Milliradians\)\s*\n([-\d.]+),([-\d.]+)",
            content)
        if m2:
            result["theta"] = {"x_mrad": float(m2.group(1)),
                               "y_mrad": float(m2.group(2))}

        # Surface positions (all occurrences)
        surf_rows = []
        for line in content.splitlines():
            m3 = re.match(
                r"X = ([\d.]+), Y = ([\d.]+), Z_afm = ([\d.]+), "
                r"Z_cam = ([+-]?[\d.]+), Z = ([\d.]+)", line.strip())
            if m3:
                surf_rows.append({
                    "X_mm": float(m3.group(1)), "Y_mm": float(m3.group(2)),
                    "Z_afm_um": float(m3.group(3)),
                    "Z_cam_delta_um": float(m3.group(4)),
                    "Z_corrected_um": float(m3.group(5))})
        if surf_rows:
            result["surfaces_pre"]  = pd.DataFrame(surf_rows[:4])
            if len(surf_rows) > 4:
                result["surfaces_post"] = pd.DataFrame(surf_rows[4:8])

        # Motor corrections
        motor_rows = []
        for line in content.splitlines():
            m4 = re.match(
                r"(\d+),([\d. ]+mm),([\d. ]+mm),([\d.]+), ?([-\d.]+),([\d.]+),([\d.]+)",
                line.strip())
            if m4:
                motor_rows.append({
                    "Motor": int(m4.group(1)),
                    "Start_um": float(m4.group(4)),
                    "Delta_um": float(m4.group(5)),
                    "Target_um": float(m4.group(6)),
                    "End_um": float(m4.group(7))})
        if motor_rows:
            result["motors"] = pd.DataFrame(motor_rows)

        # Final residual tilt result
        m5 = re.search(
            r"Residual Tilt \[um\], Allowed Spec \[um\], Result\s*\n"
            r"([\d.]+),([\d.]+),(\w+)", content)
        if m5:
            result["tilt_result"] = {"residual_um": float(m5.group(1)),
                                     "spec_um":     float(m5.group(2)),
                                     "result":      m5.group(3)}

        return result if result else None
    except Exception as e:
        warn_skip("AutotiltData", str(e)); return None


def load_recipe(runfolder: Path) -> "dict | None":
    path = find_file(runfolder, ["EffectiveRecipe.xml"])
    if path is None:
        warn_skip("EffectiveRecipe.xml"); return None
    try:
        root = ET.parse(path).getroot()
        info = {"name":        root.attrib.get("Name", "N/A"),
                "version":     root.attrib.get("Version", "N/A"),
                "description": root.attrib.get("Description", "N/A")}

        reagents, seen_r = [], set()
        for r in root.findall(".//Reagents"):
            k = (r.attrib.get("Name",""), r.attrib.get("Type",""))
            if k not in seen_r:
                seen_r.add(k); reagents.append(r.attrib)
        info["reagents"] = reagents

        protocols, seen_p = [], set()
        for p in root.iter("Protocol"):
            n = p.attrib.get("Name","")
            if n and n not in seen_p:
                seen_p.add(n)
                protocols.append({"name": n,
                                   "type": p.attrib.get("ProtocolType","")})
        info["protocols"] = protocols

        primers_el = root.find("Primers")
        if primers_el is not None:
            info["primers"] = dict(primers_el.attrib)
            info["default_primers"] = {}
            for child in primers_el:
                if child.text:
                    info["default_primers"][child.tag] = child.text.strip()

        temps, seen_t = [], set()
        for el in root.findall(".//SetThermalZoneTemp"):
            if el.attrib.get("Enable","true").lower() != "true":
                continue
            key = (el.attrib.get("Zone",""), el.attrib.get("Temperature",""))
            if key not in seen_t:
                seen_t.add(key)
                temps.append({"Zone": el.attrib.get("Zone",""),
                               "Temp_C": el.attrib.get("Temperature",""),
                               "RampRate_C_s": el.attrib.get("RampRate","")})
        info["temperatures"] = temps

        for el in root.findall(".//FocusModelGeneration"):
            info["focus_model_freq"] = el.attrib.get("CycleFrequency","N/A")
            break

        # Human-readable recipe steps: one row per unique ChemistryDefinition
        # grouped under the protocol they belong to, deduped across protocols
        skip_step_names = {"Start", "Safe State", "Anti Prime"}
        recipe_steps = []
        seen_steps = set()
        for proto in root.findall(".//Protocol"):
            pname = proto.attrib.get("Name","")
            ptype = proto.attrib.get("ProtocolType","")
            if ptype in ("SafeState","PreFlowCheck","PostFlowCheck",
                         "CustomPrimerOneTransfer","CustomPrimerTwoTransfer"):
                continue
            for cd in proto.findall(".//ChemistryDefinition"):
                cdname   = cd.attrib.get("Name","")
                cdver    = cd.attrib.get("Version","")
                if not cdname or cdname in skip_step_names:
                    continue
                key = (pname, cdname)
                if key in seen_steps:
                    continue
                seen_steps.add(key)
                # Reagents
                reagents = []
                for el in cd.iter():
                    rn = el.attrib.get("ReagentName","")
                    if rn and rn not in reagents:
                        reagents.append(rn)
                # Temps (unique non-zero)
                temps = []
                for el in cd.findall(".//SetThermalZoneTemp"):
                    if el.attrib.get("Enable","true").lower() == "true":
                        t = el.attrib.get("Temperature","")
                        if t and t not in temps:
                            temps.append(t + "°C")
                # Waits (ms -> s, unique)
                waits = []
                for el in cd.findall(".//Wait"):
                    dur = el.attrib.get("Duration","")
                    if dur:
                        s = str(int(dur)//1000) + "s"
                        if s not in waits:
                            waits.append(s)
                recipe_steps.append({
                    "Protocol":  pname,
                    "Step":      cdname,
                    "Version":   cdver if cdver else "—",
                    "Reagents":  ", ".join(reagents[:8]) if reagents else "—",
                    "Temps":     ", ".join(temps) if temps else "—",
                    "Waits":     ", ".join(waits[:4]) if waits else "—",
                })
        if recipe_steps:
            info["recipe_steps"] = recipe_steps

        return info
    except Exception as e:
        warn_skip("EffectiveRecipe.xml", str(e)); return None


# ═══════════════════════════════════════════════════════════════════════════════
# Plotting helpers
# ═══════════════════════════════════════════════════════════════════════════════

def fig_to_rl(fig, dpi: int = 150) -> RLImage:
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight",
                facecolor="white")
    buf.seek(0)
    plt.close(fig)
    return RLImage(buf)


def clean_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


def embed(story, img, aspect=0.38):
    """Scale img to full page width and append to story."""
    if img is None:
        return
    w = PAGE_W - 2 * MARGIN
    img.drawWidth  = w
    img.drawHeight = w * aspect
    story.append(img)
    story.append(Spacer(1, 0.3 * cm))


# ═══════════════════════════════════════════════════════════════════════════════
# Plot functions
# ═══════════════════════════════════════════════════════════════════════════════

# ── Primary metrics ───────────────────────────────────────────────────────────

def plot_primary_metrics(df: pd.DataFrame) -> "RLImage | None":
    df = df.dropna(subset=["Value"]).reset_index(drop=True)
    n  = len(df)
    if n == 0:
        return None

    fig, axes = plt.subplots(1, n, figsize=(min(3.2 * n, 13), 3.0))
    if n == 1:
        axes = [axes]

    for ax, (_, row) in zip(axes, df.iterrows()):
        metric, value, unit = row["Metric"], row["Value"], row["Unit"]
        col = C["blue"]
        if "Q30" in metric:
            col = C["green"] if value >= 80 else C["red"]
        elif "Loading" in metric:
            col = C["green"] if 80 <= value <= 120 else C["orange"]

        ax.barh([0], [value], color=col, alpha=0.85, height=0.5)
        if "Q30" in metric:
            ax.axvline(80, color=C["red"], linestyle="--", linewidth=1.2,
                       label="80% min spec")
            ax.legend(fontsize=7, loc="lower right")
        ax.set_yticks([])
        ax.set_xlabel(f"{unit}  ({'higher is better' if 'Q30' in metric or 'Yield' in metric or 'Reads' in metric else 'should be ~100%'})",
                      fontsize=7.5)
        ax.set_title(f"{metric}\n{value:g} {unit}", fontsize=8.5,
                     fontweight="bold", pad=3)
        clean_ax(ax)
        ax.spines["left"].set_visible(False)

    fig.suptitle("Primary Analysis Metrics  —  run output quality at a glance",
                 fontsize=10, fontweight="bold")
    fig.tight_layout()
    return fig_to_rl(fig)


# ── XY stage ─────────────────────────────────────────────────────────────────

def plot_xy_stage(df: pd.DataFrame) -> "RLImage | None":
    cycles = df["Cycle"].values.astype(float)
    ex     = df["AveragePositionErrorXum"].values.astype(float)
    ey     = df["AveragePositionErrorYum"].values.astype(float)
    spec   = 2.0

    fig = plt.figure(figsize=(13, 4.8))
    gs  = gridspec.GridSpec(2, 3, figure=fig,
                            width_ratios=[3, 3, 1.2], hspace=0.50, wspace=0.35)

    for row_i, (err, ylabel, col) in enumerate([
        (ex, "Error X [µm]\nlateral stage deviation from target", C["blue"]),
        (ey, "Error Y [µm]\nlateral stage deviation from target", C["orange"]),
    ]):
        ax_ts = fig.add_subplot(gs[row_i, 0:2])
        ax_ts.fill_between(cycles, err, 0, alpha=0.2, color=col)
        ax_ts.plot(cycles, err, color=col, linewidth=1.2)
        ax_ts.axhline( spec, color=C["red"], linestyle="--", linewidth=1.0,
                       label=f"±{spec} µm spec")
        ax_ts.axhline(-spec, color=C["red"], linestyle="--", linewidth=1.0)
        ax_ts.axhline(0, color=C["grey"], linewidth=0.5)
        ax_ts.set_ylabel(ylabel, fontsize=7.5)
        ax_ts.set_xlabel("Cycle  [sequencing cycle number]", fontsize=7.5)
        ax_ts.set_xlim(cycles[0], cycles[-1])
        ax_ts.legend(fontsize=7)
        clean_ax(ax_ts)

        ax_d = fig.add_subplot(gs[row_i, 2])
        valid = err[np.isfinite(err)].astype(float)
        if len(valid) > 1:
            ax_d.hist(valid, bins=20, color=col, alpha=0.8,
                      orientation="horizontal", edgecolor="white")
        ax_d.axhline( spec, color=C["red"], linestyle="--", linewidth=1.0)
        ax_d.axhline(-spec, color=C["red"], linestyle="--", linewidth=1.0)
        ax_d.axhline(0, color=C["grey"], linewidth=0.5)
        ax_d.set_xlabel("Count  [number of cycles in bin]", fontsize=7.5)
        ax_d.tick_params(labelleft=False)
        clean_ax(ax_d)

    fig.suptitle("XY Stage Position Errors per Cycle  —  optical stage mechanical accuracy",
                 fontsize=10, fontweight="bold")
    return fig_to_rl(fig, dpi=150)


# ── Autofocus: position error scatter ────────────────────────────────────────

def plot_af_scatter(df: pd.DataFrame) -> "RLImage | None":
    """Per-tile X/Y position error scatter, coloured by swath."""
    surfaces = sorted(df["Surface"].dropna().unique())
    swaths   = sorted(df["Swath"].dropna().unique()) if "Swath" in df.columns else [1]
    ncols    = len(surfaces)

    fig, axes = plt.subplots(1, ncols, figsize=(8 * ncols, 6.5))
    if ncols == 1:
        axes = [axes]

    for ax, surface in zip(axes, surfaces):
        sub = df[df["Surface"] == surface]
        for si, sw in enumerate(swaths):
            ssub = sub[sub["Swath"] == sw] if "Swath" in sub.columns else sub
            if ssub.empty:
                continue
            ax.scatter(ssub["PositionErrorX_um"], ssub["PositionErrorY_um"],
                       s=5, alpha=0.25,
                       color=SWATH_COLORS[si % len(SWATH_COLORS)],
                       label=f"Swath {int(sw)}", rasterized=True)
        ax.axhline(0, color=C["grey"], linewidth=0.5)
        ax.axvline(0, color=C["grey"], linewidth=0.5)
        circ = plt.Circle((0, 0), 2.0, color=C["red"], fill=False,
                           linestyle="--", linewidth=1.0)
        ax.add_patch(circ)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Position Error X  [µm — tile X miss-registration]", fontsize=8)
        ax.set_ylabel("Position Error Y  [µm — tile Y miss-registration]", fontsize=8)
        ax.set_title(f"{surface} Surface  (n={len(sub):,})", fontsize=9,
                     fontweight="bold")
        ax.legend(fontsize=7, markerscale=2, loc="upper right")
        clean_ax(ax)

    fig.suptitle(
        "Per-Tile Position Errors  \u2014  dashed circle = 2 \u00b5m spec radius\n"
        "Swath = one imaging strip along the flow cell; tiles are sub-regions within each swath",
        fontsize=10, fontweight="bold")
    try:
        fig.tight_layout()
    except Exception:
        pass
    return fig_to_rl(fig, dpi=150)


# ── Autofocus: position errors over cycles ───────────────────────────────────

def plot_af_errors_cycle(df: pd.DataFrame) -> "RLImage | None":
    surfaces = sorted(df["Surface"].dropna().unique())
    swaths   = sorted(df["Swath"].dropna().unique()) if "Swath" in df.columns else [1]
    lanes    = sorted(df["Lane"].dropna().unique()) if "Lane" in df.columns else [1]

    n_rows, n_cols = 2, len(surfaces)
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(7 * n_cols, 5.5),
                             sharex=True, sharey="row")
    if n_cols == 1:
        axes = axes.reshape(2, 1)

    for ci, surface in enumerate(surfaces):
        sub = df[df["Surface"] == surface]
        for ri, (ecol, ylabel) in enumerate([
            ("PositionErrorX_um",
             "Pos. Error X  [µm — X deviation per tile, mean ± std]"),
            ("PositionErrorY_um",
             "Pos. Error Y  [µm — Y deviation per tile, mean ± std]"),
        ]):
            ax = axes[ri, ci]
            for si, sw in enumerate(swaths):
                ssub = sub[sub["Swath"] == sw] if "Swath" in sub.columns else sub
                grp  = ssub.groupby("Cycle")[ecol].agg(["mean","std"])
                if grp.empty:
                    continue
                col = SWATH_COLORS[si % len(SWATH_COLORS)]
                ax.plot(grp.index, grp["mean"], color=col, linewidth=1.2,
                        label=f"Swath {int(sw)}")
                ax.fill_between(grp.index,
                                grp["mean"] - grp["std"],
                                grp["mean"] + grp["std"],
                                color=col, alpha=0.15)
            ax.axhline(0, color=C["grey"], linewidth=0.5)
            ax.axhline( 2, color=C["red"], linestyle="--", linewidth=0.8)
            ax.axhline(-2, color=C["red"], linestyle="--", linewidth=0.8)
            ax.set_ylabel(ylabel, fontsize=8)
            if ri == n_rows - 1:
                ax.set_xlabel("Cycle  [sequencing cycle number]", fontsize=8)
            ax.set_title(f"{surface} Surface", fontsize=9, fontweight="bold")
            if ri == 0 and ci == n_cols - 1:
                ax.legend(fontsize=7, loc="upper right")
            clean_ax(ax)

    fig.suptitle(
        "Position Errors over Cycles  —  mean ± std per swath "
        "(dashed = ±2 µm spec)",
        fontsize=10, fontweight="bold")
    fig.tight_layout()
    return fig_to_rl(fig, dpi=150)


# ── Autofocus: seeding Z ──────────────────────────────────────────────────────

def plot_af_seeding_z(df: pd.DataFrame) -> "RLImage | None":
    if "seeding_z_um" not in df.columns:
        return None
    sub = df[df["seeding_z_um"].notna()].copy()
    if sub.empty:
        return None

    surfaces = sorted(sub["Surface"].dropna().unique())
    swaths   = sorted(sub["Swath"].dropna().unique()) if "Swath" in sub.columns else [1]
    lanes    = sorted(sub["Lane"].dropna().unique()) if "Lane" in sub.columns else [1]

    fig, axes = plt.subplots(1, len(surfaces),
                             figsize=(7 * len(surfaces), 4.0),
                             sharex=True)
    if len(surfaces) == 1:
        axes = [axes]

    for ax, surface in zip(axes, surfaces):
        ssub = sub[sub["Surface"] == surface]
        for li, lane in enumerate(lanes):
            lsub = ssub[ssub["Lane"] == lane] if "Lane" in ssub.columns else ssub
            for si, sw in enumerate(swaths):
                wsub = lsub[lsub["Swath"] == sw] if "Swath" in lsub.columns else lsub
                grp  = wsub.groupby("Cycle")["seeding_z_um"].mean()
                if grp.empty:
                    continue
                col   = SWATH_COLORS[si % len(SWATH_COLORS)]
                ls_   = ["-","--",":","-."][li % 4]
                lbl   = (f"L{int(lane)} Sw{int(sw)}" if len(lanes) > 1
                         else f"Swath {int(sw)}")
                ax.plot(grp.index, grp.values, color=col, linestyle=ls_,
                        linewidth=1.3, alpha=0.9, label=lbl)

        ax.set_xlabel("Cycle  [sequencing cycle number]", fontsize=8)
        ax.set_ylabel(
            "Seeding Z Target  [µm — flow cell surface height used for focusing]",
            fontsize=8)
        ax.set_title(f"{surface} Surface", fontsize=9, fontweight="bold")
        ncol_ = max(1, len(swaths) // 4 + 1)
        ax.legend(fontsize=7, ncol=ncol_, loc="best")
        clean_ax(ax)

    fig.suptitle(
        "Seeding Z Focus Position per Cycle  —  "
        "gradual drift is normal; abrupt changes or diverging swaths are not",
        fontsize=10, fontweight="bold")
    fig.tight_layout()
    return fig_to_rl(fig, dpi=150)


# ── Autofocus: spot intensity & holding stability ─────────────────────────────

def plot_af_optical(df: pd.DataFrame) -> "RLImage | None":
    int_avail  = [c for c in ("left_avg_spot_inten","right_avg_spot_inten")
                  if c in df.columns]
    stab_avail = [c for c in ("HoldingStabilityBlue_nm","HoldingStabilityGreen_nm")
                  if c in df.columns]
    n = (1 if int_avail else 0) + (1 if stab_avail else 0)
    if n == 0:
        return None

    fig, axes = plt.subplots(1, n, figsize=(6.5 * n, 4.0))
    if n == 1:
        axes = [axes]
    idx = 0

    if int_avail:
        ax = axes[idx]; idx += 1
        labels = {"left_avg_spot_inten":  ("Left spot",  C["blue"]),
                  "right_avg_spot_inten": ("Right spot", C["orange"])}
        for col, (lbl, clr) in labels.items():
            if col not in df.columns:
                continue
            grp = df.groupby("Cycle")[col].mean()
            ax.plot(grp.index, grp.values, color=clr, linewidth=1.3, label=lbl)
        ax.set_xlabel("Cycle  [sequencing cycle number]", fontsize=8)
        ax.set_ylabel(
            "Mean Spot Intensity  [ADU — autofocus laser reflected signal strength]",
            fontsize=8)
        ax.set_title("Autofocus Spot Intensity", fontsize=9, fontweight="bold")
        ax.legend(fontsize=8)
        clean_ax(ax)

    if stab_avail:
        ax = axes[idx]; idx += 1
        stab_labels = {
            "HoldingStabilityBlue_nm":  ("Blue laser",  C["blue"]),
            "HoldingStabilityGreen_nm": ("Green laser", C["green"])}
        for col, (lbl, clr) in stab_labels.items():
            if col not in df.columns:
                continue
            for surface in sorted(df["Surface"].dropna().unique()):
                sub = df[df["Surface"] == surface]
                grp = sub.groupby("Cycle")[col].median()
                ls_ = "-" if surface == "Top" else "--"
                ax.plot(grp.index, grp.values, color=clr, linestyle=ls_,
                        linewidth=1.2, alpha=0.9, label=f"{lbl} ({surface})")
        ax.axhline(30, color=C["red"], linestyle="--", linewidth=0.9,
                   label="30 nm guideline")
        ax.set_xlabel("Cycle  [sequencing cycle number]", fontsize=8)
        ax.set_ylabel(
            "Holding Stability  [nm — stage vibration amplitude during exposure]",
            fontsize=8)
        ax.set_title("Stage Holding Stability", fontsize=9, fontweight="bold")
        ax.legend(fontsize=7, ncol=2)
        clean_ax(ax)

    fig.suptitle(
        "Optical & Mechanical Health  —  declining intensity or rising stability "
        "may signal consumable depletion or hardware issues",
        fontsize=10, fontweight="bold")
    fig.tight_layout()
    return fig_to_rl(fig, dpi=150)


# ── Autofocus: gain ───────────────────────────────────────────────────────────

def plot_af_gain(df: pd.DataFrame) -> "RLImage | None":
    if "af_gain" not in df.columns:
        return None
    surfaces = sorted(df["Surface"].dropna().unique())
    fig, ax  = plt.subplots(figsize=(10, 3.5))
    for surface in surfaces:
        sub = df[df["Surface"] == surface]
        grp = sub.groupby("Cycle")["af_gain"].agg(["mean","min","max"])
        col = SURFACE_COLORS.get(surface, C["blue"])
        ax.plot(grp.index, grp["mean"], color=col, linewidth=1.3, label=surface)
        ax.fill_between(grp.index, grp["min"], grp["max"],
                        color=col, alpha=0.15)
    ax.set_xlabel("Cycle  [sequencing cycle number]", fontsize=8)
    ax.set_ylabel(
        "AF Gain  [instrument-internal focus model sensitivity — step changes = recalibration]",
        fontsize=8)
    ax.set_title("Autofocus Gain per Cycle  —  shaded band = min/max across tiles",
                 fontsize=9, fontweight="bold")
    ax.legend(fontsize=8)
    clean_ax(ax)
    fig.tight_layout()
    return fig_to_rl(fig, dpi=150)


# ── Autofocus: StallGuard ─────────────────────────────────────────────────────

def plot_af_stallguard(df: pd.DataFrame) -> "RLImage | None":
    avail = [c for c in ("MinimumStallguardY","MinimumStallguardX")
             if c in df.columns]
    if not avail:
        return None
    fig, axes = plt.subplots(1, len(avail), figsize=(6.5 * len(avail), 3.5))
    if len(avail) == 1:
        axes = [axes]
    labels = {"MinimumStallguardY": ("StallGuard Y", C["orange"]),
              "MinimumStallguardX": ("StallGuard X", C["blue"])}
    for ax, col in zip(axes, avail):
        lbl, clr = labels[col]
        real = df[df[col].notna()]
        if real.empty:
            ax.text(0.5, 0.5, "No events recorded",
                    transform=ax.transAxes, ha="center", color=C["grey"])
            continue
        grp = real.groupby("Cycle")[col].min()
        ax.plot(grp.index, grp.values, color=clr, linewidth=1.2)
        ax.set_xlabel("Cycle  [sequencing cycle number]", fontsize=8)
        ax.set_ylabel(
            f"{lbl}  [motor load value — lower = higher friction / stall risk]",
            fontsize=8)
        ax.set_title(f"Min {lbl}  —  persistently low values may precede motor stall",
                     fontsize=9, fontweight="bold")
        clean_ax(ax)
    fig.tight_layout()
    return fig_to_rl(fig, dpi=150)


# ── AutoCenter ────────────────────────────────────────────────────────────────

def plot_autocenter(df: pd.DataFrame) -> "RLImage | None":
    spec = 0.1
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))
    for ax, col, ylabel in zip(
        axes,
        ["Difference X in mm", "Difference Y in mm"],
        ["Alignment Diff. X  [mm — deviation from expected lane X centre]",
         "Alignment Diff. Y  [mm — deviation from expected lane Y centre]"]
    ):
        if col not in df.columns:
            ax.text(0.5, 0.5, "Column not available",
                    transform=ax.transAxes, ha="center"); continue
        vals       = df[col].values.astype(float)
        facecolors = ["#2ca02c" if abs(v) <= spec else "#d62728" for v in vals]
        ax.bar(np.arange(len(vals)), vals, color=facecolors, alpha=0.85,
               edgecolor="white")
        ax.axhline( spec, color=C["red"], linestyle="--", linewidth=1.0,
                    label=f"±{spec} mm spec")
        ax.axhline(-spec, color=C["red"], linestyle="--", linewidth=1.0)
        ax.axhline(0, color=C["grey"], linewidth=0.5)
        labels = df.get("Surface", pd.Series([""] * len(df))).values
        ax.set_xticks(np.arange(len(vals)))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.set_xlabel("Measurement  [flow cell surface]", fontsize=8)
        ax.set_title(col.replace(" in mm",""), fontsize=9, fontweight="bold")
        ax.legend(fontsize=8)
        clean_ax(ax)
    fig.suptitle(
        "AutoCenter  —  flow cell lane registration accuracy "
        "(green = pass, red = fail)",
        fontsize=10, fontweight="bold")
    fig.tight_layout()
    return fig_to_rl(fig, dpi=150)


# ── AutoTilt ──────────────────────────────────────────────────────────────────

def plot_autotilt(autotilt: dict) -> "RLImage | None":
    has_surf   = "surfaces_pre" in autotilt
    has_motors = "motors" in autotilt
    if not has_surf and not has_motors:
        return None

    ncols = sum([has_surf, has_motors])
    fig   = plt.figure(figsize=(6.5 * ncols, 4.8))
    gs    = gridspec.GridSpec(1, ncols, figure=fig, wspace=0.38)
    ci    = 0
    if has_surf:
        ax = fig.add_subplot(gs[0, ci]); ci += 1
        # Only plot whichever surface datasets actually exist
        datasets = []
        if autotilt.get("surfaces_pre") is not None and not autotilt["surfaces_pre"].empty:
            datasets.append(("Before correction", autotilt["surfaces_pre"], "o", 220))
        if autotilt.get("surfaces_post") is not None and not autotilt["surfaces_post"].empty:
            datasets.append(("After correction", autotilt["surfaces_post"], "^", 220))

        sc = None
        all_z = pd.concat([d[1]["Z_corrected_um"] for d in datasets]).values.astype(float)
        vmin_, vmax_ = all_z.min() - 5, all_z.max() + 5

        for label, df_s, marker, sz in datasets:
            zvals = df_s["Z_corrected_um"].values.astype(float)
            sc = ax.scatter(df_s["X_mm"], df_s["Y_mm"],
                            c=zvals, cmap="RdYlGn", s=sz, marker=marker,
                            edgecolors="grey", linewidth=0.5,
                            vmin=vmin_, vmax=vmax_,
                            label=label, zorder=3)
            for _, row in df_s.iterrows():
                ax.annotate(f"{row['Z_corrected_um']:.1f}",
                            (row["X_mm"], row["Y_mm"]),
                            fontsize=7, ha="center", va="bottom",
                            xytext=(0, 9), textcoords="offset points")

        if sc is not None:
            cbar = fig.colorbar(sc, ax=ax, shrink=0.8)
            cbar.set_label("Z Corrected  [µm — absolute surface height]", fontsize=8)

        has_both = len(datasets) == 2
        title_suffix = "  —  ○ before, △ after correction" if has_both else "  —  before tilt correction"
        ax.set_xlabel("X  [mm — flow cell X coordinate]", fontsize=8)
        ax.set_ylabel("Y  [mm — flow cell Y coordinate]", fontsize=8)
        ax.set_title("Surface Z Map" + title_suffix, fontsize=9, fontweight="bold")
        if has_both:
            ax.legend(fontsize=7)
        clean_ax(ax)


    if has_motors:
        ax  = fig.add_subplot(gs[0, ci]); ci += 1
        mot = autotilt["motors"]
        x   = np.arange(len(mot))
        ax.bar(x - 0.2, mot["Start_um"], 0.35, label="Before",
               color=C["blue"], alpha=0.8)
        ax.bar(x + 0.2, mot["End_um"],   0.35, label="After",
               color=C["green"], alpha=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([f"Motor {int(m)}" for m in mot["Motor"]], fontsize=8)
        ax.set_ylabel(
            "Tilt Motor Position  [µm — actuator displacement from home]",
            fontsize=8)
        ax.set_xlabel("Actuator  [levelling motor ID]", fontsize=8)
        ax.set_title("Tilt Motor Positions Before/After Correction",
                     fontsize=9, fontweight="bold")
        ax.legend(fontsize=8)
        clean_ax(ax)

    fig.suptitle("AutoTilt — Flow Cell Levelling", fontsize=10, fontweight="bold")
    fig.tight_layout()
    return fig_to_rl(fig, dpi=150)


# ═══════════════════════════════════════════════════════════════════════════════
# ReportLab helpers
# ═══════════════════════════════════════════════════════════════════════════════

def get_styles():
    base = getSampleStyleSheet()
    return {
        "title":    ParagraphStyle("T",  parent=base["Title"],   fontSize=22,
                                   textColor=RL_DARK, alignment=TA_CENTER,
                                   spaceAfter=4),
        "sub":      ParagraphStyle("S",  parent=base["Normal"],  fontSize=12,
                                   textColor=RL_ACCENT, alignment=TA_CENTER,
                                   spaceAfter=8),
        "section":  ParagraphStyle("H",  parent=base["Heading1"], fontSize=12,
                                   textColor=RL_ACCENT, spaceBefore=10,
                                   spaceAfter=4),
        "body":     ParagraphStyle("B",  parent=base["Normal"],  fontSize=9,
                                   leading=13),
        "warn":     ParagraphStyle("W",  parent=base["Normal"],  fontSize=9,
                                   textColor=RL_FAIL),
        "metric_v": ParagraphStyle("MV", parent=base["Normal"],  fontSize=16,
                                   fontName="Helvetica-Bold", textColor=RL_ACCENT,
                                   alignment=TA_CENTER),
        "metric_l": ParagraphStyle("ML", parent=base["Normal"],  fontSize=8,
                                   textColor=colors.grey, alignment=TA_CENTER),
    }


def metric_box(value: str, label: str, st: dict) -> Table:
    t = Table([[Paragraph(value, st["metric_v"])],
               [Paragraph(label, st["metric_l"])]],
              colWidths=[4.5*cm], rowHeights=[1.0*cm, 0.6*cm])
    t.setStyle(TableStyle([
        ("BOX",           (0,0),(-1,-1), 0.5, RL_GREY),
        ("BACKGROUND",    (0,0),(-1,-1), RL_LIGHT_BG),
        ("ALIGN",         (0,0),(-1,-1), "CENTER"),
        ("VALIGN",        (0,0),(-1,-1), "MIDDLE"),
        ("TOPPADDING",    (0,0),(-1,-1), 4),
        ("BOTTOMPADDING", (0,0),(-1,-1), 4),
    ]))
    return t


def df_table(df: pd.DataFrame, col_widths=None, stripe=True) -> Table:
    data = [list(df.columns)] + df.astype(str).values.tolist()
    if col_widths is None:
        avail = PAGE_W - 2 * MARGIN
        col_widths = [avail / len(df.columns)] * len(df.columns)
    t = Table(data, colWidths=col_widths, repeatRows=1)
    style = [
        ("BACKGROUND",     (0,0),(-1,0),  RL_ACCENT),
        ("TEXTCOLOR",      (0,0),(-1,0),  colors.white),
        ("FONTNAME",       (0,0),(-1,0),  "Helvetica-Bold"),
        ("FONTSIZE",       (0,0),(-1,-1), 8),
        ("ALIGN",          (0,0),(-1,-1), "CENTER"),
        ("VALIGN",         (0,0),(-1,-1), "MIDDLE"),
        ("GRID",           (0,0),(-1,-1), 0.3, RL_GREY),
        ("TOPPADDING",     (0,0),(-1,-1), 3),
        ("BOTTOMPADDING",  (0,0),(-1,-1), 3),
    ]
    if stripe:
        style.append(("ROWBACKGROUNDS", (0,1),(-1,-1),
                      [colors.white, colors.HexColor("#eef3f8")]))
    t.setStyle(TableStyle(style))
    return t


def info_tbl(rows: list, key_width=4.5*cm) -> Table:
    t = Table(rows, colWidths=[key_width, PAGE_W - 2*MARGIN - key_width])
    t.setStyle(TableStyle([
        ("FONTNAME",       (0,0),(0,-1), "Helvetica-Bold"),
        ("FONTSIZE",       (0,0),(-1,-1), 9),
        ("TEXTCOLOR",      (0,0),(0,-1), RL_ACCENT),
        ("GRID",           (0,0),(-1,-1), 0.3, RL_GREY),
        ("ROWBACKGROUNDS", (0,0),(-1,-1),
         [colors.white, colors.HexColor("#eef3f8")]),
        ("TOPPADDING",     (0,0),(-1,-1), 4),
        ("BOTTOMPADDING",  (0,0),(-1,-1), 4),
        ("VALIGN",         (0,0),(-1,-1), "TOP"),
    ]))
    return t


def _header_footer(cv, doc, subtitle: str, run_id: str):
    cv.saveState()
    w, h = A4
    cv.setFillColor(RL_ACCENT)
    cv.rect(0, h-1.8*cm, w, 1.8*cm, fill=True, stroke=False)
    cv.setFillColor(colors.white)
    cv.setFont("Helvetica-Bold", 11)
    cv.drawString(MARGIN, h-1.2*cm, "Illumina Run QC Report")
    cv.setFont("Helvetica", 9)
    cv.drawRightString(w-MARGIN, h-1.2*cm, subtitle)
    cv.setFillColor(RL_GREY)
    cv.rect(0, 0, w, 0.9*cm, fill=True, stroke=False)
    cv.setFillColor(colors.white)
    cv.setFont("Helvetica", 8)
    cv.drawString(MARGIN, 0.3*cm, f"Run: {run_id}")
    cv.drawCentredString(w/2, 0.3*cm,
        f"Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    cv.drawRightString(w-MARGIN, 0.3*cm, f"Page {cv.getPageNumber()}")
    cv.restoreState()


def sec(title, st):
    return [Paragraph(title, st["section"]),
            HRFlowable(width="100%", thickness=0.5, color=RL_GREY),
            Spacer(1, 0.3*cm)]


def sp(n=1):
    return Spacer(1, n * 0.4*cm)


# ═══════════════════════════════════════════════════════════════════════════════
# Report builder
# ═══════════════════════════════════════════════════════════════════════════════

def build_report(runfolder: Path, output_pdf: Path):
    run_id = runfolder.name

    primary    = load_primary_metrics(runfolder)
    xy_stage   = load_xy_stage(runfolder)
    autofocus  = load_autofocus(runfolder)
    autocenter = load_autocenter(runfolder)
    autotilt_  = load_autotilt(runfolder)
    recipe     = load_recipe(runfolder)
    run_info   = load_run_info(runfolder)

    instrument, serial = detect_instrument(runfolder, recipe, run_info)
    st = get_styles()

    def on_page(cv, doc):
        _header_footer(cv, doc, instrument, run_id)

    doc = SimpleDocTemplate(
        str(output_pdf), pagesize=A4,
        leftMargin=MARGIN, rightMargin=MARGIN,
        topMargin=2.2*cm, bottomMargin=1.4*cm,
        title=f"Illumina QC — {run_id}")

    story = []

    # ── PAGE 1: OVERVIEW ──────────────────────────────────────────────────────
    story += [sp(2),
              Paragraph("Run QC Report", st["title"]),
              Paragraph(f"Instrument: {instrument}  |  Serial: {serial}" if instrument != "N/A" else f"Serial: {serial}", st["sub"]),
              sp(),
              HRFlowable(width="100%", thickness=1, color=RL_ACCENT),
              sp()]

    run_date = "N/A"
    m = re.match(r"(\d{2,4}-?\d{2}-?\d{2})", run_id)
    if m:
        raw = m.group(1).replace("-","")
        try:
            run_date = (datetime.strptime(raw, "%y%m%d")
                        if len(raw) == 6
                        else datetime.strptime(raw, "%Y%m%d")
                        ).strftime("%Y-%m-%d")
        except Exception:
            run_date = m.group(1)

    flowcell = run_info.get("flowcell","N/A") if run_info else "N/A"
    parts_id = run_id.split("_")
    run_num  = parts_id[2] if len(parts_id) > 2 else "N/A"

    # Parse flowcell ID from runfolder name (4th field) if not from RunInfo.xml
    fc_from_folder = parts_id[3] if len(parts_id) >= 4 else ""
    if flowcell == "N/A" and fc_from_folder:
        flowcell = fc_from_folder

    # Flowcell type hint from ID prefix (NextSeq 2000 P-series)
    FLOWCELL_TYPES = {
        "AAG": "P1 (~100M reads)",
        "AAH": "P2 (~400M reads)",
        "AAF": "P3 (~1.2B reads)",
        "AAE": "P3 (~1.2B reads)",
        "20B": "NovaSeq X 10B",
        "25B": "NovaSeq X 25B",
    }
    fc_type = ""
    if flowcell and flowcell != "N/A":
        fc_prefix = flowcell[:3].upper()
        fc_type   = FLOWCELL_TYPES.get(fc_prefix, "")
    fc_display = f"{flowcell}  ({fc_type})" if fc_type else flowcell

    # Detection sources note for overview page
    sources_used = []
    if run_info and run_info.get("instrument_type"):
        sources_used.append("RunInfo.xml (instrument model)")
    if run_info and run_info.get("serial"):
        sources_used.append("RunInfo.xml (serial)")
    else:
        sources_used.append("runfolder name (serial)")
    if run_info and run_info.get("flowcell"):
        sources_used.append("RunInfo.xml (flow cell)")
    elif fc_from_folder:
        sources_used.append("runfolder name (flow cell)")
    detection_note = "Sources: " + "; ".join(sources_used)

    info_rows = [
        ["Run Folder",  run_id],
        ["Run Date",    run_date],
        ["Instrument",  instrument if instrument != "N/A" else "N/A  (RunInfo.xml not found)"],
        ["Serial No.",  serial],
        ["Run Number",  run_num],
        ["Flow Cell ID",    fc_display],
        ["Detection Note",  detection_note],
        ["Recipe",      recipe["name"] if recipe else "N/A"],
        ["Description", recipe["description"] if recipe else "N/A"],
    ]
    if run_info and run_info.get("reads"):
        read_str = "  |  ".join(
            f"{'Index' if r['indexed'] else 'Read'}: {r['cycles']} cy"
            for r in run_info["reads"])
        info_rows.append(["Read Structure", read_str])
    if recipe and recipe.get("focus_model_freq"):
        info_rows.append(["Focus Model Update",
                          f"Every {recipe['focus_model_freq']} cycles"])

    story.append(info_tbl(info_rows))
    story.append(sp(1.5))

    # Key metric boxes
    if primary is not None and not primary.empty:
        story += sec("Key Quality Metrics", st)
        boxes = []
        for mname, label in [
            ("≥ Q30",                   "% Bases ≥ Q30"),
            ("Total Yield",             "Yield (Gbp)"),
            ("Total Reads PF",          "Reads PF (M)"),
            ("% Loading Concentration", "% Loading Conc."),
        ]:
            row_ = primary[primary["Metric"] == mname]
            if row_.empty:
                continue
            val_  = row_.iloc[0]["Value"]
            unit_ = str(row_.iloc[0]["Unit"]).strip()
            disp  = f"{val_:g} {unit_}".strip() if not pd.isna(val_) else "N/A"
            boxes.append(metric_box(disp, label, st))
        if boxes:
            pad  = [Spacer(4.8*cm,1)] * (4 - len(boxes))
            row_t = Table([boxes + pad], colWidths=[4.8*cm]*4,
                          rowHeights=[1.8*cm])
            row_t.setStyle(TableStyle([("ALIGN",(0,0),(-1,-1),"CENTER")]))
            story.append(row_t)
            story.append(sp())

    # AutoTilt result badge
    if autotilt_ and "tilt_result" in autotilt_:
        tr   = autotilt_["tilt_result"]
        bcol = RL_PASS if tr["result"] == "PASS" else RL_FAIL
        badge = Table(
            [["AutoTilt",
              Paragraph(tr["result"],
                        ParagraphStyle("b", parent=st["body"],
                                       textColor=bcol,
                                       fontName="Helvetica-Bold",
                                       alignment=TA_CENTER)),
              f"Residual: {tr['residual_um']:.4f} µm  "
              f"(spec ≤ {tr['spec_um']:.1f} µm)"]],
            colWidths=[3.5*cm, 2.5*cm, PAGE_W-2*MARGIN-6*cm])
        badge.setStyle(TableStyle([
            ("BOX",           (0,0),(-1,-1), 0.5, RL_GREY),
            ("BACKGROUND",    (0,0),(-1,-1), RL_LIGHT_BG),
            ("FONTNAME",      (0,0),(0,-1),  "Helvetica-Bold"),
            ("FONTSIZE",      (0,0),(-1,-1), 9),
            ("TEXTCOLOR",     (0,0),(0,-1),  RL_ACCENT),
            ("ALIGN",         (0,0),(-1,-1), "CENTER"),
            ("VALIGN",        (0,0),(-1,-1), "MIDDLE"),
            ("TOPPADDING",    (0,0),(-1,-1), 6),
            ("BOTTOMPADDING", (0,0),(-1,-1), 6),
        ]))
        story.append(badge)

    story.append(PageBreak())

    # ── PAGE 2: PRIMARY METRICS ───────────────────────────────────────────────
    story += sec("Primary Analysis Metrics", st)
    story.append(Paragraph(
        "Overall run output metrics. Each metric is plotted independently "
        "on its own scale so that percentage values and raw counts are never "
        "mixed on a single axis. Q30 ≥ 80% is the standard minimum threshold.",
        st["body"]))
    story.append(sp())

    if primary is not None:
        embed(story, plot_primary_metrics(primary), aspect=0.26)
        cw3 = [(PAGE_W-2*MARGIN)/3]*3
        story.append(df_table(primary.round(2), col_widths=cw3))
    else:
        story.append(Paragraph("⚠ PrimaryAnalysisMetrics.csv not found.", st["warn"]))

    story.append(PageBreak())

    # ── PAGE 3: XY STAGE ─────────────────────────────────────────────────────
    story += sec("XY Stage Position Errors", st)
    story.append(Paragraph(
        "Average per-cycle deviation of the XY imaging stage from its "
        "programmed target position. Errors beyond ±2 µm (dashed) can cause "
        "image misregistration and degrade base-call accuracy. "
        "The right panels show the error distribution across all cycles.",
        st["body"]))
    story.append(sp())

    if xy_stage is not None:
        embed(story, plot_xy_stage(xy_stage), aspect=0.50)
        stats = (xy_stage[["AveragePositionErrorXum","AveragePositionErrorYum"]]
                 .describe().round(4).reset_index())
        stats.columns = ["Statistic","Avg Error X (µm)","Avg Error Y (µm)"]
        story.append(Paragraph("Summary Statistics", st["body"]))
        story.append(sp(0.3))
        story.append(df_table(stats, col_widths=[(PAGE_W-2*MARGIN)/3]*3))
    else:
        story.append(Paragraph("⚠ XYStagePositionErrorsReport.csv not found.", st["warn"]))

    story.append(PageBreak())

    # ── PAGE 4: AF SCATTER ────────────────────────────────────────────────────
    story += sec("Autofocus — Per-Tile Position Errors (Scatter)", st)
    story.append(Paragraph(
        "Every dot is one tile measurement; colours identify swaths. "
        "The dashed circle marks the 2 µm radius specification. "
        "Asymmetric clusters, swath-specific outliers, or a large spread "
        "along one axis indicate lane or optical alignment issues.",
        st["body"]))
    story.append(sp())

    if autofocus is not None:
        embed(story, plot_af_scatter(autofocus), aspect=0.65)
    else:
        story.append(Paragraph("⚠ AutofocusReport.csv not found.", st["warn"]))

    story.append(PageBreak())

    # ── PAGE 5: AF ERRORS OVER CYCLES ────────────────────────────────────────
    story += sec("Autofocus — Position Errors over Cycles", st)
    story.append(Paragraph(
        "Mean ± std per swath per cycle. Systematic drift or sudden "
        "step-changes across many cycles indicate hardware events "
        "(thermal expansion, reagent delivery problems, focus model updates).",
        st["body"]))
    story.append(sp())

    if autofocus is not None:
        embed(story, plot_af_errors_cycle(autofocus), aspect=0.58)
    else:
        story.append(Paragraph("⚠ AutofocusReport.csv not found.", st["warn"]))

    story.append(PageBreak())

    # ── PAGE 6: SEEDING Z ─────────────────────────────────────────────────────
    story += sec("Autofocus — Seeding Z Focus Position", st)
    story.append(Paragraph(
        "The Z position targeted during the seeding phase per cycle. "
        "A gradual downward trend is normal and reflects reagent evaporation "
        "compressing the flow cell. Abrupt jumps, divergence between swaths, "
        "or differences between Top and Bottom surfaces are abnormal.",
        st["body"]))
    story.append(sp())

    if autofocus is not None:
        img_z = plot_af_seeding_z(autofocus)
        embed(story, img_z, aspect=0.40) if img_z else story.append(
            Paragraph("⚠ Seeding Z column not available.", st["warn"]))
    else:
        story.append(Paragraph("⚠ AutofocusReport.csv not found.", st["warn"]))

    story.append(PageBreak())

    # ── PAGE 7: OPTICAL HEALTH + GAIN + STALLGUARD ───────────────────────────
    story += sec("Autofocus — Optical Health", st)
    story.append(Paragraph(
        "Left/right spot intensities reflect autofocus laser signal strength. "
        "Stage holding stability (nm) is the vibration amplitude during imaging; "
        "values above 30 nm may degrade image quality.",
        st["body"]))
    story.append(sp())

    if autofocus is not None:
        embed(story, plot_af_optical(autofocus), aspect=0.42)
        story += sec("Autofocus — Gain", st)
        story.append(Paragraph(
            "Autofocus gain reflects the sensitivity of the focus model. "
            "Discrete jumps occur when the instrument recalibrates (every "
            f"{recipe.get('focus_model_freq','N/A') if recipe else 'N/A'} cycles). "
            "Continuous drift is unusual and warrants investigation.",
            st["body"]))
        story.append(sp())
        embed(story, plot_af_gain(autofocus), aspect=0.36)
        story += sec("Motor StallGuard", st)
        story.append(Paragraph(
            "Motor StallGuard measures the Y-axis motor load. "
            "Values near 99999 mean no load was detected (good). "
            "Persistent low values indicate friction or blockage risk.",
            st["body"]))
        story.append(sp())
        embed(story, plot_af_stallguard(autofocus), aspect=0.36)
    else:
        story.append(Paragraph("⚠ AutofocusReport.csv not found.", st["warn"]))

    story.append(PageBreak())

    # ── PAGE 8: AUTOCENTER ────────────────────────────────────────────────────
    story += sec("AutoCenter Alignment", st)
    story.append(Paragraph(
        "Measured deviations from expected flow-cell lane positions. "
        "Green bars are within ±0.1 mm spec; red bars fail and may cause "
        "edge-tile image quality or de-multiplexing problems.",
        st["body"]))
    story.append(sp())

    if autocenter is not None:
        embed(story, plot_autocenter(autocenter), aspect=0.40)
        disp = autocenter[[c for c in autocenter.columns if c != "Timestamp"]].round(4)
        cw_ac = [(PAGE_W-2*MARGIN)/len(disp.columns)]*len(disp.columns)
        story.append(Paragraph("Measurement Table", st["body"]))
        story.append(sp(0.3))
        story.append(df_table(disp, col_widths=cw_ac))
    else:
        story.append(Paragraph("⚠ AutoCenterData CSV not found.", st["warn"]))

    story.append(PageBreak())

    # ── PAGE 9: AUTOTILT ──────────────────────────────────────────────────────
    story += sec("AutoTilt — Flow Cell Levelling", st)
    story.append(Paragraph(
        "Three motorised actuators correct the physical tilt of the flow cell. "
        "The Z map shows surface heights before (circle) and after (triangle) "
        "correction. Residual tilt ≤ 5 µm is required to pass.",
        st["body"]))
    story.append(sp())

    if autotilt_ is not None:
        if "tilt_result" in autotilt_:
            tr = autotilt_["tilt_result"]
            fc = "green" if tr["result"] == "PASS" else "red"
            story.append(Paragraph(
                f"Residual: <b>{tr['residual_um']:.4f} µm</b>  "
                f"(spec ≤ {tr['spec_um']:.1f} µm) — "
                f"<font color='{fc}'><b>{tr['result']}</b></font>",
                st["body"]))
            story.append(sp())
        if "theta" in autotilt_:
            th = autotilt_["theta"]
            story.append(Paragraph(
                f"Tilt angles — ThetaX: {th['x_mrad']} mrad  |  "
                f"ThetaY: {th['y_mrad']} mrad",
                st["body"]))
            story.append(sp())
        embed(story, plot_autotilt(autotilt_), aspect=0.45)
        if "surfaces_pre" in autotilt_:
            story.append(Paragraph("Surface Z Positions (before correction)", st["body"]))
            story.append(sp(0.3))
            sf = autotilt_["surfaces_pre"].round(3)
            story.append(df_table(sf,
                col_widths=[(PAGE_W-2*MARGIN)/len(sf.columns)]*len(sf.columns)))
            story.append(sp())
        if "motors" in autotilt_:
            story.append(Paragraph("Tilt Motor Corrections", st["body"]))
            story.append(sp(0.3))
            mot = autotilt_["motors"].round(3)
            story.append(df_table(mot,
                col_widths=[(PAGE_W-2*MARGIN)/len(mot.columns)]*len(mot.columns)))
    else:
        story.append(Paragraph("⚠ AutotiltData CSV not found.", st["warn"]))

    story.append(PageBreak())

    # ── PAGE 10: RECIPE ───────────────────────────────────────────────────────
    story += sec("Effective Recipe Details", st)
    story.append(Paragraph(
        "The effective recipe is the resolved chemistry protocol executed by "
        "the instrument. Key parameters, thermal settings, primer assignments, "
        "and reagent map versions are listed below.",
        st["body"]))
    story.append(sp())

    if recipe is not None:
        basic = [
            ["Recipe Name",  recipe.get("name","N/A")],
            ["Version",      recipe.get("version","N/A")],
            ["Description",  recipe.get("description","N/A")],
        ]
        if recipe.get("focus_model_freq"):
            basic.append(["Focus Model Update",
                          f"Every {recipe['focus_model_freq']} cycles"])
        if recipe.get("default_primers"):
            for k, v in recipe["default_primers"].items():
                lbl = k.replace("Default","").replace("Primer","").strip()
                basic.append([f"Primer ({lbl})", v])
        if recipe.get("primers"):
            for k, v in recipe["primers"].items():
                basic.append([f"Primer Slot ({k})", v])
        story.append(info_tbl(basic))
        story.append(sp())

        if recipe.get("reagents"):
            story.append(Paragraph(
                "Reagent Maps  —  chemistry module versions used in this run",
                st["body"]))
            story.append(sp(0.3))
            rdf = pd.DataFrame(recipe["reagents"])
            story.append(df_table(rdf, col_widths=[(PAGE_W-2*MARGIN)/len(rdf.columns)]*len(rdf.columns)))
            story.append(sp())

        if recipe.get("temperatures"):
            story.append(Paragraph(
                "Thermal Zone Settings  —  unique enabled temperature setpoints",
                st["body"]))
            story.append(sp(0.3))
            tdf = pd.DataFrame(recipe["temperatures"])
            story.append(df_table(tdf, col_widths=[(PAGE_W-2*MARGIN)/len(tdf.columns)]*len(tdf.columns)))
            story.append(sp())

        if recipe.get("recipe_steps"):
            rs = pd.DataFrame(recipe["recipe_steps"])
            story.append(Paragraph(
                f"Chemistry Steps  —  {len(rs)} unique steps across all protocols "
                "(Start, Safe State and maintenance steps excluded). "
                "Reagent abbreviations are Illumina internal names; "
                "Temps and Waits are key incubation setpoints per step.",
                st["body"]))
            story.append(sp(0.3))
            # Split into two tables to keep columns readable
            # Table 1: Protocol, Step, Version, Temps, Waits
            t1 = rs[["Protocol","Step","Version","Temps","Waits"]].copy()
            avail1 = PAGE_W - 2*MARGIN
            cw1 = [avail1*0.20, avail1*0.22, avail1*0.09, avail1*0.28, avail1*0.21]
            story.append(df_table(t1, col_widths=cw1))
            story.append(sp(0.5))
            # Table 2: Step + Reagents (wider column for reagent list)
            t2 = rs[["Step","Reagents"]].copy()
            avail2 = PAGE_W - 2*MARGIN
            cw2 = [avail2*0.25, avail2*0.75]
            story.append(Paragraph("Reagents used per step:", st["body"]))
            story.append(sp(0.3))
            story.append(df_table(t2, col_widths=cw2))
    else:
        story.append(Paragraph("⚠ EffectiveRecipe.xml not found.", st["warn"]))

    # ── Build ─────────────────────────────────────────────────────────────────
    doc.build(story, onFirstPage=on_page, onLaterPages=on_page)
    print(f"[OK] Report written to: {output_pdf}")


# ═══════════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Generate a PDF QC report from an Illumina runfolder.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument("runfolder",  type=Path, help="Path to Illumina runfolder")
    parser.add_argument("output_pdf", type=Path, help="Output PDF path")
    args = parser.parse_args()

    if not args.runfolder.is_dir():
        print(f"ERROR: '{args.runfolder}' is not a directory.", file=sys.stderr)
        sys.exit(1)

    args.output_pdf.parent.mkdir(parents=True, exist_ok=True)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        build_report(args.runfolder, args.output_pdf)

    for w in caught:
        if issubclass(w.category, UserWarning):
            print(f"WARNING: {w.message}", file=sys.stderr)


if __name__ == "__main__":
    main()