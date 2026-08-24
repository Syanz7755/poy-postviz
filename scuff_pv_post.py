"""
scuff_pv_post.py -- Post-processing tool for SCUFF-EM near-field heat transfer results.

Subcommands:
    collect     Merge SRFlux output files from subtask directories into canonical TSV datasets.
    integrate   Run scuff-integrate on a collected dataset to produce integrated PVMST output.
    plot        Generate cell-map and quiver plots from a collected dataset.
    run-all     Auto-discover all configs and run the full pipeline (collect + plot + SIFlux).
    slice       Extract 1D slices, export data to TSV, and plot quantities along x or y axis.
    slice-plot  Generate slice plots from previously exported data (no re-computation).

Usage:
    # One-line full pipeline:
    python scuff_pv_post.py run-all --task-dir ./results-batches/ord --out-dir ./results-batches [--last]

    # Individual steps:
    python scuff_pv_post.py collect --task-dir ./results-batches/ord \
        --group-name ad1 --out ./results-batches/pv_out_ad1 [--last]

    python scuff_pv_post.py plot --dataset ./results-batches/pv_out_ad1 \
        --planes 0.75,1.0,1.5 --out ./results-batches/pv_out_ad1/figures

    # Slice analysis (export + plot):
    python scuff_pv_post.py slice --dataset ./results-batches/ord_processed/srflux/ad1 \
        --slice-axis y --slice-value 0.01 --temperature 300 \
        --quantity Sx,Sy,absSxy --planes 0.75 --out ./test_slice

    # Slice export only (skip plots):
    python scuff_pv_post.py slice --dataset ... --export-only --out ./test_slice

    # Re-plot from exported slice data:
    python scuff_pv_post.py slice-plot --export-dir ./test_slice --out ./test_slice_figures
"""
from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

__version__ = "0.1.0"

# ---------------------------------------------------------------------------
# Constants / regex
# ---------------------------------------------------------------------------

RUN_HEADER_RE = re.compile(r"^#\s*scuff-neq run on\s+")
SUBTASK_DIR_RE = re.compile(r"^subtask-(.+)-(\d+)$")

FALLBACK_COL_MAP = {
    0: "transform_id",
    1: "omega",
    2: "x",
    3: "y",
    4: "z",
    5: "source_id",
    6: "Sx_flux",
    7: "Sy_flux",
    8: "Sz_flux",
    9: "Mxx_flux",
    10: "Mxy_flux",
    11: "Mxz_flux",
    12: "Myx_flux",
    13: "Myy_flux",
    14: "Myz_flux",
    15: "Mzx_flux",
    16: "Mzy_flux",
    17: "Mzz_flux",
}

# Canonical output column order for merged_srflux.tsv
MERGED_TSV_COLS = [
    "omega", "transform_id", "source_id", "point_id",
    "x", "y", "z",
    "Sx_flux", "Sy_flux", "Sz_flux",
    "subtask_index", "source_file",
    "Mxx_flux", "Mxy_flux", "Mxz_flux",
    "Myx_flux", "Myy_flux", "Myz_flux",
    "Mzx_flux", "Mzy_flux", "Mzz_flux",
]

SRFLUX_OUT_COLS = [
    "transform_id", "omega", "x", "y", "z", "source_id",
    "Sx_flux", "Sy_flux", "Sz_flux",
    "Mxx_flux", "Mxy_flux", "Mxz_flux",
    "Myx_flux", "Myy_flux", "Myz_flux",
    "Mzx_flux", "Mzy_flux", "Mzz_flux",
]

SRFLUX_HEADER = """\
# merged SRFlux for scuff-integrate
# data file columns:
# 1 transform tag
# 2 omega
# 3, 4, 5 x,y,z (coordinates of eval point)
# 6 sourceObject
#  7  8  9 Px,    Py,    Pz
# 10 11 12 Txx,   Txy,   Txz
# 13 14 15 Tyx,   Tyy,   Tyz
# 16 17 18 Tzx,   Tzy,   Tzz
"""

# Bose-Einstein physical constants
HBAR = 1.0545718e-34   # J·s
K_B = 1.380649e-23     # J/K

# Geometry constants (from gen_evpoints.py)
W_GEOM = 1.0
GROUP_DELTA = {"a": 0.6, "b": 0.2, "c": 0.1, "d": 0.6, "e": 0.2, "f": 0.1}

# Name translation from header column names to canonical names
HEADER_TO_CANONICAL = {
    "transform tag": "transform_id",
    "transform_tag": "transform_id",
    "omega": "omega",
    "x": "x",
    "y": "y",
    "z": "z",
    "sourceobject": "source_id",
    "source object": "source_id",
    "px": "Sx_flux",
    "py": "Sy_flux",
    "pz": "Sz_flux",
    "txx": "Mxx_flux",
    "txy": "Mxy_flux",
    "txz": "Mxz_flux",
    "tyx": "Myx_flux",
    "tyy": "Myy_flux",
    "tyz": "Myz_flux",
    "tzx": "Mzx_flux",
    "tzy": "Mzy_flux",
    "tzz": "Mzz_flux",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def is_nan_value(s: str) -> bool:
    return s.strip().lower() in ("nan", "-nan", "+nan")


def make_point_id(x: float, y: float, z: float, tol: float = 1e-9) -> str:
    key = f"{round(x / tol)},{round(y / tol)},{round(z / tol)}"
    h = hashlib.md5(key.encode()).hexdigest()[:12]
    return f"p_{h}"


def write_text_file(path: Path, content: str) -> None:
    path.write_text(content.rstrip() + "\n", encoding="utf-8", newline="\n")


def parse_subtask_name(dirname: str) -> tuple[str, int] | None:
    m = SUBTASK_DIR_RE.match(dirname)
    if m:
        return m.group(1), int(m.group(2))
    return None


# ---------------------------------------------------------------------------
# Subtask discovery
# ---------------------------------------------------------------------------

def discover_subtasks(task_dir: Path, group_name: str) -> list[tuple[Path, int]]:
    """Find subtask directories matching subtask-{group_name}-N, return sorted by N."""
    matches: list[tuple[Path, int]] = []
    for p in task_dir.iterdir():
        if p.is_dir():
            parsed = parse_subtask_name(p.name)
            if parsed is not None and parsed[0] == group_name:
                matches.append((p, parsed[1]))
    matches.sort(key=lambda t: t[1])

    # Check for non-contiguous indices
    if matches:
        indices = [m[1] for m in matches]
        expected = list(range(indices[0], indices[0] + len(indices)))
        if indices != expected:
            print(f"WARNING: subtask indices are non-contiguous: {indices}", file=sys.stderr)

    return matches


# ---------------------------------------------------------------------------
# SRFlux parsing
# ---------------------------------------------------------------------------

def split_runs(lines: list[str]) -> list[list[str]]:
    """Split file lines into runs delimited by '# scuff-neq run on' headers."""
    runs: list[list[str]] = []
    current_run: list[str] = []
    for line in lines:
        if RUN_HEADER_RE.match(line):
            if current_run:
                runs.append(current_run)
            current_run = [line]
        else:
            current_run.append(line)
    if current_run:
        runs.append(current_run)
    return runs


def extract_header_and_data(lines: list[str], last_only: bool) -> tuple[list[str], list[str]]:
    """Return (header_lines, data_lines) from SRFlux content.

    header_lines are the comment lines from the selected run.
    data_lines are the non-comment, non-empty stripped lines.
    """
    if last_only:
        runs = split_runs(lines)
        if not runs:
            return [], []
        target_lines = runs[-1]
    else:
        target_lines = lines

    header_lines: list[str] = []
    data_lines: list[str] = []
    for line in target_lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            header_lines.append(stripped)
        else:
            data_lines.append(stripped)

    return header_lines, data_lines


def parse_column_mapping(header_lines: list[str]) -> tuple[dict[int, str], str, str | None]:
    """Parse column mapping from header lines.

    Returns (column_map, method, warning).
    method is 'header' if parsing succeeded, 'fallback' otherwise.
    """
    # Try to parse header lines like:
    # # 1 transform tag          (single col, multi-word name)
    # # 2 omega                  (single col, single name)
    # # 3, 4, 5 x,y,z (...)     (multi col, comma-separated names)
    # #  7  8  9 Px,    Py,    Pz (multi col, comma-separated names)

    col_map: dict[int, str] = {}
    warning: str | None = None

    for line in header_lines:
        # Match: # <nums_part> <names_part>
        m = re.match(r"^\s*#\s*((?:\d+[\s,]*)+)\s+(.+?)\s*$", line)
        if not m:
            continue

        nums_str = m.group(1).strip()
        names_str = m.group(2).strip()

        # Parse the column numbers
        nums = re.findall(r"\d+", nums_str)

        # Strip parenthetical suffixes from names
        names_str_clean = re.sub(r"\s*\(.*?\)\s*$", "", names_str).strip()

        if len(nums) == 1:
            # Single column: entire name string is one name (may be multi-word)
            names = [names_str_clean]
        elif "," in names_str_clean:
            # Multi-column with comma-separated names
            names = [n.strip() for n in names_str_clean.split(",") if n.strip()]
        else:
            # Multi-column with space-separated names
            names = names_str_clean.split()

        if len(nums) != len(names):
            warning = f"Mismatch: {len(nums)} column numbers but {len(names)} names in: {line}"
            continue

        for num_str, name in zip(nums, names):
            col_idx = int(num_str) - 1  # Convert 1-based to 0-based
            canonical = HEADER_TO_CANONICAL.get(name.lower().strip())
            if canonical:
                col_map[col_idx] = canonical
            else:
                col_map[col_idx] = name.lower().strip()

    if not col_map:
        return FALLBACK_COL_MAP.copy(), "fallback", "No column descriptions found in header"

    # Check we got all 18 columns
    if len(col_map) < 18:
        # Fill in any missing from fallback
        for idx, name in FALLBACK_COL_MAP.items():
            if idx not in col_map:
                col_map[idx] = name
        if warning is None:
            warning = f"Only {len(col_map)} columns parsed from header, filled rest from fallback"

    return col_map, "header", warning


def parse_srflux_rows(
    data_lines: list[str],
    col_map: dict[int, str],
) -> list[dict[str, Any]]:
    """Parse SRFlux data lines into row dicts using column mapping."""
    rows: list[dict[str, Any]] = []
    for line in data_lines:
        parts = line.split()
        if len(parts) < 18:
            continue

        row: dict[str, Any] = {}
        for col_idx, canonical_name in col_map.items():
            if col_idx >= len(parts):
                continue
            val_str = parts[col_idx]
            if is_nan_value(val_str):
                row[canonical_name] = "nan"
            elif canonical_name in ("transform_id", "source_id"):
                row[canonical_name] = int(float(val_str))
            else:
                row[canonical_name] = float(val_str)

        rows.append(row)

    return rows


# ---------------------------------------------------------------------------
# Collect command
# ---------------------------------------------------------------------------

def cmd_collect(args: argparse.Namespace) -> None:
    task_dir = Path(args.task_dir).resolve()
    group_name = args.group_name
    filebase = args.filebase
    out_dir = Path(args.out).resolve()
    coord_tol = args.coord_tol
    dup_policy = args.duplicate_policy
    last_only = args.last

    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Discovery ---
    subtasks = discover_subtasks(task_dir, group_name)
    if not subtasks:
        print(f"ERROR: no subtask-{group_name}-* directories found in {task_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Discovered {len(subtasks)} subtask directories:")
    for path, idx in subtasks:
        print(f"  {path.name} (index {idx})")

    matched_names = [p.name for p, _ in subtasks]

    # --- Parse each subtask ---
    all_rows: list[dict[str, Any]] = []
    raw_headers: dict[str, list[str]] = {}
    col_map: dict[int, str] = {}
    col_method = "fallback"
    col_warning: str | None = None

    for subtask_path, subtask_idx in subtasks:
        srflux_path = subtask_path / f"{filebase}.SRFlux"
        if not srflux_path.exists():
            print(f"WARNING: {srflux_path} not found, skipping", file=sys.stderr)
            continue

        lines = srflux_path.read_text(encoding="utf-8").splitlines()
        header_lines, data_lines = extract_header_and_data(lines, last_only)

        raw_headers[subtask_path.name] = header_lines

        # Parse column mapping from first subtask's header
        if not col_map:
            col_map, col_method, col_warning = parse_column_mapping(header_lines)
            print(f"Column mapping method: {col_method}")
            if col_warning:
                print(f"Column mapping warning: {col_warning}")

        rows = parse_srflux_rows(data_lines, col_map)

        # Annotate rows with subtask metadata
        for row in rows:
            row["subtask_index"] = subtask_idx
            row["source_file"] = f"{subtask_path.name}/{filebase}.SRFlux"

        all_rows.extend(rows)
        print(f"  {subtask_path.name}: {len(rows)} data rows")

    if not all_rows:
        print("ERROR: no data rows parsed from any subtask", file=sys.stderr)
        sys.exit(1)

    # --- Generate point_ids ---
    point_cache: dict[tuple[float, float, float], str] = {}
    for row in all_rows:
        coord_key = (row["x"], row["y"], row["z"])
        if coord_key not in point_cache:
            point_cache[coord_key] = make_point_id(row["x"], row["y"], row["z"], coord_tol)
        row["point_id"] = point_cache[coord_key]

    # --- Deduplication ---
    primary_key_cols = ["omega", "transform_id", "source_id", "point_id"]
    numeric_cols = [
        "Sx_flux", "Sy_flux", "Sz_flux",
        "Mxx_flux", "Mxy_flux", "Mxz_flux",
        "Myx_flux", "Myy_flux", "Myz_flux",
        "Mzx_flux", "Mzy_flux", "Mzz_flux",
    ]

    grouped: dict[tuple, list[dict[str, Any]]] = defaultdict(list)
    for row in all_rows:
        pk = tuple(row[c] for c in primary_key_cols)
        grouped[pk].append(row)

    deduped_rows: list[dict[str, Any]] = []
    n_identical = 0
    n_averaged = 0

    for pk, group in grouped.items():
        if len(group) == 1:
            deduped_rows.append(group[0])
            continue

        # Check if all numeric values are identical
        all_identical = True
        for col in numeric_cols:
            vals = set()
            for r in group:
                v = r.get(col)
                if isinstance(v, (int, float)):
                    vals.add(v)
                else:
                    vals.add(str(v))
            if len(vals) > 1:
                all_identical = False
                break

        if all_identical:
            deduped_rows.append(group[0])
            n_identical += len(group) - 1
        elif dup_policy == "error":
            print(f"ERROR: duplicate primary key {pk} with differing values", file=sys.stderr)
            sys.exit(1)
        else:  # mean
            merged = dict(group[0])
            for col in numeric_cols:
                vals = [r[col] for r in group if isinstance(r.get(col), (int, float))]
                if vals:
                    merged[col] = sum(vals) / len(vals)
            deduped_rows.append(merged)
            n_averaged += 1

    print(f"Deduplication: {n_identical} identical duplicates removed, {n_averaged} averaged")

    # Sort by (omega, transform_id, source_id, z, x, y)
    deduped_rows.sort(key=lambda r: (r["omega"], r["transform_id"], r["source_id"], r["z"], r["x"], r["y"]))

    # --- Build unique points ---
    unique_points: dict[str, dict[str, Any]] = {}
    for row in deduped_rows:
        pid = row["point_id"]
        if pid not in unique_points:
            z = row["z"]
            plane_id = f"z_{z:.6f}"
            unique_points[pid] = {"point_id": pid, "x": row["x"], "y": row["y"], "z": z, "plane_id": plane_id}

    # Sort points by (z, x, y) for stable output
    sorted_points = sorted(unique_points.values(), key=lambda p: (p["z"], p["x"], p["y"]))

    # --- Build unique omegas ---
    omega_seen: dict[float, dict[str, Any]] = {}
    for row in deduped_rows:
        w = row["omega"]
        if w not in omega_seen:
            omega_seen[w] = {"omega": w, "subtask_index": row["subtask_index"], "source_file": row["source_file"]}

    sorted_omegas = sorted(omega_seen.values(), key=lambda d: d["omega"])
    for i, entry in enumerate(sorted_omegas, start=1):
        entry["omega_id"] = f"w{i:06d}"

    # --- Missing point detection ---
    omega_point_counts: dict[float, set[str]] = defaultdict(set)
    for row in deduped_rows:
        omega_point_counts[row["omega"]].add(row["point_id"])

    majority_count = Counter(len(pts) for pts in omega_point_counts.values()).most_common(1)[0][0]
    missing_report: list[str] = []
    for w in sorted(omega_point_counts):
        n = len(omega_point_counts[w])
        if n != majority_count:
            missing_report.append(f"omega={w}: {n} points (expected {majority_count})")

    if missing_report:
        print(f"WARNING: {len(missing_report)} omegas have mismatched point counts")
        for line in missing_report:
            print(f"  {line}")

    # --- Write merged_srflux.tsv ---
    tsv_lines = ["\t".join(MERGED_TSV_COLS)]
    for row in deduped_rows:
        parts: list[str] = []
        for col in MERGED_TSV_COLS:
            val = row[col]
            if col in ("transform_id", "source_id", "subtask_index"):
                parts.append(str(val))
            elif col == "source_file":
                parts.append(str(val))
            elif col == "point_id":
                parts.append(str(val))
            elif isinstance(val, str):  # nan
                parts.append(val)
            else:
                parts.append(f"{val:.10e}")
        tsv_lines.append("\t".join(parts))

    write_text_file(out_dir / "merged_srflux.tsv", "\n".join(tsv_lines))

    # --- Write points.tsv ---
    pts_lines = ["point_id\tx\ty\tz\tplane_id"]
    for p in sorted_points:
        pts_lines.append(f"{p['point_id']}\t{p['x']}\t{p['y']}\t{p['z']}\t{p['plane_id']}")
    write_text_file(out_dir / "points.tsv", "\n".join(pts_lines))

    # --- Write omegas.tsv ---
    omega_lines = ["omega_id\tomega\tsubtask_index\tsource_file"]
    for entry in sorted_omegas:
        omega_lines.append(f"{entry['omega_id']}\t{entry['omega']}\t{entry['subtask_index']}\t{entry['source_file']}")
    write_text_file(out_dir / "omegas.tsv", "\n".join(omega_lines))

    # --- Write merged.SRFlux ---
    srflux_lines = [SRFLUX_HEADER.rstrip()]
    for row in deduped_rows:
        parts = []
        for col in SRFLUX_OUT_COLS:
            val = row[col]
            if isinstance(val, str):  # nan
                parts.append(val)
            else:
                parts.append(f"{val:.10e}")
        srflux_lines.append(" ".join(parts))
    write_text_file(out_dir / "merged.SRFlux", "\n".join(srflux_lines))

    # --- Write raw headers ---
    raw_hdr_dir = out_dir / "raw_headers"
    raw_hdr_dir.mkdir(parents=True, exist_ok=True)
    for subtask_name, hdr_lines in raw_headers.items():
        fname = f"{subtask_name}.SRFlux.header.txt"
        write_text_file(raw_hdr_dir / fname, "\n".join(hdr_lines))

    # --- Write manifest.yaml ---
    manifest = {
        "dataset_version": 1,
        "task_dir": str(task_dir),
        "group_name": group_name,
        "filebase": filebase,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "software": {
            "project_name": "scuff_pv_post",
            "project_version": __version__,
        },
        "discovery": {
            "matched_subtasks": matched_names,
            "n_subtasks": len(subtasks),
        },
        "column_mapping": {
            "method": col_method,
            "column_map_file": None,
            "warning": col_warning,
        },
        "coordinate": {
            "coord_tol": coord_tol,
            "point_id_method": "rounded_coordinate_hash",
        },
        "merge": {
            "primary_key": primary_key_cols,
            "duplicate_policy": dup_policy,
            "n_rows_merged": len(deduped_rows),
            "n_points": len(sorted_points),
            "n_omegas": len(sorted_omegas),
        },
        "validation": {
            "status": "pass" if not missing_report else "warn",
            "reports": {
                "duplicate_report": "duplicate_report.tsv",
                "missing_report": "missing_report.tsv",
            },
        },
    }

    (out_dir / "manifest.yaml").write_text(
        yaml.dump(manifest, default_flow_style=False, sort_keys=False),
        encoding="utf-8",
    )

    print(f"\nCollect complete:")
    print(f"  {len(deduped_rows)} rows merged")
    print(f"  {len(sorted_points)} unique points")
    print(f"  {len(sorted_omegas)} unique omegas")
    print(f"  Output: {out_dir}")


# ---------------------------------------------------------------------------
# PVMST column mapping
# ---------------------------------------------------------------------------

PVMST_TSV_COLS = [
    "transform_id", "source_id", "point_id",
    "x", "y", "z",
    "Sx", "Sy", "Sz",
    "Mxx", "Mxy", "Mxz",
    "Myx", "Myy", "Myz",
    "Mzx", "Mzy", "Mzz",
]

# Header-to-canonical mapping for PVMST files
PVMST_HEADER_TO_CANONICAL = {
    "transform tag": "transform_id",
    "transform_tag": "transform_id",
    "sourceobject": "source_id",
    "source object": "source_id",
    "x": "x",
    "y": "y",
    "z": "z",
    "px": "Sx",
    "py": "Sy",
    "pz": "Sz",
    "txx": "Mxx",
    "txy": "Mxy",
    "txz": "Mxz",
    "tyx": "Myx",
    "tyy": "Myy",
    "tyz": "Myz",
    "tzx": "Mzx",
    "tzy": "Mzy",
    "tzz": "Mzz",
}


# ---------------------------------------------------------------------------
# PVMST parsing
# ---------------------------------------------------------------------------

def parse_pvmst_file(path: Path) -> list[dict[str, Any]]:
    """Parse a .PVMST file produced by scuff-integrate.

    Returns a list of row dicts with canonical column names.
    """
    lines = path.read_text(encoding="utf-8").splitlines()

    header_lines: list[str] = []
    data_lines: list[str] = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            header_lines.append(stripped)
        else:
            data_lines.append(stripped)

    # Parse column mapping from header
    col_map: dict[int, str] = {}
    for line in header_lines:
        m = re.match(r"^\s*#\s*((?:\d+[\s,]*)+)\s+(.+?)\s*$", line)
        if not m:
            continue

        nums_str = m.group(1).strip()
        names_str = m.group(2).strip()

        nums = re.findall(r"\d+", nums_str)

        # Strip parenthetical suffixes
        names_str_clean = re.sub(r"\s*\(.*?\)\s*$", "", names_str).strip()

        if len(nums) == 1:
            names = [names_str_clean]
        elif "," in names_str_clean:
            names = [n.strip() for n in names_str_clean.split(",") if n.strip()]
        else:
            names = names_str_clean.split()

        if len(nums) != len(names):
            continue

        for num_str, name in zip(nums, names):
            col_idx = int(num_str) - 1  # 1-based to 0-based
            canonical = PVMST_HEADER_TO_CANONICAL.get(name.lower().strip())
            if canonical:
                col_map[col_idx] = canonical
            else:
                col_map[col_idx] = name.lower().strip()

    # Fallback if header parsing failed
    if not col_map:
        col_map = {
            0: "transform_id",
            1: "source_id",
            2: "x",
            3: "y",
            4: "z",
            5: "Sx",
            6: "Sy",
            7: "Sz",
            8: "Mxx",
            9: "Mxy",
            10: "Mxz",
            11: "Myx",
            12: "Myy",
            13: "Myz",
            14: "Mzx",
            15: "Mzy",
            16: "Mzz",
        }

    # Parse data rows
    rows: list[dict[str, Any]] = []
    for line in data_lines:
        parts = line.split()
        if len(parts) < 17:
            continue

        row: dict[str, Any] = {}
        for col_idx, canonical_name in col_map.items():
            if col_idx >= len(parts):
                continue
            val_str = parts[col_idx]
            if is_nan_value(val_str):
                row[canonical_name] = "nan"
            elif canonical_name in ("transform_id", "source_id"):
                row[canonical_name] = int(float(val_str))
            else:
                row[canonical_name] = float(val_str)

        rows.append(row)

    return rows


# ---------------------------------------------------------------------------
# Integrate command
# ---------------------------------------------------------------------------

def cmd_integrate(args: argparse.Namespace) -> None:
    dataset_dir = Path(args.dataset).resolve()
    temperature_file = Path(args.temperature_file).resolve()
    backend = args.backend
    scuff_integrate_path = args.scuff_integrate_path

    # Check merged.SRFlux exists
    srflux_path = dataset_dir / "merged.SRFlux"
    if not srflux_path.exists():
        print(f"ERROR: merged.SRFlux not found in {dataset_dir}", file=sys.stderr)
        sys.exit(1)

    # Check temperature file exists
    if not temperature_file.exists():
        print(f"ERROR: temperature file not found: {temperature_file}", file=sys.stderr)
        sys.exit(1)

    # Check that scuff-integrate is available
    if not shutil.which(scuff_integrate_path):
        print(
            f"ERROR: '{scuff_integrate_path}' not found in PATH.\n"
            f"  Install SCUFF-EM or use --scuff-integrate-path to specify the binary location.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Build command
    cmd = [
        scuff_integrate_path,
        "--fluxfile", str(srflux_path),
        "--temperature-file", str(temperature_file),
    ]

    print(f"Running: {' '.join(cmd)}")

    # Run scuff-integrate
    result = subprocess.run(
        cmd,
        cwd=str(dataset_dir),
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        # Save stderr for diagnostics
        stderr_path = dataset_dir / "scuff-integrate.stderr.txt"
        write_text_file(stderr_path, result.stderr)
        print(f"ERROR: scuff-integrate returned non-zero exit code {result.returncode}", file=sys.stderr)
        print(f"  stderr saved to {stderr_path}", file=sys.stderr)
        if result.stderr.strip():
            print(f"  stderr: {result.stderr[:500]}", file=sys.stderr)

        # Write failure report
        report = {
            "status": "failed",
            "backend": backend,
            "command": " ".join(cmd),
            "return_code": result.returncode,
            "merged_srflux": srflux_path.name,
            "temperature_file": temperature_file.name,
            "output_pvmst": None,
            "output_integrated": None,
            "error_stderr_file": str(stderr_path.name),
        }
        (dataset_dir / "integration_report.yaml").write_text(
            yaml.dump(report, default_flow_style=False, sort_keys=False),
            encoding="utf-8",
        )
        sys.exit(1)

    print("scuff-integrate completed successfully")

    # Look for the PVMST output file
    pvmst_path = dataset_dir / "merged.PVMST"
    if not pvmst_path.exists():
        # Try alternate naming patterns
        candidates = list(dataset_dir.glob("*.PVMST"))
        if candidates:
            pvmst_path = candidates[0]
            print(f"  Found PVMST file: {pvmst_path.name}")
        else:
            print("ERROR: no .PVMST file generated by scuff-integrate", file=sys.stderr)
            report = {
                "status": "failed",
                "backend": backend,
                "command": " ".join(cmd),
                "return_code": 0,
                "merged_srflux": srflux_path.name,
                "temperature_file": temperature_file.name,
                "output_pvmst": None,
                "output_integrated": None,
                "error": "No .PVMST file generated",
            }
            (dataset_dir / "integration_report.yaml").write_text(
                yaml.dump(report, default_flow_style=False, sort_keys=False),
                encoding="utf-8",
            )
            sys.exit(1)

    # Parse the PVMST file
    print(f"Parsing {pvmst_path.name}...")
    rows = parse_pvmst_file(pvmst_path)
    print(f"  Parsed {len(rows)} rows")

    # Generate point_ids for rows that don't have them
    for row in rows:
        if "point_id" not in row:
            row["point_id"] = make_point_id(row["x"], row["y"], row["z"])

    # Write pvmst_integrated.tsv
    tsv_lines = ["\t".join(PVMST_TSV_COLS)]
    for row in rows:
        parts: list[str] = []
        for col in PVMST_TSV_COLS:
            val = row.get(col, "nan")
            if col in ("transform_id", "source_id"):
                parts.append(str(val))
            elif col == "point_id":
                parts.append(str(val))
            elif isinstance(val, str):  # nan
                parts.append(val)
            else:
                parts.append(f"{val:.10e}")
        tsv_lines.append("\t".join(parts))

    integrated_path = dataset_dir / "pvmst_integrated.tsv"
    write_text_file(integrated_path, "\n".join(tsv_lines))
    print(f"  Wrote {integrated_path}")

    # Write success report
    report = {
        "status": "success",
        "backend": backend,
        "command": " ".join(cmd),
        "return_code": 0,
        "merged_srflux": srflux_path.name,
        "temperature_file": temperature_file.name,
        "output_pvmst": pvmst_path.name,
        "output_integrated": integrated_path.name,
    }
    report_path = dataset_dir / "integration_report.yaml"
    report_path.write_text(
        yaml.dump(report, default_flow_style=False, sort_keys=False),
        encoding="utf-8",
    )

    print(f"\nIntegration complete:")
    print(f"  PVMST: {pvmst_path.name}")
    print(f"  TSV:   {integrated_path.name}")
    print(f"  Rows:  {len(rows)}")


# ---------------------------------------------------------------------------
# Plot command
# ---------------------------------------------------------------------------

COLOR_BY_OPTIONS = [
    "absSxy", "Sz", "Sx", "Sy", "absS", "log_absSxy", "log_absS",
]

SIGNED_COLOR_BY = {"Sx", "Sy", "Sz"}


def axis_cell_edges(coords: np.ndarray) -> np.ndarray:
    """Return nonuniform cell boundaries from sorted cell-center coordinates."""
    coords = np.asarray(coords, dtype=float)
    if coords.size == 0:
        raise ValueError("cannot build cell edges for an empty coordinate axis")
    if coords.size == 1:
        return np.array([coords[0] - 0.5, coords[0] + 0.5], dtype=float)

    mids = (coords[:-1] + coords[1:]) / 2
    first_width = coords[1] - coords[0]
    last_width = coords[-1] - coords[-2]
    return np.concatenate([[coords[0] - first_width / 2], mids, [coords[-1] + last_width / 2]])


def can_use_rectilinear_cellmap(n_x: int, n_y: int, n_points: int) -> bool:
    """Avoid treating fully scattered points as an enormous rectilinear grid."""
    if n_x < 2 or n_y < 2:
        return False
    grid_cells = n_x * n_y
    return grid_cells <= max(4 * n_points, 250_000)


def compute_color_value(row: dict[str, Any], color_by: str) -> float:
    """Compute the scalar color value for a data row."""
    Sx = row.get("Sx", 0.0) or 0.0
    Sy = row.get("Sy", 0.0) or 0.0
    Sz = row.get("Sz", 0.0) or 0.0

    if color_by == "absSxy":
        return float(np.sqrt(Sx**2 + Sy**2))
    elif color_by == "Sz":
        return float(Sz)
    elif color_by == "Sx":
        return float(Sx)
    elif color_by == "Sy":
        return float(Sy)
    elif color_by == "absS":
        return float(np.sqrt(Sx**2 + Sy**2 + Sz**2))
    elif color_by == "log_absSxy":
        return float(np.log10(np.sqrt(Sx**2 + Sy**2) + 1e-30))
    elif color_by == "log_absS":
        return float(np.log10(np.sqrt(Sx**2 + Sy**2 + Sz**2) + 1e-30))
    else:
        return float(np.sqrt(Sx**2 + Sy**2))


def load_data(dataset_dir: Path) -> "pd.DataFrame":
    """Load integrated or frequency-summed data from a dataset directory.

    Prefers pvmst_integrated.tsv (from integrate stage).
    Falls back to summing merged_srflux.tsv across omegas.
    """
    integrated_path = dataset_dir / "pvmst_integrated.tsv"
    if integrated_path.exists():
        df = pd.read_csv(integrated_path, sep="\t")
        # Ensure numeric columns
        for col in ["Sx", "Sy", "Sz", "x", "y", "z"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        return df

    # Fallback: sum merged_srflux.tsv across omegas
    srflux_path = dataset_dir / "merged_srflux.tsv"
    if not srflux_path.exists():
        print(f"ERROR: neither pvmst_integrated.tsv nor merged_srflux.tsv found in {dataset_dir}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(srflux_path, sep="\t")
    flux_cols = ["Sx_flux", "Sy_flux", "Sz_flux"]
    for col in flux_cols + ["x", "y", "z", "omega"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Group by (transform_id, source_id, point_id, x, y, z) and sum flux across omegas
    group_cols = ["transform_id", "source_id", "point_id", "x", "y", "z"]
    agg_dict = {col: "sum" for col in flux_cols}
    integrated = df.groupby(group_cols, as_index=False).agg(agg_dict)

    # Rename flux columns to match pvmst_integrated.tsv naming
    integrated = integrated.rename(columns={
        "Sx_flux": "Sx",
        "Sy_flux": "Sy",
        "Sz_flux": "Sz",
    })

    return integrated


def _make_fig_for_aspect(x_vals, y_vals, target_height=7.0, max_width=14.0, min_height=4.0):
    """Create a figure whose size matches the data aspect ratio, keeping colorbar height ≈ axes height."""
    import matplotlib
    import matplotlib.pyplot as plt

    x_range = max(x_vals) - min(x_vals) if len(x_vals) > 1 else 1.0
    y_range = max(y_vals) - min(y_vals) if len(y_vals) > 1 else 1.0
    data_aspect = x_range / y_range  # width / height

    # Reserve fraction for colorbar on the right
    cbar_frac = 0.06
    axes_width_frac = 1.0 - cbar_frac - 0.05  # small margin

    # target: axes height ≈ target_height
    axes_height = target_height
    axes_width = axes_height * data_aspect

    # cap width
    if axes_width > max_width * axes_width_frac:
        axes_width = max_width * axes_width_frac
        axes_height = axes_width / data_aspect

    # floor height
    if axes_height < min_height:
        axes_height = min_height
        axes_width = axes_height * data_aspect

    fig_width = axes_width / axes_width_frac
    fig_height = axes_height + 0.8  # room for title + xlabel

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    return fig, ax


def _filter_abnormal_vectors(
    u: np.ndarray,
    v: np.ndarray,
    candidate_mask: np.ndarray,
    stderr_multiplier: float | None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Remove vectors whose log-magnitude exceeds mean + k * standard error."""
    filtered_mask = np.asarray(candidate_mask, dtype=bool).copy()
    magnitude = np.hypot(u, v)
    positive_mask = filtered_mask & np.isfinite(magnitude) & (magnitude > 0)
    log_magnitudes = np.log(magnitude[positive_mask])
    n_values = int(log_magnitudes.size)

    stats: dict[str, Any] = {
        "discard_abnormal": stderr_multiplier,
        "n_vector_values": n_values,
        "n_vectors_discarded": 0,
        "log_magnitude_mean": None,
        "log_magnitude_stderr": None,
        "log_magnitude_threshold": None,
    }
    if stderr_multiplier is None or n_values == 0:
        return filtered_mask, stats

    log_mean = float(np.mean(log_magnitudes))
    # The sample standard error is undefined for one value. Zero keeps that
    # sole vector because rejection uses a strict greater-than comparison.
    log_stderr = (
        float(np.std(log_magnitudes, ddof=1) / np.sqrt(n_values))
        if n_values > 1
        else 0.0
    )
    threshold = log_mean + stderr_multiplier * log_stderr
    abnormal_mask = np.zeros_like(filtered_mask)
    abnormal_mask[positive_mask] = np.log(magnitude[positive_mask]) > threshold
    filtered_mask[abnormal_mask] = False

    stats.update({
        "n_vectors_discarded": int(np.count_nonzero(abnormal_mask)),
        "log_magnitude_mean": log_mean,
        "log_magnitude_stderr": log_stderr,
        "log_magnitude_threshold": threshold,
    })
    return filtered_mask, stats


def _uniform_vector_components(
    u: np.ndarray,
    v: np.ndarray,
    drawable_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return vector components normalized to unit length where drawable."""
    uniform_u = np.asarray(u, dtype=float).copy()
    uniform_v = np.asarray(v, dtype=float).copy()
    magnitude = np.hypot(uniform_u, uniform_v)
    normalize_mask = np.asarray(drawable_mask, dtype=bool) & (magnitude > 0)
    uniform_u[normalize_mask] /= magnitude[normalize_mask]
    uniform_v[normalize_mask] /= magnitude[normalize_mask]
    return uniform_u, uniform_v


def cmd_plot(args: argparse.Namespace) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from plot_integrated_flux import (
        ADAPTIVE_CELL_DISTANCE_THRESHOLD,
        SPARSE_OCCUPANCY_THRESHOLD,
        plot_threshold_cellmap,
    )

    dataset_dir = Path(args.dataset).resolve()
    out_root = Path(args.out).resolve()

    # Parse planes
    plane_values = [float(p.strip()) for p in args.planes.split(",")]
    color_by_list = [c.strip() for c in args.color_by.split(",")]
    stride = args.stride
    quiver_norm = getattr(args, "quiver_normalized", True)
    quiver_dense = getattr(args, "quiver_dense", True)
    discard_abnormal = getattr(args, "discard_abnormal", None)
    quiver_uniform = getattr(args, "quiver_uniform", True)
    if discard_abnormal is not None and discard_abnormal < 0:
        raise ValueError("--discard-abnormal must be non-negative")
    plane_tol = 1e-6

    # Load manifest for group_name
    manifest_path = dataset_dir / "manifest.yaml"
    group_name = "unknown"
    if manifest_path.exists():
        with open(manifest_path, encoding="utf-8") as f:
            manifest = yaml.safe_load(f)
        group_name = manifest.get("group_name", "unknown")

    # Batch-aware output: {out_root}/{batch_name}/{group_name}/figures/cellmap/
    batch_name = args.batch_name or datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = out_root / batch_name / group_name / "figures" / "cellmap"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print(f"Loading data from {dataset_dir}...")
    df = load_data(dataset_dir)
    print(f"  Loaded {len(df)} rows")

    # Ensure Sx, Sy, Sz are numeric (may have NaN from string 'nan')
    for col in ["Sx", "Sy", "Sz"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Prepare central-period data if requested
    only_central = getattr(args, "only_central", False)
    heatmap_comparison = getattr(args, "heatmap_comparison", False)
    plot_domains: list[tuple[str, pd.DataFrame]] = [("", df)]
    if only_central:
        period = get_period_from_config(group_name)
        central_mask = (df["x"] >= period) & (df["x"] < 2 * period)
        central_df = df[central_mask].copy()
        n_central = len(central_df)
        print(f"  Central period (x in [{period:.3f}, {2*period:.3f})): {n_central}/{len(df)} rows")
        plot_domains.append(("_central", central_df))

    # Build manifest entries
    manifest_entries: list[dict[str, Any]] = []
    warnings: list[str] = []
    total_plots = 0

    for domain_suffix, plot_df, plane_z in (
        (suffix, domain_df, z)
        for suffix, domain_df in plot_domains
        for z in plane_values
    ):
        plane_data = plot_df[np.abs(plot_df["z"] - plane_z) < plane_tol].copy()
        print(f"\n  z = {plane_z}: {len(plane_data)} points")

        if len(plane_data) == 0:
            warnings.append(f"NO_POINTS at z={plane_z}")
            print(f"    WARNING: no points at z={plane_z}")
            continue

        # Check if all flux values are NaN
        all_nan = True
        for col in ["Sx", "Sy", "Sz"]:
            if col in plane_data.columns and plane_data[col].notna().any():
                all_nan = False
                break

        if all_nan:
            warnings.append(f"ALL_VALUES_NAN at z={plane_z}: no valid Poynting vector data")
            print(f"    WARNING: all flux values are NaN at z={plane_z}")

        # Extract unique grid coordinates
        x_unique = sorted(plane_data["x"].unique())
        y_unique = sorted(plane_data["y"].unique())

        # Build index maps
        x_index = {v: i for i, v in enumerate(x_unique)}
        y_index = {v: i for i, v in enumerate(y_unique)}

        # Nonuniform rectilinear grids are valid cellmaps: dense regions get
        # smaller cells, sparse regions get larger cells, and missing cells stay
        # blank. Truly scattered point clouds still fall back to scatter.
        grid_cells = len(x_unique) * len(y_unique)
        is_cellmap_grid = can_use_rectilinear_cellmap(len(x_unique), len(y_unique), len(plane_data))
        occupancy = len(plane_data) / grid_cells if grid_cells else 0.0

        use_threshold_cells = (not is_cellmap_grid) or occupancy < SPARSE_OCCUPANCY_THRESHOLD
        print(
            f"    Grid: {len(x_unique)} x {len(y_unique)}, "
            f"rectilinear_cellmap={is_cellmap_grid}, occupancy={occupancy:.3f}, "
            f"cellmap_method={'threshold_rectangles' if use_threshold_cells else 'pcolormesh'}"
        )

        for color_by in color_by_list:
            # Compute color values
            values_raw = plane_data.apply(lambda row: compute_color_value(row.to_dict(), color_by), axis=1).values
            n_finite = int(np.sum(np.isfinite(values_raw)))

            if is_cellmap_grid and not use_threshold_cells:
                # Build 2D grids
                value_grid = np.full((len(y_unique), len(x_unique)), np.nan)
                Sx_grid = np.full_like(value_grid, np.nan)
                Sy_grid = np.full_like(value_grid, np.nan)

                for idx, (_, row) in enumerate(plane_data.iterrows()):
                    ix = x_index[row["x"]]
                    iy = y_index[row["y"]]
                    value_grid[iy, ix] = values_raw[idx]
                    if "Sx" in row and np.isfinite(row.get("Sx", np.nan)):
                        Sx_grid[iy, ix] = row["Sx"]
                    if "Sy" in row and np.isfinite(row.get("Sy", np.nan)):
                        Sy_grid[iy, ix] = row["Sy"]

                # Compute nonuniform edges from local coordinate midpoints.
                x_edges = axis_cell_edges(np.array(x_unique, dtype=float))
                y_edges = axis_cell_edges(np.array(y_unique, dtype=float))

                # Determine colormap and norms
                if color_by in SIGNED_COLOR_BY:
                    cmap = "RdBu_r"
                    finite_vals = value_grid[np.isfinite(value_grid)]
                    if len(finite_vals) > 0:
                        vmax = float(np.nanpercentile(np.abs(finite_vals), 99))
                        if vmax == 0:
                            vmax = 1.0
                    else:
                        vmax = 1.0
                    vmin = -vmax
                    linear_norm = None  # uses vmin/vmax
                    log_norm = None     # signed → no log
                else:
                    cmap = "viridis"
                    pos_vals = value_grid[(value_grid > 0) & np.isfinite(value_grid)]
                    if len(pos_vals) > 0:
                        log_vmin = float(np.nanpercentile(pos_vals, 1))
                        log_vmax = float(np.nanpercentile(pos_vals, 99))
                        if log_vmin <= 0:
                            log_vmin = log_vmax * 1e-6
                    else:
                        log_vmin, log_vmax = 1e-6, 1.0
                    linear_norm = None
                    log_norm = LogNorm(vmin=log_vmin, vmax=log_vmax)
                    vmin, vmax = None, None

                # Build norm variants: [(label, pcolormesh_kw)]
                norm_variants = [("linear", dict(vmin=vmin, vmax=vmax))]
                if log_norm is not None:
                    norm_variants.append(("log", dict(norm=log_norm)))

                for norm_name, norm_kw in norm_variants:
                    # --- Cellmap only ---
                    fig, axes = _make_fig_for_aspect(x_unique, y_unique)
                    ax = axes
                    im = ax.pcolormesh(
                        x_edges, y_edges, value_grid,
                        cmap=cmap, shading="flat", **norm_kw,
                    )
                    plt.colorbar(im, ax=ax, label=f"{color_by} ({norm_name})")
                    ax.set_xlabel("x")
                    ax.set_ylabel("y")
                    ax.set_title(f"{group_name} z={plane_z:.3f}  ({color_by}, {norm_name})")
                    ax.set_aspect("equal")
                    if all_nan:
                        ax.text(0.5, 0.5, "All values NaN", transform=ax.transAxes,
                                ha="center", va="center", fontsize=16, color="gray")

                    fname = f"{batch_name}_{group_name}_z{plane_z:.3f}_{color_by}_{norm_name}_cellmap{domain_suffix}.png"
                    fig.savefig(out_dir / fname, dpi=150, bbox_inches="tight")
                    plt.close(fig)
                    total_plots += 1

                    manifest_entries.append({
                        "filename": fname, "plane_z": plane_z,
                        "color_by": color_by, "norm": norm_name, "mode": "cellmap",
                        "n_points": len(plane_data), "n_finite": n_finite,
                        "grid_cells": grid_cells, "grid_occupancy": occupancy,
                        "only_central": bool(domain_suffix),
                    })

                    # --- Cellmap + quiver ---
                    fig, axes = _make_fig_for_aspect(x_unique, y_unique)
                    ax = axes
                    im = ax.pcolormesh(
                        x_edges, y_edges, value_grid,
                        cmap=cmap, shading="flat", **norm_kw,
                    )
                    plt.colorbar(im, ax=ax, label=f"{color_by} ({norm_name})")
                    ax.set_xlabel("x")
                    ax.set_ylabel("y")
                    ax.set_title(f"{group_name} z={plane_z:.3f}  ({color_by}, {norm_name} + quiver)")
                    ax.set_aspect("equal")

                    # Quiver overlay
                    X, Y = np.meshgrid(x_unique, y_unique)
                    candidate_mask = (
                        np.isfinite(Sx_grid) & np.isfinite(Sy_grid) & np.isfinite(value_grid)
                    )
                    mask = candidate_mask
                    skip = (slice(None, None, stride), slice(None, None, stride))
                    Xs, Ys = X[skip], Y[skip]
                    Ux, Uy = Sx_grid[skip], Sy_grid[skip]
                    m, vector_filter_stats = _filter_abnormal_vectors(
                        Ux, Uy, mask[skip], discard_abnormal,
                    )
                    quiver_ux, quiver_uy = (
                        _uniform_vector_components(Ux, Uy, m)
                        if quiver_uniform
                        else (Ux, Uy)
                    )
                    if np.any(m):
                        ax.quiver(
                            Xs[m], Ys[m], quiver_ux[m], quiver_uy[m],
                            scale=None, alpha=0.7, color="k", width=0.003,
                        )

                    if all_nan:
                        ax.text(0.5, 0.5, "All values NaN", transform=ax.transAxes,
                                ha="center", va="center", fontsize=16, color="gray")

                    fname = f"{batch_name}_{group_name}_z{plane_z:.3f}_{color_by}_{norm_name}_cellmap_quiver{domain_suffix}.png"
                    fig.savefig(out_dir / fname, dpi=150, bbox_inches="tight")
                    plt.close(fig)
                    total_plots += 1

                    manifest_entries.append({
                        "filename": fname, "plane_z": plane_z,
                        "color_by": color_by, "norm": norm_name, "mode": "cellmap_quiver",
                        "n_points": len(plane_data), "n_finite": n_finite,
                        "grid_cells": grid_cells, "grid_occupancy": occupancy,
                        "only_central": bool(domain_suffix),
                        **vector_filter_stats,
                    })

                    # --- Dense physical quiver (2x density) ---
                    if quiver_dense and np.any(mask):
                        dense_stride = max(stride // 2, 1)
                        skip_d = (slice(None, None, dense_stride), slice(None, None, dense_stride))
                        Xd, Yd = X[skip_d], Y[skip_d]
                        Ud, Vd = Sx_grid[skip_d], Sy_grid[skip_d]
                        md = mask[skip_d]
                        md, dense_vector_filter_stats = _filter_abnormal_vectors(
                            Ud, Vd, md, discard_abnormal,
                        )
                        quiver_ud, quiver_vd = (
                            _uniform_vector_components(Ud, Vd, md)
                            if quiver_uniform
                            else (Ud, Vd)
                        )

                        fig, axes = _make_fig_for_aspect(x_unique, y_unique)
                        ax = axes
                        im = ax.pcolormesh(
                            x_edges, y_edges, value_grid,
                            cmap=cmap, shading="flat", **norm_kw,
                        )
                        plt.colorbar(im, ax=ax, label=f"{color_by} ({norm_name})")
                        ax.set_xlabel("x")
                        ax.set_ylabel("y")
                        ax.set_title(f"{group_name} z={plane_z:.3f}  ({color_by}, {norm_name} + quiver2x)")
                        ax.set_aspect("equal")
                        if np.any(md):
                            ax.quiver(
                                Xd[md], Yd[md], quiver_ud[md], quiver_vd[md],
                                scale=None, alpha=0.7, color="k", width=0.003,
                            )
                        if all_nan:
                            ax.text(0.5, 0.5, "All values NaN", transform=ax.transAxes,
                                    ha="center", va="center", fontsize=16, color="gray")

                        fname = f"{batch_name}_{group_name}_z{plane_z:.3f}_{color_by}_{norm_name}_cellmap_quiver2x{domain_suffix}.png"
                        fig.savefig(out_dir / fname, dpi=150, bbox_inches="tight")
                        plt.close(fig)
                        total_plots += 1

                        manifest_entries.append({
                            "filename": fname, "plane_z": plane_z,
                            "color_by": color_by, "norm": norm_name, "mode": "cellmap_quiver2x",
                            "n_points": len(plane_data), "n_finite": n_finite,
                            "only_central": bool(domain_suffix),
                            **dense_vector_filter_stats,
                        })

                    # --- Normalized quiver (uniform arrow length) ---
                    if quiver_norm and np.any(mask):
                        fig, axes = _make_fig_for_aspect(x_unique, y_unique)
                        ax = axes
                        im = ax.pcolormesh(
                            x_edges, y_edges, value_grid,
                            cmap=cmap, shading="flat", **norm_kw,
                        )
                        plt.colorbar(im, ax=ax, label=f"{color_by} ({norm_name})")
                        ax.set_xlabel("x")
                        ax.set_ylabel("y")
                        ax.set_title(f"{group_name} z={plane_z:.3f}  ({color_by}, {norm_name} + normquiver)")
                        ax.set_aspect("equal")

                        Uxn = Ux.copy()
                        Uyn = Uy.copy()
                        mag = np.sqrt(Uxn**2 + Uyn**2)
                        ok = (mag > 0) & m
                        Uxn[ok] /= mag[ok]
                        Uyn[ok] /= mag[ok]
                        if np.any(ok):
                            ax.quiver(
                                Xs[ok], Ys[ok], Uxn[ok], Uyn[ok],
                                scale=150, alpha=0.7, color="k", width=0.003,
                            )

                        if all_nan:
                            ax.text(0.5, 0.5, "All values NaN", transform=ax.transAxes,
                                    ha="center", va="center", fontsize=16, color="gray")

                        fname = f"{batch_name}_{group_name}_z{plane_z:.3f}_{color_by}_{norm_name}_cellmap_normquiver{domain_suffix}.png"
                        fig.savefig(out_dir / fname, dpi=150, bbox_inches="tight")
                        plt.close(fig)
                        total_plots += 1

                        manifest_entries.append({
                            "filename": fname, "plane_z": plane_z,
                            "color_by": color_by, "norm": norm_name, "mode": "cellmap_normquiver",
                            "n_points": len(plane_data), "n_finite": n_finite,
                            "only_central": bool(domain_suffix),
                            **vector_filter_stats,
                        })

                        # --- Dense normalized quiver (2x density) ---
                        if quiver_dense:
                            dense_stride = max(stride // 2, 1)
                            skip_d = (slice(None, None, dense_stride), slice(None, None, dense_stride))
                            Xd, Yd = X[skip_d], Y[skip_d]
                            Ud, Vd = Sx_grid[skip_d], Sy_grid[skip_d]
                            md = mask[skip_d]
                            md, dense_vector_filter_stats = _filter_abnormal_vectors(
                                Ud, Vd, md, discard_abnormal,
                            )
                            Udn = Ud.copy()
                            Vdn = Vd.copy()
                            mag_d = np.sqrt(Udn**2 + Vdn**2)
                            ok_d = (mag_d > 0) & md
                            Udn[ok_d] /= mag_d[ok_d]
                            Vdn[ok_d] /= mag_d[ok_d]

                            fig, axes = _make_fig_for_aspect(x_unique, y_unique)
                            ax = axes
                            im = ax.pcolormesh(
                                x_edges, y_edges, value_grid,
                                cmap=cmap, shading="flat", **norm_kw,
                            )
                            plt.colorbar(im, ax=ax, label=f"{color_by} ({norm_name})")
                            ax.set_xlabel("x")
                            ax.set_ylabel("y")
                            ax.set_title(f"{group_name} z={plane_z:.3f}  ({color_by}, {norm_name} + normquiver2x)")
                            ax.set_aspect("equal")
                            if np.any(ok_d):
                                ax.quiver(
                                    Xd[ok_d], Yd[ok_d], Udn[ok_d], Vdn[ok_d],
                                    scale=150, alpha=0.7, color="k", width=0.003,
                                )
                            if all_nan:
                                ax.text(0.5, 0.5, "All values NaN", transform=ax.transAxes,
                                        ha="center", va="center", fontsize=16, color="gray")

                            fname = f"{batch_name}_{group_name}_z{plane_z:.3f}_{color_by}_{norm_name}_cellmap_normquiver2x{domain_suffix}.png"
                            fig.savefig(out_dir / fname, dpi=150, bbox_inches="tight")
                            plt.close(fig)
                            total_plots += 1

                            manifest_entries.append({
                                "filename": fname, "plane_z": plane_z,
                                "color_by": color_by, "norm": norm_name, "mode": "cellmap_normquiver2x",
                                "n_points": len(plane_data), "n_finite": n_finite,
                                "only_central": bool(domain_suffix),
                                **dense_vector_filter_stats,
                            })

            else:
                # Sparse/nonrectangular fallback: axis-aligned cells bridge
                # only to same-row/same-column neighbors within the threshold.
                # Large gaps stay blank with rectangular boundaries.
                if color_by in SIGNED_COLOR_BY:
                    cmap = "RdBu_r"
                    finite_vals = values_raw[np.isfinite(values_raw)]
                    if len(finite_vals) > 0:
                        vmax = float(np.nanpercentile(np.abs(finite_vals), 99))
                        if vmax == 0:
                            vmax = 1.0
                    else:
                        vmax = 1.0
                    vmin = -vmax
                    scatter_norm_variants = [("linear", dict(vmin=vmin, vmax=vmax))]
                else:
                    cmap = "viridis"
                    pos_vals = values_raw[(values_raw > 0) & np.isfinite(values_raw)]
                    if len(pos_vals) > 0:
                        lv = float(np.nanpercentile(pos_vals, 1))
                        uv = float(np.nanpercentile(pos_vals, 99))
                        if lv <= 0:
                            lv = uv * 1e-6
                    else:
                        lv, uv = 1e-6, 1.0
                    scatter_norm_variants = [
                        ("linear", dict(vmin=None, vmax=None)),
                        ("log", dict(norm=LogNorm(vmin=lv, vmax=uv))),
                    ]

                for norm_name, norm_kw in scatter_norm_variants:
                    fig, axes = _make_fig_for_aspect(x_unique, y_unique)
                    ax = axes
                    im = plot_threshold_cellmap(
                        ax,
                        plane_data["x"].to_numpy(dtype=float),
                        plane_data["y"].to_numpy(dtype=float),
                        values_raw,
                        cmap=cmap,
                        norm_kw=norm_kw,
                    )
                    plt.colorbar(im, ax=ax, label=f"{color_by} ({norm_name})")
                    ax.set_xlabel("x")
                    ax.set_ylabel("y")
                    ax.set_title(f"{group_name} z={plane_z:.3f}  ({color_by}, {norm_name} threshold cells)")
                    ax.set_aspect("equal")
                    if all_nan:
                        ax.text(0.5, 0.5, "All values NaN", transform=ax.transAxes,
                                ha="center", va="center", fontsize=16, color="gray")

                    fname = f"{batch_name}_{group_name}_z{plane_z:.3f}_{color_by}_{norm_name}_thresholdcells{domain_suffix}.png"
                    fig.savefig(out_dir / fname, dpi=150, bbox_inches="tight")
                    plt.close(fig)
                    total_plots += 1

                    manifest_entries.append({
                        "filename": fname, "plane_z": plane_z,
                        "color_by": color_by, "norm": norm_name, "mode": "threshold_rectangles",
                        "n_points": len(plane_data), "n_finite": n_finite,
                        "grid_cells": grid_cells, "grid_occupancy": occupancy,
                        "adaptive_cell_distance_threshold": ADAPTIVE_CELL_DISTANCE_THRESHOLD,
                        "only_central": bool(domain_suffix),
                    })

                    fig, axes = _make_fig_for_aspect(x_unique, y_unique)
                    ax = axes
                    im = plot_threshold_cellmap(
                        ax,
                        plane_data["x"].to_numpy(dtype=float),
                        plane_data["y"].to_numpy(dtype=float),
                        values_raw,
                        cmap=cmap,
                        norm_kw=norm_kw,
                    )
                    plt.colorbar(im, ax=ax, label=f"{color_by} ({norm_name})")
                    ax.set_xlabel("x")
                    ax.set_ylabel("y")
                    ax.set_title(f"{group_name} z={plane_z:.3f}  ({color_by}, {norm_name} threshold cells + quiver)")
                    ax.set_aspect("equal")

                    step = max(int(stride), 1)
                    qx = plane_data["x"].to_numpy(dtype=float)[::step]
                    qy = plane_data["y"].to_numpy(dtype=float)[::step]
                    qv = values_raw[::step]
                    qu = plane_data["Sx"].to_numpy(dtype=float)[::step]
                    qvv = plane_data["Sy"].to_numpy(dtype=float)[::step]
                    candidate = np.isfinite(qx) & np.isfinite(qy) & np.isfinite(qv) & np.isfinite(qu) & np.isfinite(qvv)
                    candidate, vector_filter_stats = _filter_abnormal_vectors(
                        qu, qvv, candidate, discard_abnormal,
                    )
                    quiver_u, quiver_v = (
                        _uniform_vector_components(qu, qvv, candidate)
                        if quiver_uniform
                        else (qu, qvv)
                    )
                    if np.any(candidate):
                        ax.quiver(
                            qx[candidate], qy[candidate], quiver_u[candidate], quiver_v[candidate],
                            scale=None, alpha=0.7, color="k", width=0.003,
                        )

                    if all_nan:
                        ax.text(0.5, 0.5, "All values NaN", transform=ax.transAxes,
                                ha="center", va="center", fontsize=16, color="gray")

                    fname = f"{batch_name}_{group_name}_z{plane_z:.3f}_{color_by}_{norm_name}_thresholdcells_quiver{domain_suffix}.png"
                    fig.savefig(out_dir / fname, dpi=150, bbox_inches="tight")
                    plt.close(fig)
                    total_plots += 1

                    manifest_entries.append({
                        "filename": fname, "plane_z": plane_z,
                        "color_by": color_by, "norm": norm_name, "mode": "threshold_rectangles_quiver",
                        "n_points": len(plane_data), "n_finite": n_finite,
                        "grid_cells": grid_cells, "grid_occupancy": occupancy,
                        "adaptive_cell_distance_threshold": ADAPTIVE_CELL_DISTANCE_THRESHOLD,
                        "only_central": bool(domain_suffix),
                        **vector_filter_stats,
                    })

        # --- Heatmap comparison (side-by-side: quiver | no quiver) ---
        if heatmap_comparison and is_cellmap_grid and not use_threshold_cells and len(plane_data) > 0 and not all_nan:
            # Build absSxy grid from raw Sx, Sy
            Sx_g = np.full((len(y_unique), len(x_unique)), np.nan)
            Sy_g = np.full_like(Sx_g, np.nan)
            for _, row in plane_data.iterrows():
                ix = x_index[row["x"]]
                iy = y_index[row["y"]]
                sx_v = row.get("Sx", np.nan)
                sy_v = row.get("Sy", np.nan)
                Sx_g[iy, ix] = sx_v if np.isfinite(sx_v) else np.nan
                Sy_g[iy, ix] = sy_v if np.isfinite(sy_v) else np.nan
            absSxy_g = np.sqrt(Sx_g**2 + Sy_g**2)

            # Log norm for absSxy
            pos_g = absSxy_g[(absSxy_g > 0) & np.isfinite(absSxy_g)]
            if len(pos_g) > 0:
                lv = float(np.nanpercentile(pos_g, 1))
                uv = float(np.nanpercentile(pos_g, 99))
                lv = lv if lv > 0 else uv * 1e-6
            else:
                lv, uv = 1e-6, 1.0
            cmp_log_norm = LogNorm(vmin=lv, vmax=uv)

            fig, (ax_l, ax_r) = plt.subplots(
                1, 2, figsize=(18, 7),
                gridspec_kw=dict(width_ratios=[1, 1], wspace=0.3),
            )

            # Shared data for both panels: heatmap (log absSxy)
            im = ax_l.pcolormesh(x_edges, y_edges, absSxy_g, cmap="afmhot",
                                 shading="flat", norm=cmp_log_norm)
            ax_r.pcolormesh(x_edges, y_edges, absSxy_g, cmap="afmhot",
                            shading="flat", norm=cmp_log_norm)

            # Left: quiver overlay (log-scale arrows via normalization)
            Xg, Yg = np.meshgrid(x_unique, y_unique)
            skip_cmp = (slice(None, None, stride), slice(None, None, stride))
            Xc, Yc = Xg[skip_cmp], Yg[skip_cmp]
            Uc, Vc = Sx_g[skip_cmp], Sy_g[skip_cmp]
            mag_c = np.sqrt(Uc**2 + Vc**2)
            ok_c = (mag_c > 0) & np.isfinite(Uc) & np.isfinite(Vc)
            ok_c, comparison_filter_stats = _filter_abnormal_vectors(
                Uc, Vc, ok_c, discard_abnormal,
            )
            quiver_uc, quiver_vc = (
                _uniform_vector_components(Uc, Vc, ok_c)
                if quiver_uniform
                else (Uc, Vc)
            )
            if np.any(ok_c):
                ax_l.quiver(Xc[ok_c], Yc[ok_c], quiver_uc[ok_c], quiver_vc[ok_c],
                            scale=None, alpha=0.7, color="k", width=0.003)

            plt.colorbar(im, ax=ax_l, label="absSxy (log)")
            ax_l.set_xlabel("x"); ax_l.set_ylabel("y")
            ax_l.set_title(f"{group_name} z={plane_z:.3f} absSxy (log) + quiver")
            ax_l.set_aspect("equal")

            plt.colorbar(im, ax=ax_r, label="absSxy (log)")
            ax_r.set_xlabel("x"); ax_r.set_ylabel("y")
            ax_r.set_title(f"{group_name} z={plane_z:.3f} absSxy (log)")
            ax_r.set_aspect("equal")

            fname = f"{batch_name}_{group_name}_z{plane_z:.3f}_absSxy_log_heatmap_comparison{domain_suffix}.png"
            fig.savefig(out_dir / fname, dpi=150, bbox_inches="tight")
            plt.close(fig)
            total_plots += 1

            manifest_entries.append({
                "filename": fname, "plane_z": plane_z,
                "color_by": "absSxy", "norm": "log", "mode": "heatmap_comparison",
                "n_points": len(plane_data),
                "n_finite": int(np.sum(np.isfinite(absSxy_g))),
                "only_central": bool(domain_suffix),
                **comparison_filter_stats,
            })

    # Write plot manifest
    plot_manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": str(dataset_dir),
        "planes": plane_values,
        "mode": args.mode,
        "color_by": args.color_by,
        "stride": stride,
        "discard_abnormal": discard_abnormal,
        "quiver_uniform": quiver_uniform,
        "n_plots": total_plots,
        "plots": manifest_entries,
        "warnings": warnings,
    }
    manifest_out = out_dir / "plot_manifest.yaml"
    manifest_out.write_text(
        yaml.dump(plot_manifest, default_flow_style=False, sort_keys=False),
        encoding="utf-8",
    )

    print(f"\nPlot complete:")
    print(f"  {total_plots} figures generated")
    print(f"  Output: {out_dir}")
    if warnings:
        print(f"  Warnings:")
        for w in warnings:
            print(f"    {w}")


# ---------------------------------------------------------------------------
# Auto-discovery
# ---------------------------------------------------------------------------

def discover_config_names(task_dir: Path) -> list[str]:
    """Find unique config group names from subtask-* directories."""
    names: set[str] = set()
    for p in task_dir.iterdir():
        if p.is_dir():
            parsed = parse_subtask_name(p.name)
            if parsed is not None:
                names.add(parsed[0])
    return sorted(names)


# ---------------------------------------------------------------------------
# Run-all command
# ---------------------------------------------------------------------------

def cmd_run_all(args: argparse.Namespace) -> None:
    task_dir = Path(args.task_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    planes = args.planes
    color_by = args.color_by
    stride = args.stride
    last_only = args.last
    only_central = getattr(args, "only_central", False)
    heatmap_comparison = getattr(args, "heatmap_comparison", False)
    discard_abnormal = getattr(args, "discard_abnormal", None)
    quiver_uniform = getattr(args, "quiver_uniform", True)

    if not task_dir.is_dir():
        print(f"ERROR: task directory not found: {task_dir}", file=sys.stderr)
        sys.exit(1)

    # Auto-discover configs
    config_names = discover_config_names(task_dir)
    if not config_names:
        print(f"ERROR: no subtask-* directories found in {task_dir}", file=sys.stderr)
        sys.exit(1)

    # --- Batch timestamp directory (top-level) ---
    timestr = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    batch_dir = out_dir / timestr
    batch_dir.mkdir(parents=True, exist_ok=True)

    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Input:  {task_dir}")
    print(f"Output: {out_dir}")
    print(f"Batch:  {batch_dir}")
    print(f"Configs: {', '.join(config_names)}")
    print()

    # --- Phase 1: SRFlux collect + plot (spatial data → cellmap/quiver) ---
    print("--- Phase 1: SRFlux (spatially-resolved) ---")

    for cfg in config_names:
        case_dir = batch_dir / cfg
        data_dir = case_dir / "data"
        data_dir.mkdir(parents=True, exist_ok=True)

        # Collect
        print(f"=== [{cfg}] Collect ===")
        collect_args = argparse.Namespace(
            task_dir=str(task_dir),
            group_name=cfg,
            filebase="task",
            out=str(data_dir),
            coord_tol=1e-9,
            duplicate_policy="error",
            last=last_only,
        )
        cmd_collect(collect_args)

        # Plot (spatial: cellmap + quiver) — now uses batch-aware paths
        print(f"=== [{cfg}] Plot ===")
        plot_args = argparse.Namespace(
            dataset=str(data_dir),
            planes=planes,
            mode="cellmap-quiver",
            color_by=color_by,
            vector="Sx,Sy",
            source_mode="sum",
            stride=stride,
            quiver_normalized=args.quiver_normalized,
            quiver_dense=args.quiver_dense,
            only_central=only_central,
            heatmap_comparison=heatmap_comparison,
            discard_abnormal=discard_abnormal,
            quiver_uniform=quiver_uniform,
            batch_name=timestr,
            out=str(out_dir),
        )
        cmd_plot(plot_args)
        print()

    # --- Phase 2: SIFlux analysis (scalar PAbs/PRad vs omega) ---
    siflux_dir = batch_dir / "siflux"
    siflux_dir.mkdir(parents=True, exist_ok=True)
    print("--- Phase 2: SIFlux (spatially-integrated, scalar) ---")

    postproc_script = Path(__file__).parent / "postproc_scuffem.py"
    if not postproc_script.exists():
        print(f"WARNING: postproc_scuffem.py not found at {postproc_script}, skipping SIFlux analysis",
              file=sys.stderr)
    else:
        cmd = [sys.executable, str(postproc_script), str(task_dir),
               "--out-dir", str(siflux_dir)]
        if last_only:
            cmd.append("--last")
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=False)
        if result.returncode != 0:
            print(f"WARNING: postproc_scuffem.py exited with code {result.returncode}", file=sys.stderr)

    print()
    print("=" * 60)
    print("Run-all complete.")
    first_cfg = config_names[0] if config_names else "?"
    print(f"  Batch:     {batch_dir}/")
    print(f"    SRFlux (spatial): {batch_dir}/{first_cfg}/data + {batch_dir}/{first_cfg}/figures/cellmap/")
    print(f"    SIFlux (scalar):  {siflux_dir}/")
    print(f"  Configs: {', '.join(config_names)}")


# ---------------------------------------------------------------------------
# Slice helpers
# ---------------------------------------------------------------------------

def get_period_from_config(config_name: str) -> float:
    """Extract period p from config name (e.g., 'ad1' -> group 'a' -> p=3.2)."""
    group = config_name[0]
    if group not in GROUP_DELTA:
        print(f"ERROR: unknown group '{group}' from config '{config_name}'", file=sys.stderr)
        sys.exit(1)
    return 2.0 * (W_GEOM + GROUP_DELTA[group])


def bose_einstein_weight(omega: float, temperature: float) -> float:
    """Compute Bose-Einstein occupation number n(omega, T)."""
    if omega <= 0:
        return 0.0
    x = HBAR * omega / (K_B * temperature)
    if x > 500:
        return 0.0
    if x < 1e-15:
        return 1.0 / x if x > 0 else 0.0
    return 1.0 / (np.exp(x) - 1.0)


SLICE_QUANTITY_OPTIONS = [
    "Sx", "Sy", "Sz", "absSxy", "absS", "log_absSxy", "log_absS",
    "angle_xy", "sin_angle_xy", "cos_angle_xy",
]

QUANTITY_HELP = {
    "Sx": "Poynting vector x-component",
    "Sy": "Poynting vector y-component",
    "Sz": "Poynting vector z-component",
    "absSxy": "|S_xy| = sqrt(Sx^2 + Sy^2)",
    "absS": "|S| = sqrt(Sx^2 + Sy^2 + Sz^2)",
    "log_absSxy": "log10(|S_xy|)",
    "log_absS": "log10(|S|)",
    "angle_xy": "atan2(Sy, Sx) in radians — Poynting direction angle in xy-plane",
    "sin_angle_xy": "sin(atan2(Sy, Sx)) = Sy / |S_xy|",
    "cos_angle_xy": "cos(atan2(Sy, Sx)) = Sx / |S_xy|",
}


def compute_quantity_vectorized(Sx: np.ndarray, Sy: np.ndarray, Sz: np.ndarray, quantity: str) -> np.ndarray:
    """Vectorized quantity computation for arrays."""
    if quantity == "Sx":
        return Sx
    elif quantity == "Sy":
        return Sy
    elif quantity == "Sz":
        return Sz
    elif quantity == "absSxy":
        return np.sqrt(Sx**2 + Sy**2)
    elif quantity == "absS":
        return np.sqrt(Sx**2 + Sy**2 + Sz**2)
    elif quantity == "log_absSxy":
        return np.log10(np.sqrt(Sx**2 + Sy**2) + 1e-30)
    elif quantity == "log_absS":
        return np.log10(np.sqrt(Sx**2 + Sy**2 + Sz**2) + 1e-30)
    elif quantity == "angle_xy":
        return np.arctan2(Sy, Sx)
    elif quantity == "sin_angle_xy":
        mag = np.sqrt(Sx**2 + Sy**2)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(mag > 1e-30, Sy / mag, 0.0)
    elif quantity == "cos_angle_xy":
        mag = np.sqrt(Sx**2 + Sy**2)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(mag > 1e-30, Sx / mag, 0.0)
    else:
        return np.sqrt(Sx**2 + Sy**2)


def load_srflux_with_omega(dataset_dir: Path) -> pd.DataFrame:
    """Load merged_srflux.tsv preserving per-omega rows."""
    srflux_path = dataset_dir / "merged_srflux.tsv"
    if not srflux_path.exists():
        print(f"ERROR: merged_srflux.tsv not found in {dataset_dir}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(srflux_path, sep="\t")
    for col in ["omega", "x", "y", "z", "Sx_flux", "Sy_flux", "Sz_flux"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def filter_slice_points(
    df: pd.DataFrame,
    slice_axis: str,
    slice_value: float,
    tolerance: float,
    period: float | None = None,
) -> pd.DataFrame:
    """Select rows whose coordinates match the slice specification."""
    if slice_axis == "y":
        mask = np.abs(df["y"].values - slice_value) < tolerance
        return df[mask].copy()
    elif slice_axis == "x":
        if period is None:
            print("ERROR: period is required for x-slice", file=sys.stderr)
            sys.exit(1)
        x_mod = np.mod(slice_value, period)
        x_fold = np.mod(df["x"].values, period)
        mask = np.abs(x_fold - x_mod) < tolerance
        return df[mask].copy()
    else:
        print(f"ERROR: invalid slice_axis '{slice_axis}'", file=sys.stderr)
        sys.exit(1)


def compute_weighted_sum(df: pd.DataFrame, temperature: float) -> pd.DataFrame:
    """Apply Bose-Einstein weighting and sum flux across omegas for each spatial point."""
    omega_unique = df["omega"].dropna().unique()
    weight_map = {w: bose_einstein_weight(w, temperature) for w in omega_unique}

    df = df.copy()
    df["weight"] = df["omega"].map(weight_map)
    df["Sx_w"] = df["Sx_flux"].fillna(0) * df["weight"]
    df["Sy_w"] = df["Sy_flux"].fillna(0) * df["weight"]
    df["Sz_w"] = df["Sz_flux"].fillna(0) * df["weight"]

    grouped = df.groupby(["x", "y", "z"], as_index=False).agg(
        {"Sx_w": "sum", "Sy_w": "sum", "Sz_w": "sum"}
    )
    grouped = grouped.rename(columns={"Sx_w": "Sx", "Sy_w": "Sy", "Sz_w": "Sz"})
    return grouped


def check_periodicity(
    df: pd.DataFrame,
    slice_axis: str,
    slice_value: float,
    tolerance: float,
    period: float,
    quantity: str,
    temperature: float,
) -> dict:
    """Compare values at base position vs period-shifted positions.

    Returns dict with per-omega and summed comparison results.
    """
    result: dict[str, Any] = {"period_offsets": [], "per_omega": {}, "summed": {}}

    if slice_axis == "y":
        # Pick a reference x in the first period from the filtered data
        df_base = df[(df["x"] >= 0) & (df["x"] < period)]
        if df_base.empty:
            return result
        x_ref = float(np.median(df_base["x"].values))

        offsets = [x_ref + k * period for k in range(3)]
        result["period_offsets"] = offsets

        # Per-omega comparison
        omega_unique = sorted(df["omega"].dropna().unique())
        for omega_val in omega_unique:
            df_w = df[df["omega"] == omega_val]
            vals = []
            for x_target in offsets:
                subset = df_w[np.abs(df_w["x"].values - x_target) < tolerance]
                if subset.empty:
                    vals.append(np.nan)
                else:
                    q = compute_quantity_vectorized(
                        subset["Sx_flux"].fillna(0).values,
                        subset["Sy_flux"].fillna(0).values,
                        subset["Sz_flux"].fillna(0).values,
                        quantity,
                    )
                    vals.append(float(np.nanmean(q)))
            ref = vals[0] if np.isfinite(vals[0]) else 0.0
            eps = max(abs(ref), 1e-30)
            rel_diffs = [0.0] + [abs(v - ref) / eps if np.isfinite(v) else np.nan for v in vals[1:]]
            result["per_omega"][omega_val] = {
                "quantity_values": vals,
                "relative_diffs": rel_diffs,
            }

        # Summed comparison
        df_sum = compute_weighted_sum(df, temperature)
        vals_sum = []
        for x_target in offsets:
            subset = df_sum[np.abs(df_sum["x"].values - x_target) < tolerance]
            if subset.empty:
                vals_sum.append(np.nan)
            else:
                q = compute_quantity_vectorized(
                    subset["Sx"].values, subset["Sy"].values, subset["Sz"].values, quantity
                )
                vals_sum.append(float(np.nanmean(q)))
        ref_sum = vals_sum[0] if np.isfinite(vals_sum[0]) else 0.0
        eps_sum = max(abs(ref_sum), 1e-30)
        result["summed"] = {
            "quantity_values": vals_sum,
            "relative_diffs": [0.0] + [abs(v - ref_sum) / eps_sum if np.isfinite(v) else np.nan for v in vals_sum[1:]],
        }

    elif slice_axis == "x":
        x_mod = np.mod(slice_value, period)
        offsets = [x_mod + k * period for k in range(3)]
        result["period_offsets"] = offsets

        omega_unique = sorted(df["omega"].dropna().unique())
        for omega_val in omega_unique:
            df_w = df[df["omega"] == omega_val]
            vals = []
            for x_target in offsets:
                subset = df_w[np.abs(df_w["x"].values - x_target) < tolerance]
                if subset.empty:
                    vals.append(np.nan)
                else:
                    q = compute_quantity_vectorized(
                        subset["Sx_flux"].fillna(0).values,
                        subset["Sy_flux"].fillna(0).values,
                        subset["Sz_flux"].fillna(0).values,
                        quantity,
                    )
                    vals.append(float(np.nanmean(q)))
            ref = vals[0] if np.isfinite(vals[0]) else 0.0
            eps = max(abs(ref), 1e-30)
            rel_diffs = [0.0] + [abs(v - ref) / eps if np.isfinite(v) else np.nan for v in vals[1:]]
            result["per_omega"][omega_val] = {
                "quantity_values": vals,
                "relative_diffs": rel_diffs,
            }

        df_sum = compute_weighted_sum(df, temperature)
        vals_sum = []
        for x_target in offsets:
            subset = df_sum[np.abs(df_sum["x"].values - x_target) < tolerance]
            if subset.empty:
                vals_sum.append(np.nan)
            else:
                q = compute_quantity_vectorized(
                    subset["Sx"].values, subset["Sy"].values, subset["Sz"].values, quantity
                )
                vals_sum.append(float(np.nanmean(q)))
        ref_sum = vals_sum[0] if np.isfinite(vals_sum[0]) else 0.0
        eps_sum = max(abs(ref_sum), 1e-30)
        result["summed"] = {
            "quantity_values": vals_sum,
            "relative_diffs": [0.0] + [abs(v - ref_sum) / eps_sum if np.isfinite(v) else np.nan for v in vals_sum[1:]],
        }

    return result


# ---------------------------------------------------------------------------
# Slice export
# ---------------------------------------------------------------------------

SLICE_RAW_COLS = ["omega", "x", "y", "z", "Sx_flux", "Sy_flux", "Sz_flux", "weight", "Sx_w", "Sy_w", "Sz_w"]
SLICE_WEIGHTED_COLS = ["x", "y", "z", "Sx", "Sy", "Sz"]
SLICE_PER_OMEGA_COLS = ["omega", "x", "y", "z", "Sx_w", "Sy_w", "Sz_w"]


def slice_export_data(
    dataset_dir: Path,
    config_name: str,
    slice_axis: str,
    slice_value: float,
    temperature: float,
    tolerance: float,
    plane_values: list[float],
    quantities: list[str],
    out_dir: Path,
) -> dict[str, Any]:
    """Process slice data and export to TSV + YAML.

    Returns metadata dict for downstream use.
    """
    period = get_period_from_config(config_name)

    if slice_axis == "x":
        folded_value = np.mod(slice_value, period)
    else:
        folded_value = slice_value

    # Load and filter
    df = load_srflux_with_omega(dataset_dir)
    df = df[df["z"].isin(plane_values)]
    df_slice = filter_slice_points(df, slice_axis, folded_value, tolerance, period)

    if df_slice.empty:
        print("WARNING: no points match the slice criteria", file=sys.stderr)
        return {"status": "empty", "n_rows": 0}

    # Compute weights and weighted flux
    omega_unique = df_slice["omega"].dropna().unique()
    weight_map = {w: bose_einstein_weight(w, temperature) for w in omega_unique}

    df_slice = df_slice.copy()
    df_slice["weight"] = df_slice["omega"].map(weight_map)
    df_slice["Sx_w"] = df_slice["Sx_flux"].fillna(0) * df_slice["weight"]
    df_slice["Sy_w"] = df_slice["Sy_flux"].fillna(0) * df_slice["weight"]
    df_slice["Sz_w"] = df_slice["Sz_flux"].fillna(0) * df_slice["weight"]

    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Write slice_raw.tsv ---
    raw_lines = ["\t".join(SLICE_RAW_COLS)]
    for _, row in df_slice.iterrows():
        parts = []
        for col in SLICE_RAW_COLS:
            val = row[col]
            if col == "omega":
                parts.append(f"{val:.10e}")
            elif isinstance(val, (int, float)):
                parts.append(f"{val:.10e}")
            else:
                parts.append(str(val))
        raw_lines.append("\t".join(parts))
    write_text_file(out_dir / "slice_raw.tsv", "\n".join(raw_lines))

    # --- Write slice_weighted.tsv ---
    df_sum = compute_weighted_sum(df_slice, temperature)
    weighted_lines = ["\t".join(SLICE_WEIGHTED_COLS)]
    for _, row in df_sum.iterrows():
        parts = []
        for col in SLICE_WEIGHTED_COLS:
            val = row[col]
            if isinstance(val, (int, float)):
                parts.append(f"{val:.10e}")
            else:
                parts.append(str(val))
        weighted_lines.append("\t".join(parts))
    write_text_file(out_dir / "slice_weighted.tsv", "\n".join(weighted_lines))

    # --- Write slice_per_omega.tsv ---
    per_omega_data: dict[float, dict[str, Any]] = {}
    omega_sorted = sorted(omega_unique)
    per_omega_lines = ["\t".join(SLICE_PER_OMEGA_COLS)]
    for omega_val in omega_sorted:
        df_w = df_slice[df_slice["omega"] == omega_val]
        if df_w.empty:
            continue
        for _, row in df_w.iterrows():
            parts = [f"{omega_val:.10e}"]
            for col in ["x", "y", "z", "Sx_w", "Sy_w", "Sz_w"]:
                parts.append(f"{row[col]:.10e}")
            per_omega_lines.append("\t".join(parts))
    write_text_file(out_dir / "slice_per_omega.tsv", "\n".join(per_omega_lines))

    # --- Write slice_meta.yaml ---
    meta = {
        "status": "ok",
        "config_name": config_name,
        "period": float(period),
        "slice_axis": slice_axis,
        "slice_value": float(slice_value),
        "folded_value": float(folded_value),
        "temperature": float(temperature),
        "tolerance": float(tolerance),
        "plane_values": [float(v) for v in plane_values],
        "quantities": quantities,
        "n_omegas": int(len(omega_sorted)),
        "n_raw_rows": int(len(df_slice)),
        "n_weighted_points": int(len(df_sum)),
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    (out_dir / "slice_meta.yaml").write_text(
        yaml.dump(meta, default_flow_style=False, sort_keys=False),
        encoding="utf-8",
    )

    print(f"Slice export: {len(df_slice)} raw rows, {len(df_sum)} weighted points, {len(omega_sorted)} omegas")
    print(f"  Output: {out_dir}")

    return meta


def config_from_dataset_path(dataset_dir: Path) -> str:
    """Derive config name from dataset directory path (e.g., '.../srflux/ad1' -> 'ad1')."""
    name = dataset_dir.name
    if name.startswith("pv_out_"):
        name = name[len("pv_out_"):]
    return name


def validate_quantities(quantities: list[str]) -> list[str]:
    """Validate quantity names, print help and exit if unknown."""
    valid = set(SLICE_QUANTITY_OPTIONS)
    for q in quantities:
        if q not in valid:
            print(f"ERROR: unknown quantity '{q}'", file=sys.stderr)
            print(f"Available quantities:", file=sys.stderr)
            for name in SLICE_QUANTITY_OPTIONS:
                print(f"  {name:20s} — {QUANTITY_HELP[name]}", file=sys.stderr)
            sys.exit(1)
    return quantities


def list_quantities_and_exit():
    """Print all available quantities with descriptions and exit."""
    print("Available quantities for --quantity:")
    for name in SLICE_QUANTITY_OPTIONS:
        print(f"  {name:20s} — {QUANTITY_HELP[name]}")
    sys.exit(0)


def _load_slice_exported(export_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None, dict[str, Any]]:
    """Load exported slice data from TSV + YAML files.

    Returns (df_raw, df_weighted, df_per_omega_or_None, meta).
    """
    meta_path = export_dir / "slice_meta.yaml"
    if not meta_path.exists():
        print(f"ERROR: slice_meta.yaml not found in {export_dir}", file=sys.stderr)
        sys.exit(1)
    with open(meta_path, encoding="utf-8") as f:
        meta = yaml.safe_load(f)

    raw_path = export_dir / "slice_raw.tsv"
    if not raw_path.exists():
        print(f"ERROR: slice_raw.tsv not found in {export_dir}", file=sys.stderr)
        sys.exit(1)
    df_raw = pd.read_csv(raw_path, sep="\t")
    for col in ["omega", "x", "y", "z", "Sx_flux", "Sy_flux", "Sz_flux", "weight", "Sx_w", "Sy_w", "Sz_w"]:
        if col in df_raw.columns:
            df_raw[col] = pd.to_numeric(df_raw[col], errors="coerce")

    weighted_path = export_dir / "slice_weighted.tsv"
    if not weighted_path.exists():
        print(f"ERROR: slice_weighted.tsv not found in {export_dir}", file=sys.stderr)
        sys.exit(1)
    df_weighted = pd.read_csv(weighted_path, sep="\t")
    for col in ["x", "y", "z", "Sx", "Sy", "Sz"]:
        if col in df_weighted.columns:
            df_weighted[col] = pd.to_numeric(df_weighted[col], errors="coerce")

    df_per_omega = None
    per_omega_path = export_dir / "slice_per_omega.tsv"
    if per_omega_path.exists():
        df_per_omega = pd.read_csv(per_omega_path, sep="\t")
        for col in ["omega", "x", "y", "z", "Sx_w", "Sy_w", "Sz_w"]:
            if col in df_per_omega.columns:
                df_per_omega[col] = pd.to_numeric(df_per_omega[col], errors="coerce")

    return df_raw, df_weighted, df_per_omega, meta


def _plot_slice_from_export(
    df_raw: pd.DataFrame,
    df_weighted: pd.DataFrame,
    df_per_omega: pd.DataFrame | None,
    meta: dict[str, Any],
    quantities: list[str],
    per_omega: bool,
    check_periodicity_flag: bool,
    out_dir: Path,
) -> None:
    """Generate slice plots from exported data."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)

    slice_axis = meta["slice_axis"]
    slice_value = meta["slice_value"]
    folded_value = meta["folded_value"]
    temperature = meta["temperature"]
    tolerance = meta["tolerance"]
    period = meta["period"]
    plane_values = meta["plane_values"]

    perp_axis = "x" if slice_axis == "y" else "y"

    # --- Plot 1: Quantity vs perpendicular axis (weighted sum) ---
    print(f"Weighted sum: {len(df_weighted)} unique spatial points")
    for qty in quantities:
        vals = compute_quantity_vectorized(
            df_weighted["Sx"].values, df_weighted["Sy"].values, df_weighted["Sz"].values, qty
        )
        sort_idx = np.argsort(df_weighted[perp_axis].values)
        x_sorted = df_weighted[perp_axis].values[sort_idx]
        v_sorted = vals[sort_idx]

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(x_sorted, v_sorted, linewidth=1.0, alpha=0.8, label=qty)
        ax.set_xlabel(perp_axis)
        ax.set_ylabel(qty)
        ax.set_title(f"{qty} vs {perp_axis}  ({slice_axis}={slice_value}, T={temperature}K)")
        ax.legend()
        fig.tight_layout()

        fname = f"slice_{qty}_vs_{perp_axis}.png"
        fig.savefig(out_dir / fname, dpi=150)
        plt.close(fig)
        print(f"  Saved {fname}")

    # --- Plot 2: Per-omega slices (optional) ---
    if per_omega and df_per_omega is not None and not df_per_omega.empty:
        omega_unique = sorted(df_per_omega["omega"].dropna().unique())
        for qty in quantities:
            fig, ax = plt.subplots(figsize=(10, 5))
            for omega_val in omega_unique:
                df_w = df_per_omega[df_per_omega["omega"] == omega_val]
                if df_w.empty:
                    continue
                v = compute_quantity_vectorized(
                    df_w["Sx_w"].fillna(0).values,
                    df_w["Sy_w"].fillna(0).values,
                    df_w["Sz_w"].fillna(0).values,
                    qty,
                )
                sort_idx = np.argsort(df_w[perp_axis].values)
                ax.plot(
                    df_w[perp_axis].values[sort_idx], v[sort_idx],
                    linewidth=0.8, alpha=0.7, label=f"ω={omega_val:.3e}",
                )
            ax.set_xlabel(perp_axis)
            ax.set_ylabel(qty)
            ax.set_title(f"{qty} vs {perp_axis} per ω  ({slice_axis}={slice_value})")
            ax.legend(fontsize=7, ncol=2)
            fig.tight_layout()

            fname = f"slice_{qty}_vs_{perp_axis}_per_omega.png"
            fig.savefig(out_dir / fname, dpi=150)
            plt.close(fig)
            print(f"  Saved {fname}")

    # --- Periodicity check (optional) ---
    if check_periodicity_flag:
        for qty in quantities:
            result = check_periodicity(
                df_raw, slice_axis, folded_value, tolerance, period, qty, temperature
            )
            print(f"\n--- Periodicity check for {qty} ---")
            print(f"Period offsets: {result['period_offsets']}")
            if result["summed"]:
                print(f"Summed (T={temperature}K): values={result['summed']['quantity_values']}")
                print(f"  relative diffs: {result['summed']['relative_diffs']}")
            for omega_val, data in result.get("per_omega", {}).items():
                print(f"  ω={omega_val:.3e}: values={data['quantity_values']}")
                print(f"    rel_diffs={data['relative_diffs']}")


def cmd_slice(args: argparse.Namespace) -> None:
    """Extract 1D slices from collected SRFlux data, export, and plot."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Validate required args (relaxed for --list-quantities)
    missing = []
    if not args.dataset:
        missing.append("--dataset")
    if not args.slice_axis:
        missing.append("--slice-axis")
    if args.slice_value is None:
        missing.append("--slice-value")
    if not args.out:
        missing.append("--out")
    if missing:
        print(f"ERROR: missing required arguments: {', '.join(missing)}", file=sys.stderr)
        print("Usage: scuff_pv_post.py slice --dataset DIR --slice-axis {x,y} --slice-value VAL --out DIR", file=sys.stderr)
        sys.exit(1)

    dataset_dir = Path(args.dataset).resolve()
    out_root = Path(args.out).resolve()

    # Parse config and period
    config_name = args.config if args.config else config_from_dataset_path(dataset_dir)
    period = get_period_from_config(config_name)
    print(f"Config: {config_name}, period = {period}")

    # Parse quantities
    quantities = [q.strip() for q in args.quantity.split(",")]
    validate_quantities(quantities)

    # Parse planes
    plane_values = [float(p.strip()) for p in args.planes.split(",")]

    slice_axis = args.slice_axis
    slice_value = args.slice_value
    tolerance = args.tolerance
    temperature = args.temperature

    # Period-fold slice_value for x-axis display
    if slice_axis == "x":
        folded_value = np.mod(slice_value, period)
        print(f"x-slice: {slice_value} folded to {folded_value:.4f} (mod period {period})")
    else:
        folded_value = slice_value
        print(f"y-slice: {slice_value}")

    # Batch-aware output: {out_root}/{batch_name}/{config_name}/figures/slice/
    batch_name = args.batch_name or datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = out_root / batch_name / config_name / "figures" / "slice"
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Export step ---
    meta = slice_export_data(
        dataset_dir=dataset_dir,
        config_name=config_name,
        slice_axis=slice_axis,
        slice_value=slice_value,
        temperature=temperature,
        tolerance=tolerance,
        plane_values=plane_values,
        quantities=quantities,
        out_dir=out_dir,
    )

    if meta.get("status") == "empty":
        return

    if getattr(args, "export_only", False):
        print("Export-only mode: skipping plots.")
        return

    # --- Load exported data and plot ---
    df_raw, df_weighted, df_per_omega, meta = _load_slice_exported(out_dir)

    _plot_slice_from_export(
        df_raw=df_raw,
        df_weighted=df_weighted,
        df_per_omega=df_per_omega,
        meta=meta,
        quantities=quantities,
        per_omega=getattr(args, "per_omega", False),
        check_periodicity_flag=getattr(args, "check_periodicity", False),
        out_dir=out_dir,
    )

    # --- Plot 3: Evpoint scatter with slice line (only in full slice mode) ---
    # Load full dataset for evpoint positions
    df_full = load_srflux_with_omega(dataset_dir)
    df_full = df_full[df_full["z"].isin(plane_values)]

    for z_val in plane_values:
        df_z = df_full[df_full["z"] == z_val]
        if df_z.empty:
            continue
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.scatter(df_z["x"], df_z["y"], s=1, alpha=0.3, c="gray", label="evpoints")

        if slice_axis == "y":
            ax.axhline(y=folded_value, color="red", linewidth=1.5, linestyle="--",
                        label=f"y={folded_value}")
            for k in range(1, 3):
                ax.axhline(y=folded_value, color="red", linewidth=0.5, linestyle=":",
                           alpha=0.3)
        else:
            for k in range(3):
                xk = folded_value + k * period
                ax.axvline(x=xk, color="red", linewidth=1.5 if k == 0 else 0.8,
                           linestyle="--" if k == 0 else ":", alpha=1.0 if k == 0 else 0.4,
                           label=f"x={xk:.2f}" if k == 0 else None)
            if slice_axis == "x" and abs(slice_value - folded_value) > 0.01:
                ax.axvline(x=slice_value, color="blue", linewidth=1.0, linestyle="-.",
                           label=f"x={slice_value} (unshifted)")

        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(f"Evpoints at z={z_val}, slice: {slice_axis}={slice_value}")
        ax.legend(fontsize=8)
        ax.set_aspect("equal")
        fig.tight_layout()

        fname = f"evpoints_z{z_val}_slice_{slice_axis}.png"
        fig.savefig(out_dir / fname, dpi=150)
        plt.close(fig)
        print(f"  Saved {fname}")

    print(f"\nDone. Outputs in {out_dir}")


def cmd_slice_plot(args: argparse.Namespace) -> None:
    """Generate slice plots from previously exported data."""
    export_dir = Path(args.export_dir).resolve()
    out_root = Path(args.out).resolve()

    if not export_dir.is_dir():
        print(f"ERROR: export directory not found: {export_dir}", file=sys.stderr)
        sys.exit(1)

    df_raw, df_weighted, df_per_omega, meta = _load_slice_exported(export_dir)

    if meta.get("status") == "empty":
        print("Export data is empty, nothing to plot.")
        return

    # Use quantities from meta or override
    if args.quantity:
        quantities = [q.strip() for q in args.quantity.split(",")]
        validate_quantities(quantities)
    else:
        quantities = meta["quantities"]

    # Batch-aware output: {out_root}/{batch_name}/{config}/figures/slice/
    config_name = meta.get("config_name", "unknown")
    batch_name = args.batch_name or datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = out_root / batch_name / config_name / "figures" / "slice"
    out_dir.mkdir(parents=True, exist_ok=True)

    _plot_slice_from_export(
        df_raw=df_raw,
        df_weighted=df_weighted,
        df_per_omega=df_per_omega,
        meta=meta,
        quantities=quantities,
        per_omega=getattr(args, "per_omega", False),
        check_periodicity_flag=getattr(args, "check_periodicity", False),
        out_dir=out_dir,
    )

    print(f"\nDone. Outputs in {out_dir}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="scuff_pv_post",
        description="Post-processing tool for SCUFF-EM near-field heat transfer results.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # collect subcommand
    p_collect = subparsers.add_parser("collect", help="Merge SRFlux files into canonical TSV datasets.")
    p_collect.add_argument("--task-dir", required=True, help="Directory containing subtask-* subdirectories.")
    p_collect.add_argument("--group-name", required=True, help="Config group name (e.g., 'c1', 'ad4').")
    p_collect.add_argument("--filebase", default="task", help="Base name for SRFlux file (default: 'task').")
    p_collect.add_argument("--out", required=True, help="Output directory for canonical dataset.")
    p_collect.add_argument("--coord-tol", type=float, default=1e-9, help="Tolerance for point_id coordinate hashing.")
    p_collect.add_argument("--duplicate-policy", choices=["error", "mean"], default="error",
                           help="Policy for handling duplicates: 'error' or 'mean'.")
    p_collect.add_argument("--last", action="store_true",
                           help="Only process data after the last '# scuff-neq run on' header.")

    # integrate subcommand
    p_integrate = subparsers.add_parser("integrate", help="Run scuff-integrate on a collected dataset.")
    p_integrate.add_argument("--dataset", required=True, help="Path to canonical dataset directory (output of collect).")
    p_integrate.add_argument("--temperature-file", required=True, help="Path to TemperatureFile for scuff-integrate.")
    p_integrate.add_argument("--backend", default="scuff-integrate", help="Integration backend name (default: scuff-integrate).")
    p_integrate.add_argument("--scuff-integrate-path", default="scuff-integrate",
                             help="Path to scuff-integrate binary (default: scuff-integrate).")

    # plot subcommand
    p_plot = subparsers.add_parser("plot", help="Generate cell-map and quiver plots from a collected dataset.")
    p_plot.add_argument("--dataset", required=True, help="Path to canonical dataset directory (output of collect).")
    p_plot.add_argument("--planes", required=True, help="Comma-separated z-plane values (e.g., '0.75,1.0,1.5').")
    p_plot.add_argument("--mode", default="cellmap-quiver",
                        choices=["cellmap-quiver", "cellmap", "scatter"],
                        help="Plot mode (default: cellmap-quiver).")
    p_plot.add_argument("--color-by", default="absSxy",
                        help=f"Scalar for color. Options: {','.join(COLOR_BY_OPTIONS)} (default: absSxy).")
    p_plot.add_argument("--vector", default="Sx,Sy",
                        help="Vector components for quiver arrows (default: 'Sx,Sy').")
    p_plot.add_argument("--source-mode", default="sum", choices=["sum", "separate"],
                        help="Source aggregation mode (default: sum).")
    p_plot.add_argument("--stride", type=int, default=3,
                        help="Quiver arrow stride (default: 3).")
    p_plot.add_argument("--discard-abnormal", type=float, default=None,
                        help="Discard vectors whose log-magnitude exceeds mean + N * stderr (default: disabled).")
    p_plot.add_argument("--no-discard-abnormal", dest="discard_abnormal", action="store_const", const=None,
                        help="Disable filtering and reproduce the previous quiver behavior.")
    p_plot.add_argument("--quiver-uniform", dest="quiver_uniform", action="store_true", default=True,
                        help="Draw all quiver arrows with equal length (default).")
    p_plot.add_argument("--quiver-physical", dest="quiver_uniform", action="store_false",
                        help="Draw quiver lengths proportional to vector magnitude.")
    p_plot.add_argument("--only-central", action="store_true",
                        help="Supplement with central-period-only plots (k=1, x in [p, 2p)). Requires group_name in manifest.yaml.")
    p_plot.add_argument("--heatmap-comparison", action="store_true",
                        help="Generate side-by-side heatmap comparison: left with (log) quivers, right without quivers.")
    p_plot.add_argument("--quiver-normalized", action="store_true", default=True,
                        help="Generate additional quiver plots with uniform arrow length (default: true).")
    p_plot.add_argument("--no-quiver-normalized", dest="quiver_normalized", action="store_false",
                        help="Disable normalized quiver plots.")
    p_plot.add_argument("--quiver-dense", action="store_true", default=True,
                        help="Generate additional quiver plots with 2x density (stride/2) (default: true).")
    p_plot.add_argument("--no-quiver-dense", dest="quiver_dense", action="store_false",
                        help="Disable dense quiver plots.")
    p_plot.add_argument("--batch-name", type=str, default=None,
                        help="Batch name for output subdirectory (default: auto-timestamp).")
    p_plot.add_argument("--out", required=True, help="Output root directory for figures.")

    # run-all subcommand
    p_runall = subparsers.add_parser("run-all",
        help="Auto-discover all configs and run full pipeline (collect + plot + SIFlux).")
    p_runall.add_argument("--task-dir", required=True,
                          help="Directory containing subtask-* subdirectories.")
    p_runall.add_argument("--out-dir", required=True,
                          help="Output root directory (pv_out_{config} created inside).")
    p_runall.add_argument("--planes", default="0.75,1.0,1.5",
                          help="Comma-separated z-plane values (default: '0.75,1.0,1.5').")
    p_runall.add_argument("--color-by", default="absSxy,Sz",
                          help=f"Scalar(s) for color (default: 'absSxy,Sz'). Options: {','.join(COLOR_BY_OPTIONS)}")
    p_runall.add_argument("--stride", type=int, default=3,
                          help="Quiver arrow stride (default: 3).")
    p_runall.add_argument("--discard-abnormal", type=float, default=None,
                          help="Discard vectors whose log-magnitude exceeds mean + N * stderr (default: disabled).")
    p_runall.add_argument("--no-discard-abnormal", dest="discard_abnormal", action="store_const", const=None,
                          help="Disable filtering and reproduce the previous quiver behavior.")
    p_runall.add_argument("--quiver-uniform", dest="quiver_uniform", action="store_true", default=True,
                          help="Draw all quiver arrows with equal length (default).")
    p_runall.add_argument("--quiver-physical", dest="quiver_uniform", action="store_false",
                          help="Draw quiver lengths proportional to vector magnitude.")
    p_runall.add_argument("--quiver-normalized", action="store_true", default=True,
                          help="Generate additional quiver plots with uniform arrow length (default: true).")
    p_runall.add_argument("--no-quiver-normalized", dest="quiver_normalized", action="store_false",
                          help="Disable normalized quiver plots.")
    p_runall.add_argument("--quiver-dense", action="store_true", default=True,
                          help="Generate additional quiver plots with 2x density (stride/2) (default: true).")
    p_runall.add_argument("--no-quiver-dense", dest="quiver_dense", action="store_false",
                          help="Disable dense quiver plots.")
    p_runall.add_argument("--last", action="store_true",
                          help="Only process data after the last '# scuff-neq run on' header.")
    p_runall.add_argument("--only-central", action="store_true",
                          help="Supplement SRFlux plots with central-period-only versions (x in [p, 2p)).")
    p_runall.add_argument("--heatmap-comparison", action="store_true",
                          help="Generate side-by-side heatmap comparison (quiver | no quiver) for each z-plane.")

    # slice subcommand
    p_slice = subparsers.add_parser("slice",
        help="Extract 1D slices and plot quantities along x or y axis.")
    p_slice.add_argument("--dataset",
                         help="Path to collected dataset directory (contains merged_srflux.tsv).")
    p_slice.add_argument("--config", default=None,
                         help="Config name for period lookup (e.g., 'ad1'). Auto-derived from dataset path if omitted.")
    p_slice.add_argument("--slice-axis", choices=["x", "y"],
                         help="Axis to slice along: 'x' or 'y'.")
    p_slice.add_argument("--slice-value", type=float,
                         help="Coordinate value for the slice (auto-folded by period for x-axis).")
    p_slice.add_argument("--temperature", type=float, default=300.0,
                         help="Temperature in K for Bose-Einstein weighting (default: 300).")
    p_slice.add_argument("--tolerance", type=float, default=0.05,
                         help="Matching tolerance for slice point selection (default: 0.05).")
    p_slice.add_argument("--quantity", default="Sx,Sy,absSxy",
                         help=f"Comma-separated quantities to plot (default: 'Sx,Sy,absSxy'). "
                              f"Use --list-quantities to see all options.")
    p_slice.add_argument("--planes", default="0.75",
                         help="Comma-separated z-plane values (default: '0.75').")
    p_slice.add_argument("--per-omega", action="store_true",
                         help="Also plot per-omega slices (not just Bose-Einstein weighted sum).")
    p_slice.add_argument("--check-periodicity", action="store_true",
                         help="Run periodicity consistency check across the 3 simulated periods.")
    p_slice.add_argument("--list-quantities", action="store_true",
                         help="Print all available quantities and exit.")
    p_slice.add_argument("--export-only", action="store_true",
                         help="Export data to TSV only, skip plot generation.")
    p_slice.add_argument("--batch-name", type=str, default=None,
                         help="Batch name for output subdirectory (default: auto-timestamp).")
    p_slice.add_argument("--out",
                         help="Output root directory for exported data and figures.")

    # slice-plot subcommand
    p_slice_plot = subparsers.add_parser("slice-plot",
        help="Generate slice plots from previously exported data.")
    p_slice_plot.add_argument("--export-dir", required=True,
                              help="Directory containing slice export files (slice_raw.tsv, etc.).")
    p_slice_plot.add_argument("--batch-name", type=str, default=None,
                              help="Batch name for output subdirectory (default: auto-timestamp).")
    p_slice_plot.add_argument("--out", required=True,
                              help="Output root directory for figures.")
    p_slice_plot.add_argument("--quantity", default=None,
                              help="Override quantities to plot (default: use values from slice_meta.yaml).")
    p_slice_plot.add_argument("--per-omega", action="store_true",
                              help="Also plot per-omega slices.")
    p_slice_plot.add_argument("--check-periodicity", action="store_true",
                              help="Run periodicity consistency check.")

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    # Handle --list-quantities before anything else
    if getattr(args, "list_quantities", False):
        list_quantities_and_exit()

    if args.command == "collect":
        cmd_collect(args)
    elif args.command == "integrate":
        cmd_integrate(args)
    elif args.command == "plot":
        cmd_plot(args)
    elif args.command == "run-all":
        cmd_run_all(args)
    elif args.command == "slice":
        cmd_slice(args)
    elif args.command == "slice-plot":
        cmd_slice_plot(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
