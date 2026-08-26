"""Generate Python-integrated Poynting flux diagrams from collected SRFlux data.

This script expects a dataset directory produced by:

    python scuff_pv_post.py collect ...

It numerically integrates Sx/Sy/Sz over omega using numpy.trapezoid, then writes
cellmap + fixed-scale quiver plots. All-NaN point histories remain NaN so masked
regions stay blank in the figures.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


FLUX_COLS = ["Sx_flux", "Sy_flux", "Sz_flux"]
SIGNED_COLOR_BY = {"Sx", "Sy", "Sz"}
SPARSE_OCCUPANCY_THRESHOLD = 0.55
ADAPTIVE_CELL_DISTANCE_THRESHOLD = 0.1
ADAPTIVE_CELL_HOLE_GAP_FACTOR = 2.25
ADAPTIVE_CELL_GAP_PAD = 1.05
FIGURE_DPI = 340
LABEL_FONT_SIZE = 18
TICK_FONT_SIZE = 16
SUPPORTED_FIGURE_FORMATS = {"png", "pdf", "svg", "eps"}


def apply_matplotlib_theme() -> None:
    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "mathtext.fallback": "stix",
            "axes.labelsize": LABEL_FONT_SIZE,
            "axes.titlesize": LABEL_FONT_SIZE,
            "xtick.labelsize": TICK_FONT_SIZE,
            "ytick.labelsize": TICK_FONT_SIZE,
        }
    )


def add_matched_colorbar(fig: Any, ax: Any, artist: Any, label: str) -> Any:
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.08)
    cbar = fig.colorbar(artist, cax=cax)
    cbar.set_label(label, fontsize=LABEL_FONT_SIZE)
    cbar.ax.tick_params(labelsize=TICK_FONT_SIZE)
    return cbar


@dataclass(frozen=True)
class CellmapRender:
    artist: Any
    method: str
    grid_cells: int
    occupancy: float
    color_scale: str
    x_threshold: float | None = None
    y_threshold: float | None = None


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


def _threshold_bounds_1d(coords: np.ndarray, threshold: float) -> tuple[np.ndarray, np.ndarray]:
    coords = np.asarray(coords, dtype=float)
    order = np.argsort(coords, kind="mergesort")
    sorted_coords = coords[order]
    left = sorted_coords - threshold / 2
    right = sorted_coords + threshold / 2

    if len(sorted_coords) > 1:
        gaps = np.diff(sorted_coords)
        close_idx = np.flatnonzero(gaps <= threshold)
        mids = (sorted_coords[:-1] + sorted_coords[1:]) / 2
        right[close_idx] = mids[close_idx]
        left[close_idx + 1] = mids[close_idx]

    left_out = np.empty_like(left)
    right_out = np.empty_like(right)
    left_out[order] = left
    right_out[order] = right
    return left_out, right_out


def adaptive_axis_threshold(
    gaps_or_coords: np.ndarray,
    base_threshold: float = ADAPTIVE_CELL_DISTANCE_THRESHOLD,
    hole_gap_factor: float = ADAPTIVE_CELL_HOLE_GAP_FACTOR,
    gap_pad: float = ADAPTIVE_CELL_GAP_PAD,
) -> float:
    """Choose a near-neighbor threshold without crossing large blank gaps."""
    gaps = np.asarray(gaps_or_coords, dtype=float)
    gaps = gaps[np.isfinite(gaps) & (gaps > 0)]
    if gaps.size == 0:
        return float(base_threshold)

    typical_gap = float(np.nanpercentile(gaps, 25))
    hole_gap_floor = max(float(base_threshold) * hole_gap_factor, typical_gap * hole_gap_factor)
    near_gaps = gaps[gaps <= hole_gap_floor + 1e-12]
    if near_gaps.size == 0:
        near_gaps = np.array([float(np.min(gaps))])
    return float(max(base_threshold, float(np.max(near_gaps)) * gap_pad))


def adaptive_cell_thresholds(
    x: np.ndarray,
    y: np.ndarray,
    base_threshold: float = ADAPTIVE_CELL_DISTANCE_THRESHOLD,
) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x_gaps: list[float] = []
    y_gaps: list[float] = []

    for y_val in np.unique(y[np.isfinite(y)]):
        xs = np.array(sorted(np.unique(x[y == y_val])), dtype=float)
        if len(xs) > 1:
            x_gaps.extend(np.diff(xs).tolist())

    for x_val in np.unique(x[np.isfinite(x)]):
        ys = np.array(sorted(np.unique(y[x == x_val])), dtype=float)
        if len(ys) > 1:
            y_gaps.extend(np.diff(ys).tolist())

    return (
        adaptive_axis_threshold(np.asarray(x_gaps), base_threshold=base_threshold),
        adaptive_axis_threshold(np.asarray(y_gaps), base_threshold=base_threshold),
    )


def threshold_rectangles_from_points(
    x: np.ndarray,
    y: np.ndarray,
    x_threshold: float = ADAPTIVE_CELL_DISTANCE_THRESHOLD,
    y_threshold: float | None = None,
) -> list[np.ndarray]:
    """Build axis-aligned cells, bridging only near same-row/column neighbors."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    y_threshold = x_threshold if y_threshold is None else y_threshold
    n = len(x)
    x_left = x - x_threshold / 2
    x_right = x + x_threshold / 2
    y_bottom = y - y_threshold / 2
    y_top = y + y_threshold / 2

    for y_val in np.unique(y):
        idx = np.flatnonzero(y == y_val)
        x_left[idx], x_right[idx] = _threshold_bounds_1d(x[idx], x_threshold)

    for x_val in np.unique(x):
        idx = np.flatnonzero(x == x_val)
        y_bottom[idx], y_top[idx] = _threshold_bounds_1d(y[idx], y_threshold)

    polygons: list[np.ndarray] = []
    for i in range(n):
        polygons.append(
            np.array(
                [
                    [x_left[i], y_bottom[i]],
                    [x_right[i], y_bottom[i]],
                    [x_right[i], y_top[i]],
                    [x_left[i], y_top[i]],
                ],
                dtype=float,
            )
        )
    return polygons


def plot_threshold_cellmap(
    ax: Any,
    x: np.ndarray,
    y: np.ndarray,
    values: np.ndarray,
    cmap: str,
    norm_kw: dict[str, Any],
    threshold: float = ADAPTIVE_CELL_DISTANCE_THRESHOLD,
    x_threshold: float | None = None,
    y_threshold: float | None = None,
) -> Any:
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize

    finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(values)
    x_f, y_f, v_f = x[finite], y[finite], values[finite]
    x_threshold = threshold if x_threshold is None else x_threshold
    y_threshold = threshold if y_threshold is None else y_threshold
    polygons = threshold_rectangles_from_points(x_f, y_f, x_threshold, y_threshold)

    norm = norm_kw.get("norm")
    clim_kw = {k: v for k, v in norm_kw.items() if k in ("vmin", "vmax")}
    if norm is None and clim_kw:
        norm = Normalize(**clim_kw)
    collection = PolyCollection(polygons, array=v_f, cmap=cmap, norm=norm, edgecolors="none", linewidths=0)
    if norm is None and clim_kw:
        collection.set_clim(clim_kw.get("vmin"), clim_kw.get("vmax"))
    ax.add_collection(collection)
    if polygons:
        all_vertices = np.vstack(polygons)
        ax.set_xlim(float(np.min(all_vertices[:, 0])), float(np.max(all_vertices[:, 0])))
        ax.set_ylim(float(np.min(all_vertices[:, 1])), float(np.max(all_vertices[:, 1])))
    return collection


def read_group_name(dataset_dir: Path) -> str:
    manifest_path = dataset_dir / "manifest.yaml"
    if not manifest_path.exists():
        return dataset_dir.parent.name
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    return str(manifest.get("group_name") or dataset_dir.parent.name)


def integrate_values(omega: np.ndarray, values: np.ndarray) -> float:
    valid = np.isfinite(omega) & np.isfinite(values)
    if not np.any(valid):
        return np.nan
    omega_v = omega[valid]
    values_v = values[valid]
    if len(omega_v) == 1:
        return float(values_v[0])
    order = np.argsort(omega_v)
    return float(np.trapezoid(values_v[order], omega_v[order]))


def integrate_dataset(dataset_dir: Path) -> pd.DataFrame:
    srflux_path = dataset_dir / "merged_srflux.tsv"
    if not srflux_path.exists():
        raise FileNotFoundError(f"missing collected data: {srflux_path}")

    df = pd.read_csv(srflux_path, sep="\t")
    numeric_cols = ["omega", "x", "y", "z", *FLUX_COLS]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    group_cols = ["transform_id", "source_id", "point_id", "x", "y", "z"]
    rows: list[dict[str, Any]] = []
    for key, group in df.groupby(group_cols, sort=False, dropna=False):
        row = dict(zip(group_cols, key))
        omega = group["omega"].to_numpy(dtype=float)
        for src_col, dst_col in zip(FLUX_COLS, ["Sx", "Sy", "Sz"]):
            row[dst_col] = integrate_values(omega, group[src_col].to_numpy(dtype=float))
        rows.append(row)

    return pd.DataFrame(rows)


def color_values(df: pd.DataFrame, color_by: str) -> np.ndarray:
    sx = df["Sx"].to_numpy(dtype=float)
    sy = df["Sy"].to_numpy(dtype=float)
    sz = df["Sz"].to_numpy(dtype=float)
    if color_by == "Sx":
        return sx
    if color_by == "Sy":
        return sy
    if color_by == "Sz":
        return sz
    if color_by == "absS":
        return np.sqrt(sx**2 + sy**2 + sz**2)
    return np.sqrt(sx**2 + sy**2)


def grid_plane(plane_df: pd.DataFrame, values: np.ndarray) -> dict[str, Any]:
    x_unique = np.array(sorted(plane_df["x"].unique()), dtype=float)
    y_unique = np.array(sorted(plane_df["y"].unique()), dtype=float)
    x_index = {v: i for i, v in enumerate(x_unique)}
    y_index = {v: i for i, v in enumerate(y_unique)}

    shape = (len(y_unique), len(x_unique))
    value_grid = np.full(shape, np.nan)
    sx_grid = np.full(shape, np.nan)
    sy_grid = np.full(shape, np.nan)

    for value, row in zip(values, plane_df.itertuples(index=False)):
        ix = x_index[getattr(row, "x")]
        iy = y_index[getattr(row, "y")]
        value_grid[iy, ix] = value
        sx_grid[iy, ix] = getattr(row, "Sx")
        sy_grid[iy, ix] = getattr(row, "Sy")

    dx = np.min(np.diff(x_unique)) if len(x_unique) > 1 else 1.0
    dy = np.min(np.diff(y_unique)) if len(y_unique) > 1 else 1.0
    x_edges = axis_cell_edges(x_unique)
    y_edges = axis_cell_edges(y_unique)

    return {
        "x_unique": x_unique,
        "y_unique": y_unique,
        "x_edges": x_edges,
        "y_edges": y_edges,
        "value_grid": value_grid,
        "sx_grid": sx_grid,
        "sy_grid": sy_grid,
        "cell_size": float(min(abs(dx), abs(dy))),
    }


def derive_fixed_quiver_scale(df: pd.DataFrame, target_cell_fraction: float = 0.35) -> float:
    mag = np.hypot(df["Sx"].to_numpy(dtype=float), df["Sy"].to_numpy(dtype=float))
    mag = mag[np.isfinite(mag) & (mag > 0)]
    if mag.size == 0:
        return 1.0

    xs = np.array(sorted(df["x"].dropna().unique()), dtype=float)
    ys = np.array(sorted(df["y"].dropna().unique()), dtype=float)
    dx = np.median(np.diff(xs)) if len(xs) > 1 else 1.0
    dy = np.median(np.diff(ys)) if len(ys) > 1 else 1.0
    cell_size = float(min(abs(dx), abs(dy)))
    p95 = float(np.nanpercentile(mag, 95))
    return p95 / max(target_cell_fraction * cell_size, 1e-30)


def figure_size(x_vals: np.ndarray, y_vals: np.ndarray) -> tuple[float, float]:
    x_range = max(float(np.ptp(x_vals)), 1.0)
    y_range = max(float(np.ptp(y_vals)), 1.0)
    width = min(max(7.0 * x_range / y_range, 6.0), 15.0)
    return width, 7.8


def color_norm_kwargs(values: np.ndarray, color_by: str, color_scale: str = "log") -> dict[str, Any]:
    if color_scale == "log":
        from matplotlib.colors import LogNorm, SymLogNorm

        finite = values[np.isfinite(values)]
        if finite.size == 0:
            return {}

        if color_by in SIGNED_COLOR_BY:
            abs_finite = np.abs(finite)
            positive = abs_finite[abs_finite > 0]
            if positive.size == 0:
                return {}
            vmax = float(np.nanpercentile(abs_finite, 99))
            linthresh = float(max(np.nanpercentile(positive, 5), np.nanmin(positive)))
            if vmax <= linthresh:
                vmax = linthresh * 10
            return {"norm": SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax, base=10)}

        positive = finite[finite > 0]
        if positive.size == 0:
            return {}
        vmin = float(max(np.nanpercentile(positive, 1), np.nanmin(positive)))
        vmax = float(np.nanpercentile(positive, 99))
        if vmax <= vmin:
            vmax = vmin * 10
        return {"norm": LogNorm(vmin=vmin, vmax=vmax)}

    if color_by not in SIGNED_COLOR_BY:
        return {}

    finite = values[np.isfinite(values)]
    vmax = float(np.nanpercentile(np.abs(finite), 99)) if finite.size else 1.0
    if vmax == 0:
        vmax = 1.0
    return {"vmin": -vmax, "vmax": vmax}


def draw_coloring_layer(
    ax: Any,
    plane_df: pd.DataFrame,
    values: np.ndarray,
    grid: dict[str, Any],
    color_by: str,
    cmap: str = "rainbow",
    color_scale: str = "log",
) -> CellmapRender:
    value_grid = grid["value_grid"]
    grid_cells = len(grid["x_unique"]) * len(grid["y_unique"])
    occupancy = len(plane_df) / grid_cells if grid_cells else 0.0
    use_threshold_cells = occupancy < SPARSE_OCCUPANCY_THRESHOLD
    norm_kw = color_norm_kwargs(values, color_by, color_scale=color_scale)

    if use_threshold_cells:
        x_threshold, y_threshold = adaptive_cell_thresholds(
            plane_df["x"].to_numpy(dtype=float),
            plane_df["y"].to_numpy(dtype=float),
        )
        artist = plot_threshold_cellmap(
            ax,
            plane_df["x"].to_numpy(dtype=float),
            plane_df["y"].to_numpy(dtype=float),
            values,
            cmap=cmap,
            norm_kw=norm_kw,
            x_threshold=x_threshold,
            y_threshold=y_threshold,
        )
        return CellmapRender(
            artist=artist,
            method="threshold_rectangles",
            grid_cells=int(grid_cells),
            occupancy=float(occupancy),
            color_scale=color_scale,
            x_threshold=x_threshold,
            y_threshold=y_threshold,
        )

    artist = ax.pcolormesh(
        grid["x_edges"],
        grid["y_edges"],
        value_grid,
        cmap=cmap,
        shading="flat",
        **norm_kw,
    )
    return CellmapRender(
        artist=artist,
        method="rectilinear_pcolormesh",
        grid_cells=int(grid_cells),
        occupancy=float(occupancy),
        color_scale=color_scale,
    )


def draw_quiver_layer(
    ax: Any,
    plane_df: pd.DataFrame,
    values: np.ndarray,
    grid: dict[str, Any],
    cellmap: CellmapRender,
    stride: int,
    quiver_scale: float,
) -> bool:
    stride = max(int(stride), 1)
    if cellmap.method == "threshold_rectangles":
        qx = plane_df["x"].to_numpy(dtype=float)[::stride]
        qy = plane_df["y"].to_numpy(dtype=float)[::stride]
        sx = plane_df["Sx"].to_numpy(dtype=float)[::stride]
        sy = plane_df["Sy"].to_numpy(dtype=float)[::stride]
        qv = values[::stride]
        mask = np.isfinite(qx) & np.isfinite(qy) & np.isfinite(sx) & np.isfinite(sy) & np.isfinite(qv)
    else:
        xg, yg = np.meshgrid(grid["x_unique"], grid["y_unique"])
        skip = (slice(None, None, stride), slice(None, None, stride))
        qx = xg[skip]
        qy = yg[skip]
        sx = grid["sx_grid"][skip]
        sy = grid["sy_grid"][skip]
        mask = np.isfinite(sx) & np.isfinite(sy) & np.isfinite(grid["value_grid"][skip])

    if not np.any(mask):
        return False

    ax.quiver(
        qx[mask],
        qy[mask],
        sx[mask],
        sy[mask],
        angles="xy",
        scale_units="xy",
        scale=quiver_scale,
        color="k",
        width=0.003,
        alpha=0.75,
    )
    return True


def parse_figure_formats(raw: str) -> list[str]:
    formats = [item.strip().lower().lstrip(".") for item in raw.split(",") if item.strip()]
    if not formats:
        raise ValueError("at least one figure format is required")
    unsupported = [item for item in formats if item not in SUPPORTED_FIGURE_FORMATS]
    if unsupported:
        raise ValueError(f"unsupported figure format(s): {', '.join(unsupported)}")
    return list(dict.fromkeys(formats))


def plot_integrated(
    integrated: pd.DataFrame,
    out_dir: Path,
    group_name: str,
    batch_name: str,
    planes: list[float],
    color_by_list: list[str],
    stride: int,
    quiver_scale: float,
    include_quiver: bool = True,
    color_scale: str = "log",
    colormap: str = "rainbow",
    figure_format: str = "png",
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    apply_matplotlib_theme()
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    figure_formats = parse_figure_formats(figure_format)
    manifest_entries: list[dict[str, Any]] = []
    plane_tol = 1e-6

    for plane_z in planes:
        plane_df = integrated[np.abs(integrated["z"] - plane_z) < plane_tol].copy()
        if plane_df.empty:
            continue

        for color_by in color_by_list:
            values = color_values(plane_df, color_by)
            grid = grid_plane(plane_df, values)

            fig, ax = plt.subplots(figsize=figure_size(grid["x_unique"], grid["y_unique"]))
            cellmap = draw_coloring_layer(
                ax,
                plane_df,
                values,
                grid,
                color_by,
                cmap=colormap,
                color_scale=color_scale,
            )
            add_matched_colorbar(fig, ax, cellmap.artist, label=f"integrated {color_by}")

            quiver_drawn = False
            if include_quiver:
                quiver_drawn = draw_quiver_layer(ax, plane_df, values, grid, cellmap, stride, quiver_scale)

            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.tick_params(labelsize=TICK_FONT_SIZE)
            ax.set_aspect("equal")
            layer_label = " + fixed-scale quiver" if include_quiver else ""
            scale_label = "" if color_scale == "linear" else f" ({color_scale} color)"
            ax.set_title(f"{group_name} z={plane_z:.3f} integrated {color_by}{scale_label}{layer_label}")

            suffix = "quiver1x_fixedscale" if include_quiver else "coloronly"
            scale_part = "" if color_scale == "linear" else f"_{color_scale}"
            figure_files = {}
            for item_format in figure_formats:
                fname = f"{batch_name}_{group_name}_z{plane_z:.3f}_{color_by}_integrated_{colormap}{scale_part}_{suffix}.{item_format}"
                save_kwargs: dict[str, Any] = {"bbox_inches": "tight"}
                if item_format == "png":
                    save_kwargs["dpi"] = FIGURE_DPI
                fig.savefig(out_dir / fname, **save_kwargs)
                figure_files[item_format] = fname
            plt.close(fig)

            manifest_entries.append(
                {
                    "filename": figure_files[figure_formats[0]],
                    "figure_formats": figure_formats,
                    "figures": figure_files,
                    "plane_z": plane_z,
                    "color_by": color_by,
                    "colormap": colormap,
                    "color_scale": color_scale,
                    "quiver_stride": stride,
                    "quiver_scale": quiver_scale,
                    "include_quiver": include_quiver,
                    "quiver_drawn": quiver_drawn,
                    "n_points": int(len(plane_df)),
                    "n_finite": int(np.isfinite(values).sum()),
                    "cellmap_method": cellmap.method,
                    "adaptive_cell_base_distance_threshold": ADAPTIVE_CELL_DISTANCE_THRESHOLD if cellmap.x_threshold else None,
                    "adaptive_cell_x_distance_threshold": cellmap.x_threshold,
                    "adaptive_cell_y_distance_threshold": cellmap.y_threshold,
                    "grid_cells": cellmap.grid_cells,
                    "grid_occupancy": cellmap.occupancy,
                }
            )

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "group_name": group_name,
        "batch_name": batch_name,
        "integration": "numpy.trapezoid over omega",
        "colormap": colormap,
        "color_scale": color_scale,
        "quiver_stride": stride,
        "quiver_scale": quiver_scale,
        "include_quiver": include_quiver,
        "figure_format": figure_formats[0],
        "figure_formats": figure_formats,
        "vector_export": any(item in {"pdf", "svg", "eps"} for item in figure_formats),
        "plots": manifest_entries,
    }
    (out_dir / "integrated_flux_manifest.yaml").write_text(
        yaml.dump(manifest, default_flow_style=False, sort_keys=False),
        encoding="utf-8",
    )


def run_integrated_flux(
    dataset_dir: Path,
    out_dir: Path,
    batch_name: str | None = None,
    planes: list[float] | None = None,
    color_by_list: list[str] | None = None,
    stride: int = 1,
    quiver_scale: float | None = None,
    include_quiver: bool = True,
    color_scale: str = "log",
    colormap: str = "rainbow",
    figure_format: str = "png",
) -> dict[str, Any]:
    """Integrate one collected dataset and write TSV + figures.

    This is the module entry point used by batch_postprocess.py. The CLI below
    is intentionally thin so ordered workflows can call this function directly.
    """
    dataset_dir = Path(dataset_dir).resolve()
    out_dir = Path(out_dir).resolve()
    batch_name = batch_name or datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    group_name = read_group_name(dataset_dir)
    planes = planes or [0.75, 1.0, 1.5]
    color_by_list = color_by_list or ["absSxy", "Sz"]

    integrated = integrate_dataset(dataset_dir)
    fixed_scale = float(quiver_scale) if quiver_scale is not None else derive_fixed_quiver_scale(integrated)

    out_dir.mkdir(parents=True, exist_ok=True)
    integrated_path = out_dir / f"{group_name}_poynting_integrated_python.tsv"
    integrated.to_csv(integrated_path, sep="\t", index=False)
    plot_integrated(
        integrated,
        out_dir,
        group_name,
        batch_name,
        planes,
        color_by_list,
        stride,
        fixed_scale,
        include_quiver=include_quiver,
        color_scale=color_scale,
        colormap=colormap,
        figure_format=figure_format,
    )

    figure_paths = []
    for item_format in parse_figure_formats(figure_format):
        figure_paths.extend(sorted(out_dir.glob(f"{batch_name}_{group_name}_*_integrated_{colormap}_*.{item_format}")))
    return {
        "group_name": group_name,
        "integrated_path": integrated_path,
        "out_dir": out_dir,
        "n_figures": len(figure_paths),
        "quiver_scale": fixed_scale,
        "include_quiver": include_quiver,
        "color_scale": color_scale,
        "colormap": colormap,
        "figure_formats": parse_figure_formats(figure_format),
        "figures": figure_paths,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True, help="Collected dataset directory containing merged_srflux.tsv.")
    parser.add_argument("--out", required=True, help="Output directory for integrated TSV and figures.")
    parser.add_argument("--batch-name", default=None, help="Output filename prefix. Defaults to current UTC timestamp.")
    parser.add_argument("--planes", default="0.75,1.0,1.5", help="Comma-separated z planes.")
    parser.add_argument("--color-by", default="absSxy,Sz", help="Comma-separated: absSxy,absS,Sx,Sy,Sz.")
    parser.add_argument("--stride", type=int, default=1, help="Quiver stride. Use 1 for x1 density.")
    parser.add_argument("--quiver-scale", type=float, default=None, help="Fixed matplotlib quiver scale. Auto-derived if omitted.")
    parser.add_argument("--no-quiver", dest="include_quiver", action="store_false", help="Draw coloring only.")
    parser.add_argument("--color-scale", choices=["linear", "log"], default="log", help="Color normalization.")
    parser.add_argument("--colormap", default="rainbow", help="Matplotlib colormap name.")
    parser.add_argument("--figure-format", default="png", help="Figure format(s), comma-separated. Default: png.")
    parser.set_defaults(include_quiver=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    dataset_dir = Path(args.dataset).resolve()
    out_dir = Path(args.out).resolve()
    batch_name = args.batch_name or datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    planes = [float(p.strip()) for p in args.planes.split(",") if p.strip()]
    color_by_list = [c.strip() for c in args.color_by.split(",") if c.strip()]

    result = run_integrated_flux(
        dataset_dir=dataset_dir,
        out_dir=out_dir,
        batch_name=batch_name,
        planes=planes,
        color_by_list=color_by_list,
        stride=args.stride,
        quiver_scale=args.quiver_scale,
        include_quiver=args.include_quiver,
        color_scale=args.color_scale,
        colormap=args.colormap,
        figure_format=args.figure_format,
    )

    print(f"Wrote integrated data: {result['integrated_path']}")
    print(f"Wrote {result['n_figures']} figures to: {result['out_dir']}")
    print(f"Fixed quiver scale: {result['quiver_scale']:.10e}")


if __name__ == "__main__":
    main()
