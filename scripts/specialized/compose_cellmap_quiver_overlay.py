"""Compose cellmap and fixed-length quiver as separate image layers."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from PIL import Image

from plot_integrated_flux import (
    TICK_FONT_SIZE,
    add_matched_colorbar,
    apply_matplotlib_theme,
    color_values,
    draw_coloring_layer,
    figure_size,
    grid_plane,
)


DEFAULT_DPI = 340
DEFAULT_ARROW_LENGTH_FACTOR = 0.45


def min_positive_spacing(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    unique = np.array(sorted(np.unique(values[np.isfinite(values)])), dtype=float)
    diffs = np.diff(unique)
    diffs = diffs[diffs > 1e-12]
    if diffs.size == 0:
        return 1.0
    return float(np.min(diffs))


def fixed_quiver_components(sx: np.ndarray, sy: np.ndarray, arrow_length: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mag = np.hypot(sx, sy)
    mask = np.isfinite(sx) & np.isfinite(sy) & np.isfinite(mag) & (mag > 0)
    u = np.full_like(sx, np.nan, dtype=float)
    v = np.full_like(sy, np.nan, dtype=float)
    u[mask] = sx[mask] / mag[mask] * arrow_length
    v[mask] = sy[mask] / mag[mask] * arrow_length
    return u, v, mask


def save_layer_composite(base_path: Path, quiver_path: Path, out_path: Path, opacity: float) -> None:
    base = Image.open(base_path).convert("RGBA")
    quiver = Image.open(quiver_path).convert("RGBA")
    if quiver.size != base.size:
        quiver = quiver.resize(base.size, Image.Resampling.LANCZOS)
    if opacity < 1:
        alpha = quiver.getchannel("A")
        alpha = alpha.point(lambda p: int(p * opacity))
        quiver.putalpha(alpha)
    composed = Image.alpha_composite(base, quiver)
    composed.save(out_path)


def render_overlay(
    plane_df: pd.DataFrame,
    out_dir: Path,
    group_name: str,
    batch_name: str,
    plane_z: float,
    color_by: str,
    colormap: str,
    color_scale: str,
    arrow_length_factor: float,
    dpi: int,
    quiver_opacity: float,
) -> dict[str, Any]:
    import matplotlib

    matplotlib.use("Agg")
    apply_matplotlib_theme()
    import matplotlib.pyplot as plt

    values = color_values(plane_df, color_by)
    grid = grid_plane(plane_df, values)
    x = plane_df["x"].to_numpy(dtype=float)
    y = plane_df["y"].to_numpy(dtype=float)
    sx = plane_df["Sx"].to_numpy(dtype=float)
    sy = plane_df["Sy"].to_numpy(dtype=float)
    pixel_size = float(min(min_positive_spacing(x), min_positive_spacing(y)))
    arrow_length = pixel_size * arrow_length_factor
    u, v, quiver_mask = fixed_quiver_components(sx, sy, arrow_length)
    quiver_mask &= np.isfinite(x) & np.isfinite(y) & np.isfinite(values)

    layer_dir = out_dir / "_layers" / group_name
    layer_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    base_name = f"{batch_name}_{group_name}_z{plane_z:.3f}_{color_by}_cellmap_layer.png"
    quiver_name = f"{batch_name}_{group_name}_z{plane_z:.3f}_Sxy_quiver_layer.png"
    overlay_name = f"{batch_name}_{group_name}_z{plane_z:.3f}_{color_by}_cellmap_plus_fixedquiver_overlay.png"
    base_path = layer_dir / base_name
    quiver_path = layer_dir / quiver_name
    overlay_path = out_dir / overlay_name

    fig, ax = plt.subplots(figsize=figure_size(grid["x_unique"], grid["y_unique"]))
    cellmap = draw_coloring_layer(ax, plane_df, values, grid, color_by, cmap=colormap, color_scale=color_scale)
    add_matched_colorbar(fig, ax, cellmap.artist, label=f"integrated {color_by}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.tick_params(labelsize=TICK_FONT_SIZE)
    ax.set_aspect("equal")
    scale_label = "" if color_scale == "linear" else f" ({color_scale} color)"
    ax.set_title(f"{group_name} z={plane_z:.3f} integrated {color_by}{scale_label}")
    fig.canvas.draw()
    ax_position = ax.get_position()
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    figsize = fig.get_size_inches()
    fig.savefig(base_path, dpi=dpi)
    plt.close(fig)

    fig_q = plt.figure(figsize=figsize)
    fig_q.patch.set_alpha(0)
    ax_q = fig_q.add_axes(ax_position)
    ax_q.patch.set_alpha(0)
    if np.any(quiver_mask):
        ax_q.quiver(
            x[quiver_mask],
            y[quiver_mask],
            u[quiver_mask],
            v[quiver_mask],
            angles="xy",
            scale_units="xy",
            scale=1,
            color="k",
            width=0.0024,
            headwidth=3.4,
            headlength=4.5,
            headaxislength=3.9,
            pivot="middle",
            alpha=1.0,
        )
    ax_q.set_xlim(xlim)
    ax_q.set_ylim(ylim)
    ax_q.set_aspect("equal")
    ax_q.axis("off")
    fig_q.savefig(quiver_path, dpi=dpi, transparent=True)
    plt.close(fig_q)

    save_layer_composite(base_path, quiver_path, overlay_path, quiver_opacity)
    return {
        "filename": overlay_name,
        "cellmap_layer": str(base_path.relative_to(out_dir)),
        "quiver_layer": str(quiver_path.relative_to(out_dir)),
        "plane_z": plane_z,
        "color_by": color_by,
        "n_evpoints": int(len(plane_df)),
        "n_quivers": int(np.count_nonzero(quiver_mask)),
        "pixel_size": pixel_size,
        "arrow_length": arrow_length,
        "arrow_length_factor": arrow_length_factor,
        "quiver_opacity": quiver_opacity,
        "cellmap_method": cellmap.method,
    }


def render_group(
    integrated_path: Path,
    out_dir: Path,
    group_name: str,
    batch_name: str,
    planes: list[float],
    color_by_list: list[str],
    colormap: str,
    color_scale: str,
    arrow_length_factor: float,
    dpi: int,
    quiver_opacity: float,
) -> dict[str, Any]:
    integrated = pd.read_csv(integrated_path, sep="\t")
    for col in ["x", "y", "z", "Sx", "Sy", "Sz"]:
        integrated[col] = pd.to_numeric(integrated[col], errors="coerce")

    group_out = out_dir / group_name
    group_out.mkdir(parents=True, exist_ok=True)
    plots: list[dict[str, Any]] = []
    plane_tol = 1e-6
    for plane_z in planes:
        plane_df = integrated[np.abs(integrated["z"] - plane_z) < plane_tol].copy()
        if plane_df.empty:
            continue
        for color_by in color_by_list:
            plots.append(
                render_overlay(
                    plane_df=plane_df,
                    out_dir=group_out,
                    group_name=group_name,
                    batch_name=batch_name,
                    plane_z=plane_z,
                    color_by=color_by,
                    colormap=colormap,
                    color_scale=color_scale,
                    arrow_length_factor=arrow_length_factor,
                    dpi=dpi,
                    quiver_opacity=quiver_opacity,
                )
            )

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "group_name": group_name,
        "integrated_path": str(integrated_path),
        "render_method": "separate_cellmap_and_quiver_layers_alpha_composite",
        "drawn_together": False,
        "plots": plots,
    }
    (group_out / "cellmap_quiver_overlay_manifest.yaml").write_text(
        yaml.dump(manifest, default_flow_style=False, sort_keys=False),
        encoding="utf-8",
    )
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--configs", default="cd1,cd4")
    parser.add_argument("--planes", default="0.75,1.0,1.5")
    parser.add_argument("--color-by", default="absSxy,Sz")
    parser.add_argument("--batch-name", default=None)
    parser.add_argument("--colormap", default="rainbow")
    parser.add_argument("--color-scale", choices=["linear", "log"], default="log")
    parser.add_argument("--arrow-length-factor", type=float, default=DEFAULT_ARROW_LENGTH_FACTOR)
    parser.add_argument("--quiver-opacity", type=float, default=0.9)
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    batch_dir = Path(args.batch_dir).resolve()
    out_dir = Path(args.out).resolve()
    configs = [item.strip() for item in args.configs.split(",") if item.strip()]
    planes = [float(item.strip()) for item in args.planes.split(",") if item.strip()]
    color_by_list = [item.strip() for item in args.color_by.split(",") if item.strip()]
    batch_name = args.batch_name or batch_dir.name
    out_dir.mkdir(parents=True, exist_ok=True)

    groups: list[dict[str, Any]] = []
    for config in configs:
        integrated_path = batch_dir / config / "figures" / "integrated_flux" / f"{config}_poynting_integrated_python.tsv"
        groups.append(
            render_group(
                integrated_path=integrated_path,
                out_dir=out_dir,
                group_name=config,
                batch_name=batch_name,
                planes=planes,
                color_by_list=color_by_list,
                colormap=args.colormap,
                color_scale=args.color_scale,
                arrow_length_factor=args.arrow_length_factor,
                dpi=args.dpi,
                quiver_opacity=args.quiver_opacity,
            )
        )

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "batch_dir": str(batch_dir),
        "out_dir": str(out_dir),
        "render_method": "separate_cellmap_and_quiver_layers_alpha_composite",
        "drawn_together": False,
        "configs": configs,
        "planes": planes,
        "color_by": color_by_list,
        "colormap": args.colormap,
        "color_scale": args.color_scale,
        "arrow_length_factor": args.arrow_length_factor,
        "quiver_opacity": args.quiver_opacity,
        "groups": groups,
    }
    (out_dir / "cellmap_quiver_overlay_manifest.yaml").write_text(
        yaml.dump(manifest, default_flow_style=False, sort_keys=False),
        encoding="utf-8",
    )
    print(f"Wrote cellmap/quiver overlays to: {out_dir}")


if __name__ == "__main__":
    main()
