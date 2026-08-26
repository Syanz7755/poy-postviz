"""Draw fixed-length Sxy quivers at every integrated evpoint.

This script is independent of the cellmap plotting path. It reads integrated
Poynting TSV files and writes quiver-only figures.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from plot_integrated_flux import TICK_FONT_SIZE, apply_matplotlib_theme


DEFAULT_DPI = 340
DEFAULT_ARROW_LENGTH_FACTOR = 0.35


def min_positive_spacing(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    unique = np.array(sorted(np.unique(values[np.isfinite(values)])), dtype=float)
    diffs = np.diff(unique)
    diffs = diffs[diffs > 1e-12]
    if diffs.size == 0:
        return 1.0
    return float(np.min(diffs))


def figure_size(x_vals: np.ndarray, y_vals: np.ndarray) -> tuple[float, float]:
    x_range = max(float(np.ptp(x_vals)), 1.0)
    y_range = max(float(np.ptp(y_vals)), 1.0)
    width = min(max(7.0 * x_range / y_range, 6.0), 15.0)
    return width, 7.8


def fixed_quiver_components(
    sx: np.ndarray,
    sy: np.ndarray,
    arrow_length: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mag = np.hypot(sx, sy)
    mask = np.isfinite(sx) & np.isfinite(sy) & np.isfinite(mag) & (mag > 0)
    u = np.full_like(sx, np.nan, dtype=float)
    v = np.full_like(sy, np.nan, dtype=float)
    u[mask] = sx[mask] / mag[mask] * arrow_length
    v[mask] = sy[mask] / mag[mask] * arrow_length
    return u, v, mask


def plot_group_quivers(
    integrated_path: Path,
    out_dir: Path,
    group_name: str,
    batch_name: str,
    planes: list[float],
    arrow_length: float | None = None,
    arrow_length_factor: float = DEFAULT_ARROW_LENGTH_FACTOR,
    dpi: int = DEFAULT_DPI,
) -> dict[str, Any]:
    import matplotlib

    matplotlib.use("Agg")
    apply_matplotlib_theme()
    import matplotlib.pyplot as plt

    integrated = pd.read_csv(integrated_path, sep="\t")
    for col in ["x", "y", "z", "Sx", "Sy"]:
        integrated[col] = pd.to_numeric(integrated[col], errors="coerce")

    out_dir.mkdir(parents=True, exist_ok=True)
    plane_tol = 1e-6
    plots: list[dict[str, Any]] = []

    for plane_z in planes:
        plane_df = integrated[np.abs(integrated["z"] - plane_z) < plane_tol].copy()
        if plane_df.empty:
            continue

        x = plane_df["x"].to_numpy(dtype=float)
        y = plane_df["y"].to_numpy(dtype=float)
        sx = plane_df["Sx"].to_numpy(dtype=float)
        sy = plane_df["Sy"].to_numpy(dtype=float)
        pixel_size = min(min_positive_spacing(x), min_positive_spacing(y))
        fixed_length = float(arrow_length) if arrow_length is not None else float(pixel_size * arrow_length_factor)
        u, v, mask = fixed_quiver_components(sx, sy, fixed_length)
        finite_xy = np.isfinite(x) & np.isfinite(y)
        mask = mask & finite_xy

        fig, ax = plt.subplots(figsize=figure_size(x[finite_xy], y[finite_xy]))
        if np.any(mask):
            ax.quiver(
                x[mask],
                y[mask],
                u[mask],
                v[mask],
                angles="xy",
                scale_units="xy",
                scale=1,
                color="k",
                width=0.0022,
                headwidth=3.2,
                headlength=4.2,
                headaxislength=3.7,
                pivot="middle",
                alpha=0.82,
            )

        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.tick_params(labelsize=TICK_FONT_SIZE)
        ax.set_aspect("equal")
        ax.set_title(f"{group_name} z={plane_z:.3f} fixed-length Sxy quiver")
        fname = f"{batch_name}_{group_name}_z{plane_z:.3f}_Sxy_fixedlength_evpoint_quiver.png"
        fig.savefig(out_dir / fname, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

        plots.append(
            {
                "filename": fname,
                "plane_z": plane_z,
                "n_evpoints": int(len(plane_df)),
                "n_quivers": int(np.count_nonzero(mask)),
                "pixel_size": pixel_size,
                "arrow_length": fixed_length,
                "arrow_length_factor": None if arrow_length is not None else arrow_length_factor,
                "arrow_length_mode": "absolute" if arrow_length is not None else "pixel_size_factor",
            }
        )

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "group_name": group_name,
        "batch_name": batch_name,
        "integrated_path": str(integrated_path),
        "render_method": "fixed_length_evpoint_quiver",
        "cellmap_independent": True,
        "plots": plots,
    }
    (out_dir / "fixed_evpoint_quiver_manifest.yaml").write_text(
        yaml.dump(manifest, default_flow_style=False, sort_keys=False),
        encoding="utf-8",
    )
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-dir", required=True, help="Batch directory containing config/figures/integrated_flux.")
    parser.add_argument("--out", required=True, help="Output directory for quiver-only figures.")
    parser.add_argument("--configs", default="cd1,cd4", help="Comma-separated configs.")
    parser.add_argument("--planes", default="0.75,1.0,1.5", help="Comma-separated z planes.")
    parser.add_argument("--batch-name", default=None, help="Filename prefix. Defaults to batch directory name.")
    parser.add_argument("--arrow-length", type=float, default=None, help="Absolute data-coordinate arrow length.")
    parser.add_argument("--arrow-length-factor", type=float, default=DEFAULT_ARROW_LENGTH_FACTOR, help="Arrow length as a factor of per-plane min spacing.")
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    batch_dir = Path(args.batch_dir).resolve()
    out_root = Path(args.out).resolve()
    configs = [item.strip() for item in args.configs.split(",") if item.strip()]
    planes = [float(item.strip()) for item in args.planes.split(",") if item.strip()]
    batch_name = args.batch_name or batch_dir.name
    results: list[dict[str, Any]] = []

    for config in configs:
        integrated_path = batch_dir / config / "figures" / "integrated_flux" / f"{config}_poynting_integrated_python.tsv"
        config_out = out_root / config
        results.append(
            plot_group_quivers(
                integrated_path=integrated_path,
                out_dir=config_out,
                group_name=config,
                batch_name=batch_name,
                planes=planes,
                arrow_length=args.arrow_length,
                arrow_length_factor=args.arrow_length_factor,
                dpi=args.dpi,
            )
        )

    out_root.mkdir(parents=True, exist_ok=True)
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "batch_dir": str(batch_dir),
        "out_dir": str(out_root),
        "render_method": "fixed_length_evpoint_quiver",
        "cellmap_independent": True,
        "configs": configs,
        "planes": planes,
        "arrow_length": args.arrow_length,
        "arrow_length_factor": args.arrow_length_factor,
        "groups": results,
    }
    (out_root / "fixed_evpoint_quiver_manifest.yaml").write_text(
        yaml.dump(manifest, default_flow_style=False, sort_keys=False),
        encoding="utf-8",
    )
    print(f"Wrote fixed-length evpoint quivers to: {out_root}")


if __name__ == "__main__":
    main()
