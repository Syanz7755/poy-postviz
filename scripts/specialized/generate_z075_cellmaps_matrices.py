"""Generate z=0.75 Poynting cellmaps and CSV value matrices.

This is a task-specific exporter for integrated Poynting-vector data produced
by ``plot_integrated_flux.py``. It deliberately does not interpolate. The CSV
matrices are full x/y rectangular matrices; coordinates without an evaluation
point are written as NaN.
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

from plot_integrated_flux import (
    ADAPTIVE_CELL_DISTANCE_THRESHOLD,
    LABEL_FONT_SIZE,
    SPARSE_OCCUPANCY_THRESHOLD,
    TICK_FONT_SIZE,
    adaptive_cell_thresholds,
    add_matched_colorbar,
    apply_matplotlib_theme,
    plot_threshold_cellmap,
)


DEFAULT_CONFIGS = ["ad1", "bd1", "cd1", "cd4", "dd1", "dd05", "ed1", "fd1"]
FIELDS = ["absSxy", "Sx", "Sy", "Sz", "direction_deg"]
SUPPORTED_FIGURE_FORMATS = {"png", "pdf", "svg", "eps"}
FIELD_LABELS = {
    "absSxy": "|Sxy|",
    "Sx": "Sx",
    "Sy": "Sy",
    "Sz": "Sz",
    "direction_deg": "direction [deg]",
}
ABSSXY_VMIN_PERCENTILE = 1.0
ABSSXY_VMAX_PERCENTILE = 92.0
SIGNED_SCALE_STEP = 0.2
SIGNED_FIXED_LIMIT = 0.3
BLUE_GRAY_RED_CYCLIC_CMAP = "blue_gray_red_gray_cyclic"
BLUE_GRAY_RED_CYCLIC_GRAY = "#bfbfbf"
BLUE_GRAY_RED_CYCLIC_BLUE = "#2626ff"
BLUE_GRAY_RED_CYCLIC_RED = "#ff2626"
TWILIGHT_SHIFTED_QUARTER_CMAP = "twilight_shifted_0p25"
TWILIGHT_SHIFT_FRACTION = 0.25
TWILIGHT_SHIFTED_THREE_QUARTER_CMAP = "twilight_shifted_0p75"
TWILIGHT_THREE_QUARTER_SHIFT_FRACTION = 0.75


@dataclass(frozen=True)
class MatrixResult:
    matrix: pd.DataFrame
    x_values: np.ndarray
    y_values: np.ndarray
    finite_count: int
    nan_count: int


def axis_cell_edges(coords: np.ndarray) -> np.ndarray:
    coords = np.asarray(coords, dtype=float)
    if coords.size == 0:
        raise ValueError("cannot build cell edges for empty coordinate axis")
    if coords.size == 1:
        return np.array([coords[0] - 0.5, coords[0] + 0.5], dtype=float)
    mids = (coords[:-1] + coords[1:]) / 2.0
    first_width = coords[1] - coords[0]
    last_width = coords[-1] - coords[-2]
    return np.concatenate([[coords[0] - first_width / 2.0], mids, [coords[-1] + last_width / 2.0]])


def format_coord(value: float) -> str:
    text = f"{float(value):.12g}"
    return "0" if text == "-0" else text


def read_integrated(
    batch_dir: Path | None,
    config: str,
    source_files: dict[str, Path] | None = None,
) -> tuple[pd.DataFrame, Path]:
    if source_files and config in source_files:
        path = source_files[config]
    elif batch_dir is not None:
        path = batch_dir / config / "figures" / "integrated_flux" / f"{config}_poynting_integrated_python.tsv"
    else:
        raise ValueError(f"no integrated-data source configured for {config}")
    if not path.exists():
        raise FileNotFoundError(f"missing integrated data for {config}: {path}")
    df = pd.read_csv(path, sep="\t")
    for col in ["x", "y", "z", "Sx", "Sy", "Sz"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["absSxy"] = np.hypot(df["Sx"], df["Sy"])
    df["direction_deg"] = np.mod(np.degrees(np.arctan2(df["Sy"], df["Sx"])), 360.0)
    zero_vec = np.isclose(df["absSxy"], 0.0, rtol=0.0, atol=0.0) | ~np.isfinite(df["absSxy"])
    df.loc[zero_vec, "direction_deg"] = np.nan
    return df, path


def build_matrix(plane_df: pd.DataFrame, field: str) -> MatrixResult:
    collapsed = plane_df.groupby(["y", "x"], as_index=False, sort=True)[field].mean(numeric_only=True)
    matrix = collapsed.pivot(index="y", columns="x", values=field).sort_index(ascending=True)
    x_values = matrix.columns.to_numpy(dtype=float)
    y_values = matrix.index.to_numpy(dtype=float)
    matrix.columns = [format_coord(v) for v in x_values]
    matrix.insert(0, "y\\x", [format_coord(v) for v in y_values])
    values = matrix.drop(columns=["y\\x"]).to_numpy(dtype=float)
    finite_count = int(np.isfinite(values).sum())
    nan_count = int(np.isnan(values).sum())
    return MatrixResult(matrix=matrix, x_values=x_values, y_values=y_values, finite_count=finite_count, nan_count=nan_count)


def color_settings(
    field: str,
    values: np.ndarray,
    signed_limit: float | None = None,
    requested_cmap: str | None = None,
) -> dict[str, Any]:
    import matplotlib.colors as mcolors

    finite = values[np.isfinite(values)]
    if field == "absSxy":
        positive = finite[finite > 0]
        if positive.size == 0:
            return {
                "cmap": requested_cmap or "turbo",
                "cmap_name": requested_cmap or "turbo",
                "norm": None,
                "scale": "log",
                "vmin": None,
                "vmax": None,
                "clip": None,
            }
        vmin = float(max(np.nanpercentile(positive, ABSSXY_VMIN_PERCENTILE), np.nanmin(positive)))
        vmax = float(np.nanpercentile(positive, ABSSXY_VMAX_PERCENTILE))
        if vmax <= vmin:
            vmax = float(np.nanmax(positive))
        if vmax <= vmin:
            vmax = vmin * 10.0
        return {
            "cmap": requested_cmap or "turbo",
            "cmap_name": requested_cmap or "turbo",
            "norm": mcolors.LogNorm(vmin=vmin, vmax=vmax, clip=True),
            "scale": "log",
            "vmin": vmin,
            "vmax": vmax,
            "clip": f"{ABSSXY_VMIN_PERCENTILE:g}-{ABSSXY_VMAX_PERCENTILE:g} percentile",
        }

    if field == "direction_deg":
        cmap_name = requested_cmap or "hsv"
        if cmap_name == BLUE_GRAY_RED_CYCLIC_CMAP:
            cmap_spec: str | mcolors.Colormap = mcolors.LinearSegmentedColormap.from_list(
                BLUE_GRAY_RED_CYCLIC_CMAP,
                [
                    (0.00, BLUE_GRAY_RED_CYCLIC_BLUE),
                    (0.25, BLUE_GRAY_RED_CYCLIC_GRAY),
                    (0.50, BLUE_GRAY_RED_CYCLIC_RED),
                    (0.75, BLUE_GRAY_RED_CYCLIC_GRAY),
                    (1.00, BLUE_GRAY_RED_CYCLIC_BLUE),
                ],
                N=257,
            )
        elif cmap_name in {
            TWILIGHT_SHIFTED_QUARTER_CMAP,
            TWILIGHT_SHIFTED_THREE_QUARTER_CMAP,
        }:
            import matplotlib.pyplot as plt

            base_cmap = plt.get_cmap("twilight")
            cyclic_positions = np.linspace(0.0, 1.0, 257)
            shift_fraction = (
                TWILIGHT_SHIFT_FRACTION
                if cmap_name == TWILIGHT_SHIFTED_QUARTER_CMAP
                else TWILIGHT_THREE_QUARTER_SHIFT_FRACTION
            )
            shifted_positions = np.mod(
                cyclic_positions + shift_fraction, 1.0
            )
            cmap_spec = mcolors.ListedColormap(
                base_cmap(shifted_positions),
                name=TWILIGHT_SHIFTED_QUARTER_CMAP,
            )
        else:
            cmap_spec = cmap_name
        return {
            "cmap": cmap_spec,
            "cmap_name": cmap_name,
            "norm": mcolors.Normalize(vmin=0.0, vmax=360.0),
            "scale": "linear",
            "vmin": 0.0,
            "vmax": 360.0,
            "clip": None,
        }

    max_abs = SIGNED_FIXED_LIMIT
    ticks = [-SIGNED_FIXED_LIMIT, 0.0, SIGNED_FIXED_LIMIT]
    signed_cmap = requested_cmap or mcolors.LinearSegmentedColormap.from_list(
        "signed_saturated_red_lightgrey_blue",
        ["#d7191c", "#f0f0f0", "#2c7bb6"],
        N=256,
    )
    return {
        "cmap": signed_cmap,
        "cmap_name": requested_cmap or "signed_saturated_red_lightgrey_blue",
        "norm": mcolors.TwoSlopeNorm(vmin=-max_abs, vcenter=0.0, vmax=max_abs),
        "scale": "linear",
        "vmin": -max_abs,
        "vmax": max_abs,
        "clip": None,
        "ticks": [float(t) for t in ticks],
    }


def plot_cellmap(
    out_path: Path,
    config: str,
    z_value: float,
    field: str,
    plane_df: pd.DataFrame,
    field_values: np.ndarray,
    result: MatrixResult,
    dpi: int,
    signed_limit: float | None = None,
    requested_cmap: str | None = None,
) -> dict[str, Any]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matrix_values = result.matrix.drop(columns=["y\\x"]).to_numpy(dtype=float)
    settings = color_settings(
        field,
        field_values,
        signed_limit=signed_limit,
        requested_cmap=requested_cmap,
    )
    norm_kw = {"norm": settings["norm"]} if settings["norm"] is not None else {}

    x_range = max(float(np.ptp(result.x_values)), 1.0)
    y_range = max(float(np.ptp(result.y_values)), 1.0)
    width = min(max(7.0 * x_range / y_range, 6.5), 14.0)
    apply_matplotlib_theme()
    fig, ax = plt.subplots(figsize=(width, 7.5))
    cmap_spec = settings["cmap"]
    cmap = (plt.get_cmap(cmap_spec) if isinstance(cmap_spec, str) else cmap_spec).copy()
    cmap.set_bad(color="white", alpha=0.0)

    grid_cells = int(len(result.x_values) * len(result.y_values))
    occupancy = float(len(plane_df) / grid_cells) if grid_cells else 0.0
    if occupancy < SPARSE_OCCUPANCY_THRESHOLD:
        x_threshold, y_threshold = adaptive_cell_thresholds(
            plane_df["x"].to_numpy(dtype=float),
            plane_df["y"].to_numpy(dtype=float),
        )
        artist = plot_threshold_cellmap(
            ax,
            plane_df["x"].to_numpy(dtype=float),
            plane_df["y"].to_numpy(dtype=float),
            field_values,
            cmap=cmap,
            norm_kw=norm_kw,
            x_threshold=x_threshold,
            y_threshold=y_threshold,
        )
        cellmap_method = "threshold_rectangles"
    else:
        x_edges = axis_cell_edges(result.x_values)
        y_edges = axis_cell_edges(result.y_values)
        artist = ax.pcolormesh(x_edges, y_edges, matrix_values, cmap=cmap, shading="flat", **norm_kw)
        x_threshold = None
        y_threshold = None
        cellmap_method = "rectilinear_pcolormesh"

    cbar = add_matched_colorbar(fig, ax, artist, label=FIELD_LABELS[field])
    if field in {"Sx", "Sy", "Sz"} and settings.get("ticks"):
        cbar.set_ticks(settings["ticks"])
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.tick_params(labelsize=TICK_FONT_SIZE)
    scale_text = "log" if settings["scale"] == "log" else "linear"
    ax.set_title(f"{config} z={z_value:.3f} {FIELD_LABELS[field]} ({scale_text})")
    save_kwargs: dict[str, Any] = {"bbox_inches": "tight"}
    if out_path.suffix.lower() == ".png":
        save_kwargs["dpi"] = dpi
    fig.savefig(out_path, **save_kwargs)
    plt.close(fig)

    entry = {
        "field": field,
        "figure": str(out_path),
        "figure_format": out_path.suffix.lower().lstrip("."),
        "colormap": settings["cmap_name"],
        "scale": settings["scale"],
        "vmin": settings["vmin"],
        "vmax": settings["vmax"],
        "clip": settings["clip"],
        "ticks": settings.get("ticks"),
        "finite_count": result.finite_count,
        "nan_count": result.nan_count,
        "shape": [int(len(result.y_values)), int(len(result.x_values))],
        "cellmap_method": cellmap_method,
        "adaptive_cell_base_distance_threshold": ADAPTIVE_CELL_DISTANCE_THRESHOLD if x_threshold is not None else None,
        "adaptive_cell_x_distance_threshold": x_threshold,
        "adaptive_cell_y_distance_threshold": y_threshold,
        "grid_cells": grid_cells,
        "grid_occupancy": occupancy,
    }
    if out_path.suffix.lower() == ".png":
        entry["png"] = str(out_path)
    return entry


def parse_figure_formats(raw: str) -> list[str]:
    formats = [item.strip().lower().lstrip(".") for item in raw.split(",") if item.strip()]
    if not formats:
        raise ValueError("at least one figure format is required")
    unsupported = [item for item in formats if item not in SUPPORTED_FIGURE_FORMATS]
    if unsupported:
        raise ValueError(f"unsupported figure format(s): {', '.join(unsupported)}")
    return list(dict.fromkeys(formats))


def generate(
    batch_dir: Path | None,
    out_dir: Path,
    configs: list[str],
    z_value: float,
    dpi: int,
    figure_format: str,
    fields: list[str] | None = None,
    field_colormaps: dict[str, list[str]] | None = None,
    source_files: dict[str, Path] | None = None,
    write_matrices: bool = True,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    figure_formats = parse_figure_formats(figure_format)
    selected_fields = fields or list(FIELDS)
    invalid_fields = [field for field in selected_fields if field not in FIELDS]
    if invalid_fields:
        raise ValueError(f"unsupported field(s): {', '.join(invalid_fields)}")
    field_colormaps = field_colormaps or {}
    invalid_cmap_fields = [field for field in field_colormaps if field not in FIELDS]
    if invalid_cmap_fields:
        raise ValueError(f"colormaps configured for unsupported field(s): {', '.join(invalid_cmap_fields)}")
    manifest: dict[str, Any] = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_batch_dir": str(batch_dir) if batch_dir is not None else None,
        "source_files": {key: str(value) for key, value in (source_files or {}).items()},
        "out_dir": str(out_dir),
        "z": z_value,
        "figure_format": figure_formats[0],
        "figure_formats": figure_formats,
        "vector_export": any(item in {"pdf", "svg", "eps"} for item in figure_formats),
        "method": "full rectangular x/y matrix from integrated evpoints; absent coordinates are NaN; no interpolation",
        "direction_definition": "direction_deg = mod(degrees(atan2(Sy, Sx)), 360), range [0, 360)",
        "selected_fields": selected_fields,
        "field_colormaps": field_colormaps,
        "write_matrices": write_matrices,
        "fields": {},
        "configs": [],
    }
    if any(
        BLUE_GRAY_RED_CYCLIC_CMAP in colormaps
        for colormaps in field_colormaps.values()
    ):
        manifest["custom_colormaps"] = {
            BLUE_GRAY_RED_CYCLIC_CMAP: {
                "type": "cyclic",
                "control_points": [
                    {"position": 0.00, "angle_deg": 0.0, "color": BLUE_GRAY_RED_CYCLIC_BLUE},
                    {"position": 0.25, "angle_deg": 90.0, "color": BLUE_GRAY_RED_CYCLIC_GRAY},
                    {"position": 0.50, "angle_deg": 180.0, "color": BLUE_GRAY_RED_CYCLIC_RED},
                    {"position": 0.75, "angle_deg": 270.0, "color": BLUE_GRAY_RED_CYCLIC_GRAY},
                    {"position": 1.00, "angle_deg": 360.0, "color": BLUE_GRAY_RED_CYCLIC_BLUE},
                ],
                "masked_background": "white",
            }
        }
    if any(
        TWILIGHT_SHIFTED_QUARTER_CMAP in colormaps
        for colormaps in field_colormaps.values()
    ):
        manifest.setdefault("custom_colormaps", {})[
            TWILIGHT_SHIFTED_QUARTER_CMAP
        ] = {
            "type": "cyclic",
            "base_colormap": "twilight",
            "shift_fraction": TWILIGHT_SHIFT_FRACTION,
            "mapping": "cmap_new(x) = twilight((x + 0.25) mod 1)",
            "angle_shift_deg": 90.0,
            "masked_background": "white",
        }
    if any(
        TWILIGHT_SHIFTED_THREE_QUARTER_CMAP in colormaps
        for colormaps in field_colormaps.values()
    ):
        manifest.setdefault("custom_colormaps", {})[
            TWILIGHT_SHIFTED_THREE_QUARTER_CMAP
        ] = {
            "type": "cyclic",
            "base_colormap": "twilight",
            "shift_fraction": TWILIGHT_THREE_QUARTER_SHIFT_FRACTION,
            "mapping": "cmap_new(x) = twilight((x + 0.75) mod 1)",
            "angle_shift_deg": 270.0,
            "masked_background": "white",
        }
    for config in configs:
        df, source_path = read_integrated(batch_dir, config, source_files=source_files)
        plane_df = df[np.abs(df["z"] - z_value) < 1e-6].copy()
        if plane_df.empty:
            raise ValueError(f"no points for config={config} at z={z_value}")

        config_dir = out_dir / config
        matrix_dir = config_dir / "csv_matrices"
        figure_dir = config_dir / "cellmaps"
        if write_matrices:
            matrix_dir.mkdir(parents=True, exist_ok=True)
        figure_dir.mkdir(parents=True, exist_ok=True)

        config_entry: dict[str, Any] = {
            "config": config,
            "source_file": str(source_path),
            "source_rows_at_z": int(len(plane_df)),
            "fields": [],
        }
        sx_sy_values = plane_df[["Sx", "Sy"]].to_numpy(dtype=float)
        sx_sy_finite = sx_sy_values[np.isfinite(sx_sy_values)]
        sx_sy_signed_limit = float(np.nanmax(np.abs(sx_sy_finite))) if sx_sy_finite.size else 1.0
        for field in selected_fields:
            result = build_matrix(plane_df, field)
            csv_path = matrix_dir / f"{config}_z{z_value:.3f}_{field}_matrix.csv" if write_matrices else None
            if csv_path is not None:
                result.matrix.to_csv(csv_path, index=False, na_rep="NaN")

            requested_colormaps: list[str | None] = field_colormaps.get(field, [None])
            for requested_cmap in requested_colormaps:
                plot_entries = []
                cmap_suffix = f"_{requested_cmap}" if requested_cmap else ""
                for item_format in figure_formats:
                    figure_path = (
                        figure_dir
                        / f"{config}_z{z_value:.3f}_{field}_cellmap{cmap_suffix}.{item_format}"
                    )
                    plot_entry = plot_cellmap(
                        figure_path,
                        config,
                        z_value,
                        field,
                        plane_df,
                        plane_df[field].to_numpy(dtype=float),
                        result,
                        dpi,
                        signed_limit=sx_sy_signed_limit if field in {"Sx", "Sy"} else None,
                        requested_cmap=requested_cmap,
                    )
                    plot_entries.append(plot_entry)

                field_entry = dict(plot_entries[0])
                field_entry["csv_matrix"] = str(csv_path) if csv_path is not None else None
                field_entry["figure_formats"] = figure_formats
                field_entry["figures"] = {
                    entry["figure_format"]: entry["figure"] for entry in plot_entries
                }
                for entry in plot_entries:
                    field_entry[entry["figure_format"]] = entry["figure"]
                config_entry["fields"].append(field_entry)
        manifest["configs"].append(config_entry)

    manifest["fields"] = {
        "absSxy": {
            "definition": "sqrt(Sx^2 + Sy^2)",
            "scale": "log",
            "colormaps": field_colormaps.get("absSxy", ["turbo"]),
            "clip": f"{ABSSXY_VMIN_PERCENTILE:g}-{ABSSXY_VMAX_PERCENTILE:g} percentile",
        },
        "Sx": {
            "definition": "signed integrated Sx",
            "scale": "linear",
            "colormap": "signed_saturated_red_lightgrey_blue",
            "zero_color": "light grey #f0f0f0",
            "scaler": "fixed symmetric range [-0.3, 0.3]",
            "ticks": [-0.3, 0.0, 0.3],
        },
        "Sy": {
            "definition": "signed integrated Sy",
            "scale": "linear",
            "colormap": "signed_saturated_red_lightgrey_blue",
            "zero_color": "light grey #f0f0f0",
            "scaler": "fixed symmetric range [-0.3, 0.3]",
            "ticks": [-0.3, 0.0, 0.3],
        },
        "Sz": {
            "definition": "signed integrated Sz",
            "scale": "linear",
            "colormap": "signed_saturated_red_lightgrey_blue",
            "zero_color": "light grey #f0f0f0",
            "scaler": "fixed symmetric range [-0.3, 0.3]",
            "ticks": [-0.3, 0.0, 0.3],
        },
        "direction_deg": {
            "definition": "mod(degrees(atan2(Sy, Sx)), 360)",
            "range": "[0, 360)",
            "scale": "linear",
            "colormaps": field_colormaps.get("direction_deg", ["hsv"]),
            "colormap_type": "cyclic/periodic",
        },
    }
    (out_dir / "poy_z075_cellmaps_matrices_manifest.yaml").write_text(
        yaml.dump(manifest, sort_keys=False, default_flow_style=False),
        encoding="utf-8",
    )
    return manifest


def parse_field_colormaps(specs: list[str]) -> dict[str, list[str]]:
    parsed: dict[str, list[str]] = {}
    for spec in specs:
        if "=" not in spec:
            raise ValueError(f"invalid --field-colormap value {spec!r}; expected FIELD=CMAP[,CMAP]")
        field, raw_cmaps = spec.split("=", 1)
        field = field.strip()
        cmaps = [item.strip() for item in raw_cmaps.split(",") if item.strip()]
        if not field or not cmaps:
            raise ValueError(f"invalid --field-colormap value {spec!r}; expected FIELD=CMAP[,CMAP]")
        parsed[field] = list(dict.fromkeys(cmaps))
    return parsed


def parse_source_files(specs: list[str]) -> dict[str, Path]:
    parsed: dict[str, Path] = {}
    for spec in specs:
        if "=" not in spec:
            raise ValueError(f"invalid --source value {spec!r}; expected CONFIG=PATH")
        config, raw_path = spec.split("=", 1)
        config = config.strip()
        raw_path = raw_path.strip()
        if not config or not raw_path:
            raise ValueError(f"invalid --source value {spec!r}; expected CONFIG=PATH")
        parsed[config] = Path(raw_path).resolve()
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--batch-dir",
        default=None,
        help="Batch directory containing config/figures/integrated_flux TSVs.",
    )
    parser.add_argument(
        "--source",
        action="append",
        default=[],
        metavar="CONFIG=PATH",
        help="Explicit integrated TSV for a config; repeat for sources from different batches.",
    )
    parser.add_argument("--out", required=True, help="Output directory under analysis.")
    parser.add_argument("--configs", default=",".join(DEFAULT_CONFIGS), help="Comma-separated config names.")
    parser.add_argument(
        "--fields",
        default=",".join(FIELDS),
        help="Comma-separated fields to export.",
    )
    parser.add_argument(
        "--field-colormap",
        action="append",
        default=[],
        metavar="FIELD=CMAP[,CMAP]",
        help="Override one field with one or more colormaps; repeat for additional fields.",
    )
    parser.add_argument("--z", type=float, default=0.75)
    parser.add_argument("--dpi", type=int, default=340)
    parser.add_argument(
        "--figure-format",
        default="png",
        help="Cellmap figure format(s), comma-separated. Use pdf/svg/eps for vector output. Default: png.",
    )
    parser.add_argument(
        "--skip-matrices",
        action="store_true",
        help="Do not duplicate CSV matrices when only figure variants are requested.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configs = [item.strip() for item in args.configs.split(",") if item.strip()]
    fields = [item.strip() for item in args.fields.split(",") if item.strip()]
    source_files = parse_source_files(args.source)
    batch_dir = Path(args.batch_dir).resolve() if args.batch_dir else None
    if batch_dir is None and any(config not in source_files for config in configs):
        missing = [config for config in configs if config not in source_files]
        raise ValueError(
            "--batch-dir is required unless --source is supplied for every config; "
            f"missing: {', '.join(missing)}"
        )
    manifest = generate(
        batch_dir,
        Path(args.out).resolve(),
        configs,
        args.z,
        args.dpi,
        args.figure_format,
        fields=fields,
        field_colormaps=parse_field_colormaps(args.field_colormap),
        source_files=source_files,
        write_matrices=not args.skip_matrices,
    )
    n_configs = len(manifest["configs"])
    field_variants = [field for item in manifest["configs"] for field in item["fields"]]
    n_figures = sum(len(field["figures"]) for field in field_variants)
    n_matrices = sum(1 for field in field_variants if field["csv_matrix"])
    print(
        f"Wrote {n_figures} cellmaps and {n_matrices} CSV matrices "
        f"for {n_configs} configs to: {manifest['out_dir']}"
    )


if __name__ == "__main__":
    main()
