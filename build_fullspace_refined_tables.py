"""Build full rectangular coordinate-value tables from integrated evpoints.

This script is intentionally independent of the plotting scripts. It writes
tabular coordinate-value matrices; it does not draw figures.
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


FIELD_COLS = ["Sx", "Sy", "Sz"]
DEFAULT_SQUARE_TOLERANCE = 0.25
DEFAULT_STEP_TOLERANCE = 1e-5
DEFAULT_MAX_INTERPOLATE_STEPS = 2
DEFAULT_SNAP_TOL_FACTOR = 0.05


@dataclass(frozen=True)
class PlaneResult:
    table: pd.DataFrame
    z: float
    pixel_size: float
    nx: int
    ny: int
    original_points: int
    interpolated_points: int
    nan_points: int
    interpolated_rectangles: int


def point_key(x: float, y: float) -> tuple[float, float]:
    return (round(float(x), 12), round(float(y), 12))


def min_positive_spacing(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    unique = np.array(sorted(np.unique(values[np.isfinite(values)])), dtype=float)
    diffs = np.diff(unique)
    diffs = diffs[diffs > 1e-12]
    if diffs.size == 0:
        return 1.0
    return float(np.min(diffs))


def grid_axis(min_value: float, max_value: float, step: float) -> np.ndarray:
    span = float(max_value - min_value)
    if span <= 1e-12:
        return np.array([float(min_value)])
    n_steps = int(round(span / step))
    if abs(span / step - n_steps) > DEFAULT_STEP_TOLERANCE:
        n_steps = int(np.ceil(span / step))
    axis = float(min_value) + np.arange(n_steps + 1, dtype=float) * step
    if axis[-1] < max_value - step * DEFAULT_STEP_TOLERANCE:
        axis = np.append(axis, float(max_value))
    else:
        axis[-1] = float(max_value)
    return axis


def snap_values_to_axis(values: np.ndarray, axis: np.ndarray, tolerance: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    axis = np.asarray(axis, dtype=float)
    snapped = values.copy()
    insert_at = np.searchsorted(axis, values)
    for i, value in enumerate(values):
        candidates: list[int] = []
        if insert_at[i] < len(axis):
            candidates.append(int(insert_at[i]))
        if insert_at[i] > 0:
            candidates.append(int(insert_at[i] - 1))
        if not candidates:
            continue
        nearest = min(candidates, key=lambda idx: abs(axis[idx] - value))
        if abs(axis[nearest] - value) <= tolerance:
            snapped[i] = axis[nearest]
    return snapped


def segment_points(start: float, stop: float, step: float) -> np.ndarray:
    span = float(stop - start)
    if abs(span) <= 1e-12:
        return np.array([float(start)])
    n_steps = int(round(abs(span) / step))
    if n_steps < 1 or abs(abs(span) / step - n_steps) > DEFAULT_STEP_TOLERANCE:
        n_steps = int(np.ceil(abs(span) / step))
    return np.linspace(float(start), float(stop), n_steps + 1)


def consecutive_values_by_axis(df: pd.DataFrame) -> tuple[dict[float, list[float]], dict[float, list[float]]]:
    xs_by_y: dict[float, list[float]] = {}
    ys_by_x: dict[float, list[float]] = {}
    for y, group in df.groupby("y", sort=True):
        xs_by_y[round(float(y), 12)] = sorted({round(float(x), 12) for x in group["x"]})
    for x, group in df.groupby("x", sort=True):
        ys_by_x[round(float(x), 12)] = sorted({round(float(y), 12) for y in group["y"]})
    return xs_by_y, ys_by_x


def pair_set(axis_values: dict[float, list[float]]) -> set[tuple[float, float, float]]:
    pairs: set[tuple[float, float, float]] = set()
    for fixed, values in axis_values.items():
        for start, stop in zip(values[:-1], values[1:]):
            pairs.add((fixed, start, stop))
    return pairs


def valid_small_rectangle(
    width: float,
    height: float,
    pixel_size: float,
    max_interpolate_steps: int,
    square_tolerance: float,
) -> bool:
    if pixel_size <= 0:
        return False
    x_steps = int(round(abs(width) / pixel_size))
    y_steps = int(round(abs(height) / pixel_size))
    if x_steps <= 1 and y_steps <= 1:
        return False
    if x_steps > max_interpolate_steps or y_steps > max_interpolate_steps:
        return False
    scale = max(abs(width), abs(height))
    if scale <= 1e-12:
        return False
    return abs(abs(width) - abs(height)) / scale <= square_tolerance


def bilinear(v00: np.ndarray, v10: np.ndarray, v01: np.ndarray, v11: np.ndarray, tx: float, ty: float) -> np.ndarray:
    return (
        (1 - tx) * (1 - ty) * v00
        + tx * (1 - ty) * v10
        + (1 - tx) * ty * v01
        + tx * ty * v11
    )


def first_round_interpolate(
    averaged: pd.DataFrame,
    pixel_size: float,
    max_interpolate_steps: int,
    square_tolerance: float,
) -> tuple[dict[tuple[float, float], tuple[np.ndarray, str]], int]:
    values_by_point = {
        point_key(row.x, row.y): np.array([row.Sx, row.Sy, row.Sz], dtype=float)
        for row in averaged.itertuples(index=False)
    }
    out: dict[tuple[float, float], tuple[np.ndarray, str]] = {
        key: (value, "original") for key, value in values_by_point.items()
    }
    xs_by_y, ys_by_x = consecutive_values_by_axis(averaged)
    x_pairs = pair_set(xs_by_y)
    y_pairs = pair_set(ys_by_x)
    rectangle_count = 0

    for y0, xs_row in xs_by_y.items():
        for x0, x1 in zip(xs_row[:-1], xs_row[1:]):
            y_values = ys_by_x.get(x0, [])
            for y0_col, y1 in zip(y_values[:-1], y_values[1:]):
                if abs(y0_col - y0) > 1e-9:
                    continue
                if (x1, y0, y1) not in y_pairs:
                    continue
                if (y1, x0, x1) not in x_pairs:
                    continue
                if not valid_small_rectangle(x1 - x0, y1 - y0, pixel_size, max_interpolate_steps, square_tolerance):
                    continue

                corner_keys = [point_key(x0, y0), point_key(x1, y0), point_key(x0, y1), point_key(x1, y1)]
                if any(key not in values_by_point for key in corner_keys):
                    continue
                v00, v10, v01, v11 = [values_by_point[key] for key in corner_keys]
                rectangle_count += 1

                for x in segment_points(x0, x1, pixel_size):
                    tx = 0.0 if abs(x1 - x0) <= 1e-12 else float((x - x0) / (x1 - x0))
                    for y in segment_points(y0, y1, pixel_size):
                        key = point_key(x, y)
                        if key in values_by_point:
                            continue
                        ty = 0.0 if abs(y1 - y0) <= 1e-12 else float((y - y0) / (y1 - y0))
                        out[key] = (bilinear(v00, v10, v01, v11, tx, ty), "interpolated_round1")

    return out, rectangle_count


def build_plane_table(
    plane_df: pd.DataFrame,
    plane_z: float,
    group_name: str,
    max_interpolate_steps: int,
    square_tolerance: float,
    snap_tol_factor: float,
) -> PlaneResult:
    averaged = (
        plane_df.groupby(["x", "y"], as_index=False, sort=True)[FIELD_COLS]
        .mean(numeric_only=True)
        .dropna(how="all", subset=FIELD_COLS)
    )
    x_values = averaged["x"].to_numpy(dtype=float)
    y_values = averaged["y"].to_numpy(dtype=float)
    pixel_size = float(min(min_positive_spacing(x_values), min_positive_spacing(y_values)))
    snap_tolerance = pixel_size * snap_tol_factor
    raw_x_grid = grid_axis(float(np.min(x_values)), float(np.max(x_values)), pixel_size)
    raw_y_grid = grid_axis(float(np.min(y_values)), float(np.max(y_values)), pixel_size)
    averaged = averaged.copy()
    averaged["x"] = snap_values_to_axis(x_values, raw_x_grid, snap_tolerance)
    averaged["y"] = snap_values_to_axis(y_values, raw_y_grid, snap_tolerance)
    averaged = averaged.groupby(["x", "y"], as_index=False, sort=True)[FIELD_COLS].mean(numeric_only=True)
    x_values = averaged["x"].to_numpy(dtype=float)
    y_values = averaged["y"].to_numpy(dtype=float)

    first_round, rectangle_count = first_round_interpolate(
        averaged,
        pixel_size=pixel_size,
        max_interpolate_steps=max_interpolate_steps,
        square_tolerance=square_tolerance,
    )

    x_grid = grid_axis(float(np.min(x_values)), float(np.max(x_values)), pixel_size)
    y_grid = grid_axis(float(np.min(y_values)), float(np.max(y_values)), pixel_size)
    rows: list[dict[str, Any]] = []
    for y in y_grid:
        for x in x_grid:
            key = point_key(x, y)
            if key in first_round:
                values, stage = first_round[key]
            else:
                values = np.array([np.nan, np.nan, np.nan], dtype=float)
                stage = "nan_round2_fullspace"
            sx, sy, sz = values
            abs_sxy = float(np.hypot(sx, sy)) if np.isfinite(sx) and np.isfinite(sy) else np.nan
            if np.isfinite(abs_sxy) and abs_sxy > 0:
                sxy_unit_x = float(sx / abs_sxy)
                sxy_unit_y = float(sy / abs_sxy)
            else:
                sxy_unit_x = np.nan
                sxy_unit_y = np.nan
            rows.append(
                {
                    "group_name": group_name,
                    "z": plane_z,
                    "x": float(x),
                    "y": float(y),
                    "Sx": sx,
                    "Sy": sy,
                    "Sz": sz,
                    "absSxy": abs_sxy,
                    "Sxy_unit_x": sxy_unit_x,
                    "Sxy_unit_y": sxy_unit_y,
                    "value_stage": stage,
                }
            )

    table = pd.DataFrame(rows)
    return PlaneResult(
        table=table,
        z=plane_z,
        pixel_size=pixel_size,
        nx=len(x_grid),
        ny=len(y_grid),
        original_points=int((table["value_stage"] == "original").sum()),
        interpolated_points=int((table["value_stage"] == "interpolated_round1").sum()),
        nan_points=int((table["value_stage"] == "nan_round2_fullspace").sum()),
        interpolated_rectangles=int(rectangle_count),
    )


def write_matrix_tables(group_table: pd.DataFrame, out_dir: Path, group_name: str, fields: list[str]) -> list[str]:
    matrix_dir = out_dir / "matrices"
    matrix_dir.mkdir(parents=True, exist_ok=True)
    paths: list[str] = []
    for z, plane_df in group_table.groupby("z", sort=True):
        for field in fields:
            matrix = plane_df.pivot(index="y", columns="x", values=field).sort_index(ascending=True)
            matrix.insert(0, "y", matrix.index)
            out_path = matrix_dir / f"{group_name}_z{float(z):.3f}_{field}_fullspace_matrix.tsv"
            matrix.to_csv(out_path, sep="\t", index=False)
            paths.append(str(out_path.relative_to(out_dir)))
    return paths


def build_group(
    integrated_path: Path,
    out_dir: Path,
    group_name: str,
    planes: list[float],
    max_interpolate_steps: int,
    square_tolerance: float,
    snap_tol_factor: float,
    matrix_fields: list[str],
) -> dict[str, Any]:
    integrated = pd.read_csv(integrated_path, sep="\t")
    for col in ["x", "y", "z", *FIELD_COLS]:
        integrated[col] = pd.to_numeric(integrated[col], errors="coerce")

    out_dir.mkdir(parents=True, exist_ok=True)
    plane_results: list[PlaneResult] = []
    plane_tol = 1e-6
    for plane_z in planes:
        plane_df = integrated[np.abs(integrated["z"] - plane_z) < plane_tol].copy()
        if plane_df.empty:
            continue
        plane_results.append(
            build_plane_table(
                plane_df,
                plane_z=plane_z,
                group_name=group_name,
                max_interpolate_steps=max_interpolate_steps,
                square_tolerance=square_tolerance,
                snap_tol_factor=snap_tol_factor,
            )
        )

    group_table = pd.concat([result.table for result in plane_results], ignore_index=True)
    long_path = out_dir / f"{group_name}_fullspace_coordinate_value_table.tsv"
    group_table.to_csv(long_path, sep="\t", index=False)
    matrix_paths = write_matrix_tables(group_table, out_dir, group_name, matrix_fields)

    return {
        "group_name": group_name,
        "integrated_path": str(integrated_path),
        "max_interpolate_steps": max_interpolate_steps,
        "refine_enabled": max_interpolate_steps > 1,
        "long_table": str(long_path),
        "matrix_tables": matrix_paths,
        "planes": [
            {
                "z": result.z,
                "pixel_size": result.pixel_size,
                "nx": result.nx,
                "ny": result.ny,
                "total_points": int(len(result.table)),
                "original_points": result.original_points,
                "interpolated_round1_points": result.interpolated_points,
                "nan_round2_fullspace_points": result.nan_points,
                "interpolated_rectangles": result.interpolated_rectangles,
            }
            for result in plane_results
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-dir", required=True, help="Batch directory containing cd*/figures/integrated_flux.")
    parser.add_argument("--out", required=True, help="New output directory for full-space tables.")
    parser.add_argument("--configs", default="cd1,cd4", help="Comma-separated configs.")
    parser.add_argument("--planes", default="0.75,1.0,1.5", help="Comma-separated z planes.")
    parser.add_argument("--max-interpolate-steps", type=int, default=DEFAULT_MAX_INTERPOLATE_STEPS)
    parser.add_argument("--no-refine-configs", default="", help="Comma-separated configs that should skip first-round interpolation.")
    parser.add_argument("--square-tolerance", type=float, default=DEFAULT_SQUARE_TOLERANCE)
    parser.add_argument("--snap-tol-factor", type=float, default=DEFAULT_SNAP_TOL_FACTOR)
    parser.add_argument(
        "--matrix-fields",
        default="absSxy,Sxy_unit_x,Sxy_unit_y,Sx,Sy,Sz",
        help="Comma-separated matrix fields.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    batch_dir = Path(args.batch_dir).resolve()
    out_dir = Path(args.out).resolve()
    configs = [item.strip() for item in args.configs.split(",") if item.strip()]
    no_refine_configs = {item.strip() for item in args.no_refine_configs.split(",") if item.strip()}
    planes = [float(item.strip()) for item in args.planes.split(",") if item.strip()]
    matrix_fields = [item.strip() for item in args.matrix_fields.split(",") if item.strip()]

    results = []
    for config in configs:
        integrated_path = batch_dir / config / "figures" / "integrated_flux" / f"{config}_poynting_integrated_python.tsv"
        group_out = out_dir / config
        config_max_interpolate_steps = 1 if config in no_refine_configs else args.max_interpolate_steps
        results.append(
            build_group(
                integrated_path=integrated_path,
                out_dir=group_out,
                group_name=config,
                planes=planes,
                max_interpolate_steps=config_max_interpolate_steps,
                square_tolerance=args.square_tolerance,
                snap_tol_factor=args.snap_tol_factor,
                matrix_fields=matrix_fields,
            )
        )

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "batch_dir": str(batch_dir),
        "out_dir": str(out_dir),
        "method": "round1 local small-rectangle bilinear interpolation, round2 full rectangular grid with NaN",
        "pixel_size": "minimum positive spacing across x and y for each group/plane",
        "round2_nan_rule": "coordinates absent after round1 are written with NaN values",
        "max_interpolate_steps": args.max_interpolate_steps,
        "no_refine_configs": sorted(no_refine_configs),
        "square_tolerance": args.square_tolerance,
        "snap_tol_factor": args.snap_tol_factor,
        "matrix_fields": matrix_fields,
        "groups": results,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "fullspace_refined_tables_manifest.yaml").write_text(
        yaml.dump(manifest, default_flow_style=False, sort_keys=False),
        encoding="utf-8",
    )
    print(f"Wrote full-space refined tables to: {out_dir}")


if __name__ == "__main__":
    main()
