"""Export diagrams for every existing evpoint text file without regenerating points."""

from __future__ import annotations

import argparse
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from gen_evpoints import parse_float_list, plot_evpoints


GROUP_RE = re.compile(r"^# Evaluation points for group (?P<group>.+), d = (?P<d>[-+0-9.eE]+)$")
GEOMETRY_RE = re.compile(
    r"^# W = (?P<W>[-+0-9.eE]+), H = (?P<H>[-+0-9.eE]+), delta = (?P<delta>[-+0-9.eE]+)$"
)
SUPPORTED_FORMATS = {"png", "svg", "pdf"}


def parse_formats(raw: str) -> list[str]:
    formats = [item.strip().lower().lstrip(".") for item in raw.split(",") if item.strip()]
    formats = list(dict.fromkeys(formats))
    if not formats:
        raise ValueError("at least one output format is required")
    unsupported = [item for item in formats if item not in SUPPORTED_FORMATS]
    if unsupported:
        raise ValueError(f"unsupported output format(s): {', '.join(unsupported)}")
    return formats


def read_evpoint_file(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    header_lines = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            header_lines.append(line.strip())

    group_match = next((GROUP_RE.match(line) for line in header_lines if GROUP_RE.match(line)), None)
    geometry_match = next(
        (GEOMETRY_RE.match(line) for line in header_lines if GEOMETRY_RE.match(line)),
        None,
    )
    if group_match is None or geometry_match is None:
        raise ValueError(f"missing group/geometry metadata in {path}")

    values = np.loadtxt(path, comments="#", dtype=float)
    if values.ndim == 1:
        values = values.reshape(1, -1)
    if values.shape[1] != 3:
        raise ValueError(f"expected x/y/z columns in {path}, got shape {values.shape}")

    z_point_data = []
    for z_value in sorted(np.unique(values[:, 2])):
        layer = values[np.isclose(values[:, 2], z_value, rtol=0.0, atol=1e-12)]
        z_point_data.append(
            {
                "z": float(z_value),
                "points": [tuple(row) for row in layer[:, :2]],
            }
        )

    metadata = {
        "group": group_match.group("group"),
        "d": float(group_match.group("d")),
        "W": float(geometry_match.group("W")),
        "H": float(geometry_match.group("H")),
        "delta": float(geometry_match.group("delta")),
    }
    return metadata, z_point_data


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path("evpoints"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--formats", default="png,svg,pdf")
    parser.add_argument("--plot-z", default="0.75")
    parser.add_argument("--point-size", type=float, default=0.2)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--show-legend", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    out_root = args.out.resolve()
    figure_dir = out_root / "figures" / "evpoints"
    figure_dir.mkdir(parents=True, exist_ok=True)
    formats = parse_formats(args.formats)
    plot_z = parse_float_list(args.plot_z)

    input_paths = sorted(input_dir.glob("*.txt"), key=lambda path: path.name.lower())
    if not input_paths:
        raise FileNotFoundError(f"no evpoint text files found in {input_dir}")

    rows = []
    for input_path in input_paths:
        metadata, z_point_data = read_evpoint_file(input_path)
        output_stub = figure_dir / f"{input_path.stem}.png"
        written = plot_evpoints(
            output_stub,
            z_point_data,
            metadata["group"],
            metadata["H"],
            metadata["delta"],
            metadata["d"],
            plot_z=plot_z,
            point_size=args.point_size,
            hide_legend=not args.show_legend,
            dpi=args.dpi,
            plot_formats=formats,
        )
        plotted_z = {
            round(float(value), 12) for value in plot_z
        } if plot_z is not None else None
        plotted_points = sum(
            len(item["points"])
            for item in z_point_data
            if plotted_z is None or round(float(item["z"]), 12) in plotted_z
        )
        rows.append(
            {
                "input": str(input_path),
                **metadata,
                "available_z": [float(item["z"]) for item in z_point_data],
                "plot_z": plot_z,
                "plotted_points": plotted_points,
                "point_size": args.point_size,
                "formats": formats,
                "outputs": {path.suffix.lstrip("."): str(path) for path in written},
            }
        )

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "input_dir": str(input_dir),
        "out_dir": str(out_root),
        "method": "plot existing evpoint coordinates; no point regeneration",
        "plot_z": plot_z,
        "point_size": args.point_size,
        "hide_legend": not args.show_legend,
        "dpi": args.dpi,
        "formats": formats,
        "vector_formats": [item for item in formats if item in {"svg", "pdf"}],
        "input_count": len(input_paths),
        "figure_count": sum(len(row["outputs"]) for row in rows),
        "diagrams": rows,
    }
    manifest_path = out_root / "evpoint_diagram_manifest.yaml"
    manifest_path.write_text(
        yaml.dump(manifest, default_flow_style=False, sort_keys=False),
        encoding="utf-8",
    )
    print(
        f"Wrote {manifest['figure_count']} figures for {manifest['input_count']} "
        f"evpoint inputs to: {out_root}"
    )


if __name__ == "__main__":
    main()
