"""Regenerate evpoints/high106-2_d1.0.txt and write PNG/SVG/PDF scatter diagrams.

One-off driver that reuses helpers from gen_evpoints.py without modifying the
shared groups dict. Keeps the geometry header from the existing
high106-2_d1.0.txt (H=10, delta=0.6, d=1.0, denser_near=0.1,
denser_min_ratio=sqrt(6)) but caps z=0.75 at 9000 points so the regenerated
file satisfies the "< 10000" rule. The current high106-2 file only contains a
z=0.75 layer, so we generate only z=0.75 to preserve coverage.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import yaml

from gen_evpoints import (
    generate_evpoints_adaptive_for_z,
    plot_evpoints,
    write_evpoints,
)


# Geometry (from existing high106-2_d1.0.txt header).
H = 10
DELTA = 0.6
D = 1.0
DENSER_NEAR = 0.1
DENSER_MIN_RATIO = 2.449489742783178  # sqrt(6); preserves prior refinement ratio

# Cap below 10000.
MAX_POINTS_Z075 = 9000

# Only z=0.75 (matches the existing file's coverage).
Z_LAYERS = [0.75]

GROUP = "high106-2"
DISPLACEMENT_LABEL = "d1.0"
OUTPUT_EVPOINTS = Path(f"evpoints/{GROUP}_{DISPLACEMENT_LABEL}.txt")
PLOT_BATCH_NAME = f"high106-2-z075-{MAX_POINTS_Z075}-{datetime.now().strftime('%Y%m%d_%H%M%S')}"
PLOT_OUT_ROOT = Path("evpoints_plots") / PLOT_BATCH_NAME
PLOT_Z = [0.75]
PLOT_POINT_SIZE = 1.0
PLOT_HIDE_LEGEND = True
PLOT_DPI = 300
PLOT_FORMATS = ["png", "svg", "pdf"]


def main() -> None:
    OUTPUT_EVPOINTS.parent.mkdir(parents=True, exist_ok=True)

    z_point_data = []
    for z in Z_LAYERS:
        points, used_spacing, n_iter = generate_evpoints_adaptive_for_z(
            group=GROUP,
            H=H,
            delta=DELTA,
            d=D,
            z=z,
            max_points=MAX_POINTS_Z075,
            denser_near=DENSER_NEAR,
            denser_min_ratio=DENSER_MIN_RATIO,
        )
        z_point_data.append(
            {
                "z": z,
                "points": points,
                "spacing": used_spacing,
                "adaptive_iter": n_iter,
                "max_points": MAX_POINTS_Z075,
            }
        )

    # Defensive guard: ensure cap holds.
    n_points = len(z_point_data[0]["points"])
    if n_points >= 10000:
        raise SystemExit(
            f"Regenerated file has {n_points} points at z=0.75; cap of 10000 violated."
        )
    if n_points > MAX_POINTS_Z075:
        raise SystemExit(
            f"Adaptive spacing failed to honor max_points={MAX_POINTS_Z075}; "
            f"got {n_points}."
        )

    write_evpoints(
        str(OUTPUT_EVPOINTS),
        z_point_data,
        GROUP,
        H,
        DELTA,
        D,
        denser_near=DENSER_NEAR,
        denser_min_ratio=DENSER_MIN_RATIO,
    )

    for item in z_point_data:
        print(
            f"z={item['z']}: "
            f"{len(item['points'])} points, "
            f"max={item['max_points']}, "
            f"spacing={item['spacing']:.6f}, "
            f"adaptive_iter={item['adaptive_iter']}"
        )
    print(f"wrote {OUTPUT_EVPOINTS}")

    figure_dir = PLOT_OUT_ROOT / "all" / "figures" / "evpoints"
    figure_dir.mkdir(parents=True, exist_ok=True)

    output_stub = figure_dir / f"{GROUP}_{DISPLACEMENT_LABEL}.png"
    written = plot_evpoints(
        str(output_stub),
        z_point_data,
        GROUP,
        H,
        DELTA,
        D,
        plot_z=PLOT_Z,
        point_size=PLOT_POINT_SIZE,
        hide_legend=PLOT_HIDE_LEGEND,
        dpi=PLOT_DPI,
        plot_formats=PLOT_FORMATS,
    )

    plotted_points = sum(len(item["points"]) for item in z_point_data)
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "method": "regenerate evpoints via gen_evpoints helpers; plot scatter diagrams",
        "input": str(OUTPUT_EVPOINTS),
        "group": GROUP,
        "d": D,
        "H": H,
        "delta": DELTA,
        "denser_near": DENSER_NEAR,
        "denser_min_ratio": DENSER_MIN_RATIO,
        "max_points_z075": MAX_POINTS_Z075,
        "available_z": [item["z"] for item in z_point_data],
        "plot_z": PLOT_Z,
        "plotted_points": plotted_points,
        "point_size": PLOT_POINT_SIZE,
        "hide_legend": PLOT_HIDE_LEGEND,
        "dpi": PLOT_DPI,
        "formats": PLOT_FORMATS,
        "outputs": {path.suffix.lstrip("."): str(path) for path in written},
    }
    manifest_path = PLOT_OUT_ROOT / "diagram_manifest.yaml"
    manifest_path.write_text(
        yaml.dump(manifest, default_flow_style=False, sort_keys=False),
        encoding="utf-8",
    )

    print(f"wrote {len(written)} diagram files to: {figure_dir}")
    for path in written:
        print(f"  - {path}")
    print(f"manifest: {manifest_path}")


if __name__ == "__main__":
    main()