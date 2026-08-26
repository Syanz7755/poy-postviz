"""
Generate evaluation points for scuff-neq in vacuum regions.

Usage:
  python gen_evpoints.py
  python gen_evpoints.py [--denser-near 0.2] [--denser-min-ratio 2]
  python gen_evpoints.py --with-plot [--out-dir DIR] [--batch-name NAME]

Output:
  evpoints/{group}_d{d}.txt                        (always, for scuff-neq input)
  {out-dir}/{batch}/all/figures/evpoints/{group}_d{d}.png  (only --with-plot)
"""

import os
import argparse
from pathlib import Path

import numpy as np

# Parameters
W = 1.0
base_min_spacing = 0.08
min_edge_dist = 0.01
outexpand = 1.0
rep = 3
default_denser_near = 0.2
default_denser_min_ratio = 2.0

# Each z has its own maximum point threshold.
# Example:
#   z_values = [0.75, 1.0]
#   max_points_per_z = [10000, 2000]
z_values = [0.75,1.0,1.5]
max_points_per_z = [8000, 800, 800]
max_points_per_z_overrides = {
    ("c", 1.0): [13000, 800, 800],
    ("c", 4.0): [11000, 800, 800],
}

# Adaptive spacing parameters
spacing_growth = 1.05
max_adaptive_iter = 200

# Group definitions
groups = {
    "a": {"H": 5, "delta": 0.6, "displacements": [1.0, 2.0, 4.0, 5.0, 10.0]},
    "b": {"H": 5, "delta": 0.2, "displacements": [1.0, 2.0, 4.0, 5.0, 10.0]},
    "c": {"H": 5, "delta": 0.1, "displacements": [1.0, 2.0, 4.0, 5.0, 10.0]},
    "d": {"H": 2, "delta": 0.6, "displacements": [0.5, 1.0, 2.0, 5.0]},
    "e": {"H": 2, "delta": 0.2, "displacements": [0.5, 1.0, 2.0, 5.0]},
    "f": {"H": 2, "delta": 0.1, "displacements": [0.5, 1.0, 2.0, 5.0]},
}


def point_in_rect(x, y, x0, y0, w, h, margin=0):
    """Check if point (x, y) is inside rectangle with margin.
    Rectangle is defined by bottom-mid-point (x0, y0) and size (w, h).
    """
    x_min = x0 - w / 2 - margin
    x_max = x0 + w / 2 + margin
    y_min = y0 - margin
    y_max = y0 + h + margin
    return x_min <= x <= x_max and y_min <= y <= y_max


def solid_rectangles(H, delta, d):
    """Return solid rectangles as (bottom-middle x/y, width, height)."""
    p = 2 * (W + delta)
    rectangles = []

    for n in range(rep):
        rectangles.extend(
            [
                (n * p + W / 4, d, W / 2, H),
                ((n + 1) * p - W / 4, d, W / 2, H),
                ((2 * n + 1) * p / 2, 0, W, H),
            ]
        )

    return rectangles


def point_in_any_rect(x, y, rectangles, margin):
    """Return whether a point lies in any rectangle with the given margin."""
    return any(
        point_in_rect(x, y, x0, y0, width, height, margin)
        for x0, y0, width, height in rectangles
    )


def validate_denser_options(denser_near, denser_min_ratio):
    if denser_near < 0:
        raise ValueError(f"denser_near must be >= 0, got {denser_near}")
    if denser_min_ratio < 1:
        raise ValueError(
            f"denser_min_ratio must be >= 1, got {denser_min_ratio}"
        )


def generate_evpoints(
    group,
    H,
    delta,
    d,
    min_spacing,
    denser_near=default_denser_near,
    denser_min_ratio=default_denser_min_ratio,
):
    """Generate xy evaluation points for one configuration.

    A regular grid covers the full vacuum domain.  When denser_near is
    positive, a second grid with spacing min_spacing / denser_min_ratio is
    added in the vacuum-side padding around removed solid regions.
    """
    del group  # Reserved for future group-specific geometry.
    validate_denser_options(denser_near, denser_min_ratio)

    x_max = 6 * (W + delta)
    y_max = H + d
    rectangles = solid_rectangles(H, delta, d)

    x_coords = np.arange(min_spacing / 2 - outexpand, x_max + outexpand, min_spacing)
    y_coords = np.arange(min_spacing / 2, y_max, min_spacing)

    vacuum_points = []

    for x in x_coords:
        for y in y_coords:
            if not point_in_any_rect(x, y, rectangles, min_edge_dist):
                vacuum_points.append((x, y))

    if denser_near == 0 or denser_min_ratio == 1:
        return vacuum_points

    dense_spacing = min_spacing / denser_min_ratio
    # Use the coarse-grid origin so integer refinement ratios retain all
    # coarse points and add regularly interleaved samples.
    dense_x_coords = np.arange(
        min_spacing / 2 - outexpand,
        x_max + outexpand,
        dense_spacing,
    )
    dense_y_coords = np.arange(min_spacing / 2, y_max, dense_spacing)
    dense_points = []
    near_margin = min_edge_dist + denser_near

    for x in dense_x_coords:
        for y in dense_y_coords:
            if point_in_any_rect(x, y, rectangles, min_edge_dist):
                continue
            if point_in_any_rect(x, y, rectangles, near_margin):
                dense_points.append((x, y))

    # Rounded keys remove points shared by nested coarse and dense grids while
    # preserving stable output ordering.
    points_by_key = {
        (round(float(x), 12), round(float(y), 12)): (x, y)
        for x, y in vacuum_points
    }
    for x, y in dense_points:
        points_by_key.setdefault(
            (round(float(x), 12), round(float(y), 12)),
            (x, y),
        )

    return sorted(points_by_key.values(), key=lambda point: (point[0], point[1]))


def generate_evpoints_adaptive_for_z(
    group,
    H,
    delta,
    d,
    z,
    max_points,
    denser_near=default_denser_near,
    denser_min_ratio=default_denser_min_ratio,
):
    """Adapt spacing until the combined regular and dense grid fits the limit."""
    min_spacing = base_min_spacing

    for i in range(max_adaptive_iter):
        points = generate_evpoints(
            group,
            H,
            delta,
            d,
            min_spacing,
            denser_near=denser_near,
            denser_min_ratio=denser_min_ratio,
        )
        n_points = len(points)

        if n_points <= max_points:
            return points, min_spacing, i

        min_spacing *= spacing_growth

    raise RuntimeError(
        f"Adaptive spacing failed: group={group}, d={d}, z={z}, "
        f"last spacing={min_spacing:.6f}, points={n_points}, max_points={max_points}"
    )


def write_evpoints(
    filename,
    z_point_data,
    group,
    H,
    delta,
    d,
    denser_near,
    denser_min_ratio,
):
    """Write xyz evaluation points to txt file."""
    with open(filename, "w") as f:
        f.write(f"# Evaluation points for group {group}, d = {d}\n")
        f.write(f"# W = {W}, H = {H}, delta = {delta}\n")
        f.write(f"# p = {2 * (W + delta)}\n")
        f.write(
            f"# min_edge_dist = {min_edge_dist}, denser_near = {denser_near}, "
            f"denser_min_ratio = {denser_min_ratio}\n"
        )
        f.write("# Format: x y z\n")

        for item in z_point_data:
            z = item["z"]
            points = item["points"]
            spacing = item["spacing"]
            max_points = item["max_points"]
            dense_spacing = (
                f"{spacing / denser_min_ratio:.6f}"
                if denser_near > 0 and denser_min_ratio > 1
                else "disabled"
            )

            f.write(
                f"# z = {z}, xy_points = {len(points)}, "
                f"max_points = {max_points}, adaptive_spacing = {spacing:.6f}, "
                f"dense_spacing = {dense_spacing}\n"
            )

            for x, y in points:
                f.write(f"{x:.4f} {y:.4f} {z:.4f}\n")


def parse_float_list(raw):
    """Parse a comma-separated float list."""
    if raw is None:
        return None
    values = []
    for item in raw.split(","):
        item = item.strip()
        if item:
            values.append(float(item))
    return values


def plot_evpoints(
    plot_filename,
    z_point_data,
    group,
    H,
    delta,
    d,
    plot_z=None,
    point_size=1.0,
    hide_legend=False,
    dpi=300,
    plot_formats=None,
):
    """Plot xy distribution of evaluation points for all z layers."""
    import matplotlib.pyplot as plt

    p = 2 * (W + delta)
    plot_z_values = None if plot_z is None else {round(float(z), 12) for z in plot_z}

    fig, ax = plt.subplots(figsize=(8, 5))
    plotted_points = 0

    for item in z_point_data:
        z = item["z"]
        points = item["points"]
        if plot_z_values is not None and round(float(z), 12) not in plot_z_values:
            continue

        if len(points) == 0:
            continue

        arr = np.array(points)
        ax.scatter(arr[:, 0], arr[:, 1], s=point_size, label=f"z={z}, n={len(points)}")
        plotted_points += len(points)

    # Draw upper fins
    for n in range(rep):
        rects = [
            (n * p + W / 4 - W / 4, d, W / 2, H),
            ((n + 1) * p - W / 4 - W / 4, d, W / 2, H),
        ]

        for rx, ry, rw, rh in rects:
            ax.add_patch(
                plt.Rectangle(
                    (rx, ry),
                    rw,
                    rh,
                    fill=False,
                    linewidth=1.0,
                )
            )

    # Draw lower fins
    for n in range(rep):
        rx = (2 * n + 1) * p / 2 - W / 2
        ry = 0

        ax.add_patch(
            plt.Rectangle(
                (rx, ry),
                W,
                H,
                fill=False,
                linewidth=1.0,
            )
        )

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(
        f"group {group}, d={d}, H={H}, delta={delta}, plotted points={plotted_points}"
    )
    ax.grid(True, linewidth=0.3)
    if not hide_legend:
        ax.legend(markerscale=3, fontsize=8)

    fig.tight_layout()
    plot_path = Path(plot_filename)
    if plot_formats is None:
        formats = [plot_path.suffix.lower().lstrip(".") or "png"]
    else:
        formats = [str(item).lower().lstrip(".") for item in plot_formats]
    formats = list(dict.fromkeys(formats))
    unsupported = [item for item in formats if item not in {"png", "svg", "pdf"}]
    if unsupported:
        raise ValueError(f"unsupported plot format(s): {', '.join(unsupported)}")
    plot_stem = plot_path.stem if plot_path.suffix else plot_path.name
    written = []
    for plot_format in formats:
        output_path = plot_path.parent / f"{plot_stem}.{plot_format}"
        save_kwargs = {"bbox_inches": "tight"}
        if plot_format == "png":
            save_kwargs["dpi"] = dpi
        fig.savefig(output_path, format=plot_format, **save_kwargs)
        written.append(output_path)
    plt.close(fig)
    return written


def main():
    from datetime import datetime
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--with-plot",
        action="store_true",
        help="Generate scatter plots of x-y evaluation-point distribution.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("evpoints_plots"),
        help="Output root directory for plots (default: evpoints_plots/).",
    )
    parser.add_argument(
        "--batch-name",
        type=str,
        default=None,
        help="Batch name for plot subdirectory (default: auto-timestamp).",
    )
    parser.add_argument(
        "--denser-near",
        type=float,
        default=default_denser_near,
        help=(
            "Width of the denser vacuum-side padding around removed solids "
            "(default: 0.2; use 0 to disable)."
        ),
    )
    parser.add_argument(
        "--denser-min-ratio",
        type=float,
        default=default_denser_min_ratio,
        help=(
            "Minimum linear refinement ratio in the denser padding; dense "
            "spacing is adaptive_spacing / ratio (default: 2)."
        ),
    )
    parser.add_argument(
        "--plot-z",
        default=None,
        help="Comma-separated z values to include in plots only, e.g. 0.75.",
    )
    parser.add_argument(
        "--plot-point-size",
        type=float,
        default=1.0,
        help="Scatter marker size for evpoint plots.",
    )
    parser.add_argument(
        "--hide-legend",
        action="store_true",
        help="Do not draw legends on evpoint plots.",
    )
    parser.add_argument(
        "--skip-write-evpoints",
        action="store_true",
        help="Generate plots without writing evpoints/*.txt files.",
    )
    args = parser.parse_args()

    validate_denser_options(args.denser_near, args.denser_min_ratio)
    plot_z = parse_float_list(args.plot_z)

    point_limits = [max_points_per_z, *max_points_per_z_overrides.values()]
    if any(len(limits) != len(z_values) for limits in point_limits):
        raise ValueError(
            "z_values and every max-points-per-z list must have the same length."
        )

    os.makedirs("evpoints", exist_ok=True)

    plot_out_dir = None
    if args.with_plot:
        batch_name = args.batch_name or datetime.now().strftime("%Y%m%d_%H%M%S")
        plot_out_dir = args.out_dir.resolve() / batch_name / "all" / "figures" / "evpoints"
        plot_out_dir.mkdir(parents=True, exist_ok=True)

    count = 0

    for group, params in groups.items():
        H = params["H"]
        delta = params["delta"]

        for d in params["displacements"]:
            z_point_data = []
            config_max_points = max_points_per_z_overrides.get(
                (group, float(d)),
                max_points_per_z,
            )

            for z, max_points in zip(z_values, config_max_points):
                points, used_spacing, n_iter = generate_evpoints_adaptive_for_z(
                    group=group,
                    H=H,
                    delta=delta,
                    d=d,
                    z=z,
                    max_points=max_points,
                    denser_near=args.denser_near,
                    denser_min_ratio=args.denser_min_ratio,
                )

                z_point_data.append(
                    {
                        "z": z,
                        "points": points,
                        "spacing": used_spacing,
                        "adaptive_iter": n_iter,
                        "max_points": max_points,
                    }
                )

            txt_filename = f"evpoints/{group}_d{d:.1f}.txt"
            if not args.skip_write_evpoints:
                write_evpoints(
                    txt_filename,
                    z_point_data,
                    group,
                    H,
                    delta,
                    d,
                    denser_near=args.denser_near,
                    denser_min_ratio=args.denser_min_ratio,
                )

            if args.with_plot and plot_out_dir is not None:
                plot_filename = str(plot_out_dir / f"{group}_d{d:.1f}.png")
                plot_evpoints(
                    plot_filename,
                    z_point_data,
                    group,
                    H,
                    delta,
                    d,
                    plot_z=plot_z,
                    point_size=args.plot_point_size,
                    hide_legend=args.hide_legend,
                )

            total_points = sum(len(item["points"]) for item in z_point_data)

            print(f"{txt_filename}: total points = {total_points}")

            for item in z_point_data:
                print(
                    f"  z={item['z']}: "
                    f"{len(item['points'])} points, "
                    f"max={item['max_points']}, "
                    f"spacing={item['spacing']:.4f}, "
                    f"adaptive_iter={item['adaptive_iter']}"
                )

            if args.with_plot:
                print(f"  plot saved to {plot_out_dir / f'{group}_d{d:.1f}.png'}")

            count += 1

    print(f"\nTotal: {count} evpoint files generated")


if __name__ == "__main__":
    main()
