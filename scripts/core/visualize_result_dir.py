"""Visualize one collected poy result directory with one required argument."""
from __future__ import annotations

import argparse
from pathlib import Path

try:
    from .plot_integrated_flux import parse_figure_formats, run_integrated_flux
except ImportError:  # direct execution: python scripts/core/visualize_result_dir.py
    from plot_integrated_flux import parse_figure_formats, run_integrated_flux


def find_datasets(result_dir: Path) -> list[Path]:
    """Return collected dataset directories below *result_dir*.

    The accepted layouts are ``result_dir/merged_srflux.tsv``,
    ``result_dir/data/merged_srflux.tsv``, or a batch directory containing
    one or more ``*/data/merged_srflux.tsv`` files.
    """
    candidates: list[Path] = []
    if (result_dir / "merged_srflux.tsv").is_file():
        candidates.append(result_dir)
    if (result_dir / "data" / "merged_srflux.tsv").is_file():
        candidates.append(result_dir / "data")
    candidates.extend(
        path.parent
        for path in sorted(result_dir.glob("*/data/merged_srflux.tsv"))
    )
    candidates.extend(
        path.parent
        for path in sorted(result_dir.glob("*/*/data/merged_srflux.tsv"))
    )
    return list(dict.fromkeys(candidates))


def parse_csv(raw: str) -> list[str]:
    return [item.strip() for item in raw.split(",") if item.strip()]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Visualize collected poy results. Only --result-dir is required."
    )
    parser.add_argument(
        "--result-dir", required=True,
        help="Collected result directory, a dataset directory, or a batch directory.",
    )
    parser.add_argument(
        "--out-dir", default=None,
        help="Figure output directory. Default: <result-dir>/figures/integrated_flux.",
    )
    parser.add_argument("--batch-name", default=None)
    parser.add_argument("--planes", default="0.75,1.0,1.5")
    parser.add_argument("--color-by", default="absSxy,Sz")
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--quiver-scale", type=float, default=None)
    parser.add_argument("--no-quiver", dest="include_quiver", action="store_false")
    parser.add_argument("--color-scale", choices=["linear", "log"], default="log")
    parser.add_argument("--colormap", default="rainbow")
    parser.add_argument("--figure-format", default="png")
    parser.set_defaults(include_quiver=True)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    result_dir = Path(args.result_dir).resolve()
    if not result_dir.is_dir():
        raise SystemExit(f"result directory not found: {result_dir}")

    datasets = find_datasets(result_dir)
    if not datasets:
        raise SystemExit(
            f"no collected dataset found below {result_dir}; "
            "expected merged_srflux.tsv"
        )

    requested_out = Path(args.out_dir).resolve() if args.out_dir else None
    planes = [float(item) for item in parse_csv(args.planes)]
    color_by = parse_csv(args.color_by)
    formats = parse_figure_formats(args.figure_format)

    for dataset in datasets:
        if requested_out is None:
            out_dir = dataset / "figures" / "integrated_flux"
        elif len(datasets) == 1:
            out_dir = requested_out
        else:
            out_dir = requested_out / dataset.parent.name

        result = run_integrated_flux(
            dataset_dir=dataset,
            out_dir=out_dir,
            batch_name=args.batch_name,
            planes=planes,
            color_by_list=color_by,
            stride=args.stride,
            quiver_scale=args.quiver_scale,
            include_quiver=args.include_quiver,
            color_scale=args.color_scale,
            colormap=args.colormap,
            figure_format=", ".join(formats),
        )
        print(f"{dataset}: wrote {result['n_figures']} figures to {out_dir}")


if __name__ == "__main__":
    main()
