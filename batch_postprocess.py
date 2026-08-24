"""Ordered batch post-processing for SCUFF-EM results-batches directories.

The module entry points are:

    collect_srflux_batch(...)
    python_integrated_flux_batch(...)
    siflux_batch(...)
    run_python_flux_pipeline(...)

The CLI wraps the same functions so manual and scripted use follow the same
ordered stages.
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

import postproc_scuffem
import scuff_pv_post
from plot_integrated_flux import run_integrated_flux


@dataclass(frozen=True)
class BatchLayout:
    task_dir: Path
    out_root: Path
    batch_name: str

    @property
    def batch_dir(self) -> Path:
        return self.out_root / self.batch_name

    def config_data_dir(self, config: str) -> Path:
        return self.batch_dir / config / "data"

    def integrated_flux_dir(self, config: str) -> Path:
        return self.batch_dir / config / "figures" / "integrated_flux"

    @property
    def siflux_out_dir(self) -> Path:
        return self.batch_dir / "siflux"


def parse_csv_floats(value: str) -> list[float]:
    return [float(part.strip()) for part in value.split(",") if part.strip()]


def parse_csv_strings(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def make_layout(task_dir: str | Path, out_root: str | Path, batch_name: str | None = None) -> BatchLayout:
    task_path = Path(task_dir).resolve()
    out_path = Path(out_root).resolve()
    name = batch_name or datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    return BatchLayout(task_dir=task_path, out_root=out_path, batch_name=name)


def discover_configs(task_dir: str | Path) -> list[str]:
    task_path = Path(task_dir).resolve()
    if not task_path.is_dir():
        raise FileNotFoundError(f"task directory not found: {task_path}")
    configs = scuff_pv_post.discover_config_names(task_path)
    if not configs:
        raise RuntimeError(f"no subtask-* directories found in {task_path}")
    return configs


def write_batch_manifest(layout: BatchLayout, stages: list[dict[str, Any]]) -> Path:
    layout.batch_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "task_dir": str(layout.task_dir),
        "out_root": str(layout.out_root),
        "batch_name": layout.batch_name,
        "stages": stages,
    }
    path = layout.batch_dir / "batch_postprocess_manifest.yaml"
    path.write_text(yaml.dump(manifest, default_flow_style=False, sort_keys=False), encoding="utf-8")
    return path


def collect_srflux_config(layout: BatchLayout, config: str, last: bool = True) -> Path:
    data_dir = layout.config_data_dir(config)
    data_dir.mkdir(parents=True, exist_ok=True)
    args = argparse.Namespace(
        task_dir=str(layout.task_dir),
        group_name=config,
        filebase="task",
        out=str(data_dir),
        coord_tol=1e-9,
        duplicate_policy="error",
        last=last,
    )
    scuff_pv_post.cmd_collect(args)
    return data_dir


def collect_srflux_batch(
    task_dir: str | Path,
    out_root: str | Path,
    batch_name: str | None = None,
    configs: list[str] | None = None,
    last: bool = True,
) -> dict[str, Any]:
    layout = make_layout(task_dir, out_root, batch_name)
    config_names = configs or discover_configs(layout.task_dir)
    data_dirs = {config: collect_srflux_config(layout, config, last=last) for config in config_names}
    manifest = write_batch_manifest(
        layout,
        [{"name": "collect_srflux", "configs": config_names, "last": last}],
    )
    return {"layout": layout, "configs": config_names, "data_dirs": data_dirs, "manifest": manifest}


def python_integrated_flux_batch(
    out_root: str | Path,
    batch_name: str,
    configs: list[str] | None = None,
    planes: list[float] | None = None,
    color_by_list: list[str] | None = None,
    stride: int = 1,
    quiver_scale: float | None = None,
    include_quiver: bool = True,
    color_scale: str = "log",
    colormap: str = "rainbow",
    figure_format: str = "png",
) -> dict[str, Any]:
    out_path = Path(out_root).resolve()
    batch_dir = out_path / batch_name
    if not batch_dir.is_dir():
        raise FileNotFoundError(f"batch directory not found: {batch_dir}")

    config_names = configs or sorted(
        p.name for p in batch_dir.iterdir()
        if p.is_dir() and (p / "data" / "merged_srflux.tsv").exists()
    )
    if not config_names:
        raise RuntimeError(f"no collected config data found in {batch_dir}")

    results = {}
    for config in config_names:
        dataset_dir = batch_dir / config / "data"
        fig_dir = batch_dir / config / "figures" / "integrated_flux"
        results[config] = run_integrated_flux(
            dataset_dir=dataset_dir,
            out_dir=fig_dir,
            batch_name=batch_name,
            planes=planes,
            color_by_list=color_by_list,
            stride=stride,
            quiver_scale=quiver_scale,
            include_quiver=include_quiver,
            color_scale=color_scale,
            colormap=colormap,
            figure_format=figure_format,
        )
    return {"batch_dir": batch_dir, "configs": config_names, "results": results}


def siflux_batch(
    task_dir: str | Path,
    out_root: str | Path,
    batch_name: str,
    last: bool = True,
    freq_per_subset: int = 2,
) -> dict[str, Any]:
    layout = make_layout(task_dir, out_root, batch_name)
    siflux_dir = layout.siflux_out_dir
    siflux_dir.mkdir(parents=True, exist_ok=True)

    args = argparse.Namespace(
        results_dir=layout.task_dir,
        out_dir=siflux_dir,
        batch_name="analysis",
        last=last,
        freq_per_subset=freq_per_subset,
    )
    postproc_scuffem.main_from_args(args)
    return {"out_dir": siflux_dir / "analysis"}


def run_python_flux_pipeline(
    task_dir: str | Path,
    out_root: str | Path,
    batch_name: str | None = None,
    configs: list[str] | None = None,
    last: bool = True,
    planes: list[float] | None = None,
    color_by_list: list[str] | None = None,
    stride: int = 1,
    quiver_scale: float | None = None,
    include_quiver: bool = True,
    color_scale: str = "log",
    colormap: str = "rainbow",
    figure_format: str = "png",
    include_siflux: bool = True,
    freq_per_subset: int = 2,
) -> dict[str, Any]:
    layout = make_layout(task_dir, out_root, batch_name)
    config_names = configs or discover_configs(layout.task_dir)

    stages: list[dict[str, Any]] = []
    collect_result = collect_srflux_batch(
        task_dir=layout.task_dir,
        out_root=layout.out_root,
        batch_name=layout.batch_name,
        configs=config_names,
        last=last,
    )
    stages.append({"name": "collect_srflux", "configs": config_names, "last": last})

    flux_result = python_integrated_flux_batch(
        out_root=layout.out_root,
        batch_name=layout.batch_name,
        configs=config_names,
        planes=planes,
        color_by_list=color_by_list,
        stride=stride,
        quiver_scale=quiver_scale,
        include_quiver=include_quiver,
        color_scale=color_scale,
        colormap=colormap,
        figure_format=figure_format,
    )
    stages.append({
        "name": "python_integrated_flux",
        "configs": config_names,
        "planes": planes or [0.75, 1.0, 1.5],
        "color_by": color_by_list or ["absSxy", "Sz"],
        "stride": stride,
        "quiver_scale": quiver_scale,
        "include_quiver": include_quiver,
        "color_scale": color_scale,
        "colormap": colormap,
        "figure_format": figure_format,
    })

    siflux_result = None
    if include_siflux:
        siflux_result = siflux_batch(
            task_dir=layout.task_dir,
            out_root=layout.out_root,
            batch_name=layout.batch_name,
            last=last,
            freq_per_subset=freq_per_subset,
        )
        stages.append({"name": "siflux", "last": last, "freq_per_subset": freq_per_subset})

    manifest = write_batch_manifest(layout, stages)
    return {
        "layout": layout,
        "configs": config_names,
        "collect": collect_result,
        "python_integrated_flux": flux_result,
        "siflux": siflux_result,
        "manifest": manifest,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Ordered results-batches post-processing.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_common(p: argparse.ArgumentParser) -> None:
        p.add_argument("--task-dir", required=True, help="Input directory containing subtask-{config}-{index} folders.")
        p.add_argument("--out-root", required=True, help="Output root; batch_name is created under this directory.")
        p.add_argument("--batch-name", default=None, help="Stable output batch name. Defaults to UTC timestamp.")
        p.add_argument("--configs", default=None, help="Comma-separated configs. Default: auto-discover.")
        p.add_argument("--last", action="store_true", default=True, help="Use only the last scuff-neq run block.")
        p.add_argument("--all-runs", dest="last", action="store_false", help="Process all run blocks.")

    p_run = subparsers.add_parser("run-python-flux", help="collect SRFlux, integrate over omega in Python, plot flux diagrams.")
    add_common(p_run)
    p_run.add_argument("--planes", default="0.75,1.0,1.5", help="Comma-separated z planes.")
    p_run.add_argument("--color-by", default="absSxy,Sz", help="Comma-separated scalar fields.")
    p_run.add_argument("--stride", type=int, default=1, help="Quiver stride; 1 means x1 density.")
    p_run.add_argument("--quiver-scale", type=float, default=None, help="Fixed quiver scale; auto-derived per config if omitted.")
    p_run.add_argument("--no-quiver", dest="include_quiver", action="store_false", help="Draw coloring only.")
    p_run.add_argument("--color-scale", choices=["linear", "log"], default="log", help="Color normalization.")
    p_run.add_argument("--colormap", default="rainbow", help="Matplotlib colormap name.")
    p_run.add_argument("--figure-format", default="png", help="Figure format(s), comma-separated. Default: png.")
    p_run.set_defaults(include_quiver=True)
    p_run.add_argument("--skip-siflux", action="store_true", help="Skip SIFlux validity/recompute analysis.")
    p_run.add_argument("--freq-per-subset", type=int, default=2, help="SIFlux recompute frequency subset size.")

    p_collect = subparsers.add_parser("collect", help="Only collect SRFlux datasets.")
    add_common(p_collect)

    p_flux = subparsers.add_parser("python-integrated-flux", help="Only run Python integration + plotting on an existing batch.")
    p_flux.add_argument("--out-root", required=True)
    p_flux.add_argument("--batch-name", required=True)
    p_flux.add_argument("--configs", default=None)
    p_flux.add_argument("--planes", default="0.75,1.0,1.5")
    p_flux.add_argument("--color-by", default="absSxy,Sz")
    p_flux.add_argument("--stride", type=int, default=1)
    p_flux.add_argument("--quiver-scale", type=float, default=None)
    p_flux.add_argument("--no-quiver", dest="include_quiver", action="store_false", help="Draw coloring only.")
    p_flux.add_argument("--color-scale", choices=["linear", "log"], default="log", help="Color normalization.")
    p_flux.add_argument("--colormap", default="rainbow", help="Matplotlib colormap name.")
    p_flux.add_argument("--figure-format", default="png", help="Figure format(s), comma-separated. Default: png.")
    p_flux.set_defaults(include_quiver=True)

    p_siflux = subparsers.add_parser("siflux", help="Only run SIFlux validity/recompute analysis.")
    add_common(p_siflux)
    p_siflux.add_argument("--freq-per-subset", type=int, default=2)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    configs = parse_csv_strings(args.configs) if getattr(args, "configs", None) else None

    if args.command == "run-python-flux":
        result = run_python_flux_pipeline(
            task_dir=args.task_dir,
            out_root=args.out_root,
            batch_name=args.batch_name,
            configs=configs,
            last=args.last,
            planes=parse_csv_floats(args.planes),
            color_by_list=parse_csv_strings(args.color_by),
            stride=args.stride,
            quiver_scale=args.quiver_scale,
            include_quiver=args.include_quiver,
            color_scale=args.color_scale,
            colormap=args.colormap,
            figure_format=args.figure_format,
            include_siflux=not args.skip_siflux,
            freq_per_subset=args.freq_per_subset,
        )
        print(f"Batch complete: {result['layout'].batch_dir}")
        print(f"Manifest: {result['manifest']}")
    elif args.command == "collect":
        result = collect_srflux_batch(args.task_dir, args.out_root, args.batch_name, configs, args.last)
        print(f"Collected configs: {', '.join(result['configs'])}")
        print(f"Manifest: {result['manifest']}")
    elif args.command == "python-integrated-flux":
        result = python_integrated_flux_batch(
            out_root=args.out_root,
            batch_name=args.batch_name,
            configs=configs,
            planes=parse_csv_floats(args.planes),
            color_by_list=parse_csv_strings(args.color_by),
            stride=args.stride,
            quiver_scale=args.quiver_scale,
            include_quiver=args.include_quiver,
            color_scale=args.color_scale,
            colormap=args.colormap,
            figure_format=args.figure_format,
        )
        print(f"Python integrated flux complete: {result['batch_dir']}")
    elif args.command == "siflux":
        if args.batch_name is None:
            print("ERROR: --batch-name is required for siflux output placement", file=sys.stderr)
            sys.exit(1)
        result = siflux_batch(args.task_dir, args.out_root, args.batch_name, args.last, args.freq_per_subset)
        print(f"SIFlux output: {result['out_dir']}")
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
