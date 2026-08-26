"""
Post-process scuff-neq results: check NaN validity, generate recompute tasks,
merge valid data, and plot PAbs/PRad vs omega.

Usage:
    python postproc_scuffem.py <results_dir> [--last] [--freq-per-subset=N]
"""
from pathlib import Path
import argparse
import re
import shutil
from datetime import datetime
from collections import defaultdict


ROOT = Path(__file__).resolve().parents[2]
MATERIAL_DIR = ROOT / "inputs" / "materials"

RUN_HEADER_RE = re.compile(r"^#\s*scuff-neq run on\s+")
SUBTASK_DIR_RE = re.compile(r"^subtask-(.+)-(\d+)$")

SCUFFGEO_TEMPLATE = """\
OBJECT Lower
  MESHFILE mshs/{mesh_low}
  MATERIAL FILE_{material_file}
ENDOBJECT

OBJECT Upper
  MESHFILE mshs/{mesh_up}
  MATERIAL FILE_{material_file}
  DISPLACED 0 {displacement} 0
ENDOBJECT
"""

SBATCH_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err

mkdir -p logs

module load scuff-em

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_DYNAMIC=false

scuff-neq \\
  --geometry task.scuffgeo \\
  --epfile evpoints \\
  --omegafile omegas \\
  --emtpft \\
  --sourceobject Lower \\
  --destobject Upper
"""

SUBMIT_ALL_TEMPLATE = """\
#!/bin/bash
set -euo pipefail

for dir in recompute-*; do
  [ -d "$dir" ] || continue
  (
    cd "$dir"
    sbatch run.sbatch
  )
done
"""


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Post-process scuff-neq results: check NaN, generate recompute tasks, merge valid data."
    )
    parser.add_argument(
        "results_dir",
        type=Path,
        help="Path to batch results directory (e.g., results-batches/batch02-refined)",
    )
    parser.add_argument(
        "--last",
        action="store_true",
        help="Only analyze content after the LAST '# scuff-neq run on' header.",
    )
    parser.add_argument(
        "--freq-per-subset",
        type=int,
        default=2,
        help="Max frequencies per recompute subset (default: 2).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: <results_dir>/../postproc_<timestamp>).",
    )
    parser.add_argument(
        "--batch-name",
        type=str,
        default=None,
        help="Batch directory name (e.g., '20260630_143000'). Auto-generated if omitted.",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def write_text_file(path: Path, content: str) -> None:
    path.write_text(content.rstrip() + "\n", encoding="utf-8", newline="\n")


def is_nan_value(s: str) -> bool:
    return s.strip().lower() in ("nan", "-nan", "+nan")


# ---------------------------------------------------------------------------
# Subtask discovery
# ---------------------------------------------------------------------------

def parse_subtask_name(dirname: str) -> tuple[str, int] | None:
    m = SUBTASK_DIR_RE.match(dirname)
    if m:
        return m.group(1), int(m.group(2))
    return None


def discover_subtasks(results_dir: Path) -> list[Path]:
    subtasks = []
    for p in sorted(results_dir.iterdir()):
        if p.is_dir() and p.name.startswith("subtask-"):
            if parse_subtask_name(p.name) is not None:
                subtasks.append(p)
    return subtasks


# ---------------------------------------------------------------------------
# File parsing with --last support
# ---------------------------------------------------------------------------

def split_runs(lines: list[str]) -> list[list[str]]:
    runs: list[list[str]] = []
    current_run: list[str] = []
    for line in lines:
        if RUN_HEADER_RE.match(line):
            if current_run:
                runs.append(current_run)
            current_run = [line]
        else:
            current_run.append(line)
    if current_run:
        runs.append(current_run)
    return runs


def extract_data_lines(lines: list[str], last_only: bool) -> list[str]:
    if last_only:
        runs = split_runs(lines)
        if not runs:
            return []
        target_lines = runs[-1]
    else:
        target_lines = lines
    data_lines = []
    for line in target_lines:
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            data_lines.append(stripped)
    return data_lines


# ---------------------------------------------------------------------------
# SIFlux parsing & validity
# ---------------------------------------------------------------------------

def check_siflux_validity(data_lines: list[str]) -> tuple[list[dict], list[dict]]:
    """SIFlux: invalid if PAbs (col 3) OR PRad (col 4) is NaN."""
    valid, invalid = [], []
    for line in data_lines:
        parts = line.split()
        if len(parts) < 5:
            continue
        omega = float(parts[1])
        row = {"omega": omega, "source_dest": parts[2], "values": parts}
        if is_nan_value(parts[3]) or is_nan_value(parts[4]):
            invalid.append(row)
        else:
            valid.append(row)
    return valid, invalid


# ---------------------------------------------------------------------------
# SRFlux parsing & validity
# ---------------------------------------------------------------------------

def check_srflux_validity(data_lines: list[str]) -> tuple[list[dict], list[dict]]:
    """SRFlux: invalid if ALL of cols 6-17 (Px..Tzz) are NaN."""
    valid, invalid = [], []
    for line in data_lines:
        parts = line.split()
        if len(parts) < 18:
            continue
        omega = float(parts[1])
        flux_cols = parts[6:18]
        row = {"omega": omega, "values": parts}
        if all(is_nan_value(c) for c in flux_cols):
            invalid.append(row)
        else:
            valid.append(row)
    return valid, invalid


# ---------------------------------------------------------------------------
# Per-subtask processing
# ---------------------------------------------------------------------------

def process_subtask(subtask_dir: Path, last_only: bool) -> dict:
    name = subtask_dir.name
    config_id, subset_idx = parse_subtask_name(name)
    result = {
        "subtask_name": name,
        "config_id": config_id,
        "subset_idx": subset_idx,
        "siflux_valid": [],
        "siflux_invalid": [],
        "srflux_valid": [],
        "srflux_invalid": [],
    }
    siflux_path = subtask_dir / "task.SIFlux.EMTPFT"
    if siflux_path.exists():
        lines = siflux_path.read_text(encoding="utf-8").splitlines()
        data_lines = extract_data_lines(lines, last_only)
        result["siflux_valid"], result["siflux_invalid"] = check_siflux_validity(data_lines)
    srflux_path = subtask_dir / "task.SRFlux"
    if srflux_path.exists():
        lines = srflux_path.read_text(encoding="utf-8").splitlines()
        data_lines = extract_data_lines(lines, last_only)
        result["srflux_valid"], result["srflux_invalid"] = check_srflux_validity(data_lines)
    return result


# ---------------------------------------------------------------------------
# Invalid points collection & output
# ---------------------------------------------------------------------------

def collect_invalid_points(results: list[dict]) -> set[tuple[str, float]]:
    pts: set[tuple[str, float]] = set()
    for r in results:
        for row in r["siflux_invalid"]:
            pts.add((r["subtask_name"], row["omega"]))
        for row in r["srflux_invalid"]:
            pts.add((r["subtask_name"], row["omega"]))
    return pts


def write_invalid_points(
    invalid_points: set[tuple[str, float]], batch_dir: Path, case_ids: list[str]
) -> None:
    sorted_pts = sorted(invalid_points, key=lambda x: (x[0], x[1]))

    # Aggregate by config for per-case metadata
    by_config: dict[str, list] = defaultdict(list)
    for name, omega in sorted_pts:
        cid, _ = parse_subtask_name(name)
        by_config[cid].append((name, omega))

    for cid in case_ids:
        metadata_dir = batch_dir / cid / "metadata"
        metadata_dir.mkdir(parents=True, exist_ok=True)

        # invalid_points.txt per case
        case_pts = sorted(by_config.get(cid, []), key=lambda x: x[1])
        lines = [f"{name}  {omega:.6e}" for name, omega in case_pts]
        write_text_file(metadata_dir / "invalid_points.txt", "\n".join(lines))

        # summary.txt per case
        pairs = case_pts
        omegas = sorted(set(o for _, o in pairs))
        summary = [
            "Invalid points summary",
            f"Total invalid (subtask, omega) pairs: {len(sorted_pts)}",
            "",
            f"Config {cid}: {len(pairs)} invalid pairs, {len(omegas)} unique omegas",
        ]
        for o in omegas:
            summary.append(f"  omega = {o:.6e}")
        write_text_file(metadata_dir / "summary.txt", "\n".join(summary))


# ---------------------------------------------------------------------------
# Config info extraction from existing scuffgeo
# ---------------------------------------------------------------------------

def extract_config_info(results_dir: Path, subtask_name: str) -> dict:
    scuffgeo_path = results_dir / subtask_name / "task.scuffgeo"
    if not scuffgeo_path.exists():
        raise FileNotFoundError(f"task.scuffgeo not found: {scuffgeo_path}")
    content = scuffgeo_path.read_text(encoding="utf-8")
    mesh_files = re.findall(r"MESHFILE\s+mshs/(\S+)", content)
    displaced_match = re.search(r"DISPLACED\s+0\s+(\S+)\s+0", content)
    if len(mesh_files) < 2 or not displaced_match:
        raise ValueError(f"Could not parse scuffgeo: {scuffgeo_path}")
    return {
        "mesh_low": mesh_files[0],
        "mesh_up": mesh_files[1],
        "displacement": displaced_match.group(1),
    }


# ---------------------------------------------------------------------------
# Material file selection
# ---------------------------------------------------------------------------

def select_material_file() -> tuple[str, str]:
    """Return (material_content, display_name).

    Priority: 2k+ > 53p > 27p.  Raises FileNotFoundError if none.
    """
    # 1. 2k+ point file
    high_res = MATERIAL_DIR / "aligned_SiO2-Franta-300C.txt"
    if high_res.exists():
        content = high_res.read_text(encoding="utf-8")
        n_data = sum(1 for ln in content.splitlines() if ln.strip() and not ln.strip().startswith("#"))
        if n_data >= 2000:
            return content, f"silica-{n_data}p.dat"

    # 2. 53-point (generate from real+img)
    real_f = ROOT / "aligned-SiO2_Franta-300C-sparse_53 - real.txt"
    img_f = ROOT / "aligned-SiO2_Franta-300C-sparse_53 - img.txt"
    if real_f.exists() and img_f.exists():
        omegas, reals, imags = [], [], []
        for line in real_f.read_text(encoding="utf-8").splitlines():
            p = line.strip().split()
            if len(p) >= 2:
                omegas.append(float(p[0]))
                reals.append(float(p[1]))
        for line in img_f.read_text(encoding="utf-8").splitlines():
            p = line.strip().split()
            if len(p) >= 2:
                imags.append(float(p[1]))
        out = ["# Silica (Franta 300C) - 53 points", "# omega(rad/s)  Re(eps)  Im(eps)"]
        for w, re_val, im_val in zip(omegas, reals, imags):
            out.append(f"{w:.6e}  {re_val:.6e}  {im_val:.6e}")
        return "\n".join(out), "silica-53p.dat"

    # 3. 27-point
    p27 = MATERIAL_DIR / "silica-27p.txt"
    if not p27.exists():
        p27 = MATERIAL_DIR / "silica-27p.dat"
    if p27.exists():
        return p27.read_text(encoding="utf-8"), "silica-27p.dat"

    raise FileNotFoundError("No suitable material file found.")


# ---------------------------------------------------------------------------
# Recompute task generation
# ---------------------------------------------------------------------------

def group_invalid_by_config(invalid_points: set[tuple[str, float]]) -> dict[str, set[float]]:
    by_cfg: dict[str, set[float]] = defaultdict(set)
    for name, omega in invalid_points:
        cid, _ = parse_subtask_name(name)
        by_cfg[cid].add(omega)
    return dict(by_cfg)


def split_into_subsets(omegas: list[float], max_n: int) -> list[list[float]]:
    omegas_sorted = sorted(omegas)
    return [omegas_sorted[i:i + max_n] for i in range(0, len(omegas_sorted), max_n)]


def generate_recompute_tasks(
    invalid_points: set[tuple[str, float]],
    results_dir: Path,
    batch_dir: Path,
    freq_per_subset: int,
    case_ids: list[str],
) -> list[Path]:
    by_config = group_invalid_by_config(invalid_points)
    created_dirs: list[Path] = []

    mat_content, mat_name = select_material_file()

    for config_id in case_ids:
        recompute_dir = batch_dir / config_id / "recompute"
        recompute_dir.mkdir(parents=True, exist_ok=True)

        # Write material file at recompute level for reference
        write_text_file(recompute_dir / mat_name, mat_content)

        if config_id not in by_config:
            continue

        omegas = sorted(by_config[config_id])
        subsets = split_into_subsets(omegas, freq_per_subset)

        # Pick a representative original subtask
        representative = next(
            name for name, _ in invalid_points if parse_subtask_name(name)[0] == config_id
        )
        config_info = extract_config_info(results_dir, representative)

        for subset_idx, omega_subset in enumerate(subsets, start=1):
            # Use "recompute-" prefix to distinguish from original subtask-* dirs
            rc_name = f"recompute-{config_id}-{subset_idx}"
            rc_dir = recompute_dir / rc_name
            rc_dir.mkdir(parents=True, exist_ok=True)

            # omegas
            write_text_file(rc_dir / "omegas", "\n".join(f"{w}" for w in omega_subset))

            # material file (copy into subtask dir, matching existing pattern)
            write_text_file(rc_dir / mat_name, mat_content)

            # task.scuffgeo
            write_text_file(
                rc_dir / "task.scuffgeo",
                SCUFFGEO_TEMPLATE.format(
                    mesh_low=config_info["mesh_low"],
                    mesh_up=config_info["mesh_up"],
                    material_file=mat_name,
                    displacement=config_info["displacement"],
                ),
            )

            # evpoints
            evpoints_src = results_dir / representative / "evpoints"
            if evpoints_src.exists():
                shutil.copyfile(evpoints_src, rc_dir / "evpoints")

            # mshs/
            mshs_src = results_dir / representative / "mshs"
            if mshs_src.exists():
                mshs_dst = rc_dir / "mshs"
                if mshs_dst.exists():
                    shutil.rmtree(mshs_dst)
                shutil.copytree(mshs_src, mshs_dst)

            # run.sbatch
            write_text_file(
                rc_dir / "run.sbatch",
                SBATCH_TEMPLATE.format(job_name=f"r{config_id}{subset_idx}"),
            )

            created_dirs.append(rc_dir)

        write_text_file(recompute_dir / "submit_all.sh", SUBMIT_ALL_TEMPLATE)
    return created_dirs


# ---------------------------------------------------------------------------
# Valid data merging
# ---------------------------------------------------------------------------

SIFLUX_HEADER = "subtask_name\tconfig_id\tomega\tsource_dest\tPAbs\tPRad\tXForce\tYForce\tZForce\tXTorque\tYTorque\tZTorque"


def _siflux_entry_line(entry: dict) -> str:
    return (
        f"{entry['subtask_name']}\t{entry['config_id']}\t{entry['omega']:.6e}\t"
        f"{entry['source_dest']}\t{entry['pabs']}\t{entry['prad']}\t"
        f"{entry['xforce']}\t{entry['yforce']}\t{entry['zforce']}\t"
        f"{entry['xtorque']}\t{entry['ytorque']}\t{entry['ztorque']}"
    )


def merge_valid_siflux(
    results: list[dict], batch_dir: Path, case_ids: list[str]
) -> dict[str, list[dict]]:
    all_valid: list[dict] = []
    by_config: dict[str, list[dict]] = defaultdict(list)

    for r in results:
        cid = r["config_id"]
        for row in r["siflux_valid"]:
            entry = {
                "subtask_name": r["subtask_name"],
                "config_id": cid,
                "omega": row["omega"],
                "source_dest": row["source_dest"],
                "pabs": row["values"][3],
                "prad": row["values"][4],
                "xforce": row["values"][5],
                "yforce": row["values"][6],
                "zforce": row["values"][7],
                "xtorque": row["values"][8],
                "ytorque": row["values"][9],
                "ztorque": row["values"][10],
            }
            all_valid.append(entry)
            by_config[cid].append(entry)

    # per-case data/ directory
    for cid in case_ids:
        data_dir = batch_dir / cid / "data"
        data_dir.mkdir(parents=True, exist_ok=True)
        entries = by_config.get(cid, [])
        lines = [SIFLUX_HEADER] + [
            _siflux_entry_line(e)
            for e in sorted(entries, key=lambda e: (e["omega"], e["source_dest"]))
        ]
        write_text_file(data_dir / f"{cid}_siflux.tsv", "\n".join(lines))

    # merged_siflux.tsv at batch level
    lines = [SIFLUX_HEADER] + [
        _siflux_entry_line(e)
        for e in sorted(all_valid, key=lambda e: (e["config_id"], e["omega"], e["source_dest"]))
    ]
    write_text_file(batch_dir / "merged_siflux.tsv", "\n".join(lines))

    return dict(by_config)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def generate_plots(
    by_config: dict[str, list[dict]], batch_dir: Path, case_ids: list[str]
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for cid in case_ids:
        fig_dir = batch_dir / cid / "figures" / "siflux"
        fig_dir.mkdir(parents=True, exist_ok=True)

        entries = by_config.get(cid, [])
        by_sd: dict[str, list[dict]] = defaultdict(list)
        for e in entries:
            by_sd[e["source_dest"]].append(e)

        for quantity, col_key in [("PAbs", "pabs"), ("PRad", "prad")]:
            fig, ax = plt.subplots(figsize=(8, 5))
            for sd in sorted(by_sd):
                rows = sorted(by_sd[sd], key=lambda e: e["omega"])
                omegas = [e["omega"] for e in rows]
                vals = []
                for e in rows:
                    try:
                        vals.append(float(e[col_key]))
                    except ValueError:
                        vals.append(float("nan"))
                ax.plot(omegas, vals, "o-", label=f"sd={sd}", markersize=3)
            ax.set_xlabel("omega (rad/s)")
            ax.set_ylabel(quantity)
            ax.set_title(f"Config {cid}: {quantity} vs omega")
            ax.legend(fontsize=8)
            ax.grid(True, linewidth=0.3)
            fig.tight_layout()
            fig.savefig(fig_dir / f"{cid}_{quantity}.png", dpi=150)
            plt.close(fig)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def write_report(
    results: list[dict],
    invalid_points: set[tuple[str, float]],
    created_dirs: list[Path],
    batch_dir: Path,
    args: argparse.Namespace,
) -> None:
    lines = [
        "Post-processing report",
        "=" * 60,
        f"Results directory: {args.results_dir}",
        f"--last: {args.last}",
        f"--freq-per-subset: {args.freq_per_subset}",
        f"Timestamp: {datetime.now().isoformat()}",
        "",
        "Subtask summary:",
        "-" * 40,
    ]
    tsiv, tsii, tsrv, tsri = 0, 0, 0, 0
    for r in results:
        nv, ni = len(r["siflux_valid"]), len(r["siflux_invalid"])
        sv, si = len(r["srflux_valid"]), len(r["srflux_invalid"])
        tsiv += nv; tsii += ni; tsrv += sv; tsri += si
        lines.append(f"  {r['subtask_name']}: SIFlux valid={nv} invalid={ni}, SRFlux valid={sv} invalid={si}")

    lines += [
        "",
        "Totals:",
        f"  SIFlux: {tsiv} valid, {tsii} invalid",
        f"  SRFlux: {tsrv} valid, {tsri} invalid",
        f"  Unique invalid (subtask, omega) pairs: {len(invalid_points)}",
        "",
        f"Recompute tasks generated: {len(created_dirs)} directories",
    ]
    if created_dirs:
        lines.append("Recompute directories:")
        for d in created_dirs:
            lines.append(f"  {d.name}")

    write_text_file(batch_dir / "report.txt", "\n".join(lines))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main_from_args(args: argparse.Namespace) -> None:
    results_dir = args.results_dir.resolve()
    if not results_dir.is_dir():
        print(f"Error: results directory not found: {results_dir}")
        return

    if args.out_dir is not None:
        output_dir = args.out_dir.resolve()
    else:
        output_dir = results_dir.parent / f"postproc_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # batch_dir: {out_dir}/{batch_name}/ (or timestamp if not specified)
    if args.batch_name:
        batch_dir = output_dir / args.batch_name
    else:
        batch_dir = output_dir / datetime.now().strftime("%Y%m%d_%H%M%S")
    batch_dir.mkdir(parents=True, exist_ok=True)

    subtasks = discover_subtasks(results_dir)
    if not subtasks:
        print(f"No subtask directories found in {results_dir}")
        return
    print(f"Found {len(subtasks)} subtask directories")

    results = []
    for sd in subtasks:
        print(f"Processing {sd.name}...")
        results.append(process_subtask(sd, last_only=args.last))

    invalid_points = collect_invalid_points(results)
    print(f"Found {len(invalid_points)} invalid (subtask, omega) pairs")

    # Determine case_ids from results for per-case directory structure
    case_ids = sorted(set(r["config_id"] for r in results))

    write_invalid_points(invalid_points, batch_dir, case_ids)

    created_dirs: list[Path] = []
    if invalid_points:
        created_dirs = generate_recompute_tasks(
            invalid_points, results_dir, batch_dir, args.freq_per_subset, case_ids
        )
        print(f"Generated {len(created_dirs)} recompute subtask directories")

    if any(r["siflux_valid"] for r in results):
        by_config = merge_valid_siflux(results, batch_dir, case_ids)
        generate_plots(by_config, batch_dir, case_ids)
        print("Generated valid data merge and plots")
    else:
        print("No valid SIFlux data found; skipping merge and plots")

    write_report(results, invalid_points, created_dirs, batch_dir, args)
    print(f"\nOutput written to: {batch_dir}")
    print(f"Report: {batch_dir / 'report.txt'}")


def main() -> None:
    main_from_args(parse_args())


if __name__ == "__main__":
    main()
