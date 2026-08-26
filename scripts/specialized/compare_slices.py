"""Compare slice exports for agreement analysis."""
import argparse
import sys
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

BASE = Path(__file__).parent / "slice_runs"

SLICES = {
    "y=0.01":       BASE / "ad1_y0.01",
    "y=d-0.01":     BASE / "ad1_y_d-0.01",
    "x=(W+d)/2":    BASE / "ad1_x_half",
    "x=(W+d)*1.5":  BASE / "ad1_x_1.5p",
    "y=H+0.01":     BASE / "ad1_y_H+0.01",
    "y=H+d/2":      BASE / "ad1_y_H+d_2",
}


def load_weighted(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path / "slice_weighted.tsv", sep="\t")
    for c in ["x", "y", "z", "Sx", "Sy", "Sz"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load_meta(path: Path) -> dict:
    with open(path / "slice_meta.yaml", encoding="utf-8") as f:
        return yaml.safe_load(f)


def norm_xy(Sx, Sy):
    return np.sqrt(Sx**2 + Sy**2)


def compare_series(x1, v1, x2, v2, label1, label2, qty):
    """Compare two 1D series on potentially different x grids."""
    # Interpolate onto common grid
    x_common = np.union1d(x1, x2)
    x_common.sort()

    v1_interp = np.interp(x_common, x1, v1, left=np.nan, right=np.nan)
    v2_interp = np.interp(x_common, x2, v2, left=np.nan, right=np.nan)

    mask = np.isfinite(v1_interp) & np.isfinite(v2_interp)
    if mask.sum() < 2:
        return {"status": "insufficient_overlap", "n_common": int(mask.sum())}

    v1c = v1_interp[mask]
    v2c = v2_interp[mask]

    diff = v1c - v2c
    abs_diff = np.abs(diff)
    max_val = max(np.max(np.abs(v1c)), np.max(np.abs(v2c)), 1e-30)

    return {
        "status": "ok",
        "n_common": int(mask.sum()),
        "max_abs_diff": float(np.max(abs_diff)),
        "mean_abs_diff": float(np.mean(abs_diff)),
        "max_rel_diff": float(np.max(abs_diff / np.maximum(np.abs(v1c), 1e-30))),
        "mean_rel_diff": float(np.mean(abs_diff / np.maximum(np.abs(v1c), 1e-30))),
        "norm_diff_vs_peak": float(np.max(abs_diff) / max_val),
        "correlation": float(np.corrcoef(v1c, v2c)[0, 1]) if len(v1c) > 1 and np.std(v1c) > 0 and np.std(v2c) > 0 else np.nan,
    }


def main():
    parser = argparse.ArgumentParser(description="Compare slice exports for agreement analysis.")
    parser.add_argument("--batch-name", type=str, default=None,
                        help="Batch name for output subdirectory (default: auto-timestamp).")
    parser.add_argument("--out", type=Path, default=None,
                        help="Output root directory (default: slice_runs/).")
    args = parser.parse_args()

    out_root = args.out.resolve() if args.out else BASE
    batch_name = args.batch_name or datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = out_root / batch_name / "all" / "figures" / "compare"
    out_dir.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Load all slices
    data = {}
    for name, path in SLICES.items():
        if not path.exists():
            print(f"SKIP: {name} — {path} not found")
            continue
        try:
            df = load_weighted(path)
            meta = load_meta(path)
            data[name] = {"df": df, "meta": meta}
            print(f"Loaded {name}: {len(df)} points, axis={meta['slice_axis']}")
        except Exception as e:
            print(f"SKIP: {name} — {e}")

    if len(data) < 2:
        print("Not enough slices to compare.")
        return

    # Determine perpendicular axis for each slice
    for name, d in data.items():
        meta = d["meta"]
        d["perp_axis"] = "x" if meta["slice_axis"] == "y" else "y"
        d["perp_col"] = d["perp_axis"]

    # Define comparison pairs with physical meaning
    pairs = [
        ("y=0.01", "y=d-0.01", "Lower gap: y=0.01 vs y=d-0.01"),
        ("x=(W+d)/2", "x=(W+d)*1.5", "Fin cross-section: x=center vs x=gap"),
        ("y=H+0.01", "y=H+d/2", "Upper region: y=H+0.01 vs y=H+d/2"),
        ("y=0.01", "y=H+0.01", "Across structure: lower gap vs upper gap"),
    ]

    quantities = ["Sx", "Sy", "absSxy"]
    report_lines = []

    report_lines.append("=" * 70)
    report_lines.append("SLICE AGREEMENT REPORT — ad1 (W=1, H=5, delta=0.6, d=1, p=3.2)")
    report_lines.append("=" * 70)
    report_lines.append("")

    for s1_name, s2_name, desc in pairs:
        if s1_name not in data or s2_name not in data:
            report_lines.append(f"--- {desc} ---")
            report_lines.append(f"  SKIPPED (data not available)")
            report_lines.append("")
            continue

        d1 = data[s1_name]
        d2 = data[s2_name]
        df1 = d1["df"]
        df2 = d2["df"]
        meta1 = d1["meta"]
        meta2 = d2["meta"]

        report_lines.append(f"--- {desc} ---")
        report_lines.append(f"  {s1_name}: axis={meta1['slice_axis']}, value={meta1['folded_value']}, {len(df1)} pts")
        report_lines.append(f"  {s2_name}: axis={meta2['slice_axis']}, value={meta2['folded_value']}, {len(df2)} pts")

        # Both must have same slice axis for direct comparison
        if meta1["slice_axis"] != meta2["slice_axis"]:
            report_lines.append(f"  NOTE: different slice axes ({meta1['slice_axis']} vs {meta2['slice_axis']}), comparing as cross-sections")
            # For cross-axis comparison, compare magnitude distributions
            for qty in quantities:
                v1 = df1[qty].values if qty in df1.columns else norm_xy(df1["Sx"].values, df1["Sy"].values)
                v2 = df2[qty].values if qty in df2.columns else norm_xy(df2["Sx"].values, df2["Sy"].values)
                report_lines.append(f"  {qty}: {s1_name} mean={np.nanmean(np.abs(v1)):.4e}, {s2_name} mean={np.nanmean(np.abs(v2)):.4e}")
            report_lines.append("")
            continue

        # Same axis: compare along perpendicular coordinate
        perp = d1["perp_axis"]
        x1 = df1[perp].values
        x2 = df2[perp].values

        for qty in quantities:
            if qty == "absSxy":
                v1 = norm_xy(df1["Sx"].values, df1["Sy"].values)
                v2 = norm_xy(df2["Sx"].values, df2["Sy"].values)
            else:
                v1 = df1[qty].values
                v2 = df2[qty].values

            result = compare_series(x1, v1, x2, v2, s1_name, s2_name, qty)

            if result["status"] == "ok":
                report_lines.append(
                    f"  {qty}: corr={result['correlation']:.4f}, "
                    f"max_rel_diff={result['max_rel_diff']:.2e}, "
                    f"norm_diff_vs_peak={result['norm_diff_vs_peak']:.2e}, "
                    f"common_pts={result['n_common']}"
                )
            else:
                report_lines.append(f"  {qty}: {result['status']}")

            # Overlay plot
            fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True,
                                     gridspec_kw={"height_ratios": [3, 1]})

            sort1 = np.argsort(x1)
            sort2 = np.argsort(x2)

            axes[0].plot(x1[sort1], v1[sort1], "o-", ms=3, label=s1_name, alpha=0.8)
            axes[0].plot(x2[sort2], v2[sort2], "s--", ms=3, label=s2_name, alpha=0.8)
            axes[0].set_ylabel(qty)
            axes[0].set_title(f"{qty}: {s1_name} vs {s2_name}")
            axes[0].legend(fontsize=8)

            # Difference on common grid
            x_common = np.union1d(x1, x2)
            x_common.sort()
            v1i = np.interp(x_common, x1, v1)
            v2i = np.interp(x_common, x2, v2)
            diff = v1i - v2i
            axes[1].plot(x_common, diff, "k.-", ms=2)
            axes[1].axhline(0, color="gray", ls=":", lw=0.5)
            axes[1].set_ylabel(f"diff ({s1_name} - {s2_name})")
            axes[1].set_xlabel(perp)

            fig.tight_layout()
            def safe_name(s):
                return s.replace("=","").replace("+","p").replace("/","_").replace("(","").replace(")","").replace("*","x")
            fname = f"compare_{safe_name(s1_name)}_{safe_name(s2_name)}_{qty}.png"
            fig.savefig(out_dir / fname, dpi=150)
            plt.close(fig)

        report_lines.append("")

    # Summary table
    report_lines.append("=" * 70)
    report_lines.append("SUMMARY")
    report_lines.append("=" * 70)
    report_lines.append("")
    report_lines.append("Geometry: W=1.0, H=5, delta=0.6, d=1.0, period=3.2")
    report_lines.append("Config: ad1 (group a), z=0.75, T=300K")
    report_lines.append("")
    report_lines.append("Expected agreements:")
    report_lines.append("  y=0.01 vs y=d-0.01: mirror symmetry across gap center (y=d/2=0.5)")
    report_lines.append("  x=(W+d)/2 vs x=(W+d)*1.5: fin-center vs gap-center cross-sections")
    report_lines.append("  y=H+0.01 vs y=H+d/2: both above fins, different heights")
    report_lines.append("  y=0.01 vs y=H+0.01: lower gap vs upper gap (structure-dependent)")

    report_text = "\n".join(report_lines)
    (out_dir / "agreement_report.txt").write_text(report_text, encoding="utf-8")
    print(report_text)
    print(f"\nReport saved to {out_dir / 'agreement_report.txt'}")
    print(f"Overlay plots in {out_dir}")


if __name__ == "__main__":
    main()
