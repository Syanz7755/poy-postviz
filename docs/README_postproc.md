# SCUFF-EM Results-Batch Post-Processing

This directory contains a small ordered pipeline for SCUFF-EM `scuff-neq`
results stored as:

```text
results-batches/<raw-batch>/
  subtask-<config>-1/
    task.SRFlux
    task.SIFlux.EMTPFT
    task.scuffgeo
    evpoints
    omegas
  subtask-<config>-2/
  ...
```

Use `batch_postprocess.py` for normal batch work. The older scripts remain
available as lower-level tools.

## Recommended Command

For a new raw batch, run the ordered Python-integrated flux pipeline:

```bash
python batch_postprocess.py run-python-flux \
  --task-dir results-batches/ord3-results \
  --out-root results-batches/ord3_processed \
  --batch-name 20260708_python_integrated_flux \
  --planes 0.75,1.0,1.5 \
  --color-by absSxy,Sz \
  --stride 1
```

Default behavior:

1. Discover configs from `subtask-{config}-{index}` directories.
2. Collect each config's `task.SRFlux` subtasks into canonical TSV files.
3. Integrate `Sx/Sy/Sz` over omega in Python with `numpy.trapezoid`.
4. Plot integrated flux diagrams with `rainbow`, log-scaled coloring, fixed-scale quivers, and 1x quiver density.
5. Run SIFlux validity/recompute analysis unless `--skip-siflux` is passed.

Use `--no-quiver` when you want coloring-only figures while debugging cellmap
geometry.

The input batch is read-only. Outputs are written under:

```text
{out-root}/{batch-name}/
```

## Ordered Modules

### `batch_postprocess.py`

Top-level orchestration module. It exposes callable functions and a CLI.

Useful functions:

```python
from batch_postprocess import (
    collect_srflux_batch,
    python_integrated_flux_batch,
    siflux_batch,
    run_python_flux_pipeline,
)
```

CLI stages:

```bash
# Full ordered workflow
python batch_postprocess.py run-python-flux --task-dir ... --out-root ...

# Stage 1 only: combine SRFlux subtasks
python batch_postprocess.py collect --task-dir ... --out-root ... --batch-name ...

# Stage 2 only: Python integrate + plot an already collected batch
python batch_postprocess.py python-integrated-flux --out-root ... --batch-name ...

# Stage 3 only: SIFlux validity/recompute analysis
python batch_postprocess.py siflux --task-dir ... --out-root ... --batch-name ...
```

### `scuff_pv_post.py`

Lower-level SRFlux tool. It still provides:

- `collect`: merge raw `task.SRFlux` files into `merged_srflux.tsv`
- `plot`: legacy summed/integrated cellmap plots
- `integrate`: wrapper around `scuff-integrate`
- `slice`: 1D cross-section extraction and plots
- `run-all`: older all-in-one workflow

Use this directly when you need legacy plotting, `scuff-integrate`, or slices.

### `plot_integrated_flux.py`

Python integration and flux-diagram module. It expects a collected dataset
containing `merged_srflux.tsv`.

It now exposes:

```python
from plot_integrated_flux import run_integrated_flux
```

Direct CLI use:

```bash
python plot_integrated_flux.py \
  --dataset results-batches/ord3_processed/20260708_python_integrated_flux/cd1/data \
  --out results-batches/ord3_processed/20260708_python_integrated_flux/cd1/figures/integrated_flux \
  --batch-name 20260708_python_integrated_flux \
  --planes 0.75,1.0,1.5 \
  --color-by absSxy,Sz \
  --colormap rainbow \
  --color-scale log \
  --stride 1 \
  --no-quiver
```

### `postproc_scuffem.py`

SIFlux validity and recompute-task analysis. It can still be run directly:

```bash
python postproc_scuffem.py results-batches/ord3-results \
  --last \
  --out-dir results-batches/ord3_processed/20260708_python_integrated_flux/siflux \
  --batch-name analysis
```

It also has `main_from_args(args)` so `batch_postprocess.py` can call it as a
module.

## Output Layout

For `batch_postprocess.py run-python-flux`, outputs look like:

```text
results-batches/ord3_processed/
  20260708_python_integrated_flux/
    batch_postprocess_manifest.yaml
    cd1/
      data/
        merged_srflux.tsv
        merged.SRFlux
        points.tsv
        omegas.tsv
        manifest.yaml
        raw_headers/
      figures/
        integrated_flux/
          cd1_poynting_integrated_python.tsv
          integrated_flux_manifest.yaml
          20260708_python_integrated_flux_cd1_z0.750_absSxy_integrated_rainbow_log_quiver1x_fixedscale.png
          20260708_python_integrated_flux_cd1_z0.750_absSxy_integrated_rainbow_log_coloronly.png
          20260708_python_integrated_flux_cd1_z0.750_Sz_integrated_rainbow_log_quiver1x_fixedscale.png
          ...
    cd4/
      data/
      figures/integrated_flux/
    siflux/
      analysis/
        report.txt
        merged_siflux.tsv
        <config>/figures/siflux/
        <config>/metadata/
        <config>/recompute/
```

## Important Options

| Option | Applies to | Default | Meaning |
|---|---|---:|---|
| `--last` | collect/SIFlux | on in `batch_postprocess.py` | Use only the last `# scuff-neq run on` block. |
| `--all-runs` | collect/SIFlux | off | Process every run block. |
| `--configs cd1,cd4` | batch stages | auto | Limit processing to selected configs. |
| `--planes 0.75,1.0,1.5` | plotting | `0.75,1.0,1.5` | z planes to plot. |
| `--color-by absSxy,Sz` | plotting | `absSxy,Sz` | scalar fields for color maps. |
| `--color-scale log` | plotting | `log` | Color normalization. Signed fields use symmetric log scaling. |
| `--colormap rainbow` | plotting | `rainbow` | Matplotlib colormap for coloring. |
| `--stride 1` | plotting | `1` | Quiver downsampling; `1` means x1 density. |
| `--quiver-scale <float>` | plotting | auto | Fixed Matplotlib quiver scale. Omit for one auto scale per config. |
| `--no-quiver` | plotting | off | Draw coloring-only figures with `_coloronly.png` suffix. |
| `--skip-siflux` | full pipeline | off | Skip SIFlux validity/recompute analysis. |

## Data Semantics

- `merged_srflux.tsv` is frequency-resolved SRFlux after subtask combination.
- `*_poynting_integrated_python.tsv` is numerically integrated over omega.
- Integration uses trapezoidal integration of finite values per spatial point.
- Points whose full frequency history is NaN remain NaN and plot as blank cells.
- Quiver arrows use integrated `Sx,Sy`.
- `rainbow` is the default colormap for the Python-integrated flux figures.

## Cellmap Coloring on Nonuniform Evaluation Points

The plotting code supports nonuniform rectilinear evaluation grids, such as
fine point spacing in narrow gaps and coarser spacing elsewhere.

Layering:

- `draw_coloring_layer(...)` draws only the scalar cellmap.
- `draw_quiver_layer(...)` draws only the vector arrows.
- `plot_integrated(...)` composes those layers. Pass `include_quiver=False`
  or CLI `--no-quiver` to suppress the quiver layer.

Coloring rule:

- High-occupancy rectilinear planes use `pcolormesh` with nonuniform midpoint
  cell edges.
- Low-occupancy sparse planes use thresholded axis-aligned rectangles.
- The base near-neighbor threshold is `0.1`.
- For coarser planes, the x/y thresholds adapt from same-row/same-column gap
  distributions. This lets sparse regions use larger cells without treating
  every sparse point as an island.
- A cell extends to the midpoint between same-row/same-column neighbors only
  when that x or y gap is below the adaptive threshold for that axis.
- Larger x/y gaps stay blank, with rectangular blank-region boundaries.
- No interpolation or smoothing is used for cellmap coloring.

This behavior is implemented in both:

- `plot_integrated_flux.py` for Python-integrated flux diagrams
- `scuff_pv_post.py plot` for legacy SRFlux cellmaps

For sparse planes, manifests record `cellmap_method: threshold_rectangles`,
`adaptive_cell_base_distance_threshold: 0.1`, and the resolved
`adaptive_cell_x_distance_threshold` / `adaptive_cell_y_distance_threshold`.
For dense rectilinear planes, manifests record
`cellmap_method: rectilinear_pcolormesh` or the legacy `cellmap` mode.

## Legacy One-Line Workflow

The previous all-in-one workflow still works:

```bash
python scuff_pv_post.py run-all \
  --task-dir results-batches/ord \
  --out-dir results-batches/ord_processed \
  --last \
  --only-central
```

Prefer `batch_postprocess.py run-python-flux` when you want the ordered,
module-friendly workflow and Python frequency integration.
