# Poynting Flux Post-Processing Workspace

This workspace contains SCUFF-EM setup files, raw `results-batches`, and
post-processing scripts for spatially resolved Poynting-vector analysis.

Start here for results-batch processing:

- [README_postproc.md](README_postproc.md): ordered post-processing workflow
- [README_subtasks.md](README_subtasks.md): subtask layout and naming rules
- [README_scuffneq.md](README_scuffneq.md): SCUFF-EM run notes

## Main Workflow

Use `batch_postprocess.py` for new batch results:

```bash
python batch_postprocess.py run-python-flux \
  --task-dir results-batches/ord3-results \
  --out-root results-batches/ord3_processed \
  --batch-name 20260708_python_integrated_flux \
  --planes 0.75,1.0,1.5 \
  --color-by absSxy,Sz \
  --colormap rainbow \
  --color-scale log \
  --stride 1
```

This runs the ordered stages:

1. combine per-frequency subtasks into `merged_srflux.tsv`
2. integrate `Sx/Sy/Sz` over frequency in Python
3. generate `rainbow`, log-scaled flux diagrams with fixed-scale x1-density quivers
4. run SIFlux validity/recompute analysis

Outputs are written under `{out-root}/{batch-name}/`; raw input directories are
not modified.

## Input and output directory management

The post-processing tools accept input and output locations explicitly, so raw
SCUFF-EM task directories can be kept separate from derived results. The main
batch workflow uses `--task-dir` for the input task directory and `--out-root`
for the output root:

```bash
python batch_postprocess.py run-python-flux \
  --task-dir D:/scuffem/tasks/ord3-results \
  --out-root D:/analysis/poy \
  --batch-name ord3_case1
```

Other utilities use the same convention with `--out-dir` or `--out` when they
produce a single dataset or figure set. Check each subcommand's `--help` output
for its exact input and output options. Generated results should be stored
under a dedicated output root and are excluded from Git by `.gitignore`.

## Merge existing outputs and visualize them

Collected outputs from earlier batches can be merged into a new directory. The
inputs are read-only; duplicate rows are removed using all TSV columns, and a
`merge_manifest.yaml` records the sources:

```bash
python batch_postprocess.py merge \
  --input-dir results-batches/batch_a \
  --input-dir results-batches/batch_b \
  --out-dir results-batches/merged_cases
```

The merged directory keeps the normal `<config>/data/merged_srflux.tsv` layout,
so it can be passed directly to the existing integration/visualization stage
by giving it a batch directory and a batch name:

```bash
python batch_postprocess.py python-integrated-flux \
  --out-root results-batches \
  --batch-name merged_cases \
  --configs config_name \
  --out-dir results-batches/merged_cases/config_name/figures/integrated_flux
```

For a single collected dataset, the lower-level interface is explicit about
both locations:

```bash
python plot_integrated_flux.py \
  --dataset results-batches/merged_cases/config_name/data \
  --out results-batches/merged_cases/config_name/figures/integrated_flux
```

Scripts tied to one geometry, evpoint layout, or diagnostic result are under
`scripts/specialized/`. The root-level Python modules are the reusable batch
processing and plotting entry points.
