# Poynting Flux Post-Processing Workspace

This workspace contains SCUFF-EM setup files, raw `results-batches`, and
post-processing scripts for spatially resolved Poynting-vector analysis.

Start here for results-batch processing:

- [docs/README_postproc.md](docs/README_postproc.md): ordered post-processing workflow
- [docs/README_subtasks.md](docs/README_subtasks.md): subtask layout and naming rules
- [docs/README_scuffneq.md](docs/README_scuffneq.md): SCUFF-EM run notes

## Main Workflow

Use `scripts/core/batch_postprocess.py` for new batch results:

```bash
python -m scripts.core.batch_postprocess run-python-flux \
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
python -m scripts.core.batch_postprocess run-python-flux \
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
python -m scripts.core.batch_postprocess merge \
  --input-dir results-batches/batch_a \
  --input-dir results-batches/batch_b \
  --out-dir results-batches/merged_cases
```

The merged directory keeps the normal `<config>/data/merged_srflux.tsv` layout,
so it can be passed directly to the existing integration/visualization stage
by giving it a batch directory and a batch name:

```bash
python -m scripts.core.batch_postprocess python-integrated-flux \
  --out-root results-batches \
  --batch-name merged_cases \
  --configs config_name \
  --out-dir results-batches/merged_cases/config_name/figures/integrated_flux
```

To visualize a prepared result directory, only the result directory is
required. The script automatically recognizes a dataset directory containing
`merged_srflux.tsv`, a `data/` directory, or a batch directory containing
multiple config datasets. Figures are written below the result directory by
default:

```bash
python scripts/core/visualize_result_dir.py \
  --result-dir results-batches/merged_cases/config_name/data
```

Optional plotting arguments are forwarded to the existing integration and
plotting implementation, for example:

```bash
python scripts/core/visualize_result_dir.py \
  --result-dir results-batches/merged_cases/config_name/data \
  --out-dir figures/my_case \
  --planes 0.75,1.0,1.5 \
  --color-by absSxy,Sz \
  --no-quiver
```

The equivalent module form is `python -m scripts.core.visualize_result_dir`.

Reusable processing modules are under `scripts/core/`; scripts tied to one
geometry, evpoint layout, or diagnostic result are under `scripts/specialized/`.
Input notes are under `docs/`. The root contains only project metadata, the
main README, and material files needed by legacy SCUFF-EM setup workflows.
