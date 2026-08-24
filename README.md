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
