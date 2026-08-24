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
