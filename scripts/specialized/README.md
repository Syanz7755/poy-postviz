# Specialized poy scripts

These scripts are retained for historical or case-specific SCUFF-EM inputs,
evpoint layouts, and figure assembly. They are not part of the general batch
pipeline. `scuffem_setup/` contains geometry/material/task generators; the
other files contain case-specific diagnostics and figure assembly. Run them
from the project root when they refer to root-relative materials or result
directories, and inspect each script's `--help` output before use.

General workflows remain in `batch_postprocess.py`, `postproc_scuffem.py`,
`scuff_pv_post.py`, and `plot_integrated_flux.py`.
