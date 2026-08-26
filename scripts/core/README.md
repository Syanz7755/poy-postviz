# Core workflow modules

These modules implement the supported reusable workflow:

- `batch_postprocess.py`: collect, merge, integrate, and orchestrate results;
- `plot_integrated_flux.py`: integrate a collected dataset and create figures;
- `postproc_scuffem.py`: validate SCUFF-EM output and prepare recomputations;
- `scuff_pv_post.py`: collect spatial-flux datasets.

Run them from the project root with `python -m scripts.core.<module>` so the
package imports and project-relative paths remain stable.
