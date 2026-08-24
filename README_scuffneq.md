# scuff-neq Input Files

## Directory Structure

```
poy/
├── materials/
│   └── silica-27p.txt          # 27-point silica dielectric function
├── mshs/                       # Place your .msh files here
│   ├── a3-low.msh, a3-up.msh
│   ├── b3-low.msh, b3-up.msh
│   ├── c3-low.msh, c3-up.msh
│   ├── d3-low.msh, d3-up.msh
│   ├── e3-low.msh, e3-up.msh
│   └── f3-low.msh, f3-up.msh
├── scuffgeo/                   # Generated .scuffgeo files
│   ├── a_d1.0.scuffgeo ... a_d10.0.scuffgeo
│   ├── b_d1.0.scuffgeo ... b_d10.0.scuffgeo
│   ├── c_d1.0.scuffgeo ... c_d10.0.scuffgeo
│   ├── d_d0.5.scuffgeo ... d_d5.0.scuffgeo
│   ├── e_d0.5.scuffgeo ... e_d5.0.scuffgeo
│   └── f_d0.5.scuffgeo ... f_d5.0.scuffgeo
├── evpoints/                   # Generated evaluation point files
│   ├── a_d1.0.txt ... a_d10.0.txt
│   ├── b_d1.0.txt ... b_d10.0.txt
│   ├── c_d1.0.txt ... c_d10.0.txt
│   ├── d_d0.5.txt ... d_d5.0.txt
│   ├── e_d0.5.txt ... e_d5.0.txt
│   └── f_d0.5.txt ... f_d5.0.txt
├── gen_material.py             # Script to generate material file
├── gen_scuffgeo.py             # Script to generate .scuffgeo files
└── gen_evpoints.py             # Script to generate evaluation points
```

## Geometry Parameters

| Group | H   | delta | p = 2(W+delta) | Displacements        |
|-------|-----|-------|----------------|----------------------|
| a     | 5   | 0.6   | 3.2            | 1.0, 2.0, 4.0, 5.0, 10.0 |
| b     | 5   | 0.2   | 2.4            | 1.0, 2.0, 4.0, 5.0, 10.0 |
| c     | 5   | 0.1   | 2.2            | 1.0, 2.0, 4.0, 5.0, 10.0 |
| d     | 2   | 0.6   | 3.2            | 0.5, 1.0, 2.0, 5.0       |
| e     | 2   | 0.2   | 2.4            | 0.5, 1.0, 2.0, 5.0       |
| f     | 2   | 0.1   | 2.2            | 0.5, 1.0, 2.0, 5.0       |

- W = 1.0 (all groups)
- Evaluation planes: z = 0.75, 1.0, 1.5
- Base grid spacing: 0.08; adaptive growth factor: 1.05
- Solid-edge exclusion margin: 0.01
- Default near-solid refinement: width 0.2 at 2x linear resolution
- Default maximum points per z plane: 8000, 800, 800
- C1 maximum points per z plane: 13000, 800, 800
- C4 maximum points per z plane: 11000, 800, 800

The maximum is a hard upper bound, not a target. The adaptive generator counts
both base-grid and near-solid points, then increases both spacings together
until the combined point count fits the applicable limit.

## Usage

### 1. Place mesh files
Copy your .msh files to `./mshs/` with names: `{group}3-low.msh`, `{group}3-up.msh`

### 2. Run scuff-neq
```bash
# Example for group a, displacement 1.0
scuff-neq --geometry scuffgeo/a_d1.0.scuffgeo \
          --EPFile evpoints/a_d1.0.txt \
          --OmegaFile omega.txt \
          --power \
          --OMIT_SELFTERMS

# Batch run for all configurations
for group in a b c d e f; do
  for d in $(ls scuffgeo/${group}_d*.scuffgeo | sed 's/.*_d//;s/\.scuffgeo//'); do
    scuff-neq --geometry scuffgeo/${group}_d${d}.scuffgeo \
              --EPFile evpoints/${group}_d${d}.txt \
              --OmegaFile omega.txt \
              --power \
              --OMIT_SELFTERMS \
              --OutFile results/${group}_d${d}.out
  done
done
```

### 3. Key options
- `--OMIT_SELFTERMS`: Reduces computation (omit self-term contributions)
- `--power`: Compute Poynting vector
- `--EPFile`: Evaluation points file
- `--OmegaFile`: Frequency list (create your own)

## Regenerating Files

To regenerate with different parameters, edit the Python scripts and run:
```bash
python gen_material.py
python gen_scuffgeo.py
python gen_evpoints.py --denser-near 0.2 --denser-min-ratio 2
```

`--denser-near 0` disables near-solid refinement. The dense spacing is the
adaptive base spacing divided by `--denser-min-ratio`.
