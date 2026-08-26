#!/bin/bash
#SBATCH --job-name=snq_%g_%d
#SBATCH --output=logs/%g_d%j.out
#SBATCH --error=logs/%g_d%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=8G

# Usage: sbatch --export=g=a,d=1.0 run_single.sh

module load scuff-em

GROUP=${g}
DISP=${d}

GEOFILE="scuffgeo/${GROUP}_d${DISP}.scuffgeo"
EPFILE="evpoints/${GROUP}_d${DISP}.txt"
OUTFILE="results/${GROUP}_d${DISP}.out"

mkdir -p logs results

echo "Running: group=${GROUP}, displacement=${DISP}"
echo "Geometry: ${GEOFILE}"
echo "EPFile: ${EPFILE}"
echo "Output: ${OUTFILE}"

scuff-neq \
    --geometry ${GEOFILE} \
    --EPFile ${EPFILE} \
    --OmegaFile omega.txt \
    --power \
    --OMIT_SELFTERMS \
    --OutFile ${OUTFILE}

echo "Done: ${GROUP} d=${DISP}"
