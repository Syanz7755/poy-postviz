#!/bin/bash
# Submit all scuff-neq jobs
# Usage: bash submit_all.sh

mkdir -p logs results

# Create omega.txt if not exists (angular frequencies in rad/s)
if [ ! -f omega.txt ]; then
    echo "Creating omega.txt from material file..."
    grep -v "^#" materials/silica-27p.txt | awk '{print $1}' > omega.txt
    echo "Created omega.txt with $(wc -l < omega.txt) frequencies"
fi

# Group a, b, c: displacements 1.0, 2.0, 4.0, 5.0, 10.0
for group in a b c; do
    for d in 1.0 2.0 4.0 5.0 10.0; do
        echo "Submitting: group=${group}, d=${d}"
        sbatch --export=g=${group},d=${d} run_single.sh
    done
done

# Group d, e, f: displacements 0.5, 1.0, 2.0, 5.0
for group in d e f; do
    for d in 0.5 1.0 2.0 5.0; do
        echo "Submitting: group=${group}, d=${d}"
        sbatch --export=g=${group},d=${d} run_single.sh
    done
done

echo "All 27 jobs submitted"
