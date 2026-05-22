#!/bin/bash

## Change the directory to the location your annotated TSV files are saved
SAMPLE_DIR="/staging/biology/r12455009/test_vep/small_variant/VEP_output/clinvar_20251109.cleaned"


# FILES=($(find "$SAMPLE_DIR" -name '*.vep.tsv'))
# FILES=($(find "$SAMPLE_DIR" -name '*.vep.mane_plus_clinical.tsv'))
FILES=($(find "$SAMPLE_DIR" -name '*.vep*.tsv'))
NUM_SAMPLES=${#FILES[@]}

if [ $NUM_SAMPLES -eq 0 ]; then
    echo "No file found!"
    exit 1
fi

echo "Found $NUM_SAMPLES files. Submitting Job Array..."
sbatch --array=0-$((NUM_SAMPLES-1)) ./run_acmg.sh "${FILES[@]}"
