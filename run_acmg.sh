#!/bin/bash
#SBATCH -p ngs26G
#SBATCH -c 4
#SBATCH --mem=26g
#SBATCH -A MST109178
#SBATCH -J ACMG
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL,END

WKDIR=$PWD
LOGDIR=$WKDIR/logs
OUTDIR=$WKDIR/ACMG_output

mkdir -p $LOGDIR
mkdir -p $OUTDIR

# log file print info setting
source $WKDIR/utils/job_utils.sh
set -euo pipefail


# Arguments
FILE=("$@")
SAMPLE_INPUT=${FILE[$SLURM_ARRAY_TASK_ID]}
SAMPLE_ID=$(basename "$SAMPLE_INPUT" .tsv)
SAMPLE_OUTPUT="${SAMPLE_ID}.ACMG.tsv"
configfile="$WKDIR/utils/config.json"

TIME=`date +%Y%m%d%H%M`
logfile=${LOGDIR}/${TIME}_ACMG_${SAMPLE_ID}.log


# call function from job_utils.sh to initialize log file
start_job

module load Anaconda/Anaconda3
conda activate acmg_rule


# Run python
python ./main.py \
    --config $configfile \
    --input $SAMPLE_INPUT \
    --output $OUTDIR/$SAMPLE_OUTPUT \
    --sampleID $SAMPLE_ID 


rm -f $WKDIR/slurm-*.out
