#!/bin/bash
#SBATCH -p ngs26G
#SBATCH -c 4
#SBATCH --mem=26g
#SBATCH -A MST109178
#SBATCH -J ACMG
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL,END

# log file print info setting
source /staging/biology/r12455009/test_vep/small_variant/VEP_script/utils/job_utils.sh
set -euo pipefail

WKDIR=$PWD
LOGDIR=$WKDIR/logs
OUTDIR=$WKDIR/ACMG_output

mkdir -p $LOGDIR
mkdir -p $OUTDIR

# Arguments
FILE=("$@")
SAMPLE_INPUT=${FILE[$SLURM_ARRAY_TASK_ID]}
SAMPLE_ID=$(basename "$SAMPLE_INPUT" .tsv)
SAMPLE_OUTPUT="test_${SAMPLE_ID}.tsv"
configfile="$WKDIR/utils/config.json"

TIME=`date +%Y%m%d%H%M`
logfile=${LOGDIR}/${TIME}_ACMG_${SAMPLE_ID}.log


# call function from job_utils.sh to initialize log file
start_job

module load biology
module load BCFtools/1.18
module load Anaconda/Anaconda3
conda activate acmg_rule


# Run python
python ./main.py \
    --config $configfile \
    --input $SAMPLE_INPUT \
    --output $OUTDIR/$SAMPLE_OUTPUT \
    --sampleID $SAMPLE_ID 


# Output result
echo "$(date '+%Y-%m-%d %H:%M:%S') - Converting ACMG output TSV to VCF"

SAMPLE_VCF=${SAMPLE_INPUT%.tsv}.vcf.gz
ACMG_TSV=$OUTDIR/$SAMPLE_OUTPUT
CSQ_TSV=$OUTDIR/${SAMPLE_ID}_CSQ.tsv.gz

# Prepare ACMG TSV for bcftools
sed -i '1s/CHROM/#CHROM/' $ACMG_TSV
bgzip -f $ACMG_TSV
tabix -s1 -b2 -e2 ${ACMG_TSV}.gz

# Extract CSQ as a real file (process substitution not supported by bcftools)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\n' $SAMPLE_VCF | \
    bgzip -c -f > $CSQ_TSV
tabix -s1 -b2 -e2 $CSQ_TSV

# strip CSQ → add ACMG → add CSQ back
bcftools annotate -x INFO/CSQ $SAMPLE_VCF | \
    bcftools annotate \
        -a ${ACMG_TSV}.gz \
        -h $WKDIR/utils/header.txt \
        -c CHROM,POS,REF,ALT,.LoF_gene,.Inheritance,.MOI_description,.MOI_source,.ClinVar_GOF,.HGMD_GOF,.Pathogenicity_class,.ACMG_rules | \
    bcftools annotate \
        -a $CSQ_TSV \
        -h <(zgrep -m1 '##INFO=<ID=CSQ' $SAMPLE_VCF) \
        -c CHROM,POS,REF,ALT,INFO/CSQ \
        -Oz -o $OUTDIR/${SAMPLE_ID}.ACMG.vcf.gz

bcftools index -t -f $OUTDIR/${SAMPLE_ID}.ACMG.vcf.gz

# Cleanup
rm -f ${ACMG_TSV}.gz ${ACMG_TSV}.gz.tbi $CSQ_TSV ${CSQ_TSV}.tbi
rm -f $WKDIR/slurm-*.out

#### Notes ####
# convert tsv to vcf with vep annotation
# 另存一個沒有 CSQ 的 vcf
# 將 ACMG 的結果加上去
# 最後再將 CSQ 加回去
# 有沒有 chr 要確認
# tsv 要壓縮