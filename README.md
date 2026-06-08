# ACMG Rule Implementation

This repository automates small variant classification according to ACMG guidelines using the GRCh38 reference genome. It integrates VEP annotations and multiple evidence sources to determine pathogenicity.


## Environment Setup
### 1. Create the environment
Use the `environment.yml` file to install validated dependencies, including `polars`.

```bash
ml biology Anaconda/Anaconda3
conda env create -f environment.yml
```

### 2. Test the environment
Verify that the environment was created correctly and all dependencies are functional.

```bash
# Verify the presence of the environment
conda info --envs | grep acmg_rule

# Test key dependencies without manual activation
conda run -n acmg_rule python -c "import polars; print('Environment Test: SUCCESS (polars ' + polars.__version__ + ')')"
```


## Usage
### Input
Update the `SAMPLE_DIR` variable in `acmg_batch_submitter.sh` to point to the directory containing VEP-annotated TSV files.
> [!NOTE] 
> The pipeline will concurrently process all VEP-annotated TSV files found in this directory.

### Running pipeline
Execute the following command to start the analysis:

```bash
bash acmg_batch_submitter.sh
```

The script will automatically create an `ACMG_output` subdirectory within the current working directory to store all results.


## Configuration & Parameters
You can customize parameter cutoffs in `config.json`.

### Default mode
Set a parameter to `null` to use predefined values.

**Predefined Default Values:**
| Parameter | Default value | Note |
| --- | ---: | --- |
| AD_AF | 0.00001441287 | Allele frequency threshold for AD, XL, YL, or unspecified inheritance modes |
| AR_AF | 0.0001 | Allele frequency threshold for AR | 
| PVS1_LOEUF | 0.755 | follow [varsome](https://varsome.com/about/resources/germline-implementation/#pvs1) |
| PP3_CADD_SUP | 25.3 | follow [ClinGen guideline](https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/pdf/main.pdf) |
| PP3_CADD_MOD | 28.1 | follow [ClinGen guideline](https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/pdf/main.pdf) |
| BP4_CADD_SUP | 22.7 | follow [ClinGen guideline](https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/pdf/main.pdf) |
| BP4_CADD_MOD | 17.3 | follow [ClinGen guideline](https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/pdf/main.pdf) |
| BP4_CADD_STR | 0.15 | follow [ClinGen guideline](https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/pdf/main.pdf) |
| PP3_DANN | 0.999 | follow [varsome](https://varsome.com/about/resources/germline-implementation/#insilicopredictions) |
| BP4_DANN_SUP | 0.974 | follow [varsome](https://varsome.com/about/resources/germline-implementation/#insilicopredictions) |
| BP4_DANN_MOD | 0.915 | follow [varsome](https://varsome.com/about/resources/germline-implementation/#insilicopredictions) |
| BP4_DANN_STR | 0.478 | follow [varsome](https://varsome.com/about/resources/germline-implementation/#insilicopredictions) |
| P_SPLICEAI | 0.2 | SpliceAI prediction score threshold for pathogenic |
| B_SPLICEAI | 0.1 | SpliceAI prediction score threshold for benign |

### Custom Mode
Enter a specific numeric value to override defaults.

For example, to override allele frequency settings:

```json
"acmg_params": {
    "AD_AF": 0.00001,
    "AR_AF": 0.0005
}
```

> [!Important] 
> Ensure the file maintains valid JSON syntax when switching between `null` and specific values.


## Output description
All results will be saved in `ACMG_output`.

### Result files
- `{sample}.vep.ACMG.tsv`: VEP annotated TSV with pathogenic classfication and ACMG rules.
- `{sample}.vep.mane_plus_clinical.ACMG.tsv`: Only variants located in 65 MANE Plus Clinical transcripts were included, pathogenic classfication and ACMG rules were added to the annotated TSV.

### Key Output Columns
| Column | Description |
| :--- | :--- |
| Pathogenicity_class | Final classification (e.g., Pathogenic, Likely_pathogenic, Uncertain_significance, Likely_benign, Benign). |
| ACMG_rules | A list of criteria met by the variant (e.g., PVS1,PM2,PP3). |


## Databases for pathogenicity evaluation
Detailed versions and descriptions for all databases can be found in the VEP annotation result documentations.
Key sources used for ACMG scoring include:
- ClinVar, version 20251109
- CADD, version 1.7, using dbNSFP version 4.9a
- DANN, using dbNSFP version 4.9a
- HGMD
- LOEUF, based on gnomAD v2.1.1, liftover from GRCh37
- DVD, version 9.2
- gnomAD genomes coverage, version 3.0.1
- gnomAD genomes, version 4.1
- gnomAD exomes, version 4.1
- RepeatMasker, download from UCSC Table Browser
- SpliceAI SNV, version 1.3
- SpliceAI indel, version 1.3
- MitoMap, version 20230621

<!--
## Rules
### Implemented rule notes
#### PS3
We only apply the PS3 criterion to GoF variants, as these are more likely to be overlooked by in silico prediction tools.

### Unimplemented rules
- PS2
- PS4
- PM3
- PM6
- PP1
- PP2
- PP4
- BS3
- BS4
- BP1
- BP2
- BP5
--->