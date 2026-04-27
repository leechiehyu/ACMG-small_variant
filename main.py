import sys
import time
import json
import logging
import argparse
import polars as pl
from modules.data_loader import prepare_full_data
from modules.aggregation import run_acmg_pipeline


OUTPUT_COLUMNS = [
    # --- Variant coordinates ---
    "CHROM", "POS", "REF", "ALT",
    # --- Gene / transcript ---
    "LoF_gene", "Inheritance", "MOI_description", "MOI_source", "ClinVar_GOF", "HGMD_GOF",
    # --- ACMG result ---
    "Pathogenicity_class", "ACMG_rules",
]


def parse_args():
    """
    Recieve the arguments from the bash file
    """
    parser = argparse.ArgumentParser(description="ACMG Rule Pipeline")
    parser.add_argument("--config",   type=str, default="config.json", help="Path to config.json")
    parser.add_argument("--input",    type=str, required=True,         help="Path to input TSV")
    parser.add_argument("--output",   type=str, required=True,         help="Path to output TSV")
    parser.add_argument("--sampleID", type=str, required=True,         help="Sample identifier")
    return parser.parse_args()


def setup_logging():
    """
    Setup logging configuration
    """
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)s] %(asctime)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.StreamHandler()],
    )


def main():
    args = parse_args()
    setup_logging()

    start_time = time.time()
    logging.info(f"Starting ACMG Pipeline for {args.sampleID}")

    # 1. Load config
    try:
        with open(args.config, "r", encoding="utf-8") as f:
            config = json.load(f)

    except Exception as e:
        logging.error(f"Failed to load config: {str(e)}")
        sys.exit(1)

    # 2. Resolve hotspot path based on sampleID suffix
    db_paths = config.get("database_path", {})
    is_clinical = args.sampleID.endswith(".mane_plus_clinical")
    hotspot_path = db_paths.get("P_HOTSPOT_CLINICAL") if is_clinical else db_paths.get("P_HOTSPOT_SELECT")

    # 3. Load and annotate input data
    ## using data_loader.py
    try:
        logging.info(f"Loading data from: {args.input}")

        lf = prepare_full_data(
            input_path=args.input,
            lof_genelist_path=db_paths.get("LOF_GENELIST"),
            moi_db_path=db_paths.get("MOI"),
            gof_db_path=db_paths.get("GOF"),
            hotspot_path=hotspot_path,
            variant_path=db_paths.get("VAR_EXP")
        )

    except Exception as e:
        logging.error(f"Data loading failed: {str(e)}")
        sys.exit(1)

    # 4. Run ACMG classification
    ## using aggregation.py
    try:
        logging.info("Constructing ACMG pipeline...")

        # tmp_cols = ["_tmp_aa_change", "chrom_rank", "original_order", "_raw_CHROM"]
        # schema_names = set(lf.collect_schema().names())

        final_lf = (
            run_acmg_pipeline(lf, **config.get("acmg_params", {}))
            .sort("original_order")
            .with_columns(pl.col("_raw_CHROM").alias("CHROM"))  # restore original CHROM
            # .drop([c for c in tmp_cols if c in schema_names])
            .select(OUTPUT_COLUMNS)
        )
        
        logging.info("Pipeline constructed successfully.")

    except Exception as e:
        logging.error(f"ACMG calculation failed: {str(e)}")
        sys.exit(1)

    # 5. Write output
    try:
        logging.info(f"Saving pipeline temporary results to: {args.output}")

        final_lf.sink_csv(args.output, separator="\t", null_value=".")

    except Exception as e:
        logging.error(f"Output write failed: {e}")
        sys.exit(1)

    # 6. Summary
    minutes, seconds = divmod(time.time() - start_time, 60)

    print("-" * 40)
    print("PIPELINE SUMMARY")
    print(f"SUCCESS : {args.sampleID.split('.vep')[0]}")
    print(f"MANE    : {'MANE plus clinical' if is_clinical else 'MANE select'}")
    print(f"Time    : {int(minutes)} min {seconds:.2f} sec")
    print("-" * 40, "\n")


if __name__ == "__main__":
    main()
