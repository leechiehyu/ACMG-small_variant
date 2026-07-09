# Rescore ACMG classification
#
# Usage:
#   python rescore.py --input {sample}.ACMG.tsv --output {sample}.ACMG.rescore.tsv
# or
#   python rescore.py -i {sample}.ACMG.tsv -o {sample}.ACMG.rescore.tsv
#
# ============================================================
# User-defined scores — modify here
# ============================================================

# Point weight per evidence strength
STRENGTH_SCORES = {
    "very_strong":  8,
    "strong":       4,
    "moderate":     2,
    "supporting":   1,
}

# Classification cutoffs (inclusive)
# Pathogenic  : score >= P_CUTOFF
# Likely path.: P_CUTOFF > score >= LP_CUTOFF
# Likely ben. : LB_CUTOFF >= score > B_CUTOFF
# Benign      : score <= B_CUTOFF
CUTOFFS = {
    "Pathogenic":        10,
    "Likely_pathogenic":  6,
    "Likely_benign":     -1,
    "Benign":            -4,
}

# Rule → evidence strength
# Pathogenic rules (positive) and Benign rules (negative)
RULE_STRENGTH = {
    # Pathogenic
    "PVS1":            ("pathogenic", "very_strong"),
    "PVS1_STRONG":     ("pathogenic", "strong"),
    "PS1":             ("pathogenic", "strong"),
    "PS3":             ("pathogenic", "strong"),
    "PM1":             ("pathogenic", "moderate"),
    "PM2":             ("pathogenic", "moderate"),
    "PM4":             ("pathogenic", "moderate"),
    "PM5":             ("pathogenic", "moderate"),
    "PP3":             ("pathogenic", "supporting"),
    "PP3_MODERATE":    ("pathogenic", "moderate"),
    "PP5":             ("pathogenic", "supporting"),
    "PP5_MODERATE":    ("pathogenic", "moderate"),
    "PP5_STRONG":      ("pathogenic", "strong"),
    "PP5_VERYSTRONG":  ("pathogenic", "very_strong"),
    # Benign
    "BA1":             ("benign", "very_strong"),
    "BS1":             ("benign", "strong"),
    "BS2":             ("benign", "strong"),
    "BP3":             ("benign", "supporting"),
    "BP4":             ("benign", "supporting"),
    "BP4_MODERATE":    ("benign", "moderate"),
    "BP4_STRONG":      ("benign", "strong"),
    "BP6":             ("benign", "supporting"),
    "BP6_STRONG":      ("benign", "strong"),
    "BP6_VERYSTRONG":  ("benign", "very_strong"),
    "BP7":             ("benign", "supporting"),
}

# ============================================================
# Script logic — no need to modify below
# ============================================================

import argparse
import polars as pl


def build_weights(rule_strength: dict, strength_scores: dict) -> dict[str, int]:
    weights = {}
    for rule, (direction, strength) in rule_strength.items():
        score = strength_scores[strength]
        weights[rule] = score if direction == "pathogenic" else -score
    return weights


def rescore(lf: pl.LazyFrame, weights: dict[str, int], cutoffs: dict) -> pl.LazyFrame:
    score_expr = pl.lit(0)
    for rule, score in weights.items():
        score_expr = score_expr + (
            pl.col("ACMG_rules")
            .str.contains(r"(?:^|,)" + rule + r"(?:,|$)")
            .cast(pl.Int32) * score
        )

    p_cutoff  = cutoffs["Pathogenic"]
    lp_cutoff = cutoffs["Likely_pathogenic"]
    lb_cutoff = cutoffs["Likely_benign"]
    b_cutoff  = cutoffs["Benign"]

    return (
        lf
        .with_columns(score_expr.alias("custom_score"))
        .with_columns(
            pl.when(pl.col("custom_score") >= p_cutoff) .then(pl.lit("Pathogenic"))
            .when(pl.col("custom_score") >= lp_cutoff)  .then(pl.lit("Likely pathogenic"))
            .when(pl.col("custom_score") <= b_cutoff)   .then(pl.lit("Benign"))
            .when(pl.col("custom_score") <= lb_cutoff)  .then(pl.lit("Likely benign"))
            .otherwise(pl.lit("Uncertain significance"))
            .alias("Custom_class")
        )
        # .drop("custom_score")
        .drop("Custom_class")
    )


def main():
    parser = argparse.ArgumentParser(description="Custom ACMG Rescoring")
    parser.add_argument("--input",  "-i", type=str, required=True, help="Path to {sample}.ACMG.tsv")
    parser.add_argument("--output", "-o", type=str, required=True, help="Path to output TSV")
    args = parser.parse_args()

    weights = build_weights(RULE_STRENGTH, STRENGTH_SCORES)
    lf = pl.scan_csv(
        args.input, 
        separator = "\t", 
        null_values = ".",
        infer_schema_length = 0,
        )
    rescore(lf, weights, CUTOFFS).sink_csv(args.output, separator="\t", null_value=".")
    print(f"Done: {args.output}")


if __name__ == "__main__":
    main()