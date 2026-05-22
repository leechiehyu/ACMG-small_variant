# Classify the pathogenicity of the variants based on `acmg_rules.py`
import inspect
import polars as pl
from modules import acmg_rules

# ================================= #
# =========== Constants =========== #
# ================================= #
# Configuration mapping for user-defined parameters in config.json to function parameters in acmg_rules.py
### "user-defined parameter name in config.json": "function parameter name in acmg_rules.py"
CONFIG_MAP = {
    "AD_AF":        "ad_cutoff",
    "AR_AF":        "ar_cutoff",
    "PVS1_LOEUF":   "pvs1_loeuf_cutoff",
    "PP3_CADD_SUP": "p_cadd_sup_cutoff",
    "PP3_CADD_MOD": "p_cadd_m_cutoff",
    "BP4_CADD_SUP": "b_cadd_sup_cutoff",
    "BP4_CADD_MOD": "b_cadd_m_cutoff",
    "BP4_CADD_STR": "b_cadd_s_cutoff",
    "PP3_DANN":     "p_dann_cutoff",
    "BP4_DANN_SUP": "b_dann_sup_cutoff",
    "BP4_DANN_MOD": "b_dann_m_cutoff",
    "BP4_DANN_STR": "b_dann_s_cutoff",
    "P_SPLICEAI":   "p_spliceai_cutoff",
    "B_SPLICEAI":   "b_spliceai_cutoff",
}

# Disable dictionary
### Disable map: if key is True → disable all listed targets
DISABLE_MAP = {
    "RAW_BA1":         ["RAW_BS1", "RAW_BS2"], 
    "RAW_PVS1_STRONG": ["RAW_PVS1", "RAW_PP3", "RAW_PP3_MODERATE", 
                        "RAW_BP4", "RAW_BP4_MODERATE", "RAW_BP4_STRONG", "RAW_PM4"],
    "RAW_PVS1":        ["RAW_PP3", "RAW_PP3_MODERATE", 
                        "RAW_BP4", "RAW_BP4_MODERATE", "RAW_BP4_STRONG", "RAW_PM4"], 
    "RAW_PM1":         ["RAW_BP3"],
    "RAW_PM2":         ["RAW_BS1", "RAW_BS2"],
    "RAW_PM4":         ["RAW_PP3", "RAW_PP3_MODERATE"],
}

# Define point weights for each rule
## Weight: Very Strong(8), Strong(4), Moderate(2), Supporting(1)
P_WEIGHT_FAMILIES = {
    "PVS1_grp": {"PVS1": 8, "PVS1_STRONG": 4},  # PVS1 had been disabled by PVS1_STRONG in the step 2
    "PS1_grp": {"PS1": 4},
    "PS3_grp": {"PS3": 4},
    "PM1_grp": {"PM1": 2},
    "PM2_grp": {"PM2": 2},
    "PM4_grp": {"PM4": 2},
    "PM5_grp": {"PM5": 2},
    "PP3_grp": {"PP3": 1, "PP3_MODERATE": 2},
    "PP5_grp": {"PP5": 1, "PP5_MODERATE": 2, "PP5_VERYSTRONG": 8},
}

B_WEIGHT_FAMILIES = {
    "BA1_grp": {"BA1": -8},
    "BS1_grp": {"BS1": -4},
    "BS2_grp": {"BS2": -4},
    "BP3_grp": {"BP3": -1},
    "BP4_grp": {"BP4": -1, "BP4_MODERATE": -2, "BP4_STRONG": -4},
    "BP6_grp": {"BP6": -1, "BP6_STRONG": -4, "BP6_VERYSTRONG": -8},
    "BP7_grp": {"BP7": -1},
}



# ================================= #
# ============ Helper ============= #
# ================================= #
def _build_group_exprs(families: dict[str, dict[str, int]], descending: bool):
    """
    For step 3 point-based scoring
    For each group return:
      - one score expression (max or min across members)
      - one label expression (name of the active rule with highest/lowest weight)
    """
    score_exprs, label_exprs = [], []
 
    agg_fn   = pl.max_horizontal if not descending else pl.min_horizontal
    sort_rev = descending  # ascending for P (pick highest last), descending for B
 
    for labels in families.values():
        # Score: pick the group's extreme score
        score_exprs.append(
            agg_fn([pl.col(name).cast(pl.Int8) * w for name, w in labels.items()]).fill_null(0)
        )
        # Label: iterate from least to most extreme so the last `.when` wins
        sorted_items = sorted(labels.items(), key=lambda x: x[1], reverse=sort_rev)
        expr = pl.lit(None)
        for name, _ in sorted_items:
            expr = pl.when(pl.col(name)).then(pl.lit(name)).otherwise(expr)
        label_exprs.append(expr)
 
    return score_exprs, label_exprs



# ================================= #
# ======== Pipeline Steps ========= #
# ================================= #
def step1_calculate_raw_labels(lf: pl.LazyFrame, **user_config) -> pl.LazyFrame:
    """
    Call every check_* function in acmg_rules, passing only the params each function accepts
    user_config: dictionary of user-defined parameters (e.g., ad, ar)
    """
    expressions = []
    
    for name in dir(acmg_rules):
        if not name.startswith("check_"):
            continue
 
        func = getattr(acmg_rules, name)
        func_params = inspect.signature(func).parameters.keys()
 
        matched_args = {
            func_param: user_config[user_key]
            for user_key, func_param in CONFIG_MAP.items()
            if user_key in user_config
            and user_config[user_key] is not None
            and func_param in func_params
        }
 
        col_name = name.replace("check_", "raw_").upper()
        expressions.append(func(**matched_args).alias(col_name))
 
    return lf.with_columns(expressions)


def step2_apply_disables(lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Conflict Resolution (disable rules based on other rules)
    Define the final ACMG classification of each variant
    """
    # find all involved labels and initialize them (raw_ -> final name)
    raw_cols = [c for c in lf.collect_schema().names() if c.upper().startswith("RAW_")]

    # Pre-index: target → list of inhibitors
    inhibitor_index: dict[str, list[str]] = {}
    for trigger, targets in DISABLE_MAP.items():
        for t in targets:
            inhibitor_index.setdefault(t, []).append(trigger)
 
    expressions = []
    for raw_col in raw_cols:
        final_name = raw_col.replace("RAW_", "")
        inhibitors = inhibitor_index.get(raw_col, [])
 
        if inhibitors:
            disabled = pl.any_horizontal(pl.col(i) for i in inhibitors)
            expr = pl.when(disabled).then(pl.lit(False)).otherwise(pl.col(raw_col))
        else:
            expr = pl.col(raw_col)
 
        expressions.append(expr.alias(final_name))
 
    return lf.with_columns(expressions)


##### POINT-BASED #####
def step3_scoring(lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Point-based scoring and final pathogenicity classification
    """
    p_score_exprs, p_label_exprs = _build_group_exprs(P_WEIGHT_FAMILIES, descending=False)
    b_score_exprs, b_label_exprs = _build_group_exprs(B_WEIGHT_FAMILIES, descending=True)
 
    all_label_exprs = p_label_exprs + b_label_exprs

    # Scoring and classification
    return (
        lf.with_columns(
            pl.sum_horizontal(p_score_exprs).alias("total_p_score"),
            pl.sum_horizontal(b_score_exprs).alias("total_b_score"),
            pl.concat_list(all_label_exprs)
              .list.drop_nulls()
              .list.join(",")
              .map_elements(lambda x: "." if x == "" else x, return_dtype=pl.String)
              .alias("ACMG_rules"),
        )
        .with_columns(
            (pl.col("total_p_score") + pl.col("total_b_score")).alias("final_score")
        )
        .with_columns(
            pl.when(pl.col("final_score") <= -4)              .then(pl.lit("Benign"))
            .when(pl.col("final_score").is_between(-3, -1))   .then(pl.lit("Likely_benign"))
            .when(pl.col("final_score").is_between( 0,  5))   .then(pl.lit("Uncertain_significance"))
            .when(pl.col("final_score").is_between( 6,  9))   .then(pl.lit("Likely_pathogenic"))
            .when(pl.col("final_score") >= 10)                .then(pl.lit("Pathogenic"))
            .otherwise(pl.lit("Uncertain_significance"))
            .alias("Pathogenicity_class")
        )
    )



# ================================= #
# ========= Entry point =========== #
# ================================= #
def run_acmg_pipeline(lf: pl.LazyFrame, **user_config) -> pl.LazyFrame:
    """
    ACMG classification pipeline
    Aggregates all steps: raw label calculation, conflict resolution, scoring and final classification
    """
    input_cols = lf.collect_schema().names()

    return (
        lf.pipe(step1_calculate_raw_labels, **user_config)
          .pipe(step2_apply_disables)
          .pipe(step3_scoring)
          .select(input_cols + ["Pathogenicity_class", "ACMG_rules"])
    )

