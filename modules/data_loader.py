# Deal with all data loading and merging
import polars as pl

# ================================= #
# =========== Constants =========== #
# ================================= #
_VARIANT_SCHEMA: dict[str, type] = {
    "POS":                          pl.Int32,
    "DISTANCE":                     pl.Int32,
    "gnomADe_AF":                   pl.Float32,
    "gnomADe_AFR_AF":               pl.Float32,
    "gnomADe_AMR_AF":               pl.Float32,
    "gnomADe_ASJ_AF":               pl.Float32,
    "gnomADe_EAS_AF":               pl.Float32,
    "gnomADe_FIN_AF":               pl.Float32,
    "gnomADe_MID_AF":               pl.Float32,
    "gnomADe_NFE_AF":               pl.Float32,
    "gnomADe_REMAINING_AF":         pl.Float32,
    "gnomADe_SAS_AF":               pl.Float32,
    "gnomADg_AF":                   pl.Float32,
    "gnomADg_AFR_AF":               pl.Float32,
    "gnomADg_AMI_AF":               pl.Float32,
    "gnomADg_AMR_AF":               pl.Float32,
    "gnomADg_ASJ_AF":               pl.Float32,
    "gnomADg_EAS_AF":               pl.Float32,
    "gnomADg_FIN_AF":               pl.Float32,
    "gnomADg_MID_AF":               pl.Float32,
    "gnomADg_NFE_AF":               pl.Float32,
    "gnomADg_REMAINING_AF":         pl.Float32,
    "gnomADg_SAS_AF":               pl.Float32,
    "ada_score":                    pl.Float64,
    "rf_score":                     pl.Float64,
    "CADD_phred":                   pl.Float64,
    "CADD_raw":                     pl.Float64,
    "DANN_rankscore":               pl.Float64,
    "DANN_score":                   pl.Float64,
    "GERP___RS":                    pl.Float64,
    "LRT_score":                    pl.Float64,
    "M_CAP_score":                  pl.Float64,
    "MetaLR_score":                 pl.Float64,
    "MetaSVM_score":                pl.Float64,
    "fathmm_MKL_coding_score":      pl.Float64,
    "SpliceAI_pred_DP_AG":          pl.Int32,
    "SpliceAI_pred_DP_AL":          pl.Int32,
    "SpliceAI_pred_DP_DG":          pl.Int32,
    "SpliceAI_pred_DP_DL":          pl.Int32,
    "SpliceAI_pred_DS_AG":          pl.Float64,
    "SpliceAI_pred_DS_AL":          pl.Float64,
    "SpliceAI_pred_DS_DG":          pl.Float64,
    "SpliceAI_pred_DS_DL":          pl.Float64,
    "LOEUF":                        pl.Float64,
    "PrimateAI":                    pl.Float64,
    "MES_NCSS_downstream_acceptor": pl.Float64,
    "MES_NCSS_downstream_donor":    pl.Float64,
    "MES_NCSS_upstream_acceptor":   pl.Float64,
    "MES_NCSS_upstream_donor":      pl.Float64,
    "MES_SWA_acceptor_alt":         pl.Float64,
    "MES_SWA_acceptor_diff":        pl.Float64,
    "MES_SWA_acceptor_ref":         pl.Float64,
    "MES_SWA_acceptor_ref_comp":    pl.Float64,
    "MES_SWA_donor_alt":            pl.Float64,
    "MES_SWA_donor_diff":           pl.Float64,
    "MES_SWA_donor_ref":            pl.Float64,
    "MES_SWA_donor_ref_comp":       pl.Float64,
    "MaxEntScan_alt":               pl.Float64,
    "MaxEntScan_diff":              pl.Float64,
    "MaxEntScan_ref":               pl.Float64,
    "pHaplo":                       pl.Float64,
    "pTriplo":                      pl.Float64,
    "LoFtool":                      pl.Float64,
    "pLI_gene_value":               pl.Float64,
    "TWB1490_SNV_AF":               pl.Float64,
    "TWB_official_SNV_AF":          pl.Float64,
    "gnomAD_exome_AN":              pl.Int64,
    "gnomAD_exome_AF":              pl.Float64,
    "gnomAD_exome_nhomalt":         pl.Int64,
    "gnomAD_exome_AN_eas":          pl.Int64,
    "gnomAD_exome_AF_eas":          pl.Float64,
    "gnomAD_exome_nhomalt_eas":     pl.Int64,
    "gnomAD_genome_AN":             pl.Int64,
    "gnomAD_genome_AF":             pl.Float64,
    "gnomAD_genome_nhomalt":        pl.Int64,
    "gnomAD_genome_AN_eas":         pl.Int64,
    "gnomAD_genome_AF_eas":         pl.Float64,
    "gnomAD_genome_nhomalt_eas":    pl.Int64,
    "TWB_mtDNA_AF_het_vaf05":       pl.Float64,
    "gnomAD_mtDNA_AF_hom_eas":      pl.Float64,
    "gnomAD_mtDNA_AF_het_eas":      pl.Float64,
}

_CHROM_RANK: dict[str, int] = {
    **{str(i): i for i in range(1, 23)},
    "X": 23, "Y": 24, "MT": 25,
}

_HOTSPOT_CONSEQUENCES = ["inframe_insertion", "inframe_deletion", "missense_variant"]



# ================================= #
# ======== Shared Helpers ========= #
# ================================= #
def _standardize_chrom(col: str = "CHROM") -> pl.Expr:
    """
    Standardize chromosome names: 
    1. remove 'chr' prefix 
    2. normalize 'M' to 'MT'
    """
    return (
        pl.col(col)
        .str.replace(r"^(?i)chr", "")
        .str.replace(r"^M$", "MT")
    )

def _scan_tsv(path: str, schema: dict, **kwargs) -> pl.LazyFrame:
    """
    Wrapper around pl.scan_csv for TSV files with standardized settings
    """
    return pl.scan_csv(
        path,
        separator = '\t',
        null_values = ["."],
        schema_overrides = schema,
        infer_schema_length = 0,
        **kwargs,
    )

def _add_chrom_rank(lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Add a sorting key column based on chromosome order
     (For attach_pathogenic_hotspot function)
     """
    return lf.with_columns(
        pl.col("CHROM").replace(_CHROM_RANK, default=99).cast(pl.Int8).alias("chrom_rank")
    )

def _join_by_id_then_symbol(
    input_lf: pl.LazyFrame,
    db: pl.LazyFrame,
    cols: list[str],
    rename: dict[str, str],
) -> pl.LazyFrame:
    """
    Join `db` to `input_lf` twice — first on HGNC_ID, then on SYMBOL —
    and coalesce the results into final column names defined by `rename`.
 
    `cols` must be the non-key columns to carry over (excluding HGNC_ID / SYMBOL).
    """
    lf = input_lf.join(
        db.select(["HGNC_ID"] + cols).rename({c: f"{c}_by_id" for c in cols}),
        on="HGNC_ID", how="left",
    ).join(
        db.select(["SYMBOL"] + cols).rename({c: f"{c}_by_sym" for c in cols}),
        on="SYMBOL", how="left",
    )
 
    tmp_cols = [f"{c}_by_id" for c in cols] + [f"{c}_by_sym" for c in cols]
 
    return (
        lf.with_columns([
            pl.coalesce([f"{src}_by_id", f"{src}_by_sym"]).alias(dst)
            for src, dst in rename.items()
        ])
        .drop(tmp_cols)
    )



# ================================= #
# ============ Loaders ============ #
# ================================= #
def load_input_variants(file_path: str) -> pl.LazyFrame:
    """
    Read sample variants from TSV; 
    all unlisted columns default to String
    """
    # Catch all column names in the header
    headers = pl.read_csv(
        file_path, 
        separator = '\t', 
        n_rows = 0,
        infer_schema_length = 0
    ).columns
    # Set all columns to String type initially
    variant_schema = {col: pl.String for col in headers}
    # Define the data type for input file columns
    variant_schema.update(_VARIANT_SCHEMA)    

    return (
        pl.scan_csv(
            file_path,
            separator = '\t',
            null_values = ["."],
            schema_overrides = variant_schema,
            infer_schema_length = 0,
        )
        .with_columns(
            pl.col("CHROM").alias("_raw_CHROM"),    # Save original CHROM
            pl.col("HGVSp").str.replace(r"^.*:", "").alias("_tmp_aa_change"),
        )
        .with_columns(
            _standardize_chrom("CHROM"),
        )
        .with_row_index("original_order")
    )


def attach_lof_gene_info(input_lf: pl.LazyFrame, lof_genelist_path: str) -> pl.LazyFrame:
    """
    Join LoF gene list; 
    adds boolean column `LoF_gene`
    Rule(s): PVS1
    """
    lof_schema = {
        "SYMBOL": pl.String,
        "HGNC_ID": pl.String
    }

    lof_db = (
        _scan_tsv(
            lof_genelist_path, 
            schema = {"SYMBOL": pl.String, "HGNC_ID": pl.String},
            truncate_ragged_lines=True,
        )
        .select([
            pl.col("SYMBOL"),
            pl.col("HGNC_ID"),
            pl.lit(True).alias("is_lof")
        ])
    )

    return (
        _join_by_id_then_symbol(
            input_lf, lof_db,
            cols=["is_lof"],
            rename={"is_lof": "LoF_gene"},
        )
        .with_columns(pl.col("LoF_gene").fill_null(False))
    )


def attach_moi_info(input_lf: pl.LazyFrame, moi_db_path: str) -> pl.LazyFrame:
    """
    Join MOI database; 
    adds `Inheritance`, `MOI_description`, `MOI_source`
    Rule(s): PM2, BS2
    """
    moi_schema = {
        "SYMBOL": pl.String,
        "HGNC_ID": pl.String,
        "inheritance": pl.String,
        "description": pl.String,
        "source": pl.String
    }

    moi_db = (
        _scan_tsv(
            moi_db_path, 
            schema=moi_schema, 
            truncate_ragged_lines=True
        )
        .with_columns(
            pl.col("description").str.replace_all(";", "|")
        )
        .unique(subset=["HGNC_ID", "SYMBOL"])
    )

    # Keep the original columns of moi_db
    carry_cols = ["inheritance", "description", "source"]
    rename_cols = {
        "inheritance": "Inheritance",
        "description": "MOI_description",
        "source": "MOI_source"
    }
    
    return (
        _join_by_id_then_symbol(
            input_lf, moi_db, 
            cols=carry_cols, 
            rename=rename_cols)
    )


def attach_gof_info(input_lf: pl.LazyFrame, gof_db_path: str) -> pl.LazyFrame:
    """
    Join GoF database on CHROM/POS/REF/ALT; 
    adds `ClinVar_GOF`, `HGMD_GOF`
    Rule(s): PS3
    """
    gof_schema = {
        "CHROM": pl.String,
        "POS": pl.Int32,
        "REF": pl.String,
        "ALT": pl.String,
        "ClinVar_GOF": pl.String,
        "HGMD_GOF": pl.String
    }

    gof_db = (
        _scan_tsv(gof_db_path, schema = gof_schema)
        .with_columns(_standardize_chrom("CHROM"))
        .group_by(["CHROM", "POS", "REF", "ALT"])
        .agg(
            pl.col("ClinVar_GOF").drop_nulls().unique().str.concat("|"),
            pl.col("HGMD_GOF").drop_nulls().unique().str.concat("|"),
        )
    )

    return input_lf.join(gof_db, on=["CHROM", "POS", "REF", "ALT"], how="left")


def attach_pathogenic_hotspot(input_lf: pl.LazyFrame, hotspot_path: str) -> pl.LazyFrame:
    """
    Count pathogenic ClinVar variants within a ±25 bp window for each
    inframe-indel / missense position; result stored in `patho_count`
    Rule(s): PM1
    """
    patho_schema = {
        "CHROM": pl.String,
        "POS": pl.Int32
    }
    csq_pattern = "|".join(_HOTSPOT_CONSEQUENCES)

    # ClinVar missense & inframe indel pathogenic DB
    clinvar_ref = (
        _scan_tsv(hotspot_path, schema = patho_schema)
        .with_columns(_standardize_chrom("CHROM"))

        ## consider the "position"
        # .unique(subset=["CHROM", "POS"])

        ## consider the "variant number"
        .group_by(["CHROM", "POS"])
        .agg(pl.len().cast(pl.Int32).alias("is_patho"))

        ## flag these positions as "clinvar"
        .with_columns(pl.lit("clinvar").alias("point_type"))
        .pipe(_add_chrom_rank)
    )

    # extract only inframe indel and missense variants from the input_lf
    input_coords = (
        input_lf
        .filter(pl.col("Consequence").str.contains(csq_pattern))
        .select(
            pl.col("CHROM"), pl.col("POS"),
            pl.lit(0).alias("is_patho"),
            pl.lit("sample").alias("point_type"),
        )
        .unique()
        .pipe(_add_chrom_rank)
    )
 

    # combine input and pathogenicDB and count pathogenic hotspot
    density_df = (
        pl.concat([input_coords, clinvar_ref], how = "diagonal")
        .group_by(["CHROM", "POS", "chrom_rank"])
        .agg(
            pl.col("is_patho").sum().fill_null(0),
            pl.when((pl.col("point_type") == "sample").any())
              .then(pl.lit("sample"))
              .otherwise(pl.lit("clinvar"))
              .alias("point_type"),
        )
        .sort(["chrom_rank", "POS"])
        .rolling(
            index_column = "POS", 
            period = "51i", 
            offset = "-25i", 
            by = "CHROM"
        )
        .agg(
            pl.col("is_patho").sum().fill_null(0).alias("patho_count"),
            pl.col("point_type").first(),
        )
        .filter(pl.col("point_type") == "sample")
        .select(["CHROM", "POS", "patho_count"])
    )

    # join result to the original input df
    return (
        input_lf.join(density_df, on=["CHROM", "POS"], how="left")
        .with_columns(
            pl.when(
                (pl.col("Consequence").str.contains(csq_pattern)) & 
                (pl.col("CHROM") != "MT")
            )
            .then(pl.col("patho_count"))
            .otherwise(None)
            .alias("patho_count")
        )
    )


def attach_pathogenicDB_variant_exception(input_lf: pl.LazyFrame, variant_path: str) -> pl.LazyFrame:
    """
    Join pathogenicDB on CHROM/POS/REF/ALT; 
    adds `Variant`, `Exception`
    Rule(s): PP5, BA1, BP6
    """
    patho_schema = {
        "CHROM": pl.String,
        "POS": pl.Int32,
        "REF": pl.String,
        "ALT": pl.String,
        "Variant": pl.String,
        "Exception": pl.String
    }

    patho_db = (
        _scan_tsv(
            variant_path, 
            schema = patho_schema, 
            truncate_ragged_lines = True
        )
        .with_columns(_standardize_chrom("CHROM"))
        .unique(subset=["CHROM", "POS", "REF", "ALT"])
    )

    return input_lf.join(patho_db, on=["CHROM", "POS", "REF", "ALT"], how="left")



# ================================= #
# ========= Entry point =========== #
# ================================= #
def prepare_full_data(
    input_path: str, 
    lof_genelist_path: str,
    moi_db_path: str, 
    gof_db_path: str, 
    hotspot_path: str, 
    variant_path: str
) -> pl.LazyFrame:
    """
    Load input variants and attach all annotation databases
    """
    return (
        # load case input and standardize the chromosome names
        load_input_variants(input_path)
        # attach LoF info
        .pipe(attach_lof_gene_info, lof_genelist_path = lof_genelist_path)
        # attach MOI info
        .pipe(attach_moi_info, moi_db_path = moi_db_path)
        # attach GoF info
        .pipe(attach_gof_info, gof_db_path = gof_db_path)
        # attach pathogenic hotspot info
        .pipe(attach_pathogenic_hotspot, hotspot_path = hotspot_path)
        # attach pathogenicDB variant exception info
        .pipe(attach_pathogenicDB_variant_exception, variant_path = variant_path)
    )
