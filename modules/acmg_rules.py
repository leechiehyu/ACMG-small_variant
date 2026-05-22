# Define all criteria of ACMG rules
import operator
import polars as pl

# ================================= #
# =========== Constants =========== #
# ================================= #
### For PVS1
LOF_CONSEQUENCES = [
    "stop_gained", "frameshift_variant",
    "splice_donor_variant", "splice_acceptor_variant",
]

### For PM4, BP3
INFRAME_CONSEQUENCES = ["inframe_deletion", "inframe_insertion"]

### For PP3, BP4, BP7
SPLICE_AI_COLS = [
    "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL",
]

### For _moi_ad()
INHERITANCE_AD = ["AD", "XL", "YL"]



# ================================= #
# ======== Shared Helpers ========= #
# ================================= #
# -------------------------- #
# gnomAD coverage preprocess #
# -------------------------- #
### for PM2, BA1, BS1, BS2
def _gnomAD_min_cov():
    return (
        pl.col("gnomAD_genome_cov")
        .str.split("&")
        .list.eval(pl.element().cast(pl.Float64, strict=False))
        .list.min()
    )

# -------------- #
# gnomAD variant #
# -------------- #
### for PM2, BA1, BS1, BS2
def _gnomad_qc(mode: str):
    """
    Exist in gnomAD pass QC: 
    FILTER is PASS, AN > 2000, coverage > 20
    """
    return (
        (pl.col(f"gnomAD_{mode}_FILTER") == "PASS") & 
        (pl.col(f"gnomAD_{mode}_AN") > 2000) & 
        (_gnomAD_min_cov() > 20)
    ).fill_null(False)

def _gnomad_absent(mode: str):
    """
    Absent in gnomAD pass QC: 
    AN is null, coverage > 20
    """
    return (
        (pl.col(f"gnomAD_{mode}_AN").is_null()) & 
        (_gnomAD_min_cov() > 20)
    )

def _gnomad_af(mode: str, cutoff: float, op=operator.le):
    """
    check gnomAD af cutoff
    """
    af = pl.col(f"gnomAD_{mode}_AF")
    return op(af, cutoff).fill_null(False)

### for PM2, BS2
def _gnomad_nhomalt(mode: str, op=operator.le):
    """
    check gnomAD nhomalt cutoff
    """
    nhomalt = pl.col(f"gnomAD_{mode}_nhomalt")
    return op(nhomalt, 2).fill_null(False)

    #### lt(a, b) -> a < b
    #### le(a, b) -> a <= b
    #### gt(a, b) -> a > b
    #### ge(a, b) -> a >= b

# ------------------- #
# Mode of Inheritance #
# ------------------- #
### for PM2, BS2
def _moi_ad():
    """ inheritance AD """
    return pl.col("Inheritance").is_in(INHERITANCE_AD).fill_null(True)

def _moi_ar():
    """ inheritance AR """
    return (pl.col("Inheritance") == "AR").fill_null(False)

# -------------------- #
# In-silico prediction #
# -------------------- #
### for PP3, BP4
def _in_silico_score(col: str, cutoff: float, op=operator.ge):
    """ Generic in-silico predictor cutoff check """
    return op(pl.col(col), cutoff)

def _cadd_or_dann(cadd_cutoff: float, dann_cutoff: float, op=operator.ge):
    ### keeping the first non-null value between CADD and DANN
    return pl.coalesce([
        _in_silico_score("CADD_phred",     cadd_cutoff, op),
        _in_silico_score("DANN_rankscore", dann_cutoff, op),
    ])

### For BP4, BP7
def _spliceAI_all(cutoff: float):
    """ SpliceAI every DS <= cutoff """
    return pl.all_horizontal(pl.col(c).fill_null(1) <= cutoff for c in SPLICE_AI_COLS)



# ================================= #
# ========== ACMG Rules =========== #
# ================================= #
############################
########### PVS1 ###########
############################
# --- PVS1 sub rules --- #
def _pvs1_lof_csq():
    """ frame change in vep ensembl Effect """
    return pl.col("Consequence").fill_null("").str.contains("|".join(LOF_CONSEQUENCES))

def _pvs1_loeuf(cutoff: float): 
    """ LOEUF cutoff """
    return pl.col("LOEUF").fill_null(1) <= cutoff

def _pvs1_lof_gene():
    """ gene locate in LoF gene """
    return pl.col("LoF_gene")

def _pvs1_nmd_escaping():
    """ NMD escape """
    return (pl.col("NMD") == "NMD_escaping_variant").fill_null(False)

# --- PVS1 main rule --- #
def check_pvs1(pvs1_loeuf_cutoff: float = 0.755):
    return (
        _pvs1_lof_csq() & 
        # _pvs1_loeuf(pvs1_loeuf_cutoff) & 
        # _pvs1_lof_gene()
        (_pvs1_loeuf(pvs1_loeuf_cutoff) | _pvs1_lof_gene())
    ).fill_null(False)

def check_pvs1_strong(pvs1_loeuf_cutoff: float = 0.755):
    """
    PVS1 downgraded to Strong when variant is NMD_escaping_variant.
    """
    return (
        check_pvs1(pvs1_loeuf_cutoff) & 
        _pvs1_nmd_escaping()
    ).fill_null(False)


#############################
############ PS1 ############
#############################
# --- PS1 sub rules --- #
def _ps1_vep_transcript():
    """ same transcript """
    return (
        (pl.col("Feature") == pl.col("CLN_VEP_Feature")) & 
        pl.col("Feature").is_not_null() & 
        pl.col("CLN_VEP_Feature").is_not_null()
    )

def _ps1_vep_aachange():
    """ same amino acid change """
    return (
        (pl.col("_tmp_aa_change") == pl.col("CLN_VEP_AAchange")) & 
        pl.col("_tmp_aa_change").is_not_null() &
        pl.col("CLN_VEP_AAchange").is_not_null()
    )

# --- PS1 main rule --- #
def check_ps1():
    """
    Same transcript + same amino acid change as a known pathogenic variant
    """
    return (
        _ps1_vep_transcript() & 
        _ps1_vep_aachange()
    ).fill_null(False)


#############################
############ PS3 ############
#############################
# --- PS3 sub rules --- #
def _ps3_gofclinvar():
    """ GoF knowed by predicted ClinVar """
    return pl.col("ClinVar_GOF").is_not_null()

def _ps3_gofhgmd():
    """ GOF knowed by predicted HGMD """
    return pl.col("HGMD_GOF").is_not_null()

# --- PS3 main rule --- #
def check_ps3():
    """
    GoF evidence from ClinVar or HGMD
    """
    return (
        _ps3_gofclinvar() | 
        _ps3_gofhgmd()
    ).fill_null(False)


#############################
############ PM1 ############
#############################
# --- PM1 rule --- #
def check_pm1():
    """
    >= 4 pathogenic variants within 25 bp
    """
    return (pl.col("patho_count") >= 4).fill_null(False)


#############################
############ PM2 ############
#############################
# --- PM2 main rule --- #
def check_pm2(ad_cutoff: float = 0.00001441287, ar_cutoff: float = 0.0001):
    """
    Low / absent allele frequency in gnomAD
    """
    def _criteria(use_exome: bool, use_genome: bool):

        e_ad   = _gnomad_af("exome",  ad_cutoff)
        g_ad   = _gnomad_af("genome", ad_cutoff)
        e_ar   = _gnomad_af("exome",  ar_cutoff)
        g_ar   = _gnomad_af("genome", ar_cutoff)
        e_nhom = _gnomad_nhomalt("exome")
        g_nhom = _gnomad_nhomalt("genome")
 
        if use_exome and use_genome:
            ad = _moi_ad() & e_ad & g_ad
            ar = _moi_ar() & ((e_ar & g_ar) | (e_nhom & g_nhom))
        elif use_exome: 
            ad = _moi_ad() & e_ad
            ar = _moi_ar() & (e_ar | e_nhom)
        else:  # genome only
            ad = _moi_ad() & g_ad
            ar = _moi_ar() & (g_ar | g_nhom)
        return ad | ar
 
    qc_e, qc_g = _gnomad_qc("exome"), _gnomad_qc("genome")
    ab_e, ab_g = _gnomad_absent("exome"), _gnomad_absent("genome")
 
    return (
        (qc_e & qc_g  & _criteria(True,  True))  |  # both exome and genome have data
        (qc_e & ab_g  & _criteria(True,  False)) |  # only exome have data and genome is absent
        (ab_e & qc_g  & _criteria(False, True))  |  # only genome have data and exome is absent
        (ab_e & ab_g)                               # both exome and genome are absent
    ).fill_null(False)


#############################
############ PM4 ############
#############################
# --- PM4 sub rules --- #
def _pm4_csq_inframe():
    """ inframe change """
    return pl.col("Consequence").fill_null("").str.contains("|".join(INFRAME_CONSEQUENCES))

def _pm4_repeatMasker():
    """ This variant isn't located in a repeat region, which is annotated by rmsk database """
    return pl.col("RepeatMasker").is_null()

def _pm4_csq_stoploss():
    """ stop lost """
    return pl.col("Consequence").fill_null("").str.contains("stop_lost")

# --- PM4 main rule --- #
def check_pm4():
    """
    In-frame indel outside repeat region, or stop-loss
    """
    return (
        (_pm4_csq_inframe() & _pm4_repeatMasker()) |
        _pm4_csq_stoploss()
    ).fill_null(False)


#############################
############ PM5 ############
#############################
# --- PM5 sub rules --- #
def _pm5_transcript():
    """ same transcript """
    return (
        (pl.col("Feature") == pl.col("CLN_VEP_Feature")) & 
        (pl.col("Feature").is_not_null()) & 
        pl.col("CLN_VEP_Feature").is_not_null()
    )

def _pm5_aapos():
    """ same position of aa change """
    return (
        (pl.col("Protein_position") == pl.col("CLN_VEP_Protein_position")) & 
        (pl.col("Protein_position").is_not_null()) & 
        pl.col("CLN_VEP_Protein_position").is_not_null()
    )

# --- PM5 main rule --- #
def check_pm5():
    """
    Different pathogenic missense at the same amino acid position
    """
    return (
        _pm5_transcript() & 
        _pm5_aapos()
    ).fill_null(False)


#############################
############ PP3 ############
#############################
# --- PP3 sub rules --- #
def _pp3_spliceAI_any(cutoff: float):
    """ SpliceAI any DS cutoff >= 0.2 """
    return pl.any_horizontal(pl.col(c).fill_null(0) >= cutoff for c in SPLICE_AI_COLS)

# --- PP3 main rule --- #
def check_pp3(
    p_cadd_sup_cutoff: float = 25.3, 
    p_dann_cutoff: float = 0.999, 
    p_spliceai_cutoff: float = 0.2
):
    return (
        _cadd_or_dann(p_cadd_sup_cutoff, p_dann_cutoff, operator.ge) | 
        _pp3_spliceAI_any(p_spliceai_cutoff)
    ).fill_null(False)

def check_pp3_moderate(
    p_cadd_m_cutoff: float = 28.1, 
    p_dann_cutoff: float = 0.999, 
    p_spliceai_cutoff: float = 0.2
):
    return (
        _cadd_or_dann(p_cadd_m_cutoff, p_dann_cutoff, operator.ge) | 
        _pp3_spliceAI_any(p_spliceai_cutoff)
    ).fill_null(False)


#############################
############ PP5 ############
#############################
# --- PP5 sub rules --- #
def _pp5_clinvar_clnsig():
    """ ClinVar pathogenic/likely_pathogenic """
    return pl.col("ClinVar_CLNSIG").str.contains("Pathogenic|Likely_pathogenic").fill_null(False)

def _pp5_clinvar_clnrevstat():
    """ ClinVar review status """
    PP5_STATUSES = [
        "no_assertion_criteria_provided", 
        "criteria_provided&_single_submitter", 
        "criteria_provided&_multiple_submitters&_no_conflicts", 
        "reviewed_by_expert_panel", 
        "practice_guideline"
    ]
    return pl.col("ClinVar_CLNREVSTAT").is_in(PP5_STATUSES).fill_null(False)

def _pp5_supporting():
    """ By ClinVar PLP 0 or 1 star """
    PLP_01_STAR = [
        "ClinVar_P0", 
        "ClinVar_P1"
    ]
    return pl.col("Variant").is_in(PLP_01_STAR).fill_null(False)

def _pp5_verystrong():
    """ By clinvar PLP 3 or 4 star, DVD, MitoMap """
    PLP_34_STAR = [
        "ClinVar_P3", 
        "ClinVar_P4",
        "DeafnessVD",
        "MitoMap"
    ]
    return pl.col("Variant").is_in(PLP_34_STAR).fill_null(False)

def _pp5_dvd():
    """ DeafnessVD CLNSIG P/LP """
    return pl.col("DVD_SNV_Variant_Classification").str.contains("Pathogenic|Likely_pathogenic").fill_null(False)

# --- PP5 main rule --- #
def check_pp5():
    """
    ClinVar P/LP with 0-1 star review status
    """
    return (
        _pp5_clinvar_clnsig() & 
        _pp5_clinvar_clnrevstat() & 
        _pp5_supporting()
    ).fill_null(False)

def check_pp5_moderate():
    """
    ClinVar P/LP with 2 star review status
    """
    return pl.col("Variant").fill_null("") == "ClinVar_P2"

def check_pp5_verystrong():
    """
    ClinVar P/LP 3-4 star, DeafnessVD, or MitoMap
    """
    return (
        _pp5_verystrong() |
        _pp5_dvd()
    ).fill_null(False)


#############################
############ BA1 ############
#############################
# --- BA1 sub rules --- #
def _ba1_exception():
    """ clingen BA1 exception list """
    return pl.col("Exception").is_not_null()

# --- BA1 main rule --- #
def check_ba1():
    """
    High allele frequency (≥0.05) in gnomAD, excluding ClinGen exceptions
    """
    high_af_e = _gnomad_qc("exome")  & _gnomad_af("exome",  0.05, operator.ge)
    high_af_g = _gnomad_qc("genome") & _gnomad_af("genome", 0.05, operator.ge)

    return (
        (high_af_e | high_af_g) & 
        ~_ba1_exception()
    ).fill_null(False)


###########################
########### BS1 ###########
###########################
# --- BS1 main rule --- #
def check_bs1():
    """
    Allele frequency >= 0.01 in gnomAD
    """
    high_af_e = _gnomad_qc("exome")  & _gnomad_af("exome",  0.01, operator.ge)
    high_af_g = _gnomad_qc("genome") & _gnomad_af("genome", 0.01, operator.ge)

    return (high_af_e | high_af_g).fill_null(False)


###########################
########### BS2 ###########
###########################
# --- BS2 main rule --- #
def check_bs2(ad_cutoff: float = 0.00001441287, ar_cutoff: float = 0.0001):
    """
    Observed in healthy adults at frequency above disease expectation
    """
    def _criteria(mode: str):
        return (
            (_moi_ad() & _gnomad_af(mode, ad_cutoff, operator.gt)) |
            (_moi_ar() & (_gnomad_af(mode, ar_cutoff, operator.gt) | _gnomad_nhomalt(mode, operator.gt)))
        )

    ### either gnomAD exome or genome have data meets the criteria
    cand_e = _gnomad_qc("exome")  & _criteria("exome")
    cand_g = _gnomad_qc("genome") & _criteria("genome")

    return (cand_e | cand_g).fill_null(False)


#############################
############ BP3 ############
#############################
# --- BP3 sub rules --- #
def _bp3_csq():
    """ inframe change """
    return pl.col("Consequence").fill_null("").str.contains("|".join(INFRAME_CONSEQUENCES))

def _bp3_indel():
    """ indel variant """
    return pl.col("VARIANT_CLASS").is_in(['deletion', 'insertion']).fill_null(False)

def _bp3_repeatMasker():
    """ variant located in a repeat region, which is annotated by rmsk database """
    return pl.col("RepeatMasker").is_not_null()

def _bp3_mt():
    """ mitochondrial variants don't apply BP3 """
    return (pl.col("CHROM") == "MT").fill_null(False)

# --- BP3 main rule --- #
def check_bp3():
    """
    In-frame indel in repeat region
    (disabled for mitochondrial variants)
    """
    return (
        ~_bp3_mt() &     # Rule BP3 is disabled for mitochondrial variants, in line with ClinGen Guidelines.
        _bp3_csq() & 
        _bp3_indel() & 
        _bp3_repeatMasker()
    ).fill_null(False)


#############################
############ BP4 ############
#############################
# --- BP4 main rule --- #
def check_bp4(
    b_cadd_sup_cutoff: float = 22.7, 
    b_dann_sup_cutoff: float = 0.974, 
    b_spliceai_cutoff: float = 0.1
):
    return (
        _cadd_or_dann(b_cadd_sup_cutoff, b_dann_sup_cutoff, operator.le) & 
        _spliceAI_all(b_spliceai_cutoff)
    ).fill_null(False)

def check_bp4_moderate(
    b_cadd_m_cutoff: float = 17.3, 
    b_dann_m_cutoff: float = 0.915, 
    b_spliceai_cutoff: float = 0.1
):
    return (
        _cadd_or_dann(b_cadd_m_cutoff, b_dann_m_cutoff, operator.le) & 
        _spliceAI_all(b_spliceai_cutoff)
    ).fill_null(False)

def check_bp4_strong(
    b_cadd_s_cutoff: float = 0.15, 
    b_dann_s_cutoff: float = 0.478, 
    b_spliceai_cutoff: float = 0.1
):
    return (
        _cadd_or_dann(b_cadd_s_cutoff, b_dann_s_cutoff, operator.le) &
        _spliceAI_all(b_spliceai_cutoff)
    ).fill_null(False)


#############################
############ BP6 ############
#############################
# --- BP6 sub rules --- #
def _bp6_clinvar_clinsig():
    """ ClinVar benign/likely_benign """
    return pl.col("ClinVar_CLNSIG").str.contains("Benign|Likely_benign").fill_null(False)

def _bp6_clinvar_clnrevstat():
    """ ClinVar review status """
    BP6_STATUSES = [
        "no_assertion_criteria_provided", 
        "criteria_provided&_single_submitter"
    ]
    return pl.col("ClinVar_CLNREVSTAT").is_in(BP6_STATUSES).fill_null(False)

def _bp6_supporting():
    """ By ClinVar BLB 0 or 1 star """
    BLB_01_STAR = [
        "ClinVar_B0", 
        "ClinVar_B1"
    ]
    return pl.col("Variant").is_in(BLB_01_STAR).fill_null(False)

# --- BP6 main rule --- #
def check_bp6():
    """
    ClinVar B/LB with 0-1 star review status
    """
    return (
        _bp6_clinvar_clinsig() & 
        _bp6_clinvar_clnrevstat() & 
        _bp6_supporting()
    ).fill_null(False)

def check_bp6_strong():
    """
    ClinVar B/LB with 2 star review status
    """
    return pl.col("Variant").fill_null("") == "ClinVar_B2"

def check_bp6_verystrong():
    """
    ClinVar B/LB with 3-4 star review status
    """
    BLB_34_STAR = [
        "ClinVar_B3", 
        "ClinVar_B4"
    ]
    return pl.col("Variant").is_in(BLB_34_STAR).fill_null(False)


#############################
############ BP7 ############
#############################
# --- BP7 sub rules --- #
def _bp7_synonymous():
    """ vep synonymous variant """
    return (pl.col("Consequence").str.contains("synonymous_variant")).fill_null(False)

# --- BP7 main rule --- #
def check_bp7(b_spliceai_cutoff: float = 0.1):
    """
    Synonymous variant with no predicted splice impact
    """
    return (
        _bp7_synonymous() & 
        _spliceAI_all(b_spliceai_cutoff)
    ).fill_null(False)
