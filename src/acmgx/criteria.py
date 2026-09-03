"""ACMG/AMP criterion implementations, current ClinGen SVI guidance.

Every criterion is individually switchable. That is the whole point: the
experiment turns criteria off and measures what happens, so the engine cannot
bake any of them in.

Deliberate departures from the 2015 paper, all following ClinGen SVI:

* PM2 is applied at Supporting strength, never Moderate.
* PVS1 emits a strength from the ClinGen decision tree, not a flat VeryStrong.
* PP3/BP4 use one calibrated predictor with graded strengths. Stacking several
  correlated predictors, as ANNOVAR-style majority voting does, is not allowed.
* PP5/BP6 exist only so the blinding experiment has something to remove; they
  are deprecated for real use.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional

# ---------------------------------------------------------------- thresholds
# REVEL, from Pejaver et al. 2022 (AJHG), Table 2. These are score intervals
# calibrated to evidence strength; they are not tool defaults.
REVEL_PP3_STRONG      = 0.932    # >= 0.932
REVEL_PP3_MODERATE    = 0.773    # [0.773, 0.932)
REVEL_PP3_SUPPORTING  = 0.644    # [0.644, 0.773)
REVEL_BP4_SUPPORTING  = 0.290    # (0.183, 0.290]
REVEL_BP4_MODERATE    = 0.183    # (0.016, 0.183]
REVEL_BP4_STRONG      = 0.016    # (0.003, 0.016]
REVEL_BP4_VERYSTRONG  = 0.003    # <= 0.003

# gnomAD v4 joint grpmax filtering AF (95% CI lower bound).
BA1_FAF = 0.05          # ACMG 2015 stand-alone benign
BS1_FAF = 0.001         # generic; VCEPs set this per gene-disease
PM2_FAF = 0.0001        # "absent or extremely rare"

LOF_CONSEQUENCES = {"nonsense", "frameshift_variant",
                    "splice_donor_variant", "splice_acceptor_variant"}
SPLICE_CONSEQUENCES = {"splice_donor_variant", "splice_acceptor_variant",
                       "splice_site_variant"}

# ClinGen dosage: 3 = sufficient evidence for haploinsufficiency,
# 30 = gene associated with an autosomal recessive phenotype (LoF is the
# mechanism), 40 = dosage sensitivity unlikely.
HI_LOF_ESTABLISHED = {"3", "30"}
HI_LOF_REFUTED = {"40"}


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: Optional[str] = None
    consequence: str = "unknown"
    revel: Optional[float] = None
    faf95: Optional[float] = None
    clinvar_tier: Optional[str] = None       # aggregate ClinVar classification
    clinvar_revstat: str = ""
    aa_pos: Optional[int] = None
    aa_ref: Optional[str] = None
    aa_alt: Optional[str] = None
    hi_score: Optional[str] = None
    mane: Optional[dict] = None
    ps1_hit: bool = False                    # same AA change, known P/LP
    pm5_hit: bool = False                    # different AA change at residue, known P/LP


# ---------------------------------------------------------------- PVS1
def nmd_escapes(v: Variant) -> Optional[bool]:
    """True if a PTC here escapes nonsense-mediated decay.

    Escape if the variant lies in the last exon, or within the last 50 nt of the
    penultimate exon. Returns None when the transcript structure is unavailable.
    """
    t = v.mane
    if not t or t["n_exons"] < 1:
        return None
    exons = t["exons"]                # already in transcription order
    if len(exons) == 1:
        return True
    last, penult = exons[-1], exons[-2]
    if last[0] <= v.pos <= last[1]:
        return True
    if penult[0] <= v.pos <= penult[1]:
        dist = (penult[1] - v.pos) if t["strand"] == "+" else (v.pos - penult[0])
        return dist < 50
    return False


def pvs1(v: Variant) -> Optional[str]:
    if v.consequence not in LOF_CONSEQUENCES:
        return None
    if v.hi_score in HI_LOF_REFUTED:
        return None                          # LoF is not the disease mechanism

    established = v.hi_score in HI_LOF_ESTABLISHED
    if v.consequence in ("nonsense", "frameshift_variant"):
        escapes = nmd_escapes(v)
        if escapes is True:
            base = "PVS1_Strong"             # NMD-escaping PTC, truncates protein
        elif escapes is False:
            base = "PVS1"
        else:
            base = "PVS1_Strong"             # unknown structure -> be conservative
    else:                                    # canonical splice site
        base = "PVS1_Strong"                 # skipping vs retention not modelled

    if not established:
        # gene mechanism not curated by ClinGen: drop one strength level
        base = {"PVS1": "PVS1_Strong", "PVS1_Strong": "PVS1_Moderate"}[base]
    return base


# ---------------------------------------------------------------- frequency
def ba1(v: Variant) -> Optional[str]:
    return "BA1" if v.faf95 is not None and v.faf95 >= BA1_FAF else None

def bs1(v: Variant) -> Optional[str]:
    if v.faf95 is None:
        return None
    return "BS1" if BS1_FAF <= v.faf95 < BA1_FAF else None

def pm2(v: Variant) -> Optional[str]:
    if v.faf95 is None:
        return None
    return "PM2_Supporting" if v.faf95 < PM2_FAF else None


# ---------------------------------------------------------------- in silico
def pp3(v: Variant) -> Optional[str]:
    if v.revel is None or v.consequence != "missense_variant":
        return None
    r = v.revel
    if r >= REVEL_PP3_STRONG:     return "PP3_Strong"
    if r >= REVEL_PP3_MODERATE:   return "PP3_Moderate"
    if r >= REVEL_PP3_SUPPORTING: return "PP3"
    return None

def bp4(v: Variant) -> Optional[str]:
    if v.revel is None or v.consequence != "missense_variant":
        return None
    r = v.revel
    if r <= REVEL_BP4_VERYSTRONG: return "BP4_Very Strong"
    if r <= REVEL_BP4_STRONG:     return "BP4_Strong"
    if r <= REVEL_BP4_MODERATE:   return "BP4_Moderate"
    if r <= REVEL_BP4_SUPPORTING: return "BP4"
    return None


# ---------------------------------------------------------------- protein
def pm4(v: Variant) -> Optional[str]:
    """In-frame indel or stop-loss changing protein length."""
    if v.consequence in ("inframe_deletion", "inframe_insertion", "stop_lost"):
        return "PM4"
    return None

def bp7(v: Variant) -> Optional[str]:
    """Synonymous with no predicted splice impact.

    Splice impact is approximated by distance to the nearest exon boundary,
    because SpliceAI was not run (see docs/limitations.md). The approximation is
    conservative: a synonymous variant within 3 nt of a boundary does not get BP7.
    """
    if v.consequence != "synonymous_variant":
        return None
    t = v.mane
    if t:
        for s, e in t["exons"]:
            if s <= v.pos <= e:
                if min(v.pos - s, e - v.pos) < 3:
                    return None
    return "BP7"


# ------------------------------------------------------ ClinVar-derived
def ps1(v: Variant) -> Optional[str]:
    return "PS1" if v.ps1_hit else None

def pm5(v: Variant) -> Optional[str]:
    return "PM5" if v.pm5_hit else None

def pp5(v: Variant) -> Optional[str]:
    if v.clinvar_tier in ("P", "LP") and "no_assertion" not in (v.clinvar_revstat or ""):
        return "PP5"
    return None

def bp6(v: Variant) -> Optional[str]:
    if v.clinvar_tier in ("B", "LB") and "no_assertion" not in (v.clinvar_revstat or ""):
        return "BP6"
    return None


# ---------------------------------------------------------------- registry
CRITERIA = {
    "PVS1": pvs1, "BA1": ba1, "BS1": bs1, "PM2": pm2,
    "PP3": pp3, "BP4": bp4, "PM4": pm4, "BP7": bp7,
    "PS1": ps1, "PM5": pm5, "PP5": pp5, "BP6": bp6,
}

# Criteria this engine does not attempt. Declared, never silently skipped.
NOT_IMPLEMENTED = {
    "PM1": "requires curated hotspot/domain annotation not assembled here",
    "PP2": "requires gene-level missense constraint (gnomAD constraint not retrieved)",
    "BP1": "requires gene-level missense constraint",
    "BP3": "requires repeat-region annotation",
    "PS2": "de novo status - not available for a database benchmark set",
    "PM6": "de novo without confirmed parentage - not available",
    "PS3": "functional assay", "BS3": "functional assay",
    "PS4": "case-control prevalence",
    "PM3": "in trans observations in patients", "BP2": "in cis / in trans observations",
    "PP1": "co-segregation", "BS4": "lack of segregation",
    "PP4": "phenotype specificity", "BS2": "observed in healthy individuals",
    "BP5": "alternate molecular cause",
}

CLINVAR_DERIVED = {"PS1", "PM5", "PP5", "BP6"}
DEPRECATED = {"PP5", "BP6"}


def evaluate(v: Variant, disabled: set[str] | None = None) -> list[str]:
    """Fire every enabled criterion and return the codes that apply."""
    disabled = disabled or set()
    out = []
    for name, fn in CRITERIA.items():
        if name in disabled:
            continue
        code = fn(v)
        if code:
            out.append(code)
    return out
