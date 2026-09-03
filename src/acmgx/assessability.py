"""Which ACMG/AMP criteria can an automated system actually assess?

This split is the basis of the ceiling analysis (Arm B'). A criterion is
machine-assessable if it can be decided from the variant plus public sequence
and population resources. It is not machine-assessable if deciding it requires
wet-lab work, family data, or a clinical cohort the machine does not have.

The distinction is deliberately generous to automation: PM1 and PP2/BP1 are
counted as assessable even though their usual implementations lean on curated
variant databases, because a defensible version can be built from constraint
and domain annotation alone.
"""

# Decidable from sequence, gene models, population frequency and calibrated predictors.
MACHINE_ASSESSABLE = {
    "PVS1",  # null variant + NMD + LoF mechanism
    "PM2", "BA1", "BS1",  # gnomAD grpmax filtering AF
    "PP3", "BP4",  # calibrated in-silico predictor
    "PM4", "BP3",  # protein length change / in-frame indel in repeat
    "BP7",  # synonymous with no predicted splice effect
    "PM1",  # hotspot / functional domain
    "PP2", "BP1",  # gene-level missense constraint
}

# Decidable, but only by reading a curated assertion database. These are the
# criteria the blinding experiment removes.
KNOWLEDGE_DERIVED = {
    "PS1", "PM5",  # same / different amino acid change previously called pathogenic
    "PP5", "BP6",  # "reputable source" - deprecated by ClinGen SVI
}

# Require evidence no automated pipeline has: functional assays, pedigrees,
# case-control cohorts, confirmed parentage, or clinical detail.
NOT_MACHINE_ASSESSABLE = {
    "PS2", "PM6",   # de novo, with and without confirmed parentage
    "PS3", "BS3",   # functional assay
    "PS4",          # prevalence in affecteds vs controls
    "PM3", "BP2",   # in trans / in cis observations in patients
    "PP1", "BS4",   # co-segregation
    "PP4",          # phenotype specificity for the gene
    "BS2",          # observed in healthy individuals
    "BP5",          # alternate molecular cause identified
}

ALL_KNOWN = MACHINE_ASSESSABLE | KNOWLEDGE_DERIVED | NOT_MACHINE_ASSESSABLE


def bucket(base_code: str) -> str:
    """One of 'assessable', 'knowledge', 'unassessable', 'unknown'."""
    if base_code in MACHINE_ASSESSABLE:
        return "assessable"
    if base_code in KNOWLEDGE_DERIVED:
        return "knowledge"
    if base_code in NOT_MACHINE_ASSESSABLE:
        return "unassessable"
    return "unknown"
