#!/usr/bin/env python3
"""Tests for the ACMG engine.

These check the things that would silently corrupt a result rather than crash:
combining-rule boundaries, strength-modifier parsing, the NMD rule PVS1 depends
on, and the calibrated predictor intervals.

Run: python tests/test_acmgx.py
"""
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from acmgx.combine import combine_2015, combine_points, points, parse_code
from acmgx.criteria import Variant, evaluate, pvs1, pp3, bp4, pm2, ba1, bs1, nmd_escapes
from acmgx.assessability import bucket

fails = []
def check(name, got, want):
    if got != want:
        fails.append(f"{name}: got {got!r}, want {want!r}")

# ---------------------------------------------------------------- parsing
check("parse plain", parse_code("PM2"), ("PM2", "Moderate", "P"))
check("parse downgrade", parse_code("PM2_Supporting"), ("PM2", "Supporting", "P"))
check("parse upgrade", parse_code("PM3_Strong"), ("PM3", "Strong", "P"))
check("parse spaced", parse_code("PM3_Very Strong"), ("PM3", "VeryStrong", "P"))
check("parse benign", parse_code("BP4_Moderate"), ("BP4", "Moderate", "B"))
check("parse standalone", parse_code("BA1"), ("BA1", "StandAlone", "B"))

# -------------------------------------------------------- 2015 combining
# Richards 2015 Table 5, pathogenic rules
check("P: PVS1+PS", combine_2015(["PVS1", "PS3"]), "P")
check("P: PVS1+2PM", combine_2015(["PVS1", "PM2", "PM4"]), "P")
check("P: PVS1+PM+PP", combine_2015(["PVS1", "PM2", "PP3"]), "P")
check("P: PVS1+2PP", combine_2015(["PVS1", "PP3", "PP4"]), "P")
check("P: 2PS", combine_2015(["PS1", "PS3"]), "P")
check("P: PS+3PM", combine_2015(["PS3", "PM1", "PM2", "PM4"]), "P")
# likely pathogenic
check("LP: PVS1+PM", combine_2015(["PVS1", "PM2"]), "LP")
check("LP: PS+PM", combine_2015(["PS3", "PM2"]), "LP")
check("LP: 3PM", combine_2015(["PM1", "PM2", "PM4"]), "LP")
check("LP: PM+4PP", combine_2015(["PM2", "PP1", "PP2", "PP3", "PP4"]), "LP")
# benign
check("B: BA1", combine_2015(["BA1"]), "B")
check("B: 2BS", combine_2015(["BS1", "BS2"]), "B")
check("LB: BS+BP", combine_2015(["BS1", "BP4"]), "LB")
check("LB: 2BP", combine_2015(["BP4", "BP7"]), "LB")
# uncertain, and conflicting evidence resolves to VUS not to a side
check("VUS: single PM", combine_2015(["PM2"]), "VUS")
check("VUS: nothing", combine_2015([]), "VUS")
check("VUS: conflicting", combine_2015(["PVS1", "PS3", "BA1"]), "VUS")
# a lone PVS1 is not pathogenic on its own under the 2015 rules
check("VUS: lone PVS1", combine_2015(["PVS1"]), "VUS")

# ------------------------------------------------------------- points
check("points PVS1", points(["PVS1"]), 8)
check("points strong", points(["PS3"]), 4)
check("points modifier", points(["PM3_Strong"]), 4)
check("points benign sign", points(["BA1"]), -8)
check("points mixed", points(["PVS1", "PM2_Supporting", "PP3"]), 10)
check("tier at 10", combine_points(["PVS1", "PM2_Supporting", "PP3"]), ("P", 10))
check("tier at 6", combine_points(["PS3", "PM2"]), ("LP", 6))
check("tier at 5", combine_points(["PS3", "PP3"]), ("VUS", 5))
check("tier at -7", combine_points(["BS1", "BS2"]), ("B", -8))

# ------------------------------------------------------ REVEL intervals
# Pejaver et al. 2022 boundaries; the interval edges are where mistakes hide.
def mis(revel):
    return Variant("chr1", 1, "A", "G", gene="X", consequence="missense_variant", revel=revel)
check("PP3 strong at 0.932", pp3(mis(0.932)), "PP3_Strong")
check("PP3 moderate below", pp3(mis(0.9319)), "PP3_Moderate")
check("PP3 moderate at 0.773", pp3(mis(0.773)), "PP3_Moderate")
check("PP3 supporting at 0.644", pp3(mis(0.644)), "PP3")
check("PP3 none below 0.644", pp3(mis(0.6439)), None)
check("BP4 none at 0.291", bp4(mis(0.291)), None)
check("BP4 supporting at 0.290", bp4(mis(0.290)), "BP4")
check("BP4 moderate at 0.183", bp4(mis(0.183)), "BP4_Moderate")
check("BP4 strong at 0.016", bp4(mis(0.016)), "BP4_Strong")
check("BP4 very strong at 0.003", bp4(mis(0.003)), "BP4_Very Strong")
# a nonsense variant is not scored by a missense predictor
check("PP3 not on nonsense",
      pp3(Variant("chr1", 1, "A", "T", consequence="nonsense", revel=0.99)), None)

# -------------------------------------------------------------- frequency
def freq(f):
    return Variant("chr1", 1, "A", "G", consequence="missense_variant", faf95=f)
check("BA1 at 5%", ba1(freq(0.05)), "BA1")
check("BA1 not below", ba1(freq(0.049)), None)
check("BS1 in band", bs1(freq(0.002)), "BS1")
check("BS1 not above BA1", bs1(freq(0.06)), None)
check("PM2 when rare", pm2(freq(0.00001)), "PM2_Supporting")
check("PM2 not when common", pm2(freq(0.01)), None)
check("PM2 needs data", pm2(Variant("chr1", 1, "A", "G")), None)

# ------------------------------------------------------------------ NMD
# three-exon gene on the + strand: exons 1-100, 201-300, 401-500
plus = {"n_exons": 3, "strand": "+", "exons": [(1, 100), (201, 300), (401, 500)]}
def at(pos, tx=plus):
    return Variant("chr1", pos, "C", "T", consequence="nonsense", mane=tx)
check("NMD: first exon", nmd_escapes(at(50)), False)
check("NMD: last exon escapes", nmd_escapes(at(450)), True)
check("NMD: penultimate, far from end", nmd_escapes(at(210)), False)
check("NMD: penultimate, last 50nt", nmd_escapes(at(295)), True)
# minus strand: the "end" of the penultimate exon is its lower coordinate
minus = {"n_exons": 3, "strand": "-", "exons": [(401, 500), (201, 300), (1, 100)]}
check("NMD: minus last exon", nmd_escapes(at(50, minus)), True)
check("NMD: minus penultimate near end", nmd_escapes(at(205, minus)), True)
check("NMD: minus penultimate far", nmd_escapes(at(290, minus)), False)
check("NMD: single exon", nmd_escapes(at(50, {"n_exons": 1, "strand": "+", "exons": [(1, 500)]})), True)
check("NMD: unknown structure", nmd_escapes(Variant("chr1", 1, "C", "T")), None)

# ------------------------------------------------------------------ PVS1
def lof(cons, hi, pos=50, tx=plus):
    return Variant("chr1", pos, "C", "T", consequence=cons, hi_score=hi, mane=tx)
check("PVS1 nonsense, NMD, HI known", pvs1(lof("nonsense", "3")), "PVS1")
check("PVS1 nonsense, NMD, HI unknown", pvs1(lof("nonsense", None)), "PVS1_Strong")
check("PVS1 refuted by dosage", pvs1(lof("nonsense", "40")), None)
check("PVS1 escaping NMD downgrades", pvs1(lof("nonsense", "3", pos=450)), "PVS1_Strong")
check("PVS1 splice", pvs1(lof("splice_donor_variant", "3")), "PVS1_Strong")
check("PVS1 not on missense", pvs1(lof("missense_variant", "3")), None)
check("PVS1 AR gene counts as established", pvs1(lof("frameshift_variant", "30")), "PVS1")

# --------------------------------------------------------------- toggles
v = Variant("chr1", 1, "A", "G", gene="X", consequence="missense_variant",
            revel=0.95, faf95=0.0, clinvar_tier="P", clinvar_revstat="criteria_provided",
            ps1_hit=True, pm5_hit=True)
full = set(evaluate(v))
svi = set(evaluate(v, {"PP5", "BP6"}))
blind = set(evaluate(v, {"PP5", "BP6", "PS1", "PM5"}))
check("full fires PP5", "PP5" in full, True)
check("svi drops PP5", "PP5" in svi, False)
check("svi keeps PS1", "PS1" in svi, True)
check("blind drops PS1/PM5", {"PS1", "PM5"} & blind, set())
check("blinding only removes", blind <= svi <= full, True)

# ------------------------------------------------------- assessability
check("bucket PVS1", bucket("PVS1"), "assessable")
check("bucket PP5", bucket("PP5"), "knowledge")
check("bucket PS3", bucket("PS3"), "unassessable")
check("bucket unknown", bucket("PZ9"), "unknown")

if fails:
    print(f"FAILED {len(fails)}:")
    for f in fails:
        print("  " + f)
    sys.exit(1)
print("all acmgx tests pass")
