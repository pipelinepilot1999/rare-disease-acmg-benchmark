"""ACMG/AMP evidence combination.

Two combining rules are implemented and always reported together:

* ``combine_2015``  - the categorical rules of Richards et al. 2015.
* ``combine_points`` - the Bayesian point system of Tavtigian et al. 2018/2020,
  which ClinGen VCEP specifications have been migrating to.

They disagree on a small but interesting set of evidence combinations, and
reporting only one hides that.
"""
from __future__ import annotations

from collections import Counter
from typing import Iterable

# ---------------------------------------------------------------- strengths

VERY_STRONG, STRONG, MODERATE, SUPPORTING = "VeryStrong", "Strong", "Moderate", "Supporting"
STAND_ALONE = "StandAlone"

_PREFIX_STRENGTH = {
    "PVS": VERY_STRONG, "PS": STRONG, "PM": MODERATE, "PP": SUPPORTING,
    "BA": STAND_ALONE, "BS": STRONG, "BP": SUPPORTING,
}

_MODIFIER_ALIASES = {
    "very strong": VERY_STRONG, "verystrong": VERY_STRONG, "very_strong": VERY_STRONG,
    "strong": STRONG, "moderate": MODERATE, "supporting": SUPPORTING,
    "stand alone": STAND_ALONE, "standalone": STAND_ALONE,
}

# Signed points. Benign evidence is the mirror of pathogenic evidence.
_POINTS = {VERY_STRONG: 8, STRONG: 4, MODERATE: 2, SUPPORTING: 1, STAND_ALONE: 8}

TIERS = ("P", "LP", "VUS", "LB", "B")


def parse_code(code: str) -> tuple[str, str, str]:
    """``'PM3_Strong'`` -> ``('PM3', 'Strong', 'P')``.

    Returns (base criterion, applied strength, direction) where direction is
    ``'P'`` for pathogenic-leaning evidence and ``'B'`` for benign-leaning.
    """
    code = code.strip()
    base, _, modifier = code.partition("_")
    prefix = "".join(ch for ch in base if ch.isalpha())
    strength = _PREFIX_STRENGTH.get(prefix)
    if strength is None:                      # unknown criterion
        raise ValueError(f"unrecognised ACMG code: {code!r}")
    if modifier:
        strength = _MODIFIER_ALIASES.get(modifier.strip().lower().replace("_", " "), strength)
    direction = "B" if prefix.startswith("B") else "P"
    return base, strength, direction


# ------------------------------------------------------------- 2015 rules

def combine_2015(codes: Iterable[str]) -> str:
    """Categorical combining rules from Richards et al. 2015, Tables 5 and 6."""
    n: Counter[str] = Counter()
    for code in codes:
        base, strength, direction = parse_code(code)
        if direction == "P":
            n["PVS" if strength == VERY_STRONG else
              "PS" if strength == STRONG else
              "PM" if strength == MODERATE else "PP"] += 1
        else:
            n["BA" if strength == STAND_ALONE else
              "BS" if strength in (VERY_STRONG, STRONG) else
              "BM" if strength == MODERATE else "BP"] += 1

    pvs, ps, pm, pp = n["PVS"], n["PS"], n["PM"], n["PP"]
    ba, bs, bm, bp = n["BA"], n["BS"], n["BM"], n["BP"]
    # 2015 has no benign-moderate; fold it in as the nearest available bucket.
    bp += bm

    pathogenic = (
        (pvs >= 1 and (ps >= 1 or pm >= 2 or (pm == 1 and pp >= 1) or pp >= 2))
        or ps >= 2
        or (ps == 1 and (pm >= 3 or (pm == 2 and pp >= 2) or (pm == 1 and pp >= 4)))
    )
    likely_pathogenic = (
        (pvs == 1 and pm == 1)
        or (ps == 1 and 1 <= pm <= 2)
        or (ps == 1 and pp >= 2)
        or pm >= 3
        or (pm == 2 and pp >= 2)
        or (pm == 1 and pp >= 4)
    )
    benign = ba >= 1 or bs >= 2
    likely_benign = (bs == 1 and bp >= 1) or bp >= 2

    path_side = "P" if pathogenic else "LP" if likely_pathogenic else None
    ben_side = "B" if benign else "LB" if likely_benign else None

    if path_side and ben_side:
        return "VUS"          # criteria in both directions -> conflicting
    return path_side or ben_side or "VUS"


# ------------------------------------------------------------ points rules

def points(codes: Iterable[str]) -> int:
    """Signed Tavtigian points total."""
    total = 0
    for code in codes:
        _, strength, direction = parse_code(code)
        value = _POINTS[strength]
        total += value if direction == "P" else -value
    return total


def combine_points(codes: Iterable[str]) -> tuple[str, int]:
    """Point total mapped to a tier, per Tavtigian et al."""
    p = points(codes)
    if p >= 10:
        tier = "P"
    elif p >= 6:
        tier = "LP"
    elif p >= 0:
        tier = "VUS"
    elif p >= -6:
        tier = "LB"
    else:
        tier = "B"
    return tier, p
