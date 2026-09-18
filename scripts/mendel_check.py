#!/usr/bin/env python3
"""Confirm the trio relationship from genotypes, not from filenames.

Reads 'proband<TAB>father<TAB>mother' genotype triples on stdin.
"""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import _paths
import sys, os

def gt(t):
    t = t.replace("|", "/")
    if t in ("./.", ".", ""):
        return None
    try:
        return sorted(int(x) for x in t.split("/"))
    except ValueError:
        return None

total = err = 0
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    if len(f) < 3:
        continue
    child, father, mother = gt(f[0]), gt(f[1]), gt(f[2])
    if not child or not father or not mother or len(child) != 2:
        continue
    total += 1
    # a child allele pair is consistent if one allele can come from each parent
    ok = ((child[0] in father and child[1] in mother) or
          (child[1] in father and child[0] in mother))
    if not ok:
        err += 1

rate = err / max(total, 1)
print(f"[mendel] informative={total} violations={err} rate={rate:.4%}")
out = _paths.s(_paths.WORK / "mendel.txt")
with open(out, "w") as fh:
    fh.write(f"{total}\t{err}\t{rate:.6f}\n")
if total == 0:
    sys.exit("ASSERTION FAILED: no informative trios parsed")
if rate > 0.02:
    sys.exit(f"ASSERTION FAILED: Mendelian violation rate {rate:.2%} > 2% "
             f"- the three samples may not be a trio")
print("[mendel] OK")
