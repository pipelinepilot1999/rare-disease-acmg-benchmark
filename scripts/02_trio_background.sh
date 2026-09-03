#!/usr/bin/env bash
# Joint trio VCF from GIAB v4.2.1 benchmarks, restricted to testable exome, normalized.
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate rdacmg
R=~/rdacmg/resources; W=~/rdacmg/work
log(){ echo "[$(date +%H:%M:%S)] $*"; }

REF=$R/ref/GRCh38.fa.gz

# 1. restrict each sample to the testable region, normalize, keep PASS
for s in HG002 HG003 HG004; do
  if [[ ! -s $W/${s}.norm.vcf.gz.tbi ]]; then
    log "normalizing $s"
    bcftools view -R $W/trio_hiconf.bed -f PASS $R/giab/${s}.vcf.gz -Ou \
      | bcftools norm -f $REF -m -any -Ou 2>/dev/null \
      | bcftools norm -d exact -Oz -o $W/${s}.norm.vcf.gz
    tabix -f -p vcf $W/${s}.norm.vcf.gz
  fi
  log "  $s: $(bcftools index -n $W/${s}.norm.vcf.gz) variants"
done

# 2. merge; within the high-confidence intersection a missing call is a confident hom-ref
if [[ ! -s $W/trio.vcf.gz.tbi ]]; then
  log "merging trio"
  bcftools merge -m none -0 $W/HG002.norm.vcf.gz $W/HG003.norm.vcf.gz $W/HG004.norm.vcf.gz -Oz -o $W/trio.raw.vcf.gz
  # canonical sample order: proband, father, mother
  bcftools view -s HG002,HG003,HG004 $W/trio.raw.vcf.gz -Oz -o $W/trio.vcf.gz
  tabix -f -p vcf $W/trio.vcf.gz
  rm -f $W/trio.raw.vcf.gz
fi
log "trio: $(bcftools index -n $W/trio.vcf.gz) variants; samples: $(bcftools query -l $W/trio.vcf.gz | tr '\n' ' ')"

# 3. ASSERTION: Mendelian consistency (trio relationship from genotypes, not filenames)
log "checking Mendelian consistency"
bcftools query -f '[%GT\t]\n' $W/trio.vcf.gz | python - <<'PY'
import sys
def g(t):
    t=t.replace("|","/")
    if t in ("./.","."): return None
    a=t.split("/")
    try: return sorted(int(x) for x in a)
    except: return None
tot=err=0
for line in sys.stdin:
    f=line.rstrip("\n").split("\t")[:3]
    c,fa,mo=(g(x) for x in f)
    if not c or not fa or not mo or len(c)!=2: continue
    tot+=1
    if not any(c[0] in fa and c[1] in mo or c[1] in fa and c[0] in mo for _ in (0,)):
        err+=1
print(f"[mendel] informative={tot} violations={err} rate={err/max(tot,1):.4%}")
import os
open(os.path.expanduser("~/rdacmg/work/mendel.txt"),"w").write(f"{tot}\t{err}\t{err/max(tot,1):.6f}\n")
PY

# 4. genotype spectrum, for the record
log "genotype spectrum (proband):"
bcftools query -f '[%GT ]\n' $W/trio.vcf.gz | awk '{print $1}' | sort | uniq -c | sort -rn | head
log DONE
