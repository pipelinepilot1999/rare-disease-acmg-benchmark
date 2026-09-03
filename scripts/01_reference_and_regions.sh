#!/usr/bin/env bash
# Reference prep + exome region definition (GENCODE v50 basic CDS +/- 20bp).
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate rdacmg
R=~/rdacmg/resources; W=~/rdacmg/work; mkdir -p $W
log(){ echo "[$(date +%H:%M:%S)] $*"; }

# --- reference: GIAB ships plain gzip; bgzip it so faidx works
if [[ ! -s $R/ref/GRCh38.fa.gz.fai ]]; then
  log "recompressing reference with bgzip"
  zcat $R/ref/GRCh38.fa.gz | bgzip -@ 8 -c > $R/ref/GRCh38.bgz.fa.gz
  mv $R/ref/GRCh38.bgz.fa.gz $R/ref/GRCh38.fa.gz
  samtools faidx $R/ref/GRCh38.fa.gz
fi
log "reference contigs: $(wc -l < $R/ref/GRCh38.fa.gz.fai)"

# --- exome regions: CDS +/- 20bp, autosomes + X, merged
if [[ ! -s $W/exome_regions.bed ]]; then
  log "building exome region set from GENCODE v50"
  zcat $R/gencode/gencode.v50.basic.annotation.gtf.gz \
   | awk -F'\t' 'BEGIN{OFS="\t"} $3=="CDS" && $1 ~ /^chr([0-9]+|X)$/ {
        s=$4-1-20; if(s<0)s=0; print $1, s, $5+20 }' \
   | sort -k1,1 -k2,2n \
   | bedtools merge -i - 2>/dev/null > $W/exome_regions.bed \
   || { zcat $R/gencode/gencode.v50.basic.annotation.gtf.gz \
        | awk -F'\t' 'BEGIN{OFS="\t"} $3=="CDS" && $1 ~ /^chr([0-9]+|X)$/ {
             s=$4-1-20; if(s<0)s=0; print $1, s, $5+20 }' \
        | sort -k1,1 -k2,2n \
        | awk 'BEGIN{OFS="\t"}{if(c==$1 && $2<=e){if($3>e)e=$3} else {if(c)print c,s,e; c=$1;s=$2;e=$3}} END{if(c)print c,s,e}' \
        > $W/exome_regions.bed; }
fi
log "exome regions: $(wc -l < $W/exome_regions.bed) intervals, $(awk '{n+=$3-$2}END{printf "%.1f Mb",n/1e6}' $W/exome_regions.bed)"

# --- trio high-confidence intersection
if [[ ! -s $W/trio_hiconf.bed ]]; then
  log "intersecting the three GIAB high-confidence BEDs"
  python - <<'PY'
import os
W=os.path.expanduser("~/rdacmg/work"); R=os.path.expanduser("~/rdacmg/resources")
def load(p):
    d={}
    for line in open(p):
        f=line.split()
        if len(f)<3 or not f[0].startswith("chr"): continue
        d.setdefault(f[0],[]).append((int(f[1]),int(f[2])))
    return {k:sorted(v) for k,v in d.items()}
def inter(a,b):
    out={}
    for c in set(a)&set(b):
        A,B=a[c],b[c]; i=j=0; r=[]
        while i<len(A) and j<len(B):
            s=max(A[i][0],B[j][0]); e=min(A[i][1],B[j][1])
            if s<e: r.append((s,e))
            if A[i][1]<B[j][1]: i+=1
            else: j+=1
        if r: out[c]=r
    return out
cur=load(f"{R}/giab/HG002.bed")
for s in ("HG003","HG004"): cur=inter(cur,load(f"{R}/giab/{s}.bed"))
ex=load(f"{W}/exome_regions.bed")
cur=inter(cur,ex)
tot=0
with open(f"{W}/trio_hiconf.bed","w") as fh:
    for c in sorted(cur, key=lambda x:(len(x),x)):
        for s,e in cur[c]:
            fh.write(f"{c}\t{s}\t{e}\n"); tot+=e-s
print(f"  trio hi-conf INTERSECT exome = {tot/1e6:.1f} Mb")
PY
fi
log "trio_hiconf: $(wc -l < $W/trio_hiconf.bed) intervals, $(awk '{n+=$3-$2}END{printf "%.1f Mb",n/1e6}' $W/trio_hiconf.bed)"

# --- how much of the exome is LOST to GIAB's difficult-region exclusions (a reported limitation)
python - <<'PY'
import os
W=os.path.expanduser("~/rdacmg/work")
def bp(p): return sum(int(l.split()[2])-int(l.split()[1]) for l in open(p) if l.startswith("chr"))
e=bp(f"{W}/exome_regions.bed"); t=bp(f"{W}/trio_hiconf.bed")
print(f"[exome coverage] GENCODE CDS+/-20 = {e/1e6:.1f} Mb ; testable (GIAB hi-conf) = {t/1e6:.1f} Mb "
      f"; EXCLUDED = {(e-t)/1e6:.1f} Mb ({100*(e-t)/e:.1f}%)")
open(f"{W}/exome_coverage.txt","w").write(f"{e}\t{t}\t{(e-t)/e:.4f}\n")
PY
log DONE
