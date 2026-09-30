#!/usr/bin/env bash
# Section 10 assertion: ClinVar blinding must be verified by behaviour.
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate rdacmg
E=~/rdacmg/resources/exomiser/exomiser-cli-15.1.0
R=~/rdacmg/repo; O=~/rdacmg/blindtest; mkdir -p $O

CID=PMID_37054711_Family_10_Subject_13
VCF=~/rdacmg/cases/${CID}.pos.vcf.gz
HPO='HP:0011344,HP:0001270,HP:0000750,HP:0010864,HP:0002059,HP:0000280,HP:0000158,HP:0000154,HP:0005280,HP:0001252'
POS=1321072

# whitelist ON (data dir with clinvar db) vs OFF (same dir, flag off)
cat > $O/wl_on.properties <<EOF
exomiser.data-directory=${RDACMG_HOME:-$HOME/rdacmg}/resources/exomiser/data
exomiser.hg38.data-version=2512
exomiser.hg38.use-clinvar-white-list=true
exomiser.phenotype.data-version=2512
logging.level.org.monarchinitiative.exomiser=WARN
EOF
sed 's/white-list=true/white-list=false/' $O/wl_on.properties > $O/wl_off.properties

for c in wl_on wl_off; do
  python $R/scripts/make_job.py "$VCF" "$HPO" "${CID}.${c}" "$O" 1 > $O/${c}.yml
  java -Xmx16g -Dspring.config.location=file:$O/${c}.properties \
    -jar $E/exomiser-cli-15.1.0.jar analyse --analysis $O/${c}.yml > $O/${c}.log 2>&1
done

echo "=== causal variant (chr1:$POS, INTS11) under each condition"
printf "%-8s %-6s %-10s %-10s %-10s %-10s %s\n" cond rank combined pheno genevar variant whitelist
for c in wl_on wl_off; do
  awk -F'\t' -v p=$POS 'NR>1 && $16==p {printf "%-8s %-6s %-10s %-10s %-10s %-10s %s\n", "'"$c"'", $1, $7, $8, $9, $10, $12}' $O/${CID}.${c}.variants.tsv
done
echo
echo "=== candidate set size"
for c in wl_on wl_off; do echo -n "$c: "; tail -n +2 $O/${CID}.${c}.variants.tsv | wc -l; done
