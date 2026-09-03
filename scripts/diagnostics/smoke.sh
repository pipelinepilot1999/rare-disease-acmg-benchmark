#!/usr/bin/env bash
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate rdacmg
E=~/rdacmg/resources/exomiser/exomiser-cli-15.1.0
R=~/rdacmg/repo
O=~/rdacmg/smoke; mkdir -p $O
CASE=$(python -c "
import json;m=json.load(open('$R/data/manifest.json'))
c=[x for x in m if x['model']=='denovo'][0]
print(c['case_id'], c['vcf_pos'], ','.join(c['hpo']), c['gene'], c['truth'][0]['chrom']+':'+str(c['truth'][0]['pos']))
")
read -r CID VCF HPO GENE LOC <<< "$CASE"
echo "smoke case: $CID gene=$GENE at $LOC"
echo "$CID $GENE $LOC" > $O/case.txt

for cond in informed clinvar_blind; do
  python $R/scripts/make_job.py "$VCF" "$HPO" "${CID}.${cond}" "$O" 1 > $O/${cond}.yml
  echo "--- running $cond"
  /usr/bin/time -f "  wall=%es maxRSS=%MkB" java -Xmx24g \
    -Dspring.config.location=file:$R/config/${cond}.properties \
    -jar $E/exomiser-cli-15.1.0.jar analyse \
    --analysis $O/${cond}.yml \
    > $O/${cond}.log 2>&1 || { echo "FAILED $cond"; tail -25 $O/${cond}.log; exit 1; }
  tail -2 $O/${cond}.log
  ls -la $O/${CID}.${cond}* 2>/dev/null | head -3
done
