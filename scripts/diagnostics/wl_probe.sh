#!/usr/bin/env bash
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate rdacmg
E=$HOME/rdacmg/resources/exomiser/exomiser-cli-15.1.0
O=$HOME/rdacmg/blindtest
cat > $O/probe.properties <<EOF
exomiser.data-directory=$HOME/rdacmg/resources/exomiser/data
exomiser.hg38.data-version=2512
exomiser.hg38.use-clinvar-white-list=true
exomiser.phenotype.data-version=2512
logging.level.org.monarchinitiative.exomiser=INFO
EOF
java -Xmx16g -Dspring.config.location=file:$O/probe.properties \
  -jar $E/exomiser-cli-15.1.0.jar analyse --analysis $O/wl_on.yml > $O/probe.log 2>&1
echo "=== whitelist / clinvar lines:"
grep -iE "whitelist|white-list|white list|clinvar" $O/probe.log | head -12
echo "=== all INFO lines mentioning variants loaded:"
grep -iE "variants|filter" $O/probe.log | head -14
