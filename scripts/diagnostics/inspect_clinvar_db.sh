#!/usr/bin/env bash
set -uo pipefail
H2=$HOME/rdacmg/resources/exomiser/exomiser-cli-15.1.0/lib/h2-2.2.224.jar
cp -f $HOME/rdacmg/resources/exomiser/data/2512_hg38/2512_hg38_clinvar.mv.db /tmp/cvq.mv.db
chmod u+w /tmp/cvq.mv.db
q(){ java -cp $H2 org.h2.tools.Shell -url "jdbc:h2:file:/tmp/cvq" -user sa -password "" -sql "$1" 2>&1; }
echo "=== tables"
q "SELECT TABLE_SCHEMA, TABLE_NAME FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA <> 'INFORMATION_SCHEMA';" | head -8
echo "=== columns of CLINVAR"
q "SELECT COLUMN_NAME, DATA_TYPE FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_NAME='CLINVAR';" | head -20
echo "=== row count"
q "SELECT COUNT(*) FROM CLINVAR;" | head -4
echo "=== the INTS11 variant chr1:1321072"
q "SELECT * FROM CLINVAR WHERE CONTIG=1 AND START=1321072;" | head -6
