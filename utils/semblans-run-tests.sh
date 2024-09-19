#!/usr/bin/env bash
####################################################
# Suppress printing of error messages
# exec 2>/dev/null
# Stop on first error
set -o errexit
# Set trap on ERR to be inherited by shell functions
set -o errtrace
# Trap errors
trap 'echo Error at line: $LINENO' ERR
####################################################

mkdir -p semblans-tests-workdir || exit 1
cd semblans-tests-workdir || exit 1

DIR_EXAMPLES="/usr/local/share/semblans-examples"
DIR_ENSEMBLR="/usr/local/share/semblans-ensembl-refs"
DIR_KRAKENDB="/usr/local/share/semblans-kraken2-dbs"

semblans preprocess \
--verbose \
--threads 8 \
--ram 16 \
--retain \
--kraken-db ${DIR_KRAKENDB}/16S_Silva138,${DIR_KRAKENDB}/minikraken_8GB_2020-03-12 \
-1 ${DIR_EXAMPLES}/ChloroSubSet_1.fq \
-2 ${DIR_EXAMPLES}/ChloroSubSet_2.fq

semblans all \
--verbose \
--threads 8 \
--ram 16 \
--retain \
--kraken-db ${DIR_KRAKENDB}/16S_Silva138,${DIR_KRAKENDB}/minikraken_8GB_2020-03-12 \
--prefix assembly-prefix \
--ref-proteome ${DIR_ENSEMBLR}/ensembl_plant.pep.all.fa \
-1 ${DIR_EXAMPLES}/ChloroSubSet_1.fq \
-2 ${DIR_EXAMPLES}/ChloroSubSet_2.fq

semblans assemble \
--verbose \
--threads 8 \
--ram 16 \
--retain \
--prefix assembly-prefix \
-1 ${DIR_EXAMPLES}/ChloroSubSet_1.fq \
-2 ${DIR_EXAMPLES}/ChloroSubSet_2.fq
