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

curl -O "ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/16S_Silva138_20200326.tgz" || exit 1
tar xzf 16S_Silva138_20200326.tgz
mv 16S_SILVA138_k2db 16S_Silva138
rm -f 16S_Silva138_20200326.tgz

curl -O "ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz" || exit 1
tar xzf minikraken_8GB_202003.tgz
mv minikraken_8GB_20200312 minikraken_8GB_2020-03-12
rm -f minikraken_8GB_202003.tgz
