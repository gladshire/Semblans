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

curl -OJL "https://uofi.box.com/shared/static/0rlp6q0u5uvc161mzbdr3c0xoiti63sk" || \
curl -OL "https://www.dropbox.com/scl/fi/n49jm9i1yscrfrsj1dnq8/ensembl_animals.pep.all.fa.gz?rlkey=yemush6bm36wr4fu7dpe8h5e0&st=98kgb83l&dl=1" || \
exit 1

curl -OJL "https://uofi.box.com/shared/static/lvg7x2qrxvg8hfcmue9xv4y9t1rgfb48" || \
curl -OL "https://www.dropbox.com/scl/fi/hbvtnd9wsiwt8k7gakcmq/ensembl_plant.pep.all.fa.gz?rlkey=8cp9sn5wrt9uu4vc2pmg8xvin&st=1gesuqq0&dl=1" || \
exit 1

curl -OJL "https://uofi.box.com/shared/static/qc4nmun4apb0pik5fqxukn4qvn3wm943" || \
curl -OL "https://www.dropbox.com/scl/fi/8as6tci331utcrqrl7txz/ensembl_fungi.pep.all.fa.gz?rlkey=eyhsv35lelnao5xgd7s9dy51e&st=oc9dz9s6&dl=1" || \
exit 1

gunzip ./*.gz
