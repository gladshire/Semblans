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

docker buildx build --tag semblans:latest .
