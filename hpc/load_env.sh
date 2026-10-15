#!/usr/bin/env bash

set -euo pipefail


# Example environment configuration.
# Adapt module names and versions to your local HPC system.
# Additional system libraries may be required depending on local
# package compilation requirements.

module load gcc
module load openblas
module load jags

export TMPDIR="${TMPDIR:-/tmp}"
export RENV_CONFIG_PAK_ENABLED=TRUE