#!/bin/bash
# Resolve SpecHLA paths for both conda and development installs.
#
# Priority chain:
#   1. SPECHLA_DB env var (user override)
#   2. $CONDA_PREFIX/share/spechla/db (conda install)
#   3. Relative to this script (development install)
#
# Source this from other bash scripts:
#   source "$(dirname $(realpath $0))/spechla_env.sh"   # from script/
#   source "$(dirname $(realpath $0))/../spechla_env.sh" # from script/whole/

# Guard against repeated sourcing (prevents PATH pollution in nested scripts)
if [ -n "${_SPECHLA_ENV_SOURCED:-}" ]; then
    return 0 2>/dev/null || true
fi
_SPECHLA_ENV_SOURCED=1

if [ -z "${SPECHLA_DB:-}" ] || [ ! -d "$SPECHLA_DB" ]; then
    if [ -n "${CONDA_PREFIX:-}" ] && [ -d "$CONDA_PREFIX/share/spechla/db" ]; then
        SPECHLA_DB="$CONDA_PREFIX/share/spechla/db"
    else
        _spechla_root=$(cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/.." && pwd)
        SPECHLA_DB="$_spechla_root/db"
    fi
fi

if [ -z "${SPECHLA_SCRIPT:-}" ] || [ ! -d "$SPECHLA_SCRIPT" ]; then
    if [ -n "${CONDA_PREFIX:-}" ] && [ -d "$CONDA_PREFIX/share/spechla/script" ]; then
        SPECHLA_SCRIPT="$CONDA_PREFIX/share/spechla/script"
    else
        SPECHLA_SCRIPT=$(cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")" && pwd)
    fi
fi

# In dev mode, add the repo's bin/ to PATH so vendored binaries are found.
# In conda mode, conda already provides these tools on PATH.
if [ -z "${CONDA_PREFIX:-}" ] || [ ! -d "$CONDA_PREFIX/share/spechla" ]; then
    _spechla_bin=$(cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/../bin" && pwd)
    if [ -d "$_spechla_bin" ]; then
        export PATH="$_spechla_bin:$_spechla_bin/SpecHap/build:$_spechla_bin/extractHairs/build:$_spechla_bin/fermikit/fermi.kit:$PATH"
    fi
fi

export SPECHLA_DB
export SPECHLA_SCRIPT
