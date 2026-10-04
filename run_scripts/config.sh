#!/bin/bash
# Default settings of the run scripts. Change them on the command line,
# e.g. CHANNEL=Xn0n ./run_scripts/run_nucleus.sh

# Collision system and photon flux
: "${CHANNEL:=An0n}"               # An0n | Xn0n | PL(AnAn)
: "${FLUX_MODEL:=EFF}"             # EFF | PL | WS | TABLE (with FLUX_FILE=<file>)
: "${TARGET:=AA}"                  # AA (Pb+Pb, 5.36 TeV) | pA (p+Pb, 8.16 TeV, run_proton.sh)

# Fragmentation
: "${FRAG_TYPE:=KniehlKramer}"     # BCFY | KniehlKramer | HymnD
: "${HYMND_FILE:=input/HymnD/prompt-D0-1-109_0000.dat}"   # central HymnD member
: "${SCALE_FACTOR:=1.0}"           # fragmentation scale Q = SCALE_FACTOR * m_T
: "${SCALE_FACTORS:=0.5 1.0 2.0}"  # scale factors of run_scale_variation.sh

# Kinematic grid
: "${Y_VALS:=-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0}"
: "${PT_MIN:=0.1}"
: "${PT_STEP:=0.2}"
: "${PT_MAX:=12.0}"
: "${PT_VALS:=$(echo $(seq $PT_MIN $PT_STEP $PT_MAX))}"

# VEGAS calls per point
: "${CALLS:=2e5}"

# Folder for all results
: "${OUTPUT_ROOT:=output}"

# Parallel processes
: "${CORES:=$(( $(nproc) / 2 ))}"

export OUTPUT_ROOT CHANNEL FLUX_MODEL TARGET FRAG_TYPE HYMND_FILE SCALE_FACTOR CALLS
[[ -n "${FLUX_FILE:-}" ]] && export FLUX_FILE

# Tags used in folder and file names
channel_tag=$(echo "$CHANNEL" | tr -d '() ')
# "" for EFF, "_PL" etc. for the other fluxes
flux_tag=$([[ "$FLUX_MODEL" == "EFF" ]] && echo "" || echo "_${FLUX_MODEL}")
# "_G1" if GAMMA_AA_ONE is set
g1_tag=$([[ -n "${GAMMA_AA_ONE:-}" ]] && echo "_G1" || echo "")
