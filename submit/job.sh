#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/gpfs01/star/pwg/prozorov/dijets/pythia-jets"

# Arguments from condor:
#   $1 = ptHatMin
#   $2 = ptHatMax
#   $3 = nEvents
#   $4 = seed  (generated in the submit file)

PTMIN="${1:?ptHatMin missing}"
PTMAX="${2:?ptHatMax missing}"
NEVT="${3:?nEvents missing}"
SEED="${4:?seed missing}"


# Best-effort Cluster/Proc detection (don’t die if absent)
CLUSTER="${CLUSTER_ID:-${CLUSTER:-0}}"
PROC="${PROC_ID:-${PROC:-0}}"



# Optional: echo for quick debugging
echo "[`date`] Starting makeTree with: ptHatMin=$PTMIN ptHatMax=$PTMAX nEvents=$NEVT seed=$SEED"
echo "Hostname: $(hostname)"
echo "Cluster/Process: ${CLUSTER}/${PROC}"

# If apptainer isn't in PATH on the worker, try singularity as a fallback.
APPTAINER_BIN="$(command -v apptainer || true)"
if [[ -z "$APPTAINER_BIN" ]]; then
  APPTAINER_BIN="$(command -v singularity || true)"
fi
if [[ -z "$APPTAINER_BIN" ]]; then
  echo "Error: neither 'apptainer' nor 'singularity' found in PATH on the worker node." >&2
  exit 127
fi

IMG=$WORKDIR/rivet-pythia.sif

EXECUTABLE=$WORKDIR/makeTree

if [[ ! -x "$EXECUTABLE" ]]; then
  echo "Error: executable '$EXECUTABLE' not found " >&2
  exit 1
fi

OUTDIR=$WORKDIR/output
mkdir -p "$OUTDIR"

PREFIX="pp200_"$CLUSTER"_"$PROC

SEED="0"
# Run the job
# Note: we bind /gpfs01 because your inputs/outputs live there.
"$APPTAINER_BIN" exec -B /gpfs01 "$IMG" \
  "$EXECUTABLE" "$PTMIN" "$PTMAX" "$NEVT" "$SEED" $OUTDIR/$PREFIX

echo "[`date`] Finished."
