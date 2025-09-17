#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/gpfs01/star/pwg/prozorov/dijets/pythia-jets/submit"
SUBMIT=$WORKDIR/condor.submit
LIST=$WORKDIR/ptHatBins.list

cd $WORKDIR

mkdir -p $WORKDIR/log
# cleaning up old logs
rm -f $WORKDIR/log/pythia.*.log
rm -f $WORKDIR/log/pythia.*.err
rm -f $WORKDIR/log/pythia.*.out

# sanity: show how many lines we'll submit
N=$(grep -v '^\s*#' "$LIST" | grep -v '^\s*$' | wc -l | tr -d ' ')
echo "Submitting $N jobs from $LIST …"

i=0
grep -v '^\s*#' "$LIST" | grep -v '^\s*$' | while read -r PTMIN PTMAX NEVT _; do
  i=$((i+1))
  # Make a deterministic-per-line but unique seed (fits int32)
  # Change formula if you prefer purely random:
  RAW=$(( $(date +%s) + i*1117 ))
  SEED=$(( 1 + ( RAW % 900000000 ) ))

  echo "  [$i/$N] ptHat: $PTMIN..$PTMAX, nEvents: $NEVT, seed: $SEED"

  condor_submit \
    -append "arguments = ${PTMIN} ${PTMAX} ${NEVT} ${SEED}" \
    "$SUBMIT" >/dev/null
done


./condor_control.sh

# merge TTrees with hadd

echo "Merging TTrees..."

TREEDIR="/gpfs01/star/pwg/prozorov/dijets/pythia-jets/output"

i=0
grep -v '^\s*#' "$LIST" | grep -v '^\s*$' | while read -r PTMIN PTMAX NEVT _; do
  i=$((i+1))
  apptainer exec -B /gpfs01 /gpfs01/star/pwg/prozorov/dijets/pythia-jets/rivet-pythia.sif /usr/local/root/bin/hadd -f -j -k\
   $TREEDIR/sum_pp200_ptHat_${PTMIN}_${PTMAX}.root $TREEDIR/pp200_*_*_pThat_${PTMIN}_${PTMAX}.root
done

echo "Merging TTrees done."

echo "Running anaTrees..."
cd /gpfs01/star/pwg/prozorov/dijets/pythia-jets

apptainer exec -B /gpfs01 /gpfs01/star/pwg/prozorov/dijets/pythia-jets/rivet-pythia.sif\
 /usr/local/root/bin/root -l -b -q  anaTrees/anaTrees.cpp+
