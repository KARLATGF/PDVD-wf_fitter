#!/bin/bash
# Fit every channel of one pulser run, in parallel, then merge:
#   ./run_fit.sh <run> [n_jobs]      (default n_jobs = number of cores - 2)
# The input is magnify_files/magnify-<run>-<n>.root (next to this script).
# Outputs in data_files/ (root, txt, log) next to this script. See README.md.
cd "$(dirname "$0")" || exit 1
MAGNIFY_DIR="magnify_files"
RUN=$1
if [[ ! "$RUN" =~ ^[0-9]+$ ]]; then
    echo "Usage: $0 <run> [n_jobs]     (input: $MAGNIFY_DIR/magnify-<run>-<n>.root)"
    exit 1
fi
FILES=( "$MAGNIFY_DIR"/magnify-"$RUN"-*.root )
if [ ! -f "${FILES[0]}" ]; then
    echo "No $MAGNIFY_DIR/magnify-$RUN-*.root found"
    exit 1
fi
if (( ${#FILES[@]} > 1 )); then
    echo "More than one magnify file for run $RUN in $MAGNIFY_DIR (${FILES[*]}); keep only the one to fit"
    exit 1
fi
FILE=$(realpath "${FILES[0]}")
NJOBS=${2:-$(( $(nproc) - 2 ))}
(( NJOBS < 1 )) && NJOBS=1
NCH=6144
LOG_DIR="data_files/log/run_${RUN}"
mkdir -p "$LOG_DIR"

# Compile once (ACLiC, optimised); the jobs then only load the library.
root -l -b -q -e 'gSystem->SetFlagsOpt("-O2")' -e '.L wf_fitter_avg.C+' >| "$LOG_DIR/compile.log" 2>&1 \
    || { echo "Compilation failed, see $LOG_DIR/compile.log"; exit 1; }

rm -rf "data_files/parts/run_${RUN}" plots/fitted/*/Run_${RUN}
start=$(date +%s)
CHUNK=$(( (NCH + NJOBS - 1) / NJOBS ))
echo "Run $RUN ($FILE): $NJOBS jobs of ~$CHUNK channels"
for (( first=0; first<NCH; first+=CHUNK )); do
    last=$(( first + CHUNK - 1 )); (( last >= NCH )) && last=$(( NCH - 1 ))
    root -l -b -q -e '.L wf_fitter_avg.C+' -e "wf_fitter_avg(\"$FILE\", $RUN, $first, $last)" \
        >| "$LOG_DIR/ch${first}-${last}.log" 2>&1 &
done
wait
echo "Fits done in $(( $(date +%s) - start )) s; merging"
root -l -b -q -e 'gSystem->SetFlagsOpt("-O2")' -e ".x merge_run.C+($RUN)" 2>&1 | tee "$LOG_DIR/merge.log" | grep -v "^Info in"
