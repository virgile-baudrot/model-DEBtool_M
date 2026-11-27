#!/bin/bash
# Rerun the CT36 EFSA simulations of the 9 species in GNU Octave, then draw the figures in R.
# Usage: bash run_all_CT36EFSA.sh [outroot] [scenarios]
#   outroot   default: zoo_collection/CT36EFSA_results
#   scenarios default: "acute full"  (acute = report section 3.3, full = APPENDIX C)
# Environment: OCTAVE (default: octave), JOBS (parallel Octave processes, default 3),
#              NO_PLOT=1 to skip the R figures.
# Takes about 30 min per scenario with 3 processes.
set -e
here=$(cd "$(dirname "$0")" && pwd)
out=${1:-$here/../CT36EFSA_results}
scenarios=${2:-"acute full"}
octave=${OCTAVE:-octave}
jobs=${JOBS:-3}
export QT_QPA_PLATFORM=offscreen
mkdir -p "$out"
species="Perdix_perdix Sorex_araneus Turdus_philomelos Motacilla_flava Alauda_arvensis Podarcis_muralis Crocidura_russula Rana_temporaria Ichthyosaura_alpestris"
for sc in $scenarios; do
  printf '%s\n' $species | xargs -P "$jobs" -I{} sh -c \
    "$octave --no-gui -q '$here/run_octave_CT36EFSA.m' {} '$out' $sc > '$out/log_{}_$sc.txt' 2>&1"
  [ "$NO_PLOT" = 1 ] || for s in $species; do Rscript "$here/plot_CT36EFSA.R" "$s" "$out" "$sc"; done
done
