#!/bin/bash
# Window-length test of the acute scenario (not in the report; see ../WINDOW_LENGTH_CT36EFSA.md).
# Reruns Rana and Ichthyosaura with a 7-day window and Podarcis with 4, 8 and 16-week windows.
# The window length d is changed in a temporary copy of the species folder; the repository
# files are not modified. Figures are copied here as <variant>_<figure>.png.
# Usage: bash run_window_length_CT36EFSA.sh   (OCTAVE=path/to/octave, default: octave)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
octave=${OCTAVE:-octave}
export QT_QPA_PLATFORM=offscreen
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

variant() {   # variant <name> <species> <sed expression>...
  local name=$1 sp=$2; shift 2
  mkdir -p "$tmp/$name"
  cp -r "$root/zoo_collection/$sp" "$tmp/$name/"
  local f="$tmp/$name/$sp/predict_${sp}_CT36EFSA.m"
  for e in "$@"; do sed -i "$e" "$f"; done
  CT36_SPDIR="$tmp/$name/$sp" "$octave" --no-gui -q "$root/zoo_collection/CT36EFSA_pipeline/run_octave_CT36EFSA.m" \
    "$sp" "$tmp/out_$name" acute > "$tmp/log_$name.txt" 2>&1
  Rscript "$root/zoo_collection/CT36EFSA_pipeline/plot_CT36EFSA.R" "$sp" "$tmp/out_$name" acute >> "$tmp/log_$name.txt" 2>&1
  for p in growth_combined growth_combined_relative combined_dynamics combined_dynamics_relative; do
    [ -f "$tmp/out_$name/$sp/${sp}_$p.png" ] && cp "$tmp/out_$name/$sp/${sp}_$p.png" "$here/${name}_$p.png"
  done
  grep "runs ok" "$tmp/log_$name.txt"
}

variant rana_7d Rana_temporaria 's/^  d = 7\*3\*6;/  d = 7;/' &
variant ichthyosaura_7d Ichthyosaura_alpestris 's/^  d_nest= 7\*4\*12;/  d_nest= 7;/' 's/^  d = 7\*4\*12;/  d = 7;/' &
variant podarcis_4w Podarcis_muralis 's/^  d=7;/  d=28;/' &
wait
variant podarcis_8w Podarcis_muralis 's/^  d=7;/  d=56;/' &
variant podarcis_16w Podarcis_muralis 's/^  d=7;/  d=112;/' &
wait
