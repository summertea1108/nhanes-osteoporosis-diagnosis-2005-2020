#!/bin/bash
# Download every public NHANES file used in this analysis into ./nhanes_data
# Usage:  bash download_nhanes.sh        (needs internet; on OSC run it on a login node)
set -euo pipefail
DATA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/nhanes_data"
mkdir -p "$DATA_DIR"
BASE="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"
fetch() {  # $1 = file name without .xpt, $2 = cycle start year
  local dest="$DATA_DIR/$1.xpt"
  [[ -s "$dest" ]] && { echo "have $1"; return; }
  echo "get  $1"; curl -fsSL --retry 3 -o "$dest" "$BASE/$2/DataFiles/$1.xpt"
}
# 2005-2006 (D), 2007-2008 (E), 2009-2010 (F), 2013-2014 (H)
for pair in D:2005 E:2007 F:2009 H:2013; do
  s=${pair%%:*}; y=${pair##*:}
  for f in DEMO DXXFEM OSQ RXQ_RX HIQ HUQ BMX GHB DIQ BPQ BPX TCHOL; do fetch "${f}_${s}" "$y"; done
done
# 2017-March 2020 prepandemic file (P_ prefix; oscillometric blood pressure BPXO)
for f in DEMO DXXFEM OSQ RXQ_RX HIQ HUQ BMX GHB DIQ BPQ BPXO TCHOL; do fetch "P_${f}" 2017; done
echo "Done: $(ls "$DATA_DIR"/*.xpt | wc -l) files in $DATA_DIR"
