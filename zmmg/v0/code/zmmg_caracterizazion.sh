#!/usr/bin/env bash
set -euo pipefail

# Argomenti: versione obbligatoria, commento opzionale
if [ $# -lt 1 ]; then
  echo "Usage: $0 <versione> [commento opzionale]"
  exit 1
fi
VERSION="$1"; shift || true
COMMENT="${*:-}"

BASE_DIR="/gwpool/users/fballo/HiggsAnalysis/code/boostedZgamma"
REPO_BASE="/gwpool/users/fballo/HiggsAnalysis/code/versions"
MACRO_DIR="${BASE_DIR}/code_zmmg"
HIST_BASE="${BASE_DIR}/Ntuples_zmmg/Histogram"
PLOT_BASE="${BASE_DIR}/Plots/zmmg"

echo "[1/4] Eseguo producer zmmg..."
root -l -q -b "${MACRO_DIR}/histogram_maker_plotter_zmumug.C"

echo "[2/4] Eseguo plotter zmmg..."
root -l -q -b "${MACRO_DIR}/plotter_ZmmG.C"

DEST_BASE="${REPO_BASE}/zmmg/${VERSION}"
mkdir -p "${DEST_BASE}"
echo "[3/4] Raccolgo output in ${DEST_BASE}"

# Histograms
mkdir -p "${DEST_BASE}/hist"
shopt -s nullglob
hist_files=( "${HIST_BASE}"/*.root )
if [ ${#hist_files[@]} -gt 0 ]; then
  echo "  [hist] Copio ${#hist_files[@]} file"
  cp -a "${hist_files[@]}" "${DEST_BASE}/hist/"
else
  echo "  [hist] Nessun file ROOT in ${HIST_BASE}"
fi

# Plots
mkdir -p "${DEST_BASE}/plot"
plot_files=( "${PLOT_BASE}"/*.png "${PLOT_BASE}"/*.pdf )
plot_count=0
for f in "${plot_files[@]}"; do
  if [ -f "$f" ]; then plot_count=$((plot_count+1)); fi
done
if [ $plot_count -gt 0 ]; then
  echo "  [plot] Copio ${plot_count} file"
  cp -a ${PLOT_BASE}/*.png ${PLOT_BASE}/*.pdf "${DEST_BASE}/plot/" 2>/dev/null || true
else
  echo "  [plot] Nessun plot in ${PLOT_BASE}"
fi

# META
TIMESTAMP="$(TZ=Europe/Rome date +'%Y-%m-%dT%H:%M:%S%z')"
{
  echo "timestamp: ${TIMESTAMP}"
  echo "version: ${VERSION}"
  echo "comment: ${COMMENT}"
} > "${DEST_BASE}/META.txt"


echo "[4/4] Salvo snapshot del codice usato in questa versione..."
CODE_SNAP_DIR="${DEST_BASE}/code"
mkdir -p "${CODE_SNAP_DIR}"
cp -a "${MACRO_DIR}/histogram_maker_plotter_zmumug.C" "${CODE_SNAP_DIR}/" 2>/dev/null || true
cp -a "${MACRO_DIR}/plotter_ZmmG.C" "${CODE_SNAP_DIR}/" 2>/dev/null || true
# cp -a "${MACRO_DIR}/CMS_lumi.C" "${CODE_SNAP_DIR}/" 2>/dev/null || true
# cp -a "${MACRO_DIR}/CMS_lumi.h" "${CODE_SNAP_DIR}/" 2>/dev/null || true
# cp -a "${MACRO_DIR}/CustomColors.h" "${CODE_SNAP_DIR}/" 2>/dev/null || true
cp -a "$0" "${CODE_SNAP_DIR}/zmmg_caracterizazion.sh" 2>/dev/null || true
( cd "${DEST_BASE}" && tar -czf code_snapshot.tar.gz code && mv code_snapshot.tar.gz "${CODE_SNAP_DIR}/" )

echo "Fatto. Tutto in ${DEST_BASE}/{hist,plot,code}"
