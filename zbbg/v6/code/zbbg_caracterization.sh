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
MACRO_DIR="${BASE_DIR}/code_zbbg"
HIST_BASE="${BASE_DIR}/Ntuples_zbbg/Histogram/GParT"
PLOT_BASE="${BASE_DIR}/Plots/zbbg/GParT"

# 1) Producer
echo "[1/5] Eseguo producer..."
root -l -q -b "${MACRO_DIR}/histogram_maker_plotter_categorization.C"

# 2) Hadd
echo "[2/5] Eseguo hadd..."
# Esegui lo script hadd.csh con csh usando il percorso assoluto e verifica che esista
if [ -f "${HIST_BASE}/hadd.csh" ]; then
  bash -c "cd /gwpool/users/fballo/HiggsAnalysis/code/boostedZgamma/Ntuples_zbbg/Histogram/GParT && source hadd.csh"
else
  echo "  [hadd] File non trovato: ${HIST_BASE}/hadd.csh" >&2
  exit 1
fi

# 3) Plotter
echo "[3/5] Eseguo plotter..."
root -l -q -b "${MACRO_DIR}/plotter_categorization.C"

# 4) Raccogli artefatti in versions/zbbg/<versione>/{hist,plot,code}/CAT*/
DEST_BASE="${REPO_BASE}/zbbg/${VERSION}"
mkdir -p "${DEST_BASE}"

echo "[4/5] Raccolgo output in ${DEST_BASE}"
for kind in hist plot; do
  for CAT in CAT0 CAT1 CAT2 CAT3; do
    src=""
    dst="${DEST_BASE}/${kind}/${CAT}"
    mkdir -p "${dst}"
    shopt -s nullglob
    case "${kind}" in
      hist)
        src="${HIST_BASE}/${CAT}"
        files=( "${src}"/*.root )
        ;;
      plot)
        src="${PLOT_BASE}/${CAT}"
        files=( "${src}"/*.png "${src}"/*.pdf )
        ;;
    esac
    if [ ${#files[@]} -eq 0 ]; then
      echo "  [${kind}/${CAT}] Nessun file da copiare da ${src}"
      continue
    fi
    echo "  [${kind}/${CAT}] Copio ${#files[@]} file"
    cp -a "${files[@]}" "${dst}/"
  done
done

# Scrivi META.txt con timestamp e commento (se fornito) nella cartella versione
TIMESTAMP="$(TZ=Europe/Rome date +'%Y-%m-%dT%H:%M:%S%z')"
{
  echo "timestamp: ${TIMESTAMP}"
  echo "version: ${VERSION}"
  echo "comment: ${COMMENT}"
} > "${DEST_BASE}/META.txt"

# 5) Salva snapshot dei codici usati in questa versione
echo "[5/5] Salvo snapshot del codice usato in questa versione..."
CODE_SNAP_DIR="${DEST_BASE}/code"
mkdir -p "${CODE_SNAP_DIR}"

# Salva macro e script correlati (lista minimale e chiara)
cp -a "${MACRO_DIR}/histogram_maker_plotter_categorization.C" "${CODE_SNAP_DIR}/" 2>/dev/null || true
cp -a "${MACRO_DIR}/plotter_categorization.C" "${CODE_SNAP_DIR}/" 2>/dev/null || true
# cp -a "${MACRO_DIR}/CMS_lumi.C" "${CODE_SNAP_DIR}/" 2>/dev/null || true
# cp -a "${MACRO_DIR}/CMS_lumi.h" "${CODE_SNAP_DIR}/" 2>/dev/null || true
# cp -a "${MACRO_DIR}/CustomColors.h" "${CODE_SNAP_DIR}/" 2>/dev/null || true
# cp -a "${BASE_DIR}/fit_code/prepareDatacardsZbb.py" "${CODE_SNAP_DIR}/" 2>/dev/null || true
# cp -a "${BASE_DIR}/fit_code/makeFitDistributionsZbb.py" "${CODE_SNAP_DIR}/" 2>/dev/null || true
# cp -a "${BASE_DIR}/fit_code/runWSCreation.sh" "${CODE_SNAP_DIR}/" 2>/dev/null || true
cp -a "$0" "${CODE_SNAP_DIR}/zbbg_caracterization.sh" 2>/dev/null || true

# Crea un archivio tar.gz del codice per portabilità
( cd "${DEST_BASE}" && tar -czf code_snapshot.tar.gz code )

echo "Fatto. Tutto in ${DEST_BASE}/{hist,plot,code}/CAT{n}"