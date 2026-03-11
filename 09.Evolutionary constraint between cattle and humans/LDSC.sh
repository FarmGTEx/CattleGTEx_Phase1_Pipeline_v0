#!/usr/bin/env bash
set -euo pipefail


PROCESSED_DIR="/faststorage/project/farmgtex/zhudi/human_GWAS_ldsc"
BFILE_PREFIX="/home/dizhu/1000G_EUR_Phase3_plink/1000G.EUR.QC."
WLD_PREFIX="/home/dizhu/weights_hm3_no_hla/weights."
FRQ_PREFIX="/home/dizhu/1000G_Phase3_frq/1000G.EUR.QC."

BED_DIR="/faststorage/project/farmgtex/zhudi/run_ldsc/beds"

WORK_ANNOT="work/annot"    
WORK_L2="work/ldscores"    
MUNGED_DIR="gwas_munged"   
RESULTS_DIR="ldsc_results" 

mkdir -p "${WORK_ANNOT}" "${WORK_L2}" "${MUNGED_DIR}" "${RESULTS_DIR}"

mapfile -t BEDS < <(ls "${BED_DIR}"/*.bed | xargs -n1 basename | sed 's/\.bed$//' | sort)
if [[ ${#BEDS[@]} -eq 0 ]]; then
  echo "[ERROR] no such file BED：${BED_DIR}/*.bed" >&2
  exit 1
fi
echo "[INFO] BEDs: ${BEDS[*]}"

mapfile -t GWAS_FILES < <(ls "${PROCESSED_DIR}"/*.tsv.gz 2>/dev/null || true)
if [[ ${#GWAS_FILES[@]} -eq 0 ]]; then
  echo "[ERROR] no such file GWAS：${PROCESSED_DIR}/*.tsv.gz" >&2
  exit 1
fi
echo "[INFO] Processed GWAS files: ${#GWAS_FILES[@]}"

############################
# 1)
############################
echo "[STEP] make_annot.py -> ${WORK_ANNOT}"
for bed in "${BEDS[@]}"; do
  for c in {1..22}; do
    out_annot="${WORK_ANNOT}/${bed}.chr${c}.annot.gz"
    if [[ ! -s "${out_annot}" ]]; then
      python /home/dizhu/tools/ldsc/make_annot.py \
        --bed-file "${BED_DIR}/${bed}.bed" \
        --bimfile  "${BFILE_PREFIX}${c}.bim" \
        --annot-file "${out_annot}"
    fi
  done
done

############################
# 2) 
############################
echo "[STEP] ldsc.py --l2 -> ${WORK_L2}"
for bed in "${BEDS[@]}"; do
  for c in {1..22}; do
    out_pref="${WORK_L2}/${bed}.chr${c}"
    if [[ ! -s "${out_pref}.l2.ldscore.gz" ]]; then
      python /home/dizhu/tools/ldsc/ldsc.py \
        --l2 \
        --bfile "${BFILE_PREFIX}${c}" \
        --ld-wind-cm 1 \
	--thin-annot \
        --annot "${WORK_ANNOT}/${bed}.chr${c}.annot.gz" \
        --out   "${out_pref}"
    fi
  done
done

############################
# 3) 
############################
REF_LD_CHR=""

for bed in "${BEDS[@]}"; do
  if [[ -z "${REF_LD_CHR}" ]]; then
    REF_LD_CHR="${WORK_L2}/${bed}.chr"
  else
    REF_LD_CHR="${REF_LD_CHR},${WORK_L2}/${bed}.chr"
  fi
done
echo "[INFO] --ref-ld-chr = ${REF_LD_CHR}"

############################
# 4) 
############################
echo "[STEP] munge_sumstats.py -> ${MUNGED_DIR}"
for f in "${GWAS_FILES[@]}"; do
  trait=$(basename "${f}" .tsv.gz)
  out_pref="${MUNGED_DIR}/${trait}"
  if [[ ! -s "${out_pref}.sumstats.gz" ]]; then
    python /home/dizhu/tools/ldsc/munge_sumstats.py \
      --sumstats "${f}" \
      --out "${out_pref}" \
      --a1 alt \
      --a2 ref \
      --chr chr \
      --bp pos \
      --p P \
      --n-col N \
      --signed-sumstats beta_meta_hq,0 \
      --merge-alleles w_hm3.snplist
  fi
done

############################
# 5) 
############################
echo "[STEP] ldsc.py --h2 -> ${RESULTS_DIR}"
for sf in "${MUNGED_DIR}"/*.sumstats.gz; do
  trait=$(basename "${sf}" .sumstats.gz)
  out="${RESULTS_DIR}/${trait}__multiBED"

  CMD=(python /home/dizhu/tools/ldsc/ldsc.py
    --h2 "${sf}"
    --ref-ld-chr "${REF_LD_CHR}"
    --w-ld-chr   "${WLD_PREFIX}"
    --overlap-annot
    --print-coefficients
    --out "${out}"
  )

  echo "[RUN] ${trait}"
  "${CMD[@]}"
done


echo "[DONE] ：
- Munge ： ${MUNGED_DIR}/*.sumstats.gz
- LDSC：  ${RESULTS_DIR}/*__multiBED.log"
