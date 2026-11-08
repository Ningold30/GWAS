#!/bin/bash

#########################################
# Parse arguments
#########################################

def_data_fld="41270"
def_UKBB_data="/reference/data/UKBB_500k/feb2024_pheno_data/ukb677634.tab"
def_output="$PWD"
def_name="endo_polyps"

pheno=()
quant_mode=0
quant_threshold=""
quant_direction="ge"   # ge = >=threshold (cases), le = <=threshold (cases)

while getopts ":p:d:u:o:n:qt:r:" opt; do
  case $opt in
    p) pheno+=("$OPTARG") ;;
    d) data_fld="$OPTARG" ;;
    u) UKBB_data="$OPTARG" ;;
    o) output="$OPTARG" ;;
    n) name="$OPTARG" ;;
    q) quant_mode=1 ;;
    t) quant_threshold="$OPTARG" ;;
    r) quant_direction="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done

data_fld="${data_fld:-$def_data_fld}"
UKBB_data="${UKBB_data:-$def_UKBB_data}"
output="${output:-$def_output}"
name="${name:-$def_name}"

shift $((OPTIND -1))

#########################################
# Validate arguments
#########################################

if [[ "$quant_mode" -eq 1 ]]; then
  # In quantitative mode, -p codes are not needed
  # -d must be provided (or uses default, but warn user)
  echo "Running in QUANTITATIVE/PRE-BINARISED mode"
  echo "data field: $data_fld"
  echo "UKBB source: $UKBB_data"

  if [[ -n "$quant_threshold" ]]; then
    echo "threshold: values ${quant_direction} ${quant_threshold} = cases"
  else
    echo "no threshold set: using values as-is (field must already be 0/1 or numeric)"
  fi

  pheno_label="${data_fld}"

else
  # Standard ICD code mode: -p is required
  if [[ ${#pheno[@]} -eq 0 ]]; then
    echo "Error: at least one phenotype code required (-p), or use -q for quantitative mode" >&2
    exit 1
  fi

  # Build a label for filenames from all pheno codes (strip trailing '=' if used for exact match)
  pheno_label=$(printf "%s_" "${pheno[@]//=/}" | sed 's/_$//')

  echo "phenotype(s): ${pheno[*]}"
  echo "data field: $data_fld"
  echo "UKBB source: $UKBB_data"
  echo "output label: ${pheno_label}"
fi

#########################################
# Build grep pattern (ICD/code mode only)
#########################################

build_pattern() {
  local patterns=()
  for code in "$@"; do
    if [[ "$code" == *= ]]; then
      local c="${code%=}"
      patterns+=("\"${c}\"")   # exact match wrapped in quotes
    else
      patterns+=("\"${code}")  # prefix match - opening quote + code, anything after
    fi
  done
  local IFS="|"
  echo "${patterns[*]}"
}

if [[ "$quant_mode" -eq 0 ]]; then
  GREP_PATTERN=$(build_pattern "${pheno[@]}")
  echo "grep pattern: $GREP_PATTERN"
fi

#########################################
# Set up working directories
#########################################

cd "${output}"
tmp_working="${output}/tmp_working"
mkdir -p "${tmp_working}"

# Base name for all output files
base="${name}_${pheno_label}_datafield_${data_fld}"

#########################################
# Find column range for the data field
#########################################

range_1=$(head -1 "${UKBB_data}" | tr '\t' '\n' | grep -n "${data_fld}" | cut -d":" -f1 | head -1)
range_2=$(head -1 "${UKBB_data}" | tr '\t' '\n' | grep -n "${data_fld}" | cut -d":" -f1 | tail -1)

echo "searching columns $range_1 to $range_2"

#########################################
# Extract cases and controls
#########################################

if [[ "$quant_mode" -eq 1 ]]; then

  # -------------------------------------------------------
  # QUANTITATIVE / PRE-BINARISED MODE
  #
  # For fields like 23070 that already contain 0/1 values
  # (or any numeric value), extract all non-NA rows and
  # either use the values directly or apply a threshold.
  # -------------------------------------------------------

  echo "extracting non-missing values from field ${data_fld}"

  # Pull IID + the field columns; keep only rows with at least one non-NA value
  # (excludes the header line and rows where every instance column is empty/"")
  cut -f1,${range_1}-${range_2} "${UKBB_data}" \
    | awk 'NR==1 {next}                          # skip header
           {
             val = ""
             for (i=2; i<=NF; i++) {
               if ($i != "" && $i != "NA") { val = $i; break }
             }
             if (val != "") print $1, val
           }' \
    > "${tmp_working}/${base}_quant_raw.tmp"

  echo "$(wc -l < "${tmp_working}/${base}_quant_raw.tmp") samples with non-missing values"

  if [[ -n "$quant_threshold" ]]; then
    # Apply threshold to split into cases (1) and controls (0)
    awk -v thr="$quant_threshold" -v dir="$quant_direction" \
      '{
         if      (dir == "ge" && $2 >= thr) label = 1
         else if (dir == "le" && $2 <= thr) label = 1
         else if (dir == "gt" && $2 >  thr) label = 1
         else if (dir == "lt" && $2 <  thr) label = 1
         else if (dir == "eq" && $2 == thr) label = 1
         else label = 0
         print $1, label
       }' \
      "${tmp_working}/${base}_quant_raw.tmp" \
      > "${tmp_working}/${base}_case_control.tmp"

    n_cases=$(awk '$2==1' "${tmp_working}/${base}_case_control.tmp" | wc -l)
    n_controls=$(awk '$2==0' "${tmp_working}/${base}_case_control.tmp" | wc -l)
    echo "${n_cases} cases (value ${quant_direction} ${quant_threshold})"
    echo "${n_controls} controls"

  else
    # No threshold: treat the raw numeric value as the phenotype directly
    # (works for pre-binarised 0/1 fields like 23070)
    awk '{print $1, $2}' \
      "${tmp_working}/${base}_quant_raw.tmp" \
      > "${tmp_working}/${base}_case_control.tmp"

    n_cases=$(awk '$2==1' "${tmp_working}/${base}_case_control.tmp" | wc -l)
    n_controls=$(awk '$2==0' "${tmp_working}/${base}_case_control.tmp" | wc -l)
    echo "${n_cases} value=1 samples"
    echo "${n_controls} value=0 samples"
  fi

else

  # -------------------------------------------------------
  # STANDARD ICD / CODE MODE (original behaviour)
  # -------------------------------------------------------

  cut -f1,${range_1}-${range_2} "${UKBB_data}" \
    | grep -E "$GREP_PATTERN" \
    > "${tmp_working}/${base}_cases.tmp"

  echo "$(wc -l < "${tmp_working}/${base}_cases.tmp") cases found"

  # Controls: samples that don't match any of the pheno codes in this field
  cut -f1,${range_1}-${range_2} "${UKBB_data}" \
    | grep -v -E "$GREP_PATTERN" \
    | grep -v "^f.eid" \
    > "${tmp_working}/${base}_controls.tmp"

  echo "$(wc -l < "${tmp_working}/${base}_controls.tmp") controls found"

  # Check for overlaps
  echo "checking for overlapping samples between cases and controls:"
  awk 'NR==FNR {a[$1]; next} $1 in a {count++} END {print count+0, "overlapping samples"}' \
    "${tmp_working}/${base}_cases.tmp" \
    "${tmp_working}/${base}_controls.tmp"

  # Merge: cases get 1, controls get 0
  echo "merging cases and controls"
  awk -v name="${name}" \
    'BEGIN {print "IID", name}
     {print $1, "1"}' \
    "${tmp_working}/${base}_cases.tmp" > "${tmp_working}/${base}_case_control.tmp"

  awk '{print $1, "0"}' \
    "${tmp_working}/${base}_controls.tmp" >> "${tmp_working}/${base}_case_control.tmp"

fi

#########################################
# Filter for White British
#########################################

echo "filtering for White British and removing exclude samples"
awk 'NR==FNR {a[$2]; next} !($1 in a) {print}' \
  /reference/data/UKBB_500k/versions/lab_stuartma/exclude_samples/430k_whiteBrit_default.exclude.list \
  "${tmp_working}/${base}_case_control.tmp" \
  > "${tmp_working}/${base}_case_control_WB.tmp"

echo "$(wc -l < "${tmp_working}/${base}_case_control_WB.tmp") samples after WB filter"

#########################################
# Add covariates and write final outputs
#########################################

COVARIATES=/working/lab_tracyo/nathanI/reference_panel/UKBB_all_covariates.txt

echo "adding covariates"

# Covariate file
awk 'NR==FNR {a[$2]=$0; next} $1 in a {print $0, a[$1]}' \
  "${COVARIATES}" \
  "${tmp_working}/${base}_case_control_WB.tmp" \
  | cut -d " " -f3-4,12,14-34,38,39-41 \
  > "${output}/${base}.covariates"

# Phenotype file (FID IID PHENO format)
awk 'NR==FNR {a[$2]=$1; next} $1 in a {print a[$1], $0}' \
  "${COVARIATES}" \
  "${tmp_working}/${base}_case_control_WB.tmp" \
  > "${output}/${base}.phenotype"

#########################################
# Sex split
#########################################

echo "splitting by sex"
qsub -N splt -v output=${output},name=${name},base=${base} ~/Scripts/UKBB_phenotyping/split_by_sex.sh

timeout 5m bash -c '
while qstat -u nathanI | grep -q splt; do
    echo "sex split job is still running..."
    sleep 5
done
'

echo "$(wc -l < "${output}/${base}.phenotype.men") males"
echo "$(wc -l < "${output}/${base}.phenotype.women") females"
echo "done"
