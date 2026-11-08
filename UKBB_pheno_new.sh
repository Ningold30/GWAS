#!/bin/bash

#########################################
# Parse arguments
#########################################

def_data_fld="41270"
def_UKBB_data="/reference/data/UKBB_500k/feb2024_pheno_data/ukb677634.tab"
def_output="$PWD"
def_name="endo_polyps"

pheno=()

while getopts ":p:d:u:o:n:" opt; do
  case $opt in
    p) pheno+=("$OPTARG") ;;
    d) data_fld="$OPTARG" ;;
    u) UKBB_data="$OPTARG" ;;
    o) output="$OPTARG" ;;
    n) name="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done

data_fld="${data_fld:-$def_data_fld}"
UKBB_data="${UKBB_data:-$def_UKBB_data}"
output="${output:-$def_output}"
name="${name:-$def_name}"

shift $((OPTIND -1))

if [[ ${#pheno[@]} -eq 0 ]]; then
  echo "Error: at least one phenotype code required (-p)" >&2
  exit 1
fi

# Build a label for filenames from all pheno codes (strip trailing '=' if used for exact match)
pheno_label=$(printf "%s_" "${pheno[@]//=/}" | sed 's/_$//')

echo "phenotype(s): ${pheno[*]}"
echo "data field: $data_fld"
echo "UKBB source: $UKBB_data"
echo "output label: ${pheno_label}"

#########################################
# Build grep pattern
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

GREP_PATTERN=$(build_pattern "${pheno[@]}")
echo "grep pattern: $GREP_PATTERN"

#########################################
# Set up working directories
#########################################

cd "${output}"
tmp_working="${output}/tmp_working"
mkdir -p "${tmp_working}"

# Base name for all output files
base="${name}_${pheno_label}_datafield_${data_fld}"

#########################################
# Extract cases
#########################################

range_1=$(head -1 "${UKBB_data}" | tr '\t' '\n' | grep -n "${data_fld}" | cut -d":" -f1 | head -1)
range_2=$(head -1 "${UKBB_data}" | tr '\t' '\n' | grep -n "${data_fld}" | cut -d":" -f1 | tail -1)

echo "searching columns $range_1 to $range_2"

cut -f1,${range_1}-${range_2} "${UKBB_data}" \
  | grep -E "$GREP_PATTERN" \
  > "${tmp_working}/${base}_cases.tmp"

echo "$(wc -l < "${tmp_working}/${base}_cases.tmp") cases found"

#########################################
# Extract controls
#########################################

# Controls: samples that don't match any of the pheno codes in this field
cut -f1,${range_1}-${range_2} "${UKBB_data}" \
  | grep -v -E "$GREP_PATTERN" \
  | grep -v "^f.eid" \
  > "${tmp_working}/${base}_controls.tmp"

echo "$(wc -l < "${tmp_working}/${base}_controls.tmp") controls found"

#########################################
# Check for overlaps
#########################################

echo "checking for overlapping samples between cases and controls:"
awk 'NR==FNR {a[$1]; next} $1 in a {count++} END {print count+0, "overlapping samples"}' \
  "${tmp_working}/${base}_cases.tmp" \
  "${tmp_working}/${base}_controls.tmp"

#########################################
# Merge cases and controls
#########################################

echo "merging cases and controls"

# Cases get label 1, controls get label 0
awk -v name="${name}" \
  'BEGIN {print "IID", name}
   {print $1, "1"}' \
  "${tmp_working}/${base}_cases.tmp" > "${tmp_working}/${base}_case_control.tmp"

awk '{print $1, "0"}' \
  "${tmp_working}/${base}_controls.tmp" >> "${tmp_working}/${base}_case_control.tmp"

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













