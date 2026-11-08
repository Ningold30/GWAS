#!/bin/bash

#########################################
#Parse arguments
#########################################

# define default args
def_pheno="N840"
def_data_fld="41270"
def_UKBB_data="/reference/data/UKBB_500k/feb2024_pheno_data/ukb677634.tab"
def_output="$PWD"
def_name="endo_polyps"
# parse args
while getopts ":p:d:u:o:n:" opt; do
  case $opt in
    p) pheno="$OPTARG"
    ;;
    d) data_fld="$OPTARG"
    ;;
    u) UKBB_data="$OPTARG"
    ;;
    o) output="$OPTARG"
    ;;
    n) name="$OPTARG"
    ;;
    \?) echo "Invalid option: -$OPTARG" >&2
    ;;
  esac
done

# Set default values if not provided
pheno="${pheno:-$def_pheno}"
data_fld="${data_fld:-$def_data_fld}"
UKBB_data="${UKBB_data:-$def_UKBB_data}"
output="${output:-$def_output}"
name="${name:-$def_name}"
# Shift option index so $1 now refers to the first non-option argument
shift $((OPTIND -1))


# echo parsed flags
echo "phenotype is: $pheno"
echo "data field is: $data_fld"
echo "using $UKBB_data as ukbb source data"

#######################################
#Pick out phenotype
#######################################

cd ${output}
echo "using directory $PWD to save phenotype file"
tmp_working=${output}/tmp_working
mkdir -p ${tmp_working}

#cut out data field
range_1=$(head -1 ${UKBB_data} | tr '\t' '\n' | grep -n ${data_fld} | cut -d":" -f1 | head -1)
range_2=$(head -1 ${UKBB_data} | tr '\t' '\n' | grep -n ${data_fld} | cut -d":" -f1 | tail -1)

#Qsub job to get controls
qsub -N ${pheno} -v range_1=${range_1},range_2=${range_2},pheno=${pheno},data_fld=${data_fld},UKBB_data=${UKBB_data},tmp_working=${tmp_working},name=${name} -e /working/lab_tracyo/nathanI/GWAS/endo_polyps/UKBB/scripts/picking_out_controls.err -o /working/lab_tracyo/nathanI/GWAS/endo_polyps/UKBB/scripts/picking_out_controls.out  -l mem=4gb,walltime=02:00:00,ncpus=4 ~/Scripts/UKBB_phenotyping/UKBB_controls.sh
echo "qsub submitted to extract controls"
qstat -u nathanI

echo "extracting cases"
echo "searching for phenotype between columns $range_1 and $range_2 "

cut -f1,${range_1}-${range_2} ${UKBB_data} | grep -w ${pheno} > ${tmp_working}/UKBB_${name}_${pheno}_datafield_${data_fld}.tmp
timeout 5m bash -c '
while qstat -u nathanI | grep -q ${pheno} ; do     
   echo "controls job is still running..."
   sleep 5
done
'
wc -l ${tmp_working}/UKBB_${name}_${pheno}_datafield_${data_fld}.tmp
wc -l ${tmp_working}/UKBB_${name}_${pheno}_datafield_${data_fld}_controls.tmp

########################################
#adding covariates and filtering
########################################

echo -e "checking there are no overlapping samples between cases and controls \noverlapping samples:"

awk 'NR=FNR {a[$1];next} $1 in a {print c++} END {print c}' ${tmp_working}/UKBB_${name}_${pheno}_datafield_${data_fld}.tmp ${tmp_working}/UKBB_${name}_${pheno}_datafield_${data_fld}_controls.tmp

echo "merging cases and controls into one file"

awk -v name=${name} '{if (NR==1) print "IID",name ; else print $1,"1"}' ${tmp_working}/UKBB_${name}_${pheno}_datafield_${data_fld}.tmp | cat - ${tmp_working}/UKBB_${name}_${pheno}_datafield_${data_fld}_controls.tmp | awk -v name=${name} '{if ($2 == name || $2 == 1) print $0; else print $1, "0"}' | grep -v "f.eid" >  ${tmp_working}/UKBB_${name}_${pheno}_datafield_${data_fld}_case_control.tmp

echo 'merged saved here $tmp_working/UKBB_$name_$pheno_datafield_$data_fld_case_control.tmp'

#filter for WB
echo 'filtering for white british'
awk 'NR==FNR {a[$2];next} (!($1 in a)){print $0}' /reference/data/UKBB_500k/versions/lab_stuartma/exclude_samples/430k_whiteBrit_default.exclude.list ${tmp_working}/UKBB_${name}_${pheno}_datafield_${data_fld}_case_control.tmp > ${tmp_working}/UKBB_${name}_${pheno}_datafield_${data_fld}_case_control_WB.tmp


#add covariates
awk 'NR==FNR {a[$2]=$0;next} $1 in a {print $0,a[$1]}' /working/lab_tracyo/nathanI/reference_panel/UKBB_all_covariates.txt ${tmp_working}/UKBB_${name}_${pheno}_datafield_${data_fld}_case_control_WB.tmp | cut -d " " -f3-4,12,14-34,38,39-41 > ${output}/UKBB_${name}_${pheno}_datafield_${data_fld}.covariates

echo 'covariate file saved here $output/UKBB_$name_$pheno_datafield_$data_fld.covariates'

awk 'NR==FNR {a[$2]=$1;next} $1 in a {print a[$1],$0}' /working/lab_tracyo/nathanI/reference_panel/UKBB_all_covariates.txt ${tmp_working}/UKBB_${name}_${pheno}_datafield_${data_fld}_case_control_WB.tmp > ${output}/UKBB_${name}_${pheno}_datafield_${data_fld}.phenotype

echo 'phenotype file saved here $output/UKBB_$name_$pheno_datafield_$data_fld.phenotype'

#merge phenotype files
echo 'merging phenotpe files'
qsub -N mrge -v output=${output},name=${name},pheno=${pheno} ~/Scripts/UKBB_phenotyping/merging_pheno.sh

timeout 5m bash -c '
while qstat -u nathanI | grep -q mrge ; do     
   echo "mergin results..."
   sleep 10
done
'
#split by sex
echo 'spliting phenotpe files by sex'
qsub -N splt -v output=${output},name=${name} ~/Scripts/UKBB_phenotyping/split_by_sex.sh

qstat -u nathanI

echo 'done'








