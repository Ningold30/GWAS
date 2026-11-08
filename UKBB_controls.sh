#!/bin/bash

cut -f1,${range_1}-${range_2} ${UKBB_data} | grep -vw ${pheno} > ${tmp_working}/UKBB_${name}_${pheno}_datafield_${data_fld}_controls.tmp

echo " $pheno $data_fld controls saved here $tmp_working/UKBB_$name_$pheno_datafield_$data_fld_controls.tmp" 
