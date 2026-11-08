#PBS -l mem=8gb,walltime=04:00:00,ncpus=4
#PBS -N UKBB_pheno
#PBS -e /working/lab_tracyo/nathanI/GWAS/endo_polyps/UKBB/scripts/merging_pheno.err
#PBS -o /working/lab_tracyo/nathanI/GWAS/endo_polyps/UKBB/scripts/merging_pheno.out

cd ${output}

cat UKBB_${name}*.phenotype | sort -nrk3 | awk '!seen[$2]++' | sort -g > ${name}_UKBB_WB.phenotype.merged

cat UKBB_${name}*.covariates | awk '!seen[$2]++' | awk 'NR==FNR {a[$2];next} $2 in a {print $0}' ${name}_UKBB_WB.phenotype.merged - > ${name}_UKBB_WB.covariates.merged
