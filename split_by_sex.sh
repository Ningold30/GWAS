#PBS -l mem=4gb,walltime=01:00:00,ncpus=4
#PBS -N UKBB_pheno
#PBS -e /working/lab_tracyo/nathanI/GWAS/endo_polyps/UKBB/scripts/merging_pheno.err
#PBS -o /working/lab_tracyo/nathanI/GWAS/endo_polyps/UKBB/scripts/merging_pheno.out

cd ${output}

#1 = women 0 = men

awk '{if (NR==1 || $24 ==1) print $0}' ${output}/${base}.covariates > ${base}.covariates.women

awk '{if (NR==1 || $24 ==0) print $0}' ${output}/${base}.covariates > ${base}.covariates.men

awk 'NR==FNR {a[$2];next} $2 in a {print $0}' ${base}.covariates.men ${base}.phenotype > ${base}.phenotype.men

awk 'NR==FNR {a[$2];next} $2 in a {print $0}' ${base}.covariates.women ${base}.phenotype > ${base}.phenotype.women
