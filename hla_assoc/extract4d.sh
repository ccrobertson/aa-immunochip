#!/bin/bash


echo -e "FID\nIID\nPAT\nMAT\nSEX\nPHENOTYPE" > derived_data/allele_list.4d
awk --posix '$2~/^HLA_[^_]+_[0-9]{3}[0-9]$/ {print $2"_P"}' ../hla_imputation/imputation_results/AA_IMPUTED.bim >> derived_data/allele_list.4d
perl ~/utility_scripts/extractCols.pl --file derived_data/genotypes_recode.raw --out derived_data/genotypes_recode.raw.4d --cols derived_data/allele_list.4d --cols-file

