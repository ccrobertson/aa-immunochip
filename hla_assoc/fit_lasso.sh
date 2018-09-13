#!/bin/bash



plink1.9 --bfile ../hla_imputation/imputation_results/AA_IMPUTED --assoc

plink1.9 --bfile ../hla_imputation/imputation_results/AA_IMPUTED --lasso 0.5