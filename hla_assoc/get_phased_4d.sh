#!/bin/bash

awk --posix 'NR <=5 || $2~/^HLA_[^_]+_([0-9]{3})[0-9]$/' ../hla_imputation/imputation_results/AA_IMPUTED.bgl.phased
