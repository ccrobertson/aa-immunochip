#!/bin/bash




### Get race codes
cd /m/jdrfdn_scratch/users/ccr5ju/AA_Immuno/processed
awk -F "," '{print $5, $6, $18, $20}' /m/jdrfdn_scratch/users/projects/IMCHIP/T1DGC.2011.03_Golden_Pedigree/T1DGC.2011.03_Golden_Pedigree/Resources/T1DGC.2011.03_Resources.csv > T1DGC+resources_ID+RaceCode.txt
#Update t1dgc famids in race code data
awk 'NR==FNR {a[$2]=$1} NR!=FNR && $2 in a {$1=a[$2]; print $0} NR!=FNR && !($2 in a) {print}' ../hla_imputation/mkreference/familyrs.fam T1DGC+resources_ID+RaceCode.txt > T1DGC+resources_ID+RaceCode_updated+t1dgc+famids.txt



### Race distribution in reference panel
cd /m/jdrfdn_scratch/users/ccr5ju/AA_Immuno/hla_imputation
awk 'NR==FNR {a[$1"\t"$2]=$5"\t"$6; next} NR!=FNR && ($1"\t"$2 in a) {print $0, a[$1"\t"$2]}' mkreference/t1dgc_refpanel.fam ../processed/T1DGC+resources_ID+RaceCode_updated+t1dgc+famids.txt > t1dgc_refpanel_RaceCodes.txt

#R code:
df = read.table("t1dgc_refpanel_RaceCodes.txt")
names(df) <- c("FAMID", "IID", "Race1", "Race2","Sex","Pheno")

table(df$Race1)
#1=American Indian/ Alaskan Native 
#2=Asian 
#3=Native Hawaiian or Other Pacific Islander 
#4=Black or African American 
#5=White or Caucasian 
#6=Unknown

table(df$Sex)

table(df$Pheno)

table(df$Race1, df$Pheno)