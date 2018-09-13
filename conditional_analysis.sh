#!/bin/bash

#bash conditional_analysis.sh AA_icaa_hg19_noIBD 21 0.05

if [ $# == 3 ]; then
	study=$1
	chr=$2
	threshold=$3
else
	echo "Usage: conditional_analysis.sh <study> <chromosome> <significance threshold, alpha>"
	exit
fi

path=/m/jdrfdn_scratch/users/ccr5ju/AA_Immuno
outfile=$path/conditional_analysis/${study}.chr${chr}.conditional_analysis.txt
echo -e "CHR\tSNP\tcondition\tP" > $outfile


results=$path/assoc_results/plink.logistic.2PC.assoc.logistic.adjusted
condition_snps=condition_snp_list_chr${chr}

start_snp=$(awk -v x=$chr '$1==x {print $2}' $results | head -n1) #select top snp on specified chromosome
start_pval=$(awk -v x=$chr '$1==x {print $3}' $results | head -n1) #select pval associated with top snp on specified chromosome
echo -e "$chr\t$start_snp\t.\t$start_pval" >> $outfile

echo "start SNP: $start_snp"
newsnp=$start_snp
count=1
while [ $newsnp != "." ] ; do
	if [ $count -eq 1 ] 
	then
		echo $newsnp > $condition_snps  #create new file
	else
		echo $newsnp >> $condition_snps  #write next snp to file
	fi
	count+=1
	echo -e "Conditioning on $newsnp"
	file=${study}.chr${chr}.$newsnp
	plink2 --bfile $path/${study} --chr $chr --logistic sex --covar $path/pc2.txt --adjust --condition-list $condition_snps --out $file  
	awk -v x=$newsnp -v t=$threshold 'NR==2 {if ($3 < t) print $1, $2, x, $3}' ${file}.assoc.logistic.adjusted >> $outfile  		#if significant, print results to file
	
	#generate locuszoom plot
	#/m/cphg-expr1/cphg-expr1/locuszoom/bin/locuszoom --epacts results_2PCs_epactsformat.txt --markercol MARKER_ID --pvalcol PVALUE --chr 6 --start 29500000 --end 33500000 --build hg19 --pop AFR --source 1000G_March2012 --plotonly
	
	
	newsnp=$(awk -v t=$threshold 'NR==2 {if ($3 < t) {print $2} else {print "."}}'  ${file}.assoc.logistic.adjusted) 	#if significant, set newsnp to condition on
	pvalue=$(awk 'NR==2 {print $3}'  ${file}.assoc.logistic.adjusted)  
	echo -e
	echo !!!!!!!!!!!!!!! after conditioning, $newsnp is significant with a p-value of $pvalue !!!!!!!!!!!!!!!
	echo -e
	echo -e
done

echo "done"
