R_OPTS=--no-restore --no-save
DAT=derived_data
TAB=tables
S=scripts

LOCAL=$local/AA_Immunochip/tables_and_figures 


all: $(TAB)/table1.doc $(TAB)/table1b.doc $(TAB)/hwe_plots.4d.pdf $(TAB)/table2.doc $(TAB)/table3.doc

$(TAB)/table1.doc: $(DAT)/genotypes_recode.raw.4d $(S)/multiallelic_assoc.R ../data_original/pc2.txt
	R CMD BATCH $(R_OPTS) '--args $< $@ ../data_original/pc2.txt' $(S)/multiallelic_assoc.R multiallelic_assoc.Rout
	scp $@ $(LOCAL) #Move target to desktop
	
$(TAB)/table1b.doc: $(DAT)/genotypes_recode.raw.4d $(S)/biallelic_assoc.R ../data_original/pc2.txt
	R CMD BATCH $(R_OPTS) '--args $< $@ ../data_original/pc2.txt' $(S)/biallelic_assoc.R biallelic_assoc.Rout
	scp $@ $(LOCAL) #Move target to desktop
		
$(TAB)/hwe_plots.4d.pdf: $(DAT)/genotypes_recode.raw.4d $(S)/hwe.R 
	R CMD BATCH '--args $< $@' $(S)/hwe.R hwe.Rout
	scp $@ $(LOCAL) #Move target to desktop
	
$(DAT)/genotypes_recode.raw: ../hla_imputation/imputation_results/AA_IMPUTED.bed
	plink2 -bfile ../hla_imputation/imputation_results/AA_IMPUTED --recode A --out $(DAT)/genotypes_recode

$(DAT)/genotypes_recode.raw.4d: $(DAT)/genotypes_recode.raw ../hla_imputation/imputation_results/AA_IMPUTED.bim
	bash $(S)/extract4d.sh

$(DAT)/phased_4d.txt: ../hla_imputation/imputation_results/AA_IMPUTED.bgl.phased
	bash $(S)/get_phased_4d.sh > $@

$(DAT)/haplotypes.4d.txt: $(DAT)/phased_4d.txt ../hla_imputation/imputation_results/AA_IMPUTED.fam
	perl $(S)/get_haplotypes.pl --in $< --out $@ --fam ../hla_imputation/imputation_results/AA_IMPUTED.fam

$(TAB)/table2.doc: $(S)/haplotype_assoc.R $(DAT)/haplotypes.4d.txt
	R CMD BATCH $(R_OPTS) $(S)/haplotype_assoc.R	
	scp $@ $(LOCAL) #Move target to desktop

$(TAB)/table3.doc: $(S)/multiallelic_conditional_assoc.R $(DAT)/genotypes_recode.raw.4d ../data_original/pc2.txt 
	R CMD BATCH $(R_OPTS) $(S)/multiallelic_conditional_assoc.R 
	scp $@ $(LOCAL)
	
	
#to automatically recreate a target each time, add .FORCE to its prerequisites ( and uncomment the target .FORCE below)    
#.FORCE:	
 
#to generate all targets run make with the -B option 
#make -B
