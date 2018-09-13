options(stringsAsFactors=FALSE)
setwd("/m/jdrfdn_scratch/users/ccr5ju/AA_Immuno/hla_assoc")
#install.packages("lmtest", repos="https://cloud.r-project.org")
library(lmtest)
library(rtf)

genofile = "derived_data/genotypes_recode.raw.4d"
rtffilename = "tables/table3.doc"
phenofile ="../data_original/pc2.txt"

rtffile <- RTF(rtffilename)  # this can be an .rtf or a .doc
geno = read.table(genofile, comment.char="~", header=TRUE)
pheno = read.table(phenofile, comment.char="~", header=TRUE)
covariates = names(pheno)[c(-1,-2)]
df = merge(geno, pheno, by="IID")
df$t1d = df$PHENOTYPE-1

getCondPvals = function(genes, condition_genes) {
	addParagraph(rtffile, paste("Conditioning on:", paste(condition_genes,collapes="\t"),"\n", sep=" "))
	for (k in 1:length(genes)) {
		gene = genes[k]
		
		condition_alleles = NULL
		for(f in 1:length(condition_genes)) {
			condition_alleles = c(condition_alleles, names(geno)[grep(condition_genes[f], names(geno))])
		}	
		predictor_alleles = names(geno)[grep(gene, names(geno))]
		
		#Quality filtering (exclude subjects whose alleles do not add to 2 for a given locus)
		X = geno[,c(predictor_alleles, condition_alleles)]; rownames(X) <- geno$IID
		exclude = rownames(X)[!rowSums(X)==2*(length(condition_genes)+1)]
		df_filt = df[!df$IID%in%exclude,]
		
		#Collapse rare alleles
		collapseAlleles = function(alleles, input.df) {
			var_count = NULL 
			for (i in 1:length(alleles)) {
				var_count[i] = sum(input.df[,alleles[i]])
			}
			common_alleles = alleles[var_count>=30]
			rare_alleles = alleles[var_count<30]
			
			if (length(rare_alleles)<=1) {
				collapsed = common_alleles		
			} else if (length(rare_alleles)>1) {
				
				groups = paste(unlist(lapply(strsplit(rare_alleles, split="_"),"[[", 1)), 
						unlist(lapply(strsplit(rare_alleles, split="_"),"[[", 2)),
						substring(unlist(lapply(strsplit(rare_alleles, split="_"),"[[", 3)), 1, 2),
						unlist(lapply(strsplit(rare_alleles, split="_"),"[[", 4)), sep="_")
				common_groups = NULL
				for (i in 1:length(unique(groups))) {
					group = unique(groups)[i]
					alleles_in_group = rare_alleles[groups==group]
					if (length(alleles_in_group)>1) {
						input.df[,group] = rowSums(input.df[,alleles_in_group])
						if (sum(input.df[,group]) >= 30) {
							common_groups = c(common_groups, group)	
						}
					}		
				}
				collapsed = c(common_alleles, common_groups) 
			}
			input.df[,collapsed]
		}

		predictor_alleles_collapsed = names(collapseAlleles(predictor_alleles, df_filt))
		condition_alleles_collapsed = names(collapseAlleles(condition_alleles, df_filt))
		all_alleles_collapsed = c(predictor_alleles_collapsed, condition_alleles_collapsed)	
		df_filt_new = data.frame(df_filt[,c("IID", "t1d", "PC1", "PC2")], collapseAlleles(predictor_alleles, df_filt), collapseAlleles(condition_alleles, df_filt))
		
		#Remove subjects with rare 2 digit alleles
		exclude2 = df_filt_new$IID[!rowSums(df_filt_new[,all_alleles_collapsed])==2*(length(condition_genes)+1)]
		common_df = df_filt_new[!df_filt_new$IID%in%exclude2,]
				
		#Get allele info
		case_count = NULL
		case_freq = NULL
		case_percent = NULL
		control_count = NULL
		control_freq = NULL
		control_percent = NULL
		for (i in 1:length(all_alleles_collapsed)) {
			case_count[i] = sum(common_df[common_df$t1d==1,all_alleles_collapsed[i]])
			control_count[i] = sum(common_df[common_df$t1d==0,all_alleles_collapsed[i]])
			case_freq[i] = sum(common_df[common_df$t1d==1,all_alleles_collapsed[i]])/(2*dim(common_df[common_df$t1d==1,])[1])
			case_percent[i] = case_freq[i]*100
			control_freq[i] = sum(common_df[common_df$t1d==0,all_alleles_collapsed[i]])/(2*dim(common_df[common_df$t1d==0,])[1])
			control_percent[i] = control_freq[i]*100
		}
		allele_info = data.frame(case_count, control_count, case_freq, control_freq); rownames(allele_info) <- all_alleles_collapsed
		#allele_info$Cases = with(allele_info, paste(case_count, " (", format(round(case_percent, 1),nsmall=1, trim=TRUE),")", sep=""))
		#allele_info$Controls = with(allele_info, paste(control_count, " (", format(round(control_percent, 1),nsmall=1, trim=TRUE),")", sep=""))
		allele_info$Cases = with(allele_info, format(round(case_percent, 1),nsmall=1, trim=TRUE))
		allele_info$Controls = with(allele_info, format(round(control_percent, 1),nsmall=1, trim=TRUE))
		
		
		##Get reference allele ### NOTE: MUST CHOOSE A REFERENCE ALLELE FOR EACH GENE INCLUDED IN THE MODEL!!!!
		allele_info$OR = (allele_info$case_freq/(1-allele_info$case_freq))/(allele_info$control_freq/(1-allele_info$control_freq))
		
		ref_alleles = NULL
		for (g in c(gene, condition_genes)) {
			sub.df = allele_info[grep(g, rownames(allele_info)),]
			min_logOR = min(abs(log(sub.df$OR)))
			ref_allele = rownames(sub.df)[abs(log(sub.df$OR))==min_logOR]
			ref_alleles = c(ref_alleles, ref_allele)
		}	
		#ref_allele = rownames(allele_info)[abs(log(OR))==min(abs(log(OR)), na.rm=TRUE) & !is.nan(OR)]
		
		
		
		#Fit model
		alleles_model = all_alleles_collapsed[!all_alleles_collapsed%in%ref_alleles]
		null_ind_vars = c(alleles_model[alleles_model%in%condition_alleles_collapsed], covariates)
		full_ind_vars = c(alleles_model, covariates)
		out.glm.null = glm(as.formula(paste("t1d~",paste(null_ind_vars,collapse="+"))), family="binomial", data=common_df)
		out.glm.full = glm(as.formula(paste("t1d~",paste(full_ind_vars,collapse="+"))), family="binomial", data=common_df)
		
		full.coeff = coef(summary(out.glm.full))[alleles_model,]
		full.ub = exp(full.coeff[,1]+1.96*full.coeff[,2])
		full.lb = exp(full.coeff[,1]-1.96*full.coeff[,2])
		full.OR = exp(full.coeff[,1])
		full.p = full.coeff[,4]
		full.df=data.frame(allele_info[alleles_model,], OR=format(round(full.OR, 2),nsmall=2), LB=format(round(full.lb,2),nsmall=2), UB=format(round(full.ub,2),nsmall=2), P=format(full.p, scientific=TRUE, digits=2))	
		full.df$OR_format = with(full.df, paste(OR, " [", LB,"-",UB,"]",sep=""))
		ref_allele.df = data.frame(allele_info[ref_allele,], OR=1.00, LB=NA, UB=NA, P=NA, OR_format="1.00 (reference)")
		out.df = rbind(full.df, ref_allele.df)
		
		rowname_strsplit = strsplit(rownames(out.df), split="_")
		out.df$Gene=unlist(lapply(rowname_strsplit,"[[", 2))
		out.df$Allele=unlist(lapply(rowname_strsplit,"[[", 3))
		
		OUT.TABLE = out.df[order(out.df$Allele, decreasing=FALSE),c("Allele", "Cases", "Controls", "OR_format", "P")]
		lrt_pvalue = format(lrtest(out.glm.full, out.glm.null)$"Pr(>Chisq)"[2], scientific=TRUE, digits=2)
		#addParagraph(rtffile, paste("Table 1-", k, " Regression coefficients for ", gene, sep=""))
		#addTable(rtffile, OUT.TABLE)
		addParagraph(rtffile, paste(gene,"\t", lrt_pvalue,"\n", sep=" "))
		print(paste(gene,"\t", lrt_pvalue,"\t sample size:",dim(common_df)[1],"\n", sep=" "))
	}
	addParagraph(rtffile, "\n")
}


getCondPvals(
		genes = c("HLA_A", "HLA_B", "HLA_C", "HLA_DRB1", "HLA_DQA1", "HLA_DPA1", "HLA_DPB1"), 
		condition_genes = c("HLA_DQB1")
) 


getCondPvals(
		genes = c("HLA_A", "HLA_B", "HLA_C", "HLA_DQA1", "HLA_DPA1", "HLA_DPB1"), 
		condition_genes = c("HLA_DQB1", "HLA_DRB1")
) 

getCondPvals(
		genes = c("HLA_A", "HLA_B", "HLA_C", "HLA_DPA1", "HLA_DPB1"), 
		condition_genes = c("HLA_DQB1", "HLA_DRB1", "HLA_DQA1")
) 

getCondPvals(
		genes = c("HLA_B", "HLA_C", "HLA_DPA1", "HLA_DPB1"), 
		condition_genes = c("HLA_DQB1", "HLA_DRB1", "HLA_DQA1", "HLA_A")
) 

getCondPvals(
		genes = c("HLA_C", "HLA_DPA1", "HLA_DPB1"), 
		condition_genes = c("HLA_DQB1", "HLA_DRB1", "HLA_DQA1", "HLA_A", "HLA_B")
) 


p_bonf = 0.05/8

done(rtffile)


#HLA-A (p=3.5x10-8), HLA-B (p=2.1x10-12) HLA-DRB1 (p=8.0x10-23) and HLA-DQA1 (p=1.7x10-12) remained significant