options(stringsAsFactors=FALSE)
setwd("/m/jdrfdn_scratch/users/ccr5ju/AA_Immuno/hla_assoc")
#install.packages("lmtest", repos="https://cloud.r-project.org", lib="/home/ccr5ju/R/x86_64-redhat-linux-gnu-library/3.3")
library(lmtest, lib.loc="/home/ccr5ju/R/x86_64-redhat-linux-gnu-library/3.3")
library(rtf)
library(ggplot2)

genofile = "derived_data/genotypes_recode.raw.4d"
rtffilename = "tables/table1.doc"
phenofile ="../data_original/pc2.txt"
hufile = "ng.2015_Hu.supptable2.reformatted.csv"



#Get data
geno = read.table(genofile, comment.char="~", header=TRUE)
pheno = read.table(phenofile, comment.char="~", header=TRUE)
covariates = c("PC1", "PC2","SEX")
genes = c("HLA_A", "HLA_B", "HLA_C", "HLA_DRB1", "HLA_DQA1", "HLA_DQB1", "HLA_DPA1", "HLA_DPB1")
data = merge(geno, pheno, by="IID")
data$t1d = data$PHENOTYPE-1

#Quality filtering (exclude subjects whose alleles do not add to 2 for a given locus)
filterData = function(gene, df) {
	alleles=names(geno)[grep(gene, names(geno))]
	X = geno[,alleles]; rownames(X) <- geno$IID
	rownames(X)[!rowSums(X)==2]
}
exclude = unique(unlist(sapply(genes, filterData, geno)))
data_filt = data[!data$IID%in%exclude,]

#Look at excluded subjects
chisq.test(table(data$IID%in%exclude, data$t1d))
chisq.test(table(data$IID%in%exclude, data$SEX))

#Collapse or remove rare alleles
collapseAlleles = function(gene, df) {
	alleles=names(df)[grep(gene, names(df))]
	var_count = NULL 
	for (i in 1:length(alleles)) {
		var_count[i] = sum(df[,alleles[i]])
	}
	common_alleles = alleles[var_count>=30]
	rare_alleles = alleles[var_count<30]
	
	if (length(rare_alleles)<=1) {
		alleles_new = common_alleles		
	} else if (length(rare_alleles)>1) {
		
		groups = paste(unlist(lapply(strsplit(rare_alleles, split="_"),"[[", 1)), 
				unlist(lapply(strsplit(rare_alleles, split="_"),"[[", 2)),
				substring(unlist(lapply(strsplit(rare_alleles, split="_"),"[[", 3)), 1, 2),
				unlist(lapply(strsplit(rare_alleles, split="_"),"[[", 4)), sep="_")
		rare_allele_dict = list()
		common_groups = NULL
		for (i in 1:length(unique(groups))) {
			group = unique(groups)[i]
			alleles_in_group = rare_alleles[groups==group]
			rare_allele_dict[[i]] = list(group=group, alleles=alleles_in_group, count=sum(df[,alleles_in_group]))
			if (length(alleles_in_group)>1) {
				df[,group] = rowSums(df[,alleles_in_group])
				if (sum(df[,group]) >= 30) {
					common_groups = c(common_groups, group)	
				}
			}		
		}
		alleles_new = c(common_alleles, common_groups) 
	}
	out = df[,alleles_new]; row.names(out) <- df$IID
	out
}
collapsed = do.call(cbind,lapply(genes, collapseAlleles, data_filt))

#Generate data frames for association analysis with alleles
all_df = data.frame(data_filt[,c("IID","t1d",covariates, names(data_filt)[grep("HLA", names(data_filt))])])
common_df = data.frame(data_filt[,c("IID", "t1d", covariates)], collapsed) 
#exclude subjects with rare alleles --> 
common_df_2 = common_df[rowSums(common_df[,c(-1,-2,-3,-4,-5)])==16,]

#Get case and control frequencies of prepped data
hu = read.csv(hufile)
hu_hla = hu[grep("HLA_",hu$SNP2HLA.variant.name),]
rownames(hu_hla) = paste(hu_hla$SNP2HLA.variant.name, "P", sep="_")
getAlleleInfo = function(df){
	#Get allele info
	alleles = names(df)[grep("HLA", names(df))]
	count = NULL
	case_count = NULL
	case_freq = NULL
	case_percent = NULL
	control_count = NULL
	control_freq = NULL
	control_percent = NULL
	for (i in 1:length(alleles)) {
		count[i] = sum(df[,alleles[i]])
		case_count[i] = sum(df[df$t1d==1,alleles[i]])
		control_count[i] = sum(df[df$t1d==0,alleles[i]])
		case_freq[i] = sum(df[df$t1d==1,alleles[i]])/(2*dim(df[df$t1d==1,])[1])
		control_freq[i] = sum(df[df$t1d==0,alleles[i]])/(2*dim(df[df$t1d==0,])[1])
	}
	allele_info = data.frame(count, case_count, control_count, case_freq, control_freq); rownames(allele_info) <- alleles
	allele_info$case_AF = with(allele_info, format(round(case_freq, 3),nsmall=1, trim=TRUE))
	allele_info$control_AF = with(allele_info, format(round(control_freq, 3),nsmall=1, trim=TRUE))
	allele_info$OR_unadjusted = round((allele_info$case_freq/(1-allele_info$case_freq))/(allele_info$control_freq/(1-allele_info$control_freq)), digits=2)
	
	#add Hu et al info
	allele_info$Hu.case_AF = hu_hla[rownames(allele_info),"Case.frequency"]
	allele_info$Hu.control_AF = hu_hla[rownames(allele_info),"Control.frequency"]
	allele_info$Hu.OR =  hu_hla[rownames(allele_info),"OR.Total."]
	allele_info
}
allele_info_all = getAlleleInfo(all_df) #calculate freqs in all alleles
allele_info_common = getAlleleInfo(common_df) #calculate freqs in common/aggregated alleles only (note, do NOT want to use common_df_2 because this will bias AF estimates) 


####Fit multiallelic model
#multiallelic analysis <-- note, DO want to use common_df_2 because otherwise, rare alleles will be interpreted as equivalence to reference allele in model 
fitMultiModel = function(gene, df, allele_info, covariates) {
	
	#get alleles
	alleles = names(df)[grep(gene, names(df))]
	
	#Get ref_allele
	#neutral ref allele
	info = allele_info[alleles,]
	ref_allele = rownames(info)[abs(log(info$OR))==min(abs(log(info$OR)), na.rm=TRUE) & !is.nan(info$OR)]
	#high-risk reference allele
	#ref_allele2 = rownames(allele_info)[log(OR)==max(log(OR), na.rm=TRUE) & !is.nan(OR)]
	
	#Fit model
	alleles_model = alleles[!alleles==ref_allele]
	out.glm.null = glm(as.formula(paste("t1d~",paste(covariates,collapse="+"))), family="binomial", data=df)
	out.glm.full = glm(as.formula(paste("t1d~",paste(c(alleles_model, covariates),collapse="+"))), family="binomial", data=df)
	
	#get model info
	full.coeff = coef(summary(out.glm.full))[alleles_model,]
	full.OR = exp(full.coeff[,1])
	full.ub = exp(full.coeff[,1]+1.96*full.coeff[,2])
	full.lb = exp(full.coeff[,1]-1.96*full.coeff[,2])
	full.p = full.coeff[,4]
		
	#format model info output
	full.df=data.frame(allele_info[alleles_model,], OR=format(round(full.OR, 2),nsmall=2), LB=format(round(full.lb,2),nsmall=2), UB=format(round(full.ub,2),nsmall=2), P=format(full.p, scientific=TRUE, digits=2))
	full.df$OR_format = with(full.df, paste(OR, " [", LB,"-",UB,"]",sep=""))
	ref_allele.df = data.frame(allele_info[ref_allele,], OR=1.00, LB=NA, UB=NA, P=NA, OR_format="1.00 (reference)")
	out.df = rbind(full.df, ref_allele.df)

	#### Testing out different statistical significance measures
#	#get full-1 p-values
#	getM1Pval = function(l) {
#		alleles_m1 = alleles_model[!alleles_model==l]
#		out.glm.fullm1 = glm(as.formula(paste("t1d~",paste(c(alleles_m1, covariates),collapse="+"))), family="binomial", data=common_df)
#		format(lrtest(out.glm.full, out.glm.fullm1)$"Pr(>Chisq)"[2], scientific=TRUE, digits=2)
#	}
#	fullm1.p = sapply(alleles_model, getM1Pval)
#	
#	#mu + beta1 - mu0
#	teststat = summary(out.glm.full)$coef["(Intercept)",1] + summary(out.glm.full)$coef[alleles_model,1] - summary(out.glm.null)$coef["(Intercept)",1]
#	
	#full.df=data.frame(allele_info[alleles_model,], OR=format(round(full.OR, 2),nsmall=2), LB=format(round(full.lb,2),nsmall=2), UB=format(round(full.ub,2),nsmall=2), P=format(full.p, scientific=TRUE, digits=2), P.m1 = fullm1.p, beta_star=teststat)	
	#full.df$OR_format = with(full.df, paste(OR, " [", LB,"-",UB,"]",sep=""))
	#ref_allele.df = data.frame(allele_info[ref_allele,], OR=1.00, LB=NA, UB=NA, P=NA, P.m1=NA, beta_star=NA, OR_format="1.00 (reference)", OR_bm)
	#out.df = rbind(full.df, ref_allele.df)
		
	rowname_strsplit = strsplit(rownames(out.df), split="_")
	out.df$Gene=unlist(lapply(rowname_strsplit,"[[", 2))
	out.df$Allele=unlist(lapply(rowname_strsplit,"[[", 3))
	
	OUT.TABLE = out.df[order(out.df$Allele, decreasing=FALSE),c("Allele", "case_AF", "control_AF", "OR", "OR_format", "P")]
	lrt_pvalue = format(lrtest(out.glm.full, out.glm.null)$"Pr(>Chisq)"[2], scientific=TRUE, digits=2)
	list(OUT.TABLE, lrt_pvalue)
}
OUT.MULT = lapply(genes, fitMultiModel, common_df_2, allele_info_common, covariates)

##multiallelic conditional analysis  ## UPDATE THIS !!!
getCondPvals = function(condition_genes) {
#	addParagraph(rtffile, paste("Conditioning on:", paste(condition_genes,collapes="\t"),"\n", sep=" "))
#	for (k in 1:length(genes)) {
#		gene = genes[k]
#		
#		condition_alleles = NULL
#		for(f in 1:length(condition_genes)) {
#			condition_alleles = c(condition_alleles, names(geno)[grep(condition_genes[f], names(geno))])
#		}	
#		predictor_alleles = names(geno)[grep(gene, names(geno))]
#		
#		#Quality filtering (exclude subjects whose alleles do not add to 2 for a given locus)
#		X = geno[,c(predictor_alleles, condition_alleles)]; rownames(X) <- geno$IID
#		exclude = rownames(X)[!rowSums(X)==2*(length(condition_genes)+1)]
#		df_filt = df[!df$IID%in%exclude,]
#		
#		#Collapse rare alleles
#		collapseAlleles = function(alleles, input.df) {
#			var_count = NULL 
#			for (i in 1:length(alleles)) {
#				var_count[i] = sum(input.df[,alleles[i]])
#			}
#			common_alleles = alleles[var_count>=30]
#			rare_alleles = alleles[var_count<30]
#			
#			if (length(rare_alleles)<=1) {
#				collapsed = common_alleles		
#			} else if (length(rare_alleles)>1) {
#				
#				groups = paste(unlist(lapply(strsplit(rare_alleles, split="_"),"[[", 1)), 
#						unlist(lapply(strsplit(rare_alleles, split="_"),"[[", 2)),
#						substring(unlist(lapply(strsplit(rare_alleles, split="_"),"[[", 3)), 1, 2),
#						unlist(lapply(strsplit(rare_alleles, split="_"),"[[", 4)), sep="_")
#				common_groups = NULL
#				for (i in 1:length(unique(groups))) {
#					group = unique(groups)[i]
#					alleles_in_group = rare_alleles[groups==group]
#					if (length(alleles_in_group)>1) {
#						input.df[,group] = rowSums(input.df[,alleles_in_group])
#						if (sum(input.df[,group]) >= 30) {
#							common_groups = c(common_groups, group)	
#						}
#					}		
#				}
#				collapsed = c(common_alleles, common_groups) 
#			}
#			input.df[,collapsed]
#		}
#		
#		predictor_alleles_collapsed = names(collapseAlleles(predictor_alleles, df_filt))
#		condition_alleles_collapsed = names(collapseAlleles(condition_alleles, df_filt))
#		all_alleles_collapsed = c(predictor_alleles_collapsed, condition_alleles_collapsed)	
#		df_filt_new = data.frame(df_filt[,c("IID", "t1d", "PC1", "PC2")], collapseAlleles(predictor_alleles, df_filt), collapseAlleles(condition_alleles, df_filt))
#		
#		#Remove subjects with rare 2 digit alleles
#		exclude2 = df_filt_new$IID[!rowSums(df_filt_new[,all_alleles_collapsed])==2*(length(condition_genes)+1)]
#		common_df = df_filt_new[!df_filt_new$IID%in%exclude2,]
#		
#		#Get allele info
#		case_count = NULL
#		case_freq = NULL
#		case_percent = NULL
#		control_count = NULL
#		control_freq = NULL
#		control_percent = NULL
#		for (i in 1:length(all_alleles_collapsed)) {
#			case_count[i] = sum(common_df[common_df$t1d==1,all_alleles_collapsed[i]])
#			control_count[i] = sum(common_df[common_df$t1d==0,all_alleles_collapsed[i]])
#			case_freq[i] = sum(common_df[common_df$t1d==1,all_alleles_collapsed[i]])/(2*dim(common_df[common_df$t1d==1,])[1])
#			case_percent[i] = case_freq[i]*100
#			control_freq[i] = sum(common_df[common_df$t1d==0,all_alleles_collapsed[i]])/(2*dim(common_df[common_df$t1d==0,])[1])
#			control_percent[i] = control_freq[i]*100
#		}
#		allele_info = data.frame(case_count, control_count, case_freq, control_freq); rownames(allele_info) <- all_alleles_collapsed
#		#allele_info$Cases = with(allele_info, paste(case_count, " (", format(round(case_percent, 1),nsmall=1, trim=TRUE),")", sep=""))
#		#allele_info$Controls = with(allele_info, paste(control_count, " (", format(round(control_percent, 1),nsmall=1, trim=TRUE),")", sep=""))
#		allele_info$Cases = with(allele_info, format(round(case_percent, 1),nsmall=1, trim=TRUE))
#		allele_info$Controls = with(allele_info, format(round(control_percent, 1),nsmall=1, trim=TRUE))
#		
#		
#		##Get reference allele ### NOTE: MUST CHOOSE A REFERENCE ALLELE FOR EACH GENE INCLUDED IN THE MODEL!!!!
#		allele_info$OR = (allele_info$case_freq/(1-allele_info$case_freq))/(allele_info$control_freq/(1-allele_info$control_freq))
#		
#		ref_alleles = NULL
#		for (g in c(gene, condition_genes)) {
#			sub.df = allele_info[grep(g, rownames(allele_info)),]
#			min_logOR = min(abs(log(sub.df$OR)))
#			ref_allele = rownames(sub.df)[abs(log(sub.df$OR))==min_logOR]
#			ref_alleles = c(ref_alleles, ref_allele)
#		}	
#		#ref_allele = rownames(allele_info)[abs(log(OR))==min(abs(log(OR)), na.rm=TRUE) & !is.nan(OR)]
#		
#		
#		
#		#Fit model
#		alleles_model = all_alleles_collapsed[!all_alleles_collapsed%in%ref_alleles]
#		null_ind_vars = c(alleles_model[alleles_model%in%condition_alleles_collapsed], covariates)
#		full_ind_vars = c(alleles_model, covariates)
#		out.glm.null = glm(as.formula(paste("t1d~",paste(null_ind_vars,collapse="+"))), family="binomial", data=common_df)
#		out.glm.full = glm(as.formula(paste("t1d~",paste(full_ind_vars,collapse="+"))), family="binomial", data=common_df)
#		
#		full.coeff = coef(summary(out.glm.full))[alleles_model,]
#		full.ub = exp(full.coeff[,1]+1.96*full.coeff[,2])
#		full.lb = exp(full.coeff[,1]-1.96*full.coeff[,2])
#		full.OR = exp(full.coeff[,1])
#		full.p = full.coeff[,4]
#		full.df=data.frame(allele_info[alleles_model,], OR=format(round(full.OR, 2),nsmall=2), LB=format(round(full.lb,2),nsmall=2), UB=format(round(full.ub,2),nsmall=2), P=format(full.p, scientific=TRUE, digits=2))	
#		full.df$OR_format = with(full.df, paste(OR, " [", LB,"-",UB,"]",sep=""))
#		ref_allele.df = data.frame(allele_info[ref_allele,], OR=1.00, LB=NA, UB=NA, P=NA, OR_format="1.00 (reference)")
#		out.df = rbind(full.df, ref_allele.df)
#		
#		rowname_strsplit = strsplit(rownames(out.df), split="_")
#		out.df$Gene=unlist(lapply(rowname_strsplit,"[[", 2))
#		out.df$Allele=unlist(lapply(rowname_strsplit,"[[", 3))
#		
#		OUT.TABLE = out.df[order(out.df$Allele, decreasing=FALSE),c("Allele", "Cases", "Controls", "OR_format", "P")]
#		lrt_pvalue = format(lrtest(out.glm.full, out.glm.null)$"Pr(>Chisq)"[2], scientific=TRUE, digits=2)
#		#addParagraph(rtffile, paste("Table 1-", k, " Regression coefficients for ", gene, sep=""))
#		#addTable(rtffile, OUT.TABLE)
#		addParagraph(rtffile, paste(gene,"\t", lrt_pvalue,"\n", sep=" "))
#		print(paste(gene,"\t", lrt_pvalue,"\t sample size:",dim(common_df)[1],"\n", sep=" "))
#	}
#	addParagraph(rtffile, "\n")
}
#getCondPvals(condition_genes = c("HLA_DQB1")) 
#getCondPvals(condition_genes = c("HLA_DQB1", "HLA_DRB1")) 
#getCondPvals(condition_genes = c("HLA_DQB1", "HLA_DRB1", "HLA_DQA1")) 
#getCondPvals(condition_genes = c("HLA_DQB1", "HLA_DRB1", "HLA_DQA1", "HLA_A")) 
#getCondPvals(condition_genes = c("HLA_DQB1", "HLA_DRB1", "HLA_DQA1", "HLA_A", "HLA_B")) 


getVarianceExplained = function(gene) {
	#Get genotype freq and OR
	k = which(genes==gene)
	model_info = OUT.MULT[[k]][[1]]
	alleles = names(common_df)[grep(gene, names(common_df))]
	genotypes = NULL
	control_count = NULL
	frq_obs = NULL
	frq_est = NULL
	OR_est = NULL
	for (i in 1:length(alleles)) {
		for (j in 1:length(alleles)) {
			allele12 = paste(sort(c(alleles[i], alleles[j])),collapse="/")
			if (!allele12 %in% genotypes) {
				genotypes = c(genotypes, allele12)
				if (i==j) {
					count = sum(common_df[common_df$t1d==0,alleles[i]]==2)
					frq = allele_info_common[alleles[i],]$control_freq*allele_info_common[alleles[j],]$control_freq
					OR = as.numeric(model_info[alleles[i],'OR'])^2
				} else {
					count = sum(common_df[common_df$t1d==0,alleles[i]]==1 & common_df[common_df$t1d==0,alleles[j]]==1)
					frq = 2*allele_info_common[alleles[i],]$control_freq*allele_info_common[alleles[j],]$control_freq
					OR = as.numeric(model_info[alleles[i],'OR'])*as.numeric(model_info[alleles[j],'OR'])
				}
				control_count = c(control_count, count) 
				frq_obs = c(frq_obs, count/sum(common_df$t1d==0))
				frq_est =  c(frq_est, frq)
				OR_est = c(OR_est, OR)
			}
		}
	}
	g_df = data.frame(genotypes, control_count, frq_obs, frq_est, OR_est)
	g_df$est_count = g_df$frq_est*sum(common_df$t1d==0)
	g_df$OR_1 = sqrt(g_df$OR_est)
	#plot(g_df$est_count, g_df$control_count) + abline(coef=c(0,1))
	
	K = 0.002  #population prevalence
	g_df$prev_g = (g_df$frq_est*K)/sum(g_df$frq_est*g_df$OR_est)
	ref_allele = rownames(model_info)[model_info$OR_format=="1.00 (reference)"]
	g_df$T = qnorm(g_df$prev_g, 0, 1, lower.tail=FALSE)
	T_ref = g_df[g_df$genotypes==paste(ref_allele,ref_allele, sep="/"),"T"] 
	g_df$shift = g_df$T-T_ref
	x_bar = sum(g_df$frq_est*g_df$shift)
	var_g = sum(g_df$frq_est*(g_df$shift-x_bar)^2)
	list(gene, var_g/(var_g+1), g_df)
}
#lapply(genes, getVarianceExplained)
getVarianceExplained("HLA_DQB1")


#create file for INDI-V
	model_info = OUT.MULT[[which(genes=="HLA_DQB1")]][[1]]
	model_info$OR_2 = as.numeric(model_info$OR)^2 
	model_info$K = 0.002
	model_info$Lambdas = 15
	write.table(model_info[,c("control_AF", "OR", "OR_2", "K","Lambdas")], file="input_for_INDI-V.csv",row.names=TRUE, col.names=TRUE, quote=FALSE, sep=",")

#alleles = names(common_df)[grep("HLA_", names(common_df))]
#genotypes = matrix(NA, nrow=dim(common_df)[1], ncol=9)
#colnames(genotypes) = c("IID", "A", "B", "C", "DPA1","DPB1","DQA1","DQB1","DRB1")
#for (i in 1:dim(common_df)[1]) {
#	g = sort(c(alleles[common_df[i,alleles]==1], rep(alleles[common_df[i,alleles]==2],2)))
#	genos = rep(NA, 8)
#	for (m in 1:8) {
#		genos[m] = paste(g[2*m-1],g[2*m], sep="/")
#	}
#	genotypes[i,] = c(common_df$IID[i], genos)
#}

#Get variance explained
getVarExplainedPerAllele = function(allele) {
	K = 0.002
	t = table(common_df[,allele], common_df$t1d)
	p = (t[2,1]+2*t[3,1])/(2*sum(common_df$t1d==0))
	q = 1-p
	OR_Aa = (t[1,1]*t[2,2])/(t[2,1]*t[1,2])  
	OR_AA = (t[1,1]*t[3,2])/(t[3,1]*t[1,2])
	k_aa = K/(q^2 + 2*p*q*OR_Aa + q^2*OR_AA) 
	P_0 = k_aa
	P_1 = k_aa*OR_Aa
	P_2 = k_aa*OR_AA
	T_0 = qnorm(P_0, 0, 1, lower.tail=FALSE)
	T_1 = qnorm(P_1, 0, 1, lower.tail=FALSE)
	T_2 = qnorm(P_2, 0, 1, lower.tail=FALSE)
	tt = list(c(p^2, T_2-T_0), c(2*p*q, T_1-T_0), c(q^2, 0))
	x_bar = tt[[1]][1]*tt[[1]][2] + tt[[2]][1]*tt[[2]][2] + tt[[3]][1]*tt[[3]][2]
	var_g = tt[[1]][1]*(tt[[1]][2]-x_bar)^2 + tt[[2]][1]*(tt[[2]][2]-x_bar)^2 + tt[[2]][1]*(tt[[2]][2]-x_bar)^2
	var_g/(var_g+1)
}
getVarExplainedPerAllele("HLA_A_0101_P")
getVarExplainedPerAllele("HLA_DQB1_0201_P")
getVarExplainedPerAllele("HLA_DRB1_0301_P")




####Fit biallelic models
#fit biallelic models (common alleles only)
fitBiModel = function(df, allele_info, covariates) {
	alleles = names(df)[grep("HLA", names(df))]
	out.glm.null = glm(as.formula(paste("t1d~",paste(covariates,collapse="+"))), family="binomial", data=df)		
	#get biallelic p-value
	getBAPval = function(l) {
		if (allele_info[l,"count"] >= 30) {
			out.glm.l = glm(as.formula(paste("t1d~",paste(c(l, covariates),collapse="+"))), family="binomial", data=df)
			beta = summary(out.glm.l)$coeff[l,"Estimate"]
			or = format(exp(beta), digits=2, nsmall=2)
			lrt_p = format(lrtest(out.glm.l, out.glm.null)$"Pr(>Chisq)"[2], scientific=TRUE)	
		} else {
			or = NA
			lrt_p = NA
		}
		c(or, lrt_p)
	}
	l.p = as.data.frame(t(sapply(alleles, getBAPval))); names(l.p) <- c("OR","P")
	#format model info output
	out.df=data.frame(allele_info[alleles,], l.p[alleles,])	
	rowname_strsplit = strsplit(rownames(out.df), split="_")
	out.df$Gene=unlist(lapply(rowname_strsplit,"[[", 2))
	out.df$Allele=unlist(lapply(rowname_strsplit,"[[", 3))
	OUT.BI = out.df[order(out.df$Gene, out.df$Allele, decreasing=FALSE),c("Gene", "Allele", "count", "control_AF", "case_AF", "OR_unadjusted","OR", "P", "Hu.control_AF", "Hu.case_AF", "Hu.OR")]	
	OUT.BI
}
OUT.BI.COMMON = fitBiModel(common_df, allele_info_common, covariates)  #note, here we use common_df instead of common_df_2 because there is no reference allele, simply testing allele vs other alleles (or allele vs not allele)
OUT.BI.ALL = fitBiModel(all_df, allele_info_all, covariates)

p_thresh = 0.05/dim(OUT.BI.COMMON)[1]


####Biallelic conditional analysis
fitBiModel = function(df, allele_info, covariates) {
	alleles_1 = names(df)[grep("HLA", names(df))]
	alleles = alleles_1[!alleles_1%in%covariates]
	out.glm.null = glm(as.formula(paste("t1d~",paste(covariates,collapse="+"))), family="binomial", data=df)		
	#get biallelic p-value
	getBAPval = function(l) {
		if (allele_info[l,"count"] >= 30) {
			out.glm.l = glm(as.formula(paste("t1d~",paste(c(l, covariates),collapse="+"))), family="binomial", data=df)
			beta = summary(out.glm.l)$coeff[l,"Estimate"]
			or = format(exp(beta), digits=2, nsmall=2)
			lrt_p = lrtest(out.glm.l, out.glm.null)$"Pr(>Chisq)"[2]	
		} else {
			or = NA
			lrt_p = NA
		}
		c(or, lrt_p)
	}
	l.p = as.data.frame(t(sapply(alleles, getBAPval))); names(l.p) <- c("OR","P")
	#format model info output
	out.df=data.frame(allele_info[alleles,], l.p[alleles,])	
	rowname_strsplit = strsplit(rownames(out.df), split="_")
	out.df$Gene=unlist(lapply(rowname_strsplit,"[[", 2))
	out.df$Allele=unlist(lapply(rowname_strsplit,"[[", 3))
	OUT.BI = out.df[order(out.df$Gene, out.df$Allele, decreasing=FALSE),c("Gene", "Allele", "count", "control_AF", "case_AF", "OR_unadjusted","OR", "P", "Hu.control_AF", "Hu.case_AF", "Hu.OR")]	
	OUT.BI
}
#initialize
p_thresh = 0.05/dim(OUT.BI.COMMON)[1]
p_min = 0
round = 0
condition_alleles = NULL
covars =  c("PC1", "PC2","SEX")
top_allele = NA
#run rounds
while(p_min<p_thresh) {
	cat(round, p_min, top_allele,"\n")
	out = fitBiModel(common_df, allele_info_common, covars)	
	p_min = min(as.numeric(out$P))
	top_allele = row.names(out)[out$P==p_min]	
	
	#update
	round = round+1
	condition_alleles = c(condition_alleles, top_allele)
	covars = c(covars, top_allele)
}
out_all = glm(as.formula(paste("t1d~",paste(covars[1:(length(covars)-1)],collapse="+"))), family="binomial", data=common_df)
OUT.COND = as.data.frame(summary(out_all)$coef[condition_alleles[1:(length(condition_alleles)-1)],c("Estimate", "Pr(>|z|)")])
OUT.COND$OR = round(exp(OUT.COND$Estimate),digits=2) 
OUT.COND$P = format(OUT.COND[,"Pr(>|z|)"], scientific=TRUE, digits=2)
rowname_strsplit = strsplit(rownames(OUT.COND), split="_")
OUT.COND$Gene=unlist(lapply(rowname_strsplit,"[[", 2))
OUT.COND$Allele=unlist(lapply(rowname_strsplit,"[[", 3))



#Write tables to file
rtffile <- RTF(rtffilename)  # this can be an .rtf or a .doc

addParagraph(rtffile, "Supplementary Table 5. Associations of T1D in African Americans with classical HLA alleles")
addTable(rtffile, OUT.BI.COMMON)
addParagraph(rtffile,"Odds ratio estimates and p-values generated by logistic regression, adjusting for 2 principal components and sex. Rare 4-digit alleles are combined into 2-digit groups. Only alleles (4-digit or 2-digit aggregates of rare 4-digit) with total allele count >30 are analyzed. \n\n")

#addParagraph(rtffile, "Supplementary Table ##. Association of T1D in African Americans with classical HLA alleles")
#addTable(rtffile, OUT.BI.ALL)
#addParagraph(rtffile,"Odds ratio estimates and p-values generated by logistic regression, adjusting for 2 principal components and sex. Only alleles with total allele count >30 are analyzed.\n\n")

addParagraph(rtffile, "Supplementary Table 6. Nine classical HLA alleles are independently associated with T1D in African Americans")
addTable(rtffile, OUT.COND[,c("Gene", "Allele", "OR", "P")])
addParagraph(rtffile,"Odds ratio estimates and p-values generated by logistic regression, adjusting for 2 principal components and sex. Only alleles with total allele count >30 are analyzed.\n\n")



#for (k in 1:length(genes)) {
#	gene = genes[k]
#	addParagraph(rtffile, paste("Supplementary Table 2-", k, ". Regression coefficients for ", gene, sep=""))
#	addTable(rtffile, OUT.MULT[[k]][[1]])
#	addParagraph(rtffile, paste("P-value for likelihood ratio test (full vs null) of model:", OUT.MULT[[k]][[2]],"\n\n", sep=" "))
#}

done(rtffile)
#scp tables/table1.doc "$local/AA_Immunochip/tables_and_figures/"








##Look at associations in Eur vs Afr
OUT.BI.COMMON$MHC_class = NA 
OUT.BI.COMMON$MHC_class[OUT.BI.COMMON$Gene %in% c("DRB1", "DQA1","DQB1", "DPA1", "DPB1")] <- "class II"
OUT.BI.COMMON$MHC_class[OUT.BI.COMMON$Gene %in% c("A", "B","C")] <- "class I"
OUT.BI.COMMON$diff = log(as.numeric(OUT.BI.COMMON$OR)) - log(as.numeric(OUT.BI.COMMON$Hu.OR)) 
ggplot(OUT.BI.COMMON, aes(diff, as.numeric(count))) + geom_point()
ggplot(OUT.BI.COMMON, aes(diff, as.numeric(control_AF))) + geom_point() + xlab("Difference in log(OR) AA")
OUT.BI.COMMON[order(OUT.BI.COMMON$diff),]


pdf("tables/afr_vs_eur.pdf")
ggplot(OUT.BI.COMMON, aes(log(as.numeric(OR)), log(-log(as.numeric(P))), colour=MHC_class)) + geom_point() + xlab("Log(OR) in African Americans") + ylab("-log10(P) in African Americans")
ggplot(OUT.BI.COMMON, aes(log(as.numeric(Hu.OR)), log(-log(as.numeric(P))), colour=MHC_class)) + geom_point() + xlab("Log(OR) in Europeans") + ylab("-log10(P) in African Americans")



ggplot(OUT.BI.COMMON, aes(log(as.numeric(Hu.OR)), log(as.numeric(OR)), colour=Gene)) + geom_point() + xlab("Log(OR) in Europeans") + ylab("Log(OR) in African Americans") + geom_abline(intercept=0, slope=1, color="red", linetype="dashed", size=1.5)
ggplot(OUT.BI.COMMON, aes(log(as.numeric(Hu.OR)), log(as.numeric(OR)), colour=MHC_class)) + geom_point()+ xlab("Log(OR) in Europeans") + ylab("Log(OR) in African Americans")
dev.off()

ggplot(OUT.BI.COMMON, aes(log(as.numeric(Hu.OR)), log(as.numeric(OR)), colour=as.numeric(control_AF))) + geom_point()
ggplot(OUT.BI.COMMON, aes(log(as.numeric(Hu.OR)), log(as.numeric(OR)), colour=as.numeric(Hu.control_AF))) + geom_point()
ggplot(OUT.BI.ALL, aes(log(as.numeric(Hu.case_AF)), log(as.numeric(case_AF)), colour=as.numeric(Hu.OR))) + geom_point()
