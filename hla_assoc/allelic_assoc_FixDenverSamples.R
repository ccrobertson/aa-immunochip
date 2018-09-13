options(stringsAsFactors=FALSE)
setwd("/data1/ccr5ju/AA_Immuno/FixDenverSamples")
library(lmtest)
library(rtf)
library(ggplot2)

genofile = "/data1/ccr5ju/AA_Immuno/hla_assoc/processed_impdata/genotypes_recode.raw.4d"
rtffilename = "allelic_assoc_tables.doc"
phenofile ="/data1/ccr5ju/AA_Immuno/raw/pc2.txt"
hufile = "/data1/ccr5ju/AA_Immuno/hla_assoc/sources/ng.2015_Hu.supptable2.reformatted.csv"


#Get data
geno = read.table(genofile, comment.char="~", header=TRUE)
pheno = read.table(phenofile, comment.char="~", header=TRUE)
covariates = c("PC1", "PC2","SEX")
genes = c("HLA_A", "HLA_B", "HLA_C", "HLA_DRB1", "HLA_DQA1", "HLA_DQB1", "HLA_DPA1", "HLA_DPB1")
data = merge(geno, pheno, by="IID")

#Fix Denver samples --> set as controls
denver_fixlist = read.table("remove_t2d_bdc.txt", header=FALSE)
names(denver_fixlist) = c("FID", "IID", "reason")
data$PHENOTYPE[data$IID%in%denver_fixlist$IID] <- 1
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
#cd /Users/Cassie/Box\ Sync/Rich\ Lab/AA_Immunochip/FixDenverSamples
#scp ccr5ju@dobby.cphg.virginia.edu:/data1/ccr5ju/AA_Immuno/FixDenverSamples/supp_tables.doc .


