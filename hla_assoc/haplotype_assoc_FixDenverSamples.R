options(stringsAsFactors=FALSE)
setwd("/data1/ccr5ju/AA_Immuno/FixDenverSamples")
library(rtf)
library(lmtest)

geno_1 = read.table("/data1/ccr5ju/AA_Immuno/hla_assoc/processed_impdata/haplotypes.4d.txt", header=TRUE, sep="\t", comment.char="~") #file with 2 rows per subject
geno = geno_1[complete.cases(geno_1),] #remove missing data
geno$IID = unlist(lapply(strsplit(geno$uniq_id, split="_"), "[[", 2))
t1d_map = unique(geno[,c("IID", "case")]); rownames(t1d_map) <- t1d_map$IID
pcs=read.table("/data1/ccr5ju/AA_Immuno/raw/pc2.txt", header=TRUE, comment.char="~"); rownames(pcs) <- pcs$IID
geno_2 = read.table("/data1/ccr5ju/AA_Immuno/hla_assoc/processed_impdata/genotypes_recode.raw.4d", comment.char="~", header=TRUE)
sex_map = geno_2[,c("IID", "SEX")]; rownames(sex_map) <- sex_map$IID
covariates = c("PC1", "PC2","SEX")

erlichfile = "/data1/ccr5ju/AA_Immuno/hla_assoc/sources/diab.2008_Erlich.table2_v2.txt"
erlich = read.table(erlichfile, sep="\t", header=TRUE)
write.table(erlich, file="diab.2008_Erlich.table2_v3.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
rownames(erlich) <- paste("H", sprintf( '%04d', erlich$DRB1), sprintf( '%04d',erlich$DQA1), sprintf( '%04d',erlich$DQB1), sep="_")


##Create haplotype matrix (Note this takes time -- best to load from file)
#haplotypes = unique(geno$haplotype)
#subjects = unique(geno$IID)
#hapMat = matrix(NA, nrow=length(subjects), ncol=length(haplotypes))
#for (i in 1:length(subjects)) {
#	sub = subjects[i]
#	for (j in 1:length(haplotypes)) {
#		hap = haplotypes[j]
#		hapMat[i, j] = sum(geno$IID==sub & geno$haplotype==hap)
#	}
#}
#rownames(hapMat) <- subjects
#colnames(hapMat) <- paste("H", gsub(":", "_", haplotypes), sep="_")
#save(hapMat, file="hapMat.RData")

#Get data
load(file="/data1/ccr5ju/AA_Immuno/hla_assoc/processed_impdata/hapMat.RData")

#Quality filtering
exclude = rownames(hapMat)[!rowSums(hapMat)==2]
hapMat_filt = hapMat[!rownames(hapMat)%in%exclude,] 

#Create data frame
df = data.frame(t1d=t1d_map[rownames(hapMat_filt),"case"]-1, pcs[rownames(hapMat_filt),c("PC1","PC2")], SEX=sex_map[rownames(hapMat_filt),"SEX"], hapMat_filt)

#Fix Denver samples --> set as controls
denver_fixlist = read.table("remove_t2d_bdc.txt", header=FALSE)
names(denver_fixlist) = c("FID", "IID", "reason")
df[row.names(df)%in%denver_fixlist$IID,"t1d"] <- 0


#Get haplotype info
getHapInfo = function(df){
	haps = names(df)[grep("H_", names(df))]
	count = NULL
	case_count = NULL
	case_freq = NULL
	case_percent = NULL
	control_count = NULL
	control_freq = NULL
	control_percent = NULL
	for (i in 1:length(haps)) {
		count[i] = sum(df[,haps[i]])
		case_count[i] = sum(df[df$t1d==1,haps[i]])
		control_count[i] = sum(df[df$t1d==0,haps[i]])
		case_freq[i] = sum(df[df$t1d==1,haps[i]])/(2*dim(df[df$t1d==1,])[1])
		control_freq[i] = sum(df[df$t1d==0,haps[i]])/(2*dim(df[df$t1d==0,])[1])
	}
	hap_info = data.frame(count, case_count, control_count, case_freq, control_freq); rownames(hap_info) <- haps
	hap_info$case_AF = with(hap_info, format(round(case_freq, 3),nsmall=1, trim=TRUE))
	hap_info$control_AF = with(hap_info, format(round(control_freq, 3),nsmall=1, trim=TRUE))
	hap_info$OR_unadjusted = round((hap_info$case_freq/(1-hap_info$case_freq))/(hap_info$control_freq/(1-hap_info$control_freq)), digits=2)
	hap_info$Erlich_case_AF = erlich[haps,"Proband"]/100
	hap_info$Erlich_control_AF = erlich[haps,"AFBAC"]/100
	hap_info$Erlich_OR = erlich[haps,"OR"]
	hap_info
}
hap_info = getHapInfo(df) #calculate freqs in all alleles

#Fit model
fitBiModel = function(df, hap_info, covariates) {
	haps_1 = names(df)[grep("H_", names(df))]
	haps = haps_1[!haps_1%in%covariates]
	out.glm.null = glm(as.formula(paste("t1d~",paste(covariates,collapse="+"))), family="binomial", data=df)
	
	getBetaPval = function(hap) {
		if (hap_info[hap,"count"] >= 30) {
			out.glm.full = glm(as.formula(paste("t1d~",paste(c(hap, covariates),collapse="+"))), family="binomial", data=df)
			beta = summary(out.glm.full)$coeff[hap,"Estimate"]
			or = format(exp(beta), digits=2, nsmall=2)
			lrt_p = format(lrtest(out.glm.full, out.glm.null)$"Pr(>Chisq)"[2], scientific=TRUE, digits=2)	
		} else {
			or = NA
			lrt_p = NA
		}
		c(or, lrt_p)
	}
	hap.p = as.data.frame(t(sapply(haps, getBetaPval))); names(hap.p) <- c("OR","P")
	#format model info output
	out.df=data.frame(hap_info[haps,], hap.p[haps,])	
	out.df$DRB1 =  sapply(strsplit(haps, split="_"),"[[", 2)
	out.df$DQA1 =  sapply(strsplit(haps, split="_"),"[[", 3)
	out.df$DQB1 =  sapply(strsplit(haps, split="_"),"[[", 4)	
	OUT.BI = out.df[order(out.df$DRB1, out.df$DQA1, out.df$DQB1, decreasing=FALSE),c("DRB1","DQA1","DQB1", "count", "control_AF", "case_AF", "OR_unadjusted","OR", "P", "Erlich_control_AF", "Erlich_case_AF", "Erlich_OR")]
	OUT.BI
} 
OUT.BI.ALL = fitBiModel(df, hap_info, covariates)
OUT.BI.COMMON = OUT.BI.ALL[OUT.BI.ALL$count>=30,]

write.table(OUT.BI.COMMON, file="haplotype_associations_common.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
write.table(OUT.BI.ALL, file="haplotype_associations_all.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


#####Conditional analysis
fitBiModel = function(df, hap_info, covariates) {
	haps_1 = names(df)[grep("H_", names(df))]
	haps = haps_1[!haps_1%in%covariates]
	out.glm.null = glm(as.formula(paste("t1d~",paste(covariates,collapse="+"))), family="binomial", data=df)
	
	getBetaPval = function(hap) {
		if (hap_info[hap,"count"] >= 30) {
			out.glm.full = glm(as.formula(paste("t1d~",paste(c(hap, covariates),collapse="+"))), family="binomial", data=df)
			beta = summary(out.glm.full)$coeff[hap,"Estimate"]
			or = format(exp(beta), digits=2, nsmall=2)
			lrt_p = lrtest(out.glm.full, out.glm.null)$"Pr(>Chisq)"[2]	
		} else {
			or = NA
			lrt_p = NA
		}
		c(or, lrt_p)
	}
	hap.p = as.data.frame(t(sapply(haps, getBetaPval))); names(hap.p) <- c("OR","P")
	#format model info output
	out.df=data.frame(hap_info[haps,], hap.p[haps,])	
	out.df$DRB1 =  sapply(strsplit(haps, split="_"),"[[", 2)
	out.df$DQA1 =  sapply(strsplit(haps, split="_"),"[[", 3)
	out.df$DQB1 =  sapply(strsplit(haps, split="_"),"[[", 4)	
	OUT.BI = out.df[order(out.df$DRB1, out.df$DQA1, out.df$DQB1, decreasing=FALSE),c("DRB1","DQA1","DQB1", "count", "control_AF", "case_AF", "OR_unadjusted","OR", "P", "Erlich_control_AF", "Erlich_case_AF", "Erlich_OR")]
	OUT.BI
} 
common_haps = row.names(hap_info[hap_info$count>30,])
common_df = df[,c("t1d","PC1","PC2","SEX",common_haps)]
p_thresh = 0.05/length(common_haps)
p_min = 0
round = 0
condition_haps = NULL
covars =  c("PC1", "PC2","SEX")
top_hap = NA
while(p_min<p_thresh) {
	cat(round, p_min, top_hap,"\n")
	out = fitBiModel(common_df, hap_info, covars)	
	p_min = min(as.numeric(out$P))
	top_hap = row.names(out)[out$P==p_min]	
	
	#update
	round = round+1
	condition_haps = c(condition_haps, top_hap)
	covars = c(covars, top_hap)
}
out_all = glm(as.formula(paste("t1d~",paste(covars[1:(length(covars)-1)],collapse="+"))), family="binomial", data=common_df)
OUT.COND = as.data.frame(summary(out_all)$coef[condition_haps[1:(length(condition_haps)-1)],c("Estimate", "Pr(>|z|)")])
OUT.COND$OR = round(exp(OUT.COND$Estimate),digits=2) 
OUT.COND$P = format(OUT.COND[,"Pr(>|z|)"], scientific=TRUE, digits=2)
rowname_strsplit = strsplit(rownames(OUT.COND), split="_")
OUT.COND$DRB1=unlist(lapply(rowname_strsplit,"[[", 2))
OUT.COND$DQA1=unlist(lapply(rowname_strsplit,"[[", 3))
OUT.COND$DQB1=unlist(lapply(rowname_strsplit,"[[", 4))



#Write tables to file
rtffile <- RTF("haplotype_assoc_tables.doc")

addParagraph(rtffile, "Supplementary Table ##. Association of T1D in African Americans with MHC Class II Haplotypes")
addTable(rtffile, OUT.BI.COMMON)
addParagraph(rtffile, "Only haplotypes with total allele count >30 are presented. Odds ratios were estimated and p-values generated by logistic regression, adjusting for 2 principal components and sex.\n\n")
#addTable(rtffile, OUT.BI.ALL)
addParagraph(rtffile, "\n\n")


addParagraph(rtffile, "Supplementary Table 7. Fourteen HLA class II haplotypes are independently associated with T1D in African Americans")
addTable(rtffile, OUT.COND[,c("DRB1","DQA1","DQB1","OR","P")])
addParagraph(rtffile, "Odds ratios and p-values generated by multivariable logistic regression of T1D risk, including 2 principal components, sex, and all 9 alleles as independent variables. Significance evaluated at a p<1.5x10-3 threshold (Bonferroni correction at family-wise error rate of 0.05 given ≤33 haplotypes tested at each round of conditional analysis).\n\n")



done(rtffile)
#cd /Users/Cassie/Box\ Sync/Rich\ Lab/AA_Immunochip/FixDenverSamples
#scp haplotype_assoc_tables.doc .










