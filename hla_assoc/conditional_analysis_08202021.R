options(stringsAsFactors=FALSE)
setwd("/data1/ccr5ju/AA_Immuno/hla_assoc")
#setwd("/data1/ccr5ju/AA_Immuno/thesis")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define conditional analysis functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
runConditionalAnalysis = function(infile, threshold=NULL) {

	### Read dosage data and add column names
	dose = read.table(infile)
	fam = read.table("../hla_imputation/imputation_results/AA_IMPUTED.fam", comment.char="~")
	names(fam) <- c("FID","IID", "FAT", "MOT", "SEX", "CASE") #case=2 control=1
	uniq_ids = fam$IID
	colnames(dose) <- c("Marker", "Allele1", "Allele2", uniq_ids)

	### Remove Denver samples
	denver_samples = read.table("../FixDenverSamples/remove_t2d_bdc.txt")
	names(denver_samples) <- c("FID", "IID", "V3")
	uniq_ids_keep = uniq_ids[!uniq_ids %in% denver_samples$IID]
	dose = dose[,c("Marker", "Allele1", "Allele2", uniq_ids_keep)]

	### Create pheno data.frame
	pcs = read.table("../raw/pc10.txt", header=TRUE, comment.char="~")
	pheno = merge(fam, pcs, by="IID")
	row.names(pheno) = pheno$IID

	### Filter for common alleles
	dose$counts = rowSums(dose[,uniq_ids_keep])
	dose$AF = dose$counts/(2*length(uniq_ids_keep))
	dose$MAF = sapply(dose$AF, function(x){min(x, 1-x)})
	dose_common = dose[dose$MAF>0.025,]

	if (is.null(threshold)) {
		threshold=0.05/nrow(dose_common)
	}

	getGLMStats = function(glm_obj, ind_var) {
		coefs = coef(summary(glm_obj))
		est = coefs[ind_var,"Estimate"]
		or = exp(est)
		p = coefs[ind_var,4]
		c(est, or, p)
	}

	runAssoc = function(x, y) { #<<<<---- INSERT COVARIATES INTO MODEL
		indVarsMat = cbind(x, cov)
		out_glm = glm(y~indVarsMat, family="binomial")
		case_freq = sum(x[y==1])/(2*length(x[y==1]))
		control_freq = sum(x[y==0])/(2*length(x[y==0]))
		freq = sum(x)/(2*length(x))
		c(case_freq, control_freq, freq, getGLMStats(out_glm, "indVarsMatx"))
	}

	runCondAssoc = function(x, y, condition_mat) { #<<<<---- INSERT COVARIATES INTO MODEL
		X_new = cbind(x, condition_mat)
		indVarsMat = cbind(X_new, cov)
		out_glm = glm(y~indVarsMat, family="binomial")
		case_freq = sum(x[y==1])/(2*length(x[y==1]))
		control_freq = sum(x[y==0])/(2*length(x[y==0]))
		freq = sum(x)/(2*length(x))
		c(case_freq, control_freq, freq, getGLMStats(out_glm, "indVarsMatx"))
	}

	getInfo = function(assoc_mat) {
		min_p = min(assoc_mat[,"P"], na.rm=TRUE)
		topVar = rownames(assoc_mat)[assoc_mat[,"P"]==min_p]
		list(topVar, assoc_mat[topVar,])
	}

	#Define input matrices
	X <- t(as.matrix(dose_common[,uniq_ids_keep])); colnames(X) <- dose_common[,1]; rownames(X) <- uniq_ids_keep
	y <- pheno[uniq_ids_keep,]$CASE-1
	cov <- as.matrix(pheno[uniq_ids_keep, c("SEX","PC1","PC2")])

	#Initialize loop
	min_p=0
	condition_list_all=c()
	condition_list_pruned=c()
	output = list()

	round=1
	while (min_p<threshold) {
		cat(round,"\n")

		#run model
		if (length((condition_list_all))==0) {
			out = t(apply(X=X, FUN=runAssoc, MARGIN=2, y))
		} else {
			out = t(apply(X=X[,!colnames(X)%in%condition_list_all], FUN=runCondAssoc, MARGIN=2, y, X[,condition_list_pruned]))
		}
		colnames(out) <- c("case_AF", "control_AF", "total_AF", "Est","OR", "P")
		#export association results for each round
		output[[round]] = list(round=round, condition_snps=condition_list_all, strongest_variant=getInfo(out)[[1]], assoc_stats=getInfo(out)[[2]], out, threshold=threshold)

		#update
		if (class(output[[round]]$assoc_stats)=="numeric"){
			min_p = output[[round]]$assoc_stats["P"]
		} else if (class(output[[round]]$assoc_stats)=="matrix") {
			min_p = output[[round]]$assoc_stats[,"P"][1]
		}
		condition_list_all = c(condition_list_all, output[[round]]$strongest_variant)
		condition_list_pruned = c(condition_list_pruned, output[[round]]$strongest_variant[1])
		round = round+1
	}
	output
}

## NOTE There was a bug in this function (it assumed there was only one top variant per round).
## It did not affect analysis in the Diabetes Care 2018 manuscript. But old versions of the file conditional_analysis_aa.out
## will have some errors. This was fixed as of August 2021 (while assembling tables for my thesis)
getSummaryTable = function(output, outfile) {
	print_out = data.frame(unlist(lapply(output, '[[', 3)), do.call("rbind", lapply(output, '[[', 4)))
	names(print_out) <- c("VAR", "case_AF", "control_AF", "total_AF", "Est","OR", "P")
	write.table(print_out, file=outfile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
	print(print_out)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run conditional analyses
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### 2-digit HLA alleles (dosages)
out.2d =runConditionalAnalysis(infile = "processed_impdata/dosages.2d")
saveRDS(out.2d, file="assoc/conditional_analysis_object.2d.Rdata")
getSummaryTable(out.2d, outfile="assoc/conditional_analysis_2d.out")

### 4-digit HLA alleles (dosages)
out.4d =runConditionalAnalysis(infile = "processed_impdata/dosages.4d")
saveRDS(out.4d, file="assoc/conditional_analysis_object.4d.Rdata")
getSummaryTable(out.4d, outfile="assoc/conditional_analysis_4d.out")


### Amino acids (dosages) - #Description of amino acid naming conventions: http://hla.alleles.org/alleles/formats.html
out.aa =runConditionalAnalysis(infile = "processed_impdata/dosages.aa")
save(out.aa, file="assoc/conditional_analysis_object.aa.RData")
getSummaryTable(out.aa, outfile="assoc/conditional_analysis_aa.out")

### Everything combined (HLA alleles, amino acids, and SNPs)
out.all =runConditionalAnalysis(infile = "processed_impdata/dosages.all")
saveRDS(out.all, file="assoc/conditional_analysis_object.all.Rdata")
getSummaryTable(out.all, outfile="assoc/conditional_analysis_all.out")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Explore AA marker diversity
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fam = read.table("../hla_imputation/imputation_results/AA_IMPUTED.fam", comment.char="~")
names(fam) <- c("FID","IID", "FAT", "MOT", "SEX", "CASE") #case=2 control=1
uniq_ids = paste(fam$FID, fam$IID, sep="_")

bim = read.table("../hla_imputation/imputation_results/AA_IMPUTED.bim", comment.char="~")
names(bim) <- c("chr", "Marker", "cm", "Position", "Allele1", "Allele2")

dose.aa = read.table("processed_impdata/dosages.aa")
colnames(dose.aa) <- c("Marker", "Allele1", "Allele2", uniq_ids)
dose.aa = merge(bim[,c("Marker","Position")], dose.aa[,c("Marker", uniq_ids)], by=c("Marker"))

dose.aa$counts = rowSums(dose.aa[,uniq_ids])
dose.aa$AF = dose.aa$counts/(2*length(uniq_ids))
dose.aa$MAF = sapply(dose.aa$AF, function(x){min(x, 1-x)})
sum(dose.aa$MAF>0.05)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cluster markers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#library(ggplot2, lib.loc="/mnt/t1/users/ccr5ju/software/R_libraries/x86_64-redhat-linux-gnu-library/3.3")
library(ggplot2)

fam = read.table("../hla_imputation/imputation_results/AA_IMPUTED.fam", comment.char="~")
names(fam) <- c("FID","IID", "FAT", "MOT", "SEX", "CASE") #case=2 control=1
uniq_ids = paste(fam$FID, fam$IID, sep="_")

bim = read.table("../hla_imputation/imputation_results/AA_IMPUTED.bim", comment.char="~")
names(bim) <- c("chr", "Marker", "cm", "Position", "Allele1", "Allele2")

dose = read.table("processed_impdata/dosages.all")
colnames(dose) <- c("Marker", "Allele1", "Allele2", uniq_ids)
dose = merge(bim[,c("Marker","Position")], dose[,c("Marker", uniq_ids)], by=c("Marker"))

doseLong = melt(dose[,c("Marker","Position", uniq_ids)], id.vars=c("Marker","Position"), variable.name="Subject", value.name="Dosage")
doseLong$Marker_ordered = factor(doseLong$Marker, levels=dose$Marker[order(dose$Position)])

png("hla_heatmap.png")
ggplot(doseLong[sample(x=seq(1,nrow(doseLong)),size=1000),]) +
	geom_tile(aes(x=Marker_ordered, y=Subject, fill=Dosage)) +
	theme(axis.text = element_blank(),
				axis.ticks = element_blank())
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot Amino Acid associations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2, lib.loc="/mnt/t1/users/ccr5ju/software/R_libraries/x86_64-redhat-linux-gnu-library/3.3")


getData = function(d) {
	varsplit = strsplit(row.names(d), split="_")
	d$gene = sapply(varsplit,"[[", 2)
	d$Gene = factor(d$gene, levels=c("A","C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"), labels=c("HLA-A","HLA-C", "HLA-B", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1"))
	d$aa_pos = sapply(varsplit,"[[", 3)
	d$aa_pos_numeric=d$aa_pos
	d$aa_pos_numeric[grep("x", d$aa_pos)] = sapply(strsplit(d$aa_pos[grep("x", d$aa_pos)],split="x"),"[[",1)
	d$aa_pos_numeric = as.numeric(d$aa_pos_numeric)
	#d$bp_pos = as.numeric(sapply(varsplit,"[[", 4))
	d
}


##### PLOTS FOR POSTER
makePlots = function(df, title, xlab) {
	axisstyle.y <- element_text(face = "bold", size = 15)
	axisstyle.x <- element_text(size = 12)
	axislabelstyle.x =element_text(size=16)
	ggplot(df, aes(x=aa_pos_numeric, y=-log10(P))) + geom_point() + facet_grid(.~Gene) +
			ggtitle(title) + ylab("-log10(p)") + xlab(xlab)+ #xlab("Amino acid position")
			theme(axis.text.y = axisstyle.y, axis.text.x = axisstyle.x, axis.title.x=axislabelstyle.x) +
			theme(aspect.ratio=1/2) +
			theme(strip.text.x = element_text(size = 12, face="bold"))
}

### Round 1 plot
pdf("assoc/amino_acid_conditional_analysis_round1.pdf")
	d = getData(as.data.frame(out.aa[[1]][[5]]))
	makePlots(d[d$gene%in%c("DRB1", "DQA1", "DQB1"),], title="", xlab="")
dev.off()

### Round 2 plot
pdf("assoc/amino_acid_conditional_analysis_round2.pdf")
	d = getData(as.data.frame(out.aa[[2]][[5]]))
	makePlots(d[d$gene%in%c("DRB1", "DQA1", "DQB1"),], title="",xlab="")
dev.off()

### Round 3 plot
pdf("assoc/amino_acid_conditional_analysis_round3.pdf")
	d = getData(as.data.frame(out.aa[[3]][[5]]))
	makePlots(d[d$gene%in%c("DRB1", "DQA1", "DQB1"),], title="", xlab="Amino acid position")
dev.off()

#scp assoc/amino_acid_conditional_analysis_round*.pdf "$local/AA_Immunochip"



##### PLOTS FOR EXPLORATORY ANALYSIS
makePlots2 = function(df, title) {
	axisstyle.y <- element_text(face = "bold", size = 8)
	axisstyle.x <- element_text(size = 8)
	ggplot(df, aes(x=aa_pos_numeric, y=-log10(P))) + geom_point() + facet_grid(.~Gene) +
			ggtitle(title) + ylab("-log10(p)") + xlab("")+ #xlab("Amino acid position")
			theme(axis.text.y = axisstyle.y, axis.text.x = axisstyle.x) +
			theme(aspect.ratio=1/2) +
			theme(strip.text.x = element_text(size = 12, face="bold"))
}

#All rounds, all genes
pdf("assoc/amino_acid_conditional_analysis_all.pdf")
for (i in 1:length(out.aa)) {
#for (i in 1:5) {
	d = getData(as.data.frame(out.aa[[i]][[5]]))
	p = makePlots2(d, title="")
	print(p)
}
dev.off()
#scp assoc/amino_acid_conditional_analysis_all.pdf "$local/AA_Immunochip"


#associations with remaining variants
dose = read.table("processed_impdata/dosages.all")
aa_condition_list = sumtable$VAR[!sumtable$VAR%in%c("AA_DQA1_52_32717206_H","AA_DQA1_54_32717212")]

fam = read.table("../hla_imputation/imputation_results/AA_IMPUTED.fam", comment.char="~")
names(fam) <- c("FID","IID", "FAT", "MOT", "SEX", "CASE") #case=2 control=1
uniq_ids = paste(fam$FID, fam$IID, sep="_")

X <- t(as.matrix(dose[,4:dim(dose)[2]])); colnames(X) <- dose[,1]; rownames(X) <- uniq_ids
y <- fam$CASE-1

runCondAssoc = function(x, y, condition_mat) {
	if (max(cor(x, condition_mat))<0.95) {
		X_new = cbind(x, condition_mat)
		out_glm = glm(y~X_new, family="binomial")
		case_freq = mean(x[y==1])
		control_freq = mean(x[y==0])
		freq = mean(x)
		out=c(case_freq, control_freq, freq, getGLMStats(out_glm, "X_newx"))
	} else {
		case_freq = mean(x[y==1])
		control_freq = mean(x[y==0])
		freq = mean(x)
		out=c(case_freq, control_freq, freq, NA, NA, NA)
	}

}
getGLMStats = function(glm_obj, ind_var) {
	coefs = coef(summary(glm_obj))
	est = coefs[ind_var,"Estimate"]
	or = exp(est)
	p = coefs[ind_var,4]
	c(est, or, p)
}
out = as.data.frame(t(apply(X=X[,!colnames(X)%in%aa_condition_list], FUN=runCondAssoc, MARGIN=2, y, X[,aa_condition_list])))
colnames(out) <- c("case_AF", "control_AF", "total_AF", "Est","OR", "P")
out$VAR = row.names(out)
write.table(out, "assoc/conditioned_on_all_aa.txt")
#sed 's/"//g' conditioned_on_all_aa.txt | awk 'NR==1 {print "SNP",$0} NR>1 {print}' > conditioned_on_all_aa_errorsRemoved.txt
#scp conditioned_on_all_aa_errorsRemoved.txt "$local/AA_Immunochip"

out_sorted = out[order(out$P),]
